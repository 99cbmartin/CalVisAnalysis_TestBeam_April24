#include <TFile.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>

TGraph* profileToGraph(TProfile* prof) {
    if (!prof) return nullptr;
    TGraph* graph = new TGraph();
    for (int i = 1; i <= prof->GetNbinsX()-24; ++i) {
        graph->SetPoint(i-1, prof->GetBinCenter(i), prof->GetBinContent(i));
    }
    return graph;
}

double getMaxValue(TGraph* graph) {
    double max = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < graph->GetN(); ++i) {
        double x, y;
        graph->GetPoint(i, x, y);
        if (y > max) max = y;
    }
    return max;
}

void drawGraphsOnCanvas(TCanvas* canvas, TGraph* graph1, TGraph* graph2, TGraph* graph3, TGraph* graph4, Color_t* colors, const char* filename, int canvasNum) {
    canvas->cd();
    
    TGraph* graphs[4] = {graph1, graph2, graph3, graph4};
    double maxValues[4] = {
        graph1 ? getMaxValue(graph1) : -std::numeric_limits<double>::infinity(),
        graph2 ? getMaxValue(graph2) : -std::numeric_limits<double>::infinity(),
        graph3 ? getMaxValue(graph3) : -std::numeric_limits<double>::infinity(),
        graph4 ? getMaxValue(graph4) : -std::numeric_limits<double>::infinity()
    };
    
    for (int i = 0; i < 4; ++i) {
        for (int j = i + 1; j < 4; ++j) {
            if (maxValues[i] < maxValues[j]) {
                std::swap(maxValues[i], maxValues[j]);
                std::swap(graphs[i], graphs[j]);
            }
        }
    }
    
    for (int i = 0; i < 4; ++i) {
        if (graphs[i]) {
            graphs[i]->SetLineColor(colors[i]);
            graphs[i]->SetLineStyle(i == 0 ? 1 : 2);  
            graphs[i]->Draw(i == 0 ? "AL" : "L SAME");
        }
    }
    
    TString canvasName = TString::Format("%s_Canvas%d.png", filename, canvasNum);
    canvas->SaveAs(canvasName);
    std::cout << "Saved canvas: " << canvasName << std::endl;
}

void processFile(const char* filename) {
    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    Color_t colors[] = {kRed, kBlue, kGreen, kMagenta};

    TGraph* graphs1[4];
    for (int i = 0; i < 4; ++i) {
        TProfile* profile = dynamic_cast<TProfile*>(file->Get(TString::Format("AveragePulse_Ch%d", i)));
        graphs1[i] = profile ? profileToGraph(profile) : nullptr;
    }

    TGraph* graphs2[2];
    for (int i = 0; i < 2; ++i) {
        int ch = (i == 0) ? 4 : 7;
        TProfile* profile = dynamic_cast<TProfile*>(file->Get(TString::Format("AveragePulse_Ch%d", ch)));
        graphs2[i] = profile ? profileToGraph(profile) : nullptr;
    }

    TGraph* graphs3[6];
    for (int i = 0; i < 4; ++i) {
        graphs3[i] = graphs1[i];
    }
    for (int i = 0; i < 2; ++i) {
        graphs3[4 + i] = graphs2[i];
    }

    TCanvas* canvas1 = new TCanvas("canvas1", "AveragePulse_Ch0-Ch3", 800, 600);
    drawGraphsOnCanvas(canvas1, graphs1[0], graphs1[1], graphs1[2], graphs1[3], colors, filename, 1);
    delete canvas1;

    TCanvas* canvas2 = new TCanvas("canvas2", "AveragePulse_Ch4, Ch7", 800, 600);
    drawGraphsOnCanvas(canvas2, graphs2[0], graphs2[1], nullptr, nullptr, colors, filename, 2);
    delete canvas2;

    TCanvas* canvas3 = new TCanvas("canvas3", "AveragePulse_Ch0-Ch3, Ch4, Ch7", 800, 600);
    drawGraphsOnCanvas(canvas3, graphs3[0], graphs3[1], graphs3[2], graphs3[3], colors, filename, 3);
    drawGraphsOnCanvas(canvas3, graphs3[4], graphs3[5], nullptr, nullptr, colors, filename, 3);  
    delete canvas3;

    for (int i = 0; i < 4; ++i) delete graphs1[i];
    for (int i = 0; i < 2; ++i) delete graphs2[i];

    file->Close();
    delete file;
}

int main() {
    const char* filename = "testbeam_plots/June/week2/run582_LG.root";
    processFile(filename);
    return 0;
}


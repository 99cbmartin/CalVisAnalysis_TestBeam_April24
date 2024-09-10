#include <TFile.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <algorithm>
#include <vector>

// Function to convert a TProfile to a TGraph
TGraph* profileToGraph(TProfile* prof) {
    if (!prof) return nullptr;
    TGraph* graph = new TGraph();
    for (int i = 1; i <= prof->GetNbinsX()-24; ++i) {
        graph->SetPoint(i-1, prof->GetBinCenter(i), prof->GetBinContent(i));
    }
    return graph;
}

// Function to get the maximum value of a TGraph
double getMaxValue(TGraph* graph) {
    double max = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < graph->GetN(); ++i) {
        double x, y;
        graph->GetPoint(i, x, y);
        if (fabs(y) > max) max =fabs(y);
    }
    return max;
}

void drawGraphsOnCanvas(TCanvas* canvas, const std::vector<TGraph*>& graphs, const std::vector<Color_t>& colors) {
    canvas->cd();
    std::vector<std::pair<TGraph*, double>> graphMaxPairs;
    for (size_t i = 0; i < graphs.size(); ++i) {
        graphMaxPairs.emplace_back(graphs[i], getMaxValue(graphs[i]));
    }

    std::sort(graphMaxPairs.begin(), graphMaxPairs.end(), [](const auto& a, const auto& b) {
        return a.second > b.second;
    });

    for (size_t i = 0; i < graphMaxPairs.size(); ++i) {
        graphMaxPairs[i].first->SetLineColor(colors[i]);
        graphMaxPairs[i].first->Draw(i == 0 ? "AL" : "L SAME");
    }
}

void processFile(const char* filename) {
    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    Color_t colors[] = {kRed, kBlue, kGreen, kMagenta, kCyan, kYellow};

    std::vector<TGraph*> graphs1, graphs2, graphs3;

    int channels1[] = {0, 1, 2, 3};
    int channels2[] = {4, 7};

    for (int ch : channels1) {
        TProfile* profile = dynamic_cast<TProfile*>(file->Get(TString::Format("AveragePulse_Ch%d", ch)));
        if (profile) {
            TGraph* graph = profileToGraph(profile);
            if (graph) {
                graphs1.push_back(graph);
                graphs3.push_back(graph);
            }
        }
    }

    for (int ch : channels2) {
        TProfile* profile = dynamic_cast<TProfile*>(file->Get(TString::Format("AveragePulse_Ch%d", ch)));
        if (profile) {
            TGraph* graph = profileToGraph(profile);
            if (graph) {
                graphs2.push_back(graph);
                graphs3.push_back(graph);
            }
        }
    }

    TCanvas* canvas1 = new TCanvas("canvas1", "AveragePulse_Ch0-Ch3", 800, 600);
    drawGraphsOnCanvas(canvas1, graphs1, std::vector<Color_t>(colors, colors + graphs1.size()));
    canvas1->SaveAs(TString::Format("%s_Canvas1.png", filename));

    TCanvas* canvas2 = new TCanvas("canvas2", "AveragePulse_Ch4, Ch7", 800, 600);
    drawGraphsOnCanvas(canvas2, graphs2, std::vector<Color_t>(colors, colors + graphs2.size()));
    canvas2->SaveAs(TString::Format("%s_Canvas2.png", filename));

    TCanvas* canvas3 = new TCanvas("canvas3", "AveragePulse_Ch0-Ch3, Ch4, Ch7", 800, 600);
    drawGraphsOnCanvas(canvas3, graphs3, std::vector<Color_t>(colors, colors + graphs3.size()));
    canvas3->SaveAs(TString::Format("%s_Canvas3.png", filename));

    delete canvas1;
    delete canvas2;
    delete canvas3;
    for (auto graph : graphs1) delete graph;
    for (auto graph : graphs2) delete graph;
    for (auto graph : graphs3) delete graph;

    file->Close();
    delete file;
}

int main() {
    const char* filename = "testbeam_plots/June/week2/run582_LG.root";
    processFile(filename);
    return 0;
}


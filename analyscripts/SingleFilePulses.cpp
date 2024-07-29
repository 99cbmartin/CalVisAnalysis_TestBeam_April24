//this code should plots multiple average pulses from different channels in a single file and order them from largest to smallest for plotting axes to be well adjusted
#include <TFile.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>

TGraph* Convert(TProfile* prof) {
    if (!prof) return nullptr;
    TGraph* graph = new TGraph();
    for (int i = 1; i <= prof->GetNbinsX()-24; ++i) {//truncate the last 24 ns since electronics produce spikes that confuse getMaxValue
        graph->SetPoint(i-1, prof->GetBinCenter(i), prof->GetBinContent(i));
    }
    return graph;
}

double getMaxValue(TGraph* graph) {
    double max = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < graph->GetN(); ++i) {
        double x, y;
        graph->GetPoint(i, x, y);
        if (fabs(y) > max) max = fabs(y);
    }
    return max;
}

void draw(const char* filename, int* channels, int numChannels) {
    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    Color_t colors[] ={kRed, kBlue, kGreen, kMagenta, kBlack, kOrange};
    if (numChannels > sizeof(colors)/sizeof(colors[0])) {
        std::cerr << "Error: More channels than colors available" << std::endl;
        return;
    }

    TGraph* graphs[numChannels];
    double maxValues[numChannels];
    for (int i = 0; i < numChannels; ++i) {
        TProfile* profile = dynamic_cast<TProfile*>(file->Get(TString::Format("AveragePulse_Ch%d", channels[i])));
        graphs[i] = profile ? Convert(profile) : nullptr;
        maxValues[i] = graphs[i] ? getMaxValue(graphs[i]) : -std::numeric_limits<double>::infinity();
    }

//    for (int i = 0; i < numChannels; ++i) {
//        for (int j = i + 1; j < numChannels; ++j) {
//            if (maxValues[i] < maxValues[j]) {
//                std::swap(maxValues[i], maxValues[j]);
//                std::swap(graphs[i], graphs[j]);
//                std::swap(channels[i], channels[j]);
 //           }
 //       }
 //   }

    TCanvas* canvas = new TCanvas("canvas", "Average Pulses", 800, 600);
    for (int i = 0; i < numChannels; ++i) {
        if (graphs[i]) {
            graphs[i]->SetLineColor(colors[i]);
            graphs[i]->SetLineStyle(i == 0 ? 1 : 2); 
            graphs[i]->SetTitle(TString::Format("Channel %d", channels[i]));
            if (i == 0) {
                graphs[i]->GetYaxis()->SetRangeUser(-1.2,0.2);
		graphs[i]->GetXaxis()->SetTitle("t-T_{trig} (ns)");
                graphs[i]->GetYaxis()->SetTitle("Amp (mV)");
                graphs[i]->SetTitle("Norm Average Pulses Back-End Channels");
                graphs[i]->SetStats(0);
		graphs[i]->Scale(1/maxValues[i]);
                graphs[i]->Draw("AL");
            } else {
		graphs[i]->Scale(1/maxValues[i]);
                graphs[i]->Draw("L SAME");
            }
        }
    }
    TString canvasName = TString::Format("%s_NormBackPulses.png", filename);
    canvas->SaveAs(canvasName);
    std::cout << "Saved canvas: " << canvasName << std::endl;

    for (int i = 0; i < numChannels; ++i) delete graphs[i];
    delete canvas;
    file->Close();
    delete file;
}

int main() {
    const char* filename = "testbeam_plots/June/week3/run526_LG.root";
    int channels[] = {4,5,6,7};  // Specify the channels to plot
    int numChannels = sizeof(channels) / sizeof(channels[0]);
    draw(filename, channels, numChannels);
    return 0;
}


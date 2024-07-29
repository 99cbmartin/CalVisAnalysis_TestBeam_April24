#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <iostream>
#include <TAxis.h>
int main() {
    TFile* file = TFile::Open("../testbeam_plots/July/week4/run284_LG.root", "READ");
    if (!file || !file->IsOpen()) {
        std::cerr << "Error: Could not open run594_LG.root file!" << std::endl;
        return 1;
    }
for (int R = 1 ; R < 5 ; R++){
    for (int i = 0; i < 8; ++i) {
        TString graphName = Form("R%dAveragePulseGraph_Ch%d",R,i);
        
        TGraph* graph = (TGraph*)file->Get(graphName);
        if (!graph) {
            std::cerr << "Error: Could not find graph " << graphName << " in file!" << std::endl;
            continue;
        }

        TCanvas* canvas = new TCanvas("canvas", "Canvas", 800, 600);
       	graph->SetLineColor(kRed);
	graph->GetXaxis()->SetTitle("t-t_{trig}} (ns)");
	graph->GetYaxis()->SetTitle("amp (mV)");
	graph->SetTitle(Form("Region %d Average Pulse Ch_%d",R,i));
	 graph->Draw("AL");

        TString pngFileName = Form("../testbeam_plots/July/week4/SDLFVU_Regions_2/R%d_Unfil%s.png",R,graphName.Data());
        canvas->Print(pngFileName);

        delete canvas;
    }
}
    file->Close();
    delete file;

    return 0;
}


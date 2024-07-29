
//This script will compare the average pulses of each SiPM to it's corresponding partner on the other end of the crystal

#include <iostream>
#include <TFile.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TString.h>
#include <filesystem>
#include <string>



namespace fs = std::filesystem;

 void ensureDirectoryExists(const std::string& path) {
     if (!fs::exists(path)) {
             fs::create_directories(path);
                 }
 }

double getMaxValue(TGraph* graph) {
    if (!graph) return -std::numeric_limits<double>::infinity();
    
    double maxVal = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < graph->GetN(); ++i) {
        double x, y;
        graph->GetPoint(i, x, y);
        if (fabs(y) > maxVal) {
            maxVal =fabs(y);
        }
    }
    return maxVal;
}


void processFile(const char* path, const char* filename, const char* ftype, const char* canvasNamePrefix, const char* outputpath) {
    TString fullFile = TString::Format("%s%s%s",path,filename,ftype);
    TFile* file = TFile::Open(fullFile);
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    auto profileToGraph = [](TProfile* prof) -> TGraph* {
        if (!prof) return nullptr;
        TGraph* graph = new TGraph();
        for (int i = 1; i <= prof->GetNbinsX()-24; ++i) {
            graph->SetPoint(i-1, prof->GetBinCenter(i), prof->GetBinContent(i));
        }
        return graph;
    };

    // Array of channel pairs for comparison
    int channelPairs[4][2] = {{0, 5}, {3, 6}, {1, 4}, {2, 7}};

    for (int i = 0; i < 4; ++i) {
        int ch1 = channelPairs[i][0];
        int ch2 = channelPairs[i][1];

        TString profName1 = TString::Format("AveragePulse_Ch%d", ch1);
        TString profName2 = TString::Format("AveragePulse_Ch%d", ch2);

        TProfile* prof1 = dynamic_cast<TProfile*>(file->Get(profName1));
        TProfile* prof2 = dynamic_cast<TProfile*>(file->Get(profName2));

        if (!prof1 || !prof2) {
            std::cerr << "Error retrieving profiles: " << profName1 << " or " << profName2 << std::endl;
            continue;
        }

        TGraph* graph1 = profileToGraph(prof1);
        TGraph* graph2 = profileToGraph(prof2);

        TCanvas* canvas = new TCanvas(TString::Format("%s_%d_vs_%d", canvasNamePrefix, ch1, ch2),
                                      TString::Format("%s: Channel %d vs Channel %d", filename, ch1, ch2), 800, 600);

        double max1 = getMaxValue(graph1);
        double max2 = getMaxValue(graph2);
	
	
	if (max1 > max2) {
            graph1->SetLineColor(kRed);
	    graph1->Scale(1/max1);
	    graph1->GetXaxis()->SetTitle("t-T_{trig} (ns)");
	    graph1->GetYaxis()->SetTitle("Amp (mV)");
	    graph1->SetTitle(Form("Average Pulse Channel %d Vs. Channel %d",channelPairs[i][0],channelPairs[i][1]));
	    graph1->SetStats(0);
	    graph1->GetYaxis()->SetRangeUser(-1.5,0.5);
            graph1->Draw("AL");
	    graph2->SetStats(0);
            graph2->SetLineColor(kBlue);
	    graph2->Scale(1/max2);
            graph2->Draw("L SAME");
        } else {
            graph2->SetLineColor(kBlue);
	    graph2->GetXaxis()->SetTitle("t-T_{trig} (ns)");
            graph2->GetYaxis()->SetTitle("Amp (mV)");
            graph2->SetTitle(Form("Average Pulse Channel %d Vs. Channel %d",channelPairs[i][0],channelPairs[i][1]));
	    graph2->SetStats(0);
	    graph2->GetYaxis()->SetRangeUser(-1.5,0.5);
	    graph2->Scale(1/max2);
            graph2->Draw("AL");
            graph1->SetLineColor(kRed);
	    graph1->SetStats(0);
	    graph1->Scale(1/max1);
            graph1->Draw("L SAME");
        }
	
       // graph1->SetLineColor(kRed);
//	graph1->SetStats(0);
  //      graph1->Draw("AL");
//	graph2->SetStats(0);
  //      graph2->SetLineColor(kBlue);
    //    graph2->Draw("L SAME");

        TString canvasName = TString::Format("%s%s_%s_%d_vs_%d.png",outputpath, canvasNamePrefix, filename, ch1, ch2);
        canvas->SaveAs(canvasName);

        delete graph1;
        delete graph2;
        delete canvas;
    }

    file->Close();
    delete file;
}

int main() {
    const char* filenames[] = {"run284_LG", "run290_LG", "run594_LG", "run582_LG"};
    const char* ftype = ".root";
    const char* canvasNamePrefix = "Norm_FrontVsBack";
    const char* path = "testbeam_plots/June/week2/";
    const char* savepath = "plots/";
    TString outputpath  = TString::Format("%s%s",path,savepath); 
    ensureDirectoryExists(outputpath.Data());
    for (const char* filename : filenames) {
        processFile(path, filename, ftype, canvasNamePrefix, outputpath);
    }

    return 0;
}


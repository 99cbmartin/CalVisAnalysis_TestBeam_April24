//this script will compare the Filtered to Unfiltered data on each SiPM for a given configuration of a crystal at a specific angle for HG and LG

#include <TProfile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TString.h>
#include <iostream>
#include <filesystem>

namespace fs = std::filesystem;

void ensureDirectoryExists(const std::string& path) {
    if (!fs::exists(path)) {
        fs::create_directories(path);
    }
}

TGraph* profileToGraph(TProfile* prof) {
    if (!prof) return nullptr;
    TGraph* graph = new TGraph();
    for (int i = 1; i <= prof->GetNbinsX()-24; ++i) {
        graph->SetPoint(i-1, prof->GetBinCenter(i), prof->GetBinContent(i));
    }
    return graph;
}

double getMaxValue(TGraph* graph) {
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

void compareFilesForChannel(const char* file1, const char* file2, int channel, const char* outputDir, const char* sourceDir, Color_t color1, Color_t color2){
    TString filename1 = TString::Format("%s%s.root", sourceDir,file1);
    TString filename2 = TString::Format("%s%s.root",sourceDir, file2);

    TFile* f1 = TFile::Open(filename1);
    TFile* f2 = TFile::Open(filename2);

    if (!f1 || f1->IsZombie() || !f2 || f2->IsZombie()) {
        std::cerr << "Error opening files: " << filename1 << " or " << filename2 << std::endl;
        return;
    }

    TString profName1 = TString::Format("AveragePulse_Ch%d", channel);
    TString profName2 = TString::Format("AveragePulse_Ch%d", channel);

    TProfile* prof1 = dynamic_cast<TProfile*>(f1->Get(profName1));
    TProfile* prof2 = dynamic_cast<TProfile*>(f2->Get(profName2));

    if (!prof1 || !prof2) {
        std::cerr << "Error retrieving profiles: " << profName1 << " or " << profName2 << std::endl;
        f1->Close();
        f2->Close();
        return;
    }

    TGraph* graph1 = profileToGraph(prof1);
    TGraph* graph2 = profileToGraph(prof2);

    if (!graph1 || !graph2) {
        std::cerr << "Error converting profiles to graphs: " << profName1 << " or " << profName2 << std::endl;
        f1->Close();
        f2->Close();
        delete graph1;
        delete graph2;
        return;
    }

    double max1 = getMaxValue(graph1);
    double max2 = getMaxValue(graph2);

    TCanvas* canvas = new TCanvas(TString::Format("AvgPulse0V90_%s_Channel%d", file1, channel),
                                  TString::Format("Average Pulse Channel %d: %s vs %s", channel, file1, file2), 800, 600);

    if (max1 > max2) {
        graph1->SetLineColor(color1); 
        graph1->SetLineStyle(1);      
	graph1->GetXaxis()->SetTitle("t-T_{trig} (ns)");
        graph1->GetYaxis()->SetTitle("Amp (mV)");
        graph1->SetTitle(Form("Average Pulse 0 Vs. 90 Channel %d",channel));
        graph1->SetStats(0);
	graph1->GetYaxis()->SetRangeUser(-1.5,0.5);
     	graph1->Scale(1/max1);
	graph1->Draw("AL");
        graph2->Scale(1/max2);
        graph2->SetLineColor(color2);
        graph2->SetLineStyle(2);      
        graph2->Draw("L SAME");
    } else {
        graph2->SetLineColor(color2);   
        graph2->SetLineStyle(2);       
        graph2->GetXaxis()->SetTitle("t-T_{trig} (ns)");
        graph2->GetYaxis()->SetTitle("Amp (mV)");
        graph2->SetTitle(Form("Average Pulse 0 Vs. 90 Channel %d",channel));
        graph2->SetStats(0);
	graph2->GetYaxis()->SetRangeUser(-1.5,0.5);
	graph2->Scale(1/max2);
	graph2->Draw("AL");
        graph1->Scale(1/max1);
        graph1->SetLineColor(color1);  
        graph1->SetLineStyle(1);      
        graph1->Draw("L SAME");
    }

    TString canvasName = TString::Format("%s/Norm_AvgPulse0V90_%s_Channel%d.png", outputDir, file1, channel);
    canvas->SaveAs(canvasName);

    std::cout << "Saved canvas: " << canvasName << std::endl;

    delete graph1;
    delete graph2;
    delete canvas;
    
    f1->Close();
    f2->Close();
}

int main() {
    const char* files1[] = {"run582_LG", "run290_LG"};
    const char* files2[] = {"run594_LG", "run284_LG"};
    const char* outputDir = "testbeam_plots/June/week2/plots";
    const char* sourceDir = "testbeam_plots/June/week2/";
    Color_t col[4] ={kGreen,kRed,kMagenta,kBlue}; 
    ensureDirectoryExists(outputDir);

    for (int channel = 0; channel < 8; ++channel) {
        compareFilesForChannel(files1[0], files2[0], channel, outputDir, sourceDir,col[1],col[3]);
        compareFilesForChannel(files1[1], files2[1], channel, outputDir, sourceDir,col[0],col[2]);
    }

    return 0;
}


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

void compareFilesForChannel(const char* file1, const char* file2,const char* file3, int channel, const char* outputDir, const char* sourceDir, Color_t color1, Color_t color2, Color_t color3){
    TString filename1 = TString::Format("%s%s.root", sourceDir,file1);
    TString filename2 = TString::Format("%s%s.root",sourceDir, file2);
    TString filename3 = TString::Format("%s%s.root",sourceDir,file3);

    TFile* f1 = TFile::Open(filename1);
    TFile* f2 = TFile::Open(filename2);
    TFile* f3 = TFile::Open(filename3);

    if (!f1 || f1->IsZombie() || !f2 || f2->IsZombie() || !f3 || f3->IsZombie()) {
        std::cerr << "Error opening files: " << filename1 << " or " << filename2 << "or" << filename3 << std::endl;
        return;
    }

    TString profName1 = TString::Format("AveragePulse_Ch%d", channel);
    TString profName2 = TString::Format("AveragePulse_Ch%d", channel);
    TString profName3 = TString::Format("AveragePulse_Ch%d",channel);

    TProfile* prof1 = dynamic_cast<TProfile*>(f1->Get(profName1));
    TProfile* prof2 = dynamic_cast<TProfile*>(f2->Get(profName2));
    TProfile* prof3 = dynamic_cast<TProfile*>(f3->Get(profName3));

    if (!prof1 || !prof2 || !prof3) {
        std::cerr << "Error retrieving profiles: " << profName1 << " or " << profName2 <<"or"<< profName3 << std::endl;
        f1->Close();
        f2->Close();
	f3->Close();
        return;
    }

    TGraph* graph1 = profileToGraph(prof1);
    TGraph* graph2 = profileToGraph(prof2);
    TGraph* graph3 = profileToGraph(prof3);
    if (!graph1 || !graph2 || !graph3) {
        std::cerr << "Error converting profiles to graphs: " << profName1 << " or " << profName2 <<"or"<< profName3 << std::endl;
        f1->Close();
        f2->Close();
	f3->Close();
        delete graph1;
        delete graph2;
	delete graph3;
        return;
    }

    double max1 = getMaxValue(graph1);
    double max2 = getMaxValue(graph2);
    double max3 = getMaxValue(graph3);

    TCanvas* canvas = new TCanvas(TString::Format("AvgPulseFilterComparisons_Channel%d", channel),
                                  TString::Format("Average Pulse Filter Comparisons Channel %d", channel), 800, 600);

    
        graph1->SetLineColor(color1);  
        graph1->SetLineStyle(1);      
	graph1->GetXaxis()->SetTitle("t-T_{trig} (ns)");
        graph1->GetYaxis()->SetTitle("Amp (mV)");
        graph1->SetTitle(Form("Average Pulse 610 nm Vs. 560 nm:  Ch_%d",channel));
        graph1->SetStats(0);
	graph1->GetYaxis()->SetRangeUser(-1.5,0.5);
	//graph1->GetYaxis()->SetRangeUser(-max1-10,10);
        graph1->Scale(1/max1);
	graph1->Draw("AL");
        graph2->Scale(1/max2);
        graph2->SetLineColor(color2);   
        graph2->SetLineStyle(2);       
        graph2->Draw("L SAME");
	graph3->Scale(1/max3);
        graph3->SetLineColor(color3);
        graph3->SetLineStyle(2);
        graph3->Draw("L SAME");
//    } else if (max2 > max1 && max2 > max3) {
//       graph2->SetLineColor(color2);   
//        graph2->SetLineStyle(2);       
//        graph2->GetXaxis()->SetTitle("t-T_{trig} (ns)");
//        graph2->GetYaxis()->SetTitle("Amp (mV)");
//        graph2->SetTitle(Form("Average Pulse 610 nm Filter Vs. 560 nm: Ch_%d",channel));
//        graph2->SetStats(0);
//	graph2->GetYaxis()->SetRangeUser(-1.5,0.5);
//	graph2->Scale(1/max2);
//	graph2->Draw("AL");
//        graph1->Scale(1/max1);
//        graph1->SetLineColor(color1);  
//        graph1->SetLineStyle(1);      
//        graph1->Draw("L SAME");
//	graph3->Scale(1/max3);
//        graph3->SetLineColor(color3);
//        graph3->SetLineStyle(2);
//        graph3->Draw("L SAME");
//    }else{
//	graph3->SetLineColor(color3);
//        graph3->SetLineStyle(2);
//        graph3->GetXaxis()->SetTitle("t-T_{trig} (ns)");
//        graph3->GetYaxis()->SetTitle("Amp (mV)");
//        graph3->SetTitle(Form("Average Pulse 610 nm Filter Vs. 560 nm: Ch_%d",channel));
//        graph3->SetStats(0);
//        graph3->GetYaxis()->SetRangeUser(-1.5,0.5);
//        graph3->Scale(1/max3);
//        graph3->Draw("AL");
//        graph1->Scale(1/max1);
//        graph1->SetLineColor(color1);
//        graph1->SetLineStyle(1);
//        graph1->Draw("L SAME");
//        graph2->Scale(1/max2);
//        graph2->SetLineColor(color2);
//        graph2->SetLineStyle(2);
//        graph2->Draw("L SAME"); 
//	 }

    TString canvasName = TString::Format("%s/Norm_AvgPulseFilterComp_Channel%d.png", outputDir, channel);
    canvas->SaveAs(canvasName);

    std::cout << "Saved canvas: " << canvasName << std::endl;

    delete graph1;
    delete graph2;
    delete graph3;
    delete canvas;
 
    
    f1->Close();
    f2->Close();
    f3->Close();
}

int main() {
    const char* files1[] = {"run594_LG"};
    const char* files2[] = {"run280_LG"};
    const char* files3[] = {"run526_LG"};
    const char* outputDir = "testbeam_plots/June/week3/plots";
    const char* sourceDir = "testbeam_plots/June/week3/";
    Color_t col[6] ={kGreen,kRed,kMagenta,kBlue, kOrange, kPink}; 
    ensureDirectoryExists(outputDir);

    for (int channel = 0; channel < 8; ++channel) {
        compareFilesForChannel(files1[0], files2[0],files3[0], channel, outputDir, sourceDir,col[3],col[4],col[5]);
        //compareFilesForChannel(files1[1], files2[1], channel, outputDir, sourceDir,col[2],col[3]);
    }

    return 0;
}


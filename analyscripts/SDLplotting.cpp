#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <TAxis.h>
#include <TMath.h>

double GetXpoint(TGraph* graph, double yValue){
int NumP = graph -> GetN();
double* xValues = graph->GetX();
double*  yValues = graph->GetY();

for (int i = 0; i < NumP-1 ; i++){
	 if ((yValues[i] <= yValue && yValues[i+1] >= yValue) || (yValues[i] >= yValue && yValues[i+1] <= yValue)) {
             double x1 = xValues[i];
             double x2 = xValues[i+1];
             double y1 = yValues[i];
             double y2 = yValues[i+1];
             //(x - x1) / (x2 - x1) = (yValue - y1) / (y2 - y1)
             double x = x1 + (yValue - y1) * (x2 - x1) / (y2 - y1);
             return x;


}

}
std::cerr << "Warning: Specified yValue is out of the range of the graph." << std::endl;
    return -1;
}


int main() {
    TFile* file = TFile::Open("../DESYSDLTGraphs_Regions.root", "READ");
    if (!file || !file->IsOpen()) {
        std::cerr << "Error: Could not open SDLTGraphs.root file!" << std::endl;
    }
    TGraph* residual[8][10];//Made an array of these for future studies of delays
    // Loop over channels and delays
    for (int ich = 0; ich < 8; ich++) {
        for (int del = 0; del <= 0; del++) {//Just Focused on 1 ns delay
      	    residual[ich][del] = new TGraph();//Form("Residual_Ch%d_del%dns",ich,del),Form("Residual_Ch%d_del%dns",ich,del));//Probably needs more arguments      
	    TString nameU = Form("R2SDL_Ch%d_del%d_Unfil", ich, del);
            TString nameF = Form("R2SDL_Ch%d_del%d_RG610", ich, del);
	    //TGraph* graphR = new TGraph();
            TGraph* graphU = (TGraph*)file->Get(nameU);
            TGraph* graphF = (TGraph*)file->Get(nameF);
            TGraph* agraphU = new TGraph();  	
	    double* xU = graphU->GetX();
	    double* yU = graphU->GetY();
	    double* yF = graphF->GetY();
	    double alignedxU[graphU->GetN()];
            if (!graphU || !graphF) {
                std::cerr << "Error: Could not find graphs " << nameU << " or " << nameF << " in file!" << std::endl;
                continue;
            }

            double maxU = graphU->GetY()[TMath::LocMax(graphU->GetN(), graphU->GetY())];
            double maxF = graphF->GetY()[TMath::LocMax(graphF->GetN(), graphF->GetY())];
	    double AtU = GetXpoint(graphU,maxU/2);
	    double AtF = GetXpoint(graphF,maxF/2);
	    double timediff = AtU-AtF;
           for (int i = 0; i < graphU->GetN()-1; i++){
		if(timediff < 0){
		alignedxU[i] = xU[i]-timediff;
		}else{
		alignedxU[i]=xU[i]+timediff;
		}
		agraphU->SetPoint(i,alignedxU[i], yU[i]);	
                residual[ich][del]->SetPoint(i,alignedxU[i],yU[i]-yF[i]);

            }	
	

            TCanvas* canvas = new TCanvas("canvas", "Canvas", 800, 600);
            canvas->cd();
	    residual[ich][del]->SetLineColor(kGreen);
	    residual[ich][del]->SetLineWidth(2);
            agraphU->SetLineColor(kBlue);
            graphU->SetLineWidth(2);
            graphF->SetLineColor(kRed);
            graphF->SetLineWidth(2);

            graphU->GetXaxis()->SetTitle("t - t_{trig} (ns)");
            graphU->GetYaxis()->SetTitle("a(t) - a(t-d)");
	    graphF->GetXaxis()->SetTitle("t - t_{MCP} (ns)");
            graphF->GetYaxis()->SetTitle("a(t) - a(t-d)");
	    agraphU->SetTitle(Form("Channel %d Delay %d ns",ich,del+1));
            graphF->SetTitle(Form("Channel %d Delay %d ns",ich,del+1));
	    agraphU->GetXaxis()->SetRangeUser(0,500);
	    graphF->GetXaxis()->SetRangeUser(0,500);	
            // Draw the graphs with the larger maximum first
            if (maxU > maxF) {
                agraphU->Draw("AL");
                graphF->Draw("L SAME");
		residual[ich][del]->Draw("L SAME");
            } else {
                graphF->Draw("AL");
                agraphU->Draw("L SAME");
		residual[ich][del]->Draw("L SAME");
            }

            TLegend* legend = new TLegend(0.75, 0.75, 0.9, 0.9);
            legend->AddEntry(agraphU, "Unfiltered", "l");
            legend->AddEntry(graphF, "Filtered", "l");
	    legend->AddEntry(residual[ich][del],"Residual","l");
            legend->Draw();

            canvas->Print(Form("../testbeam_plots/July/week4/SDLFVU_Align/SDL_Ch%d_del%d.png", ich, del));

            delete legend;
            delete canvas;
        }
    }

    file->Close();
    delete file;
}


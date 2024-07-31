#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <TAxis.h>
#include <TMath.h>
double alignmentsetup(TGraph* graph, int ich){
int NumP = graph->GetN();
double* XValues = graph ->GetX();
double* YValues = graph->GetY();
double minY = std::numeric_limits<double>::max();
double halfmax = 0.0;
double Xlim = 0.0;
int turnpoint = 0;
double x = 100000;
//double xstep = 0;
	for (int i = 0; i < NumP -1; i++){
		if (ich == 0){
//		std::cout<<"y values are: "<<YValues[i]<<std::endl;
	}	if (minY > YValues[i]){
			minY = YValues[i];
			Xlim = XValues[i];
			turnpoint = i;
		}
		halfmax = minY/2;
	}
		
	std::cout<<"The Min for channel ("<<ich<<") is: "<<minY<< " mV...and that happens at: "<<Xlim<<" ns"<<std::endl;
		

	for (int j = 0; j < turnpoint ; j++){
		
		if ((YValues[j] <= halfmax && YValues[j+1] >= halfmax) || (YValues[j] >= halfmax && YValues[j+1] <= halfmax)) {
             		double x1 = XValues[j];
             		double x2 = XValues[j+1];
             		double y1 = YValues[j];
             		double y2 = YValues[j+1];
            		//(x - x1) / (x2 - x1) = (halfmax - y1) / (y2 - y1)
                        double x = x1 + (halfmax - y1) * (x2 - x1) / (y2 - y1);
                        return x;
             
		}		
	}

 return x;
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
	    double* xF = graphF->GetX();
	    double* yF = graphF->GetY();
	    int overlap = 0;
	    double alignedxU[graphU->GetN()];
            if (!graphU || !graphF) {
                std::cerr << "Error: Could not find graphs " << nameU << " or " << nameF << " in file!" << std::endl;
                continue;
            }

	    double AtU = alignmentsetup(graphU,ich);
	    double AtF = alignmentsetup(graphF,ich);
	    double timediff = AtU-AtF;
	   std::cout<<"The time diff for channel ("<<ich<<") is: "<<timediff<<std::endl;
           for (int i = 0; i < graphU->GetN()-1; i++){
		if(timediff < 0){
		alignedxU[i] = xU[i]-timediff;
		}else{
		alignedxU[i]=xU[i]+timediff;
		}
		if(alignedxU[i] == xF[i]){
		overlap++;
		}
		agraphU->SetPoint(i,alignedxU[i], yU[i]);	
                residual[ich][del]->SetPoint(i,alignedxU[i],yU[i]-yF[i]);

            }	
	std::cout<<"overlap: "<<overlap<<std::endl;

            TCanvas* canvas = new TCanvas("canvas", "Canvas", 800, 600);
            canvas->cd();
	    residual[ich][del]->SetLineColor(kGreen);
	    residual[ich][del]->SetLineWidth(2);
            agraphU->SetLineColor(kBlue);
           agraphU->SetLineWidth(2);
            graphF->SetLineColor(kRed);
            graphF->SetLineWidth(2);

            agraphU->GetXaxis()->SetTitle("t - t_{trig} (ns)");
            agraphU->GetYaxis()->SetTitle("a(t) - a(t-d)");
	    graphF->GetXaxis()->SetTitle("t - t_{MCP} (ns)");
            graphF->GetYaxis()->SetTitle("a(t) - a(t-d)");
	    agraphU->SetTitle(Form("Channel %d Delay %d ns",ich,del+1));
            graphF->SetTitle(Form("Channel %d Delay %d ns",ich,del+1));
	    agraphU->GetXaxis()->SetRangeUser(0,500);
	    graphF->GetXaxis()->SetRangeUser(0,500);	
            // Draw the graphs with the larger maximum first
          //  if (maxU > maxF) {
           //     agraphU->Draw("AL");
            //    graphF->Draw("L SAME");
//		residual[ich][del]->Draw("L SAME");
  //          } else {
                graphF->Draw("AL");
                agraphU->Draw("L SAME");
		residual[ich][del]->Draw("L SAME");
    //        }

            TLegend* legend = new TLegend(0.75, 0.75, 0.9, 0.9);
            legend->AddEntry(agraphU, "Unfiltered", "l");
            legend->AddEntry(graphF, "Filtered", "l");
	    legend->AddEntry(residual[ich][del],"Residual","l");
            legend->Draw();
		//if(ich < 4){
            canvas->SaveAs(Form("../SDLFVU_Align/SDL_Ch%d_del%d.png", ich, del));
//	}
          delete legend;
            delete canvas;
        }
    }

    file->Close();
    delete file;
}


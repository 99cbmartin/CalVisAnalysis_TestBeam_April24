#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <TAxis.h>
#include <TMath.h>
#include <iomanip> 
#include <cmath>
double alignmentsetup(TGraph* graph, int ich){
int NumP = graph->GetN();
double* XValues = graph ->GetX();
double* YValues = graph->GetY();
double minY = std::numeric_limits<double>::max();
double halfmax = 0.0;
double Xlim = 0.0;
int turnpoint = 0;
double x = 100000;
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
		
//	std::cout<<"The Min for channel ("<<ich<<") is: "<<halfmax<< " mV...and that happens at: "<<Xlim<<" ns"<<std::endl;
		

	for (int j = 0; j < turnpoint ; j++){
		
		if ((YValues[j] <= halfmax && YValues[j+1] >= halfmax) || (YValues[j] >= halfmax && YValues[j+1] <= halfmax)) {
             		double x1 = XValues[j];
             		double x2 = XValues[j+1];
             		double y1 = YValues[j];
             		double y2 = YValues[j+1];
            		//(x - x1) / (x2 - x1) = (halfmax - y1) / (y2 - y1)
                        double x = x1 + (halfmax - y1) * (x2 - x1) / (y2 - y1);
 std::cout<<"The halfMin for channel ("<<ich<<") is: "<<halfmax<< " mV...and that happens at: "<<x<<std::endl;        
              return x;
             
		}		
	}
	//std::cout<<"The halfMin for channel ("<<ich<<") is: "<<halfmax<< " mV...and that happens at: "<<x<<std::endl;
 return x;
//std::cout<<"The halfMin for channel ("<<ich<<") is: "<<halfmax<< " mV...and that happens at: "<<x<<std::endl;
}

TGraph* residualplot(int ich, TGraph* Ugraph, TGraph* Fgraph, double* alignedXU){//must feed it the two original TGraphs and the newly aligned XArray

TGraph* residualTG = new TGraph();
int NumP = Ugraph->GetN();
//double* newX = alignedXU();//quesationable
double NewUY[NumP];// = {0.0};
double NewFY[NumP];// = {0.0};
double* oldUY = Ugraph->GetY();
double* oldFY= Fgraph->GetY();
double* oldUX = Ugraph->GetX();
double* oldFX = Fgraph->GetX();
double* residualY[NumP];


// Print table header
    std::cout << std::left << std::setw(10) << "Channel " 
              << std::setw(10) << "TargetX " 
              << std::setw(10) << "oldX 1 " 
              << std::setw(10) << "oldX 2 " 
              << std::setw(10) << "oldY 1 " 
              << std::setw(10) << "oldY 2 " 
              << std::setw(15) << "New Y Value " 
              <<std::setw(15)<<"Old F Y "
	      <<std::setw(15)<<"residual "
	      << std::endl;
    std::cout << std::string(75, '-') << std::endl;

 for (int k  = 0 ; k < NumP-1 ; k++){

//std::cout<<"in Channel "<<ich<<" the "<<k<<"th point value was: "<<oldUX[k]<<" but now is: "<<alignedXU[k]<<std::endl;
//std::cout<<"in Channel "<<ich<<" the "<<k<<"th point value was: "<<oldFX[k]<<std::endl;
double targetX = alignedXU[k];
double Ux1 = 10000;
double Ux2 = 10000;
double Uy1 = 10000;
double Uy2 = 10000;
	for (int j = 49; j < 150 ; j++){
                	if ((oldUX[j] <= targetX && oldUX[j+1] >= targetX) || (oldUX[j] >= targetX && oldUX[j+1] <= targetX)) {
                        	Ux1 = oldUX[j];
                      		Ux2 = oldUX[j+1];
                        	Uy1 = oldUY[j];
                        	Uy2 = oldUY[j+1];
                        	//(targetX - x1) / (x2 - x1) = (newY[k] - y1) / (y2 - y1)
                        	NewUY[k] = Uy1 + (targetX - Ux1)*(Uy2-Uy1)/(Ux2-Ux1); //Check this later
            			



				std::cout << std::left << std::setw(10) << ich<< " "
                          << std::setw(10) << targetX<< " "
                          << std::setw(10) << Ux1<< " "
                          << std::setw(10) << Ux2<< " "
                          << std::setw(10) << Uy1<< " "
                          << std::setw(10) << Uy2<< " "
                          << std::setw(15) << NewUY[k]<< " "
			  << std::setw(15)<<oldFY[k]<<" "
			  <<std::setw(15)<<NewUY[k]-oldFY[k]<<" "
                          << std::endl;            	
			}
	}
			if (k < 150 && k > 99){

//			std::cout<<"the channel: "<<ich<<std::endl;
	//		std::cout<<"The target is: "<<targetX<<std::endl;
	//		std::cout<<"firstly: x1= "<<Ux1<<" then x2="<<Ux2<<" their diff= "<<Ux2-Ux1<<std::endl;
	//		std::cout<<"Then: y1= "<<Uy1<<" then y2="<<Uy2<<" their diff= "<<Uy2-Uy1<<std::endl;
	//		std::cout<<"Finding the interpolated value of: "<<NewUY[k]<<std::endl;
	//		std::cout<<"The Residual: "<<NewUY[k]-oldFY[k]<<std::endl;


                        }
                       


	}
	for (int l  = 0 ; l < 140 ; l++){
		residualTG->SetPoint(l,oldUX[l],NewUY[l]-oldFY[l]);
		 //std::cout<<NewUY[l]<<" minus "<<oldFY[l]<<" gives the "<< l<<"th residual: "<<NewUY[l]-oldFY[l]<<std::endl;
		//std::cout<<NewUY[l]-oldFY[l]<<" is the value of the "<< l<<"th residual"<<std::endl;

	}	
return residualTG;
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
	    TString nameU = Form("SDL_Ch%d_del%d_Unfil", ich, del);
            TString nameF = Form("SDL_Ch%d_del%d_RG610", ich, del);
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
	   //std::cout<<"The time diff for channel ("<<ich<<") is: "<<timediff<<std::endl;
           for (int i = 0; i < graphU->GetN()-1; i++){
		if(timediff < 0){
		alignedxU[i] = xU[i]+timediff;
		}if (timediff > 0){
		alignedxU[i]=xU[i]-timediff;
		}
		if(alignedxU[i] == xF[i]){
		overlap++;
		}
		std::cout << "Original xU[" << i << "]: " << xU[i] << ", Shifted xU[" << i << "] by: "<<timediff<< " to get: " << alignedxU[i] << std::endl;

		agraphU->SetPoint(i,xU[i]-timediff, yU[i]);	
               // residual[ich][del]->SetPoint(i,alignedxU[i],yU[i]-yF[i]);

            }	
	//std::cout<<"overlap: "<<overlap<<std::endl;
	    residual[ich][del] = residualplot(ich,graphU,graphF,alignedxU);//Check here
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
  
		graphF->GetYaxis()->SetRangeUser(-1.2,0.2);
		graphF->GetXaxis()->SetRangeUser(80,140);
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
            canvas->SaveAs(Form("../SDLFVU_Align_Z/FullSDL_Ch%d_del%d.png", ich, del),"Q");
//	}
          delete legend;
            delete canvas;
        }
    }

    file->Close();
    delete file;
}


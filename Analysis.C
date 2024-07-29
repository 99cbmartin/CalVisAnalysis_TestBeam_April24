#define Analysis_cxx
#include "Analysis.h"
#include <TH1D.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TProfile.h>
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include "TMath.h"
#include <TFile.h>
#include <TString.h>
#include <filesystem>

//quick function for getting the max of a TGraph
double getMaxValue(TGraph* graph) {
    double maxVal = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < graph->GetN()-24; ++i) {//snip the last 24 ns off since there is large noise jump
        double x, y;
        graph->GetPoint(i, x, y);
        if (fabs(y) > maxVal) {
            maxVal =fabs(y);
        }
    }
    return maxVal;
}

TGraph* Convert(TProfile* prof) {//Quick function to turn a TProfile to a TGraph
    if (!prof) return nullptr;
    TGraph* graph = new TGraph();
    for (int i = 1; i <= prof->GetNbinsX(); ++i) {
        graph->SetPoint(i-1, prof->GetBinCenter(i), prof->GetBinContent(i));
    }
    return graph;
}

void Analysis::Loop(std::string outputFileName)
{
Int_t icol[8] = { kOrange+10, kOrange-3, kSpring-1, kAzure+6, kAzure-1, kViolet-1, kBlack, kRed+1};//first and last color too close
TH1D* triggerRMS = new TH1D("triggerRMS","Trigger RMS;Value;Counts",100,-1,1);
TH1D* trigger_Amp = new TH1D("trigger_Amp", "Trigger Amplitude; Value;Counts",100,400,500);
TH1D* trigger_time = new TH1D("trigger_time","Trigger Time;time(ns);Counts",100,0,10);
TH1D* h001 = new TH1D("EventRMS","Baseline Noise of Events;RMS;Counts",100,0,100);
//TCanvas* triggerAlignment = new TCanvas("triggerAlignment","Aligned Triggers",800,600);
//TCanvas* SingleEvents = new TCanvas("SingleEvents", "Aligned Single Events",800,600);
TProfile* AveragePulse[8];
TProfile* RAveragePulse[8][4];//Organized by Channel and then Regions 1-4
//TProfile* R2AveragePulse[8];
//TProfile* R3AveragePulse[8];
//TProfile* R4AveragePulse[8];
TGraph* AveragePulseGraph[8];
TGraph* AvgPGraphR[8][4];//Also Organized Regions 1-4
//TGraph* AvgPGraphR2[8];
//TGraph* AvgPGraphR3[8];
//TGraph* AvgPGraphR4[8];
TProfile* AveragePulseFront = new TProfile("AveragePulse_Front","AveragePulse_Front",1024,0,1023);
TProfile* AveragePulseBack = new TProfile("AveragePulse_Back","AveragePulse_Back",1024,0,1023);
TProfile* AveragePulseMid = new TProfile("AveragePulse_Mid","AveragePulse_Mid",1024,0,1023);//Added this to cover any leftouver distinct pulses such as Channels 5 and 6 in run526
TH1D* h_noiseRMS[8];
TH1D* h_noiseRMStotal = new TH1D("AllEvtRMS","RMS of All Events",600,0,8);
TH1D* h_eventAmp = new TH1D("EvtAmps","Amplitudes of Events",1000,0,1000);
TH1D* h_ChAmp[8];
TH1D* baselinehisto = new TH1D("AllBaseline","AllBaseline",100,0,100);
TH1D* baselineH[8];
TProfile* AvgSDL[8][10];
TProfile* RAvgSDL[8][10][4];
//TProfile* R2AvgSDL[8][10];
//TProfile* R3AvgSDL[8][10];
//TProfile* R4AvgSDL[8][10];
int Delays[10]={1,2,3,4,5,6,7,8,9,10};//Delay times in ns
double RAvgmax[8][4] ={0,0};
for (int ich =0 ; ich <8; ich++){
	for(int del=0; del < 10;del++){
		AvgSDL[ich][del] = new TProfile(Form("Average_SDL_Ch%d_del%d_ns",ich,Delays[del]),Form("Average_SDL_Ch%d_del%d_ns",ich,Delays[del]),1024,0,1023);
	 	for (int R = 1; R < 5; R++){ 
			RAvgSDL[ich][del][R] = new TProfile(Form("R%dAverage_SDL_Ch%d_del%d_ns",R,ich,Delays[del]),Form("R%dAverage_SDL_Ch%d_del%d_ns",R,ich,Delays[del]),1024,0,1023);
		}
	 //R2AvgSDL[ich][del] = new TProfile(Form("R2Average_SDL_Ch%d_del%d_ns",ich,Delays[del]),Form("R2Average_SDL_Ch%d_del%d_ns",ich,Delays[del]),1024,0,1023);
	// R3AvgSDL[ich][del] = new TProfile(Form("R3Average_SDL_Ch%d_del%d_ns",ich,Delays[del]),Form("R3Average_SDL_Ch%d_del%d_ns",ich,Delays[del]),1024,0,1023);
	// R4AvgSDL[ich][del] = new TProfile(Form("R4Average_SDL_Ch%d_del%d_ns",ich,Delays[del]),Form("R4Average_SDL_Ch%d_del%d_ns",ich,Delays[del]),1024,0,1023);
	 }
 }

//Time Stuff
Int_t samples = 1024; //Number of samples performed by the DRS system
double time[1024];
double lowertime = 0.0;
double uppertime = 0.0;
double average = 0.0;
int num = 0;
int num2 = 0;

for (int ich = 0; ich <8; ich++){
	AveragePulse[ich] = new TProfile(Form("AveragePulse_Ch%d",ich),Form("Average Pulse for Channel %d",ich),samples,0,1023);//will need to impliment this to be more modular for upper and lower bins in cases of hori_interv \=1
	for (int R = 1 ; R < 5; R++){
		RAveragePulse[ich][R] = new TProfile(Form("R%dAveragePulse_Ch%d",R,ich),Form("R%d Average Pulse for Channel %d",R,ich),samples,0,1023);
		AvgPGraphR[ich][R] = new TGraph();
	}
//R2AveragePulse[ich] = new TProfile(Form("R2AveragePulse_Ch%d",ich),Form("R2 Average Pulse for Channel %d",ich),samples,0,1023);
//R3AveragePulse[ich] = new TProfile(Form("R3AveragePulse_Ch%d",ich),Form("R3 Average Pulse for Channel %d",ich),samples,0,1023);
//R4AveragePulse[ich] = new TProfile(Form("R4AveragePulse_Ch%d",ich),Form("R4 Average Pulse for Channel %d",ich),samples,0,1023);
AveragePulseGraph[ich] = new TGraph();//Form("AveragePulseG_Ch%d",ich),Form("Average Pulse for Channel %d",ich));
//AvgPGraphR1[ich] = new TGraph();
//AvgPGraphR2[ich] = new TGraph();
//AvgPGraphR3[ich] = new TGraph();
//AvgPGraphR4[ich] = new TGraph();
auto chstr = std::to_string(ich);
h_noiseRMS[ich] = new TH1D(Form("RMS_Noise_Ch%d",ich),Form("RMS_Noise_Ch%d",ich),200,0,3);
h_ChAmp[ich] = new TH1D(Form("ChAmp%d",ich),Form("ChAmp%d",ich),300,0,1000);
baselineH[ich] = new TH1D(Form("BaselineCh%d",ich),Form("BaselineCh%d",ich),100,0,100);
}

  TFile* pulseFile = new TFile(outputFileName.c_str(),"recreate");
  TGraph * grSnglPul[8];
  TGraph * grSnglPulA[8];
  TGraph* grTrigPul;
  TGraph* grTrigPulA;
  if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {//begin event loop
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
bool isOK = true;//Quality event flag to grow with requirments to be determined [Good spot?]

//===================================\\
//              Trigger              \\
//===================================\\
//Trigger Analysis and timing alignment

//Take the average of the trigger before the step (A)
//Average of the trigger after the step (B)
//find the amplitude of the trigger time (A+B/2)
//finally take the timestamp of that achieved trigger amplitude
double trigrange[4] = {1,100,300,1000};//Some ranges for the base and plateu of the trigger
double s[4];
double trigbase=0.0;
double trigplat=0.0; 
double trigamp=0.0;//(A+B)/2, the amplitude at which I count the trigger time 
double trigger_t=0.0; 
double noise=0.0; //RMS of the same area before the step of the pulse
double trigger_amp=0.0;//A-B, the amplitude of the actual trigger pulse with respect to the baseline
//the key is the horizontal_interval
//define a time array for each event based on horizontal interval
for (int sam = 0 ;sam < samples; sam++){
time[sam] = sam * horizontal_interval;

}	
for (int n =0 ; n<4 ; n++){
	s[n] =trigrange[n]/horizontal_interval;
}

double trigbasesum=0.0;
double trigbasecount=0;
double trigplatsum=0.0;
double trigplatcount=0;
int base = s[1]-s[0];
std::vector<double> basevalues;

for (int isamp = 0; isamp < samples; isamp++){

	if (isamp > s[0] && isamp < s[1]){
		trigbasecount++;
		trigbasesum+= trigger[0][isamp];
		basevalues.push_back(trigger[0][isamp]);
	}
	if (isamp > s[2] && isamp < s[3]){
		trigplatcount++;
		trigplatsum+=trigger[0][isamp];
	}		

}

trigbase = trigbasesum/trigbasecount;
trigplat = trigplatsum/trigplatcount;
trigamp = (trigbase+trigplat)/2;
trigger_amp = fabs(trigbase -trigplat);
noise = TMath::RMS(basevalues.size(),basevalues.data());

for (int isamp = 0; isamp < samples; isamp++){
	if((trigger[0][isamp] < trigamp && trigger[0][isamp+1] >= trigamp) ||
	   (trigger[0][isamp] > trigamp && trigger[0][isamp+1] <= trigamp)){
	double t1 = isamp * horizontal_interval;
	double t2 = (isamp + 1) *horizontal_interval;
	double v1 = trigger[0][isamp];
	double v2 = trigger[0][isamp +1];
	trigger_t = t1 + (trigamp -v1)*(t2-t1)/(v2-v1);
	break;
	}
}
lowertime =(time[0]-trigger_t);
uppertime = (time[1023]);
triggerRMS->Fill(noise);
trigger_Amp->Fill(trigger_amp);
trigger_time->Fill(trigger_t);

auto evtnstr = std::to_string(event);

//===================================\\
//           Pre-Analysis	     \\
//===================================\\
 //In this section, I derive some basic results and accumulate quality event management 
 

 double channelMax[8] = {0.0};
 double noiseRMS[8];
 double noisebin = 100/horizontal_interval;//finds the final bin for 100 ns point which beyond, the pulse noise is expected to dominate less    
 double baseline[8] = {0.0};
 double amp[8][1024] = {0,0};//Just the adjusted output of SiPM's with consideration of the baseline such that pulses align with y = 0 mV
grTrigPul = new TGraph();
string trigname = "trigger_event"+evtnstr;
grTrigPul->SetName(trigname.c_str());
grTrigPulA = new TGraph();
string trignameA = "Alitrigger_event"+evtnstr;
grTrigPulA->SetName(trignameA.c_str());

for (int ich=0;ich<8; ich++){
	double baselinesum=0.0;
	double baselinecount=0.0;
	noiseRMS[ich] = TMath::RMS(noisebin,channels[ich]);	
	h_noiseRMS[ich]->Fill(noiseRMS[ich]);
	h_noiseRMStotal->Fill(noiseRMS[ich]);	
	
	grSnglPul[ich] = new TGraph();
	grSnglPulA[ich] = new TGraph();
	auto chstr = std::to_string(ich);
	string gname = "singlePulse_event"+evtnstr+"_ch"+chstr;
	grSnglPul[ich]->SetName(gname.c_str());
	string gnameA = "AlignedSinglePulses_event"+evtnstr+"_ch"+chstr;
	grSnglPulA[ich]->SetName(gnameA.c_str());
	
	
	for (int isamp=0;isamp<samples;isamp++){
		if ( isamp > 1 && isamp < 35){
			
			baselinesum += channels[ich][isamp];
			baselinecount++;
		}	
	 
		grSnglPul[ich]->SetPoint(isamp,(isamp*horizontal_interval),channels[ich][isamp]);
	  	
	}
baseline[ich] = baselinesum/baselinecount;
 
	for (int isamp=0;isamp<samples;isamp++){
	   amp[ich][isamp] = channels[ich][isamp]-baseline[ich]; 
		if (fabs(channels[ich][isamp]-baseline[ich]) > fabs(channelMax[ich])){

                        channelMax[ich] = fabs(channels[ich][isamp]-baseline[ich]);
                }
        }
//isOK = isOK && (channelMax[ich] > 3* noiseRMS[ich]);//Event Amplitude quality check, possible to use 4 or 5
isOK = isOK && noise < 1 && noise > 0.5 && trigger_amp > 440 && trigger_amp < 455;//Trigger conditions for quality events
//isOK = isOK && channelMax[ich] < 10 && noiseRMS[ich] <1.1;//SiPM conditions for quality events, Needs to be enforced channel-by-channel due to overlap of dead-channel activity in quality events
baselinehisto->Fill(baseline[ich]);
baselineH[ich]->Fill(baseline[ich]);
}//channel loop end

//===================================\\
//           Pulse Analysis          \\
//===================================\\
//Actual Analysis of Pulse Shapes and Timing

if (isOK){//Quality Events only beyond this point
	for (int ich=0;ich<8; ich++){
	 		

		
		if (channelMax[ich] > 6*noiseRMS[ich] && noiseRMS[ich] < 1){//Manual Event Quality Cuts
			AveragePulse[ich]->SetLineColor(icol[ich]);
			grSnglPulA[ich]->SetLineColor(icol[ich]);
			grSnglPulA[ich]->GetXaxis()->SetTitle("time - T_{Trig} (ns)");
			grSnglPulA[ich]->GetYaxis()->SetTitle("a (mV)");

			Float_t ampchan = *std::max_element(channels[ich],channels[ich]+1024);
			//std::cout<<"  Channel "<<ich<<" max "<<ampchan<<std::endl;

			h_eventAmp->Fill(channelMax[ich]);
			h_ChAmp[ich]->Fill(channelMax[ich]);
			double SDLRegions[8][4]=
			{ {0,75,300,600},{0,50,150,300},{0,100,300,600},{0,50,200,400},//Front-End Channels
			  {0,15,35,60},{0,0,0,0},{0,0,0,0},{0,15,60,120} //Back-End Channels (5 & 6 Often Dead)
	                };

				for (int isamp=0;isamp<samples;isamp++){
          			//	AveragePulse[ich]->Fill(isamp*horizontal_interval-trigger_t+abs(lowertime),channels[ich][isamp]-baseline[ich]);
          				AveragePulse[ich]->Fill(isamp*horizontal_interval-trigger_t+abs(lowertime),amp[ich][isamp]);
        				grSnglPulA[ich]->SetPoint(isamp,(isamp*horizontal_interval-trigger_t),(channels[ich][isamp]-baseline[ich]));
					if (ich < 4){
					AveragePulseFront->Fill(isamp*horizontal_interval-trigger_t+abs(lowertime),amp[ich][isamp]);
					}if(ich == 4 || ich == 7){//adjust this to include channels 5 & 6 *if* those channels are not dead
					AveragePulseBack->Fill(isamp*horizontal_interval-trigger_t+abs(lowertime),amp[ich][isamp]);
					}
					if (ich == 5 || ich == 6){
					AveragePulseMid->Fill(isamp*horizontal_interval-trigger_t+abs(lowertime),amp[ich][isamp]);
					}
					
				//Regions 1 through 3 cut on amplitude for SDL studies
		
				for (int R = 1; R <3 ; R++){//Maintaining the Regions 1-3 at the moment, to augment must expand the arguments of SDLRegion array
					
					if ( channelMax[ich] > SDLRegions[ich][R-1] && channelMax[ich] < SDLRegions[ich][R]){
					RAveragePulse[ich][R]->Fill(isamp*horizontal_interval-trigger_t_abs(lowertime),amp[ich][isamp]);
					}	


					/*
  					if (ich == 0 && channelMax[0] > 0 && channelMax[0] < 75){
					R1AveragePulse[ich]->Fill(isamp*horizontal_interval-trigger_t+abs(lowertime),amp[ich][isamp]);	
					}
					if (ich == 0 && channelMax[0] > 75 && channelMax[0] < 300){
					R2AveragePulse[ich]->Fill(isamp*horizontal_interval-trigger_t+abs(lowertime),amp[ich][isamp]);
                                        }
					if (ich == 0 && channelMax[0] > 300 && channelMax[0] < 600){
					R3AveragePulse[ich]->Fill(isamp*horizontal_interval-trigger_t+abs(lowertime),amp[ich][isamp]);	
                                        }
					if (ich == 1 && channelMax[1] > 0 && channelMax[1] < 50){
					R1AveragePulse[ich]->Fill(isamp*horizontal_interval-trigger_t+abs(lowertime),amp[ich][isamp]);
					}
					if (ich == 1 && channelMax[1] > 50 && channelMax[1] < 150){
					R2AveragePulse[ich]->Fill(isamp*horizontal_interval-trigger_t+abs(lowertime),amp[ich][isamp]);
                                        }
					if (ich == 1 && channelMax[1] > 150 && channelMax[1] < 300){
					R3AveragePulse[ich]->Fill(isamp*horizontal_interval-trigger_t+abs(lowertime),amp[ich][isamp]);
                                        }
					if (ich == 2  && channelMax[2] > 0 && channelMax[2] < 100){
					R1AveragePulse[ich]->Fill(isamp*horizontal_interval-trigger_t+abs(lowertime),amp[ich][isamp]);
					}
					if (ich == 2  && channelMax[2] > 100 && channelMax[2] < 300){
					R2AveragePulse[ich]->Fill(isamp*horizontal_interval-trigger_t+abs(lowertime),amp[ich][isamp]);
                                        }
					if (ich == 2  && channelMax[2] > 300 && channelMax[2] < 600){
					R3AveragePulse[ich]->Fill(isamp*horizontal_interval-trigger_t+abs(lowertime),amp[ich][isamp]);
                                        }
					if (ich == 3 && channelMax[3] > 0 && channelMax[3] < 50){
					R1AveragePulse[ich]->Fill(isamp*horizontal_interval-trigger_t+abs(lowertime),amp[ich][isamp]);
					}
					if (ich == 3 && channelMax[3] > 50 && channelMax[3] < 200){
					R2AveragePulse[ich]->Fill(isamp*horizontal_interval-trigger_t+abs(lowertime),amp[ich][isamp]);
                                        }
				        if (ich == 3 && channelMax[3] > 200 && channelMax[3] < 400){
					R3AveragePulse[ich]->Fill(isamp*horizontal_interval-trigger_t+abs(lowertime),amp[ich][isamp]);
                                        }
					if (ich == 4 && channelMax[4] > 0 && channelMax[4] < 15){
					R1AveragePulse[ich]->Fill(isamp*horizontal_interval-trigger_t+abs(lowertime),amp[ich][isamp]);	
					}
					if (ich == 4 && channelMax[4] > 15 && channelMax[4] < 35){
					R2AveragePulse[ich]->Fill(isamp*horizontal_interval-trigger_t+abs(lowertime),amp[ich][isamp]);	
                                        }
					if (ich == 4 && channelMax[4] > 35 && channelMax[4] < 60){
					R3AveragePulse[ich]->Fill(isamp*horizontal_interval-trigger_t+abs(lowertime),amp[ich][isamp]);
                                        }
					if (ich == 4 && channelMax[4] > 0 && channelMax[4] < 60){
					R4AveragePulse[ich]->Fill(isamp*horizontal_interval-trigger_t+abs(lowertime),amp[ich][isamp]);
					}
					if (ich == 7 && channelMax[7] > 0 && channelMax[7] < 15){
					R1AveragePulse[ich]->Fill(isamp*horizontal_interval-trigger_t+abs(lowertime),amp[ich][isamp]);
					}
					if (ich == 7 && channelMax[7] > 15 && channelMax[7] < 60){
					R2AveragePulse[ich]->Fill(isamp*horizontal_interval-trigger_t+abs(lowertime),amp[ich][isamp]);
                                        }
					if (ich == 7 && channelMax[7] > 60 && channelMax[7] < 120){
					R3AveragePulse[ich]->Fill(isamp*horizontal_interval-trigger_t+abs(lowertime),amp[ich][isamp]);
                                        }
					if (ich == 7 && channelMax[7] > 0 && channelMax[7] < 120){
					R4AveragePulse[ich]->Fill(isamp*horizontal_interval-trigger_t+abs(lowertime),amp[ich][isamp]);
					}*/
				}
				//if (noiseRMS[ich] > 1.5){
				//	grSnglPulA[ich]->Write();
			//}//SDL Average Pulse	
				for (int del = 0 ; del < 10; del++){
					for(int isamp=0; isamp<samples; isamp++){
					double x = isamp*horizontal_interval-trigger_t+abs(lowertime);
						if (isamp >= Delays[del]){
							AvgSDL[ich][del]->Fill(x, (amp[ich][isamp]-amp[ich][isamp-(Delays[del])]));
							for (int R = 1; R < 4; R++){
								if( channelMax[ich] > SDLRegions[ich][R-1] && channelMax[ich] < SDLRegions[ich][R]){
									RAvgSDL[ich][del][R]->Fill(x, (amp[ich][isamp]-amp[ich][isamp]-Delays[del]));
								}
							}
							


	
	/*				//Regions 1 through 3 cut on amplitude for SDL studies
					if (ich == 0 && channelMax[0] > 0 && channelMax[0] < 75){
					R1AvgSDL[ich][del]->Fill(x, (amp[ich][isamp]-amp[ich][isamp-(Delays[del])]));
                                        }
                                        if (ich == 0 && channelMax[0] > 75 && channelMax[0] < 300){
					R2AvgSDL[ich][del]->Fill(x, (amp[ich][isamp]-amp[ich][isamp-(Delays[del])]));
                                        }
                                        if (ich == 0 && channelMax[0] > 300 && channelMax[0] < 600){
					R3AvgSDL[ich][del]->Fill(x, (amp[ich][isamp]-amp[ich][isamp-(Delays[del])]));
                                        }
                                        if (ich == 1 && channelMax[1] > 0 && channelMax[1] < 50){
					R1AvgSDL[ich][del]->Fill(x, (amp[ich][isamp]-amp[ich][isamp-(Delays[del])]));
                                        }
                                        if (ich == 1 && channelMax[1] > 50 && channelMax[1] < 150){
					R2AvgSDL[ich][del]->Fill(x, (amp[ich][isamp]-amp[ich][isamp-(Delays[del])]));
                                        }
                                        if (ich == 1 && channelMax[1] > 150 && channelMax[1] < 300){
					R3AvgSDL[ich][del]->Fill(x, (amp[ich][isamp]-amp[ich][isamp-(Delays[del])]));
                                        }
                                        if (ich == 2  && channelMax[2] > 0 && channelMax[2] < 100){
					R1AvgSDL[ich][del]->Fill(x, (amp[ich][isamp]-amp[ich][isamp-(Delays[del])]));
                                        }
                                        if (ich == 2  && channelMax[2] > 100 && channelMax[2] < 300){
					R2AvgSDL[ich][del]->Fill(x, (amp[ich][isamp]-amp[ich][isamp-(Delays[del])]));
                                        }
                                        if (ich == 2  && channelMax[2] > 300 && channelMax[2] < 600){
					R3AvgSDL[ich][del]->Fill(x, (amp[ich][isamp]-amp[ich][isamp-(Delays[del])]));
                                        }
                                        if (ich == 3 && channelMax[3] > 0 && channelMax[3] < 50){
					R1AvgSDL[ich][del]->Fill(x, (amp[ich][isamp]-amp[ich][isamp-(Delays[del])]));
                                        }
                                        if (ich == 3 && channelMax[3] > 50 && channelMax[3] < 200){
					R2AvgSDL[ich][del]->Fill(x, (amp[ich][isamp]-amp[ich][isamp-(Delays[del])]));
                                        }
                                        if (ich == 3 && channelMax[3] > 200 && channelMax[3] < 400){
					R3AvgSDL[ich][del]->Fill(x, (amp[ich][isamp]-amp[ich][isamp-(Delays[del])]));
                                        }
                                        if (ich == 4 && channelMax[4] > 0 && channelMax[4] < 15){
					R1AvgSDL[ich][del]->Fill(x, (amp[ich][isamp]-amp[ich][isamp-(Delays[del])]));
                                        }
                                        if (ich == 4 && channelMax[4] > 15 && channelMax[4] < 35){
					R2AvgSDL[ich][del]->Fill(x, (amp[ich][isamp]-amp[ich][isamp-(Delays[del])]));
                                        }
                                        if (ich == 4 && channelMax[4] > 35 && channelMax[4] < 60){
					R3AvgSDL[ich][del]->Fill(x, (amp[ich][isamp]-amp[ich][isamp-(Delays[del])]));
                                        }
					if (ich == 4 && channelMax[4] > 0 && channelMax[4] < 60){
					R4AvgSDL[ich][del]->Fill(x, (amp[ich][isamp]-amp[ich][isamp-(Delays[del])]));
					}
                                        if (ich == 7 && channelMax[7] > 0 && channelMax[7] < 15){
					R1AvgSDL[ich][del]->Fill(x, (amp[ich][isamp]-amp[ich][isamp-(Delays[del])]));
                                        }
                                        if (ich == 7 && channelMax[7] > 15 && channelMax[7] < 60){
					R2AvgSDL[ich][del]->Fill(x, (amp[ich][isamp]-amp[ich][isamp-(Delays[del])]));
                                        }
                                        if (ich == 7 && channelMax[7] > 60 && channelMax[7] < 120){
					R3AvgSDL[ich][del]->Fill(x, (amp[ich][isamp]-amp[ich][isamp-(Delays[del])]));
                                        }
					if (ich == 7 && channelMax[7] > 0 && channelMax[7] < 120){
                                        R4AvgSDL[ich][del]->Fill(x, (amp[ich][isamp]-amp[ich][isamp-(Delays[del])]));
					}*/

						}
					}
				}
		}//Manual Quality events
	}//Channel Loop

	for (int isamp=0;isamp<samples;isamp++){
		grTrigPul->SetPoint(isamp,isamp*horizontal_interval,trigger[0][isamp]);
         	grTrigPulA->SetPoint(isamp,isamp*horizontal_interval-trigger_t,trigger[0][isamp]);
//		AveragePulseFront->SetPoint(isamp,(isamp*horizontal_interval-trigger_t),(channels[0][isamp]-baseline[ich]));
//		AveragePulseBack->SetPoint(
	}

}//event quality 

}//event loop

baselinehisto->Write();
AveragePulseFront->Write();
AveragePulseBack->Write();
AveragePulseMid->Write();


TFile* SDLfile = new TFile("DESYSDLTGraphs_Regions.root", "UPDATE");
    if (!SDLfile->IsOpen()) {
        std::cerr << "Error: Could not open output.root file!" << std::endl;
        return;
     }
for (int ich = 0; ich < 8; ich++){
	pulseFile->cd();
	AveragePulse[ich]->GetXaxis()->SetTitle("time - T_{Trig} (ns)");
	AveragePulse[ich]->GetYaxis()->SetTitle("a (mV)");
	AveragePulse[ich]->Write();
	h_noiseRMS[ich]->Write();
	h_ChAmp[ich]->Write();
	baselineH[ich]->Write();
	for (int del = 0 ; del < 1; del++){//There are 10 elements but only focused on the 1 ns delay
	//double SDLmax = getMaxValue(AvgSDL[ich][del]);
	//AvgSDL[ich][del]->Scale(1/SDLmax);
	TGraph* tempgraph = new TGraph();
	for (int R = 1; R < 4; R++){
	TGraph* tempgraphR = new TGraph();
	tempgraphR->SetName(Form("R%dSDL_Ch%d_del%d_RG610",R,ich,del));
	for (int i =1; i <=RAvgSDL[ich][del][R]->GetNbinsX();i++){
	tempgraphR->SetPoint(i-1,RAvgSDL[ich][del][R]->GetBinCenter(i),-1*(RAvgSDL[ich][del][R]->GetBinContent(i)));
	}
	tempgraphR->Scale(-1/RAvgmax[ich][R]);
	SDLfile->cd();
	tempgraphR->Write();
	}
//	TGraph* tempgraphR1 = new TGraph();
//	TGraph* tempgraphR2 = new TGraph();
//	TGraph* tempgraphR3 = new TGraph();
//	TGraph* tempgraphR4 = new TGraph();
	tempgraph->SetName(Form("SDL_Ch%d_del%d_RG610",ich,del));
//	tempgraphR1->SetName(Form("R1SDL_Ch%d_del%d_RG610",ich,del));
//	tempgraphR2->SetName(Form("R2SDL_Ch%d_del%d_RG610",ich,del));
//	tempgraphR3->SetName(Form("R3SDL_Ch%d_del%d_RG610",ich,del));
//	tempgraphR4->SetName(Form("R4SDL_Ch%d_del%d_RG610",ich,del));
//		for(int i =1; i <=AvgSDL[ich][del]->GetNbinsX();i++){
//		tempgraph->SetPoint(i-1,AvgSDL[ich][del]->GetBinCenter(i),-1*(AvgSDL[ich][del]->GetBinContent(i)));
//			}
//		for(int i =1; i <=R1AvgSDL[ich][del]->GetNbinsX();i++){
  //              tempgraphR1->SetPoint(i-1,R1AvgSDL[ich][del]->GetBinCenter(i),-1*(R1AvgSDL[ich][del]->GetBinContent(i)));
    //                    }
//		for(int i =1; i <=R2AvgSDL[ich][del]->GetNbinsX();i++){
  //              tempgraphR2->SetPoint(i-1,R2AvgSDL[ich][del]->GetBinCenter(i),-1*(R2AvgSDL[ich][del]->GetBinContent(i)));
    //                    }
//		for(int i =1; i <=R3AvgSDL[ich][del]->GetNbinsX();i++){
  ///              tempgraphR3->SetPoint(i-1,R3AvgSDL[ich][del]->GetBinCenter(i),-1*(R3AvgSDL[ich][del]->GetBinContent(i)));
     //                   }
//		for(int i =1; i <=R4AvgSDL[ich][del]->GetNbinsX();i++){
///		tempgraphR4->SetPoint(i-1,R4AvgSDL[ich][del]->GetBinCenter(i),-1*(R4AvgSDL[ich][del]->GetBinContent(i)));
//		}
	double SDLmax = getMaxValue(tempgraph);
//	double R1SDLmax = getMaxValue(tempgraphR1);
//	double R2SDLmax = getMaxValue(tempgraphR2);
//	double R3SDLmax = getMaxValue(tempgraphR3);
//	double R4SDLmax = getMaxValue(tempgraphR4);
        tempgraph->Scale(-1/SDLmax);
//	tempgraphR1->Scale(-1/R1SDLmax);
//	tempgraphR2->Scale(-1/R2SDLmax);
///	tempgraphR3->Scale(-1/R3SDLmax);
///	tempgraphR4->Scale(-1/R4SDLmax);
	SDLfile->cd();
	tempgraph->Write();	
//	tempgraphR1->Write();
//	tempgraphR2->Write();
//	tempgraphR3->Write();
//	tempgraphR4->Write();
	}
	
	//double Avgmax = getMaxValue(AveragePulse[ich]);
	//AveragePulse[ich]->Scale(1/Avgmax);
	AveragePulseGraph[ich]=Convert(AveragePulse[ich]);
//	AvgPGraphR[ich][R] = Convert(RAveragePulse[ich][R]);
//	AvgPGraphR2[ich] = Convert(R2AveragePulse[ich]);
//	AvgPGraphR3[ich] = Convert(R3AveragePulse[ich]);
//	AvgPGraphR4[ich] = Convert(R4AveragePulse[ich]);
        AveragePulseGraph[ich]->SetName(Form("AveragePulseGraph_Ch%d", ich));
//	AvgPGraphR1[ich]->SetName(Form("R1AveragePulseGraph_Ch%d", ich));
//	AvgPGraphR2[ich]->SetName(Form("R2AveragePulseGraph_Ch%d", ich));
//	AvgPGraphR3[ich]->SetName(Form("R3AveragePulseGraph_Ch%d", ich));
//	AvgPGraphR4[ich]->SetName(Form("R4AveragePulseGraph_Ch%d",ich));
	double Avgmax = getMaxValue(AveragePulseGraph[ich]);
//	double R1Avgmax = getMaxValue(AvgPGraphR1[ich]);
//	double R2Avgmax = getMaxValue(AvgPGraphR2[ich]);
//	double R3Avgmax = getMaxValue(AvgPGraphR3[ich]);
//	double R4Avgmax = getMaxValue(AvgPGraphR4[ich]);
        AveragePulseGraph[ich]->Scale(1/Avgmax);
//	AvgPGraphR1[ich]->Scale(1/R1Avgmax);
//	AvgPGraphR2[ich]->Scale(1/R2Avgmax);
//	AvgPGraphR3[ich]->Scale(1/R3Avgmax);
//	AvgPGraphR4[ich]->Scale(1/R3Avgmax);
	for (int R = 1; R < 5 ; R++){
		AvgPGraphR[ich][R] = Convert(RAveragePulse[ich][R]);
		AvgPGraphR[ich][R]->SetName(Form("R%dAveragePulseGraph_Ch%d",R, ich));
		double RAvgmax[ich][R] = getMaxValue(AvgPGraphR[ich][R]);
		AvgPGraphR[ich][R]->Scale(1/RAvgmax[ich][R]);
		pulseFile->cd();
		AvgPGraphR[ich][R]->Write();
	}
	//pulseFile->cd();
	AveragePulseGraph[ich]->Write();
	//AvgPGraphR1[ich]->Write();
///	AvgPGraphR2[ich]->Write();
//	AvgPGraphR3[ich]->Write();	
//	AvgPGraphR4[ich]->Write();
	//Willing to bet that a seg fault happens since I don't fill AvgPGraph for any channels other than 4 and 7. Will have to fill the rest with a single data point perhaps
}
SDLfile->Close();
pulseFile->cd();
   h_noiseRMStotal->Write();
   h_eventAmp->Write();
   trigger_Amp->Write();
   triggerRMS->Write();
   trigger_time->Write();
   pulseFile->Write();
   pulseFile->Close();
   std::cout<<"Done :)"<<std::endl;

}//End


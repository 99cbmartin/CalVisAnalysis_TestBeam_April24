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
void Analysis::Loop(std::string outputFileName)
{
Int_t icol[8] = { kOrange+10, kOrange-3, kSpring-1, kAzure+6, kAzure-1, kViolet-1, kBlack, kRed+1};//first and last color too close
TH1D* triggerRMS = new TH1D("triggerRMS","Trigger RMS;Value;Counts",100,-1,1);
TH1D* trigger_Amp = new TH1D("trigger_Amp", "Trigger Amplitude; Value;Counts",100,400,500);
TH1D* trigger_time = new TH1D("trigger_time","Trigger Time;time(ns);Counts",100,0,10);
TH1D* h001 = new TH1D("EventRMS","Baseline Noise of Events;RMS;Counts",100,0,100);
TCanvas* triggerAlignment = new TCanvas("triggerAlignment","Aligned Triggers",800,600);
TCanvas* SingleEvents = new TCanvas("SingleEvents", "Aligned Single Events",800,600);
TProfile* AveragePulse[8]; 
TH1D* h_noiseRMS[8];
TH1D* h_noiseRMStotal = new TH1D("AllEvtRMS","RMS of All Events",100,0,1.5);
TH1D* h_eventAmp = new TH1D("EvtAmps","Amplitudes of Events",900,-500,500);
TH1D* h_ChAmp[8];

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
auto chstr = std::to_string(ich);
h_noiseRMS[ich] = new TH1D(Form("RMS_Noise_Ch%d",ich),Form("RMS_Noise_Ch%d",ich),200,-1.5,1.5);
h_ChAmp[ich] = new TH1D(Form("ChAmp%d",ich),Form("ChAmp%d",ich),100,0,100);
}

  TFile* pulseFile = new TFile(outputFileName.c_str(),"recreate");
  TGraph * grSnglPul[8];
  TGraph * grSnglPulA[8];
  TGraph* grTrigPul;
  TGraph* grTrigPulA;
  if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
//begin event loop
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
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
//std::cout<<"trigger amplitudes are: "<<trigger_amp<<std::endl;
noise = TMath::RMS(basevalues.size(),basevalues.data());
for (int isamp = 0; isamp < samples; isamp++){
	if(trigger[0][isamp] < trigamp && trigger[0][isamp+1] >= trigamp ||
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
	   
		if (fabs(channels[ich][isamp]) > fabs(channelMax[ich])){

                        channelMax[ich] = fabs(channels[ich][isamp]-baseline[ich]);
                }
        }
//isOK = isOK && (channelMax[ich] > 3* noiseRMS[ich]);//Event Amplitude quality check, possible to use 4 or 5
isOK = isOK && noise < 1 && noise > 0.5 && trigger_amp > 440 && trigger_amp < 455;//Trigger conditions for quality events
isOK = isOK && channelMax[ich] < 10 && noiseRMS[ich] <1.1;//SiPM conditions for quality events
}//channel loop end
//===================================\\
//           Pulse Analysis          \\
//===================================\\
//Actual Analysis of Pulse Shapes and Timing
num2++;
if (isOK){//Quality Events only beyond this point
if (channelMax[0] > 10){
std::cout<<"The max of channel 0 is: "<<channelMax[2]<<" on event: "<<jentry<<std::endl;
}
num++;
average += channelMax[2];
for (int ich=0;ich<8; ich++){

//AveragePulse[ich]->SetLineColor(icol[ich]);
grSnglPulA[ich]->SetLineColor(icol[ich]);
grSnglPulA[ich]->GetXaxis()->SetTitle("time - T_{Trig} (ns)");
grSnglPulA[ich]->GetYaxis()->SetTitle("a (mV)");
Float_t ampchan = *std::max_element(channels[ich],channels[ich]+1024);
//std::cout<<"  Channel "<<ich<<" max "<<ampchan<<std::endl;
h_eventAmp->Fill(channelMax[ich]);
h_ChAmp[ich]->Fill(channelMax[ich]);

for (int isamp=0;isamp<samples;isamp++){
          AveragePulse[ich]->Fill(isamp*horizontal_interval-trigger_t+abs(lowertime),channels[ich][isamp]-baseline[ich]);
          grSnglPulA[ich]->SetPoint(isamp,(isamp*horizontal_interval-trigger_t),(channels[ich][isamp]-baseline[ich]));
}
if (jentry < 100){
grSnglPulA[ich]->Write();
}


}//Channel Loop
//if (channelMax[0] > 100) {num++;}
//if (channelMax[0] < 400) {num2++;}
for (int isamp=0;isamp<samples;isamp++){
	  grTrigPul->SetPoint(isamp,isamp*horizontal_interval,trigger[0][isamp]);
          grTrigPulA->SetPoint(isamp,isamp*horizontal_interval-trigger_t,trigger[0][isamp]);
}

if (jentry == 0){
triggerAlignment->cd();
grTrigPulA->Draw();
SingleEvents->cd();
//grSnglPulA->Draw();
}
//if (jentry >6) break;
//triggerAlignment->cd();
//grTrigPulA->SetLineColor(icol[jentry]);
//grTrigPulA->Draw("SAME");
//SingleEvents->cd();
//grSnglPulA->Draw("SAME");
 //grTrigPul->Write();
 //grTrigPulA->Write();
//grSnglPulA->Write();
//SingleEvents->Write();
}//event quality 
   }//event loop
std::cout<<"After cut "<<num<<"before cut: "<<num2<<std::endl;
for (int ich = 0; ich < 8; ich++){
AveragePulse[ich]->GetXaxis()->SetTitle("time - T_{Trig} (ns)");
AveragePulse[ich]->GetYaxis()->SetTitle("a (mV)");
AveragePulse[ich]->Write();
h_noiseRMS[ich]->Write();
h_ChAmp[ich]->Write();
}
h_noiseRMStotal->Write();
h_eventAmp->Write();
//   triggerAlignment->Write();
   trigger_Amp->Write();
   triggerRMS->Write();
   trigger_time->Write();
   pulseFile->Write();
   pulseFile->Close();
   std::cout<<"Wrote output file with single pulses"<<std::endl;
}


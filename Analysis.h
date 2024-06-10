//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Apr 15 02:14:05 2024 by ROOT version 6.30/02
// from TTree tree/DRS Data
// found on file: ../bc-35v-tbbox.root
//////////////////////////////////////////////////////////

#ifndef Analysis_h
#define Analysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
// Header file for the classes stored in the TTree if any.

class Analysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
  // Float_t         vertical_gain[8];
  // Float_t         vertical_offset[8];
   Float_t         horizontal_interval;
  // Double_t        horizontal_offset;
   Int_t           event;
   Double_t        trigger_time;
  // Double_t        trigger_offset;
  // Double_t        time[1024];
  // Int_t           samples;
   Float_t         channels[10][1024];
   Float_t	   trigger[2][1024];

   // List of branches
   //TBranch        *b_vertical_gain;   //!
  // TBranch        *b_vertical_offset;   //!
   TBranch        *b_horizontal_interval;   //!
 //  TBranch        *b_horizontal_offset;   //!
   TBranch        *b_event;   //!
   TBranch        *b_trigger_time;   //!
   //TBranch        *b_trigger_offset;   //!
   //TBranch        *b_time;   //!
   //TBranch        *b_samples;   //!
   TBranch        *b_channels;   //!
   TBranch	  *b_trigger; //!

   Analysis(TTree *tree=0);
   virtual ~Analysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(std::string outputFileName);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Analysis_cxx
Analysis::Analysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
std::cerr << "Error: No TTree of TChain" << std::endl;
std::exit(EXIT_FAILURE);
//      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../store/shared/calvision/run_285/outfile_LG.root");
  //    if (!f || !f->IsOpen()) {
    //     f = new TFile("../../store/shared/calvision/run_285/outfile_LG.root");
    //  }
     // f->GetObject("tree",tree);

   }
   Init(tree);
}

Analysis::~Analysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Analysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Analysis::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Analysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

  // fChain->SetBranchAddress("vertical_gain", vertical_gain, &b_vertical_gain);
   //fChain->SetBranchAddress("vertical_offset", vertical_offset, &b_vertical_offset);
   fChain->SetBranchAddress("horizontal_interval", &horizontal_interval, &b_horizontal_interval);
  // fChain->SetBranchAddress("horizontal_offset", &horizontal_offset, &b_horizontal_offset);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("trigger_time", &trigger_time, &b_trigger_time);
  // fChain->SetBranchAddress("trigger_offset", &trigger_offset, &b_trigger_offset);
  // fChain->SetBranchAddress("time", time, &b_time);
 //  fChain->SetBranchAddress("samples", &samples, &b_samples);
   fChain->SetBranchAddress("channels", channels, &b_channels);
   fChain->SetBranchAddress("trigger",trigger,&b_trigger);
   Notify();
}

Bool_t Analysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Analysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Analysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Analysis_cxx

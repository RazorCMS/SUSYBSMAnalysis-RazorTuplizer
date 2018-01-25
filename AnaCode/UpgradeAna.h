//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jan 19 18:54:53 2018 by ROOT version 6.10/05
// from TTree RazorEvents/selected miniAOD information
// found on file: razorNtuple.root
//////////////////////////////////////////////////////////

#ifndef UpgradeAna_h
#define UpgradeAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <cmath>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

#include "NtupleVariables.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include<string>
#include "TF1.h"
#include "TAxis.h"

class UpgradeAna : public NtupleVariables{
public :
  UpgradeAna(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data");
   virtual ~UpgradeAna();
   
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *);
  void     BookHistogram(const char *);

  TFile *oFile; 

  // Define histograms here
  TH1D *h_dR;
  TH1D *h_delPt;
  TH1D *h_leadingMuonPt;
  TH1D *h_leadingMuonTime;
  TH1D *h_MatchedTrackPt;
  TH1D *h_ChargedFractionR04;
  TH1D *h_ChargedFractionR04_dT;
  TH1D *h_CMSchargedIso04;
  TH1D *h_diff;
};

#endif

#ifdef UpgradeAna_cxx



void UpgradeAna::BookHistogram(const char *outFileName) {

  oFile = new TFile(outFileName, "recreate");
  h_dR=new TH1D("deltaR","deltaR",100,0.,0.1);
  h_delPt=new TH1D("deltaPt","deltaPt",100,0.,100.);
  h_leadingMuonPt=new TH1D("mu1Pt","mu1Pt",50,0.,500.);
  h_MatchedTrackPt=new TH1D("Tr1Pt","Tr1Pt",50,0.,500.);
  h_leadingMuonTime=new TH1D("mu1T","mu1T",100,2.,-2.);
  h_ChargedFractionR04=new TH1D("ChargedFrac","ChargedFrac",50,0.,1.);
  h_ChargedFractionR04_dT=new TH1D("ChargedFrac_dT","ChargedFrac_dT",50,0.,1.);
  h_CMSchargedIso04=new TH1D("ChargedFracCMS","ChargedFracCMS",50,0.,1.);
  h_diff=new TH1D("h_diff","h_diff",100,-1.,1.);
}

UpgradeAna::UpgradeAna(const TString &inputFileList, const char *outFileName, const char* dataset)
{
TChain *tree = new TChain("ntuples/RazorEvents");

  if( ! FillChain(tree, inputFileList) ) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;
  
  NtupleVariables::Init(tree);
  BookHistogram(outFileName);
  }
}

Bool_t UpgradeAna::FillChain(TChain *chain, const TString &inputFileList) {

  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if(!infile.is_open()) {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;
  while(1) {
    infile >> buffer;
    if(!infile.good()) break;
    std::cout << "Adding tree from " << buffer.c_str() << std::endl;                                                              
    chain->Add(buffer.c_str());
  }
  std::cout << "No. of Entries in this tree : " << chain->GetEntries() << std::endl;
  
  return kTRUE;
 }
UpgradeAna::~UpgradeAna()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
  oFile->cd();
  oFile->Write();
  oFile->Close();
}


Long64_t UpgradeAna::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
    if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      // Notify();
   }
   return centry;
}


#endif // #ifdef UpgradeAna_cxx

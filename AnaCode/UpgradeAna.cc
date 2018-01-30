#define UpgradeAna_cxx
#include "UpgradeAna.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <cstring>
#include<string>
#include<fstream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <cmath>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include<string>
#include "TF1.h"
#include "TAxis.h"

using namespace std;
int main(int argc, char* argv[])
{

  if(argc < 3) {
    std::cerr << "Please give 3 arguments " << "runList " << " " << "outputFileName" << " " << "dataset"<<std::endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];
  
  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if(!infile.is_open()) {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }
  TChain *treeChain = new TChain("ntuples/RazorEvents");
  std::cout << "TreeUtilities : FillChain " << std::endl;
  while(1) {
    infile >> buffer;
    if(!infile.good()) break;
    std::cout << "Adding tree from " << buffer.c_str() << std::endl;      
    treeChain->Add(buffer.c_str());
  }
  std::cout << "No. of Entries in this tree : " << treeChain->GetEntries() << std::endl;

  UpgradeAna UpgradeP2(treeChain);
  cout << "dataset " << data << " " << endl;
  UpgradeP2.Loop(outFileName);
  return 0;
}

void UpgradeAna::Loop(const char *outFileName)
{

  
  TFile *oFile = new TFile(outFileName, "recreate");

  BookHistograms();
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      float pt=-1.;
      int itrMu=-99;
      for(int mu=0;mu<nMuons;mu++){
	if(jentry==244)cout<<muonE[mu]<<" "<<muonPt[mu]<<" "<<muonEta[mu]<<" "<<muonPhi[mu]<< " "<<muon_d0[mu]<< " "<<muon_dZ[mu]<<endl;
	if(muonPt[mu]>pt){
	  itrMu=mu;
	  pt=muonPt[mu];
	}
      }
      float dRmin=99.0;
      int iTr=-99;
      float dR,dT;
    
      for(int itTrack= 0;itTrack<allTrackT->size();itTrack++){
	dR=DeltaR(muonEta[itrMu],muonPhi[itrMu],(*allTrackEta)[itTrack],(*allTrackPhi)[itTrack]);
	//cout<<dR<<endl;
	if(dR<dRmin){
	  dRmin=dR;
	  //cout<<dRmin<<endl;
	  //cout<<itTrack<<" "<<(*allTrackPt)[itTrack]<< " "<<(*allTrackPhi)[itTrack]<<endl;
	  iTr=itTrack;
	}
      }
      if(jentry==244)cout<<"====================\n";
          if(!dRmin && itrMu!=-99 && iTr!=-99){
	    float deltaPt=muonPt[itrMu]-(*allTrackPt)[iTr];
	    //h_dR->Fill(fabs(dRmin));
	    //h_delPt->Fill(fabs(deltaPt));
	    cout<<muonPt[itrMu]<<endl;
	    h_leadingMuonPt->Fill(muonPt[itrMu]);
	    h_leadingMuonTime->Fill((*allTrackdT)[iTr]);
	    //      cout<<dRmin<<" "<<itrMu<<" "<<muonPt[itrMu]<<" "<<iTr<<" "<<(*allTrackPt)[iTr]<<endl;
	    h_MatchedTrackPt->Fill((*allTrackPt)[iTr]);
	    float SumChargedPt=0.0;
	    float t_SumChargedPt=0.0;
	    dT=9999.;
	    for(int itTrack= 0;itTrack<allTrackT->size();itTrack++){
	    if(itTrack==iTr)continue;
	    dR=DeltaR(muonEta[itrMu],muonPhi[itrMu],(*allTrackEta)[itTrack],(*allTrackPhi)[itTrack]);
	    dT=(*allTrackdT)[iTr]-(*allTrackdT)[itTrack];
	    if(dR<0.4)SumChargedPt+=(*allTrackPt)[itTrack];
	    if(dR<0.4 && fabs(dT)<0.090)t_SumChargedPt+=(*allTrackPt)[itTrack];
	    }
	    h_ChargedFractionR04->Fill(SumChargedPt/(*allTrackPt)[iTr]);
	    h_ChargedFractionR04_dT->Fill(t_SumChargedPt/(*allTrackPt)[iTr]);
	    h_diff->Fill((muon_chargedIso[itrMu]+muon_pileupIso[itrMu])-SumChargedPt);
	    h_CMSchargedIso04->Fill((muon_chargedIso[itrMu]+muon_pileupIso[itrMu])/(*allTrackPt)[iTr]);
	    }


   }

   oFile->Write();
   oFile->Close();
}

void UpgradeAna::BookHistograms(){
  
  h_dR=new TH1D("deltaR","deltaR",100,0.,0.1);
  h_delPt=new TH1D("deltaPt","deltaPt",100,0.,100.);
  h_leadingMuonPt=new TH1D("mu1Pt","mu1Pt",50,0.,500.);
  h_MatchedTrackPt=new TH1D("Tr1Pt","Tr1Pt",50,0.,500.);
  h_leadingMuonTime=new TH1D("mu1T","mu1T",100,2.,-2.);
  h_ChargedFractionR04=new TH1D("ChargedFrac","ChargedFrac",50,0.,1.);
  h_ChargedFractionR04_dT=new TH1D("ChargedFrac_dT","ChargedFrac_dT",50,0.,1.);
  h_CMSchargedIso04=new TH1D("ChargedFracCMS","ChargedFracCMS",50,0.,1.);
  h_diff=new TH1D("h_diff","h_diff",100,100.,0.);

}

#define UpgradeAna_cxx
#include "UpgradeAna.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <cstring>
#include<string>
int main(int argc, char* argv[])
{

  if(argc < 3) {
    cerr << "Please give 3 arguments " << "runList " << " " << "outputFileName" << " " << "dataset"<<endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];
  UpgradeAna UpgradeP2(inputFileList, outFileName, data);
  cout << "dataset " << data << " " << endl;
  UpgradeP2.EventLoop(data);

  return 0;
}
void UpgradeAna::EventLoop(const char *data)
{

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    float pt=-1.;
    int itrMu=-99;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    for(int mu=0;mu<nMuons;mu++){
      //if(jentry<1)cout<<muonE[mu]<<" "<<muonPt[mu]<<" "<<muonEta[mu]<<" "<<muonPhi[mu]<< " "<<muon_d0[mu]<< " "<<muon_dZ[mu]<<endl;
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
    //cout<<"====================\n";
    if(!dRmin && itrMu!=-99){
      float deltaPt=muonPt[itrMu]-(*allTrackPt)[iTr];
      h_dR->Fill(fabs(dRmin));
      h_delPt->Fill(fabs(deltaPt));
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
      h_diff->Fill((muon_chargedIso[itrMu]+muon_pileupIso[itrMu])/(*allTrackPt)[iTr]-SumChargedPt/(*allTrackPt)[iTr]);
      h_CMSchargedIso04->Fill((muon_chargedIso[itrMu]+muon_pileupIso[itrMu])/(*allTrackPt)[iTr]);
    }

    //float deltaZ=muon_dZ[mu]-(*allTrackdZ)[itTrack];
    //  float deltaPt=muonPt[mu]-(*allTrackPt)[itTrack];
    // if(fabs(deltaZ)<0.01){
    //	h_delZ->Fill(fabs(deltaZ));
    //	h_delPt->Fill(fabs(deltaPt));
    //  }
      //if(jentry<1)cout<<(*allTrackdZ)[itTrack]<<endl;
	
	//DeltaR(muonEta[mu],muonPhi[mu]
  //    }

//}
     //if(jentry<10)cout<<"==============================\n";
  }
}


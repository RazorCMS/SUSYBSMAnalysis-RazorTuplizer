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
#include "TRandom3.h"
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
  //std::cout << "dataset " << data << " " << endl;
  UpgradeP2.Loop(outFileName);
  //std::cout << "done\n";
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
      //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx  Order muons according to pt  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      for(int mu=0;mu<nMuons;mu++){
	//if(jentry==244)cout<<muonE[mu]<<" "<<muonPt[mu]<<" "<<muonEta[mu]<<" "<<muonPhi[mu]<< " "<<muon_d0[mu]<< " "<<muon_dZ[mu]<<endl;
	if(muonPt[mu]>pt){
	  itrMu=mu;
	  pt=muonPt[mu];
	}
      }
      float dRmin=99.0;
      int iTr=-99;
      float dR,dT;
      //cout<<"Entering Loop\n";



      // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx find muon track xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      for(int itTrack= 0;itTrack<allTrackEta->size();itTrack++){
	dR=DeltaR(muonEta[itrMu],muonPhi[itrMu],(*allTrackEta)[itTrack],(*allTrackPhi)[itTrack]);
	//cout<<dR<<endl;
	if(dR<dRmin){
	  dRmin=dR;
	  //cout<<dRmin<<endl;
	  //cout<<itTrack<<" "<<(*allTrackPt)[itTrack]<< " "<<(*allTrackPhi)[itTrack]<<endl;
	  iTr=itTrack;
	}
      }
      double myIso=-9999;
      double myIso_dT=-9999;
      if(!dRmin && itrMu!=-99 && iTr!=-99 && muonPt[itrMu]>20. && fabs(muonEta[itrMu])<2.4){

	h_leadingMuonPt->Fill(muonPt[itrMu]);
	h_leadingMuonTime->Fill((*allTrackdT)[iTr]);
	h_MatchedTrackPt->Fill((*allTrackPt)[iTr]);
	float SumChargedPt=0.0;
	float t_SumChargedPt=0.0;
	dT=9999.;
	//cout<<"Entering Loop 2\n";
	for(int itTrack= 0;itTrack<allTrackEta->size();itTrack++){
	  if(itTrack==iTr)continue;
	  dR=DeltaR(muonEta[itrMu],muonPhi[itrMu],(*allTrackEta)[itTrack],(*allTrackPhi)[itTrack]);
	  dT=(*allTrackdT)[iTr]-(*allTrackdT)[itTrack];
	  if(dR<0.4)h_trackDz->Fill((*allTrackDz)[itTrack]); 
	  if(dR<0.4 && fabs((*allTrackDz)[itTrack])<0.5/* && fabs((*allTrackDxy)[itTrack])<0.1 /*&& fabs((*allTrackX)[itTrack]-pvX)<0.1*/){
	    
	    //if(dR<0.4 && ((*allTrackZ)[itTrack]==pvZ) && ((*allTrackY)[itTrack]==pvY) && ((*allTrackX)[itTrack]==pvX)){
	    SumChargedPt+=(*allTrackPt)[itTrack];
	    h_trackX->Fill((*allTrackX)[itTrack]-pvX);
	    h_trackY->Fill((*allTrackY)[itTrack]-pvY);
	    h_trackZ->Fill((*allTrackZ)[itTrack]-pvZ);
	    //h_trackDz->Fill((*allTrackDz)[itTrack]);
	    //h_trackDxy->Fill((*allTrackDxy)[itTrack]);
	    h_dT->Fill(dT);
	  }
	  if(dR<0.4 && fabs(dT)<0.030 && fabs((*allTrackDz)[itTrack])<0.5 /*&& fabs((*allTrackY)[itTrack]-pvY)<0.1 && fabs((*allTrackX)[itTrack]-pvX)<0.1*/)t_SumChargedPt+=(*allTrackPt)[itTrack];
	}
	//cout<<"Loop 2 done\n";
	myIso=SumChargedPt/(*allTrackPt)[iTr];
	myIso_dT=t_SumChargedPt/(*allTrackPt)[iTr];
	h_ChargedFractionR04->Fill(myIso);
	h_ChargedFractionR04_dT->Fill(myIso_dT);
	h_diff->Fill((muon_chargedIso[itrMu]-SumChargedPt));
	h_CMSchargedIso04->Fill((muon_chargedIso[itrMu])/(*allTrackPt)[iTr]);

	Int_t itGenMu=-99;
	pt=-1.;
	for(int i=0;i<nGenParticle;i++){
	  if(fabs(gParticleId[i])==13 && gParticleStatus[i]==1){
	    pt=gParticlePt[i];
	    if(fabs(muonPt[itrMu]-pt)<2. && fabs(gParticleEta[i]-muonEta[itrMu])<0.001 && fabs(gParticlePhi[i]-muonPhi[itrMu])<0.001){
	      itGenMu=i;
	      h_genMuon->Fill(gParticlePt[itGenMu]);
	      h_genMuon_vtx->Fill(nPV); 
	    if(myIso_dT<0.15 && myIso_dT!=-9999){
	      h_recoMuon_dT->Fill(gParticlePt[itGenMu]);
	      h_recoMuon_dTvtx->Fill(nPV);
	    }
	    if(myIso<0.15 && myIso!=-9999){
	      h_recoMuon->Fill(gParticlePt[itGenMu]);
	      h_recoMuon_vtx->Fill(nPV);
	    }
	    break;
	    }
	    //cout<<"ptGen: "<<pt<<" gPArtStat: "<<gParticleStatus[i]<<" EtaGen:"<<gParticleEta[i]<<" gParticlePhi: "<<gParticlePhi[i]<<" muonPt[itrMu]:"<<muonPt[itrMu]<<" muonEta:"<<muonEta[itrMu]<<" muonPhi: "<<muonPhi[itrMu]<<endl;
	  }
	}
	
      }
   

      /******************************************************************************************************************************************************************
                                                         Photon ID calculation
      *******************************************************************************************************************************************************************/
      TLorentzVector g1,g2,hcand;
      vector <float>  t1,t2;
      float tmin=100;
      int pvCand=-99;

      int itrPho=-99;
      pt=-1.;
      for(int pho=0;pho<nPhotons;pho++){ // Find leading photon=========
	//cout<<phoE[pho]<<" "<<phoPt[pho]<<" "<<phoEta[pho]<<" "<<phoPhi[pho]<<endl;
	//cout<<pho_superClusterEnergy[pho]<<" "<<phoPt[pho]<<" "<<pho_superClusterEta[pho]<<" "<<pho_superClusterPhi[pho]<<endl;
	//cout<<"=====================\n";
	
	if(phoPt[pho]>pt){
	  itrPho=pho;
	  pt=phoPt[pho];
	}
      }
      int itrPho2=-99;
      pt=-1;
      for(int pho=0;pho<nPhotons && pho!=itrPho;pho++){ // Find 2nd photon
	if(phoPt[pho]>pt){
          itrPho2=pho;
          pt=phoPt[pho];
        }
      }
      if(itrPho!=-99 && itrPho2!=-99){
	g1.SetPtEtaPhiE(phoPt[itrPho],phoEta[itrPho],phoPhi[itrPho],phoE[itrPho]);
	g2.SetPtEtaPhiE(phoPt[itrPho2],phoEta[itrPho2],phoPhi[itrPho2],phoE[itrPho2]);
	hcand=g1+g2;
	/*if(fabs(g2.Eta())<1.48 && fabs(g1.Eta())<1.48)*/{
	  //cout<<"P1: "<<pho_superClusterSeedT[itrPho]<<" "<<pho_superClusterSeedZ[itrPho]<<" "<<pho_superClusterSeedX[itrPho]<<" "<<pho_superClusterSeedY[itrPho]<<" "<<pow(pho_superClusterSeedX[itrPho],2)+pow(pho_superClusterSeedY[itrPho],2)<<endl;
	  //cout<<"P2: "<<pho_superClusterSeedT[itrPho2]<<" "<<pho_superClusterSeedZ[itrPho2]<<" "<<pho_superClusterSeedX[itrPho2]<<" "<<pho_superClusterSeedY[itrPho2]<<" "<<pow(pho_superClusterSeedX[itrPho2],2)+pow(pho_superClusterSeedY[itrPho2],2)<<endl;}
	  if(g1.Pt()>20. && g2.Pt()>20. && fabs(g2.Eta())<2.4 && fabs(g1.Eta())<2.4){
	    h_diPhoPt->Fill(hcand.Pt());
	    h_diPhoM->Fill(hcand.M());
	    h_diPhoEta->Fill(hcand.Eta());
	    h_diPhoPhi->Fill(hcand.Phi());
	
	    h_PhoEtaDiff->Fill(g1.Eta()-g2.Eta());
	    h_PhoPhiDiff->Fill(DeltaPhi(g1.Phi(),g2.Phi()));
	    float pho_t1,pho_t2;
	    float dist1 = pow((pho_superClusterSeedZ[itrPho]-genVertexZ),2.)+pow((pho_superClusterSeedX[itrPho]-genVertexX),2.)+pow((pho_superClusterSeedY[itrPho]-genVertexY),2.);
	    float dist2 = pow((pho_superClusterSeedZ[itrPho2]-genVertexZ),2.)+pow((pho_superClusterSeedX[itrPho2]-genVertexX),2.)+pow((pho_superClusterSeedY[itrPho2]-genVertexY),2.);
	    //float dist1 = pow((pho_superClusterSeedZ[itrPho]-pvAllZ[pv]),2.);
	    //float dist2 = pow((pho_superClusterSeedZ[itrPho2]-pvAllZ[pv]),2.);
	    TRandom3 randomPhotonTime1(1111);
	    
	    pho_t1=randomPhotonTime1.Gaus(genVertexT+sqrt(dist1)*0.01*pow(10,9)/(2.99792458*pow(10,8)),0.03);
	    TRandom3 randomPhotonTime2(1111);
	    pho_t2=randomPhotonTime2.Gaus(genVertexT+sqrt(dist2)*0.01*pow(10,9)/(2.99792458*pow(10,8)),0.03);
	    //cout<<pho_t1<<" "<<pho_t2<<endl;
	    //============================== finding photon vertex====================================
	    for(int pv=0;pv<nPVAll;pv++){
	      float temp1,temp2,d1,d2;
	    
	      d1 = pow((pho_superClusterSeedZ[itrPho]-pvAllZ[pv]),2.)+pow((pho_superClusterSeedX[itrPho]-pvAllX[pv]),2.)+pow((pho_superClusterSeedY[itrPho]-pvAllY[pv]),2.);
	      d2 = pow((pho_superClusterSeedZ[itrPho2]-pvAllZ[pv]),2.)+pow((pho_superClusterSeedX[itrPho2]-pvAllX[pv]),2.)+pow((pho_superClusterSeedY[itrPho2]-pvAllY[pv]),2.);
	      temp1=sqrt(d1)*0.01*pow(10,9)/(2.99792458*pow(10,8))+pho_t1;
	      temp2=sqrt(d2)*0.01*pow(10,9)/(2.99792458*pow(10,8))+pho_t2;
	      //if(fabs(temp1-temp2)<0.1)cout<<pvAllZ[pv]<<" "<<pvAllT[pv]<<" "<<temp1<< " "<<temp2<<" "<<hcand.M()<<" "<<hcand.Pt()<<" "<<hcand.Phi()<< " "<<hcand.Eta()<<endl;
	      if(fabs(temp1-temp2)<0.1 && fabs(temp1-temp2)<tmin){
		tmin=fabs(temp1-temp2);
		pvCand=pv;
	      }
	    }
	    if(pvCand!=-99 /*&& hcand.M()<140. && hcand.M()>110.*/){
	      
	      h_diPhoPVIndex->Fill(pvCand);
	      if(fabs(g1.Eta()-g2.Eta())<0.8)h_genVsRecodiff_etal0p8->Fill(genVertexZ-pvAllZ[pvCand]);
	      if(fabs(g1.Eta()-g2.Eta())>0.8)h_genVsRecodiff_etag0p8->Fill(genVertexZ-pvAllZ[pvCand]);
	      if(fabs(g1.Eta()-g2.Eta())<0.8)h_genVspvdiff_etal0p8->Fill(genVertexZ-pvZ);
	      if(fabs(g1.Eta()-g2.Eta())>0.8)h_genVspvdiff_etag0p8->Fill(genVertexZ-pvZ);
	      h_tMinVsDeta->Fill(fabs(g1.Eta()-g2.Eta()),tmin);
	      h_diPhoMwithVtx->Fill(hcand.M());
	      
	      //---------------------------------------------------------leading photon---------------------------------------------------------                                                          
	      
	      float SumChargedPt=0.0;
	      float t_SumChargedPt=0.0;
	      for(int itTrack= 0;itTrack<allTrackEta->size();itTrack++){
		dR=DeltaR(pho_superClusterEta[itrPho],pho_superClusterPhi[itrPho],(*allTrackEta)[itTrack],(*allTrackPhi)[itTrack]);
		float trk_vtxT=(*allTrackT)[itTrack]-(*allTrackdT)[itTrack];
		dT=pvAllT[pvCand]-trk_vtxT;
		if(dR<0.4)h_trackDz_pho->Fill((*allTrackDz)[itTrack]); 
		if(dR<0.4 && fabs((*allTrackDz)[itTrack])<0.1/* && fabs((*allTrackDxy)[itTrack])<0.1 /*&& fabs((*allTrackX)[itTrack]-pvX)<0.1*/){
	    
		  //if(dR<0.4 && ((*allTrackZ)[itTrack]==pvZ) && ((*allTrackY)[itTrack]==pvY) && ((*allTrackX)[itTrack]==pvX)){
		  SumChargedPt+=(*allTrackPt)[itTrack];
		  h_trackX_pho->Fill((*allTrackX)[itTrack]-pvX);
		  h_trackY_pho->Fill((*allTrackY)[itTrack]-pvY);
		  h_trackZ_pho->Fill((*allTrackZ)[itTrack]-pvZ);
		  //h_trackDz_pho->Fill((*allTrackDz)[itTrack]);
		  h_trackDxy_pho->Fill((*allTrackDxy)[itTrack]);
		  h_dT_pho->Fill(dT);
		}
		if(dR<0.4 && fabs(dT)<0.090 && fabs((*allTrackDz)[itTrack])<0.5 /*&& fabs((*allTrackY)[itTrack]-pvY)<0.1 && fabs((*allTrackX)[itTrack]-pvX)<0.1*/)t_SumChargedPt+=(*allTrackPt)[itTrack];
	      }
	      //cout<<"Loop 2 done\n";
	      myIso=SumChargedPt/phoPt[itrPho];
	      myIso_dT=t_SumChargedPt/phoPt[itrPho];
	      h_ChargedFractionR04_pho->Fill(myIso);
	      h_ChargedFractionR04_dT_pho->Fill(myIso_dT);
	      h_diff_pho->Fill((pho_pfIsoChargedHadronIso[itrPho]-SumChargedPt));
	      h_CMSchargedIso04_pho->Fill(pho_pfIsoChargedHadronIso[itrPho]/phoPt[itrPho]);



	      Int_t itGenPho=-99;
	      pt=-1.;
	      for(int i=0;i<nGenParticle;i++){
		if(fabs(gParticleId[i])==22 && gParticleStatus[i]==1){
		  pt=gParticlePt[i];
		  if(fabs(phoPt[itrPho]-pt)<2. && fabs(gParticleEta[i]-phoEta[itrPho])<0.001 && fabs(gParticlePhi[i]-phoPhi[itrPho])<0.001){
		    itGenPho=i;
		    h_genPho->Fill(gParticlePt[itGenPho]);
		    h_genPho_vtx->Fill(nPV); 
		    if(myIso_dT<0.15 && myIso_dT!=-9999){
		      h_recoPho_dT->Fill(gParticlePt[itGenPho]);
		      h_recoPho_dTvtx->Fill(nPV);
		    }
		    if(myIso<0.15 && myIso!=-9999){
		      h_recoPho->Fill(gParticlePt[itGenPho]);
		      h_recoPho_vtx->Fill(nPV);
		    }
		    break;
		  }
		}
	      }



	    

	    //---------------------------------------------------------2nd leading photon---------------------------------------------------------
	    SumChargedPt=0.0;
	    t_SumChargedPt=0.0;
	    for(int itTrack= 0;itTrack<allTrackEta->size();itTrack++){
	      dR=DeltaR(pho_superClusterEta[itrPho2],pho_superClusterPhi[itrPho2],(*allTrackEta)[itTrack],(*allTrackPhi)[itTrack]);
	      float trk_vtxT2=(*allTrackT)[itTrack]-(*allTrackdT)[itTrack];
	      dT=pvAllT[pvCand]-trk_vtxT2;
	      if(dR<0.4)h_trackDz_pho->Fill((*allTrackDz)[itTrack]); 
	      if(dR<0.4 && fabs((*allTrackDz)[itTrack])<0.1/* && fabs((*allTrackDxy)[itTrack])<0.1 /*&& fabs((*allTrackX)[itTrack]-pvX)<0.1*/){
	    
		//if(dR<0.4 && ((*allTrackZ)[itTrack]==pvZ) && ((*allTrackY)[itTrack]==pvY) && ((*allTrackX)[itTrack]==pvX)){
		SumChargedPt+=(*allTrackPt)[itTrack];
		h_trackX_pho->Fill((*allTrackX)[itTrack]-pvX);
		h_trackY_pho->Fill((*allTrackY)[itTrack]-pvY);
		h_trackZ_pho->Fill((*allTrackZ)[itTrack]-pvZ);
		//h_trackDz_pho->Fill((*allTrackDz)[itTrack]);
		h_trackDxy_pho->Fill((*allTrackDxy)[itTrack]);
		h_dT_pho->Fill(dT);
	      }
	      if(dR<0.4 && fabs(dT)<0.090 && fabs((*allTrackDz)[itTrack])<0.5 /*&& fabs((*allTrackY)[itTrack]-pvY)<0.1 && fabs((*allTrackX)[itTrack]-pvX)<0.1*/)t_SumChargedPt+=(*allTrackPt)[itTrack];
	    }
	    //cout<<"Loop 2 done\n";
	    myIso=SumChargedPt/phoPt[itrPho2];
	    myIso_dT=t_SumChargedPt/phoPt[itrPho2];
	    h_ChargedFractionR04_pho->Fill(myIso);
	    h_ChargedFractionR04_dT_pho->Fill(myIso_dT);
	    h_diff_pho->Fill((pho_pfIsoChargedHadronIso[itrPho2]-SumChargedPt));
	    h_CMSchargedIso04_pho->Fill(pho_pfIsoChargedHadronIso[itrPho2]/phoPt[itrPho2]);


	    Int_t itGenPho2=-99;
	      pt=-1.;
	      for(int i=0;i<nGenParticle;i++){
		if(fabs(gParticleId[i])==22 && gParticleStatus[i]==1 && i!=itGenPho){
		  pt=gParticlePt[i];
		  if(fabs(phoPt[itrPho2]-pt)<2. && fabs(gParticleEta[i]-phoEta[itrPho2])<0.001 && fabs(gParticlePhi[i]-phoPhi[itrPho2])<0.001){
		    itGenPho2=i;
		    h_genPho->Fill(gParticlePt[itGenPho2]);
		    h_genPho_vtx->Fill(nPV); 
		    if(myIso_dT<0.15 && myIso_dT!=-9999){
		      h_recoPho_dT->Fill(gParticlePt[itGenPho2]);
		      h_recoPho_dTvtx->Fill(nPV);
		    }
		    if(myIso<0.15 && myIso!=-9999){
		      h_recoPho->Fill(gParticlePt[itGenPho2]);
		      h_recoPho_vtx->Fill(nPV);
		    }
		    break;
		  }
		}
	      }


	    }
	  

	    h_leadingPhoPt->Fill(phoPt[itrPho]);
	    h_leadingPhoTime->Fill(pho_t1);
	    h_nleadingPhoPt->Fill(phoPt[itrPho2]);
	    h_nleadingPhoTime->Fill(pho_t2);
	  }
	}  

      }
   }
   
      //      cout<<"=====================================================\n";
      
   //h_recoMuon_dT->Sumw2();
   //h_recoMuon->Sumw2();
   //h_genMuon->Sumw2();
   oFile->Write();
   oFile->Close();
}

void UpgradeAna::BookHistograms(){
  
  h_dR=new TH1D("deltaR","deltaR",100,0.,0.1);
  h_delPt=new TH1D("deltaPt","deltaPt",100,0.,100.);
  h_dT=new TH1D("dT","dT",1000,2.,-2.);
  h_leadingMuonPt=new TH1D("mu1Pt","mu1Pt",50,0.,500.);
  h_MatchedTrackPt=new TH1D("Tr1Pt","Tr1Pt",50,0.,500.);
  h_leadingMuonTime=new TH1D("mu1T","mu1T",100,-2.,2.);
  h_ChargedFractionR04=new TH1D("ChargedFrac","ChargedFrac",50,0.,1.);
  h_ChargedFractionR04_dT=new TH1D("ChargedFrac_dT","ChargedFrac_dT",50,0.,1.);
  h_CMSchargedIso04=new TH1D("ChargedFracCMS","ChargedFracCMS",50,0.,1.);
  h_diff=new TH1D("h_diff","h_diff",100,100.,0.);
  h_trackX = new TH1D("trackX04","trackX04",100,-1.,1.);
  h_trackY = new TH1D("trackY04","trackY04",100,-1.,1.);
  h_trackZ = new TH1D("trackZ04","trackZ04",100,-1.,1.);
  h_trackDz = new TH1D("trackDz","trackDz",100,-1.,1.);
  h_trackDxy = new TH1D("trackDxy","trackDxy",100,-1.,1.);

  h_recoMuon_dT=new TH1D("recoMuon_dT","recoMuon_dT",50,0.,500.);
  h_recoMuon=new TH1D("recoMuon","recoMuon",50,0.,500.);
  h_genMuon=new TH1D("genMuon","genMuon",50,0.,500.);

  h_recoMuon_dTvtx=new TH1D("recoMuon_dTvtx","recoMuon_dTvtx",50,0.,500.);
  h_recoMuon_vtx=new TH1D("recoMuon_vtx","recoMuon_vtx",50,0.,500.);
  h_genMuon_vtx=new TH1D("genMuon_vtx","genMuon_vtx",50,0.,500.);


  //============================== photon histograms ================================================//
  t1vsPvZ=new TH2D("t1vsPVz","t1vsPVz",100,-15.0,15.0,100,-0.2,0.2);
  t2vsPvZ=new TH2D("t2vsPVz","t2vsPVz",100,-15.0,15.0,100,-0.2,0.2);
  h_dT_pho=new TH1D("dT_pho","dT_pho",1000,2.,-2.);
  h_leadingPhoPt=new TH1D("pho1Pt","pho1Pt",50,0.,500.);
  h_leadingPhoTime=new TH1D("pho1T","pho1T",100,3.,23.);

  h_nleadingPhoPt=new TH1D("pho2Pt","pho2Pt",50,0.,500.);
  h_nleadingPhoTime=new TH1D("pho2T","pho2T",100,3.,23.);
  h_ChargedFractionR04_pho=new TH1D("ChargedFrac_pho","ChargedFrac_pho",50,0.,1.);
  h_ChargedFractionR04_dT_pho=new TH1D("ChargedFrac_dT_pho","ChargedFrac_dT_pho",50,0.,1.);
  h_CMSchargedIso04_pho=new TH1D("ChargedFracCMS_pho","ChargedFracCMS_pho",50,0.,1.);
  h_diff_pho=new TH1D("h_diff_pho","h_diff_pho",100,100.,0.);
  h_trackX_pho = new TH1D("trackX04_pho","trackX04_pho",100,-1.,1.);
  h_trackY_pho = new TH1D("trackY04_pho","trackY04_pho",100,-1.,1.);
  h_trackZ_pho = new TH1D("trackZ04_pho","trackZ04_pho",100,-1.,1.);
  h_trackDz_pho = new TH1D("trackDz_pho","trackDz_pho",100,-1.,1.);
  h_trackDxy_pho = new TH1D("trackDxy_pho","trackDxy_pho",100,-1.,1.);

  h_diPhoPt = new TH1D("diphoPt","diphoPt",100,0.,1000.);
  h_diPhoM = new TH1D("diphoM","diphoM",50,0.,200.);
  h_diPhoEta = new TH1D("diphoEta","diphoEta",30,-3.,3.);
  h_diPhoPhi = new TH1D("diphoPhi","diphoPhi",30,-4.,4.);
  h_PhoEtaDiff = new TH1D("phoEtaDiff","phoEtaDiff",30,-3.,3.);
  h_PhoPhiDiff = new TH1D("phoPhiDiff","phoPhiDiff",30,-4.,4.);
  h_diPhoPVIndex = new TH1D("diPhoPVIndex","diPhoPVIndex",200,0,200);
  h_genVsRecodiff_etal0p8 = new TH1D("genVsRecoZdiff_etal0.8","genVsRecoZdiff_etal0.8",60,-30.,30.);
  h_genVsRecodiff_etag0p8 = new TH1D("genVsRecoZdiff_etag0.8","genVsRecoZdiff_etag0.8",60,-30.,30.);
  h_genVspvdiff_etal0p8 = new TH1D("genVspvZdiff_etal0.8","genVspvZdiff_etal0.8",60,-30.,30.);
  h_genVspvdiff_etag0p8 = new TH1D("genVspvZdiff_etag0.8","genVspvZdiff_etag0.8",60,-30.,30.);
  h_tMinVsDeta=new TH2D("tMinVsDeta","tMinvsDeta",15,0.,3.,50,0.,0.1);
  h_diPhoMwithVtx = new TH1D("diphoMwithVtx","diphoMwithVtx",50,0.,200.);
  
  h_recoPho_dT=new TH1D("recoPho_dT","recoPho_dT",50,0.,500.);
  h_recoPho=new TH1D("recoPho","recoPho",50,0.,500.);
  h_genPho=new TH1D("genPho","genPho",50,0.,500.);

  h_recoPho_dTvtx=new TH1D("recoPho_dTvtx","recoPho_dTvtx",50,0.,500.);
  h_recoPho_vtx=new TH1D("recoPho_vtx","recoPho_vtx",50,0.,500.);
  h_genPho_vtx=new TH1D("genPho_vtx","genPho_vtx",50,0.,500.);
  
}

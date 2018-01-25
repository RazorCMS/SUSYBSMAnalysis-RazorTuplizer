#define TPrimeAna_cxx
#include "TPrimeAna.h"
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
    cerr << "Please give 4 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" << "data type"<<endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];
  const char *isData        = argv[4];
  TPrimeAna Tprime(inputFileList, outFileName, data, isData);
  cout << "dataset " << data << " " << endl;
  Tprime.EventLoop(data, isData);

  return 0;
}
void TPrimeAna::EventLoop(const char *data,const char *isData)
{

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    /*if(t_signalType==7){           
      for(unsigned int i=0; i<t_genPartP4->size(); i++) {
	cout<<i<< " id "<<(*t_genPartID)[i]<<" status "<< (*t_genPartStatus)[i]<<" dau1 "<<(*t_genPartDau1ID)[i]<< " dau2 "<<(*t_genPartDau2ID)[i]<<" mom 1 "<<(*t_genPartMom1ID)[i]<<" mom 2 "<<(*t_genPartMom2ID)[i]<<endl;

      }
      cout<<"================================================\n";
      }*/
    int b1=-99,b2=-99,W1=-99,W2=-99,t1=-99,t2=-99,jetindex=-99;
    float eta1,eta2,phi1,phi2,dr=99,dR,dr1;
    double dphi;
    
    /*=============================================================================================
      Gen info begins
      ==============================================================================================*/
    int Zflag=0, Vflag=0, veto=1,found1=0,found2=0,doubleZ=0,search=-99,iZ=-99,itop=-99,Mom=0,confirm=0,iZ2=-99,Mom2=0,qq=0,ll=0,sig=0;
    double evtwt=t_evtwt;
    double evtwt_Up=t_evtwt  ;
    double evtwt_Down=t_evtwt  ;
    char checkData[]="T";
    if(strcmp(checkData,isData)!=0){ //strcmp is 0 when strings are equal!!
      evtwt*=pow(0.93,t_jetTopJetP4->size()); //top tagging SF, +/- 0.09 to be applied later in datacard
      evtwt_Up*=evtwt;
      evtwt_Down*=evtwt;
      //==================== apply btag SFs==============================//                                                            
	 
      evtwt_Up *= t_btagsf_bcUp;
      evtwt_Up *= t_btagsf_lUp;

      evtwt_Down *= t_btagsf_bcDown;
      evtwt_Down *= t_btagsf_lDown;
	       
      //cout<<evtwt<<" "<<evtwt_Up<<" "<<evtwt_Down<<endl;
    }  
    
    /*=============================================================================================
      Leptonic Event Selection begins
      ==========================================================================================*/
       
    if(t_signalType == 9){//signalType = 7;=>tZtZ 
      //cout<<"ok\n";
      h_rawSize->Fill(t_jetAK4P4->size());
      int lep1=0; int lep2=0;int temp=0;
      int index[50],ctr;
      vector <int> elid,muid,goodElid,goodMuid;
      ctr=0;
      double sqmu,sqe,mass,s1,s2,dPhi,st,stAK8;
      double px,py,pt;
      sqmu=0.105658*0.105658;
      sqe=0.511*0.511*0.000001;
      double htAK8,htAK4;
      TLorentzVector v1(-1000,-1000,-1000,-1000);
      h_Cutflow->Fill(0.);
      h_Cutflow_DvsMC->Fill(0.);
      h_Ele_Cutflow_DvsMC->Fill(0.);
      //cout<<"=============================================\n";
      //cout<<t_muonP4->size()<<" "<<t_elecP4->size()<<endl;
     
      if(t_muonP4->size()>=2){
	for(unsigned int i=0; i<t_muonP4->size()-1; i++) {
	  for(unsigned int j=i+1; j<t_muonP4->size(); j++) {
	    //if((*t_muCharge)[i]*(*t_muCharge)[j]<0){
	    //cout<<((*t_muonP4)[i]+(*t_muonP4)[j]).Pt()<<" "<<(*t_ZllP4)[0].Pt()<<endl;
	      if(fabs(((*t_muonP4)[i]+(*t_muonP4)[j]).Pt()-(*t_ZllP4)[0].Pt())< 0.001){
		goodMuid.push_back(i);
		goodMuid.push_back(j);
		break;
	      }
	      //}
	  }
	}
      }
      //cout<<goodMuid.size()<<endl;

      if(t_elecP4->size()>=2){
	for(unsigned int i=0; i<t_elecP4->size()-1; i++) {
	  for(unsigned int j=i+1; j<t_elecP4->size(); j++) {
	    //if((*t_elCharge)[i]*(*t_elCharge)[j]){
	      if(fabs(((*t_elecP4)[i]+(*t_elecP4)[j]).Pt()-(*t_ZllP4)[0].Pt())<0.001){
		goodElid.push_back(i);
		goodElid.push_back(j);
		break;
	      }
	      //}
	  }
	}
	}
 
      int e1=99;int e2=99;
      int mu1=99;int mu2=99;
      //===============Electron/Muon reco efficiency=============================================================//
      for(unsigned int j=0; j<t_genPartP4->size(); j++) {
	if(fabs((*t_genPartID)[j])==11 && ((*t_genPartMom1ID)[j]==23 || (*t_genPartMom2ID)[j]==23)){
	//cout<<(*t_genPartP4)[j].Pt()<<endl;
	if(e1==99)e1=j;
	else e2=j;
      }
      else if(fabs((*t_genPartID)[j])==13 && ((*t_genPartMom1ID)[j]==23 || (*t_genPartMom2ID)[j]==23)){
	if(mu1==99)mu1=j;
	else mu2=j;
      }
    }
    if(e1!=99 && e2!=99){
      if((*t_genPartP4)[e1].Pt()>(*t_genPartP4)[e2].Pt()){
	h_genEle1Pt->Fill((*t_genPartP4)[e1].Pt());
	h_genEle2Pt->Fill((*t_genPartP4)[e2].Pt());
	LepReco(e1,e2,0,elid);
      }
      else{
	h_genEle2Pt->Fill((*t_genPartP4)[e1].Pt());
	h_genEle1Pt->Fill((*t_genPartP4)[e2].Pt());
	LepReco(e2,e1,0,elid);
      }
    }
    if(mu1!=99 && mu2!=99){
      if((*t_genPartP4)[mu1].Pt()>(*t_genPartP4)[mu2].Pt()){
	h_genMu1Pt->Fill((*t_genPartP4)[mu1].Pt());
	h_genMu2Pt->Fill((*t_genPartP4)[mu2].Pt());
	LepReco(mu1,mu2,1,muid);
      }
      else{
	h_genMu2Pt->Fill((*t_genPartP4)[mu1].Pt());
	h_genMu1Pt->Fill((*t_genPartP4)[mu2].Pt());
	LepReco(mu2,mu1,1,muid);
      }
    }

    // vector <int> cleanJetAK4=selGoodAK4Jet(elid,muid);//remove any AK4 jet having dr<0.4 with a lepton
    //htAK4=HTCleanAK4(cleanJetAK4);
    vector <int>cleanJetAK4;
    for(unsigned int j=0; j<t_jetAK4P4->size(); j++) {
      cleanJetAK4.push_back(j);
    }
   
    if(goodElid.size()==2){
      h_Cutflow->Fill(11.);
      vector<int>::iterator itEle = goodElid.begin();
      eta1=(*t_elecP4)[*itEle].Eta();
      phi1=(*t_elecP4)[*itEle].Phi();
      eta2=(*t_elecP4)[*(itEle+1)].Eta();
      phi2=(*t_elecP4)[*(itEle+1)].Phi();
	 
      v1= (*t_elecP4)[*itEle]+(*t_elecP4)[*(itEle+1)];
      dR= DeltaR(eta1,phi1,eta2,phi2);
      dPhi=fabs(DeltaPhi(phi1,phi2));
	 
      //=============calculate dielectron invariant mass====================//
      s1=sqrt(((*t_elecP4)[*itEle].P()*(*t_elecP4)[*itEle].P()+sqe)*((*t_elecP4)[*(itEle+1)].P()*(*t_elecP4)[*(itEle+1)].P()+sqe));
      s2=(*t_elecP4)[*itEle].Px()*(*t_elecP4)[*(itEle+1)].Px()+(*t_elecP4)[*itEle].Py()*(*t_elecP4)[*(itEle+1)].Py()+(*t_elecP4)[*itEle].Pz()*(*t_elecP4)[*(itEle+1)].Pz();
      mass=sqrt(2.0*(sqe+(s1-s2)));
	 
      st=t_HT+(*t_elecP4)[*itEle].Pt()+(*t_elecP4)[*(itEle+1)].Pt();
	 
	 
      h_ele_leadingLepPt_Zcut->Fill((*t_elecP4)[*itEle].Pt());
      h_ele_leadingLep2Pt_Zcut->Fill((*t_elecP4)[*(itEle+1)].Pt());
      h_ele_leadingLepPhi_Zcut->Fill((*t_elecP4)[*itEle].Phi());
      h_ele_leadingLep2Phi_Zcut->Fill((*t_elecP4)[*(itEle+1)].Phi());
      h_ele_leadingLepEta_Zcut->Fill((*t_elecP4)[*itEle].Eta());
      h_ele_leadingLep2Eta_Zcut->Fill((*t_elecP4)[*(itEle+1)].Eta());
      h_ele_dR_Zcut->Fill(dR);
      h_ele_dPhi_Zcut->Fill(dPhi);
      h_ele_mass_Zcut->Fill(mass);
      h_ele_st_Zcut->Fill(st);
      h_ele_dRvsdPhiZcut->Fill(dR,dPhi);
      px=(*t_elecP4)[*itEle].Px()+(*t_elecP4)[*(itEle+1)].Px();
      py=(*t_elecP4)[*itEle].Py()+(*t_elecP4)[*(itEle+1)].Py();
      pt=sqrt(px*px+py*py);
      h_dielept->Fill(pt);
      if(dPhi<2.5){
	h_ele_leadingLepPt_dPhicut->Fill((*t_elecP4)[*itEle].Pt());
	h_ele_leadingLep2Pt_dPhicut->Fill((*t_elecP4)[*(itEle+1)].Pt());
	h_ele_leadingLepPhi_dPhicut->Fill((*t_elecP4)[*itEle].Phi());
	h_ele_leadingLep2Phi_dPhicut->Fill((*t_elecP4)[*(itEle+1)].Phi());
	h_ele_leadingLepEta_dPhicut->Fill((*t_elecP4)[*itEle].Eta());
	h_ele_leadingLep2Eta_dPhicut->Fill((*t_elecP4)[*(itEle+1)].Eta());
	h_ele_dR_dPhicut->Fill(dR);
	h_ele_dPhi_dPhicut->Fill(dPhi);
	h_ele_mass_dPhicut->Fill(mass);
	h_ele_st_dPhicut->Fill(st);
      }
	 
      if(dR>0.2){
	h_ele_leadingLepPt_dRcut->Fill((*t_elecP4)[*itEle].Pt());
	h_ele_leadingLep2Pt_dRcut->Fill((*t_elecP4)[*(itEle+1)].Pt());
	h_ele_leadingLepPhi_dRcut->Fill((*t_elecP4)[*itEle].Phi());
	h_ele_leadingLep2Phi_dRcut->Fill((*t_elecP4)[*(itEle+1)].Phi());
	h_ele_leadingLepEta_dRcut->Fill((*t_elecP4)[*itEle].Eta());
	h_ele_leadingLep2Eta_dRcut->Fill((*t_elecP4)[*(itEle+1)].Eta());
	h_ele_dR_dRcut->Fill(dR);
	h_ele_dPhi_dRcut->Fill(dPhi);
	h_ele_mass_dRcut->Fill(mass);
	h_ele_st_dRcut->Fill(st);
      }
    }
       
       
    if(goodMuid.size()==2){
      vector<int>::iterator itMu = goodMuid.begin();
      h_Cutflow->Fill(12.);
      eta1=(*t_muonP4)[*itMu].Eta();
      phi1=(*t_muonP4)[*itMu].Phi();
      eta2=(*t_muonP4)[*(itMu+1)].Eta();
      phi2=(*t_muonP4)[*(itMu+1)].Phi();
	 
      v1= (*t_muonP4)[*itMu]+(*t_muonP4)[*(itMu+1)];
      dR= DeltaR(eta1,phi1,eta2,phi2);
      dPhi=fabs(DeltaPhi(phi1,phi2));
	 
      //=============calculate dimuon invariant mass====================//
      s1=sqrt(((*t_muonP4)[*itMu].P()*(*t_muonP4)[*itMu].P()+sqmu)*((*t_muonP4)[*(itMu+1)].P()*(*t_muonP4)[*(itMu+1)].P()+sqmu));
      s2=(*t_muonP4)[*itMu].Px()*(*t_muonP4)[*(itMu+1)].Px()+(*t_muonP4)[*itMu].Py()*(*t_muonP4)[*(itMu+1)].Py()+(*t_muonP4)[*itMu].Pz()*(*t_muonP4)[*(itMu+1)].Pz();
      mass=sqrt(2.0*(sqmu+(s1-s2)));
	 
	 
      st=t_HT+(*t_muonP4)[*itMu].Pt()+(*t_muonP4)[*(itMu+1)].Pt();
	 
	 
      h_mu_leadingLepPt_Zcut->Fill((*t_muonP4)[*itMu].Pt());
      h_mu_leadingLep2Pt_Zcut->Fill((*t_muonP4)[*(itMu+1)].Pt());
      h_mu_leadingLepPhi_Zcut->Fill((*t_muonP4)[*itMu].Phi());
      h_mu_leadingLep2Phi_Zcut->Fill((*t_muonP4)[*(itMu+1)].Phi());
      h_mu_leadingLepEta_Zcut->Fill((*t_muonP4)[*itMu].Eta());
      h_mu_leadingLep2Eta_Zcut->Fill((*t_muonP4)[*(itMu+1)].Eta());
      h_mu_dR_Zcut->Fill(dR);
      h_mu_dPhi_Zcut->Fill(dPhi);
      h_mu_mass_Zcut->Fill(mass);
      h_mu_st_Zcut->Fill(st);
      h_mu_dRvsdPhiZcut->Fill(dR,dPhi);
      px=(*t_muonP4)[*itMu].Px()+(*t_muonP4)[*(itMu+1)].Px();
      py=(*t_muonP4)[*itMu].Py()+(*t_muonP4)[*(itMu+1)].Py();
      pt=sqrt(px*px+py*py);
      h_dimupt->Fill(pt);
      if(dPhi<2.5){
	h_mu_leadingLepPt_dPhicut->Fill((*t_muonP4)[*itMu].Pt());
	h_mu_leadingLep2Pt_dPhicut->Fill((*t_muonP4)[*(itMu+1)].Pt());
	h_mu_leadingLepPhi_dPhicut->Fill((*t_muonP4)[*itMu].Phi());
	h_mu_leadingLep2Phi_dPhicut->Fill((*t_muonP4)[*(itMu+1)].Phi());
	h_mu_leadingLepEta_dPhicut->Fill((*t_muonP4)[*itMu].Eta());
	h_mu_leadingLep2Eta_dPhicut->Fill((*t_muonP4)[*(itMu+1)].Eta());
	h_mu_dR_dPhicut->Fill(dR);
	h_mu_dPhi_dPhicut->Fill(dPhi);
	h_mu_mass_dPhicut->Fill(mass);
	h_mu_st_dPhicut->Fill(st);
      }
	 
      if(dR>0.2){
	h_mu_leadingLepPt_dRcut->Fill((*t_muonP4)[*itMu].Pt());
	h_mu_leadingLep2Pt_dRcut->Fill((*t_muonP4)[*(itMu+1)].Pt());
	h_mu_leadingLepPhi_dRcut->Fill((*t_muonP4)[*itMu].Phi());
	h_mu_leadingLep2Phi_dRcut->Fill((*t_muonP4)[*(itMu+1)].Phi());
	h_mu_leadingLepEta_dRcut->Fill((*t_muonP4)[*itMu].Eta());
	h_mu_leadingLep2Eta_dRcut->Fill((*t_muonP4)[*(itMu+1)].Eta());
	h_mu_dR_dRcut->Fill(dR);
	h_mu_dPhi_dRcut->Fill(dPhi);
	h_mu_mass_dRcut->Fill(mass);
	h_mu_st_dRcut->Fill(st);
      }
    }
    if(goodMuid.size()>=2 || goodElid.size()>=2)
      h_Nlep->Fill(goodMuid.size()+goodElid.size()-2);
    /*======================================================================================================
      Now select only those events with at least one Z decaying to 2 mu/ele
      ========================================================================================================*/
    //sig is incorporated because I'm looking at tZtZevent.
    int Zcheck=0;int Zidx=99;
    if((goodMuid.size()==2 || goodElid.size()==2) && (goodMuid.size()!= goodElid.size())){ //So that only one Z decays leptonically.
	 
      h_Cutflow->Fill(2.0);
      h_Cutflow_DvsMC->Fill(2.);
      h_Ele_Cutflow_DvsMC->Fill(2.);
      if(goodMuid.size()==2)h_Cutflow->Fill(13.0);
      else h_Cutflow->Fill(14.0);
      for(unsigned int j=0; j<t_metP4->size(); j++){
	h_met->Fill((*t_metP4)[j].Pt());
	phi1=v1.Phi();
	phi2=(*t_metP4)[j].Phi();
	dPhi=fabs(DeltaPhi(phi1,phi2));
	h_metDPhi->Fill(dPhi);
      }
      //for matching reconstructed Z with gen Z
      dr=99;
      for(unsigned int i=0; i<t_genPartP4->size(); i++) {
	if((*t_genPartID)[i]==23){
	  phi1=v1.Phi();
	  phi2=(*t_genPartP4)[i].Phi();
	  eta1=v1.Eta();
	  eta2=(*t_genPartP4)[i].Eta();
	  dR=DeltaR(eta1,phi1,eta2,phi2);
	  h_dRZ->Fill(dR);
	  if(dR<0.3 && dR<dr){
	    Zcheck=1;
	    Zidx=i;
	    dr=dR;
	  }
	}
      }
      vector <int>cleanJetAK8;
      for(unsigned int j=0; j<t_jetAK8P4->size(); j++) {
	cleanJetAK8.push_back(j);
      }
      
      //================= Checking Reco Top===============================================//
      int NrecoTop=t_jetTopJetP4->size();
      /*vector<int> Topid;
      for (vector<int>::iterator it = cleanJetAK8.begin() ; it != cleanJetAK8.end(); ++it){	
	if( (*t_jetAK8_SoftDropMass)[*it]>105. && (*t_jetAK8_SoftDropMass)[*it]<220. && fabs((*t_jetAK8P4)[*it].Eta())<2.4 && (*t_jetAK8_tau3)[*it]/(*t_jetAK8_tau2)[*it]<0.67 &&  (*t_jetAK8P4)[*it].Pt()>400.0){
	   
	  h_recoTopPt->Fill((*t_jetAK8P4)[*it].Pt());
	  NrecoTop++;
	  Topid.push_back(*it);
	}
	}*/
      h_NrecoTop->Fill(NrecoTop);
         
      //cout<<NrecoTop<<" "<<t_jetTopJetPt->size()<<endl;
      //================= Checking Reco Higgs===============================================//
      int NrecoH=0;
      vector<int> Hid;
      for (vector<int>::iterator it = cleanJetAK8.begin() ; it != cleanJetAK8.end(); ++it){
	
	if((*t_jetAK8_MassPruned)[*it]>105. && (*t_jetAK8_MassPruned)[*it]<135. && fabs((*t_jetAK8P4)[*it].Eta())<2.4  && (*t_jetAK8_tau2)[*it]/(*t_jetAK8_tau1)[*it]<0.6){

	  h_recoHPt->Fill((*t_jetAK8P4)[*it].Pt());
	  NrecoH++;
	  Hid.push_back(*it);
	}
      }
      h_NrecoH->Fill(NrecoH);
      //================== reco b==========================================//
	 

      vector<int> cleanBJet;
      for (vector<int>::iterator it = cleanJetAK4.begin() ; it != cleanJetAK4.end(); ++it){
	if((*t_jetAK4CSV)[*it]>0.8484 && fabs((*t_jetAK4P4)[*it].Eta())<2.4)cleanBJet.push_back(*it);
      }
      h_Nrecob->Fill(cleanBJet.size());
      vector<int>::iterator itB=cleanBJet.begin();
      if(cleanBJet.size()!=0) h_LeadingrecobPt->Fill((*t_jetAK4P4)[*itB].Pt());

      vector<int> cleanBJetT;
      vector<int> cleanBJetL;
      for (vector<int>::iterator it = cleanJetAK4.begin() ; it != cleanJetAK4.end(); ++it){
	if((*t_jetAK4CSV)[*it]>0.9535 && fabs((*t_jetAK4P4)[*it].Eta())<2.4)cleanBJetT.push_back(*it);
	if((*t_jetAK4CSV)[*it]>0.5426 && fabs((*t_jetAK4P4)[*it].Eta())<2.4)cleanBJetL.push_back(*it);
      }

      //===================== reco AK4 jets which aren't b tagged =======================================//

      vector<int> cleanNonBJet;
      int flag=0;
      for (vector<int>::iterator it = cleanJetAK4.begin() ; it != cleanJetAK4.end(); ++it){
	for (vector<int>::iterator itb = cleanBJetL.begin() ; itb != cleanBJetL.end(); ++itb){
	  if(*it==*itb){flag=1;break;}
	}
	if(!flag)cleanNonBJet.push_back(*it);
      }
	
      //==================== reco V jets ===================================//
      int NrecoV=t_jetWJetP4->size();
      /*vector<int> Vid;
      for (vector<int>::iterator it = cleanJetAK8.begin() ; it != cleanJetAK8.end(); ++it){
	if((*t_jetAK8_MassPruned)[*it]>65. && (*t_jetAK8_MassPruned)[*it]<105. && fabs((*t_jetAK8P4)[*it].Eta())<2.4 && (*t_jetAK8_tau2)[*it]/(*t_jetAK8_tau1)[*it]<0.6 && (*t_jetAK8P4)[*it].Pt()>180.0 ){
	  NrecoV++;
	  Vid.push_back(*it);
	}
	}*/
      //================================= Jets ================================================//
      
      for (vector<int>::iterator it = cleanJetAK8.begin() ; it != cleanJetAK8.end(); ++it){
	h_JetsAK8Pt->Fill((*t_jetAK8P4)[*it].Pt());
	h_JetsAK8E->Fill((*t_jetAK8P4)[*it].E());
	h_JetsAK8Eta->Fill((*t_jetAK8P4)[*it].Eta());
	h_JetsAK8Phi->Fill((*t_jetAK8P4)[*it].Phi());
      }
      h_NJetsAK8->Fill(cleanJetAK8.size());
   	
      
      if(cleanJetAK8.size()!=0)h_LeadingAK8Pt->Fill((*t_jetAK8P4)[*cleanJetAK8.begin()].Pt());
      
    
      for (vector<int>::iterator it = cleanJetAK4.begin() ; it != cleanJetAK4.end(); ++it){
	h_JetsAK4Pt->Fill((*t_jetAK4P4)[*it].Pt());
	h_JetsAK4E->Fill((*t_jetAK4P4)[*it].E());
	h_JetsAK4Eta->Fill((*t_jetAK4P4)[*it].Eta());
	h_JetsAK4Phi->Fill((*t_jetAK4P4)[*it].Phi());
      }
      h_NJetsAK4->Fill(cleanJetAK4.size());
      //	cout<<(*t_jetAK4CSV)[*it]<<endl;
      
      
      vector<int>::iterator itJ = cleanJetAK4.begin();
      if(cleanJetAK4.size()!=0)h_LeadingAK4Pt->Fill((*t_jetAK4P4)[*itJ].Pt());
      int Njets;
      if(t_jetAK8P4->size()>=1){//at least 1 fat jet
	Njets=t_jetAK8P4->size()+cleanJetAK4.size();
	if(Njets>=4) h_4jetEvent->Fill(Njets);
      }
	 
      //htAK8= HTCleanAK8(cleanJetAK8);
      h_JetsAK8HT->Fill(htAK8);
     
      h_JetsAK4HT->Fill(t_HT); // ht of AK4 jets was already calculated using HTCleanAK4
      
      //cout<<"DataMc\n";
      //=======================================================================================//
      //--------------------------Plots for Data vs MC----------------------------------------//
      //=======================================================================================//
      vector <float> chi_output;
      
      int check=1;
      int flmu1=1; int flmu2=1;
      
      vector<int>::iterator itMu = goodMuid.begin();
      if(goodMuid.size()!=0 && (*t_muonP4)[*itMu].Pt()>=30. && (*t_muonP4)[*(itMu+1)].Pt()>=20.){
	h_Cutflow_DvsMC->Fill(3.);
	if(t_HT>300.){
	  h_Cutflow_DvsMC->Fill(4.);
	  
	  if((*t_muonP4)[*itMu].Pt()>=500. && fabs((*t_muonP4)[*itMu].Eta())>1.8)flmu1=0; //remove endcap muons with high pt
	  if((*t_muonP4)[*(itMu+1)].Pt()>=500. && fabs((*t_muonP4)[*(itMu+1)].Eta())>1.8)flmu2=0;
	  
	  if(flmu1 && flmu2){

	    h_Cutflow_DvsMC->Fill(5.);
	    
	   
	    h_Cutflow_DvsMC->Fill(6.);
	    //cout<<evtwt<<endl;
	   
	    //cout<<evtwt<<" "<<evtwt_Up<<" "<<evtwt_Down<<endl;
	    //cout<<evtwt<<endl;
	    h_Cutflow_DvsMC->Fill(7.);
	    //cout<<"MU\n";
	    h_l1Pt->Fill((*t_muonP4)[*itMu].Pt(),evtwt);
	    h_l1Eta->Fill((*t_muonP4)[*itMu].Eta(),evtwt);
	    h_l1Phi->Fill((*t_muonP4)[*itMu].Phi(),evtwt);
	  
	    h_l1Pt_Up->Fill((*t_muonP4)[*itMu].Pt(),evtwt_Up);
	    h_l1Eta_Up->Fill((*t_muonP4)[*itMu].Eta(),evtwt_Up);
	    h_l1Phi_Up->Fill((*t_muonP4)[*itMu].Phi(),evtwt_Up);
	      
	    h_l1Pt_Down->Fill((*t_muonP4)[*itMu].Pt(),evtwt_Down);
	    h_l1Eta_Down->Fill((*t_muonP4)[*itMu].Eta(),evtwt_Down);
	    h_l1Phi_Down->Fill((*t_muonP4)[*itMu].Phi(),evtwt_Down);

	    h_l2Pt->Fill((*t_muonP4)[*(itMu+1)].Pt(),evtwt);
	    h_l2Eta->Fill((*t_muonP4)[*(itMu+1)].Eta(),evtwt);
	    h_l2Phi->Fill((*t_muonP4)[*(itMu+1)].Phi(),evtwt);

	    h_l2Pt_Up->Fill((*t_muonP4)[*(itMu+1)].Pt(),evtwt_Up);
	    h_l2Eta_Up->Fill((*t_muonP4)[*(itMu+1)].Eta(),evtwt_Up);
	    h_l2Phi_Up->Fill((*t_muonP4)[*(itMu+1)].Phi(),evtwt_Up);

	    h_l2Pt_Down->Fill((*t_muonP4)[*(itMu+1)].Pt(),evtwt_Down);
	    h_l2Eta_Down->Fill((*t_muonP4)[*(itMu+1)].Eta(),evtwt_Down);
	    h_l2Phi_Down->Fill((*t_muonP4)[*(itMu+1)].Phi(),evtwt_Down);

	    h_dilepMass->Fill(v1.M(),evtwt);
	    h_dilepMass_Up->Fill(v1.M(),evtwt_Up);
	    h_dilepMass_Down->Fill(v1.M(),evtwt_Down);

	    eta1=(*t_muonP4)[*itMu].Eta();
	    phi1=(*t_muonP4)[*itMu].Phi();
	    eta2=(*t_muonP4)[*(itMu+1)].Eta();
	    phi2=(*t_muonP4)[*(itMu+1)].Phi();
	    dR= DeltaR(eta1,phi1,eta2,phi2);
	    dPhi=fabs(DeltaPhi(phi1,phi2));

	    h_l1l2_dR->Fill(dR,evtwt);
	    h_l1l2_dPhi->Fill(dPhi,evtwt);
	    h_l1l2_dR_Up->Fill(dR,evtwt_Up);
	    h_l1l2_dPhi_Up->Fill(dPhi,evtwt_Up);
	    h_l1l2_dR_Down->Fill(dR,evtwt_Down);
	    h_l1l2_dPhi_Down->Fill(dPhi,evtwt_Down);

	    h_dilepPt->Fill(v1.Pt(),evtwt);
	    h_dilepEta->Fill(v1.Eta(),evtwt);
	    h_dilepPhi->Fill(v1.Phi(),evtwt);
	    h_ST->Fill(st,evtwt);
	    h_npv->Fill(ta_npv,evtwt);
	    h_dilepPt_Up->Fill(v1.Pt(),evtwt_Up);
	    h_dilepEta_Up->Fill(v1.Eta(),evtwt_Up);
	    h_dilepPhi_Up->Fill(v1.Phi(),evtwt_Up);
	    h_ST_Up->Fill(st,evtwt_Up);
	    h_npv_Up->Fill(ta_npv,evtwt_Up);
	    h_dilepPt_Down->Fill(v1.Pt(),evtwt_Down);
	    h_dilepEta_Down->Fill(v1.Eta(),evtwt_Down);
	    h_dilepPhi_Down->Fill(v1.Phi(),evtwt_Down);
	    h_ST_Down->Fill(st,evtwt_Down);
	    h_npv_Down->Fill(ta_npv,evtwt_Down);
	
	    h_NbjetL->Fill(cleanBJetL.size(),evtwt);
	    h_NbjetM->Fill(cleanBJet.size(),evtwt);
	    h_NAK4->Fill(cleanJetAK4.size(),evtwt);
	    h_NAK8->Fill(cleanJetAK8.size(),evtwt);
	    h_NbjetL_Up->Fill(cleanBJetL.size(),evtwt_Up);
	    h_NbjetM_Up->Fill(cleanBJet.size(),evtwt_Up);
	    h_NAK4_Up->Fill(cleanJetAK4.size(),evtwt_Up);
	    h_NAK8_Up->Fill(cleanJetAK8.size(),evtwt_Up);
	    h_NbjetL_Down->Fill(cleanBJetL.size(),evtwt_Down);
	    h_NbjetM_Down->Fill(cleanBJet.size(),evtwt_Down);
	    h_NAK4_Down->Fill(cleanJetAK4.size(),evtwt_Down);
	    h_NAK8_Down->Fill(cleanJetAK8.size(),evtwt_Down);
	    	      
	    for (vector<int>::iterator it = cleanJetAK4.begin() ; it != cleanJetAK4.end(); ++it){

	      h_AK4Pt->Fill((*t_jetAK4P4)[*it].Pt(),evtwt);
	      h_AK4Eta->Fill((*t_jetAK4P4)[*it].Eta(),evtwt);
	      h_AK4Phi->Fill((*t_jetAK4P4)[*it].Phi(),evtwt);

	      h_AK4Pt_Up->Fill((*t_jetAK4P4)[*it].Pt(),evtwt_Up);
	      h_AK4Eta_Up->Fill((*t_jetAK4P4)[*it].Eta(),evtwt_Up);
	      h_AK4Phi_Up->Fill((*t_jetAK4P4)[*it].Phi(),evtwt_Up);

	      h_AK4Pt_Down->Fill((*t_jetAK4P4)[*it].Pt(),evtwt_Down);
	      h_AK4Eta_Down->Fill((*t_jetAK4P4)[*it].Eta(),evtwt_Down);
	      h_AK4Phi_Down->Fill((*t_jetAK4P4)[*it].Phi(),evtwt_Down);
	    }
	    h_AK4HT->Fill(t_HT,evtwt);
	    h_AK4HT_Up->Fill(t_HT,evtwt_Up);
	    h_AK4HT_Down->Fill(t_HT,evtwt_Down);

	    if(cleanJetAK4.size()){
	      h_LeadAK4Pt->Fill((*t_jetAK4P4)[*cleanJetAK4.begin()].Pt(),evtwt);
	      h_LeadAK4Eta->Fill((*t_jetAK4P4)[*cleanJetAK4.begin()].Eta(),evtwt);
	      h_LeadAK4Phi->Fill((*t_jetAK4P4)[*cleanJetAK4.begin()].Phi(),evtwt);

	      h_LeadAK4Pt_Up->Fill((*t_jetAK4P4)[*cleanJetAK4.begin()].Pt(),evtwt_Up);
	      h_LeadAK4Eta_Up->Fill((*t_jetAK4P4)[*cleanJetAK4.begin()].Eta(),evtwt_Up);
	      h_LeadAK4Phi_Up->Fill((*t_jetAK4P4)[*cleanJetAK4.begin()].Phi(),evtwt_Up);

	      h_LeadAK4Pt_Down->Fill((*t_jetAK4P4)[*cleanJetAK4.begin()].Pt(),evtwt_Down);
	      h_LeadAK4Eta_Down->Fill((*t_jetAK4P4)[*cleanJetAK4.begin()].Eta(),evtwt_Down);
	      h_LeadAK4Phi_Down->Fill((*t_jetAK4P4)[*cleanJetAK4.begin()].Phi(),evtwt_Down);
	    }


	      
	    for (vector<int>::iterator it = cleanJetAK8.begin() ; it != cleanJetAK8.end(); ++it){
	      h_AK8Pt->Fill((*t_jetAK8P4)[*it].Pt(),evtwt);
	      h_AK8Eta->Fill((*t_jetAK8P4)[*it].Eta(),evtwt);
	      h_AK8Phi->Fill((*t_jetAK8P4)[*it].Phi(),evtwt);

	      h_AK8Pt_Up->Fill((*t_jetAK8P4)[*it].Pt(),evtwt_Up);
	      h_AK8Eta_Up->Fill((*t_jetAK8P4)[*it].Eta(),evtwt_Up);
	      h_AK8Phi_Up->Fill((*t_jetAK8P4)[*it].Phi(),evtwt_Up);

	      h_AK8Pt_Down->Fill((*t_jetAK8P4)[*it].Pt(),evtwt_Down);
	      h_AK8Eta_Down->Fill((*t_jetAK8P4)[*it].Eta(),evtwt_Down);
	      h_AK8Phi_Down->Fill((*t_jetAK8P4)[*it].Phi(),evtwt_Down);
	    }
	    h_AK8HT->Fill(htAK8,evtwt);
	    h_AK8HT_Up->Fill(htAK8,evtwt_Up);
	    h_AK8HT_Down->Fill(htAK8,evtwt_Down);
	    if(cleanJetAK8.size()){
	      h_LeadAK8Pt->Fill((*t_jetAK8P4)[*cleanJetAK8.begin()].Pt(),evtwt);
	      h_LeadAK8Eta->Fill((*t_jetAK8P4)[*cleanJetAK8.begin()].Eta(),evtwt);
	      h_LeadAK8Phi->Fill((*t_jetAK8P4)[*cleanJetAK8.begin()].Phi(),evtwt);

	      h_LeadAK8Pt_Up->Fill((*t_jetAK8P4)[*cleanJetAK8.begin()].Pt(),evtwt_Up);
	      h_LeadAK8Eta_Up->Fill((*t_jetAK8P4)[*cleanJetAK8.begin()].Eta(),evtwt_Up);
	      h_LeadAK8Phi_Up->Fill((*t_jetAK8P4)[*cleanJetAK8.begin()].Phi(),evtwt_Up);

	      h_LeadAK8Pt_Down->Fill((*t_jetAK8P4)[*cleanJetAK8.begin()].Pt(),evtwt_Down);
	      h_LeadAK8Eta_Down->Fill((*t_jetAK8P4)[*cleanJetAK8.begin()].Eta(),evtwt_Down);
	      h_LeadAK8Phi_Down->Fill((*t_jetAK8P4)[*cleanJetAK8.begin()].Phi(),evtwt_Down);
	    }
 
	    chi_output=GetChi(NrecoTop,NrecoV,v1,goodMuid,goodElid, cleanJetAK4, cleanBJetL);
	    //cout<<"Mu ok\n";
	    //cout<<chi_output[0];
	    if(chi_output[0]!=10000){
	      //cout<<"meow\n";
	      h_Mu_Chi2->Fill(chi_output[0],evtwt);
	      h_Mu_Chi2_Up->Fill(chi_output[0],evtwt_Up);
	      h_Mu_Chi2_Down->Fill(chi_output[0],evtwt_Down);
	      //cout<<"hmm\n";
	      if(NrecoV>=2 && NrecoTop>=1){
		h_Mu_Chi2_1t->Fill(chi_output[0],evtwt);
		h_Mu_Chi2_1t_Up->Fill(chi_output[0],evtwt_Up);
		h_Mu_Chi2_1t_Down->Fill(chi_output[0],evtwt_Down);

	      }
	      if(NrecoTop==1 && !NrecoV && cleanJetAK4.size()>=5){
		h_Mu_Chi2_1t0V->Fill(chi_output[0],evtwt);
		h_Mu_Chi2_1t0V_Up->Fill(chi_output[0],evtwt_Up);
		h_Mu_Chi2_1t0V_Down->Fill(chi_output[0],evtwt_Down);
	      }
	      if(!NrecoTop && NrecoV==1 && cleanJetAK4.size()>=6 ){
		h_Mu_Chi2_0t1V->Fill(chi_output[0],evtwt);
		h_Mu_Chi2_0t1V_Up->Fill(chi_output[0],evtwt_Up);
		h_Mu_Chi2_0t1V_Down->Fill(chi_output[0],evtwt_Down);
	      }
	      if(!NrecoTop && NrecoV>=2 && cleanJetAK4.size()>=4 ){
		h_Mu_Chi2_0t2V->Fill(chi_output[0],evtwt);
		h_Mu_Chi2_0t2V_Up->Fill(chi_output[0],evtwt_Up);
		h_Mu_Chi2_0t2V_Down->Fill(chi_output[0],evtwt_Down);
	      }
	      if(NrecoTop>=1 && NrecoV==1){
		h_Mu_Chi2_1t1V->Fill(chi_output[0],evtwt);
		h_Mu_Chi2_1t1V_Up->Fill(chi_output[0],evtwt_Up);
		h_Mu_Chi2_1t1V_Down->Fill(chi_output[0],evtwt_Down);
	      }
	      h_Mu_TMass_chi2_lep->Fill(chi_output[1],evtwt);
	      h_Mu_TMass_chi2_lep_Up->Fill(chi_output[1],evtwt_Up);
	      h_Mu_TMass_chi2_lep_Down->Fill(chi_output[1],evtwt_Down);
	      h_Mu_TMass_chi2_had->Fill(chi_output[2],evtwt);
	      h_Mu_TMass_chi2_had_Up->Fill(chi_output[2],evtwt_Up);
	      h_Mu_TMass_chi2_had_Down->Fill(chi_output[2],evtwt_Down);
	      h_Mu_TMass_chi2->Fill(chi_output[1],evtwt);
	      h_Mu_TMass_chi2_Up->Fill(chi_output[1],evtwt_Up);
	      h_Mu_TMass_chi2_Down->Fill(chi_output[1],evtwt_Down);
	      h_Mu_TMass_chi2->Fill(chi_output[2],evtwt);
	      h_Mu_TMass_chi2_Up->Fill(chi_output[2],evtwt_Up);
	      h_Mu_TMass_chi2_Down->Fill(chi_output[2],evtwt_Down);
	    }
	    if(chi_output[0]<20.){
	      h_Mu_Chi2_cut->Fill(chi_output[0],evtwt);
	      h_Mu_Chi2_cut_Up->Fill(chi_output[0],evtwt_Up);
	      h_Mu_Chi2_cut_Down->Fill(chi_output[0],evtwt_Down);
	      if(NrecoV>=2 && NrecoTop>=1){
		h_Mu_Chi2_1t_cut->Fill(chi_output[0],evtwt);
		h_Mu_Chi2_1t_cut_Up->Fill(chi_output[0],evtwt_Up);
		h_Mu_Chi2_1t_cut_Down->Fill(chi_output[0],evtwt_Down);
	      }
	      if(NrecoTop==1 && !NrecoV && cleanJetAK4.size()>=5){
		h_Mu_Chi2_1t0V_cut->Fill(chi_output[0],evtwt);
		h_Mu_Chi2_1t0V_cut_Up->Fill(chi_output[0],evtwt_Up);
		h_Mu_Chi2_1t0V_cut_Down->Fill(chi_output[0],evtwt_Down);

	      }
	      if(!NrecoTop && NrecoV==1 && cleanJetAK4.size()>=6 ){
		h_Mu_Chi2_0t1V_cut->Fill(chi_output[0],evtwt);
		h_Mu_Chi2_0t1V_cut_Up->Fill(chi_output[0],evtwt_Up);
		h_Mu_Chi2_0t1V_cut_Down->Fill(chi_output[0],evtwt_Down);
	      }
	      if(!NrecoTop && NrecoV>=2 && cleanJetAK4.size()>=4 ){
		h_Mu_Chi2_0t2V_cut->Fill(chi_output[0]);
		h_Mu_Chi2_0t2V_cut_Up->Fill(chi_output[0],evtwt_Up);
		h_Mu_Chi2_0t2V_cut_Down->Fill(chi_output[0],evtwt_Down);

	      }
	      if(NrecoTop>=1 && NrecoV==1){
		h_Mu_Chi2_1t1V_cut->Fill(chi_output[0]);
		h_Mu_Chi2_1t1V_cut_Up->Fill(chi_output[0],evtwt_Up);
		h_Mu_Chi2_1t1V_cut_Down->Fill(chi_output[0],evtwt_Down);
	      }
	      h_Mu_TMass_chi2_lep_cut->Fill(chi_output[1],evtwt);
	      h_Mu_TMass_chi2_lep_cut_Up->Fill(chi_output[1],evtwt_Up);
	      h_Mu_TMass_chi2_lep_cut_Down->Fill(chi_output[1],evtwt_Down);
	      h_Mu_TMass_chi2_had_cut->Fill(chi_output[2],evtwt);
	      h_Mu_TMass_chi2_had_cut_Up->Fill(chi_output[2],evtwt_Up);
	      h_Mu_TMass_chi2_had_cut_Down->Fill(chi_output[2],evtwt_Down);
	      h_Mu_TMass_chi2_cut->Fill(chi_output[1],evtwt);
	      h_Mu_TMass_chi2_cut_Up->Fill(chi_output[1],evtwt_Up);
	      h_Mu_TMass_chi2_cut_Down->Fill(chi_output[1],evtwt_Down);

	      h_Mu_TMass_chi2_cut->Fill(chi_output[2],evtwt);
	      h_Mu_TMass_chi2_cut_Up->Fill(chi_output[2],evtwt_Up);
	      h_Mu_TMass_chi2_cut_Down->Fill(chi_output[2],evtwt_Down);
	    }
	  }//if(flmu1 &&flmu2)
	}//if(t_HT>300.)
      }//if(goodMuid.size()!=0 && (*t_muonP4)[*itMu].Pt()>=30. && (*t_muonP4)[*(itMu+1)].Pt()>=20.)

    



      //x=x=x=x=x=x=x=x=x Electron data x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x//
      //cout<<check<<endl;
      vector<int>::iterator itEle = goodElid.begin();
      if(goodElid.size()!=0 && (*t_elecP4)[*itEle].Pt()>=30. && (*t_elecP4)[*(itEle+1)].Pt()>=20.){
	//if(!(*t_elecIsTight)[*itEle])cout<<(*t_elecIsTight)[*itEle]<<" "<<(*t_elecIsTight)[(*itEle+1)]<<endl;
	h_Ele_Cutflow_DvsMC->Fill(3.);
	if(t_HT>300.){
	  h_Ele_Cutflow_DvsMC->Fill(4.);
	  
	  //if((*t_elecP4)[*itEle].Pt()>=500. && fabs((*t_elecP4)[*itEle].Eta())>1.8)flmu1=0; //remove endcap muons with high pt
	  //if((*t_elecP4)[*(itEle+1)].Pt()>=500. && fabs((*t_elecP4)[*(itEle+1)].Eta())>1.8)flmu2=0;
	  
	  //if(flmu1 && flmu2){

	  h_Ele_Cutflow_DvsMC->Fill(5.);
	 
	 
	  h_Ele_Cutflow_DvsMC->Fill(6.);
	  //cout<<evtwt<<endl;
	   
	  h_Ele_Cutflow_DvsMC->Fill(7.);
	  //cout<<"El\n";
	  h_Ele_l1Pt->Fill((*t_elecP4)[*itEle].Pt(),evtwt);
	  h_Ele_l1Eta->Fill((*t_elecP4)[*itEle].Eta(),evtwt);
	  h_Ele_l1Phi->Fill((*t_elecP4)[*itEle].Phi(),evtwt);
	      
	  h_Ele_l1Pt_Down->Fill((*t_elecP4)[*itEle].Pt(),evtwt_Down);
	  h_Ele_l1Eta_Down->Fill((*t_elecP4)[*itEle].Eta(),evtwt_Down);
	  h_Ele_l1Phi_Down->Fill((*t_elecP4)[*itEle].Phi(),evtwt_Down);
	      
	  h_Ele_l1Pt_Up->Fill((*t_elecP4)[*itEle].Pt(),evtwt_Up);
	  h_Ele_l1Eta_Up->Fill((*t_elecP4)[*itEle].Eta(),evtwt_Up);
	  h_Ele_l1Phi_Up->Fill((*t_elecP4)[*itEle].Phi(),evtwt_Up);

	  h_Ele_l2Pt->Fill((*t_elecP4)[*(itEle+1)].Pt(),evtwt);
	  h_Ele_l2Eta->Fill((*t_elecP4)[*(itEle+1)].Eta(),evtwt);
	  h_Ele_l2Phi->Fill((*t_elecP4)[*(itEle+1)].Phi(),evtwt);
	      

	  h_Ele_l2Pt_Up->Fill((*t_elecP4)[*(itEle+1)].Pt(),evtwt_Up);
	  h_Ele_l2Eta_Up->Fill((*t_elecP4)[*(itEle+1)].Eta(),evtwt_Up);
	  h_Ele_l2Phi_Up->Fill((*t_elecP4)[*(itEle+1)].Phi(),evtwt_Up);

	  h_Ele_l2Pt_Down->Fill((*t_elecP4)[*(itEle+1)].Pt(),evtwt_Down);
	  h_Ele_l2Eta_Down->Fill((*t_elecP4)[*(itEle+1)].Eta(),evtwt_Down);
	  h_Ele_l2Phi_Down->Fill((*t_elecP4)[*(itEle+1)].Phi(),evtwt_Down);

	  h_Ele_dilepMass->Fill(v1.M(),evtwt);
	  h_Ele_dilepMass_Up->Fill(v1.M(),evtwt_Up);
	  h_Ele_dilepMass_Down->Fill(v1.M(),evtwt_Down);

	  eta1=(*t_elecP4)[*itEle].Eta();
	  phi1=(*t_elecP4)[*itEle].Phi();
	  eta2=(*t_elecP4)[*(itEle+1)].Eta();
	  phi2=(*t_elecP4)[*(itEle+1)].Phi();
	  dR= DeltaR(eta1,phi1,eta2,phi2);
	  dPhi=fabs(DeltaPhi(phi1,phi2));

	  h_Ele_l1l2_dR->Fill(dR,evtwt);
	  h_Ele_l1l2_dPhi->Fill(dPhi,evtwt);

	  h_Ele_l1l2_dR_Up->Fill(dR,evtwt_Up);
	  h_Ele_l1l2_dPhi_Up->Fill(dPhi,evtwt_Up);

	  h_Ele_l1l2_dR_Down->Fill(dR,evtwt_Down);
	  h_Ele_l1l2_dPhi_Down->Fill(dPhi,evtwt_Down);

	  h_Ele_dilepPt->Fill(v1.Pt(),evtwt);
	  h_Ele_dilepEta->Fill(v1.Eta(),evtwt);
	  h_Ele_dilepPhi->Fill(v1.Phi(),evtwt);
	  h_Ele_ST->Fill(st,evtwt);
	  h_Ele_npv->Fill(ta_npv,evtwt);
	
	  h_Ele_dilepPt_Up->Fill(v1.Pt(),evtwt_Up);
	  h_Ele_dilepEta_Up->Fill(v1.Eta(),evtwt_Up);
	  h_Ele_dilepPhi_Up->Fill(v1.Phi(),evtwt_Up);
	  h_Ele_ST_Up->Fill(st,evtwt_Up);
	  h_Ele_npv_Up->Fill(ta_npv,evtwt_Up);

	  h_Ele_dilepPt_Down->Fill(v1.Pt(),evtwt_Down);
	  h_Ele_dilepEta_Down->Fill(v1.Eta(),evtwt_Down);
	  h_Ele_dilepPhi_Down->Fill(v1.Phi(),evtwt_Down);
	  h_Ele_ST_Down->Fill(st,evtwt_Down);
	  h_Ele_npv_Down->Fill(ta_npv,evtwt_Down);

	  h_Ele_NbjetL->Fill(cleanBJetL.size(),evtwt);
	  h_Ele_NbjetM->Fill(cleanBJet.size(),evtwt);
	  h_Ele_NAK4->Fill(cleanJetAK4.size(),evtwt);
	  h_Ele_NAK8->Fill(cleanJetAK8.size(),evtwt);

	  h_Ele_NbjetL_Up->Fill(cleanBJetL.size(),evtwt_Up);
	  h_Ele_NbjetM_Up->Fill(cleanBJet.size(),evtwt_Up);
	  h_Ele_NAK4_Up->Fill(cleanJetAK4.size(),evtwt_Up);
	  h_Ele_NAK8_Up->Fill(cleanJetAK8.size(),evtwt_Up);

	  h_Ele_NbjetL_Down->Fill(cleanBJetL.size(),evtwt_Down);
	  h_Ele_NbjetM_Down->Fill(cleanBJet.size(),evtwt_Down);
	  h_Ele_NAK4_Down->Fill(cleanJetAK4.size(),evtwt_Down);
	  h_Ele_NAK8_Down->Fill(cleanJetAK8.size(),evtwt_Down);
	      
	  for (vector<int>::iterator it = cleanJetAK4.begin() ; it != cleanJetAK4.end(); ++it){
	    h_Ele_AK4Pt->Fill((*t_jetAK4P4)[*it].Pt(),evtwt);
	    h_Ele_AK4Eta->Fill((*t_jetAK4P4)[*it].Eta(),evtwt);
	    h_Ele_AK4Phi->Fill((*t_jetAK4P4)[*it].Phi(),evtwt);
		
	    h_Ele_AK4Pt_Up->Fill((*t_jetAK4P4)[*it].Pt(),evtwt_Up);
	    h_Ele_AK4Eta_Up->Fill((*t_jetAK4P4)[*it].Eta(),evtwt_Up);
	    h_Ele_AK4Phi_Up->Fill((*t_jetAK4P4)[*it].Phi(),evtwt_Up);

	    h_Ele_AK4Pt_Down->Fill((*t_jetAK4P4)[*it].Pt(),evtwt_Down);
	    h_Ele_AK4Eta_Down->Fill((*t_jetAK4P4)[*it].Eta(),evtwt_Down);
	    h_Ele_AK4Phi_Down->Fill((*t_jetAK4P4)[*it].Phi(),evtwt_Down);
	  }
	  h_Ele_AK4HT->Fill(t_HT,evtwt);
	  h_Ele_AK4HT_Up->Fill(t_HT,evtwt_Up);
	  h_Ele_AK4HT_Down->Fill(t_HT,evtwt_Down);
	  if(cleanJetAK4.size()){
	    h_Ele_LeadAK4Pt->Fill((*t_jetAK4P4)[*cleanJetAK4.begin()].Pt(),evtwt);
	    //if(fabs((*t_jetAK4P4)[*cleanJetAK4.begin()].Eta())>2.4)cout<<"oh\n";
	    h_Ele_LeadAK4Eta->Fill((*t_jetAK4P4)[*cleanJetAK4.begin()].Eta(),evtwt);
	    h_Ele_LeadAK4Phi->Fill((*t_jetAK4P4)[*cleanJetAK4.begin()].Phi(),evtwt);

	    h_Ele_LeadAK4Pt_Up->Fill((*t_jetAK4P4)[*cleanJetAK4.begin()].Pt(),evtwt_Up);
	    //if(fabs((*t_jetAK4P4)[*cleanJetAK4.begin()].Eta())>2.4)cout<<"oh\n";                                                        
	    h_Ele_LeadAK4Eta_Up->Fill((*t_jetAK4P4)[*cleanJetAK4.begin()].Eta(),evtwt_Up);
	    h_Ele_LeadAK4Phi_Up->Fill((*t_jetAK4P4)[*cleanJetAK4.begin()].Phi(),evtwt_Up);

	    h_Ele_LeadAK4Pt_Down->Fill((*t_jetAK4P4)[*cleanJetAK4.begin()].Pt(),evtwt_Down);
	    //if(fabs((*t_jetAK4P4)[*cleanJetAK4.begin()].Eta())>2.4)cout<<"oh\n";
	    h_Ele_LeadAK4Eta_Down->Fill((*t_jetAK4P4)[*cleanJetAK4.begin()].Eta(),evtwt_Down);
	    h_Ele_LeadAK4Phi_Down->Fill((*t_jetAK4P4)[*cleanJetAK4.begin()].Phi(),evtwt_Down);
	  }


	      
	  for (vector<int>::iterator it = cleanJetAK8.begin() ; it != cleanJetAK8.end(); ++it){
	    h_Ele_AK8Pt->Fill((*t_jetAK8P4)[*it].Pt(),evtwt);
	    h_Ele_AK8Eta->Fill((*t_jetAK8P4)[*it].Eta(),evtwt);
	    h_Ele_AK8Phi->Fill((*t_jetAK8P4)[*it].Phi(),evtwt);

	    h_Ele_AK8Pt_Up->Fill((*t_jetAK8P4)[*it].Pt(),evtwt_Up);
	    h_Ele_AK8Eta_Up->Fill((*t_jetAK8P4)[*it].Eta(),evtwt_Up);
	    h_Ele_AK8Phi_Up->Fill((*t_jetAK8P4)[*it].Phi(),evtwt_Up);

	    h_Ele_AK8Pt_Down->Fill((*t_jetAK8P4)[*it].Pt(),evtwt_Down);
	    h_Ele_AK8Eta_Down->Fill((*t_jetAK8P4)[*it].Eta(),evtwt_Down);
	    h_Ele_AK8Phi_Down->Fill((*t_jetAK8P4)[*it].Phi(),evtwt_Down);
	  }
	  h_Ele_AK8HT->Fill(htAK8,evtwt);
	  h_Ele_AK8HT_Up->Fill(htAK8,evtwt_Up);
	  h_Ele_AK8HT_Down->Fill(htAK8,evtwt_Down);
	  if(cleanJetAK8.size()){
	    h_Ele_LeadAK8Pt->Fill((*t_jetAK8P4)[*cleanJetAK8.begin()].Pt(),evtwt);
	    h_Ele_LeadAK8Eta->Fill((*t_jetAK8P4)[*cleanJetAK8.begin()].Eta(),evtwt);
	    h_Ele_LeadAK8Phi->Fill((*t_jetAK8P4)[*cleanJetAK8.begin()].Phi(),evtwt);

	    h_Ele_LeadAK8Pt_Up->Fill((*t_jetAK8P4)[*cleanJetAK8.begin()].Pt(),evtwt_Up);
	    h_Ele_LeadAK8Eta_Up->Fill((*t_jetAK8P4)[*cleanJetAK8.begin()].Eta(),evtwt_Up);
	    h_Ele_LeadAK8Phi_Up->Fill((*t_jetAK8P4)[*cleanJetAK8.begin()].Phi(),evtwt_Up);

	    h_Ele_LeadAK8Pt_Down->Fill((*t_jetAK8P4)[*cleanJetAK8.begin()].Pt(),evtwt_Down);
	    h_Ele_LeadAK8Eta_Down->Fill((*t_jetAK8P4)[*cleanJetAK8.begin()].Eta(),evtwt_Down);
	    h_Ele_LeadAK8Phi_Down->Fill((*t_jetAK8P4)[*cleanJetAK8.begin()].Phi(),evtwt_Down);
	  }

	  chi_output=GetChi(NrecoTop,NrecoV,v1,goodMuid,goodElid, cleanJetAK4, cleanBJetL);
	  //cout<<"El ok\n";
	  //cout<<chi_output[0];
	  if(chi_output[0]!=10000){
	    //cout<<"meow1\n";
	    h_Ele_Chi2->Fill(chi_output[0],evtwt);
	    h_Ele_Chi2_Up->Fill(chi_output[0],evtwt_Up);
	    h_Ele_Chi2_Down->Fill(chi_output[0],evtwt_Down);
	    //cout<<"ho\n";
	    if(NrecoV>=2 && NrecoTop>=1){
	      h_Ele_Chi2_1t->Fill(chi_output[0],evtwt);
	      h_Ele_Chi2_1t_Up->Fill(chi_output[0],evtwt_Up);
	      h_Ele_Chi2_1t_Down->Fill(chi_output[0],evtwt_Down);

	    }
	    if(NrecoTop==1 && !NrecoV && cleanJetAK4.size()>=5){
	      h_Ele_Chi2_1t0V->Fill(chi_output[0],evtwt);
	      h_Ele_Chi2_1t0V_Up->Fill(chi_output[0],evtwt_Up);
	      h_Ele_Chi2_1t0V_Down->Fill(chi_output[0],evtwt_Down);
	    }
	    if(!NrecoTop && NrecoV==1 && cleanJetAK4.size()>=6 ){
	      h_Ele_Chi2_0t1V->Fill(chi_output[0],evtwt);
	      h_Ele_Chi2_0t1V_Up->Fill(chi_output[0],evtwt_Up);
	      h_Ele_Chi2_0t1V_Down->Fill(chi_output[0],evtwt_Down);
	    }
	    if(!NrecoTop && NrecoV>=2 && cleanJetAK4.size()>=4 ){
	      h_Ele_Chi2_0t2V->Fill(chi_output[0],evtwt);
	      h_Ele_Chi2_0t2V_Up->Fill(chi_output[0],evtwt_Up);
	      h_Ele_Chi2_0t2V_Down->Fill(chi_output[0],evtwt_Down);
	    }
	    if(NrecoTop>=1 && NrecoV==1){
	      h_Ele_Chi2_1t1V->Fill(chi_output[0],evtwt);
	      h_Ele_Chi2_1t1V_Up->Fill(chi_output[0],evtwt_Up);
	      h_Ele_Chi2_1t1V_Down->Fill(chi_output[0],evtwt_Down);
	    }
	    h_Ele_TMass_chi2_lep->Fill(chi_output[1],evtwt);
	    h_Ele_TMass_chi2_lep_Up->Fill(chi_output[1],evtwt_Up);
	    h_Ele_TMass_chi2_lep_Down->Fill(chi_output[1],evtwt_Down);
	    h_Ele_TMass_chi2_had->Fill(chi_output[2],evtwt);
	    h_Ele_TMass_chi2_had_Up->Fill(chi_output[2],evtwt_Up);
	    h_Ele_TMass_chi2_had_Down->Fill(chi_output[2],evtwt_Down);
	    h_Ele_TMass_chi2->Fill(chi_output[1],evtwt);
	    h_Ele_TMass_chi2_Up->Fill(chi_output[1],evtwt_Up);
	    h_Ele_TMass_chi2_Down->Fill(chi_output[1],evtwt_Down);
	    h_Ele_TMass_chi2->Fill(chi_output[2],evtwt);
	    h_Ele_TMass_chi2_Up->Fill(chi_output[2],evtwt_Up);
	    h_Ele_TMass_chi2_Down->Fill(chi_output[2],evtwt_Down);
	  }
	  if(chi_output[0]<20.){
	    h_Ele_Chi2_cut->Fill(chi_output[0],evtwt);
	    h_Ele_Chi2_cut_Up->Fill(chi_output[0],evtwt_Up);
	    h_Ele_Chi2_cut_Down->Fill(chi_output[0],evtwt_Down);
	    if(NrecoV>=2 && NrecoTop>=1){
	      h_Ele_Chi2_1t_cut->Fill(chi_output[0],evtwt);
	      h_Ele_Chi2_1t_cut_Up->Fill(chi_output[0],evtwt_Up);
	      h_Ele_Chi2_1t_cut_Down->Fill(chi_output[0],evtwt_Down);
	    }
	    if(NrecoTop==1 && !NrecoV && cleanJetAK4.size()>=5){
	      h_Ele_Chi2_1t0V_cut->Fill(chi_output[0],evtwt);
	      h_Ele_Chi2_1t0V_cut_Up->Fill(chi_output[0],evtwt_Up);
	      h_Ele_Chi2_1t0V_cut_Down->Fill(chi_output[0],evtwt_Down);

	    }
	    if(!NrecoTop && NrecoV==1 && cleanJetAK4.size()>=6 ){
	      h_Ele_Chi2_0t1V_cut->Fill(chi_output[0],evtwt);
	      h_Ele_Chi2_0t1V_cut_Up->Fill(chi_output[0],evtwt_Up);
	      h_Ele_Chi2_0t1V_cut_Down->Fill(chi_output[0],evtwt_Down);
	    }
	    if(!NrecoTop && NrecoV>=2 && cleanJetAK4.size()>=4 ){
	      h_Ele_Chi2_0t2V_cut->Fill(chi_output[0]);
	      h_Ele_Chi2_0t2V_cut_Up->Fill(chi_output[0],evtwt_Up);
	      h_Ele_Chi2_0t2V_cut_Down->Fill(chi_output[0],evtwt_Down);

	    }
	    if(NrecoTop>=1 && NrecoV==1){
	      h_Ele_Chi2_1t1V_cut->Fill(chi_output[0]);
	      h_Ele_Chi2_1t1V_cut_Up->Fill(chi_output[0],evtwt_Up);
	      h_Ele_Chi2_1t1V_cut_Down->Fill(chi_output[0],evtwt_Down);
	    }
	    h_Ele_TMass_chi2_lep_cut->Fill(chi_output[1],evtwt);
	    h_Ele_TMass_chi2_lep_cut_Up->Fill(chi_output[1],evtwt_Up);
	    h_Ele_TMass_chi2_lep_cut_Down->Fill(chi_output[1],evtwt_Down);
	    h_Ele_TMass_chi2_had_cut->Fill(chi_output[2],evtwt);
	    h_Ele_TMass_chi2_had_cut_Up->Fill(chi_output[2],evtwt_Up);
	    h_Ele_TMass_chi2_had_cut_Down->Fill(chi_output[2],evtwt_Down);
	    h_Ele_TMass_chi2_cut->Fill(chi_output[1],evtwt);
	    h_Ele_TMass_chi2_cut_Up->Fill(chi_output[1],evtwt_Up);
	    h_Ele_TMass_chi2_cut_Down->Fill(chi_output[1],evtwt_Down);
	    h_Ele_TMass_chi2_cut->Fill(chi_output[2],evtwt);
	    h_Ele_TMass_chi2_cut_Up->Fill(chi_output[2],evtwt_Up);
	    h_Ele_TMass_chi2_cut_Down->Fill(chi_output[2],evtwt_Down);
	  }

	  //}//if(flmu1 &&flmu2)
	}//if(t_HT>300.)
      }//if(goodElid.size()!=0 && (*t_elecP4)[*itEle].Pt()>=30. && (*t_elecP4)[*(itEle+1)].Pt()>=20.)


      
      //============================xxxxxxxx===================================================//

      if(cleanJetAK4.size()>=4){
	h_Cutflow->Fill(3.0);
	if(goodMuid.size())h_Cutflow->Fill(15.0);
	else h_Cutflow->Fill(16.0);
	
	  if(t_HT>500.){
	    h_Cutflow->Fill(5.0);
	    if(goodMuid.size())h_Cutflow->Fill(19.0);
	    else h_Cutflow->Fill(20.0);
	    h_JetsAK8HT_HT->Fill(htAK8);
	    h_JetsAK4HT_HT->Fill(t_HT);
	    
	     
	    h_LeadingAK4Pt_HT->Fill((*t_jetAK4P4)[*itJ].Pt());
	    for (vector<int>::iterator it = cleanJetAK4.begin() ; it != cleanJetAK4.end(); ++it){
	      h_JetsAK4Pt_HT->Fill((*t_jetAK4P4)[*it].Pt());
	    }
	    h_NJetsAK4_HT->Fill(cleanJetAK4.size());
	    for (vector<int>::iterator it = cleanJetAK8.begin() ; it != cleanJetAK8.end(); ++it){
	      h_JetsAK8Pt_HT->Fill((*t_jetAK8P4)[*it].Pt());
	    }
	    h_NJetsAK8_HT->Fill(cleanJetAK8.size());
	      
	    if(cleanJetAK8.size()!=0)h_LeadingAK8Pt_HT->Fill((*t_jetAK8P4)[*cleanJetAK8.begin()].Pt());
	     
	    if(cleanJetAK8.size()>=1){//at least 1 fat jet
	      Njets=cleanJetAK8.size()+cleanJetAK4.size();
	      h_4jetEvent_HT->Fill(Njets);
	    }
	      
	     
	    h_Nlep_HT->Fill(goodMuid.size()+goodElid.size()-2);
	     
	    
	    int bveto=1;int TopWMix=0;
	     
	    if(cleanBJet.size()>=1 && cleanBJetL.size()>=1) {h_bL->Fill(cleanBJetL.size()); h_bM->Fill(cleanBJet.size());}
	     
	    
	     
	    h_NrecobT_HT->Fill(cleanBJetT.size());

	    h_Nrecob_HT->Fill(cleanBJet.size());

	    h_NrecobL_HT->Fill(cleanBJetL.size());

	    h_NrecoTop_HT->Fill(NrecoTop);

	    // if(cleanBJetL.size()==1 && cleanBJet.size()==1){h_Cutflow->Fill(6.0); bveto=0;} // select events with only 1 medium b
	    if(cleanBJetL.size()>1){ h_Cutflow->Fill(7.0);bveto=0;
	      h_HTcheck->Fill(t_HT);
	      h_STcheck->Fill(st);
	      h_HTAK8check->Fill(htAK8);
	      stAK8=st-t_HT+htAK8;
	      h_STAK8check->Fill(stAK8);
	      if(goodMuid.size())h_Cutflow->Fill(21.0);
	      else h_Cutflow->Fill(22.0);
	    } //if at least 2 loose b, then at least 1 medium

	    if(NrecoTop>=2){ //both reco tops matched to gen top
	      h_Nb2t->Fill(cleanBJetL.size());
	      h_NW2t->Fill(t_jetWJetP4->size());
	    }
	    if(!bveto){
	      if(goodMuid.size()!=0){
	      
		vector<int>::iterator itMu = goodMuid.begin();
		h_mu_leadingLepPt->Fill((*t_muonP4)[*itMu].Pt());
		h_mu_leadingLep2Pt->Fill((*t_muonP4)[*(itMu+1)].Pt());
		h_mu_Lep2by1Pt->Fill((*t_muonP4)[*(itMu+1)].Pt()/(*t_muonP4)[*itMu].Pt());
		h_mu_Lep2vs1Pt->Fill((*t_muonP4)[*(itMu+1)].Pt(),(*t_muonP4)[*itMu].Pt());
		h_mu_Lep2Ptfrac->Fill((*t_muonP4)[*(itMu+1)].Pt()/((*t_muonP4)[*itMu].Pt()+(*t_muonP4)[*(itMu+1)].Pt()));
		px=(*t_muonP4)[*itMu].Px()+(*t_muonP4)[*(itMu+1)].Px();
		py=(*t_muonP4)[*itMu].Py()+(*t_muonP4)[*(itMu+1)].Py();
		pt=sqrt(px*px+py*py);
	      
		h_dimupt_HT->Fill(pt);
	      }
	     
	      if(goodElid.size()!=0){
		vector<int>::iterator itEle = goodElid.begin();
		h_ele_leadingLepPt->Fill((*t_elecP4)[*itEle].Pt());
		h_ele_leadingLep2Pt->Fill((*t_elecP4)[*(itEle+1)].Pt());
		h_ele_Lep2by1Pt->Fill((*t_elecP4)[*(itEle+1)].Pt()/(*t_elecP4)[*itEle].Pt());
		h_ele_Lep2vs1Pt->Fill((*t_elecP4)[*(itEle+1)].Pt(),(*t_elecP4)[*itEle].Pt());
		h_ele_Lep2Ptfrac->Fill((*t_elecP4)[*(itEle+1)].Pt()/((*t_elecP4)[*itEle].Pt()+(*t_elecP4)[*(itEle+1)].Pt()));
		px=(*t_elecP4)[*itEle].Px()+(*t_elecP4)[*(itEle+1)].Px();
		py=(*t_elecP4)[*itEle].Py()+(*t_elecP4)[*(itEle+1)].Py();
		pt=sqrt(px*px+py*py);
		h_dielept_HT->Fill(pt);
	       
	      }
	     
	      int b;
	      int W;
	      TLorentzVector Wb,TP,AKV;
	      double chi_1,chi_2,chi_3,chi,chi_cf,chi_1_a,chi_1_b,chi_2_a,chi_2_b,chi_1_c,chi_2_d,chi_2_e,chi_2_c,chi_2_f;
	      int qj,qpj;
	      double dr1,dr2;
	      dr1=99;dr2=99;
	      TLorentzVector TP1,TP2,Wb1,Wb2,V,Wbc;
	      double dphi1;
	      double wq,wqp;
	      wq=99;wqp=99;
	      //cout<<"starting reco\n";
	      //================================ Top matching with gen ====================================================// 
	      int foundTop=0;
	      int foundWb=0;
	      for(unsigned int i=0; i<t_genPartP4->size(); i++) {
		if(fabs((*t_genPartID)[i])==6){
		  dr=99;jetindex=99;
		  for (vector<int>::iterator it = cleanJetAK8.begin() ; it != cleanJetAK8.end(); ++it){
		    if(fabs((*t_jetAK8P4)[*it].Eta())<2.4 && (*t_jetAK8_tau3)[*it]/(*t_jetAK8_tau2)[*it]<0.67 &&  (*t_jetAK8P4)[*it].Pt()>400.0){
		      phi1=(*t_jetAK8P4)[*it].Phi();
		      eta1=(*t_jetAK8P4)[*it].Eta();
		      phi2=(*t_genPartP4)[i].Phi();
		      eta2=(*t_genPartP4)[i].Eta();
		      dR=DeltaR(eta1,phi1,eta2,phi2);
		      if(dR<0.3 && dR<dr){
			dr=dR;
			jetindex=*it;
			//     topId=i;
		      }
		    }
		  }
		  if(jetindex!=99){
		    //cout<<"Case 1\n";
		    qj=99;W=99;
		    for(unsigned int q=i+1; q<t_genPartP4->size(); q++) {
		      if(fabs((*t_genPartID)[q])==5 && ((*t_genPartMom1ID)[q]==(*t_genPartID)[i] || (*t_genPartMom2ID)[q]==(*t_genPartID)[i])) qj=q;
    
		      if(fabs((*t_genPartID)[q])==24 && ((*t_genPartMom1ID)[q]==(*t_genPartID)[i] || (*t_genPartMom2ID)[q]==(*t_genPartID)[i]))W=q;		     
		    }
		    if(qj!=99 && W!=99){
		      phi1=(*t_jetAK8P4)[jetindex].Phi();
		      eta1=(*t_jetAK8P4)[jetindex].Eta();
		      phi2=(*t_genPartP4)[qj].Phi();
		      eta2=(*t_genPartP4)[qj].Eta();
		      dR=DeltaR(eta1,phi1,eta2,phi2);
		      if(dR <0.4){
			phi2=(*t_genPartP4)[W].Phi();
			eta2=(*t_genPartP4)[W].Eta();
			dr1=DeltaR(eta1,phi1,eta2,phi2);
			if(dr1<0.4){
			  for(unsigned int w=W+1; w<t_genPartP4->size(); w++) {
			    if(fabs((*t_genPartID)[w])<=5 && fabs((*t_genPartID)[w])>=1 && ((*t_genPartMom1ID)[w]==(*t_genPartID)[W] || (*t_genPartMom2ID)[w]==(*t_genPartID)[W])){ 
			      if(wq==99) wq=w;
			      else wqp=w;
			    }
			  }
			  if(wq!=99 && wqp!=99){
			    phi1=(*t_jetAK8P4)[jetindex].Phi();
			    eta1=(*t_jetAK8P4)[jetindex].Eta();
			    phi2=(*t_genPartP4)[wq].Phi();
			    eta2=(*t_genPartP4)[wq].Eta();
			    dr=DeltaR(eta1,phi1,eta2,phi2);
			    if(dr<0.4){
			      phi2=(*t_genPartP4)[wqp].Phi();
			      eta2=(*t_genPartP4)[wqp].Eta();
			      if(dr1<0.4){
				h_topMass->Fill((*t_jetAK8_SoftDropMass)[jetindex]);
				foundTop=1;
			      }
			    }
			  }
			}
		      }
		    }
		  }
		    
		  //========================== W+b matching with gen top======================================//
	       
		  else{ //iff no AK8 jet directly matches the gen top, look for W+b
		    //cout<<"case 2\n";
		    dr=99;b=99;W=99;
		    for(unsigned int w=0; w<t_genPartP4->size(); w++) {
		      if(fabs((*t_genPartID)[w])==24 && ((*t_genPartMom1ID)[w]==(*t_genPartID)[i] || (*t_genPartMom2ID)[w]==(*t_genPartID)[i]) ){
			dr=99;jetindex=99;
			for (unsigned int it = 0 ; it<t_jetWJetP4->size(); ++it){
			  phi1=(*t_jetWJetP4)[it].Phi();
			  eta1=(*t_jetWJetP4)[it].Eta();
			  phi2=(*t_genPartP4)[w].Phi();
			  eta2=(*t_genPartP4)[w].Eta();
			  dR=DeltaR(eta1,phi1,eta2,phi2);
			  if(dR<0.3 && dR<dr){
			    dr=dR;
			    W=it;
			  }
			}
			if(W!=99){
			  dr=99;b=99;qj=99;qpj=99;
			  for(unsigned int q=w+1; q<t_genPartP4->size(); q++) {
			    if(fabs((*t_genPartID)[q])<6 && fabs((*t_genPartID)[q])>=1 && ((*t_genPartMom1ID)[q]==(*t_genPartID)[w] || (*t_genPartMom2ID)[q]==(*t_genPartID)[w])){
			      if(qj==99)qj=q;
			      else qpj=q;
			    }
			  }
			  //cout<<qj<<" "<<qpj<<endl;
			  double dr1;
			  int findb=0;
			  if(qj!=99 && qpj!=99){
			    phi1=(*t_jetAK8P4)[W].Phi();
			    eta1=(*t_jetAK8P4)[W].Eta();
			    phi2=(*t_genPartP4)[qj].Phi();
			    eta2=(*t_genPartP4)[qj].Eta();
			    dR=DeltaR(eta1,phi1,eta2,phi2);
			    if(dR <0.6){
			      phi2=(*t_genPartP4)[qpj].Phi();
			      eta2=(*t_genPartP4)[qpj].Eta();
			      dr1=DeltaR(eta1,phi1,eta2,phi2);
			      if(dr1<0.6)findb=1;
			    }
			  }
			  if(findb){
			    //cout<<"cool\n";
			    for(unsigned int bj=0; bj<t_genPartP4->size(); bj++) {
			      if(fabs((*t_genPartID)[bj])==5 && ((*t_genPartMom1ID)[bj]==(*t_genPartID)[i] || (*t_genPartMom2ID)[bj]==(*t_genPartID)[i])){
			     
				for (vector<int>::iterator itb = cleanBJetL.begin() ; itb != cleanBJetL.end(); ++itb){ 
				  phi1=(*t_jetAK4P4)[*itb].Phi();
				  eta1=(*t_jetAK4P4)[*itb].Eta();
				  phi2=(*t_genPartP4)[bj].Phi();
				  eta2=(*t_genPartP4)[bj].Eta();
				  dR=DeltaR(eta1,phi1,eta2,phi2);
				  if(dR<0.4 && dR<dr){
				    dr=dR;
				    b=*itb;
				  }
				}
				break;
			      }
			    }
			    if(b!=99){
			      Wb=(*t_jetAK8P4)[W]+(*t_jetAK4P4)[b];
			      h_topMassW->Fill(Wb.M());
			      foundWb=1;
			    }
			  }//if(findb)
			}//if(W!=99){
			break;
		      }//if(fabs((*t_genPartID)[w])==24 && ((*t_genPartMom1ID)[w]==(*t_genP.....
		    }//for(unsigned int w=0; w<t_genPartP4->size(); w++) {
		  }//else{
		}//if(fabs((*t_genPartID)[i])==6){
	      }//for(unsigned int i=0; i<t_genPartP4->size(); i++) {
	      int q1;int q2;
	      //if(!foundTop && !foundWb){
	      b=99;q1=99;q2=99;
	      //cout<<"resolved\n";
	      for(unsigned int i=0; i<t_genPartP4->size(); i++) {
		if(fabs((*t_genPartID)[i])==6){
		  b=99;q1=99;q2=99;
		  for(unsigned int q=i+1; q<t_genPartP4->size(); q++) {
	      
		    if(b==99 && fabs((*t_genPartID)[q])==5 && ((*t_genPartMom1ID)[q]==(*t_genPartID)[i] || (*t_genPartMom2ID)[q]==(*t_genPartID)[i])) {
		      dr=99;b=99;
		      for (vector<int>::iterator itb = cleanBJetL.begin() ; itb != cleanBJetL.end(); ++itb){ 
			phi1=(*t_jetAK4P4)[*itb].Phi();
			 
			eta1=(*t_jetAK4P4)[*itb].Eta();
			 
			phi2=(*t_genPartP4)[q].Phi();
			eta2=(*t_genPartP4)[q].Eta();
			//cout<<"no\n";
			dR=DeltaR(eta1,phi1,eta2,phi2);
			if(dR<0.4 && dR<dr){
			  if((q1!=99 || q2!=99) && (q1!=*itb || q2!=*itb)){ 
			    dr=dR;
			    b=*itb;
			  }
			  else{
			    dr=dR;
			    b=*itb;			     
			  }
			}
		      }
		      // continue;
		    }
		     
		    if(q1==99 && q2==99 && fabs((*t_genPartID)[q])==24 && ((*t_genPartMom1ID)[q]==(*t_genPartID)[i] || (*t_genPartMom2ID)[q]==(*t_genPartID)[i])){
		      qj=99;qpj=99;	 
		      for(unsigned int wd=q+1; wd<t_genPartP4->size(); wd++) {
			if(fabs((*t_genPartID)[wd])<6 && fabs((*t_genPartID)[wd])>=1 && ((*t_genPartMom1ID)[wd]==(*t_genPartID)[q] || (*t_genPartMom2ID)[wd]==(*t_genPartID)[q])){
			  if(qj==99)qj=q;
			  else qpj=q;
			}
		      }
		      if(qj!=99 && qpj!=99){
			dr=99;
			for (vector<int>::iterator itAK1 = cleanJetAK4.begin() ;itAK1!=cleanJetAK4.end() ; ++itAK1){
			  phi1=(*t_jetAK4P4)[*itAK1].Phi();
			  eta1=(*t_jetAK4P4)[*itAK1].Eta();
			  phi2=(*t_genPartP4)[qj].Phi();
			  eta2=(*t_genPartP4)[qj].Eta();
			  dR=DeltaR(eta1,phi1,eta2,phi2);
			  if(dR<0.4 && dR<dr){
			    if(b!=99){
			      if(b!=*itAK1){dr=dR;q1=*itAK1;}
			    }
			    else{dr=dR; q1=*itAK1;}
			    //cout<<"q1:"<<q1<<endl;
			  }
			}
			if(q1!=99){
			  dr=99;
			  for (vector<int>::iterator itAK2 = cleanJetAK4.begin() ;itAK2!=cleanJetAK4.end() ; ++itAK2){
			    if(*itAK2!=q1){
			      phi1=(*t_jetAK4P4)[*itAK2].Phi();
			      eta1=(*t_jetAK4P4)[*itAK2].Eta();
			      phi2=(*t_genPartP4)[qpj].Phi();
			      eta2=(*t_genPartP4)[qpj].Eta();
			      dR=DeltaR(eta1,phi1,eta2,phi2);
			      if(dR<0.4 && dR<dr){
				if(b!=99){
				  if(b!=*itAK2){dr=dR;q2=*itAK2;}
				}
				else{dr=dR;q2=*itAK2;}
			      }
			    }
			  }
			}
		
		      }
		       
		    }
		    if(q1!=99 && q2 !=99 && b!=99){
		      //cout<<"ok\n";
		      //cout<<q1<<" "<<q2<<endl;
		      V=(*t_jetAK4P4)[q1]+(*t_jetAK4P4)[q2];
		      //cout<<"ooh\n";
		      if(V.M() <105 && V.M()>65){
			phi1=V.Phi();
			phi2=V.Eta();
			phi2=(*t_genPartP4)[q].Phi();
			eta2=(*t_genPartP4)[q].Eta();
			dR=DeltaR(eta1,phi1,eta2,phi2);
			if(dR<0.4 ){
			  Wb=V+(*t_jetAK4P4)[b];
			  phi1=Wb.Phi();
			  phi2=Wb.Eta();
			  phi2=(*t_genPartP4)[i].Phi();
			  eta2=(*t_genPartP4)[i].Eta();
			  if(dR<0.4 ){
			    h_resTop->Fill(Wb.M());
			    break;
			  }
			}
		      }
		    }
		  }

		}
	      }
	       
	      dr=99;
	      int match[]={99,99};
	      int tcand[]={99,99};
	       
	      int Vcheck=0; int Vidx=-99;


	      //===================== to calculate the chi2==============================//
	      chi=10000;chi_cf=10000;
	      chi_1=0;chi_2=0;
	      int tidx=-99;int tidx1=-99;int Wbidx=-99;int Wbidx1=-99; int bidx1=-99; int bidx=-99;
	      int ctr1=0;int ctr2=0;int ctr3=0;int ctrb=0;int ctr4=0;int ctr5=0;
	      double diff1,diff2,diff3,diff4;
	      double chi_V2[12];
	      int chi_idx[]={0,1,2,3,4,5,6,7,8,9,10,11};
	      int type[]={0,0,0,0,0,0};
	      for(int i=0;i<12;i++)chi_V2[i]=10000;
	      TLorentzVector TPfl,TPfh,TPl1,TPh1,TPl2,TPh2, TPl3,TPh3,TPl4,TPh4,Wbf1,Wbf2,AKVf,V1,V2;
	      chi_output=GetChi(NrecoTop,NrecoV,v1,goodMuid,goodElid, cleanJetAK4, cleanBJetL);
	      //cout<<chi_output[0];
	      if(chi_output[0]!=10000){
		//cout<<"meow\n";
		h_Chi2->Fill(chi_output[0]);
		if(NrecoV>=2 && NrecoTop>=1)h_Chi2_1t->Fill(chi_output[0]);
		if(NrecoTop==1 && !NrecoV && cleanJetAK4.size()>=5)h_Chi2_1t0V->Fill(chi_output[0]);
		if(!NrecoTop && NrecoV==1 && cleanJetAK4.size()>=6 )h_Chi2_0t1V->Fill(chi_output[0]);
		if(!NrecoTop && NrecoV>=2 && cleanJetAK4.size()>=4 )h_Chi2_0t2V->Fill(chi_output[0]);
		if(NrecoTop>=1 && NrecoV==1)h_Chi2_1t1V->Fill(chi_output[0]);
		h_TMass_chi2_lep->Fill(chi_output[1]);
		h_TMass_chi2_had->Fill(chi_output[2]);
		h_TMass_chi2->Fill(chi_output[1]);
		h_TMass_chi2->Fill(chi_output[2]);
	      }
	      if(chi_output[0]<20.){
		h_Chi2_cut->Fill(chi_output[0]);
		if(NrecoV>=2 && NrecoTop>=1)h_Chi2_1t_cut->Fill(chi_output[0]);
                if(NrecoTop==1 && !NrecoV && cleanJetAK4.size()>=5)h_Chi2_1t0V_cut->Fill(chi_output[0]);
                if(!NrecoTop && NrecoV==1 && cleanJetAK4.size()>=6 )h_Chi2_0t1V_cut->Fill(chi_output[0]);
                if(!NrecoTop && NrecoV>=2 && cleanJetAK4.size()>=4 )h_Chi2_0t2V_cut->Fill(chi_output[0]);
                if(NrecoTop>=1 && NrecoV==1)h_Chi2_1t1V_cut->Fill(chi_output[0]);
		h_TMass_chi2_lep_cut->Fill(chi_output[1]);
		h_TMass_chi2_had_cut->Fill(chi_output[2]);
		h_TMass_chi2_cut->Fill(chi_output[1]);
		h_TMass_chi2_cut->Fill(chi_output[2]);
	      }
	      

	      evtwt=t_evtwt;
	      evtwt_Up=t_evtwt  ;
	      evtwt_Down=t_evtwt  ;
	      vector<int>::iterator itMu = goodMuid.begin();
	      if(goodMuid.size()!=0 && (*t_muonP4)[*itMu].Pt()>=30. && (*t_muonP4)[*(itMu+1)].Pt()>=20.){
		if((*t_muonP4)[*itMu].Pt()>=500. && fabs((*t_muonP4)[*itMu].Eta())>1.8)flmu1=0; //remove endcap muons with high pt
		if((*t_muonP4)[*(itMu+1)].Pt()>=500. && fabs((*t_muonP4)[*(itMu+1)].Eta())>1.8)flmu2=0;
	  
		if(flmu1 && flmu2){
		 
		   
		  //if(evtwt_Up<evtwt_Down)cout<<evtwt<<" "<<evtwt_Up<<" "<<evtwt_Down<<endl;
		  chi_output=GetChi(NrecoTop,NrecoV,v1,goodMuid,goodElid, cleanJetAK4, cleanBJetL);
		  //cout<<"Mu ok\n";
		  //cout<<chi_output[0];
		  if(chi_output[0]!=10000){
		    //cout<<"meow\n";
		    h_UB_Mu_Chi2->Fill(chi_output[0],evtwt);
		    h_UB_Mu_Chi2_Up->Fill(chi_output[0],evtwt_Up);
		    h_UB_Mu_Chi2_Down->Fill(chi_output[0],evtwt_Down);
		    //cout<<"hmm\n";
		    if(NrecoV>=2 && NrecoTop>=1){
		      h_UB_Mu_Chi2_1t->Fill(chi_output[0],evtwt);
		      h_UB_Mu_Chi2_1t_Up->Fill(chi_output[0],evtwt_Up);
		      h_UB_Mu_Chi2_1t_Down->Fill(chi_output[0],evtwt_Down);
			
		    }
		    if(NrecoTop==1 && !NrecoV && cleanJetAK4.size()>=5){
		      h_UB_Mu_Chi2_1t0V->Fill(chi_output[0],evtwt);
		      h_UB_Mu_Chi2_1t0V_Up->Fill(chi_output[0],evtwt_Up);
		      h_UB_Mu_Chi2_1t0V_Down->Fill(chi_output[0],evtwt_Down);
		    }
		    if(!NrecoTop && NrecoV==1 && cleanJetAK4.size()>=6 ){
		      h_UB_Mu_Chi2_0t1V->Fill(chi_output[0],evtwt);
		      h_UB_Mu_Chi2_0t1V_Up->Fill(chi_output[0],evtwt_Up);
		      h_UB_Mu_Chi2_0t1V_Down->Fill(chi_output[0],evtwt_Down);
		    }
		    if(!NrecoTop && NrecoV>=2 && cleanJetAK4.size()>=4 ){
		      h_UB_Mu_Chi2_0t2V->Fill(chi_output[0],evtwt);
		      h_UB_Mu_Chi2_0t2V_Up->Fill(chi_output[0],evtwt_Up);
		      h_UB_Mu_Chi2_0t2V_Down->Fill(chi_output[0],evtwt_Down);
		    }
		    if(NrecoTop>=1 && NrecoV==1){
		      h_UB_Mu_Chi2_1t1V->Fill(chi_output[0],evtwt);
		      h_UB_Mu_Chi2_1t1V_Up->Fill(chi_output[0],evtwt_Up);
		      h_UB_Mu_Chi2_1t1V_Down->Fill(chi_output[0],evtwt_Down);
		    }
		    h_UB_Mu_TMass_chi2_lep->Fill(chi_output[1],evtwt);
		    h_UB_Mu_TMass_chi2_lep_Up->Fill(chi_output[1],evtwt_Up);
		    h_UB_Mu_TMass_chi2_lep_Down->Fill(chi_output[1],evtwt_Down);
		    h_UB_Mu_TMass_chi2_had->Fill(chi_output[2],evtwt);
		    h_UB_Mu_TMass_chi2_had_Up->Fill(chi_output[2],evtwt_Up);
		    h_UB_Mu_TMass_chi2_had_Down->Fill(chi_output[2],evtwt_Down);
		    h_UB_Mu_TMass_chi2->Fill(chi_output[1],evtwt);
		    h_UB_Mu_TMass_chi2_Up->Fill(chi_output[1],evtwt_Up);
		    h_UB_Mu_TMass_chi2_Down->Fill(chi_output[1],evtwt_Down);
		    h_UB_Mu_TMass_chi2->Fill(chi_output[2],evtwt);
		    h_UB_Mu_TMass_chi2_Up->Fill(chi_output[2],evtwt_Up);
		    h_UB_Mu_TMass_chi2_Down->Fill(chi_output[2],evtwt_Down);
		  }
		  if(chi_output[0]<20.){
		    h_UB_Mu_Chi2_cut->Fill(chi_output[0],evtwt);
		    h_UB_Mu_Chi2_cut_Up->Fill(chi_output[0],evtwt_Up);
		    h_UB_Mu_Chi2_cut_Down->Fill(chi_output[0],evtwt_Down);
		    if(NrecoV>=2 && NrecoTop>=1){
		      h_UB_Mu_Chi2_1t_cut->Fill(chi_output[0],evtwt);
		      h_UB_Mu_Chi2_1t_cut_Up->Fill(chi_output[0],evtwt_Up);
		      h_UB_Mu_Chi2_1t_cut_Down->Fill(chi_output[0],evtwt_Down);
		    }
		    if(NrecoTop==1 && !NrecoV && cleanJetAK4.size()>=5){
		      h_UB_Mu_Chi2_1t0V_cut->Fill(chi_output[0],evtwt);
		      h_UB_Mu_Chi2_1t0V_cut_Up->Fill(chi_output[0],evtwt_Up);
		      h_UB_Mu_Chi2_1t0V_cut_Down->Fill(chi_output[0],evtwt_Down);

		    }
		    if(!NrecoTop && NrecoV==1 && cleanJetAK4.size()>=6 ){
		      h_UB_Mu_Chi2_0t1V_cut->Fill(chi_output[0],evtwt);
		      h_UB_Mu_Chi2_0t1V_cut_Up->Fill(chi_output[0],evtwt_Up);
		      h_UB_Mu_Chi2_0t1V_cut_Down->Fill(chi_output[0],evtwt_Down);
		    }
		    if(!NrecoTop && NrecoV>=2 && cleanJetAK4.size()>=4 ){
		      h_UB_Mu_Chi2_0t2V_cut->Fill(chi_output[0]);
		      h_UB_Mu_Chi2_0t2V_cut_Up->Fill(chi_output[0],evtwt_Up);
		      h_UB_Mu_Chi2_0t2V_cut_Down->Fill(chi_output[0],evtwt_Down);

		    }
		    if(NrecoTop>=1 && NrecoV==1){
		      h_UB_Mu_Chi2_1t1V_cut->Fill(chi_output[0]);
		      h_UB_Mu_Chi2_1t1V_cut_Up->Fill(chi_output[0],evtwt_Up);
		      h_UB_Mu_Chi2_1t1V_cut_Down->Fill(chi_output[0],evtwt_Down);
		    }
		    h_UB_Mu_TMass_chi2_lep_cut->Fill(chi_output[1],evtwt);
		    h_UB_Mu_TMass_chi2_lep_cut_Up->Fill(chi_output[1],evtwt_Up);
		    h_UB_Mu_TMass_chi2_lep_cut_Down->Fill(chi_output[1],evtwt_Down);
		    h_UB_Mu_TMass_chi2_had_cut->Fill(chi_output[2],evtwt);
		    h_UB_Mu_TMass_chi2_had_cut_Up->Fill(chi_output[2],evtwt_Up);
		    h_UB_Mu_TMass_chi2_had_cut_Down->Fill(chi_output[2],evtwt_Down);
		    h_UB_Mu_TMass_chi2_cut->Fill(chi_output[1],evtwt);
		    h_UB_Mu_TMass_chi2_cut_Up->Fill(chi_output[1],evtwt_Up);
		    h_UB_Mu_TMass_chi2_cut_Down->Fill(chi_output[1],evtwt_Down);

		    h_UB_Mu_TMass_chi2_cut->Fill(chi_output[2],evtwt);
		    h_UB_Mu_TMass_chi2_cut_Up->Fill(chi_output[2],evtwt_Up);
		    h_UB_Mu_TMass_chi2_cut_Down->Fill(chi_output[2],evtwt_Down);
		  }
		 
		}//if(flmu1 &&flmu2)
	      }//if(goodMuid.size()!=0 && (*t_muonP4)[*itMu].Pt()>=30. && (*t_muonP4)[*(itMu+1)].Pt()>=20.)
	      //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
	      vector<int>::iterator itEle = goodElid.begin();
	      check=1;
	      if(goodElid.size()!=0 && (*t_elecP4)[*itEle].Pt()>=30. && (*t_elecP4)[*(itEle+1)].Pt()>=20.){
		//if(!(*t_elecIsTight)[*itEle])cout<<(*t_elecIsTight)[*itEle]<<" "<<(*t_elecIsTight)[(*itEle+1)]<<endl;
	
	
		chi_output=GetChi(NrecoTop,NrecoV,v1,goodMuid,goodElid, cleanJetAK4, cleanBJetL);
		//if(evtwt_Up<evtwt_Down)cout<<evtwt<<" "<<evtwt_Up<<" "<<evtwt_Down<<endl;
		//cout<<"El ok\n";
		//cout<<chi_output[0];
		if(chi_output[0]!=10000){
		  //cout<<"meow1\n";
		  h_UB_Ele_Chi2->Fill(chi_output[0],evtwt);
		  h_UB_Ele_Chi2_Up->Fill(chi_output[0],evtwt_Up);
		  h_UB_Ele_Chi2_Down->Fill(chi_output[0],evtwt_Down);
		  //cout<<"ho\n";
		  if(NrecoV>=2 && NrecoTop>=1){
		    h_UB_Ele_Chi2_1t->Fill(chi_output[0],evtwt);
		    h_UB_Ele_Chi2_1t_Up->Fill(chi_output[0],evtwt_Up);
		    h_UB_Ele_Chi2_1t_Down->Fill(chi_output[0],evtwt_Down);

		  }
		  if(NrecoTop==1 && !NrecoV && cleanJetAK4.size()>=5){
		    h_UB_Ele_Chi2_1t0V->Fill(chi_output[0],evtwt);
		    h_UB_Ele_Chi2_1t0V_Up->Fill(chi_output[0],evtwt_Up);
		    h_UB_Ele_Chi2_1t0V_Down->Fill(chi_output[0],evtwt_Down);
		  }
		  if(!NrecoTop && NrecoV==1 && cleanJetAK4.size()>=6 ){
		    h_UB_Ele_Chi2_0t1V->Fill(chi_output[0],evtwt);
		    h_UB_Ele_Chi2_0t1V_Up->Fill(chi_output[0],evtwt_Up);
		    h_UB_Ele_Chi2_0t1V_Down->Fill(chi_output[0],evtwt_Down);
		  }
		  if(!NrecoTop && NrecoV>=2 && cleanJetAK4.size()>=4 ){
		    h_UB_Ele_Chi2_0t2V->Fill(chi_output[0],evtwt);
		    h_UB_Ele_Chi2_0t2V_Up->Fill(chi_output[0],evtwt_Up);
		    h_UB_Ele_Chi2_0t2V_Down->Fill(chi_output[0],evtwt_Down);
		  }
		  if(NrecoTop>=1 && NrecoV==1){
		    h_UB_Ele_Chi2_1t1V->Fill(chi_output[0],evtwt);
		    h_UB_Ele_Chi2_1t1V_Up->Fill(chi_output[0],evtwt_Up);
		    h_UB_Ele_Chi2_1t1V_Down->Fill(chi_output[0],evtwt_Down);
		  }
		  h_UB_Ele_TMass_chi2_lep->Fill(chi_output[1],evtwt);
		  h_UB_Ele_TMass_chi2_lep_Up->Fill(chi_output[1],evtwt_Up);
		  h_UB_Ele_TMass_chi2_lep_Down->Fill(chi_output[1],evtwt_Down);
		  h_UB_Ele_TMass_chi2_had->Fill(chi_output[2],evtwt);
		  h_UB_Ele_TMass_chi2_had_Up->Fill(chi_output[2],evtwt_Up);
		  h_UB_Ele_TMass_chi2_had_Down->Fill(chi_output[2],evtwt_Down);
		  h_UB_Ele_TMass_chi2->Fill(chi_output[1],evtwt);
		  h_UB_Ele_TMass_chi2_Up->Fill(chi_output[1],evtwt_Up);
		  h_UB_Ele_TMass_chi2_Down->Fill(chi_output[1],evtwt_Down);
		  h_UB_Ele_TMass_chi2->Fill(chi_output[2],evtwt);
		  h_UB_Ele_TMass_chi2_Up->Fill(chi_output[2],evtwt_Up);
		  h_UB_Ele_TMass_chi2_Down->Fill(chi_output[2],evtwt_Down);
		}
		if(chi_output[0]<20.){
		  h_UB_Ele_Chi2_cut->Fill(chi_output[0],evtwt);
		  h_UB_Ele_Chi2_cut_Up->Fill(chi_output[0],evtwt_Up);
		  h_UB_Ele_Chi2_cut_Down->Fill(chi_output[0],evtwt_Down);
		  if(NrecoV>=2 && NrecoTop>=1){
		    h_UB_Ele_Chi2_1t_cut->Fill(chi_output[0],evtwt);
		    h_UB_Ele_Chi2_1t_cut_Up->Fill(chi_output[0],evtwt_Up);
		    h_UB_Ele_Chi2_1t_cut_Down->Fill(chi_output[0],evtwt_Down);
		  }
		  if(NrecoTop==1 && !NrecoV && cleanJetAK4.size()>=5){
		    h_UB_Ele_Chi2_1t0V_cut->Fill(chi_output[0],evtwt);
		    h_UB_Ele_Chi2_1t0V_cut_Up->Fill(chi_output[0],evtwt_Up);
		    h_UB_Ele_Chi2_1t0V_cut_Down->Fill(chi_output[0],evtwt_Down);

		  }
		  if(!NrecoTop && NrecoV==1 && cleanJetAK4.size()>=6 ){
		    h_UB_Ele_Chi2_0t1V_cut->Fill(chi_output[0],evtwt);
		    h_UB_Ele_Chi2_0t1V_cut_Up->Fill(chi_output[0],evtwt_Up);
		    h_UB_Ele_Chi2_0t1V_cut_Down->Fill(chi_output[0],evtwt_Down);
		  }
		  if(!NrecoTop && NrecoV>=2 && cleanJetAK4.size()>=4 ){
		    h_UB_Ele_Chi2_0t2V_cut->Fill(chi_output[0]);
		    h_UB_Ele_Chi2_0t2V_cut_Up->Fill(chi_output[0],evtwt_Up);
		    h_UB_Ele_Chi2_0t2V_cut_Down->Fill(chi_output[0],evtwt_Down);

		  }
		  if(NrecoTop>=1 && NrecoV==1){
		    h_UB_Ele_Chi2_1t1V_cut->Fill(chi_output[0]);
		    h_UB_Ele_Chi2_1t1V_cut_Up->Fill(chi_output[0],evtwt_Up);
		    h_UB_Ele_Chi2_1t1V_cut_Down->Fill(chi_output[0],evtwt_Down);
		  }
		  h_UB_Ele_TMass_chi2_lep_cut->Fill(chi_output[1],evtwt);
		  h_UB_Ele_TMass_chi2_lep_cut_Up->Fill(chi_output[1],evtwt_Up);
		  h_UB_Ele_TMass_chi2_lep_cut_Down->Fill(chi_output[1],evtwt_Down);
		  h_UB_Ele_TMass_chi2_had_cut->Fill(chi_output[2],evtwt);
		  h_UB_Ele_TMass_chi2_had_cut_Up->Fill(chi_output[2],evtwt_Up);
		  h_UB_Ele_TMass_chi2_had_cut_Down->Fill(chi_output[2],evtwt_Down);
		  h_UB_Ele_TMass_chi2_cut->Fill(chi_output[1],evtwt);
		  h_UB_Ele_TMass_chi2_cut_Up->Fill(chi_output[1],evtwt_Up);
		  h_UB_Ele_TMass_chi2_cut_Down->Fill(chi_output[1],evtwt_Down);
		  h_UB_Ele_TMass_chi2_cut->Fill(chi_output[2],evtwt);
		  h_UB_Ele_TMass_chi2_cut_Up->Fill(chi_output[2],evtwt_Up);
		  h_UB_Ele_TMass_chi2_cut_Down->Fill(chi_output[2],evtwt_Down);
		}

	      }//if(goodElid.size()!=0 && (*t_elecP4)[*itEle].Pt()>=30. && (*t_elecP4)[*(itEle+1)].Pt()>=20.)

	      //========================= To find Vjet Mean and sigma ===============================//
	      for(unsigned int i=0; i<t_genPartP4->size(); i++) {
		if(fabs((*t_genPartID)[i])==23){
		  dr=99;Vidx=99;
		  for (vector<int>::iterator it = cleanJetAK8.begin() ; it != cleanJetAK8.end(); ++it){
		    if(fabs((*t_jetAK8P4)[*it].Eta())<2.4 && (*t_jetAK8_tau2)[*it]/(*t_jetAK8_tau1)[*it]<0.6 && (*t_jetAK8P4)[*it].Pt()>200.0 ){
		      phi1=(*t_jetAK8P4)[*it].Phi();
		      eta1=(*t_jetAK8P4)[*it].Eta();
		      phi2=(*t_genPartP4)[i].Phi();
		      eta2=(*t_genPartP4)[i].Eta();
		      dR=DeltaR(eta1,phi1,eta2,phi2);
		      if(dR<0.3 && dR<dr){
			dr=dR;
			Vidx=*it;
		      }
		    }
		  }
		  if(Vidx!=99){
		    qj=99;qpj=99;
                    for(unsigned int q=Vidx+1; q<t_genPartP4->size(); q++) {
                      if(fabs((*t_genPartID)[q])<6 && fabs((*t_genPartID)[q])>=1 && ((*t_genPartMom1ID)[q]==(*t_genPartID)[i] || (*t_genPartMom2ID)[q]==(*t_genPartID)[i])){
			if(qj==99)qj=q;
			else qpj=q;
		      }
		    }
		    if(qj!=99 && qpj!=99){
		      phi1=(*t_jetAK8P4)[Vidx].Phi();
		      eta1=(*t_jetAK8P4)[Vidx].Eta();
		      phi2=(*t_genPartP4)[qj].Phi();
		      eta2=(*t_genPartP4)[qj].Eta();
		      dR=DeltaR(eta1,phi1,eta2,phi2);
		      if(dR <0.4){
			phi2=(*t_genPartP4)[qpj].Phi();
			eta2=(*t_genPartP4)[qpj].Eta();
			dr1=DeltaR(eta1,phi1,eta2,phi2);
			if(dr1<0.4) h_ZMass->Fill((*t_jetAK8_MassPruned)[Vidx]);
		      }
		    }
		  }
		}
	      }
	      //========================= To find Vjet Mean and sigma ===============================//
	      for(unsigned int i=0; i<t_genPartP4->size(); i++) {
		if(fabs((*t_genPartID)[i])==24){
		  dr=99;Vidx=99;
		  for (vector<int>::iterator it = cleanJetAK8.begin() ; it != cleanJetAK8.end(); ++it){
		    if(fabs((*t_jetAK8P4)[*it].Eta())<2.4 && (*t_jetAK8_tau2)[*it]/(*t_jetAK8_tau1)[*it]<0.6 && (*t_jetAK8P4)[*it].Pt()>200.){
		      phi1=(*t_jetAK8P4)[*it].Phi();
		      eta1=(*t_jetAK8P4)[*it].Eta();
		      phi2=(*t_genPartP4)[i].Phi();
		      eta2=(*t_genPartP4)[i].Eta();
		      dR=DeltaR(eta1,phi1,eta2,phi2);
		      if(dR<0.3 && dR<dr){
			dr=dR;
			Vidx=*it;
		      }
		    }
		  }
		  if(Vidx!=99){
		    qj=99;qpj=99;
		    for(unsigned int q=Vidx+1; q<t_genPartP4->size(); q++) {
		      if(fabs((*t_genPartID)[q])<6 && fabs((*t_genPartID)[q])>=1 && ((*t_genPartMom1ID)[q]==(*t_genPartID)[i] || (*t_genPartMom2ID)[q]==(*t_genPartID)[i])){
			if(qj==99)qj=q;
			else qpj=q;
		      }
		    }
		    if(qj!=99 && qpj!=99){
		      phi1=(*t_jetAK8P4)[Vidx].Phi();
		      eta1=(*t_jetAK8P4)[Vidx].Eta();
		      phi2=(*t_genPartP4)[qj].Phi();
		      eta2=(*t_genPartP4)[qj].Eta();
		      dR=DeltaR(eta1,phi1,eta2,phi2);
		      if(dR <0.4){
			phi2=(*t_genPartP4)[qpj].Phi();
			eta2=(*t_genPartP4)[qpj].Eta();
			dr1=DeltaR(eta1,phi1,eta2,phi2);
			if(dr1<0.4) h_WMass->Fill((*t_jetAK8_MassPruned)[Vidx]);
		      }
		    }
		  }
		    
		}
	      }
	      //==================================================================================================//
	      h_tvsW->Fill(NrecoTop,NrecoV);
	      h_NrecoH_bv->Fill(NrecoH);
	      h_tvsH->Fill(NrecoTop,NrecoH);
	      h_Topvsb->Fill(NrecoTop, cleanBJetL.size());
	      h_Vvsb->Fill(NrecoV, cleanBJetL.size());
	      h_TopvsAK4->Fill(NrecoTop, cleanJetAK4.size());
	      h_VvsAK4->Fill(NrecoV, cleanJetAK4.size());
	      h_tvsNonb->Fill(NrecoTop,cleanNonBJet.size());
	      h_WvsNonb->Fill(NrecoV,cleanNonBJet.size());
	      h_tvsWvsNonb->Fill(NrecoTop,NrecoV,cleanNonBJet.size());

	      dphi=99;jetindex=-99;
	    }//if(!bveto)  
	     

	  }//if(ht>500)
      }
	
    }
  
    
    } //if(signalType..)
}	 
    

	 


  //for matched top jet, plot dr vs pt.
}

//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar  9 11:36:36 2016 by ROOT version 5.34/25
// from TTree tree/tree
// found on file: test_tree.root
//////////////////////////////////////////////////////////

#ifndef TPrimeAna_H
#define TPrimeAna_H

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <cmath>
// Header file for the classes stored in the TTree if any.
#include <vector>
#include "NtupleVariables.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"
#include<string>
#include "TF1.h"
#include "TAxis.h"
// Fixed size dimensions of array or collections stored in the TTree if any.

class TPrimeAna : public NtupleVariables {
public :
  
  TPrimeAna(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data",const char *isData="F");
   virtual ~TPrimeAna();
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  void     EventLoop(const char *, const char *);
  void     BookHistogram(const char *);
  void     ZFinder(int iZ, int type, vector <int> lepid);
  void     LepReco(int e1,int e2,int type, vector <int> lepid);
  vector <TLorentzVector>* getTLV(int part);
  vector <int> sort(unsigned int lepctr,unsigned int *lep);
  vector <int> selGoodAK4Jet(vector <int> elid,vector <int> muid);
  double   HTCleanAK4(vector <int> cleanJetid);
  double   HTCleanAK8(vector <int> cleanJetid);
  vector <float> chi_calculator(int bid, int AK2, int AK3, int AK4,int AK5, int AK6,TLorentzVector v1 );
  double chi_cal(int AK1, int AK2, int AK3, int  AK4, int AK5);
 
  TFile *oFile;
  TFile *file;
  TFile *file1;
  TFile *file2;
  TFile *file3;
 
  vector <TLorentzVector> *t_elecP4;
  vector <TLorentzVector> *t_muonP4;
  vector <TLorentzVector> *t_genPartP4;
  vector <TLorentzVector> *t_jetAK4P4;
  vector <TLorentzVector> *t_jetAK8P4;
  vector <TLorentzVector> *t_metP4;
  vector <TLorentzVector> *t_jetTopJetP4;
  vector <TLorentzVector> *t_jetWJetP4;
  vector <TLorentzVector> *t_ZllP4;
//define histograms
  TH1D *h_Cutflow;
  TH1D *h_rawSize;
  TH1D *h_genDPhi_Ztop;
  TH1D *h_genTMass;
  TH1D *h_genEle1Pt;
  TH1D *h_genEle2Pt;
  TH1D *h_genMatchedEle1Pt;
  TH1D *h_genMatchedEle2Pt;
  TH1D *h_genMu1Pt;
  TH1D *h_genMu2Pt;
  TH1D *h_genMatchedMu1Pt;
  TH1D *h_genMatchedMu2Pt;
  TH1D *h_genZElPt;
  TH1D *h_genZMuPt;
  TH1D *h_genMZElPt;
  TH1D *h_genMZMuPt;

  //======== gen b======//
  TH1D *h_Nb;
  TH1D *h_LeadingbPt;
  //======== reco b======//
  TH1D *h_Nrecob;
  TH1D *h_LeadingrecobPt;
  //========gen top=====//
  TH1D *h_TopPt;
  TH1D *h_LeadingTopPt;
  TH1D *h_Ntop;
  //=======================//
  TH1D *h_recoTopPt;
  TH1D *h_NrecoTop;
  TH1D *h_TopPtdr;
  TH1D *h_Topeta;
  TH1D *h_TopPteta;
  TH1D *h_TopPtPrumass;
  TH1D *h_TopPtSDMass;
  TH1D *h_TopPtcut;
  TH1D *h_JetsAK8Pt;
  TH1D *h_JetsAK8E;
  TH1D *h_JetsAK8Phi;
  TH1D *h_JetsAK8Eta;
  TH1D *h_JetsAK8HT;
  TH1D *h_JetsAK4HT;
  TH1D *h_NJetsAK8;
  TH1D *h_LeadingAK8Pt;
  TH1D *h_JetsAK4Pt;
  TH1D *h_JetsAK4E;
  TH1D *h_JetsAK4Phi;
  TH1D *h_JetsAK4Eta;
  TH1D *h_NJetsAK4;
  TH1D *h_LeadingAK4Pt;
  TH1D *h_4jetEvent;
  
  TH1D *h_deltaRMin;
  TH2D *h2_AK8vstopPt;
  TH1D *h_JetsAK8Ptdr;
  TH1D *h_AK8toptau1;
  TH1D *h_AK8toptau2;
  TH1D *h_AK8toptau3;
  TH1D *h_AK8toptau2bytau1;
  TH1D *h_AK8toptau3bytau2;
  TH1D *h_AK8prunedmass;
  TH1D *h_AK8SoftDropMass;
  TH1D *h_AK8NSubjets;
  TH2D *h2_AK8toptau1vs2;
  TH2D *h2_AK8toptau2vs3;

  //========== w ===========//
  TH1D *h_recoWPt;
  TH1D *h_NrecoW;
  TH1D *h_WdeltaRMin;
  TH1D *h_WPt;
  TH1D *h_WPtdr;
  TH1D *h_WPtcut;
   TH1D *h_AK8Wtau1;
  TH1D *h_AK8Wtau2;
  TH1D *h_AK8Wtau3;
  TH1D *h_AK8Wtau2bytau1;
  TH1D *h_AK8Wtau3bytau2;
  TH1D *h_AK8Wprunedmass;
  TH1D *h_AK8WSoftDropMass;
  TH1D *h_AK8WNSubjets;
  TH2D *h2_AK8Wtau1vs2;
  TH2D *h2_AK8Wtau2vs3;

//======== H =============//
  TH1D *h_recoHPt;
  TH1D *h_NrecoH;
  
  //====== reco leptons =====//
 
  
  //====== reco leptons with Z mass cut =====//
  
  TH1D *h_ele_leadingLepPt_Zcut;
  TH1D *h_ele_leadingLep2Pt_Zcut;
  TH1D *h_ele_leadingLepPhi_Zcut;
  TH1D *h_ele_leadingLep2Phi_Zcut;
  TH1D *h_ele_leadingLepEta_Zcut;
  TH1D *h_ele_leadingLep2Eta_Zcut;
  TH1D *h_ele_dR_Zcut;
  TH1D *h_ele_dPhi_Zcut;
  TH1D *h_ele_mass_Zcut;
  TH1D *h_ele_st_Zcut;
  TH2D *h_ele_dRvsdPhiZcut;
  
  TH1D *h_mu_leadingLepPt_Zcut;
  TH1D *h_mu_leadingLep2Pt_Zcut;
  TH1D *h_mu_leadingLepPhi_Zcut;
  TH1D *h_mu_leadingLep2Phi_Zcut;
  TH1D *h_mu_leadingLepEta_Zcut;
  TH1D *h_mu_leadingLep2Eta_Zcut;
  TH1D *h_mu_dR_Zcut;
  TH1D *h_mu_dPhi_Zcut;
  TH1D *h_mu_mass_Zcut;
  TH1D *h_mu_st_Zcut;
  TH2D *h_mu_dRvsdPhiZcut;
  TH1D *h_dimupt;
  TH1D *h_dielept;

  TH1D *h_Nlep;
  TH1D *h_mu1_iso;
  TH1D *h_mu2_iso;
  TH1D *h_ele1_iso;
  TH1D *h_ele2_iso;
  //=========== reco leptons with Z mass cut and DPhi cut=========//
  TH1D *h_ele_leadingLepPt_dPhicut;
  TH1D *h_ele_leadingLep2Pt_dPhicut;
  TH1D *h_ele_leadingLepPhi_dPhicut;
  TH1D *h_ele_leadingLep2Phi_dPhicut;
  TH1D *h_ele_leadingLepEta_dPhicut;
  TH1D *h_ele_leadingLep2Eta_dPhicut;
  TH1D *h_ele_dR_dPhicut;
  TH1D *h_ele_dPhi_dPhicut;
  TH1D *h_ele_mass_dPhicut;
  TH1D *h_ele_st_dPhicut;
  TH1D *h_mu_leadingLepPt_dPhicut;
  TH1D *h_mu_leadingLep2Pt_dPhicut;
  TH1D *h_mu_leadingLepPhi_dPhicut;
  TH1D *h_mu_leadingLep2Phi_dPhicut;
  TH1D *h_mu_leadingLepEta_dPhicut;
  TH1D *h_mu_leadingLep2Eta_dPhicut;
  TH1D *h_mu_dR_dPhicut;
  TH1D *h_mu_dPhi_dPhicut;
  TH1D *h_mu_mass_dPhicut;
  TH1D *h_mu_st_dPhicut;

  //=========== reco leptons with Z mass cut and Dr cut=========//
  TH1D *h_ele_leadingLepPt_dRcut;
  TH1D *h_ele_leadingLep2Pt_dRcut;
  TH1D *h_ele_leadingLepPhi_dRcut;
  TH1D *h_ele_leadingLep2Phi_dRcut;
  TH1D *h_ele_leadingLepEta_dRcut;
  TH1D *h_ele_leadingLep2Eta_dRcut;
  TH1D *h_ele_dR_dRcut;
  TH1D *h_ele_dPhi_dRcut;
  TH1D *h_ele_mass_dRcut;
  TH1D *h_ele_st_dRcut;
  TH1D *h_mu_leadingLepPt_dRcut;
  TH1D *h_mu_leadingLep2Pt_dRcut;
  TH1D *h_mu_leadingLepPhi_dRcut;
  TH1D *h_mu_leadingLep2Phi_dRcut;
  TH1D *h_mu_leadingLepEta_dRcut;
  TH1D *h_mu_leadingLep2Eta_dRcut;
  TH1D *h_mu_dR_dRcut;
  TH1D *h_mu_dPhi_dRcut;
  TH1D *h_mu_mass_dRcut;
  TH1D *h_mu_st_dRcut;
  //====== gen leptons=========//
   TH1D *h_ZPt;
  TH1D *h_genEle[40];
  TH1D *h_genMu[40];
  TH2D *h_genMu_dRvsdPhiZcut;
  TH2D *h_genEle_dRvsdPhiZcut;
  TH2D *h_genMu_dRvsdPhi;
  TH2D *h_genEle_dRvsdPhi;
  TH1D *h_genEle_diPt;
  TH1D *h_genMu_diPt;

  //====== met ===============//

  TH1D *h_met;
  TH1D *h_metDPhi;

  //============HT/ST check=======//
  TH1D *h_HTcheck;
  TH1D *h_STcheck;
  TH1D *h_HTAK8check;
  TH1D *h_STAK8check;


  //======= Data vs MC ===========//
  TH1D *h_Cutflow_DvsMC;
  TH1D *h_checkTrigNoDZ;
  TH1D *h_checkTrigDZ;
  TH1D *h_l1Pt;
  TH1D *h_l1Eta;
  TH1D *h_l1Phi;
  TH1D *h_l1Pt_Up;
  TH1D *h_l1Eta_Up;
  TH1D *h_l1Phi_Up;
  TH1D *h_l1Pt_Down;
  TH1D *h_l1Eta_Down;
  TH1D *h_l1Phi_Down;
	  
  TH1D *h_l2Pt;
  TH1D *h_l2Eta;
  TH1D *h_l2Phi;
  TH1D *h_l2Pt_Up;
  TH1D *h_l2Eta_Up;
  TH1D *h_l2Phi_Up;
  TH1D *h_l2Pt_Down;
  TH1D *h_l2Eta_Down;
  TH1D *h_l2Phi_Down;

  TH1D *h_l1l2_dR;
  TH1D *h_l1l2_dPhi;
  TH1D *h_l1l2_dR_Up;
  TH1D *h_l1l2_dPhi_Up;
  TH1D *h_l1l2_dR_Down;
  TH1D *h_l1l2_dPhi_Down;

  TH1D *h_dilepPt;
  TH1D *h_dilepEta;
  TH1D *h_dilepPhi;
  TH1D *h_dilepMass;
  TH1D *h_dilepMass_NoCut;
  TH1D *h_ST;
  TH1D *h_npv;
  TH1D *h_dilepPt_Up;
  TH1D *h_dilepEta_Up;
  TH1D *h_dilepPhi_Up;
  TH1D *h_dilepMass_Up;
  TH1D *h_dilepMass_NoCut_Up;
  TH1D *h_ST_Up;
  TH1D *h_npv_Up;
  TH1D *h_dilepPt_Down;
  TH1D *h_dilepEta_Down;
  TH1D *h_dilepPhi_Down;
  TH1D *h_dilepMass_Down;
  TH1D *h_dilepMass_NoCut_Down;
  TH1D *h_ST_Down;
  TH1D *h_npv_Down;

  TH1D *h_NbjetL;
  TH1D *h_NbjetM;
  TH1D *h_NAK4;
  TH1D *h_NAK8;
  TH1D *h_NbjetL_Up;
  TH1D *h_NbjetM_Up;
  TH1D *h_NAK4_Up;
  TH1D *h_NAK8_Up;
  TH1D *h_NbjetL_Down;
  TH1D *h_NbjetM_Down;
  TH1D *h_NAK4_Down;
  TH1D *h_NAK8_Down;
  
  TH1D *h_AK8Pt;
  TH1D *h_AK8Phi;
  TH1D *h_AK8Eta;
  TH1D *h_AK8HT;
  TH1D *h_LeadAK8Pt;
  TH1D *h_LeadAK8Phi;
  TH1D *h_LeadAK8Eta;
  TH1D *h_AK8Pt_Up;
  TH1D *h_AK8Phi_Up;
  TH1D *h_AK8Eta_Up;
  TH1D *h_AK8HT_Up;
  TH1D *h_LeadAK8Pt_Up;
  TH1D *h_LeadAK8Phi_Up;
  TH1D *h_LeadAK8Eta_Up;
  TH1D *h_AK8Pt_Down;
  TH1D *h_AK8Phi_Down;
  TH1D *h_AK8Eta_Down;
  TH1D *h_AK8HT_Down;
  TH1D *h_LeadAK8Pt_Down;
  TH1D *h_LeadAK8Phi_Down;
  TH1D *h_LeadAK8Eta_Down;
  
  TH1D *h_AK4Pt;
  TH1D *h_AK4Phi;
  TH1D *h_AK4Eta;
  TH1D *h_AK4HT;
  TH1D *h_LeadAK4Pt;
  TH1D *h_LeadAK4Phi;
  TH1D *h_LeadAK4Eta;
  TH1D *h_AK4Pt_Up;
  TH1D *h_AK4Phi_Up;
  TH1D *h_AK4Eta_Up;
  TH1D *h_AK4HT_Up;
  TH1D *h_LeadAK4Pt_Up;
  TH1D *h_LeadAK4Phi_Up;
  TH1D *h_LeadAK4Eta_Up;
  TH1D *h_AK4Pt_Down;
  TH1D *h_AK4Phi_Down;
  TH1D *h_AK4Eta_Down;
  TH1D *h_AK4HT_Down;
  TH1D *h_LeadAK4Pt_Down;
  TH1D *h_LeadAK4Phi_Down;
  TH1D *h_LeadAK4Eta_Down;

  TH1D *h_Mu_Chi2;
  TH1D *h_Mu_Chi2_Up;
  TH1D *h_Mu_Chi2_Down;
  TH1D *h_Mu_Chi2_1t;
  TH1D *h_Mu_Chi2_1t_Up;
  TH1D *h_Mu_Chi2_1t_Down;
  TH1D *h_Mu_Chi2_0t1V;
  TH1D *h_Mu_Chi2_0t1V_Up;
  TH1D *h_Mu_Chi2_0t1V_Down;
  TH1D *h_Mu_Chi2_0t2V;
  TH1D *h_Mu_Chi2_0t2V_Up;
  TH1D *h_Mu_Chi2_0t2V_Down;
  TH1D *h_Mu_Chi2_1t1V;
  TH1D *h_Mu_Chi2_1t1V_Up;
  TH1D *h_Mu_Chi2_1t1V_Down;
  TH1D *h_Mu_Chi2_1t0V;
  TH1D *h_Mu_Chi2_1t0V_Up;
  TH1D *h_Mu_Chi2_1t0V_Down;
  TH1D *h_Mu_TMass_chi2_lep;
  TH1D *h_Mu_TMass_chi2_lep_Up;
  TH1D *h_Mu_TMass_chi2_lep_Down;
  TH1D *h_Mu_TMass_chi2_had;
  TH1D *h_Mu_TMass_chi2_had_Up;
  TH1D *h_Mu_TMass_chi2_had_Down;
  TH1D *h_Mu_TMass_chi2;
  TH1D *h_Mu_TMass_chi2_Up;
  TH1D *h_Mu_TMass_chi2_Down;
  
  TH1D *h_Mu_Chi2_cut;
  TH1D *h_Mu_Chi2_cut_Up;
  TH1D *h_Mu_Chi2_cut_Down;
  TH1D *h_Mu_Chi2_1t_cut;
  TH1D *h_Mu_Chi2_1t_cut_Up;
  TH1D *h_Mu_Chi2_1t_cut_Down;
  TH1D *h_Mu_Chi2_0t1V_cut;
  TH1D *h_Mu_Chi2_0t1V_cut_Up;
  TH1D *h_Mu_Chi2_0t1V_cut_Down;
  TH1D *h_Mu_Chi2_0t2V_cut;
  TH1D *h_Mu_Chi2_0t2V_cut_Up;
  TH1D *h_Mu_Chi2_0t2V_cut_Down;
  TH1D *h_Mu_Chi2_1t1V_cut;
  TH1D *h_Mu_Chi2_1t1V_cut_Up;
  TH1D *h_Mu_Chi2_1t1V_cut_Down;
  TH1D *h_Mu_Chi2_1t0V_cut;
  TH1D *h_Mu_Chi2_1t0V_cut_Up;
  TH1D *h_Mu_Chi2_1t0V_cut_Down;
  TH1D *h_Mu_TMass_chi2_lep_cut;
  TH1D *h_Mu_TMass_chi2_lep_cut_Up;
  TH1D *h_Mu_TMass_chi2_lep_cut_Down;
  TH1D *h_Mu_TMass_chi2_had_cut;
  TH1D *h_Mu_TMass_chi2_had_cut_Up;
  TH1D *h_Mu_TMass_chi2_had_cut_Down;
  TH1D *h_Mu_TMass_chi2_cut;
  TH1D *h_Mu_TMass_chi2_cut_Up;
  TH1D *h_Mu_TMass_chi2_cut_Down;

  //x=x=x=x=x=x=x=x=x=x=x=//
  TH1D *h_Ele_Cutflow_DvsMC;
  TH1D *h_Ele_checkTrigNoDZ;
  TH1D *h_Ele_checkTrigDZ;
  TH1D *h_Ele_l1Pt;
  TH1D *h_Ele_l1Eta;
  TH1D *h_Ele_l1Phi;
	  
  TH1D *h_Ele_l2Pt;
  TH1D *h_Ele_l2Eta;
  TH1D *h_Ele_l2Phi;

  TH1D *h_Ele_l1l2_dR;
  TH1D *h_Ele_l1l2_dPhi;

  TH1D *h_Ele_dilepPt;
  TH1D *h_Ele_dilepEta;
  TH1D *h_Ele_dilepPhi;
  TH1D *h_Ele_dilepMass;
  TH1D *h_Ele_dilepMass_NoCut;
  TH1D *h_Ele_ST;
  TH1D *h_Ele_npv;
  
  TH1D *h_Ele_NbjetL;
  TH1D *h_Ele_NbjetM;
  TH1D *h_Ele_NAK4;
  TH1D *h_Ele_NAK8;
  
  TH1D *h_Ele_AK8Pt;
  TH1D *h_Ele_AK8Phi;
  TH1D *h_Ele_AK8Eta;
  TH1D *h_Ele_AK8HT;
  TH1D *h_Ele_LeadAK8Pt;
  TH1D *h_Ele_LeadAK8Phi;
  TH1D *h_Ele_LeadAK8Eta;
  
  TH1D *h_Ele_AK4Pt;
  TH1D *h_Ele_AK4Phi;
  TH1D *h_Ele_AK4Eta;
  TH1D *h_Ele_AK4HT;
  TH1D *h_Ele_LeadAK4Pt;
  TH1D *h_Ele_LeadAK4Phi;
  TH1D *h_Ele_LeadAK4Eta;

  TH1D *h_Ele_l1Pt_Up;
  TH1D *h_Ele_l1Eta_Up;
  TH1D *h_Ele_l1Phi_Up;
	  
  TH1D *h_Ele_l2Pt_Up;
  TH1D *h_Ele_l2Eta_Up;
  TH1D *h_Ele_l2Phi_Up;

  TH1D *h_Ele_l1l2_dR_Up;
  TH1D *h_Ele_l1l2_dPhi_Up;

  TH1D *h_Ele_dilepPt_Up;
  TH1D *h_Ele_dilepEta_Up;
  TH1D *h_Ele_dilepPhi_Up;
  TH1D *h_Ele_dilepMass_Up;
  TH1D *h_Ele_dilepMass_NoCut_Up;
  TH1D *h_Ele_ST_Up;
  TH1D *h_Ele_npv_Up;
  
  TH1D *h_Ele_NbjetL_Up;
  TH1D *h_Ele_NbjetM_Up;
  TH1D *h_Ele_NAK4_Up;
  TH1D *h_Ele_NAK8_Up;
  
  TH1D *h_Ele_AK8Pt_Up;
  TH1D *h_Ele_AK8Phi_Up;
  TH1D *h_Ele_AK8Eta_Up;
  TH1D *h_Ele_AK8HT_Up;
  TH1D *h_Ele_LeadAK8Pt_Up;
  TH1D *h_Ele_LeadAK8Phi_Up;
  TH1D *h_Ele_LeadAK8Eta_Up;
  
  TH1D *h_Ele_AK4Pt_Up;
  TH1D *h_Ele_AK4Phi_Up;
  TH1D *h_Ele_AK4Eta_Up;
  TH1D *h_Ele_AK4HT_Up;
  TH1D *h_Ele_LeadAK4Pt_Up;
  TH1D *h_Ele_LeadAK4Phi_Up;
  TH1D *h_Ele_LeadAK4Eta_Up;

  TH1D *h_Ele_l1Pt_Down;
  TH1D *h_Ele_l1Eta_Down;
  TH1D *h_Ele_l1Phi_Down;
	  
  TH1D *h_Ele_l2Pt_Down;
  TH1D *h_Ele_l2Eta_Down;
  TH1D *h_Ele_l2Phi_Down;

  TH1D *h_Ele_l1l2_dR_Down;
  TH1D *h_Ele_l1l2_dPhi_Down;

  TH1D *h_Ele_dilepPt_Down;
  TH1D *h_Ele_dilepEta_Down;
  TH1D *h_Ele_dilepPhi_Down;
  TH1D *h_Ele_dilepMass_Down;
  TH1D *h_Ele_dilepMass_NoCut_Down;
  TH1D *h_Ele_ST_Down;
  TH1D *h_Ele_npv_Down;
  
  TH1D *h_Ele_NbjetL_Down;
  TH1D *h_Ele_NbjetM_Down;
  TH1D *h_Ele_NAK4_Down;
  TH1D *h_Ele_NAK8_Down;
  
  TH1D *h_Ele_AK8Pt_Down;
  TH1D *h_Ele_AK8Phi_Down;
  TH1D *h_Ele_AK8Eta_Down;
  TH1D *h_Ele_AK8HT_Down;
  TH1D *h_Ele_LeadAK8Pt_Down;
  TH1D *h_Ele_LeadAK8Phi_Down;
  TH1D *h_Ele_LeadAK8Eta_Down;
  
  TH1D *h_Ele_AK4Pt_Down;
  TH1D *h_Ele_AK4Phi_Down;
  TH1D *h_Ele_AK4Eta_Down;
  TH1D *h_Ele_AK4HT_Down;
  TH1D *h_Ele_LeadAK4Pt_Down;
  TH1D *h_Ele_LeadAK4Phi_Down;
  TH1D *h_Ele_LeadAK4Eta_Down;

  TH1D *h_Ele_Chi2;
  TH1D *h_Ele_Chi2_Up;
  TH1D *h_Ele_Chi2_Down;
  TH1D *h_Ele_Chi2_1t;
  TH1D *h_Ele_Chi2_1t_Up;
  TH1D *h_Ele_Chi2_1t_Down;
  TH1D *h_Ele_Chi2_0t1V;
  TH1D *h_Ele_Chi2_0t1V_Up;
  TH1D *h_Ele_Chi2_0t1V_Down;
  TH1D *h_Ele_Chi2_0t2V;
  TH1D *h_Ele_Chi2_0t2V_Up;
  TH1D *h_Ele_Chi2_0t2V_Down;
  TH1D *h_Ele_Chi2_1t1V;
  TH1D *h_Ele_Chi2_1t1V_Up;
  TH1D *h_Ele_Chi2_1t1V_Down;
  TH1D *h_Ele_Chi2_1t0V;
  TH1D *h_Ele_Chi2_1t0V_Up;
  TH1D *h_Ele_Chi2_1t0V_Down;
  TH1D *h_Ele_TMass_chi2_lep;
  TH1D *h_Ele_TMass_chi2_lep_Up;
  TH1D *h_Ele_TMass_chi2_lep_Down;
  TH1D *h_Ele_TMass_chi2_had;
  TH1D *h_Ele_TMass_chi2_had_Up;
  TH1D *h_Ele_TMass_chi2_had_Down;
  TH1D *h_Ele_TMass_chi2;
  TH1D *h_Ele_TMass_chi2_Up;
  TH1D *h_Ele_TMass_chi2_Down;

  TH1D *h_Ele_Chi2_cut;
  TH1D *h_Ele_Chi2_cut_Up;
  TH1D *h_Ele_Chi2_cut_Down;
  TH1D *h_Ele_Chi2_1t_cut;
  TH1D *h_Ele_Chi2_1t_cut_Up;
  TH1D *h_Ele_Chi2_1t_cut_Down;
  TH1D *h_Ele_Chi2_0t1V_cut;
  TH1D *h_Ele_Chi2_0t1V_cut_Up;
  TH1D *h_Ele_Chi2_0t1V_cut_Down;
  TH1D *h_Ele_Chi2_0t2V_cut;
  TH1D *h_Ele_Chi2_0t2V_cut_Up;
  TH1D *h_Ele_Chi2_0t2V_cut_Down;
  TH1D *h_Ele_Chi2_1t1V_cut;
  TH1D *h_Ele_Chi2_1t1V_cut_Up;
  TH1D *h_Ele_Chi2_1t1V_cut_Down;
  TH1D *h_Ele_Chi2_1t0V_cut;
  TH1D *h_Ele_Chi2_1t0V_cut_Up;
  TH1D *h_Ele_Chi2_1t0V_cut_Down;
  TH1D *h_Ele_TMass_chi2_lep_cut;
  TH1D *h_Ele_TMass_chi2_lep_cut_Up;
  TH1D *h_Ele_TMass_chi2_lep_cut_Down;
  TH1D *h_Ele_TMass_chi2_had_cut;
  TH1D *h_Ele_TMass_chi2_had_cut_Up;
  TH1D *h_Ele_TMass_chi2_had_cut_Down;
  TH1D *h_Ele_TMass_chi2_cut;
  TH1D *h_Ele_TMass_chi2_cut_Up;
  TH1D *h_Ele_TMass_chi2_cut_Down;

  //===== HT sel ============//
  TH1D *h_dimupt_HT;
  TH1D *h_dielept_HT;
  TH1D *h_Nlep_HT;
  

  TH1D *h_JetsAK8Pt_HT;
  TH1D *h_JetsAK8HT_HT;
  TH1D *h_JetsAK4HT_HT;
  TH1D *h_NJetsAK8_HT;
  TH1D *h_LeadingAK8Pt_HT;
  TH1D *h_JetsAK4Pt_HT;
  TH1D *h_NJetsAK4_HT;
  TH1D *h_LeadingAK4Pt_HT;
  TH1D *h_4jetEvent_HT;

  TH1D *h_NrecobT_HT;
  TH1D *h_Nrecob_HT;
  TH1D *h_NrecobL_HT;

  TH1D *h_bL;
  TH1D *h_bM;

  TH1D *h_NrecoTop_HT;
  TH2D *h_bmatrix;

  //====== B veto ===========//
  TH1D *h_ele_leadingLepPt;
  TH1D *h_ele_leadingLep2Pt;
  TH1D *h_ele_Lep2by1Pt;
  TH1D *h_ele_Lep2Ptfrac;
  TH2D *h_ele_Lep2vs1Pt;
  TH1D *h_mu_leadingLepPt;
  TH1D *h_mu_leadingLep2Pt;
  TH1D *h_mu_Lep2by1Pt;
  TH1D *h_mu_Lep2Ptfrac;
  TH2D *h_mu_Lep2vs1Pt;
  TH2D *h_tvsW;
  TH3D *h_tvsWvsNonb;
  TH1D *h_TopWMix;
  TH2D *h_tvsH;
  TH2D *h_tvsNonb;
  TH2D *h_WvsNonb;
  TH1D *h_NrecoH_bv;
  TH1D *h_1W_Dilep_dPhi;
  TH1D *h_LW_Dilep_dPhi;
  TH1D *h_1Top_Dilep_dPhi;
  TH1D *h_LTop_Dilep_dPhi;
  TH1D *h_recoTMass1b;
  TH1D *h_recoTMass2b;

  TH1D *h_recoTMass_2top;
  TH1D *h_recoTMass_1top1W;
  TH1D *h_recoTMass_2W;

  TH1D *h_TMass;
  TH1D *h_TMass_1W;
  TH1D *h_TMass_2W;
  TH1D *h_topMass;
  TH1D *h_topMassW;
  TH1D *h_resTop;
  TH1D *h_WMass;
  TH1D *h_ZMass;
  TH1D *h_genV;
  TH1D *h_gentL;
  TH1D *h_genV2;
  TH1D *h_AK4_Wb;
  TH1D *h_topTerm;
  TH1D *h_WbTerm;
  TH1D *h_AK4WbTerm;
  TH1D *h_AK4WTerm;
  TH1D *h_ZTerm;
  TH1D *h_Chi2;
  TH1D *h_Chi2_1t;
  TH1D *h_Chi2_0t1V;
  TH1D *h_Chi2_0t2V;
  TH1D *h_Chi2_1t1V;
  TH1D *h_Chi2_1t0V;
  TH1D *h_TMass_chi2_lep;
  TH1D *h_TMass_chi2_had;
  TH1D *h_TMass_chi2;

  TH1D *h_Chi2_cut;
  TH1D *h_Chi2_1t_cut;
  TH1D *h_Chi2_0t1V_cut;
  TH1D *h_Chi2_0t2V_cut;
  TH1D *h_Chi2_1t1V_cut;
  TH1D *h_Chi2_1t0V_cut;
  TH1D *h_TMass_chi2_lep_cut;
  TH1D *h_TMass_chi2_had_cut;
  TH1D *h_TMass_chi2_cut;


  TH1D *h_UB_Mu_Chi2;
  TH1D *h_UB_Mu_Chi2_Up;
  TH1D *h_UB_Mu_Chi2_Down;
  TH1D *h_UB_Mu_Chi2_1t;
  TH1D *h_UB_Mu_Chi2_1t_Up;
  TH1D *h_UB_Mu_Chi2_1t_Down;
  TH1D *h_UB_Mu_Chi2_0t1V;
  TH1D *h_UB_Mu_Chi2_0t1V_Up;
  TH1D *h_UB_Mu_Chi2_0t1V_Down;
  TH1D *h_UB_Mu_Chi2_0t2V;
  TH1D *h_UB_Mu_Chi2_0t2V_Up;
  TH1D *h_UB_Mu_Chi2_0t2V_Down;
  TH1D *h_UB_Mu_Chi2_1t1V;
  TH1D *h_UB_Mu_Chi2_1t1V_Up;
  TH1D *h_UB_Mu_Chi2_1t1V_Down;
  TH1D *h_UB_Mu_Chi2_1t0V;
  TH1D *h_UB_Mu_Chi2_1t0V_Up;
  TH1D *h_UB_Mu_Chi2_1t0V_Down;
  TH1D *h_UB_Mu_TMass_chi2_lep;
  TH1D *h_UB_Mu_TMass_chi2_lep_Up;
  TH1D *h_UB_Mu_TMass_chi2_lep_Down;
  TH1D *h_UB_Mu_TMass_chi2_had;
  TH1D *h_UB_Mu_TMass_chi2_had_Up;
  TH1D *h_UB_Mu_TMass_chi2_had_Down;
  TH1D *h_UB_Mu_TMass_chi2;
  TH1D *h_UB_Mu_TMass_chi2_Up;
  TH1D *h_UB_Mu_TMass_chi2_Down;
  
  TH1D *h_UB_Mu_Chi2_cut;
  TH1D *h_UB_Mu_Chi2_cut_Up;
  TH1D *h_UB_Mu_Chi2_cut_Down;
  TH1D *h_UB_Mu_Chi2_1t_cut;
  TH1D *h_UB_Mu_Chi2_1t_cut_Up;
  TH1D *h_UB_Mu_Chi2_1t_cut_Down;
  TH1D *h_UB_Mu_Chi2_0t1V_cut;
  TH1D *h_UB_Mu_Chi2_0t1V_cut_Up;
  TH1D *h_UB_Mu_Chi2_0t1V_cut_Down;
  TH1D *h_UB_Mu_Chi2_0t2V_cut;
  TH1D *h_UB_Mu_Chi2_0t2V_cut_Up;
  TH1D *h_UB_Mu_Chi2_0t2V_cut_Down;
  TH1D *h_UB_Mu_Chi2_1t1V_cut;
  TH1D *h_UB_Mu_Chi2_1t1V_cut_Up;
  TH1D *h_UB_Mu_Chi2_1t1V_cut_Down;
  TH1D *h_UB_Mu_Chi2_1t0V_cut;
  TH1D *h_UB_Mu_Chi2_1t0V_cut_Up;
  TH1D *h_UB_Mu_Chi2_1t0V_cut_Down;
  TH1D *h_UB_Mu_TMass_chi2_lep_cut;
  TH1D *h_UB_Mu_TMass_chi2_lep_cut_Up;
  TH1D *h_UB_Mu_TMass_chi2_lep_cut_Down;
  TH1D *h_UB_Mu_TMass_chi2_had_cut;
  TH1D *h_UB_Mu_TMass_chi2_had_cut_Up;
  TH1D *h_UB_Mu_TMass_chi2_had_cut_Down;
  TH1D *h_UB_Mu_TMass_chi2_cut;
  TH1D *h_UB_Mu_TMass_chi2_cut_Up;
  TH1D *h_UB_Mu_TMass_chi2_cut_Down;

  //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
  TH1D *h_UB_Ele_Chi2;
  TH1D *h_UB_Ele_Chi2_Up;
  TH1D *h_UB_Ele_Chi2_Down;
  TH1D *h_UB_Ele_Chi2_1t;
  TH1D *h_UB_Ele_Chi2_1t_Up;
  TH1D *h_UB_Ele_Chi2_1t_Down;
  TH1D *h_UB_Ele_Chi2_0t1V;
  TH1D *h_UB_Ele_Chi2_0t1V_Up;
  TH1D *h_UB_Ele_Chi2_0t1V_Down;
  TH1D *h_UB_Ele_Chi2_0t2V;
  TH1D *h_UB_Ele_Chi2_0t2V_Up;
  TH1D *h_UB_Ele_Chi2_0t2V_Down;
  TH1D *h_UB_Ele_Chi2_1t1V;
  TH1D *h_UB_Ele_Chi2_1t1V_Up;
  TH1D *h_UB_Ele_Chi2_1t1V_Down;
  TH1D *h_UB_Ele_Chi2_1t0V;
  TH1D *h_UB_Ele_Chi2_1t0V_Up;
  TH1D *h_UB_Ele_Chi2_1t0V_Down;
  TH1D *h_UB_Ele_TMass_chi2_lep;
  TH1D *h_UB_Ele_TMass_chi2_lep_Up;
  TH1D *h_UB_Ele_TMass_chi2_lep_Down;
  TH1D *h_UB_Ele_TMass_chi2_had;
  TH1D *h_UB_Ele_TMass_chi2_had_Up;
  TH1D *h_UB_Ele_TMass_chi2_had_Down;
  TH1D *h_UB_Ele_TMass_chi2;
  TH1D *h_UB_Ele_TMass_chi2_Up;
  TH1D *h_UB_Ele_TMass_chi2_Down;

  TH1D *h_UB_Ele_Chi2_cut;
  TH1D *h_UB_Ele_Chi2_cut_Up;
  TH1D *h_UB_Ele_Chi2_cut_Down;
  TH1D *h_UB_Ele_Chi2_1t_cut;
  TH1D *h_UB_Ele_Chi2_1t_cut_Up;
  TH1D *h_UB_Ele_Chi2_1t_cut_Down;
  TH1D *h_UB_Ele_Chi2_0t1V_cut;
  TH1D *h_UB_Ele_Chi2_0t1V_cut_Up;
  TH1D *h_UB_Ele_Chi2_0t1V_cut_Down;
  TH1D *h_UB_Ele_Chi2_0t2V_cut;
  TH1D *h_UB_Ele_Chi2_0t2V_cut_Up;
  TH1D *h_UB_Ele_Chi2_0t2V_cut_Down;
  TH1D *h_UB_Ele_Chi2_1t1V_cut;
  TH1D *h_UB_Ele_Chi2_1t1V_cut_Up;
  TH1D *h_UB_Ele_Chi2_1t1V_cut_Down;
  TH1D *h_UB_Ele_Chi2_1t0V_cut;
  TH1D *h_UB_Ele_Chi2_1t0V_cut_Up;
  TH1D *h_UB_Ele_Chi2_1t0V_cut_Down;
  TH1D *h_UB_Ele_TMass_chi2_lep_cut;
  TH1D *h_UB_Ele_TMass_chi2_lep_cut_Up;
  TH1D *h_UB_Ele_TMass_chi2_lep_cut_Down;
  TH1D *h_UB_Ele_TMass_chi2_had_cut;
  TH1D *h_UB_Ele_TMass_chi2_had_cut_Up;
  TH1D *h_UB_Ele_TMass_chi2_had_cut_Down;
  TH1D *h_UB_Ele_TMass_chi2_cut;
  TH1D *h_UB_Ele_TMass_chi2_cut_Up;
  TH1D *h_UB_Ele_TMass_chi2_cut_Down;
  
  
  TH1D *h_AKV;
  TH1D *h_mass_diff;


  
  TH1D *h_resZ;
  TH1D *h_Nos;
  TH1D *h_dRZ;
  TH1D *h_Nb2t;
  TH1D *h_NW2t;
  TH1D *h_Jeta;
  TH1D *h_Jpt;
  TH1D *h_Jphi;
  TH1D *h_dR_JZ;
  TH1D *h_nAK8;
  TH1D *h_topTag;
  TH1D *h_Vtag;
  TH2D *h_AK8vsb;
  TH2D *h_Topvsb;
  TH2D *h_Vvsb;
  TH2D *h_TopvsAK4;
  TH2D *h_VvsAK4;
  TH1D *h_Jneweta;
  TH1D *h_Jnewphi;
  TH1D *h_Jnewpt;
  TH1D *h_ratioPt;



  vector <float> GetChi(int NrecoTop,int NrecoV,TLorentzVector v1,vector <int> goodMuid,vector <int> goodElid, vector <int> cleanJetAK4, vector <int> cleanBJetL);

 };

#endif

#ifdef TPrimeAna_cxx

void TPrimeAna::BookHistogram(const char *outFileName) {

//  char hname[200], htit[200];
//  double xlow = 0.0,  xhigh = 2000.0;
//  int nbins = 2000;0

 

  oFile = new TFile(outFileName, "recreate");
  oFile->mkdir("Cutflow");
  oFile->cd("Cutflow");
  h_Cutflow=new TH1D("Cutflow","Number of Events surviving ",23,0.0,23);
  h_rawSize =new TH1D("AK4Njets","Number of AK4 Jets",21,0.0,21);
  
  oFile->mkdir("Jets");
  oFile->cd("Jets");
  h_JetsAK8Pt=new TH1D("JetsAK8Pt","P_{T} for AK8 Jets",50,0.0,2000.);
  h_JetsAK8E=new TH1D("JetsAK8E","E for AK8 Jets",50,100.0,3000.);
  h_JetsAK8Phi=new TH1D("JetsAK8Phi","#phi for AK8 Jets",32,3.2,3.2);
  h_JetsAK8Eta=new TH1D("JetsAK8Eta","#eta for AK8 Jets",40,-4.,4.);
  h_NJetsAK8=new TH1D("AK8Njets","Number of AK8 Jets",21,0.0,21);
  h_JetsAK8HT=new TH1D("JetsAK8HT","H_{T} for AK8 Jets",100,0.0,6000.);
  h_LeadingAK8Pt=new TH1D("LeadingAK8Pt","P_{T} for Leading AK8 Jet",50,0.0,2000.);
  
  h_JetsAK4Pt=new TH1D("JetsAK4Pt","P_{T} for AK4 Jets",50,0.0,2000.);
  h_JetsAK4E=new TH1D("JetsAK4E","E for AK4 Jets",50,30.0,1000.);
  h_JetsAK4Phi=new TH1D("JetsAK4Phi","#phi for AK4 Jets",32,3.2,3.2);
  h_JetsAK4Eta=new TH1D("JetsAK4Eta","#eta for AK4 Jets",40,-4.,4.);
  h_NJetsAK4=new TH1D("AK4Njets","Number of AK4 Jets",21,0.0,21);
  h_JetsAK4HT=new TH1D("JetsAK4HT","H_{T} for AK4 Jets",100,0.0,6000.);
  h_LeadingAK4Pt=new TH1D("LeadingAK4Pt","P_{T} for Leading AK4 Jet",50,0.0,2000.);

  h_4jetEvent=new TH1D("4jetEvent","Number of Events with >=4 Jets and at least one AK8 Jet",21,0.0,21);
  //==========================================================================================//
  oFile->mkdir("gen_quarks");
  oFile->cd("gen_quarks");
  //================== gen b===================================//
  h_Nb=new TH1D("Nb","Number of b quarks",10,0.0,10);
  h_LeadingbPt=new TH1D("LeadingbPt","P_{T} for leading b quark",38,0.0,2000.);
  //================== gen top==================================//
  h_Ntop=new TH1D("Ntop","Number of Tops",10,0.0,10);
  h_TopPt=new TH1D("TopPt","P_{T} for tops",38,100.0,2000.);
  h_LeadingTopPt=new TH1D("LeadingTopPt","P_{T} for leading top",38,0.0,2000.);
 
  //================== reco b===================================//
  oFile->mkdir("reco_bquark");
  oFile->cd("reco_bquark");
  
  h_Nrecob=new TH1D("Nrecob","Number of reconstructed b quarks",10,0.0,10);
  h_LeadingrecobPt=new TH1D("LeadingrecobPt","P_{T} for leading reconstructed  b quark",38,100.0,2000.);


   //============ gen info=========================================//

  oFile->mkdir("gen_info");
  oFile->cd("gen_info");
  h_genDPhi_Ztop=new TH1D("genDPhi_Ztop","#Delta#Phi between top and Z->l^{+}l^{-} ( other V decays hadronically)",32,-3.2,3.2);
  h_genTMass=new TH1D("genTMass","Gen level TPrime Mass",80,0.0,4000.);
  h_genEle1Pt=new TH1D("genEle1Pt","Gen level Leading Electron Pt",50,0.0,1000.);
  
  h_genMatchedEle1Pt=new TH1D("genMatchedEle1Pt","Gen level Leading Matched Electron Pt",50,0.0,1000.);
  h_genMu1Pt=new TH1D("genMu1Pt","Gen level Leading Muon Pt",50,0.0,1000.);
  h_genMatchedMu1Pt=new TH1D("genMatchedMu1Pt","Gen level Leading Matched Muon Pt",50,0.0,1000.);
  h_genEle2Pt=new TH1D("genEle2Pt","Gen level sub leading Electron Pt",50,0.0,1000.);
  
  h_genMatchedEle2Pt=new TH1D("genMatchedEle2Pt","Gen level sub leading Matched Electron Pt",50,0.0,1000.);
  h_genMu2Pt=new TH1D("genMu2Pt","Gen level sub leading Muon Pt",50,0.0,1000.);
  h_genMatchedMu2Pt=new TH1D("genMatchedMu2Pt","Gen level sub leading Matched Muon Pt",50,0.0,1000.);
  h_genZElPt=new TH1D("genZElPt","Gen level Z(ee) Pt",100,0.0,2000.);
  h_genZMuPt=new TH1D("genZMuPt","Gen level Z(#mu#mu) Pt",100,0.0,2000.);
  h_genMZElPt=new TH1D("genMZElPt","Gen level Matched Z(ee) Pt",100,0.0,2000.);
  h_genMZMuPt=new TH1D("genMZMuPt","Gen level Matched Z(#mu#mu) Pt",100,0.0,2000.);
  //=================== reco top================================//
  oFile->mkdir("reco_Top");
  oFile->cd("reco_Top");

  h_NrecoTop=new TH1D("NrecoTop","Number of reconstructed Tops",10,0.0,10);
  h_recoTopPt=new TH1D("recoTopPt","P_{T} for reco tops",50,0.0,2000.);
  h_deltaRMin=new TH1D("DeltaRMin","#delta R between AK8 Jets and top quarks",20,0.,1.0);
  h2_AK8vstopPt=new TH2D("AK8vstopPt","p_{T} of top quarks vs (matched AK8 jet P_{T}/top p_{T})",50,100.0,1200.0,20,0.0,2.0);
  h_JetsAK8Ptdr=new TH1D("JetsAK8Ptdr","P_{T} for AK8 Jets with a dr cut of 0.6",50,100.0,2000.);
  h_TopPtdr=new TH1D("TopPtdr","P_{T} for tops with dr cut of 0.6",38,100.0,2000.);
  h_Topeta=new TH1D("Topeta","#eta for tops after dr<0.6 cut",40,-4.0,4.0);
  h_TopPteta=new TH1D("TopPteta","P_{T} for tops after |#eta|<2.4 cut",38,100.0,2000.);
  h_TopPtcut=new TH1D("TopPtcut","P_{T} for tops with cuts",38,100.0,2000.);
  h_TopPtPrumass=new TH1D("TopPtPruMass","P_{T} for tops with pruned mass cut",38,100.0,2000.);
  h_TopPtSDMass=new TH1D("TopPtSDMass","P_{T} for tops with soft drop mass cut",38,100.0,2000.);
  h_AK8toptau1=new TH1D("AK8toptau1","#tau1 for AK8 Jets matched with Top",20,0.0,1.0);
  h_AK8toptau2=new TH1D("AK8toptau2","#tau2 for AK8 Jets matched with Top",20,0.0,1.0);
  h_AK8toptau3=new TH1D("AK8toptau3","#tau3 for AK8 Jets matched with Top",20,0.0,1.0);
  h_AK8toptau2bytau1=new TH1D("AK8toptau2bytau1","#tau_{2}/#tau_{1} for AK8 Jets matched with Top",20,0.0,1.0);
  h_AK8toptau3bytau2=new TH1D("AK8toptau3bytau2","#tau_{3}/#tau_{2} for AK8 Jets matched with Top",20,0.0,1.0);
  h_AK8prunedmass=new TH1D("AK8prunedmass","Pruned mass for AK8 Jets matched with top",12,0.0,300.);
  h_AK8SoftDropMass=new TH1D("AK8SoftDropMass","Soft mass drop for AK8 Jets matched with top",12,0.0,300.0);
  h_AK8NSubjets=new TH1D("AK8NSubjets","NSubjets for AK8 Jets matched with top",10,0.0,10);
  h2_AK8toptau1vs2=new TH2D("AK8toptau1vs2","#tau_{2} vs #tau_{1} for AK8 Jets matched with Top",16,0.0,0.8,16,0.0,0.8);
  h2_AK8toptau2vs3=new TH2D("AK8toptau2vs3","#tau_{3} vs #tau_{2} for AK8 Jets matched with Top",10,0.0,0.5,10,0.0,0.5);
  //=================== reco Higgs ================================//
  
  oFile->mkdir("reco_H");
  oFile->cd("reco_H");

  h_NrecoH=new TH1D("NrecoH","Number of reconstructed Higgs",10,0.0,10);
  h_recoHPt=new TH1D("recoHPt","P_{T} for reco Higgs",50,0.0,2000.);
  
  //==============for W======================================

  oFile->mkdir("W");
  oFile->cd("W");
  
  h_NrecoW=new TH1D("NrecoW","Number of reconstructed W",10,0.0,10);
  h_recoWPt=new TH1D("recoWPt","P_{T} for reco W",50.0,0.0,2000.);
  h_WdeltaRMin=new TH1D("WDeltaRMin","#delta R between AK8 Jets and W",20,0.,1.0);
  h_WPt=new TH1D("WPt","P_{T} for W",50,100.0,2000.);
  h_WPtdr=new TH1D("WPtdr","P_{T} for W with dr cut of 0.6",50,100.0,2000.);
  h_WPtcut=new TH1D("WPtcut","P_{T} for W with cuts",50,100.0,2000.);

  h_AK8Wtau1=new TH1D("AK8Wtau1","#tau1 for AK8 Jets matched with W",20,0.0,1.0);
  h_AK8Wtau2=new TH1D("AK8Wtau2","#tau2 for AK8 Jets matched with W",20,0.0,1.0);
  h_AK8Wtau3=new TH1D("AK8Wtau3","#tau3 for AK8 Jets matched with W",20,0.0,1.0);
  h_AK8Wtau2bytau1=new TH1D("AK8Wtau2bytau1","#tau_{2}/#tau_{1} for AK8 Jets matched with W",20,0.0,1.0);
  h_AK8Wtau3bytau2=new TH1D("AK8Wtau3bytau2","#tau_{3}/#tau_{2} for AK8 Jets matched with W",20,0.0,1.0);
  h_AK8Wprunedmass=new TH1D("AK8Wprunedmass","Pruned mass for AK8 Jets matched with W",12,0.0,300.);
  h_AK8WSoftDropMass=new TH1D("AK8WSoftDropMass","Soft mass drop for AK8 Jets matched with W",12,0.0,300.0);
  h_AK8WNSubjets=new TH1D("AK8WNSubjets","NSubjets for AK8 Jets matched with W",10,0.0,10);
  h2_AK8Wtau1vs2=new TH2D("AK8Wtau1vs2","#tau_{2} vs #tau_{1} for AK8 Jets matched with W",16,0.0,0.8,16,0.0,0.8);
  h2_AK8Wtau2vs3=new TH2D("AK8Wtau2vs3","#tau_{3} vs #tau_{2} for AK8 Jets matched with W",10,0.0,0.5,10,0.0,0.5);
  //============================= gen leptons=====================================//

  oFile->mkdir("gen_lept");
  oFile->cd("gen_lept");
  h_ZPt=new TH1D("ZPt","P_{T} for Z",50,0.0,1500.);
  h_genEle[0]=new TH1D("ele_l1Pt","P_{T} for leading electron",50,0.0,1000.);
  h_genEle[1]=new TH1D("ele_l2Pt","P_{T} for 2nd leading electron",50,0.0,1000.);
  h_genMu[0]=new TH1D("mu_l1Pt","P_{T} for leading muon",50,0.0,1000.);
  h_genMu[1]=new TH1D("mu_l2Pt","P_{T} for 2nd leading muon",50,0.0,1000.);
  h_genEle[2]=new TH1D("ele_l1Phi","#phi for leading electron",32,3.2,3.2);
  h_genEle[3]=new TH1D("ele_l1Eta","#eta for leading electron",40,-4.,4.);
  h_genEle[4]=new TH1D("ele_l2Phi","#phi for 2nd leading electron",32,3.2,3.2);
  h_genEle[5]=new TH1D("ele_l2Eta","#eta for 2nd leading electron",40,-4.,4.);
  h_genEle[6]=new TH1D("eledPhi","#Delta#phi between two leading electrons",20,0.0,4.);
  h_genEle[7]=new TH1D("eledr","#Delta R between two leading electrons",30,0.0,1.5);
  h_genEle[8]=new TH1D("elemass","Invariant mass for the two leading electrons",20,0.0,200.0);
  h_genEle[9]= new TH1D("ele_st","S_{T} for leading dielectron and AK4 jets",100,0.0,6000.);
  h_genMu[2]=new TH1D("mu_l1Phi","#phi for leading muon",32,3.2,3.2);
  h_genMu[3]=new TH1D("mu_l1Eta","#eta for leading muon",40,-4.,4.);
  h_genMu[4]=new TH1D("mu_l2Phi","#phi for 2nd leading muon",32,3.2,3.2);
  h_genMu[5]=new TH1D("mu_l2Eta","#eta for 2nd leading muon",40,-4.,4.);
  h_genMu[6]=new TH1D("mudr","#Delta R between two leading muons",30,0.0,1.5);
  h_genMu[7]=new TH1D("mudPhi","#Delta#phi between two leading muons",20,0.0,4.);
  h_genMu[8]=new TH1D("mumass","Invariant mass for the two leading muons",20,0.0,200.0);
  h_genMu[9]= new TH1D("mu_st","S_{T} for leading dimuon and AK4 jets",100,0.0,6000.);
  h_genMu_dRvsdPhi= new TH2D("genMu_dRvsdPhi","dR vs d#phi for leading dimuon pair",30,0.0,1.5,20,0.0,4.0);
  h_genEle_dRvsdPhi= new TH2D("genEle_dRvsdPhi","dR vs d#phi for leading dielectron pair",30,0.0,1.5,20,0.0,4.0);

  
  oFile->mkdir("gen_lept_Zcut");
  oFile->cd("gen_lept_Zcut");
  
  h_genEle[10]=new TH1D("ele_l1Pt","P_{T} for leading electron",50,0.0,1000.);
  h_genEle[11]=new TH1D("ele_l2Pt","P_{T} for 2nd leading electron",50,0.0,1000.);
  h_genMu[10]=new TH1D("mu_l1Pt","P_{T} for leading muon",50,0.0,1000.);
  h_genMu[11]=new TH1D("mu_l2Pt","P_{T} for 2nd leading muon",50,0.0,1000.);
  h_genEle[12]=new TH1D("ele_l1Phi","#phi for leading electron",32,3.2,3.2);
  h_genEle[13]=new TH1D("ele_l1Eta","#eta for leading electron",40,-4.,4.);
  h_genEle[14]=new TH1D("ele_l2Phi","#phi for 2nd leading electron",32,3.2,3.2);
  h_genEle[15]=new TH1D("ele_l2Eta","#eta for 2nd leading electron",40,-4.,4.);
  h_genEle[16]=new TH1D("eledPhi","#Delta#phi between two leading electrons",20,0.0,4.);
  h_genEle[17]=new TH1D("eledr","#Delta R between two leading electrons",30,0.0,1.5);
  h_genEle[18]=new TH1D("elemass","Invariant mass for the two leading electrons",20,0.0,200.0);
  h_genEle[19]= new TH1D("ele_st","S_{T} for leading dielectron and AK4 jets",100,0.0,6000.);
  h_genMu[12]=new TH1D("mu_l1Phi","#phi for leading muon",32,3.2,3.2);
  h_genMu[13]=new TH1D("mu_l1Eta","#eta for leading muon",40,-4.,4.);
  h_genMu[14]=new TH1D("mu_l2Phi","#phi for 2nd leading muon",32,3.2,3.2);
  h_genMu[15]=new TH1D("mu_l2Eta","#eta for 2nd leading muon",40,-4.,4.);
  h_genMu[16]=new TH1D("mudr","#Delta R between two leading muons",30,0.0,1.5);
  h_genMu[17]=new TH1D("mudPhi","#Delta#phi between two leading muons",20,0.0,4.);
  h_genMu[18]=new TH1D("mumass","Invariant mass for the two leading muons",20,0.0,200.0);
  h_genMu[19]= new TH1D("mu_st","S_{T} for leading dimuon and AK4 jets",100,0.0,6000.);
  h_genMu_dRvsdPhiZcut= new TH2D("genMu_dRvsdPhi","dR vs d#phi for leading dimuon pair",30,0.0,1.5,20,0.0,4.0);
  h_genEle_dRvsdPhiZcut= new TH2D("genEle_dRvsdPhi","dR vs d#phi for leading dielectron pair",30,0.0,1.5,20,0.0,4.0);
   h_genEle_diPt=new TH1D("dielePt","P_{T} for leading dielectron pair",50,0.0,1500.);
  h_genMu_diPt=new TH1D("dimuPt","P_{T} for leading dimuon pair",50,0.0,1500.);
  //=======================gen lept with dR cut=========================//
  oFile->mkdir("gen_lept_dRcut");
  oFile->cd("gen_lept_dRcut");
  h_genEle[20]=new TH1D("ele_l1Pt","P_{T} for leading electron",50,0.0,1000.);
  h_genEle[21]=new TH1D("ele_l2Pt","P_{T} for 2nd leading electron",50,0.0,1000.);
  
  h_genEle[22]=new TH1D("ele_l1Phi","#phi for leading electron",32,3.2,3.2);
  h_genEle[23]=new TH1D("ele_l1Eta","#eta for leading electron",40,-4.,4.);
  h_genEle[24]=new TH1D("ele_l2Phi","#phi for 2nd leading electron",32,3.2,3.2);
  h_genEle[25]=new TH1D("ele_l2Eta","#eta for 2nd leading electron",40,-4.,4.);
  h_genEle[26]=new TH1D("eledPhi","#Delta#phi between two leading electrons",20,0.0,4.);
  h_genEle[27]=new TH1D("eledr","#Delta R between two leading electrons",30,0.0,1.5);
  h_genEle[28]=new TH1D("elemass","Invariant mass for the two leading electrons",20,0.0,200.0);
  h_genEle[29]= new TH1D("ele_st","S_{T} for leading dielectron and AK4 jets",100,0.0,6000.);

  h_genMu[20]=new TH1D("mu_l1Pt","P_{T} for leading muon",50,0.0,1000.);
  h_genMu[21]=new TH1D("mu_l2Pt","P_{T} for 2nd leading muon",50,0.0,1000.);
  h_genMu[22]=new TH1D("mu_l1Phi","#phi for leading muon",32,-3.2,3.2);
  h_genMu[23]=new TH1D("mu_l1Eta","#eta for leading muon",40,-4.,4.);
  h_genMu[24]=new TH1D("mu_l2Phi","#phi for 2nd leading muon",32,-3.2,3.2);
  h_genMu[25]=new TH1D("mu_l2Eta","#eta for 2nd leading muon",40,-4.,4.);
  h_genMu[26]=new TH1D("mudr","#Delta R between two leading muons",30,0.0,1.5);
  h_genMu[27]=new TH1D("mudPhi","#Delta#phi between two leading muons",20,0.0,4.);
  h_genMu[28]=new TH1D("mumass","Invariant mass for the two leading muons",20,0.0,200.0);
  h_genMu[29]= new TH1D("mu_st","S_{T} for leading dimuon and AK4 jets",100,0.0,6000.);

   //=======================gen lept with dPhi cut=========================//
  oFile->mkdir("gen_lept_dPhicut");
  oFile->cd("gen_lept_dPhicut");
  h_genEle[30]=new TH1D("ele_l1Pt","P_{T} for leading electron",50,0.0,1000.);
  h_genEle[31]=new TH1D("ele_l2Pt","P_{T} for 2nd leading electron",50,0.0,1000.);
  
  h_genEle[32]=new TH1D("ele_l1Phi","#phi for leading electron",32,3.2,3.2);
  h_genEle[33]=new TH1D("ele_l1Eta","#eta for leading electron",40,-4.,4.);
  h_genEle[34]=new TH1D("ele_l2Phi","#phi for 2nd leading electron",32,3.2,3.2);
  h_genEle[35]=new TH1D("ele_l2Eta","#eta for 2nd leading electron",40,-4.,4.);
  h_genEle[36]=new TH1D("eledPhi","#Delta#phi between two leading electrons",20,0.0,4.);
  h_genEle[37]=new TH1D("eledr","#Delta R between two leading electrons",30,0.0,1.5);
  h_genEle[38]=new TH1D("elemass","Invariant mass for the two leading electrons",20,0.0,200.0);
  h_genEle[39]= new TH1D("ele_st","S_{T} for leading dielectron and AK4 jets",100,0.0,6000.);

  h_genMu[30]=new TH1D("mu_l1Pt","P_{T} for leading muon",50,0.0,1000.);
  h_genMu[31]=new TH1D("mu_l2Pt","P_{T} for 2nd leading muon",50,0.0,1000.);
  h_genMu[32]=new TH1D("mu_l1Phi","#phi for leading muon",32,-3.2,3.2);
  h_genMu[33]=new TH1D("mu_l1Eta","#eta for leading muon",40,-4.,4.);
  h_genMu[34]=new TH1D("mu_l2Phi","#phi for 2nd leading muon",32,-3.2,3.2);
  h_genMu[35]=new TH1D("mu_l2Eta","#eta for 2nd leading muon",40,-4.,4.);
  h_genMu[36]=new TH1D("mudr","#Delta R between two leading muons",30,0.0,1.5);
  h_genMu[37]=new TH1D("mudPhi","#Delta#phi between two leading muons",20,0.0,4.);
  h_genMu[38]=new TH1D("mumass","Invariant mass for the two leading muons",20,0.0,200.0);
  h_genMu[39]= new TH1D("mu_st","S_{T} for leading dimuon and AK4 jets",100,0.0,6000.);
  // ====================== for reco leptons=======================================//
  oFile->mkdir("reco_lept");
  oFile->cd("reco_lept");

  h_mu1_iso= new TH1D("mu1_iso","Iso04 for leading muon",40,0.0,0.4);
  h_mu2_iso= new TH1D("mu2_iso","Iso04 for sub leading muon",40,0.0,0.4);
  h_ele1_iso= new TH1D("ele1_iso","Iso03 for leading electron",30,0.0,0.3);
  h_ele2_iso= new TH1D("ele2_iso","Iso03 for sub leading electron",30,0.0,0.3);
  //==================== reco leptons with Z mass cut========================//
  oFile->mkdir("reco_lept_Zcut");
  oFile->cd("reco_lept_Zcut");
   
  h_ele_leadingLepPt_Zcut=new TH1D("ele_l1Pt_Zcut","P_{T} for leading electron",50,0.0,1000.);
  h_ele_leadingLep2Pt_Zcut=new TH1D("ele_l2Pt_Zcut","P_{T} for 2nd leading electron",50,0.0,1000.);
  h_mu_leadingLepPt_Zcut=new TH1D("mu_l1Pt_Zcut","P_{T} for leading muon",50,0.0,1000.);
  h_mu_leadingLep2Pt_Zcut=new TH1D("mu_l2Pt_Zcut","P_{T} for 2nd leading muon",50,0.0,1000.);
  
  h_ele_leadingLepPhi_Zcut=new TH1D("ele_l1Phi_Zcut","#phi for leading electron",32,3.2,3.2);
  h_ele_leadingLepEta_Zcut=new TH1D("ele_l1Eta_Zcut","#eta for leading electron",40,-4.,4.);
  h_ele_leadingLep2Phi_Zcut=new TH1D("ele_l2Phi_Zcut","#phi for 2nd leading electron",32,3.2,3.2);
  h_ele_leadingLep2Eta_Zcut=new TH1D("ele_l2Eta_Zcut","#eta for 2nd leading electron",40,-4.,4.);
  
  h_ele_dPhi_Zcut=new TH1D("eledPhi_Zcut","#Delta#phi between two leading electrons",20,0.0,4.);
  h_ele_dR_Zcut=new TH1D("eledr_Zcut","#Delta R between two leading electrons",30,0.0,1.5);
  h_ele_mass_Zcut=new TH1D("elemass_Zcut","Invariant mass for the two leading electrons",20,0.0,200.0);
  h_ele_st_Zcut= new TH1D("ele_st_Zcut","S_{T} for leading dielectron and AK4 jets",100,0.0,6000.);
  h_ele_dRvsdPhiZcut= new TH2D("ele_dRvsdPhi_Zcut","dR vs d#phi for leading di electron pair",30,0.0,1.5,20,0.0,4.0);
  
  h_mu_leadingLepPhi_Zcut=new TH1D("mu_l1Phi_Zcut","#phi for leading muon",32,-3.2,3.2);
  h_mu_leadingLepEta_Zcut=new TH1D("mu_l1Eta_Zcut","#eta for leading muon",40,-4.,4.);
  h_mu_leadingLep2Phi_Zcut=new TH1D("mu_l2Phi_Zcut","#phi for 2nd leading muon",32,-3.2,3.2);
  h_mu_leadingLep2Eta_Zcut=new TH1D("mu_l2Eta_Zcut","#eta for 2nd leading muon",40,-4.,4.);
  
  h_mu_dR_Zcut=new TH1D("mudr_Zcut","#Delta R between two leading muons",30,0.0,1.5);
  h_mu_dPhi_Zcut=new TH1D("mudPhi_Zcut","#Delta#phi between two leading muons",20,0.0,4.);
  h_mu_mass_Zcut=new TH1D("mumass_Zcut","Invariant mass for the two leading muons",20,0.0,200.0);
  h_mu_st_Zcut= new TH1D("mu_st_Zcut","S_{T} for leading dimuon and AK4 jets",100,0.0,6000.);
  h_mu_dRvsdPhiZcut= new TH2D("mu_dRvsdPhi_Zcut","dR vs d#phi for leading dimuon pair",30,0.0,1.5,20,0.0,4.0);
  h_dielept=new TH1D("dielePt","P_{T} for dielectron",50,0.0,1500.);
  h_dimupt=new TH1D("dimuPt","P_{T} for dimuon",50,0.0,1500.);


   h_Nlep=new TH1D("Nlep","Number of reconstructed Muons+electrons beyond the dilepton pair",10,0.0,10);
  //==================== reco leptons with dPhi mass cut========================//
  oFile->mkdir("reco_lept_dPhicut");
  oFile->cd("reco_lept_dPhicut");
   
  h_ele_leadingLepPt_dPhicut=new TH1D("ele_l1Pt_dPhicut","P_{T} for leading electron",50,0.0,1000.);
  h_ele_leadingLep2Pt_dPhicut=new TH1D("ele_l2Pt_dPhicut","P_{T} for 2nd leading electron",50,0.0,1000.);
  h_mu_leadingLepPt_dPhicut=new TH1D("mu_l1Pt_dPhicut","P_{T} for leading muon",50,0.0,1000.);
  h_mu_leadingLep2Pt_dPhicut=new TH1D("mu_l2Pt_dPhicut","P_{T} for 2nd leading muon",50,0.0,1000.);
  h_ele_leadingLepPhi_dPhicut=new TH1D("ele_l1Phi_dPhicut","#phi for leading electron",32,-3.2,3.2);
  h_ele_leadingLepEta_dPhicut=new TH1D("ele_l1Eta_dPhicut","#eta for leading electron",40,-4.,4.);
  h_ele_leadingLep2Phi_dPhicut=new TH1D("ele_l2Phi_dPhicut","#phi for 2nd leading electron",32,-3.2,3.2);
  h_ele_leadingLep2Eta_dPhicut=new TH1D("ele_l2Eta_dPhicut","#eta for 2nd leading electron",40,-4.,4.);
  h_ele_dPhi_dPhicut=new TH1D("eledPhi_dPhicut","#Delta#phi between two leading electrons",20,0.0,4.);
  h_ele_dR_dPhicut=new TH1D("eledr_dPhicut","#Delta R between two leading electrons",30,0.0,1.5);
  h_ele_mass_dPhicut=new TH1D("elemass_dPhicut","Invariant mass for the two leading electrons",20,0.0,200.0);
  h_ele_st_dPhicut= new TH1D("ele_st_dPhicut","S_{T} for leading dielectron and AK4 jets",100,0.0,6000.);
  h_mu_leadingLepPhi_dPhicut=new TH1D("mu_l1Phi_dPhicut","#phi for leading muon",32,-3.2,3.2);
  h_mu_leadingLepEta_dPhicut=new TH1D("mu_l1Eta_dPhicut","#eta for leading muon",40,-4.,4.);
  h_mu_leadingLep2Phi_dPhicut=new TH1D("mu_l2Phi_dPhicut","#phi for 2nd leading muon",32,-3.2,3.2);
  h_mu_leadingLep2Eta_dPhicut=new TH1D("mu_l2Eta_dPhicut","#eta for 2nd leading muon",40,-4.,4.);
  h_mu_dR_dPhicut=new TH1D("mudr_dPhicut","#Delta R between two leading muons",30,0.0,1.5);
  h_mu_dPhi_dPhicut=new TH1D("mudPhi_dPhicut","#Delta#phi between two leading muons",20,0.0,4.);
  h_mu_mass_dPhicut=new TH1D("mumass_dPhicut","Invariant mass for the two leading muons",20,0.0,200.0);
  h_mu_st_dPhicut= new TH1D("mu_st_dPhicut","S_{T} for leading dimuon and AK4 jets",100,0.0,6000.);


  //==================== reco leptons with dR mass cut========================//
  oFile->mkdir("reco_lept_dRcut");
  oFile->cd("reco_lept_dRcut");
   
  h_ele_leadingLepPt_dRcut=new TH1D("ele_l1Pt_dRcut","P_{T} for leading electron",50,0.0,1000.);
  h_ele_leadingLep2Pt_dRcut=new TH1D("ele_l2Pt_dRcut","P_{T} for 2nd leading electron",50,0.0,1000.);
  h_mu_leadingLepPt_dRcut=new TH1D("mu_l1Pt_dRcut","P_{T} for leading muon",50,0.0,1000.);
  h_mu_leadingLep2Pt_dRcut=new TH1D("mu_l2Pt_dRcut","P_{T} for 2nd leading muon",50,0.0,1000.);
  h_ele_leadingLepPhi_dRcut=new TH1D("ele_l1Phi_dRcut","#phi for leading electron",32,-3.2,3.2);
  h_ele_leadingLepEta_dRcut=new TH1D("ele_l1Eta_dRcut","#eta for leading electron",40,-4.,4.);
  h_ele_leadingLep2Phi_dRcut=new TH1D("ele_l2Phi_dRcut","#phi for 2nd leading electron",32,-3.2,3.2);
  h_ele_leadingLep2Eta_dRcut=new TH1D("ele_l2Eta_dRcut","#eta for 2nd leading electron",40,-4.,4.);
  h_ele_dPhi_dRcut=new TH1D("eledPhi_dRcut","#Delta#phi between two leading electrons",20,0.0,4.);
  h_ele_dR_dRcut=new TH1D("eledr_dRcut","#Delta R between two leading electrons",30,0.0,1.5);
  h_ele_mass_dRcut=new TH1D("elemass_dRcut","Invariant mass for the two leading electrons",20,0.0,200.0);
  h_ele_st_dRcut= new TH1D("ele_st_dRcut","S_{T} for leading dielectron and AK4 jets",100,0.0,6000.);
  h_mu_leadingLepPhi_dRcut=new TH1D("mu_l1Phi_dRcut","#phi for leading muon",32,3.2,3.2);
  h_mu_leadingLepEta_dRcut=new TH1D("mu_l1Eta_dRcut","#eta for leading muon",40,-4.,4.);
  h_mu_leadingLep2Phi_dRcut=new TH1D("mu_l2Phi_dRcut","#phi for 2nd leading muon",32,-3.2,3.2);
  h_mu_leadingLep2Eta_dRcut=new TH1D("mu_l2Eta_dRcut","#eta for 2nd leading muon",40,-4.,4.);
  h_mu_dR_dRcut=new TH1D("mudr_dRcut","#Delta R between two leading muons",30,0.0,1.5);
  h_mu_dPhi_dRcut=new TH1D("mudPhi_dRcut","#Delta#phi between two leading muons",20,0.0,4.);
  h_mu_mass_dRcut=new TH1D("mumass_dRcut","Invariant mass for the two leading muons",20,0.0,200.0);
  h_mu_st_dRcut= new TH1D("mu_st_dRcut","S_{T} for leading dimuon and AK4 jets",100,0.0,6000.);
  // ===============================================================================//
  oFile->mkdir("MET");
  oFile->cd("MET");
  h_met=new TH1D("met","MET",100,0.0,6000.);
  h_metDPhi=new TH1D("met_dPhi","#Delta#phi between Z->l^{+}l^{-} and MET",20,0.0,4.);

  //============================= Data vs MC =========================================//
  oFile->mkdir("DatavsMC");
  oFile->cd("DatavsMC");
  h_Cutflow_DvsMC=new TH1D("Cutflow_DvsMC","Number of Events surviving ",21,0.0,21);
  h_checkTrigNoDZ = new TH1D("checkTrigNoDZ","Trigger check without DZ",3,0.0,3);
  h_checkTrigDZ = new TH1D("checkTrigDZ","Trigger check with DZ",3,0.0,3);

  h_l1Pt = new TH1D("l1Pt","P_{T} for leading lepton",100,0.0,1000.);;
  h_l1Eta = new TH1D("l1Eta","#eta for leading lepton",40,-4.,4.);
  h_l1Phi = new TH1D("l1Phi","#phi for leading lepton",32,-3.2,3.2);
	  

  h_l1Pt_Up = new TH1D("l1Pt_Up","P_{T} for leading lepton",100,0.0,1000.);;
  h_l1Eta_Up = new TH1D("l1Eta_Up","#eta for leading lepton",40,-4.,4.);
  h_l1Phi_Up = new TH1D("l1Phi_Up","#phi for leading lepton",32,-3.2,3.2);

  h_l1Pt_Down = new TH1D("l1Pt_Down","P_{T} for leading lepton",100,0.0,1000.);;
  h_l1Eta_Down = new TH1D("l1Eta_Down","#eta for leading lepton",40,-4.,4.);
  h_l1Phi_Down = new TH1D("l1Phi_Down","#phi for leading lepton",32,-3.2,3.2);

  h_l2Pt = new TH1D("l2Pt","P_{T} for sub leading lepton",100,0.0,1000.);;
  h_l2Eta = new TH1D("l2Eta","#eta for sub leading lepton",40,-4.,4.);
  h_l2Phi = new TH1D("l2Phi","#phi for sub leading lepton",32,-3.2,3.2);

  h_l2Pt_Up = new TH1D("l2Pt_Up","P_{T} for sub leading lepton",100,0.0,1000.);;
  h_l2Eta_Up = new TH1D("l2Eta_Up","#eta for sub leading lepton",40,-4.,4.);
  h_l2Phi_Up = new TH1D("l2Phi_Up","#phi for sub leading lepton",32,-3.2,3.2);

  h_l2Pt_Down = new TH1D("l2Pt_Down","P_{T} for sub leading lepton",100,0.0,1000.);;
  h_l2Eta_Down = new TH1D("l2Eta_Down","#eta for sub leading lepton",40,-4.,4.);
  h_l2Phi_Down = new TH1D("l2Phi_Down","#phi for sub leading lepton",32,-3.2,3.2);

  h_l1l2_dR = new TH1D("l1l2dR","#Delta R between two leading leptons",80,0.0,4.);
  h_l1l2_dPhi = new TH1D("l1l2dPhi","#Delta#phi between two leading leptons",20,0.0,4.);
  h_dilepPt = new TH1D("ZPt","P_{T} for Z",50,0.0,1500.);
  h_dilepEta= new TH1D("ZEta","#eta for Z",40,-4.,4.);
  h_dilepPhi = new TH1D("ZPhi","#phi for Z",32,-3.2,3.2);
  h_dilepMass=new TH1D("dilepMAss","Invariant mass for the two leading leptons",100,0.0,200.0);
  h_dilepMass_NoCut=new TH1D("dilepMAss_NoCut","Invariant mass for the two leading leptons",75.,50.0,200.0);
  h_ST=new TH1D("STAK4","S_{T} for AK4 Jets and leading dilepton pair(ee/ #mu#mu)",120,0.0,6000.);
  h_npv=new TH1D("NPV","Number of PV",50,0.0,50);
  h_NbjetL=new TH1D("NrecobL","Number of reconstructed b quarks (CSV>0.5426)",10,0.0,10);
  h_NbjetM=new TH1D("NrecobM","Number of reconstructed b quarks (CSV>0.8484)",10,0.0,10);
  h_NAK4=new TH1D("AK4Njets","Number of AK4 Jets",21,0.0,21);
  h_NAK8=new TH1D("AK8Njets","Number of AK8 Jets",21,0.0,21);
  
  h_l1l2_dR_Up = new TH1D("l1l2dR_Up","#Delta R between two leading leptons",80,0.0,4.);
  h_l1l2_dPhi_Up = new TH1D("l1l2dPhi_Up","#Delta#phi between two leading leptons",20,0.0,4.);
  h_dilepPt_Up = new TH1D("ZPt_Up","P_{T} for Z",50,0.0,1500.);
  h_dilepEta_Up= new TH1D("ZEta_Up","#eta for Z",40,-4.,4.);
  h_dilepPhi_Up = new TH1D("ZPhi_Up","#phi for Z",32,-3.2,3.2);
  h_dilepMass_Up =new TH1D("dilepMAss_Up","Invariant mass for the two leading leptons",100,0.0,200.0);
  h_dilepMass_NoCut_Up =new TH1D("dilepMAss_NoCut_Up","Invariant mass for the two leading leptons",75.,50.0,200.0);
  h_ST_Up =new TH1D("STAK4_Up","S_{T} for AK4 Jets and leading dilepton pair(ee/ #mu#mu)",120,0.0,6000.);
  h_npv_Up =new TH1D("NPV_Up","Number of PV",50,0.0,50);
  h_NbjetL_Up =new TH1D("NrecobL_Up","Number of reconstructed b quarks (CSV>0.5426)",10,0.0,10);
  h_NbjetM_Up =new TH1D("NrecobM_Up","Number of reconstructed b quarks (CSV>0.8484)",10,0.0,10);
  h_NAK4_Up =new TH1D("AK4Njets_Up","Number of AK4 Jets",21,0.0,21);
  h_NAK8_Up =new TH1D("AK8Njets_Up","Number of AK8 Jets",21,0.0,21);

  h_l1l2_dR_Down = new TH1D("l1l2dR_Down","#Delta R between two leading leptons",80,0.0,4.);
  h_l1l2_dPhi_Down = new TH1D("l1l2dPhi_Down","#Delta#phi between two leading leptons",20,0.0,4.);
  h_dilepPt_Down = new TH1D("ZPt_Down","P_{T} for Z",50,0.0,1500.);
  h_dilepEta_Down= new TH1D("ZEta_Down","#eta for Z",40,-4.,4.);
  h_dilepPhi_Down = new TH1D("ZPhi_Down","#phi for Z",32,-3.2,3.2);
  h_dilepMass_Down =new TH1D("dilepMAss_Down","Invariant mass for the two leading leptons",100,0.0,200.0);
  h_dilepMass_NoCut_Down =new TH1D("dilepMAss_NoCut_Down","Invariant mass for the two leading leptons",75.,50.0,200.0);
  h_ST_Down =new TH1D("STAK4_Down","S_{T} for AK4 Jets and leading dilepton pair(ee/ #mu#mu)",120,0.0,6000.);
  h_npv_Down =new TH1D("NPV_Down","Number of PV",50,0.0,50);
  h_NbjetL_Down =new TH1D("NrecobL_Down","Number of reconstructed b quarks (CSV>0.5426)",10,0.0,10);
  h_NbjetM_Down =new TH1D("NrecobM_Down","Number of reconstructed b quarks (CSV>0.8484)",10,0.0,10);
  h_NAK4_Down =new TH1D("AK4Njets_Down","Number of AK4 Jets",21,0.0,21);
  h_NAK8_Down =new TH1D("AK8Njets_Down","Number of AK8 Jets",21,0.0,21);

  h_AK8Pt=new TH1D("JetsAK8Pt","P_{T} for AK8 Jets",50,0.0,2000.);
  h_AK8Phi=new TH1D("JetsAK8Phi","#phi for AK8 Jets",32,-3.2,3.2);
  h_AK8Eta=new TH1D("JetsAK8Eta","#eta for AK8 Jets",40,-4.,4.);
  h_AK8HT=new TH1D("JetsAK8HT","H_{T} for AK8 Jets",100,0.0,6000.);
  h_LeadAK8Pt=new TH1D("LeadingAK8Pt","P_{T} for Leading AK8 Jet",50,0.0,2000.);
  h_LeadAK8Phi=new TH1D("LeadingAK8Phi","#phi for AK8 Jets",32,-3.2,3.2);
  h_LeadAK8Eta=new TH1D("LeadingAK8Eta","#eta for AK8 Jets",40,-4.,4.);

  h_AK8Pt_Up = new TH1D("JetsAK8Pt_Up","P_{T} for AK8 Jets",50,0.0,2000.);
  h_AK8Phi_Up = new TH1D("JetsAK8Phi_Up","#phi for AK8 Jets",32,-3.2,3.2);
  h_AK8Eta_Up = new TH1D("JetsAK8Eta_Up","#eta for AK8 Jets",40,-4.,4.);
  h_AK8HT_Up = new TH1D("JetsAK8HT_Up","H_{T} for AK8 Jets",100,0.0,6000.);
  h_LeadAK8Pt_Up = new TH1D("LeadingAK8Pt_Up","P_{T} for Leading AK8 Jet",50,0.0,2000.);
  h_LeadAK8Phi_Up =new TH1D("LeadingAK8Phi_Up","#phi for AK8 Jets",32,-3.2,3.2);
  h_LeadAK8Eta_Up =new TH1D("LeadingAK8Eta_Up","#eta for AK8 Jets",40,-4.,4.);
  
  h_AK8Pt_Down = new TH1D("JetsAK8Pt_Down","P_{T} for AK8 Jets",50,0.0,2000.);
  h_AK8Phi_Down = new TH1D("JetsAK8Phi_Down","#phi for AK8 Jets",32,-3.2,3.2);
  h_AK8Eta_Down = new TH1D("JetsAK8Eta_Down","#eta for AK8 Jets",40,-4.,4.);
  h_AK8HT_Down = new TH1D("JetsAK8HT_Down","H_{T} for AK8 Jets",100,0.0,6000.);
  h_LeadAK8Pt_Down = new TH1D("LeadingAK8Pt_Down","P_{T} for Leading AK8 Jet",50,0.0,2000.);
  h_LeadAK8Phi_Down =new TH1D("LeadingAK8Phi_Down","#phi for AK8 Jets",32,-3.2,3.2);
  h_LeadAK8Eta_Down =new TH1D("LeadingAK8Eta_Down","#eta for AK8 Jets",40,-4.,4.);

  h_AK4Pt=new TH1D("JetsAK4Pt","P_{T} for AK4 Jets",50,0.0,2000.);
  h_AK4Phi=new TH1D("JetsAK4Phi","#phi for AK4 Jets",32,-3.2,3.2);
  h_AK4Eta=new TH1D("JetsAK4Eta","#eta for AK4 Jets",40,-4.,4.);
  h_AK4HT=new TH1D("JetsAK4HT","H_{T} for AK4 Jets",100,0.0,6000.);
  h_LeadAK4Pt=new TH1D("LeadingAK4Pt","P_{T} for Leading AK4 Jet",50,0.0,2000.);
  h_LeadAK4Phi=new TH1D("LeadingAK4Phi","#phi for AK4 Jets",32,-3.2,3.2);
  h_LeadAK4Eta=new TH1D("LeadingAK4Eta","#eta for AK4 Jets",40,-4.,4.);

  h_AK4Pt_Up=new TH1D("JetsAK4Pt_Up","P_{T} for AK4 Jets",50,0.0,2000.);
  h_AK4Phi_Up=new TH1D("JetsAK4Phi_Up","#phi for AK4 Jets",32,-3.2,3.2);
  h_AK4Eta_Up=new TH1D("JetsAK4Eta_Up","#eta for AK4 Jets",40,-4.,4.);
  h_AK4HT_Up=new TH1D("JetsAK4HT_Up","H_{T} for AK4 Jets",100,0.0,6000.);
  h_LeadAK4Pt_Up=new TH1D("LeadingAK4Pt_Up","P_{T} for Leading AK4 Jet",50,0.0,2000.);
  h_LeadAK4Phi_Up=new TH1D("LeadingAK4Phi_Up","#phi for AK4 Jets",32,-3.2,3.2);
  h_LeadAK4Eta_Up=new TH1D("LeadingAK4Eta_Up","#eta for AK4 Jets",40,-4.,4.);

  h_AK4Pt_Down=new TH1D("JetsAK4Pt_Down","P_{T} for AK4 Jets",50,0.0,2000.);
  h_AK4Phi_Down=new TH1D("JetsAK4Phi_Down","#phi for AK4 Jets",32,-3.2,3.2);
  h_AK4Eta_Down=new TH1D("JetsAK4Eta_Down","#eta for AK4 Jets",40,-4.,4.);
  h_AK4HT_Down=new TH1D("JetsAK4HT_Down","H_{T} for AK4 Jets",100,0.0,6000.);
  h_LeadAK4Pt_Down=new TH1D("LeadingAK4Pt_Down","P_{T} for Leading AK4 Jet",50,0.0,2000.);
  h_LeadAK4Phi_Down=new TH1D("LeadingAK4Phi_Down","#phi for AK4 Jets",32,-3.2,3.2);
  h_LeadAK4Eta_Down=new TH1D("LeadingAK4Eta_Down","#eta for AK4 Jets",40,-4.,4.);


  h_Mu_Chi2=new TH1D("Chi2","Chi2",60,0.,240.);
  h_Mu_Chi2_Up=new TH1D("Chi2_Up","Chi2",60,0.,240.);
  h_Mu_Chi2_Down=new TH1D("Chi2_Down","Chi2",60,0.,240.);
  h_Mu_Chi2_1t=new TH1D("Chi2_1t","Chi2",60,0.,240.);
  h_Mu_Chi2_1t_Up=new TH1D("Chi2_1t_Up","Chi2",60,0.,240.);
  h_Mu_Chi2_1t_Down=new TH1D("Chi2_1t_Down","Chi2",60,0.,240.);
  h_Mu_Chi2_0t1V=new TH1D("Chi2_0t1V","Chi2",60,0.,240.);
  h_Mu_Chi2_0t1V_Up=new TH1D("Chi2_0t1V_Up","Chi2",60,0.,240.);
  h_Mu_Chi2_0t1V_Down=new TH1D("Chi2_0t1V_Down","Chi2",60,0.,240.);
  h_Mu_Chi2_0t2V=new TH1D("Chi2_0t2V","Chi2 ",60,0.,240.);
  h_Mu_Chi2_0t2V_Up=new TH1D("Chi2_0t2V_Up","Chi2 ",60,0.,240.);
  h_Mu_Chi2_0t2V_Down=new TH1D("Chi2_0t2V_Down","Chi2 ",60,0.,240.);
  h_Mu_Chi2_1t1V=new TH1D("Chi2_1t1V","Chi2",60,0.,240.);
  h_Mu_Chi2_1t1V_Up=new TH1D("Chi2_1t1V_Up","Chi2",60,0.,240.);
  h_Mu_Chi2_1t1V_Down=new TH1D("Chi2_1t1V_Down","Chi2",60,0.,240.);
  h_Mu_Chi2_1t0V=new TH1D("Chi2_1t0V","Chi2",60,0.,240.);
  h_Mu_Chi2_1t0V_Up=new TH1D("Chi2_1t0V_Up","Chi2",60,0.,240.);
  h_Mu_Chi2_1t0V_Down=new TH1D("Chi2_1t0V_Down","Chi2",60,0.,240.);
  
  h_Mu_TMass_chi2_lep=new TH1D("TMass_chi2_lep","reconstructed TPrime Mass (leptonic)",30,0.0,3000.);
  h_Mu_TMass_chi2_lep_Up=new TH1D("TMass_chi2_lep_Up","reconstructed TPrime Mass (leptonic)",30,0.0,3000.);
  h_Mu_TMass_chi2_lep_Down=new TH1D("TMass_chi2_lep_Down","reconstructed TPrime Mass (leptonic)",30,0.0,3000.);
  h_Mu_TMass_chi2_had=new TH1D("TMass_chi2_had","reconstructed TPrime Mass (hadronic)",30,0.0,3000.);
  h_Mu_TMass_chi2_had_Up=new TH1D("TMass_chi2_had_Up","reconstructed TPrime Mass (hadronic)",30,0.0,3000.);
  h_Mu_TMass_chi2_had_Down=new TH1D("TMass_chi2_had_Down","reconstructed TPrime Mass (hadronic)",30,0.0,3000.);
  h_Mu_TMass_chi2=new TH1D("TMass_chi2","reconstructed TPrime Mass",30,0.0,3000.);
  h_Mu_TMass_chi2_Up=new TH1D("TMass_chi2_Up","reconstructed TPrime Mass",30,0.0,3000.);
  h_Mu_TMass_chi2_Down=new TH1D("TMass_chi2_Down","reconstructed TPrime Mass",30,0.0,3000.);
  h_Mu_Chi2_cut=new TH1D("Chi2_c","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Mu_Chi2_cut_Up=new TH1D("Chi2_c_Up","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Mu_Chi2_cut_Down=new TH1D("Chi2_c_Down","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Mu_Chi2_1t_cut=new TH1D("Chi2_1t_c","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Mu_Chi2_1t_cut_Up=new TH1D("Chi2_1t_c_Up","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Mu_Chi2_1t_cut_Down=new TH1D("Chi2_1t_c_Down","Chi2 ( #chi^{2}<20.)",20,0.,40.);

  h_Mu_Chi2_0t1V_cut=new TH1D("Chi2_0t1V_c","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Mu_Chi2_0t1V_cut_Up=new TH1D("Chi2_0t1V_c_Up","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Mu_Chi2_0t1V_cut_Down=new TH1D("Chi2_0t1V_c_Down","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Mu_Chi2_0t2V_cut=new TH1D("Chi2_0t2V_c","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Mu_Chi2_0t2V_cut_Up=new TH1D("Chi2_0t2V_c_Up","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Mu_Chi2_0t2V_cut_Down=new TH1D("Chi2_0t2V_c_Down","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Mu_Chi2_1t1V_cut=new TH1D("Chi2_1t1V_c","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Mu_Chi2_1t1V_cut_Up=new TH1D("Chi2_1t1V_c_Up","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Mu_Chi2_1t1V_cut_Down=new TH1D("Chi2_1t1V_c_Down","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Mu_Chi2_1t0V_cut=new TH1D("Chi2_1t0V_c","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Mu_Chi2_1t0V_cut_Up=new TH1D("Chi2_1t0V_c_Up","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Mu_Chi2_1t0V_cut_Down=new TH1D("Chi2_1t0V_c_Down","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Mu_TMass_chi2_lep_cut=new TH1D("TMass_chi2_lep_c","reconstructed TPrime Mass (leptonic) ( #chi^{2}<20.)",30,0.0,3000.);
  h_Mu_TMass_chi2_lep_cut_Up=new TH1D("TMass_chi2_lep_c_Up","reconstructed TPrime Mass (leptonic) ( #chi^{2}<20.)",30,0.0,3000.);
  h_Mu_TMass_chi2_lep_cut_Down=new TH1D("TMass_chi2_lep_c_Down","reconstructed TPrime Mass (leptonic) ( #chi^{2}<20.)",30,0.0,3000.);
  h_Mu_TMass_chi2_had_cut=new TH1D("TMass_chi2_had_c","reconstructed TPrime Mass (hadronic) ( #chi^{2}<20.)",30,0.0,3000.);
  h_Mu_TMass_chi2_had_cut_Up=new TH1D("TMass_chi2_had_c_Up","reconstructed TPrime Mass (hadronic) ( #chi^{2}<20.)",30,0.0,3000.);
  h_Mu_TMass_chi2_had_cut_Down=new TH1D("TMass_chi2_had_c_Down","reconstructed TPrime Mass (hadronic) ( #chi^{2}<20.)",30,0.0,3000.);
  h_Mu_TMass_chi2_cut=new TH1D("TMass_chi2_c","reconstructed TPrime Mass ( #chi^{2}<20.)",30,0.0,3000.);
  h_Mu_TMass_chi2_cut_Up=new TH1D("TMass_chi2_c_Up","reconstructed TPrime Mass ( #chi^{2}<20.)",30,0.0,3000.);
  h_Mu_TMass_chi2_cut_Down=new TH1D("TMass_chi2_c_Down","reconstructed TPrime Mass ( #chi^{2}<20.)",30,0.0,3000.);
  
  //x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x=x//


  h_Ele_Cutflow_DvsMC=new TH1D("Ele_Cutflow_DvsMC","Number of Events surviving ",21,0.0,21);
  h_Ele_checkTrigNoDZ = new TH1D("Ele_checkTrigNoDZ","Trigger check without DZ",3,0.0,3);
  h_Ele_checkTrigDZ = new TH1D("Ele_checkTrigDZ","Trigger check with DZ",3,0.0,3);
  h_Ele_l1Pt = new TH1D("Ele_l1Pt","P_{T} for leading lepton",100,0.0,1000.);;
  h_Ele_l1Eta = new TH1D("Ele_l1Eta","#eta for leading lepton",40,-4.,4.);
  h_Ele_l1Phi = new TH1D("Ele_l1Phi","#phi for leading lepton",32,-3.2,3.2);
	  
  h_Ele_l2Pt = new TH1D("Ele_l2Pt","P_{T} for sub leading lepton",100,0.0,1000.);;
  h_Ele_l2Eta = new TH1D("Ele_l2Eta","#eta for sub leading lepton",40,-4.,4.);
  h_Ele_l2Phi = new TH1D("Ele_l2Phi","#phi for sub leading lepton",32,-3.2,3.2);

  h_Ele_l1l2_dR = new TH1D("Ele_l1l2dR","#Delta R between two leading leptons",80,0.0,4.);
  h_Ele_l1l2_dPhi = new TH1D("Ele_l1l2dPhi","#Delta#phi between two leading leptons",20,0.0,4.);
  h_Ele_dilepPt = new TH1D("Ele_ZPt","P_{T} for Z",50,0.0,1500.);
  h_Ele_dilepEta= new TH1D("Ele_ZEta","#eta for Z",40,-4.,4.);
  h_Ele_dilepPhi = new TH1D("Ele_ZPhi","#phi for Z",32,-3.2,3.2);
  h_Ele_dilepMass=new TH1D("Ele_dilepMAss","Invariant mass for the two leading leptons",100,0.0,200.0);
  h_Ele_dilepMass_NoCut=new TH1D("Ele_dilepMAss_Nocut","Invariant mass for the two leading leptons",75,50.0,200.0);
  h_Ele_ST=new TH1D("Ele_STAK4","S_{T} for AK4 Jets and leading dilepton pair(ee/ #mu#mu)",120,0.0,6000.);
  h_Ele_npv=new TH1D("Ele_NPV","Number of PV",50,0.0,50);
  h_Ele_NbjetL=new TH1D("Ele_NrecobL","Number of reconstructed b quarks (CSV>0.5426)",10,0.0,10);
  h_Ele_NbjetM=new TH1D("Ele_NrecobM","Number of reconstructed b quarks (CSV>0.8484)",10,0.0,10);
  h_Ele_NAK4=new TH1D("Ele_AK4Njets","Number of AK4 Jets",21,0.0,21);
  h_Ele_NAK8=new TH1D("Ele_AK8Njets","Number of AK8 Jets",21,0.0,21);
  
  h_Ele_AK8Pt=new TH1D("Ele_JetsAK8Pt","P_{T} for AK8 Jets",50,0.0,2000.);
  h_Ele_AK8Phi=new TH1D("Ele_JetsAK8Phi","#phi for AK8 Jets",32,-3.2,3.2);
  h_Ele_AK8Eta=new TH1D("Ele_JetsAK8Eta","#eta for AK8 Jets",40,-4.,4.);
  h_Ele_AK8HT=new TH1D("Ele_JetsAK8HT","H_{T} for AK8 Jets",100,0.0,6000.);
  h_Ele_LeadAK8Pt=new TH1D("Ele_LeadingAK8Pt","P_{T} for Leading AK8 Jet",50,0.0,2000.);
  h_Ele_LeadAK8Phi=new TH1D("Ele_LeadingAK8Phi","#phi for AK8 Jets",32,-3.2,3.2);
  h_Ele_LeadAK8Eta=new TH1D("Ele_LeadingAK8Eta","#eta for AK8 Jets",40,-4.,4.);
  
  h_Ele_AK4Pt=new TH1D("Ele_JetsAK4Pt","P_{T} for AK4 Jets",50,0.0,2000.);
  h_Ele_AK4Phi=new TH1D("Ele_JetsAK4Phi","#phi for AK4 Jets",32,-3.2,3.2);
  h_Ele_AK4Eta=new TH1D("Ele_JetsAK4Eta","#eta for AK4 Jets",40,-4.,4.);
  h_Ele_AK4HT=new TH1D("Ele_JetsAK4HT","H_{T} for AK4 Jets",100,0.0,6000.);
  h_Ele_LeadAK4Pt=new TH1D("Ele_LeadingAK4Pt","P_{T} for Leading AK4 Jet",50,0.0,2000.);
  h_Ele_LeadAK4Phi=new TH1D("Ele_LeadingAK4Phi","#phi for AK4 Jets",32,-3.2,3.2);
  h_Ele_LeadAK4Eta=new TH1D("Ele_LeadingAK4Eta","#eta for AK4 Jets",40,-4.,4.);


  h_Ele_l1Pt_Up = new TH1D("Ele_l1Pt_Up","P_{T} for leading lepton",100,0.0,1000.);;
  h_Ele_l1Eta_Up = new TH1D("Ele_l1Eta_Up","#eta for leading lepton",40,-4.,4.);
  h_Ele_l1Phi_Up = new TH1D("Ele_l1Phi_Up","#phi for leading lepton",32,-3.2,3.2);
	  
  h_Ele_l2Pt_Up = new TH1D("Ele_l2Pt_Up","P_{T} for sub leading lepton",100,0.0,1000.);;
  h_Ele_l2Eta_Up = new TH1D("Ele_l2Eta_Up","#eta for sub leading lepton",40,-4.,4.);
  h_Ele_l2Phi_Up = new TH1D("Ele_l2Phi_Up","#phi for sub leading lepton",32,-3.2,3.2);

  h_Ele_l1l2_dR_Up = new TH1D("Ele_l1l2dR_Up","#Delta R between two leading leptons",80,0.0,4.);
  h_Ele_l1l2_dPhi_Up = new TH1D("Ele_l1l2dPhi_Up","#Delta#phi between two leading leptons",20,0.0,4.);
  h_Ele_dilepPt_Up = new TH1D("Ele_ZPt_Up","P_{T} for Z",50,0.0,1500.);
  h_Ele_dilepEta_Up = new TH1D("Ele_ZEta_Up","#eta for Z",40,-4.,4.);
  h_Ele_dilepPhi_Up = new TH1D("Ele_ZPhi_Up","#phi for Z",32,-3.2,3.2);
  h_Ele_dilepMass_Up =new TH1D("Ele_dilepMAss_Up","Invariant mass for the two leading leptons",100,0.0,200.0);
  h_Ele_dilepMass_NoCut_Up =new TH1D("Ele_dilepMAss_Nocut_Up","Invariant mass for the two leading leptons",75,50.0,200.0);
  h_Ele_ST_Up =new TH1D("Ele_STAK4_Up","S_{T} for AK4 Jets and leading dilepton pair(ee/ #mu#mu)",120,0.0,6000.);
  h_Ele_npv_Up =new TH1D("Ele_NPV_Up","Number of PV",50,0.0,50);
  h_Ele_NbjetL_Up =new TH1D("Ele_NrecobL_Up","Number of reconstructed b quarks (CSV>0.5426)",10,0.0,10);
  h_Ele_NbjetM_Up =new TH1D("Ele_NrecobM_Up","Number of reconstructed b quarks (CSV>0.8484)",10,0.0,10);
  h_Ele_NAK4_Up =new TH1D("Ele_AK4Njets_Up","Number of AK4 Jets",21,0.0,21);
  h_Ele_NAK8_Up =new TH1D("Ele_AK8Njets_Up","Number of AK8 Jets",21,0.0,21);
  
  h_Ele_AK8Pt_Up =new TH1D("Ele_JetsAK8Pt_Up","P_{T} for AK8 Jets",50,0.0,2000.);
  h_Ele_AK8Phi_Up=new TH1D("Ele_JetsAK8Phi_Up","#phi for AK8 Jets",32,-3.2,3.2);
  h_Ele_AK8Eta_Up =new TH1D("Ele_JetsAK8Eta_Up","#eta for AK8 Jets",40,-4.,4.);
  h_Ele_AK8HT_Up =new TH1D("Ele_JetsAK8HT_Up","H_{T} for AK8 Jets",100,0.0,6000.);
  h_Ele_LeadAK8Pt_Up =new TH1D("Ele_LeadingAK8Pt_Up","P_{T} for Leading AK8 Jet",50,0.0,2000.);
  h_Ele_LeadAK8Phi_Up =new TH1D("Ele_LeadingAK8Phi_Up","#phi for AK8 Jets",32,-3.2,3.2);
  h_Ele_LeadAK8Eta_Up =new TH1D("Ele_LeadingAK8Eta_Up","#eta for AK8 Jets",40,-4.,4.);
  
  h_Ele_AK4Pt_Up =new TH1D("Ele_JetsAK4Pt_Up","P_{T} for AK4 Jets",50,0.0,2000.);
  h_Ele_AK4Phi_Up=new TH1D("Ele_JetsAK4Phi_Up","#phi for AK4 Jets",32,-3.2,3.2);
  h_Ele_AK4Eta_Up=new TH1D("Ele_JetsAK4Eta_Up","#eta for AK4 Jets",40,-4.,4.);
  h_Ele_AK4HT_Up=new TH1D("Ele_JetsAK4HT_Up","H_{T} for AK4 Jets",100,0.0,6000.);
  h_Ele_LeadAK4Pt_Up=new TH1D("Ele_LeadingAK4Pt_Up","P_{T} for Leading AK4 Jet",50,0.0,2000.);
  h_Ele_LeadAK4Phi_Up=new TH1D("Ele_LeadingAK4Phi_Up","#phi for AK4 Jets",32,-3.2,3.2);
  h_Ele_LeadAK4Eta_Up=new TH1D("Ele_LeadingAK4Eta_Up","#eta for AK4 Jets",40,-4.,4.);




  
  h_Ele_l1Pt_Down = new TH1D("Ele_l1Pt_Down","P_{T} for leading lepton",100,0.0,1000.);;
  h_Ele_l1Eta_Down = new TH1D("Ele_l1Eta_Down","#eta for leading lepton",40,-4.,4.);
  h_Ele_l1Phi_Down = new TH1D("Ele_l1Phi_Down","#phi for leading lepton",32,-3.2,3.2);
	  
  h_Ele_l2Pt_Down = new TH1D("Ele_l2Pt_Down","P_{T} for sub leading lepton",100,0.0,1000.);;
  h_Ele_l2Eta_Down = new TH1D("Ele_l2Eta_Down","#eta for sub leading lepton",40,-4.,4.);
  h_Ele_l2Phi_Down = new TH1D("Ele_l2Phi_Down","#phi for sub leading lepton",32,-3.2,3.2);

  h_Ele_l1l2_dR_Down = new TH1D("Ele_l1l2dR_Down","#Delta R between two leading leptons",80,0.0,4.);
  h_Ele_l1l2_dPhi_Down = new TH1D("Ele_l1l2dPhi_Down","#Delta#phi between two leading leptons",20,0.0,4.);
  h_Ele_dilepPt_Down = new TH1D("Ele_ZPt_Down","P_{T} for Z",50,0.0,1500.);
  h_Ele_dilepEta_Down= new TH1D("Ele_ZEta_Down","#eta for Z",40,-4.,4.);
  h_Ele_dilepPhi_Down = new TH1D("Ele_ZPhi_Down","#phi for Z",32,-3.2,3.2);
  h_Ele_dilepMass_Down =new TH1D("Ele_dilepMAss_Down","Invariant mass for the two leading leptons",100,0.0,200.0);
  h_Ele_dilepMass_NoCut_Down =new TH1D("Ele_dilepMAss_Nocut_Down","Invariant mass for the two leading leptons",75,50.0,200.0);
  h_Ele_ST_Down =new TH1D("Ele_STAK4_Down","S_{T} for AK4 Jets and leading dilepton pair(ee/ #mu#mu)",120,0.0,6000.);
  h_Ele_npv_Down =new TH1D("Ele_NPV_Down","Number of PV",50,0.0,50);
  h_Ele_NbjetL_Down =new TH1D("Ele_NrecobL_Down","Number of reconstructed b quarks (CSV>0.5426)",10,0.0,10);
  h_Ele_NbjetM_Down =new TH1D("Ele_NrecobM_Down","Number of reconstructed b quarks (CSV>0.8484)",10,0.0,10);
  h_Ele_NAK4_Down =new TH1D("Ele_AK4Njets_Down","Number of AK4 Jets",21,0.0,21);
  h_Ele_NAK8_Down =new TH1D("Ele_AK8Njets_Down","Number of AK8 Jets",21,0.0,21);
  
  h_Ele_AK8Pt_Down =new TH1D("Ele_JetsAK8Pt_Down","P_{T} for AK8 Jets",50,0.0,2000.);
  h_Ele_AK8Phi_Down =new TH1D("Ele_JetsAK8Phi_Down","#phi for AK8 Jets",32,-3.2,3.2);
  h_Ele_AK8Eta_Down =new TH1D("Ele_JetsAK8Eta_Down","#eta for AK8 Jets",40,-4.,4.);
  h_Ele_AK8HT_Down =new TH1D("Ele_JetsAK8HT_Down","H_{T} for AK8 Jets",100,0.0,6000.);
  h_Ele_LeadAK8Pt_Down =new TH1D("Ele_LeadingAK8Pt_Down","P_{T} for Leading AK8 Jet",50,0.0,2000.);
  h_Ele_LeadAK8Phi_Down =new TH1D("Ele_LeadingAK8Phi_Down","#phi for AK8 Jets",32,-3.2,3.2);
  h_Ele_LeadAK8Eta_Down =new TH1D("Ele_LeadingAK8Eta_Down","#eta for AK8 Jets",40,-4.,4.);
  
  h_Ele_AK4Pt_Down =new TH1D("Ele_JetsAK4Pt_Down","P_{T} for AK4 Jets",50,0.0,2000.);
  h_Ele_AK4Phi_Down =new TH1D("Ele_JetsAK4Phi_Down","#phi for AK4 Jets",32,-3.2,3.2);
  h_Ele_AK4Eta_Down =new TH1D("Ele_JetsAK4Eta_Down","#eta for AK4 Jets",40,-4.,4.);
  h_Ele_AK4HT_Down =new TH1D("Ele_JetsAK4HT_Down","H_{T} for AK4 Jets",100,0.0,6000.);
  h_Ele_LeadAK4Pt_Down =new TH1D("Ele_LeadingAK4Pt_Down","P_{T} for Leading AK4 Jet",50,0.0,2000.);
  h_Ele_LeadAK4Phi_Down =new TH1D("Ele_LeadingAK4Phi_Down","#phi for AK4 Jets",32,-3.2,3.2);
  h_Ele_LeadAK4Eta_Down =new TH1D("Ele_LeadingAK4Eta_Down","#eta for AK4 Jets",40,-4.,4.);

  
  h_Ele_Chi2=new TH1D("Ele_Chi2","Chi2",60,0.,240.);
  h_Ele_Chi2_Up=new TH1D("Ele_Chi2_Up","Chi2",60,0.,240.);
  h_Ele_Chi2_Down=new TH1D("Ele_Chi2_Down","Chi2",60,0.,240.);
  h_Ele_Chi2_1t=new TH1D("Ele_Chi2_1t","Chi2",60,0.,240.);
  h_Ele_Chi2_1t_Up=new TH1D("Ele_Chi2_1t_Up","Chi2",60,0.,240.);
  h_Ele_Chi2_1t_Down=new TH1D("Ele_Chi2_1t_Down","Chi2",60,0.,240.);
  h_Ele_Chi2_0t1V=new TH1D("Ele_Chi2_0t1V","Chi2",60,0.,240.);
  h_Ele_Chi2_0t1V_Up=new TH1D("Ele_Chi2_0t1V_Up","Chi2",60,0.,240.);
  h_Ele_Chi2_0t1V_Down=new TH1D("Ele_Chi2_0t1V_Down","Chi2",60,0.,240.);
  h_Ele_Chi2_0t2V=new TH1D("Ele_Chi2_0t2V","Chi2 ",60,0.,240.);
  h_Ele_Chi2_0t2V_Up=new TH1D("Ele_Chi2_0t2V_Up","Chi2 ",60,0.,240.);
  h_Ele_Chi2_0t2V_Down=new TH1D("Ele_Chi2_0t2V_Down","Chi2 ",60,0.,240.);
  h_Ele_Chi2_1t1V=new TH1D("Ele_Chi2_1t1V","Chi2",60,0.,240.);
  h_Ele_Chi2_1t1V_Up=new TH1D("Ele_Chi2_1t1V_Up","Chi2",60,0.,240.);
  h_Ele_Chi2_1t1V_Down=new TH1D("Ele_Chi2_1t1V_Down","Chi2",60,0.,240.);
  h_Ele_Chi2_1t0V=new TH1D("Ele_Chi2_1t0V","Chi2",60,0.,240.);
  h_Ele_Chi2_1t0V_Up=new TH1D("Ele_Chi2_1t0V_Up","Chi2",60,0.,240.);
  h_Ele_Chi2_1t0V_Down=new TH1D("Ele_Chi2_1t0V_Down","Chi2",60,0.,240.);
  
  h_Ele_TMass_chi2_lep=new TH1D("Ele_TMass_chi2_lep","reconstructed TPrime Mass (leptonic)",30,0.0,3000.);
  h_Ele_TMass_chi2_lep_Up=new TH1D("Ele_TMass_chi2_lep_Up","reconstructed TPrime Mass (leptonic)",30,0.0,3000.);
  h_Ele_TMass_chi2_lep_Down=new TH1D("Ele_TMass_chi2_lep_Down","reconstructed TPrime Mass (leptonic)",30,0.0,3000.);
  h_Ele_TMass_chi2_had=new TH1D("Ele_TMass_chi2_had","reconstructed TPrime Mass (hadronic)",30,0.0,3000.);
  h_Ele_TMass_chi2_had_Up=new TH1D("Ele_TMass_chi2_had_Up","reconstructed TPrime Mass (hadronic)",30,0.0,3000.);
  h_Ele_TMass_chi2_had_Down=new TH1D("Ele_TMass_chi2_had_Down","reconstructed TPrime Mass (hadronic)",30,0.0,3000.);
  h_Ele_TMass_chi2=new TH1D("Ele_TMass_chi2","reconstructed TPrime Mass",30,0.0,3000.);
  h_Ele_TMass_chi2_Up=new TH1D("Ele_TMass_chi2_Up","reconstructed TPrime Mass",30,0.0,3000.);
  h_Ele_TMass_chi2_Down=new TH1D("Ele_TMass_chi2_Down","reconstructed TPrime Mass",30,0.0,3000.);
  h_Ele_Chi2_cut=new TH1D("Ele_Chi2_c","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Ele_Chi2_cut_Up=new TH1D("Ele_Chi2_c_Up","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Ele_Chi2_cut_Down=new TH1D("Ele_Chi2_c_Down","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Ele_Chi2_1t_cut=new TH1D("Ele_Chi2_1t_c","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Ele_Chi2_1t_cut_Up=new TH1D("Ele_Chi2_1t_c_Up","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Ele_Chi2_1t_cut_Down=new TH1D("Ele_Chi2_1t_c_Down","Chi2 ( #chi^{2}<20.)",20,0.,40.);

  h_Ele_Chi2_0t1V_cut=new TH1D("Ele_Chi2_0t1V_c","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Ele_Chi2_0t1V_cut_Up=new TH1D("Ele_Chi2_0t1V_c_Up","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Ele_Chi2_0t1V_cut_Down=new TH1D("Ele_Chi2_0t1V_c_Down","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Ele_Chi2_0t2V_cut=new TH1D("Ele_Chi2_0t2V_c","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Ele_Chi2_0t2V_cut_Up=new TH1D("Ele_Chi2_0t2V_c_Up","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Ele_Chi2_0t2V_cut_Down=new TH1D("Ele_Chi2_0t2V_c_Down","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Ele_Chi2_1t1V_cut=new TH1D("Ele_Chi2_1t1V_c","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Ele_Chi2_1t1V_cut_Up=new TH1D("Ele_Chi2_1t1V_c_Up","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Ele_Chi2_1t1V_cut_Down=new TH1D("Ele_Chi2_1t1V_c_Down","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Ele_Chi2_1t0V_cut=new TH1D("Ele_Chi2_1t0V_c","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Ele_Chi2_1t0V_cut_Up=new TH1D("Ele_Chi2_1t0V_c_Up","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Ele_Chi2_1t0V_cut_Down=new TH1D("Ele_Chi2_1t0V_c_Down","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Ele_TMass_chi2_lep_cut=new TH1D("Ele_TMass_chi2_lep_c","reconstructed TPrime Mass (leptonic) ( #chi^{2}<20.)",30,0.0,3000.);
  h_Ele_TMass_chi2_lep_cut_Up=new TH1D("Ele_TMass_chi2_lep_c_Up","reconstructed TPrime Mass (leptonic) ( #chi^{2}<20.)",30,0.0,3000.);
  h_Ele_TMass_chi2_lep_cut_Down=new TH1D("Ele_TMass_chi2_lep_c_Down","reconstructed TPrime Mass (leptonic) ( #chi^{2}<20.)",30,0.0,3000.);
  h_Ele_TMass_chi2_had_cut=new TH1D("Ele_TMass_chi2_had_c","reconstructed TPrime Mass (hadronic) ( #chi^{2}<20.)",30,0.0,3000.);
  h_Ele_TMass_chi2_had_cut_Up=new TH1D("Ele_TMass_chi2_had_c_Up","reconstructed TPrime Mass (hadronic) ( #chi^{2}<20.)",30,0.0,3000.);
  h_Ele_TMass_chi2_had_cut_Down=new TH1D("Ele_TMass_chi2_had_c_Down","reconstructed TPrime Mass (hadronic) ( #chi^{2}<20.)",30,0.0,3000.);
  h_Ele_TMass_chi2_cut=new TH1D("Ele_TMass_chi2_c","reconstructed TPrime Mass ( #chi^{2}<20.)",30,0.0,3000.);
  h_Ele_TMass_chi2_cut_Up=new TH1D("Ele_TMass_chi2_c_Up","reconstructed TPrime Mass ( #chi^{2}<20.)",30,0.0,3000.);
  h_Ele_TMass_chi2_cut_Down=new TH1D("Ele_TMass_chi2_c_Down","reconstructed TPrime Mass ( #chi^{2}<20.)",30,0.0,3000.);


  //============================= HT Sel =============================================//
  oFile->mkdir("HTsel");
  oFile->cd("HTsel");
  h_LeadingAK4Pt_HT=new TH1D("LeadingAK4Pt","P_{T} for Leading AK4 Jet",50,0.0,2000.);
  h_NJetsAK4_HT=new TH1D("AK4Njets","Number of AK4 Jets",21,0.0,21);
  h_JetsAK4HT_HT=new TH1D("JetsAK4HT","H_{T} for AK4 Jets",100,0.0,6000.);
  h_JetsAK4Pt_HT=new TH1D("JetsAK4Pt","P_{T} for AK4 Jets",50,0.0,2000.);

  h_LeadingAK8Pt_HT=new TH1D("LeadingAK8Pt","P_{T} for Leading AK8 Jet",50,0.0,2000.);
  h_NJetsAK8_HT=new TH1D("AK8Njets","Number of AK8 Jets",21,0.0,21);
  h_JetsAK8HT_HT=new TH1D("JetsAK8HT","H_{T} for AK8 Jets",100,0.0,6000.);
  h_JetsAK8Pt_HT=new TH1D("JetsAK8Pt","P_{T} for AK8 Jets",50,0.0,2000.);

  
 
  h_Nlep_HT=new TH1D("Nlep","Number of reconstructed Muons+electrons beyond the dilepton pair",10,0.0,10);

  h_4jetEvent_HT=new TH1D("4jetEvent","Number of Events with >=4 AK4 Jets and at least one AK8 Jet",21,0.0,21);

  h_NrecobT_HT=new TH1D("NrecobT_HT","Number of reconstructed b quarks (CSV>0.9535)",10,0.0,10);
  h_Nrecob_HT=new TH1D("Nrecob_HT","Number of reconstructed b quarks (CSV>0.8484)",10,0.0,10);
  h_NrecobL_HT=new TH1D("NrecobL_HT","Number of reconstructed b quarks (CSV>0.5426)",10,0.0,10);

  h_bM=new TH1D("Nrecob_HT_M","Number of reconstructed b quarks (CSV>0.8484) (at least 1 loose and and 1 medium)",10,0.0,10);
  h_bL=new TH1D("Nrecob_HT_L","Number of reconstructed b quarks (CSV>0.5426) (at least 1 loose and and 1 medium)",10,0.0,10);

  h_NrecoTop_HT=new TH1D("NrecoTop","Number of reconstructed Tops",10,0.0,10);

  h_bmatrix=new TH2D("bmatrix","Distribution for two leading b's",3,0.0,3.,3,0.0,3.);
   //=========================== HT/ST check ===================================//
  oFile->mkdir("Check");
  oFile->cd("Check");
  h_HTcheck=new TH1D("HTAK4","H_{T} for AK4 Jets",120,0.0,6000.);
  h_STcheck=new TH1D("STAK4","S_{T} for AK4 Jets and leading dilepton pair(ee/ #mu#mu)",120,0.0,6000.);
  h_HTAK8check=new TH1D("HTAK8","H_{T} for AK8 Jets",120,0.0,6000.);
  h_STAK8check=new TH1D("STAK8","S_{T} for AK8 Jets and leading dilepton pair(ee/ #mu#mu)",120,0.0,6000.);
 
  //=========================== b veto ================================//
  oFile->mkdir("bVeto");
  oFile->cd("bVeto");
  h_ele_leadingLepPt=new TH1D("ele_l1Pt","P_{T} for leading electron",100,0.0,1000.);
  h_ele_leadingLep2Pt=new TH1D("ele_l2Pt","P_{T} for 2nd leading electron",100,0.0,1000.);
  h_ele_Lep2by1Pt=new TH1D("ele_l2byl1Pt","Lep 2 P_{T}/Lep 1 P_{T}",20,0.0,1.);
  h_ele_Lep2vs1Pt=new TH2D("ele_l2vsl1Pt","Lep 2 P_{T} vs Lep 1 P_{T}",20,0.0,200.,80,0.0,800.);
  h_ele_Lep2Ptfrac=new TH1D("ele_l2Ptfrac","Lep 2 P_{T}/(Lep 2 P_{T}+Lep 1 P_{T})",20,0.0,1.);
  h_mu_leadingLepPt=new TH1D("mu_l1Pt","P_{T} for leading muon",100,0.0,1000.);
  h_mu_leadingLep2Pt=new TH1D("mu_l2Pt","P_{T} for 2nd leading muon",100,0.0,1000.);
  h_mu_Lep2by1Pt=new TH1D("mu_l2byl1Pt","Lep 2 P_{T}/Lep 1 P_{T}",20,0.0,1.);
  h_mu_Lep2vs1Pt=new TH2D("mu_l2vsl1Pt","Lep 2 P_{T} vs Lep 1 P_{T}",20,0.0,200.,80,0.0,800.);
  h_mu_Lep2Ptfrac=new TH1D("mu_l2Ptfrac","Lep 2 P_{T}/(Lep 2 P_{T}+Lep 1 P_{T})",20,0.0,1.);
  h_dielept_HT=new TH1D("dielePt","P_{T} for dielectron",150,0.0,1500.);
  h_dimupt_HT=new TH1D("dimuPt","P_{T} for dimuon",150,0.0,1500.);

  h_tvsW=new TH2D("TopvsW","Number of Tops vs Number of W",10,0.0,10.,10,0.0,10.);
  h_tvsNonb=new TH2D("TopvsNonb","Number of Tops vs Number of AK4 non b jets",10,0.0,10.,10,0.0,10.);
  h_WvsNonb=new TH2D("WvsNonb","Number of W vs Number of AK4 non b jets",10,0.0,10.,10,0.0,10.);
  h_tvsWvsNonb=new TH3D("TopvsWvsNonb","Number of Tops vs Number of W vs Number of AK4 Non b jets",10,0.0,10.,10,0.0,10.,10,0.0,10.);
  h_TopWMix=new TH1D("TopWMix","Number of tops also reconstructed as W (or vice-versa)",10,0.0,10);
  h_NrecoH_bv=new TH1D("NrecoH","Number of reconstructed Higgs",10,0.0,10);
 

  h_tvsH=new TH2D("TopvsH","Number of Tops vs Number of H",10,0.0,10.,10,0.0,10.);

  h_1W_Dilep_dPhi=new TH1D("1W_Dilep_dPhi","#Delta#Phi between W+b and dilepton ( single W)",32,-3.2,3.2);
  h_LW_Dilep_dPhi=new TH1D("LW_Dilep_dPhi","#Delta#Phi between W+b and dilepton (multiple W)",32,-3.2,3.2);
  h_1Top_Dilep_dPhi=new TH1D("1Top_Dilep_dPhi","#Delta#Phi between top and dilepton (1 top)",32,-3.2,3.2);
  h_LTop_Dilep_dPhi=new TH1D("LTop_Dilep_dPhi","#Delta#Phi between top and dilepton (leading top)",32,-3.2,3.2);
  h_recoTMass1b=new TH1D("recoTMass1b","reconstructed TPrime Mass with only 1 medium b events",80,0.0,4000.);
  h_recoTMass2b=new TH1D("recoTMass2b","reconstructed TPrime Mass with at least 2 loose b",80,0.0,4000.);


  h_recoTMass_2top=new TH1D("recoTMass_2top","reconstructed TPrime Mass with events having 2 tops",80,0.0,4000.);
  h_recoTMass_1top1W=new TH1D("recoTMass_1top1W","reconstructed TPrime Mass with events having 1 top and 1 W",80,0.0,4000.);
  h_recoTMass_2W=new TH1D("recoTMass_2W","reconstructed TPrime Mass with events having 2 W",80,0.0,4000.);
  h_TMass=new TH1D("TMass","reconstructed TPrime Mass with events having 2 tops matched with gen",80,0.0,4000.);
  h_TMass_1W=new TH1D("TMass_1W","reconstructed TPrime Mass with events having 1 top and 1W matched with gen",80,0.0,4000.);
  h_TMass_2W=new TH1D("TMass_2W","reconstructed TPrime Mass with events having 2W matched with gen",80,0.0,4000.);
  h_topMass=new TH1D("topMass","Mass of AK8",35,0.,350.);
  h_topMassW=new TH1D("topMassW","reconstructed Mass with  W+b matched with gen",35,0.,350.);
  h_resTop=new TH1D("resTopMass","reconstructed Mass with  AK4+AK4+AK4 matched with gen",50,100.,350.);
  h_WMass=new TH1D("WMass","reconstructed V Mass matched with gen W",20,0.,140.);
  h_ZMass=new TH1D("ZMass","reconstructed V Mass matched with gen Z",20,0.,140.);
  h_genV=new TH1D("genVPt","P_{T} for gen Vjet",50,0.0,2000.);
  h_gentL=new TH1D("gent_LPt","P_{T} for gen top",50,0.0,2000.);
  h_genV2=new TH1D("genV2Pt","P_{T} for gen Vjet",50,0.0,2000.);
  h_topTerm=new TH1D("Chi2_topTerm","Top term",30,0.0,15.);
  h_WbTerm=new TH1D("Chi2_WbTerm","W+b term",30,0.0,15.);
  h_AK4WbTerm=new TH1D("Chi2_AK4WbTerm","AK4W+b term",30,0.0,15.);
  h_AK4WTerm=new TH1D("Chi2_AK4WTerm","W term",30,0.0,15.);
  h_ZTerm=new TH1D("Chi2_ZTerm","Z term",30,0.0,15.);
  h_AK4_Wb= new TH1D("AK4_Wb","Mass of AK4 +b",35,0.,350.);
  h_AKV= new TH1D("AKV","Mass of AK4 Vjet",35,0.,350.);
  h_resZ=new TH1D("resZ_dR","#delta R between possible Z daughters",40,0.,2.0);
  h_Nos=new TH1D("Nos","Number of Events surviving ",3,0.0,3);
  h_dRZ=new TH1D("dRZ","#delta R between dilepton candidates and gen Z",20,0.,1.0);


  
  h_Chi2=new TH1D("Chi2_MC","Chi2",60,0.,240.);
  h_Chi2_1t=new TH1D("Chi2_1t_MC","Chi2",60,0.,240.);
  h_Chi2_0t1V=new TH1D("Chi2_0t1V_MC","Chi2",60,0.,240.);
  h_Chi2_0t2V=new TH1D("Chi2_0t2V_MC","Chi2 ",60,0.,240.);
  h_Chi2_1t1V=new TH1D("Chi2_1t1V_MC","Chi2",60,0.,240.);
  h_Chi2_1t0V=new TH1D("Chi2_1t0V_MC","Chi2",60,0.,240.);
  
  h_TMass_chi2_lep=new TH1D("TMass_chi2_lep_MC","reconstructed TPrime Mass (leptonic)",30,0.0,3000.);
  h_TMass_chi2_had=new TH1D("TMass_chi2_had_MC","reconstructed TPrime Mass (hadronic)",30,0.0,3000.);
  h_TMass_chi2=new TH1D("TMass_chi2_MC","reconstructed TPrime Mass",30,0.0,3000.);
  h_Chi2_cut=new TH1D("Chi2_c_MC","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Chi2_1t_cut=new TH1D("Chi2_1t_c_MC","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  
  h_Chi2_0t1V_cut=new TH1D("Chi2_0t1V_c_MC","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Chi2_0t2V_cut=new TH1D("Chi2_0t2V_c_MC","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Chi2_1t1V_cut=new TH1D("Chi2_1t1V_c_MC","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_Chi2_1t0V_cut=new TH1D("Chi2_1t0V_c_MC","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_TMass_chi2_lep_cut=new TH1D("TMass_chi2_lep_c_MC","reconstructed TPrime Mass (leptonic) ( #chi^{2}<20.)",30,0.0,3000.);
  h_TMass_chi2_had_cut=new TH1D("TMass_chi2_had_c_MC","reconstructed TPrime Mass (hadronic) ( #chi^{2}<20.)",30,0.0,3000.);
  h_TMass_chi2_cut=new TH1D("TMass_chi2_c_MC","reconstructed TPrime Mass ( #chi^{2}<20.)",30,0.0,3000.);
  

  h_UB_Mu_Chi2=new TH1D("Chi2","Chi2",60,0.,240.);
  h_UB_Mu_Chi2_Up=new TH1D("Chi2_Up","Chi2",60,0.,240.);
  h_UB_Mu_Chi2_Down=new TH1D("Chi2_Down","Chi2",60,0.,240.);
  h_UB_Mu_Chi2_1t=new TH1D("Chi2_1t","Chi2",60,0.,240.);
  h_UB_Mu_Chi2_1t_Up=new TH1D("Chi2_1t_Up","Chi2",60,0.,240.);
  h_UB_Mu_Chi2_1t_Down=new TH1D("Chi2_1t_Down","Chi2",60,0.,240.);
  h_UB_Mu_Chi2_0t1V=new TH1D("Chi2_0t1V","Chi2",60,0.,240.);
  h_UB_Mu_Chi2_0t1V_Up=new TH1D("Chi2_0t1V_Up","Chi2",60,0.,240.);
  h_UB_Mu_Chi2_0t1V_Down=new TH1D("Chi2_0t1V_Down","Chi2",60,0.,240.);
  h_UB_Mu_Chi2_0t2V=new TH1D("Chi2_0t2V","Chi2 ",60,0.,240.);
  h_UB_Mu_Chi2_0t2V_Up=new TH1D("Chi2_0t2V_Up","Chi2 ",60,0.,240.);
  h_UB_Mu_Chi2_0t2V_Down=new TH1D("Chi2_0t2V_Down","Chi2 ",60,0.,240.);
  h_UB_Mu_Chi2_1t1V=new TH1D("Chi2_1t1V","Chi2",60,0.,240.);
  h_UB_Mu_Chi2_1t1V_Up=new TH1D("Chi2_1t1V_Up","Chi2",60,0.,240.);
  h_UB_Mu_Chi2_1t1V_Down=new TH1D("Chi2_1t1V_Down","Chi2",60,0.,240.);
  h_UB_Mu_Chi2_1t0V=new TH1D("Chi2_1t0V","Chi2",60,0.,240.);
  h_UB_Mu_Chi2_1t0V_Up=new TH1D("Chi2_1t0V_Up","Chi2",60,0.,240.);
  h_UB_Mu_Chi2_1t0V_Down=new TH1D("Chi2_1t0V_Down","Chi2",60,0.,240.);
  
  h_UB_Mu_TMass_chi2_lep=new TH1D("TMass_chi2_lep","reconstructed TPrime Mass (leptonic)",30,0.0,3000.);
  h_UB_Mu_TMass_chi2_lep_Up=new TH1D("TMass_chi2_lep_Up","reconstructed TPrime Mass (leptonic)",30,0.0,3000.);
  h_UB_Mu_TMass_chi2_lep_Down=new TH1D("TMass_chi2_lep_Down","reconstructed TPrime Mass (leptonic)",30,0.0,3000.);
  h_UB_Mu_TMass_chi2_had=new TH1D("TMass_chi2_had","reconstructed TPrime Mass (hadronic)",30,0.0,3000.);
  h_UB_Mu_TMass_chi2_had_Up=new TH1D("TMass_chi2_had_Up","reconstructed TPrime Mass (hadronic)",30,0.0,3000.);
  h_UB_Mu_TMass_chi2_had_Down=new TH1D("TMass_chi2_had_Down","reconstructed TPrime Mass (hadronic)",30,0.0,3000.);
  h_UB_Mu_TMass_chi2=new TH1D("TMass_chi2","reconstructed TPrime Mass",30,0.0,3000.);
  h_UB_Mu_TMass_chi2_Up=new TH1D("TMass_chi2_Up","reconstructed TPrime Mass",30,0.0,3000.);
  h_UB_Mu_TMass_chi2_Down=new TH1D("TMass_chi2_Down","reconstructed TPrime Mass",30,0.0,3000.);
  h_UB_Mu_Chi2_cut=new TH1D("Chi2_c","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Mu_Chi2_cut_Up=new TH1D("Chi2_c_Up","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Mu_Chi2_cut_Down=new TH1D("Chi2_c_Down","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Mu_Chi2_1t_cut=new TH1D("Chi2_1t_c","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Mu_Chi2_1t_cut_Up=new TH1D("Chi2_1t_c_Up","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Mu_Chi2_1t_cut_Down=new TH1D("Chi2_1t_c_Down","Chi2 ( #chi^{2}<20.)",20,0.,40.);

  h_UB_Mu_Chi2_0t1V_cut=new TH1D("Chi2_0t1V_c","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Mu_Chi2_0t1V_cut_Up=new TH1D("Chi2_0t1V_c_Up","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Mu_Chi2_0t1V_cut_Down=new TH1D("Chi2_0t1V_c_Down","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Mu_Chi2_0t2V_cut=new TH1D("Chi2_0t2V_c","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Mu_Chi2_0t2V_cut_Up=new TH1D("Chi2_0t2V_c_Up","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Mu_Chi2_0t2V_cut_Down=new TH1D("Chi2_0t2V_c_Down","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Mu_Chi2_1t1V_cut=new TH1D("Chi2_1t1V_c","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Mu_Chi2_1t1V_cut_Up=new TH1D("Chi2_1t1V_c_Up","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Mu_Chi2_1t1V_cut_Down=new TH1D("Chi2_1t1V_c_Down","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Mu_Chi2_1t0V_cut=new TH1D("Chi2_1t0V_c","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Mu_Chi2_1t0V_cut_Up=new TH1D("Chi2_1t0V_c_Up","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Mu_Chi2_1t0V_cut_Down=new TH1D("Chi2_1t0V_c_Down","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Mu_TMass_chi2_lep_cut=new TH1D("TMass_chi2_lep_c","reconstructed TPrime Mass (leptonic) ( #chi^{2}<20.)",30,0.0,3000.);
  h_UB_Mu_TMass_chi2_lep_cut_Up=new TH1D("TMass_chi2_lep_c_Up","reconstructed TPrime Mass (leptonic) ( #chi^{2}<20.)",30,0.0,3000.);
  h_UB_Mu_TMass_chi2_lep_cut_Down=new TH1D("TMass_chi2_lep_c_Down","reconstructed TPrime Mass (leptonic) ( #chi^{2}<20.)",30,0.0,3000.);
  h_UB_Mu_TMass_chi2_had_cut=new TH1D("TMass_chi2_had_c","reconstructed TPrime Mass (hadronic) ( #chi^{2}<20.)",30,0.0,3000.);
  h_UB_Mu_TMass_chi2_had_cut_Up=new TH1D("TMass_chi2_had_c_Up","reconstructed TPrime Mass (hadronic) ( #chi^{2}<20.)",30,0.0,3000.);
  h_UB_Mu_TMass_chi2_had_cut_Down=new TH1D("TMass_chi2_had_c_Down","reconstructed TPrime Mass (hadronic) ( #chi^{2}<20.)",30,0.0,3000.);
  h_UB_Mu_TMass_chi2_cut=new TH1D("TMass_chi2_c","reconstructed TPrime Mass ( #chi^{2}<20.)",30,0.0,3000.);
  h_UB_Mu_TMass_chi2_cut_Up=new TH1D("TMass_chi2_c_Up","reconstructed TPrime Mass ( #chi^{2}<20.)",30,0.0,3000.);
  h_UB_Mu_TMass_chi2_cut_Down=new TH1D("TMass_chi2_c_Down","reconstructed TPrime Mass ( #chi^{2}<20.)",30,0.0,3000.);
  
  //xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx//
  
  h_UB_Ele_Chi2=new TH1D("Ele_Chi2","Chi2",60,0.,240.);
  h_UB_Ele_Chi2_Up=new TH1D("Ele_Chi2_Up","Chi2",60,0.,240.);
  h_UB_Ele_Chi2_Down=new TH1D("Ele_Chi2_Down","Chi2",60,0.,240.);
  h_UB_Ele_Chi2_1t=new TH1D("Ele_Chi2_1t","Chi2",60,0.,240.);
  h_UB_Ele_Chi2_1t_Up=new TH1D("Ele_Chi2_1t_Up","Chi2",60,0.,240.);
  h_UB_Ele_Chi2_1t_Down=new TH1D("Ele_Chi2_1t_Down","Chi2",60,0.,240.);
  h_UB_Ele_Chi2_0t1V=new TH1D("Ele_Chi2_0t1V","Chi2",60,0.,240.);
  h_UB_Ele_Chi2_0t1V_Up=new TH1D("Ele_Chi2_0t1V_Up","Chi2",60,0.,240.);
  h_UB_Ele_Chi2_0t1V_Down=new TH1D("Ele_Chi2_0t1V_Down","Chi2",60,0.,240.);
  h_UB_Ele_Chi2_0t2V=new TH1D("Ele_Chi2_0t2V","Chi2 ",60,0.,240.);
  h_UB_Ele_Chi2_0t2V_Up=new TH1D("Ele_Chi2_0t2V_Up","Chi2 ",60,0.,240.);
  h_UB_Ele_Chi2_0t2V_Down=new TH1D("Ele_Chi2_0t2V_Down","Chi2 ",60,0.,240.);
  h_UB_Ele_Chi2_1t1V=new TH1D("Ele_Chi2_1t1V","Chi2",60,0.,240.);
  h_UB_Ele_Chi2_1t1V_Up=new TH1D("Ele_Chi2_1t1V_Up","Chi2",60,0.,240.);
  h_UB_Ele_Chi2_1t1V_Down=new TH1D("Ele_Chi2_1t1V_Down","Chi2",60,0.,240.);
  h_UB_Ele_Chi2_1t0V=new TH1D("Ele_Chi2_1t0V","Chi2",60,0.,240.);
  h_UB_Ele_Chi2_1t0V_Up=new TH1D("Ele_Chi2_1t0V_Up","Chi2",60,0.,240.);
  h_UB_Ele_Chi2_1t0V_Down=new TH1D("Ele_Chi2_1t0V_Down","Chi2",60,0.,240.);
  
  h_UB_Ele_TMass_chi2_lep=new TH1D("Ele_TMass_chi2_lep","reconstructed TPrime Mass (leptonic)",30,0.0,3000.);
  h_UB_Ele_TMass_chi2_lep_Up=new TH1D("Ele_TMass_chi2_lep_Up","reconstructed TPrime Mass (leptonic)",30,0.0,3000.);
  h_UB_Ele_TMass_chi2_lep_Down=new TH1D("Ele_TMass_chi2_lep_Down","reconstructed TPrime Mass (leptonic)",30,0.0,3000.);
  h_UB_Ele_TMass_chi2_had=new TH1D("Ele_TMass_chi2_had","reconstructed TPrime Mass (hadronic)",30,0.0,3000.);
  h_UB_Ele_TMass_chi2_had_Up=new TH1D("Ele_TMass_chi2_had_Up","reconstructed TPrime Mass (hadronic)",30,0.0,3000.);
  h_UB_Ele_TMass_chi2_had_Down=new TH1D("Ele_TMass_chi2_had_Down","reconstructed TPrime Mass (hadronic)",30,0.0,3000.);
  h_UB_Ele_TMass_chi2=new TH1D("Ele_TMass_chi2","reconstructed TPrime Mass",30,0.0,3000.);
  h_UB_Ele_TMass_chi2_Up=new TH1D("Ele_TMass_chi2_Up","reconstructed TPrime Mass",30,0.0,3000.);
  h_UB_Ele_TMass_chi2_Down=new TH1D("Ele_TMass_chi2_Down","reconstructed TPrime Mass",30,0.0,3000.);
  h_UB_Ele_Chi2_cut=new TH1D("Ele_Chi2_c","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Ele_Chi2_cut_Up=new TH1D("Ele_Chi2_c_Up","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Ele_Chi2_cut_Down=new TH1D("Ele_Chi2_c_Down","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Ele_Chi2_1t_cut=new TH1D("Ele_Chi2_1t_c","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Ele_Chi2_1t_cut_Up=new TH1D("Ele_Chi2_1t_c_Up","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Ele_Chi2_1t_cut_Down=new TH1D("Ele_Chi2_1t_c_Down","Chi2 ( #chi^{2}<20.)",20,0.,40.);

  h_UB_Ele_Chi2_0t1V_cut=new TH1D("Ele_Chi2_0t1V_c","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Ele_Chi2_0t1V_cut_Up=new TH1D("Ele_Chi2_0t1V_c_Up","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Ele_Chi2_0t1V_cut_Down=new TH1D("Ele_Chi2_0t1V_c_Down","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Ele_Chi2_0t2V_cut=new TH1D("Ele_Chi2_0t2V_c","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Ele_Chi2_0t2V_cut_Up=new TH1D("Ele_Chi2_0t2V_c_Up","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Ele_Chi2_0t2V_cut_Down=new TH1D("Ele_Chi2_0t2V_c_Down","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Ele_Chi2_1t1V_cut=new TH1D("Ele_Chi2_1t1V_c","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Ele_Chi2_1t1V_cut_Up=new TH1D("Ele_Chi2_1t1V_c_Up","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Ele_Chi2_1t1V_cut_Down=new TH1D("Ele_Chi2_1t1V_c_Down","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Ele_Chi2_1t0V_cut=new TH1D("Ele_Chi2_1t0V_c","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Ele_Chi2_1t0V_cut_Up=new TH1D("Ele_Chi2_1t0V_c_Up","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Ele_Chi2_1t0V_cut_Down=new TH1D("Ele_Chi2_1t0V_c_Down","Chi2 ( #chi^{2}<20.)",20,0.,40.);
  h_UB_Ele_TMass_chi2_lep_cut=new TH1D("Ele_TMass_chi2_lep_c","reconstructed TPrime Mass (leptonic) ( #chi^{2}<20.)",30,0.0,3000.);
  h_UB_Ele_TMass_chi2_lep_cut_Up=new TH1D("Ele_TMass_chi2_lep_c_Up","reconstructed TPrime Mass (leptonic) ( #chi^{2}<20.)",30,0.0,3000.);
  h_UB_Ele_TMass_chi2_lep_cut_Down=new TH1D("Ele_TMass_chi2_lep_c_Down","reconstructed TPrime Mass (leptonic) ( #chi^{2}<20.)",30,0.0,3000.);
  h_UB_Ele_TMass_chi2_had_cut=new TH1D("Ele_TMass_chi2_had_c","reconstructed TPrime Mass (hadronic) ( #chi^{2}<20.)",30,0.0,3000.);
  h_UB_Ele_TMass_chi2_had_cut_Up=new TH1D("Ele_TMass_chi2_had_c_Up","reconstructed TPrime Mass (hadronic) ( #chi^{2}<20.)",30,0.0,3000.);
  h_UB_Ele_TMass_chi2_had_cut_Down=new TH1D("Ele_TMass_chi2_had_c_Down","reconstructed TPrime Mass (hadronic) ( #chi^{2}<20.)",30,0.0,3000.);
  h_UB_Ele_TMass_chi2_cut=new TH1D("Ele_TMass_chi2_c","reconstructed TPrime Mass ( #chi^{2}<20.)",30,0.0,3000.);
  h_UB_Ele_TMass_chi2_cut_Up=new TH1D("Ele_TMass_chi2_c_Up","reconstructed TPrime Mass ( #chi^{2}<20.)",30,0.0,3000.);
  h_UB_Ele_TMass_chi2_cut_Down=new TH1D("Ele_TMass_chi2_c_Down","reconstructed TPrime Mass ( #chi^{2}<20.)",30,0.0,3000.);

  h_Nb2t=new TH1D("Nrecob_2top","Number of reconstructed loose b quarks",10,0.0,10);
  h_NW2t=new TH1D("NrecoW_2top","Number of reconstructed W",10,0.0,10);
  h_Jeta=new TH1D("Jeta","#eta for AK8 jets",40,-4.0,4.0);
  h_Jphi=new TH1D("Jphi","Phi for jets",32,-3.2,3.2);
  h_Jpt=new TH1D("AK8pt","P_{T} for AK8 jets",40,0.0,2000.);
  h_dR_JZ=new TH1D("dRJZ","#delta R between dilepton candidates and AK8 jets",40,0.,4.0);
  h_nAK8=new TH1D("nAK8","Number of AK8 jets with Pt>180. and |#eta|<2.4",10,0.0,10);
   h_topTag=new TH1D("topTag","Number of top jets",10,0.0,10);
  h_Vtag=new TH1D("VTag","Number of V jets",10,0.0,10);
  h_AK8vsb=new TH2D("AK8 vs b","Number of AK8 jets vs number of loose b jets",10,0.0,10., 10,0.,10.);
  h_Topvsb=new TH2D("NrecoTopvsb","Number of reco tops vs number of loose b jets",10,0.0,10., 10,0.,10.);
  h_TopvsAK4=new TH2D("NrecoTopvsAK4","Number of reco tops vs number of AK4 jets",10,0.0,10., 10,0.,10.);
  h_Vvsb=new TH2D("NrecoVvsb","Number of reco V jets vs number of loose b jets",10,0.0,10., 10,0.,10.);
  h_VvsAK4=new TH2D("NrecoVvsAK4","Number of reco V jets vs number of AK4 jets",10,0.0,10., 10,0.,10.);
  h_Jneweta=new TH1D("newJeta","#eta for new AK8 jets",40,-4.0,4.0);
  h_Jnewphi=new TH1D("newJphi","Phi for new jets",32,-3.2,3.2);
  h_Jnewpt=new TH1D("newAK8pt","P_{T} for new AK8 jets",40,0.0,2000.);
  h_ratioPt=new TH1D("ratio","Ratio of new to old P_{T}",20,0.0,1.);
  }

TPrimeAna::TPrimeAna(const TString &inputFileList, const char *outFileName, const char* dataset, const char *isData)
{
TChain *tree = new TChain("ana/tree");

  if( ! FillChain(tree, inputFileList) ) {
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;
    char temp[]="T";

    if(strcmp(temp,isData)==0)std::cout<<"Initiating analysis on Data"<<endl;
    else std::cout<<"Initiating analysis on MC"<<endl;
  }
  
  NtupleVariables::Init(tree);
  BookHistogram(outFileName);
  }

Bool_t TPrimeAna::FillChain(TChain *chain, const TString &inputFileList) {

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
TPrimeAna::~TPrimeAna()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
  oFile->cd();
  oFile->Write();
  oFile->Close();
}


Long64_t TPrimeAna::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
  t_elecP4 = getTLV(0);
  t_muonP4 = getTLV(1);
  t_genPartP4 = getTLV(2);
  t_jetAK4P4 = getTLV(3);
  t_jetAK8P4 = getTLV(4);
  t_metP4 = getTLV(5);
  t_jetTopJetP4 = getTLV(6);
  t_jetWJetP4 = getTLV(7);
  t_ZllP4 = getTLV(8);
   if (centry < 0) return centry;
    if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      // Notify();
   }
   return centry;
}

vector <TLorentzVector>* TPrimeAna::getTLV( int part){
   
  TLorentzVector temp;
  if( part==0){
    std::vector<TLorentzVector> *elP4 = new std::vector <TLorentzVector>();
    for(unsigned int i=0; i<t_elPt->size(); i++) {
      //cout<<"el ok\n";
      temp.SetPtEtaPhiE((*t_elPt)[i],(*t_elEta)[i],(*t_elPhi)[i],(*t_elE)[i]);
      elP4->push_back(temp);
    }
    return elP4;
  }
  else if(part==1){
    std::vector<TLorentzVector> *muP4 = new std::vector <TLorentzVector>();
    for(unsigned int i=0; i<t_muPt->size(); i++) {
      //cout<<"mu ok\n";
      temp.SetPtEtaPhiE((*t_muPt)[i],(*t_muEta)[i],(*t_muPhi)[i],(*t_muE)[i]);
      muP4->push_back(temp);
    }
    return muP4;
  }
  else if(part==2){
    std::vector<TLorentzVector> *genP4 = new std::vector <TLorentzVector>();
    //cout<<"damn!"<<endl;
    for(unsigned int i=0; i<t_genPartID->size(); i++) {
      //cout<<(*t_genPartPt)[i]<<endl;
      temp.SetPtEtaPhiE((*t_genPartPt)[i],(*t_genPartEta)[i],(*t_genPartPhi)[i],(*t_genPartE)[i]);
      genP4->push_back(temp);
    }
    return genP4;
  }
  else if(part==3){
    std::vector<TLorentzVector> *AK4P4 = new std::vector <TLorentzVector>();
    for(unsigned int i=0; i<t_jetAK4Pt->size(); i++) {
      temp.SetPtEtaPhiE((*t_jetAK4Pt)[i],(*t_jetAK4Eta)[i],(*t_jetAK4Phi)[i],(*t_jetAK4E)[i]);
      AK4P4->push_back(temp);
    }
    return AK4P4;
  }
  else if(part==4){
    std::vector<TLorentzVector> *AK8P4 = new std::vector <TLorentzVector>();
    for(unsigned int i=0; i<t_jetAK8Pt->size(); i++) {
      temp.SetPtEtaPhiE((*t_jetAK8Pt)[i],(*t_jetAK8Eta)[i],(*t_jetAK8Phi)[i],(*t_jetAK8E)[i]);
      AK8P4->push_back(temp);
    }
    return AK8P4;
  }
  else if(part==5){
    std::vector<TLorentzVector> *metP4 = new std::vector <TLorentzVector>();
    for(unsigned int i=0; i<t_metPt->size(); i++) {
      temp.SetPtEtaPhiE((*t_metPt)[i],(*t_metEta)[i],(*t_metPhi)[i],(*t_metE)[i]);
      metP4->push_back(temp);
    }
    return metP4;
  }
  else if(part==6){
    std::vector<TLorentzVector> *jetTopJetP4 = new std::vector <TLorentzVector>();
    for(unsigned int i=0; i<t_jetTopJetPt->size(); i++) {
      temp.SetPtEtaPhiE((*t_jetTopJetPt)[i],(*t_jetTopJetEta)[i],(*t_jetTopJetPhi)[i],(*t_jetTopJetE)[i]);
      jetTopJetP4->push_back(temp);
    }
    return jetTopJetP4;
  }

  else if(part==7){
    std::vector<TLorentzVector> *jetWJetP4 = new std::vector <TLorentzVector>();
    for(unsigned int i=0; i<t_jetWJetPt->size(); i++) {
      temp.SetPtEtaPhiE((*t_jetWJetPt)[i],(*t_jetWJetEta)[i],(*t_jetWJetPhi)[i],(*t_jetWJetE)[i]);
      jetWJetP4->push_back(temp);
    }
    return jetWJetP4;
  }

   else if(part==8){
    std::vector<TLorentzVector> *ZllP4 = new std::vector <TLorentzVector>();
    for(unsigned int i=0; i<t_ZllPt->size(); i++) {
      temp.SetPtEtaPhiE((*t_ZllPt)[i],(*t_ZllEta)[i],(*t_ZllPhi)[i],(*t_ZllE)[i]);
      ZllP4->push_back(temp);
    }
    return ZllP4;
  }
}
double TPrimeAna::HTCleanAK4(vector <int> cleanJetid) {
  double ht = 0.0;
  for (vector<int>::iterator it = cleanJetid.begin() ; it != cleanJetid.end(); ++it){
    ht+=(*t_jetAK4P4)[*it].Pt() ;
  }
  return ht;
}
double TPrimeAna::HTCleanAK8(vector <int> cleanJetid) {
  double ht = 0.0;
  for (vector<int>::iterator it = cleanJetid.begin() ; it != cleanJetid.end(); ++it){
    ht+=(*t_jetAK8P4)[*it].Pt() ;
  }
  return ht;
}


vector <int> TPrimeAna::sort(unsigned int lepctr,unsigned int *lep )
{
  bool swapped = true;
      int j = 0;
      int tmp;
      while (swapped) {
            swapped = false;
            j++;
            for (int i = 0; i < lepctr - j; i++) {
	      if ((*t_genPartP4)[lep[i+1]].Pt() > (*t_genPartP4)[lep[i]].Pt()) {
                        tmp = lep[i];
                        lep[i] = lep[i + 1];
                        lep[i + 1] = tmp;
                        swapped = true;
                  }
            }
      }
      vector <int> lepid;
      for(unsigned int i=0;i<lepctr;i++){
	lepid.push_back(lep[i]);
      }
      return lepid;
}
vector <int> TPrimeAna::selGoodAK4Jet(vector <int > elid, vector <int> muid)
{
  vector <int> badjet,goodjet;
  float eta1,phi1,eta2,phi2,dR,mindR;
  int jetid,flag;
  for(vector<int>::iterator it = elid.begin() ; it != elid.end(); ++it){
	 eta2=(*t_elecP4)[*it].Eta();
	 phi2=(*t_elecP4)[*it].Phi();
	 jetid=-99;
	 dR=100;
	 mindR=0.4;
	 for(unsigned int j=0; j<t_jetAK4P4->size(); j++) {
	   eta1=(*t_jetAK4P4)[j].Eta();
	   phi1=(*t_jetAK4P4)[j].Phi();
	   dR=DeltaR(eta1,phi1,eta2,phi2);
	   if(dR<0.4 && dR<mindR){
	     jetid=j;
	     mindR=dR;
	   }
	 }
	 if(jetid!=-99)badjet.push_back(jetid);
       }
    for(vector<int>::iterator it = muid.begin() ; it != muid.end(); ++it){
	 eta2=(*t_muonP4)[*it].Eta();
	 phi2=(*t_muonP4)[*it].Phi();
	 jetid=-99;
	 dR=100;
	 mindR=0.4;
	 for(unsigned int j=0; j<t_jetAK4P4->size(); j++) {
	   eta1=(*t_jetAK4P4)[j].Eta();
	   phi1=(*t_jetAK4P4)[j].Phi();
	   dR=DeltaR(eta1,phi1,eta2,phi2);
	    if(dR<0.4 && dR<mindR){
	     jetid=j;
	     mindR=dR;
	   }
	   }
	if(jetid!=-99)badjet.push_back(jetid);		    
       }
   for(unsigned int j=0; j<t_jetAK4P4->size(); j++) {
     flag=0;
     for (vector<int>::iterator it = badjet.begin() ; it != badjet.end(); ++it)
       {
	 if(*it==j){flag=1;break;}
       }
     if(!flag && (*t_jetAK4P4)[j].Pt()>30. && fabs((*t_jetAK4P4)[j].Eta())<2.4) goodjet.push_back(j);
   }
   return goodjet;
}
void TPrimeAna::LepReco(int e1,int e2,int type, vector <int> lepid){
  int itEl=99; int itMu=99;
  double phi1,eta1,phi2,eta2,dR;
  double drMin=999.0; int iimin=-1;
  if(type==0){
    
    drMin=999.0; iimin=-1;
    for( vector<int>::iterator i=lepid.begin(); i!=lepid.end(); i++) {
      eta1=(*t_elecP4)[*i].Eta();
      phi1=(*t_elecP4)[*i].Phi();
      eta2=(*t_genPartP4)[e1].Eta();
      phi2=(*t_genPartP4)[e1].Phi();
      dR= DeltaR(eta1,phi1,eta2,phi2);
      if(dR<drMin /*&& ((*t_elCharge)[i]*(*t_genPartID)[e1]>0)*/){ // it isn't important to find the best matching candidate. As long as I find even one electron within dr<0.4 of the gen electron, I know that the gen will get matched to a reco and thus straightaway fill the gen electron pt
	//h_genMatchedEle1Pt->Fill((*t_genPartP4)[e1].Pt());
	//break;
	drMin=dR; iimin=*i;
      }
    }
    if(drMin<0.4 && iimin>=0)
      h_genMatchedEle1Pt->Fill((*t_genPartP4)[e1].Pt());
    drMin=999.0;iimin=-1;
    for( vector<int>::iterator i=lepid.begin(); i!=lepid.end(); i++) {
      eta1=(*t_elecP4)[*i].Eta();
      phi1=(*t_elecP4)[*i].Phi();
      eta2=(*t_genPartP4)[e2].Eta();
      phi2=(*t_genPartP4)[e2].Phi();
      dR= DeltaR(eta1,phi1,eta2,phi2);
      if(dR<drMin/* && ((*t_elCharge)[i]*(*t_genPartID)[e2]>0)*/){   
	drMin=dR; iimin=*i;
      }
    }
    if(drMin<0.4 && iimin>=0)
      h_genMatchedEle2Pt->Fill((*t_genPartP4)[e2].Pt());
  }
  else{
    for( vector<int>::iterator i=lepid.begin(); i!=lepid.end(); i++) {
      eta1=(*t_muonP4)[*i].Eta();
      phi1=(*t_muonP4)[*i].Phi();
      eta2=(*t_genPartP4)[e1].Eta();
      phi2=(*t_genPartP4)[e1].Phi();
      dR= DeltaR(eta1,phi1,eta2,phi2);
      if(dR<drMin /*&& ((*t_muCharge)[i]*(*t_genPartID)[e1]>0)*/){ // it isn't important to find the best matching candidate. As long as I find even one electron within dr<0.4 of the gen electron, I know that the gen will get matched to a reco and thus straightaway fill the gen electron pt 
	drMin=dR; iimin=*i;
      }
    }
    if(drMin<0.4 && iimin>=0)
      h_genMatchedMu1Pt->Fill((*t_genPartP4)[e1].Pt());
    drMin=999.0;iimin=-1;
    for( vector<int>::iterator i=lepid.begin(); i!=lepid.end(); i++) {
      eta1=(*t_muonP4)[*i].Eta();
      phi1=(*t_muonP4)[*i].Phi();
      eta2=(*t_genPartP4)[e2].Eta();
      phi2=(*t_genPartP4)[e2].Phi();
      dR= DeltaR(eta1,phi1,eta2,phi2);
      if(dR<drMin /*&& ((*t_muCharge)[i]*(*t_genPartID)[e2]>0)*/){   
	drMin=dR; iimin=*i;
      }
    }
    if(drMin<0.4 && iimin>=0)
      h_genMatchedMu2Pt->Fill((*t_genPartP4)[e2].Pt());
  }
}
void TPrimeAna::ZFinder(int iZ,int type,vector <int> lepid)
{
  int e1=99;int e2=99;
  int mu1=99;int mu2=99;
  double dR,phi1,eta1,phi2,eta2;
  int found1=0;
  if(type==0){
    for(unsigned int j=iZ+1; j<t_genPartP4->size(); j++) {
      if((fabs((*t_genPartID)[j])==11) && ((*t_genPartMom1ID)[j]==23 || (*t_genPartMom2ID)[j]==23)){
	if(e1==99)e1=j;
	else e2=j;
      }
    }
    //cout<<iZ<<" "<<e1<<" "<<e2<<endl; 
    found1=0;
    for( vector<int>::iterator i=lepid.begin(); i!=lepid.end(); i++) {
      //cout<<"vla\n";
      eta1=(*t_elecP4)[*i].Eta();
      phi1=(*t_elecP4)[*i].Phi();
      cout<<"vlavla\n";
      eta2=(*t_genPartP4)[e1].Eta();
      phi2=(*t_genPartP4)[e1].Phi();
      //cout<<"meh\n";
      dR= DeltaR(eta1,phi1,eta2,phi2);
      if(dR<0.4 /*&& ((*t_elCharge)[i]*(*t_genPartID)[e1]>0)*/){
 	found1=1;
	break;
      }
    }
    if(found1){
      for( vector<int>::iterator i=lepid.begin(); i!=lepid.end(); i++) {
	eta1=(*t_elecP4)[*i].Eta();
	phi1=(*t_elecP4)[*i].Phi();
	eta2=(*t_genPartP4)[e2].Eta();
	phi2=(*t_genPartP4)[e2].Phi();
	dR= DeltaR(eta1,phi1,eta2,phi2);
	if(dR<0.4 /*&& ((*t_elCharge)[i]*(*t_genPartID)[e2]>0)*/){
	  h_genMZElPt->Fill((*t_genPartP4)[iZ].Pt());
	  break;
	}
      }
    }
  }
  else if(type==1){
    for(unsigned int j=iZ+1; j<t_genPartP4->size(); j++) {
      if((fabs((*t_genPartID)[j])==13) && ((*t_genPartMom1ID)[j]==23 || (*t_genPartMom2ID)[j]==23)){
	if(mu1==99)mu1=j;
	else mu2=j;
      }
    }
    //cout<<iZ<<" "<<e1<<" "<<e2<<endl; 
    found1=0;
    for( vector<int>::iterator i=lepid.begin(); i!=lepid.end(); i++) {
      eta1=(*t_muonP4)[*i].Eta();
      phi1=(*t_muonP4)[*i].Phi();
      eta2=(*t_genPartP4)[mu1].Eta();
      phi2=(*t_genPartP4)[mu1].Phi();
      dR= DeltaR(eta1,phi1,eta2,phi2);
      if(dR<0.4 /*&& ((*t_muCharge)[i]*(*t_genPartID)[mu1]>0)*/){
	found1=1;
	break;
      }
    }
    if(found1){
      for( vector<int>::iterator i=lepid.begin(); i!=lepid.end(); i++) {
	eta1=(*t_muonP4)[*i].Eta();
	phi1=(*t_muonP4)[*i].Phi();
	eta2=(*t_genPartP4)[mu2].Eta();
	phi2=(*t_genPartP4)[mu2].Phi();
	dR= DeltaR(eta1,phi1,eta2,phi2);
	if(dR<0.4 /*&& ((*t_muCharge)[i]*(*t_genPartID)[mu2]>0)*/){
	  h_genMZMuPt->Fill((*t_genPartP4)[iZ].Pt());
	  break;
	}
      }
    }
	
  }
}

double TPrimeAna::chi_cal(int AK1, int AK2, int AK3, int  AK4, int AK5)
{
  TLorentzVector Wb1,V,W;
  float chi[3];
  for(int i=0;i<3;i++)chi[i]=10000;
  V=(*t_jetAK4P4)[AK1]+(*t_jetAK4P4)[AK2];
  //========= case 0 ==================                                                                                                      
  W=(*t_jetAK4P4)[AK3]+(*t_jetAK4P4)[AK4];
  if(W.M()<105. && W.M()>65.){
    Wb1=W+(*t_jetAK4P4)[AK5];
    chi[0]=pow(((*t_jetTopJet_SoftDropMass)[0]-179.7)/(18.4),2.)+pow((Wb1.M()-172)/(17.6),2.)+pow((V.M()-88.5)/(9.5),2.);
  }

  //========= case 1 ==================                                                                                                        
  W=(*t_jetAK4P4)[AK3]+(*t_jetAK4P4)[AK5];
  if(W.M()<105. && W.M()>65.){
    Wb1=W+(*t_jetAK4P4)[AK4];
    chi[1]=pow(((*t_jetTopJet_SoftDropMass)[0]-179.7)/(18.4),2.)+pow((Wb1.M()-172)/(17.6),2.)+pow((V.M()-88.5)/(9.5),2.);
  }


  //========= case 2 ==================                                                                                                        
  W=(*t_jetAK4P4)[AK5]+(*t_jetAK4P4)[AK4];
  if(W.M()<105. && W.M()>65.){
    Wb1=W+(*t_jetAK4P4)[AK3];
    chi[2]=pow(((*t_jetTopJet_SoftDropMass)[0]-179.7)/(18.4),2.)+pow((Wb1.M()-172)/(17.6),2.)+pow((V.M()-88.5)/(9.5),2.);
  }
  //cout<<chi[0]<<" "<<chi[1]<<" "<<chi[2]<<endl;
  if(chi[0]<chi[1] && chi[0]<chi[2]) return chi[0];
  if(chi[1]<chi[0] && chi[1]<chi[2]) return chi[1];
  if(chi[2]<chi[1] && chi[2]<chi[0]) return chi[2];
  return 10000;
}
vector <float> TPrimeAna::chi_calculator(int bid, int AK2, int AK3, int  AK4, int AK5, int AK6, TLorentzVector v1 )
{
  TLorentzVector Wb1,V,Wb2,W,V1,V2;
  float chi[10],chi_1;
  for(int i=0;i<10;i++)chi[i]=10000;
  vector <float> res;
  int chi_idx[]={0,1,2,3,4,5,6,7,8,9};
  int type[]={0,0,0,0,0,0,0,0,0,0};
  //Wb1= (*t_jetWJetP4)[0]+(*t_jetAK4P4)[bid];

    //========= case 0 ==================
  if((*t_jetAK4CSV)[AK2]>0.5426 || (*t_jetAK4CSV)[AK3]>0.5426 || (*t_jetAK4CSV)[AK4]>0.5426){
    Wb1= (*t_jetWJetP4)[0]+(*t_jetAK4P4)[bid];
    W=(*t_jetAK4P4)[AK2]+(*t_jetAK4P4)[AK3];
    V1=(*t_jetAK4P4)[AK2]+(*t_jetAK4P4)[AK4];
    V2=(*t_jetAK4P4)[AK4]+(*t_jetAK4P4)[AK3];
    if((W.M()<105. && W.M()>65.) || (V1.M()<105. && V1.M()>65.) || (V2.M()<105. && V2.M()>65.)){
      Wb2=(*t_jetAK4P4)[AK2]+(*t_jetAK4P4)[AK3]+(*t_jetAK4P4)[AK4];
      V=(*t_jetAK4P4)[AK5]+(*t_jetAK4P4)[AK6];
      chi[0]=pow((Wb1.M()-184)/(18.1),2.)+pow((Wb2.M()-171.3)/(17.7),2.)+pow((V.M()-88.5)/(9.5),2.);//+pow((v1.M()-88.5)/(9.5),2.0);//+pow((W.M()-80.8)/(9.4),2.);
    
      Wb1=V+(*t_jetAK4P4)[bid];
      chi_1=pow((Wb1.M()-171.3)/(17.7),2.)+pow((Wb2.M()-171.3)/(17.7),2.)+pow(((*t_jetWJet_MassPruned)[0]-88.5)/(9.5),2.);
      if(chi[0]<chi_1){
	chi[0]=chi_1;
	type[0]=1;
      }

      h_AK4WbTerm->Fill(pow((Wb2.M()-171.3)/(17.7),2.));
      h_AK4WTerm->Fill(pow((W.M()-80.8)/(9.4),2.));
      h_WbTerm->Fill(pow((Wb1.M()-184.)/(18.1),2.));
      h_ZTerm->Fill(pow((V.M()-88.5)/(9.5),2.));
      
      
    }
    
  }
  //========= case 1 ==================
  if((*t_jetAK4CSV)[AK2]>0.5426 || (*t_jetAK4CSV)[AK3]>0.5426 || (*t_jetAK4CSV)[AK5]>0.5426){
    Wb1= (*t_jetWJetP4)[0]+(*t_jetAK4P4)[bid];
    W=(*t_jetAK4P4)[AK2]+(*t_jetAK4P4)[AK3];
    V1=(*t_jetAK4P4)[AK2]+(*t_jetAK4P4)[AK5];
    V2=(*t_jetAK4P4)[AK5]+(*t_jetAK4P4)[AK3];
    if((W.M()<105. && W.M()>65.) || (V1.M()<105. && V1.M()>65.) || (V2.M()<105. && V2.M()>65.)){
      Wb2=W+(*t_jetAK4P4)[AK5];
      V=(*t_jetAK4P4)[AK4]+(*t_jetAK4P4)[AK6];
      chi[1]=pow((Wb1.M()-184.)/(18.1),2.)+pow((Wb2.M()-171.3)/(17.7),2.)+pow((V.M()-88.5)/(9.5),2.);//+pow((v1.M()-88.5)/(9.5),2.0);//+pow((W.M()-80.8)/(9.4),2.);
      
      Wb1=V+(*t_jetAK4P4)[bid];
      chi_1=pow((Wb1.M()-171.3)/(17.7),2.)+pow((Wb2.M()-171.3)/(17.7),2.)+pow(((*t_jetWJet_MassPruned)[0]-88.5)/(9.5),2.);
      if(chi[1]<chi_1){
	chi[1]=chi_1;
	type[1]=1;
      }

      h_AK4WbTerm->Fill(pow((Wb2.M()-171.3)/(17.7),2.));
      h_AK4WTerm->Fill(pow((W.M()-80.8)/(9.4),2.));
      h_WbTerm->Fill(pow((Wb1.M()-184.)/(18.1),2.));
      h_ZTerm->Fill(pow((V.M()-88.5)/(9.5),2.));
    }
  }
    //========= case 2 ==================
  if((*t_jetAK4CSV)[AK2]>0.5426 || (*t_jetAK4CSV)[AK3]>0.5426 || (*t_jetAK4CSV)[AK6]>0.5426){
    Wb1= (*t_jetWJetP4)[0]+(*t_jetAK4P4)[bid];
    W=(*t_jetAK4P4)[AK2]+(*t_jetAK4P4)[AK3];
    V1=(*t_jetAK4P4)[AK2]+(*t_jetAK4P4)[AK6];
    V2=(*t_jetAK4P4)[AK6]+(*t_jetAK4P4)[AK3];
    if((W.M()<105. && W.M()>65.) || (V1.M()<105. && V1.M()>65.) || (V2.M()<105. && V2.M()>65.)){
      Wb2=W+(*t_jetAK4P4)[AK6];
      V=(*t_jetAK4P4)[AK4]+(*t_jetAK4P4)[AK5];
      chi[2]=pow((Wb1.M()-184.)/(18.1),2.)+pow((Wb2.M()-171.3)/(17.7),2.)+pow((V.M()-88.5)/(9.5),2.);//+pow((v1.M()-88.5)/(9.5),2.0);//+pow((W.M()-80.8)/(9.4),2.);
      Wb1=V+(*t_jetAK4P4)[bid];
      chi_1=pow((Wb1.M()-171.3)/(17.7),2.)+pow((Wb2.M()-171.3)/(17.7),2.)+pow(((*t_jetWJet_MassPruned)[0]-88.5)/(9.5),2.);
      if(chi[2]<chi_1){
	chi[2]=chi_1;
	type[2]=1;
      }

      h_AK4WbTerm->Fill(pow((Wb2.M()-171.3)/(17.7),2.));
      h_AK4WTerm->Fill(pow((W.M()-80.8)/(9.4),2.));
      h_WbTerm->Fill(pow((Wb1.M()-184.)/(18.1),2.));
      h_ZTerm->Fill(pow((V.M()-88.5)/(9.5),2.));
    }
  }
    //========= case 3 ==================
  if((*t_jetAK4CSV)[AK2]>0.5426 || (*t_jetAK4CSV)[AK4]>0.5426 || (*t_jetAK4CSV)[AK5]>0.5426){
    Wb1= (*t_jetWJetP4)[0]+(*t_jetAK4P4)[bid];
    W=(*t_jetAK4P4)[AK2]+(*t_jetAK4P4)[AK4];
    V1=(*t_jetAK4P4)[AK2]+(*t_jetAK4P4)[AK5];
    V2=(*t_jetAK4P4)[AK4]+(*t_jetAK4P4)[AK5];
    if((W.M()<105. && W.M()>65.) || (V1.M()<105. && V1.M()>65.) || (V2.M()<105. && V2.M()>65.)){
      Wb2=W+(*t_jetAK4P4)[AK5];
      V=(*t_jetAK4P4)[AK3]+(*t_jetAK4P4)[AK6];
      chi[3]=pow((Wb1.M()-184.)/(18.1),2.)+pow((Wb2.M()-171.3)/(17.7),2.)+pow((V.M()-88.5)/(9.5),2.);//+pow((v1.M()-88.5)/(9.5),2.0);//+pow((W.M()-80.8)/(9.4),2.);
      Wb1=V+(*t_jetAK4P4)[bid];
      chi_1=pow((Wb1.M()-171.3)/(17.7),2.)+pow((Wb2.M()-171.3)/(17.7),2.)+pow(((*t_jetWJet_MassPruned)[0]-88.5)/(9.5),2.);
      if(chi[4]<chi_1){
	chi[4]=chi_1;
	type[4]=1;
      }

      h_AK4WbTerm->Fill(pow((Wb2.M()-171.3)/(17.7),2.));
      h_AK4WTerm->Fill(pow((W.M()-80.8)/(9.4),2.));
      h_WbTerm->Fill(pow((Wb1.M()-184.)/(18.1),2.));
      h_ZTerm->Fill(pow((V.M()-88.5)/(9.5),2.));
    }
  }
    //========= case 4 ==================
  if((*t_jetAK4CSV)[AK2]>0.5426 || (*t_jetAK4CSV)[AK4]>0.5426 || (*t_jetAK4CSV)[AK6]>0.5426){
    Wb1= (*t_jetWJetP4)[0]+(*t_jetAK4P4)[bid];
    W=(*t_jetAK4P4)[AK2]+(*t_jetAK4P4)[AK4];
    if(W.M()<105. && W.M()>65.){
      Wb2=W+(*t_jetAK4P4)[AK6];
      V=(*t_jetAK4P4)[AK3]+(*t_jetAK4P4)[AK5];
      chi[4]=pow((Wb1.M()-184.)/(18.1),2.)+pow((Wb2.M()-171.3)/(17.7),2.)+pow((V.M()-88.5)/(9.5),2.);//+pow((W.M()-80.8)/(9.4),2.);
      
      Wb1=V+(*t_jetAK4P4)[bid];
      chi_1=pow((Wb1.M()-171.3)/(17.7),2.)+pow((Wb2.M()-171.3)/(17.7),2.)+pow(((*t_jetWJet_MassPruned)[0]-88.5)/(9.5),2.);
      if(chi[4]<chi_1){
	chi[4]=chi_1;
	type[4]=1;
      }

      h_AK4WbTerm->Fill(pow((Wb2.M()-171.3)/(17.7),2.));
      h_AK4WTerm->Fill(pow((W.M()-80.8)/(9.4),2.));
      h_WbTerm->Fill(pow((Wb1.M()-184.)/(18.1),2.));
      h_ZTerm->Fill(pow((V.M()-88.5)/(9.5),2.));
    }
  }
    //========= case 5 ==================
  if((*t_jetAK4CSV)[AK2]>0.5426 || (*t_jetAK4CSV)[AK5]>0.5426 || (*t_jetAK4CSV)[AK6]>0.5426){
    Wb1= (*t_jetWJetP4)[0]+(*t_jetAK4P4)[bid];
    W=(*t_jetAK4P4)[AK2]+(*t_jetAK4P4)[AK5];
    V1=(*t_jetAK4P4)[AK2]+(*t_jetAK4P4)[AK4];
    V2=(*t_jetAK4P4)[AK4]+(*t_jetAK4P4)[AK3];
    if((W.M()<105. && W.M()>65.) || (V1.M()<105. && V1.M()>65.) || (V2.M()<105. && V2.M()>65.)){
      Wb2=W+(*t_jetAK4P4)[AK6];
      V=(*t_jetAK4P4)[AK3]+(*t_jetAK4P4)[AK4];
      chi[5]=pow((Wb1.M()-184.)/(18.1),2.)+pow((Wb2.M()-171.3)/(17.7),2.)+pow((V.M()-88.5)/(9.5),2.);//+pow((v1.M()-88.5)/(9.5),2.0);//+pow((W.M()-80.8)/(9.4),2.);
      Wb1=V+(*t_jetAK4P4)[bid];
      chi_1=pow((Wb1.M()-171.3)/(17.7),2.)+pow((Wb2.M()-171.3)/(17.7),2.)+pow(((*t_jetWJet_MassPruned)[0]-88.5)/(9.5),2.);
      if(chi[5]<chi_1){
	chi[5]=chi_1;
	type[5]=1;
      }

      h_AK4WbTerm->Fill(pow((Wb2.M()-171.3)/(17.7),2.));
      h_AK4WTerm->Fill(pow((W.M()-80.8)/(9.4),2.));
      h_WbTerm->Fill(pow((Wb1.M()-184.)/(18.1),2.));
      h_ZTerm->Fill(pow((V.M()-88.5)/(9.5),2.));
    }
  }
   //========= case 6 ==================
  if((*t_jetAK4CSV)[AK3]>0.5426 || (*t_jetAK4CSV)[AK4]>0.5426 || (*t_jetAK4CSV)[AK5]>0.5426){
    Wb1= (*t_jetWJetP4)[0]+(*t_jetAK4P4)[bid];
    W=(*t_jetAK4P4)[AK3]+(*t_jetAK4P4)[AK4];
    V1=(*t_jetAK4P4)[AK5]+(*t_jetAK4P4)[AK4];
    V2=(*t_jetAK4P4)[AK5]+(*t_jetAK4P4)[AK3];
    if((W.M()<105. && W.M()>65.) || (V1.M()<105. && V1.M()>65.) || (V2.M()<105. && V2.M()>65.)){
      Wb2=W+(*t_jetAK4P4)[AK5];
      V=(*t_jetAK4P4)[AK2]+(*t_jetAK4P4)[AK6];
      chi[6]=pow((Wb1.M()-184.)/(18.1),2.)+pow((Wb2.M()-171.3)/(17.7),2.)+pow((V.M()-88.5)/(9.5),2.);//+pow((v1.M()-88.5)/(9.5),2.0);//+pow((W.M()-80.8)/(9.4),2.);
      Wb1=V+(*t_jetAK4P4)[bid];
      chi_1=pow((Wb1.M()-171.3)/(17.7),2.)+pow((Wb2.M()-171.3)/(17.7),2.)+pow(((*t_jetWJet_MassPruned)[0]-88.5)/(9.5),2.);
      if(chi[6]<chi_1){
	chi[6]=chi_1;
	type[6]=1;
      }

      h_AK4WbTerm->Fill(pow((Wb2.M()-171.3)/(17.7),2.));
      h_AK4WTerm->Fill(pow((W.M()-80.8)/(9.4),2.));
      h_WbTerm->Fill(pow((Wb1.M()-184.)/(18.1),2.));
      h_ZTerm->Fill(pow((V.M()-88.5)/(9.5),2.));
    }
  }

   //========= case 7 ==================
  if((*t_jetAK4CSV)[AK3]>0.5426 || (*t_jetAK4CSV)[AK4]>0.5426 || (*t_jetAK4CSV)[AK6]>0.5426){
    Wb1= (*t_jetWJetP4)[0]+(*t_jetAK4P4)[bid];
    W=(*t_jetAK4P4)[AK3]+(*t_jetAK4P4)[AK4];
    V1=(*t_jetAK4P4)[AK6]+(*t_jetAK4P4)[AK4];
    V2=(*t_jetAK4P4)[AK6]+(*t_jetAK4P4)[AK3];
    if((W.M()<105. && W.M()>65.) || (V1.M()<105. && V1.M()>65.) || (V2.M()<105. && V2.M()>65.)){
      Wb2=W+(*t_jetAK4P4)[AK6];
      V=(*t_jetAK4P4)[AK2]+(*t_jetAK4P4)[AK5];
      chi[7]=pow((Wb1.M()-184.)/(18.1),2.)+pow((Wb2.M()-171.3)/(17.7),2.)+pow((V.M()-88.5)/(9.5),2.);//+pow((v1.M()-88.5)/(9.5),2.0);//+pow((W.M()-80.8)/(9.4),2.);
      Wb1=V+(*t_jetAK4P4)[bid];
      chi_1=pow((Wb1.M()-171.3)/(17.7),2.)+pow((Wb2.M()-171.3)/(17.7),2.)+pow(((*t_jetWJet_MassPruned)[0]-88.5)/(9.5),2.);
      if(chi[7]<chi_1){
	chi[7]=chi_1;
	type[7]=1;
      }

      h_AK4WbTerm->Fill(pow((Wb2.M()-171.3)/(17.7),2.));
      h_AK4WTerm->Fill(pow((W.M()-80.8)/(9.4),2.));
      h_WbTerm->Fill(pow((Wb1.M()-184.)/(18.1),2.));
      h_ZTerm->Fill(pow((V.M()-88.5)/(9.5),2.));
    }
  }

   //========= case 8 ==================
  if((*t_jetAK4CSV)[AK3]>0.5426 || (*t_jetAK4CSV)[AK5]>0.5426 || (*t_jetAK4CSV)[AK6]>0.5426){
    Wb1= (*t_jetWJetP4)[0]+(*t_jetAK4P4)[bid];
    W=(*t_jetAK4P4)[AK3]+(*t_jetAK4P4)[AK5];
    V1=(*t_jetAK4P4)[AK5]+(*t_jetAK4P4)[AK6];
    V2=(*t_jetAK4P4)[AK6]+(*t_jetAK4P4)[AK3];
    if((W.M()<105. && W.M()>65.) || (V1.M()<105. && V1.M()>65.) || (V2.M()<105. && V2.M()>65.)){
      Wb2=W+(*t_jetAK4P4)[AK6];
      V=(*t_jetAK4P4)[AK2]+(*t_jetAK4P4)[AK4];
      chi[8]=pow((Wb1.M()-184.)/(18.1),2.)+pow((Wb2.M()-171.3)/(17.7),2.)+pow((V.M()-88.5)/(9.5),2.);//+pow((v1.M()-88.5)/(9.5),2.0);//+pow((W.M()-80.8)/(9.4),2.);
      Wb1=V+(*t_jetAK4P4)[bid];
      chi_1=pow((Wb1.M()-171.3)/(17.7),2.)+pow((Wb2.M()-171.3)/(17.7),2.)+pow(((*t_jetWJet_MassPruned)[0]-88.5)/(9.5),2.);
      if(chi[8]<chi_1){
	chi[8]=chi_1;
	type[8]=1;
      }

      h_AK4WbTerm->Fill(pow((Wb2.M()-171.3)/(17.7),2.));
      h_AK4WTerm->Fill(pow((W.M()-80.8)/(9.4),2.));
      h_WbTerm->Fill(pow((Wb1.M()-184.)/(18.1),2.));
      h_ZTerm->Fill(pow((V.M()-88.5)/(9.5),2.));
    }
  }

   //========= case 9 ==================
  if((*t_jetAK4CSV)[AK4]>0.5426 || (*t_jetAK4CSV)[AK5]>0.5426 || (*t_jetAK4CSV)[AK6]>0.5426){
    Wb1= (*t_jetWJetP4)[0]+(*t_jetAK4P4)[bid];
    W=(*t_jetAK4P4)[AK4]+(*t_jetAK4P4)[AK5];
    V1=(*t_jetAK4P4)[AK6]+(*t_jetAK4P4)[AK4];
    V2=(*t_jetAK4P4)[AK5]+(*t_jetAK4P4)[AK6];
    if((W.M()<105. && W.M()>65.) || (V1.M()<105. && V1.M()>65.) || (V2.M()<105. && V2.M()>65.)){
      Wb2=W+(*t_jetAK4P4)[AK6];
      V=(*t_jetAK4P4)[AK2]+(*t_jetAK4P4)[AK3];
      chi[9]=pow((Wb1.M()-184.)/(18.1),2.)+pow((Wb2.M()-171.3)/(17.7),2.)+pow((V.M()-88.5)/(9.5),2.);//+pow((v1.M()-88.5)/(9.5),2.0);//+pow((W.M()-80.8)/(9.4),2.);
      Wb1=V+(*t_jetAK4P4)[bid];
      chi_1=pow((Wb1.M()-171.3)/(17.7),2.)+pow((Wb2.M()-171.3)/(17.7),2.)+pow(((*t_jetWJet_MassPruned)[0]-88.5)/(9.5),2.);
      if(chi[9]<chi_1){
	chi[9]=chi_1;
	type[9]=1;
      }

      h_AK4WbTerm->Fill(pow((Wb2.M()-171.3)/(17.7),2.));
      h_AK4WTerm->Fill(pow((W.M()-80.8)/(9.4),2.));
      h_WbTerm->Fill(pow((Wb1.M()-184.)/(18.1),2.));
      h_ZTerm->Fill(pow((V.M()-88.5)/(9.5),2.));
    }
  }
  
  double temp;
  for(int i=1;i<10;++i)
    {
        for(int j=0;j<(10-i);++j)
            if(chi[j]>chi[j+1])
            {
                temp=chi[j];
                chi[j]=chi[j+1];
                chi[j+1]=temp;
		
		temp=type[j];
		type[j]=type[j+1];
		type[j+1]=temp;

		chi_idx[j]=j+1;
		chi_idx[j+1]=j;
            }
    }

  res.push_back(chi[0]);
  res.push_back(chi_idx[0]);
  res.push_back(type[0]);
  return res;
}
vector <float> TPrimeAna::GetChi(int NrecoTop,int NrecoV,TLorentzVector v1,vector <int> goodMuid,vector <int> goodElid, vector <int> cleanJetAK4, vector <int> cleanBJetL){
  double dr=99;
  int match[]={99,99};
  int tcand[]={99,99};
	       
  int Vcheck=0; int Vidx=-99;
  vector <float> chi_output;

  //===================== to calculate the chi2==============================//
  double chi=10000;double chi_cf=10000;
  double chi_1=10000; double chi_2=10000;double chi_3=10000;
  int tidx=-99;int tidx1=-99;int Wbidx=-99;int Wbidx1=-99; int bidx1=-99; int bidx=-99;
  int ctr1=0;int ctr2=0;int ctr3=0;int ctrb=0;int ctr4=0;int ctr5=0;
  double diff1,diff2,diff3,diff4,temp;
  double chi_V2[12];
  int chi_idx[]={0,1,2,3,4,5,6,7,8,9,10,11};
  int type[]={0,0,0,0,0,0};
  for(int i=0;i<12;i++)chi_V2[i]=10000;
  TLorentzVector TPfl,TPfh,TPl1,TPh1,TPl2,TPh2, TPl3,TPh3,TPl4,TPh4,Wbf1,Wbf2,AKVf,V, V1,V2,Wb1,Wb2,Wb,Wbc;
  //cout<<"ok\n";
  if(NrecoV>=2 && NrecoTop>=1){//cout<<"1\n";
    for(unsigned int itT = 0 ; itT <t_jetTopJetP4->size(); ++itT) {
      for(unsigned int itW1 = 0 ; ; ++itW1) {
	for(unsigned int itW2 = itW1+1 ;itW2<t_jetWJetP4->size()  ; ++itW2) {
	  for (vector<int>::iterator it = cleanBJetL.begin() ; it != cleanBJetL.end(); ++it){
			 
	    Wb1=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*it];
	    Wb2=(*t_jetWJetP4)[itW2]+(*t_jetAK4P4)[*it];
			   
	    chi_1=pow(((*t_jetWJet_MassPruned)[itW1]-88.5)/(9.5),2.)+pow((Wb2.M()-184.)/(18.1),2.)+pow(((*t_jetTopJet_SoftDropMass)[itT]-179.7)/(18.4),2.);//+pow((v1.M()-88.5)/(9.5),2.0);//+pow(((*t_jetAK8_MassPruned)[*itW2]-80.84)/(9.4),2.);
	    //chi_1=pow(((*t_jetAK8_MassPruned)[*itW1]-88.3)/(10.1),2.)+pow((Wb2.M()-177.0)/(19.2),2.)+pow(((*t_jetTopJet_SoftDropMass)[itT]-177.0)/(22.),2.)+pow(((*t_jetAK8_MassPruned)[*itW2]-88.3)/(10.1),2.);
	    diff1=fabs((v1+(*t_jetTopJetP4)[itT]).M()-((*t_jetWJetP4)[itW1]+Wb2).M());
	    diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[itW1]+(*t_jetTopJetP4)[itT]).M());
	    if(diff1<diff2){TPl1=v1+(*t_jetTopJetP4)[itT];TPh1=(*t_jetWJetP4)[itW1]+Wb2;}
	    else {TPl1=v1+Wb2;TPh1=(*t_jetWJetP4)[itW1]+(*t_jetTopJetP4)[itT];}
			   
	    chi_2=pow(((*t_jetWJet_MassPruned)[itW2]-88.5)/(9.5),2.)+pow((Wb1.M()-184.)/(18.1),2.)+pow(((*t_jetTopJet_SoftDropMass)[itT]-179.7)/(18.4),2.);//+pow((v1.M()-88.5)/(9.5),2.0);//+pow(((*t_jetAK8_MassPruned)[*itW1]-80.84)/(9.4),2.);
	    //chi_2=pow(((*t_jetAK8_MassPruned)[*itW2]-88.3)/(10.1),2.)+pow((Wb1.M()-177.0)/(19.2),2.)+pow(((*t_jetTopJet_SoftDropMass)[itT]-177.0)/(22.),2.)+pow(((*t_jetAK8_MassPruned)[*itW1]-88.3)/(10.1),2.);
	    diff3=fabs((v1+(*t_jetTopJetP4)[itT]).M()-((*t_jetWJetP4)[itW1]+Wb1).M());
	    diff4=fabs((v1+Wb1).M()-((*t_jetWJetP4)[itW1]+(*t_jetTopJetP4)[itT]).M());
	    if(diff3<diff4){diff3;TPl2=v1+(*t_jetTopJetP4)[itT];TPh2=(*t_jetWJetP4)[itW1]+Wb2;}
	    else{ diff4;TPl2=v1+Wb2;TPh2=(*t_jetWJetP4)[itW1]+(*t_jetTopJetP4)[itT];}
			   
	    if(chi_1<chi_2 && chi_1<chi){
	      chi=chi_1;
	      TPfl=TPl1;TPfh=TPh1;
	    }
	    if(chi_2<chi_1 && chi_2<chi){
	      chi=chi_2;
	      TPfl=TPl2;TPfh=TPh2;
	    }
		       
	  }
	}
	ctr1++;
	if(ctr1>=t_jetWJetP4->size()-1){break;}
      }
    }
    if(chi>240.)chi=239.;
    chi_output.push_back(chi);
    chi_output.push_back(TPfl.M());
    chi_output.push_back(TPfh.M());
    return chi_output;
  }
  else if(!NrecoTop && NrecoV==1 && cleanJetAK4.size()>=6 ){
    //cout<<"2\n";
    for (vector<int>::iterator itAK1 = cleanJetAK4.begin() ; ; ++itAK1){
      for (vector<int>::iterator itAK2 = itAK1+1 ; ; ++itAK2){
	for (vector<int>::iterator itAK3 = itAK2+1 ; ; ++itAK3){
	  for (vector<int>::iterator itAK4 = itAK3+1 ; ; ++itAK4){
	    for (vector<int>::iterator itAK5 = itAK4+1 ; ; ++itAK5){
	      for (vector<int>::iterator itAK6 = itAK5+1 ;itAK6!= cleanJetAK4.end(); ++itAK6){

		vector <float> chi11,chi22,chi33,chi44,chi55,chi66;
		for(int i=1;i<=2;i++){
		  chi22.push_back(10000);
		  chi11.push_back(10000);
		  chi33.push_back(10000);
		  chi44.push_back(10000);
		  chi55.push_back(10000);
		  chi66.push_back(10000);
		}
			     
		if((*t_jetAK4CSV)[*itAK1]>0.5426 && ((*t_jetAK4CSV)[*itAK2]>0.5426 || (*t_jetAK4CSV)[*itAK3]>0.5426 || (*t_jetAK4CSV)[*itAK4]>0.5426 || (*t_jetAK4CSV)[*itAK5]>0.5426 || (*t_jetAK4CSV)[*itAK6]>0.5426)) {
			       
		  vector < float> chi1=chi_calculator(*itAK1,*itAK2,*itAK3,*itAK4,*itAK5,*itAK6, v1);
		  chi11=chi1;
		}
			     
		if((*t_jetAK4CSV)[*itAK2]>0.5426 && ((*t_jetAK4CSV)[*itAK1]>0.5426 || (*t_jetAK4CSV)[*itAK3]>0.5426 || (*t_jetAK4CSV)[*itAK4]>0.5426 || (*t_jetAK4CSV)[*itAK5]>0.5426 || (*t_jetAK4CSV)[*itAK6]>0.5426)) {
		  vector < float> chi2 =chi_calculator( *itAK2,*itAK1,*itAK3,*itAK4,*itAK5,*itAK6,v1);
		  chi22=chi2;
		}
			    
			    
			     
		if((*t_jetAK4CSV)[*itAK3]>0.5426 && ((*t_jetAK4CSV)[*itAK2]>0.5426 || (*t_jetAK4CSV)[*itAK1]>0.5426 || (*t_jetAK4CSV)[*itAK4]>0.5426 || (*t_jetAK4CSV)[*itAK5]>0.5426 || (*t_jetAK4CSV)[*itAK6]>0.5426)){
		  vector < float> chi3= chi_calculator(*itAK3,*itAK2,*itAK1,*itAK4,*itAK5,*itAK6,v1);
		  chi33=chi3;
		}
			     
		if((*t_jetAK4CSV)[*itAK4]>0.5426 && ((*t_jetAK4CSV)[*itAK2]>0.5426 || (*t_jetAK4CSV)[*itAK3]>0.5426 || (*t_jetAK4CSV)[*itAK1]>0.5426 || (*t_jetAK4CSV)[*itAK5]>0.5426 || (*t_jetAK4CSV)[*itAK6]>0.5426)){
		  vector < float>  chi4= chi_calculator(*itAK4,*itAK2,*itAK3,*itAK1,*itAK5,*itAK6,v1);
		  chi44=chi4;
		}
			     
		if((*t_jetAK4CSV)[*itAK5]>0.5426 && ((*t_jetAK4CSV)[*itAK2]>0.5426 || (*t_jetAK4CSV)[*itAK3]>0.5426 || (*t_jetAK4CSV)[*itAK4]>0.5426 || (*t_jetAK4CSV)[*itAK6]>0.5426 || (*t_jetAK4CSV)[*itAK6]>0.5426)){
		  vector < float> chi5= chi_calculator( *itAK5,*itAK2,*itAK3,*itAK4,*itAK1,*itAK6,v1);
		  chi55=chi5;
		}
			     
		if((*t_jetAK4CSV)[*itAK6]>0.5426 && ((*t_jetAK4CSV)[*itAK2]>0.5426 || (*t_jetAK4CSV)[*itAK3]>0.5426 || (*t_jetAK4CSV)[*itAK4]>0.5426 || (*t_jetAK4CSV)[*itAK5]>0.5426 || (*t_jetAK4CSV)[*itAK1]>0.5426)){
		  vector < float> chi6= chi_calculator( *itAK6,*itAK2,*itAK3,*itAK4,*itAK5,*itAK1,v1);
		  chi66=chi6;
		}

		if(chi11[0]<chi && chi11[0]<chi22[0] && chi11[0]<chi33[0] && chi11[0]<chi44[0] && chi11[0]<chi55[0] && chi11[0]<chi66[0]) {
		  chi=chi11[0];
		  Wb1= (*t_jetWJetP4)[0]+(*t_jetAK4P4)[*itAK1];
		  int temp =chi11[1];
		  switch(temp)
		    {
		    case 0 :
		      if(chi11[2]==0){
			Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4];
			V=(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK6];
			diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			break;
		      }
		      else{
			Wb1=(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK1];
			Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4];
			diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[0]+Wb2).M());
			diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[0]+Wb1).M());
			if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[0]+Wb2;}
			else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[0]+Wb1;}
			break;
		      }
		    case 1:
		      if(chi11[2]==0){
			Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK5];
			V=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK6];
			diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			break;
		      }
		      else{
			Wb1=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK1];
			Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK5];
			diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[0]+Wb2).M());
			diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[0]+Wb1).M());
			if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[0]+Wb2;}
			else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[0]+Wb1;}
			break;
		      }
		    case 2:
		      if(chi11[2]==0){
			Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK6];
			V=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5];
			diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			break;
		      }
		      else{
			Wb1=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK1];
			Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK6];
			diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[0]+Wb2).M());
			diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[0]+Wb1).M());
			if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[0]+Wb2;}
			else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[0]+Wb1;}
			break;
		      }
		    case 3 :
		      if(chi11[2]==0){
			Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5];
			V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK6];
			diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			break;
		      }
		      else{
			Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5];
			Wb1=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK1];
			diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[0]+Wb2).M());
			diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[0]+Wb1).M());
			if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[0]+Wb2;}
			else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[0]+Wb1;}
			break;
		      }
		    case 4 :
		      if(chi11[2]==0){
			Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK6];
			V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK5];
			diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			break;
		      }
		      else{
			Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK6];
			Wb1=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK1];
			diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[0]+Wb2).M());
			diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[0]+Wb1).M());
			if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[0]+Wb2;}
			else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[0]+Wb1;}
			break;
		      }
		    case 5 :
		      if(chi11[2]==0){
			Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK6];
			V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4];
			diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			break;
		      }
		      else{
			Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK6];
			Wb1=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK1];
			diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[0]+Wb2).M());
			diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[0]+Wb1).M());
			if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[0]+Wb2;}
			else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[0]+Wb1;}
			break;
		      }
		    case 6 :
		      if(chi11[2]==0){
			Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5];
			V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK6];
			diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			break;
		      }
		      else{
			Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK4];
			Wb1=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK1];
			diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[0]+Wb2).M());
			diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[0]+Wb1).M());
			if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[0]+Wb2;}
			else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[0]+Wb1;}
			break;
		      }
		    case 7 :
		      if(chi11[2]==0){
			Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK6];
			V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK5];
			diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			break;
		      }
		      else{
			Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK4];
			Wb1=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK1];
			diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[0]+Wb2).M());
			diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[0]+Wb1).M());
			if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[0]+Wb2;}
			else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[0]+Wb1;}
			break;
		      }
		    case 8 :
		      if(chi11[2]==0){
			Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK6];
			V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4];
			diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			break;
		      }
		      else{
			Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK5];
			Wb1=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK1];
			diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[0]+Wb2).M());
			diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[0]+Wb1).M());
			if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[0]+Wb2;}
			else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[0]+Wb1;}
			break;
		      }
		    case 9 :
		      if(chi11[2]==0){
			Wb2=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK6];
			V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3];
			diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			break;
		      }
		      else{
			Wb2=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK5];
			Wb1=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK1];
			diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[0]+Wb2).M());
			diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[0]+Wb1).M());
			if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[0]+Wb2;}
			else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[0]+Wb1;}

			break;
		      }
		    }
		}
		else{
		  if(chi22[0]<chi && chi22[0]<chi11[0] && chi22[0]<chi33[0] && chi22[0]<chi44[0] && chi22[0]<chi55[0] && chi22[0]<chi66[0]) {
		    chi=chi22[0];
		    Wb1= (*t_jetWJetP4)[0]+(*t_jetAK4P4)[*itAK2];
		    int temp= chi22[1];
		    switch(temp)
		      {
		      case 0 :
			if(chi22[2]==0){
			  Wb2=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4];
			  V=(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK6];
			  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			  break;
			}
			else{
			  Wb2=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4];
			  Wb1= (*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK2];
			  diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[0]+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[0]+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[0]+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[0]+Wb1;}
			  break;
			}
				     
		      case 1:
			if(chi22[2]==0){
			  Wb2=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK5];
			  V=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK6];
			  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			  break;
			}
			else{
			  Wb2=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK5];
			  Wb1= (*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK2];
			  diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[0]+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[0]+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[0]+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[0]+Wb1;}
			  break;
			}
		      case 2:
			if(chi22[2]==0){
			  Wb2=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK6];
			  V=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5];
			  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			  break;
			}
			else{
			  Wb2=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK6];
			  Wb1= (*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK2];
			  diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[0]+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[0]+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[0]+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[0]+Wb1;}
			  break;
			}
		      case 3 :
			if(chi22[2]==0){
			  Wb2=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5];
			  V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK6];
			  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			  break;
			}
			else{
			  Wb2=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5];
			  Wb1= (*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK2];
			  diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[0]+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[0]+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[0]+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[0]+Wb1;}
			  break;
			}
		      case 4 :
			if(chi22[2]==0){
			  Wb2=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK6];
			  V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK5];
			  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			  break;
			}
			else{
			  Wb2=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK6];
			  Wb1= (*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK2];
			  diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[0]+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[0]+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[0]+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[0]+Wb1;}
			  break;
			}
		      case 5 :
			if(chi22[2]==0){
			  Wb2=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK6];
			  V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4];
			  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			  break;
			}
			else{
			  Wb2=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK6];
			  Wb1= (*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK2];
			  diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[0]+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[0]+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[0]+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[0]+Wb1;}
			  break;
			}
		      case 6 :
			if(chi22[2]==0){
			  Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5];
			  V=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK6];
			  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			  break;
			}
			else{
			  Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK4];
			  Wb1= (*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK2];
			  diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[0]+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[0]+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[0]+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[0]+Wb1;}
			  break;
			}
		      case 7 :
			if(chi22[2]==0){
			  Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK6];
			  V=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK5];
			  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			  break;
			}
			else{
			  Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK4];
			  Wb1= (*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK2];
			  diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[0]+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[0]+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[0]+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[0]+Wb1;}
			  break;
			}
		      case 8 :
			if(chi22[2]==0){
			  Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK6];
			  V=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK4];
			  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			  break;
			}
			else{
			  Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK5];
			  Wb1= (*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK2];
			  diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[0]+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[0]+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[0]+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[0]+Wb1;}
			  break;
			}
		      case 9 :
			if(chi22[2]==0){
			  Wb2=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK6];
			  V=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK3];
			  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			  break;
			}
			else{
			  Wb2=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK5];
			  Wb1= (*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK2];
			  diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[0]+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[0]+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[0]+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[0]+Wb1;}
			  break;
			}
		      }
		  }
		  else{
		    if(chi33[0]<chi && chi33[0]<chi22[0] && chi33[0]<chi11[0] && chi33[0]<chi44[0] && chi33[0]<chi55[0] && chi33[0]<chi66[0]) {
		      chi=chi33[0];
		      Wb1= (*t_jetWJetP4)[0]+(*t_jetAK4P4)[*itAK3];
		      int temp=chi33[1];
		      switch(temp)
			{
			case 0 :
			  if(chi33[2]==0){
			    Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK4];
			    V=(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK6];
			    diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			    diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			    else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			    break;
			  }
			  else{
			    Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK4];
			    Wb1=(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK3];
			    diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			    diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			    else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			    break;
			  }
			case 1:
			  if(chi33[2]==0){
			    Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK5];
			    V=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK6];
			    diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			    diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			    else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			    break;
			  }
			  else{
			    Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK5];
			    Wb1=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK3];
			    diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			    diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			    else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			    break;
			  }
			case 2:
			  if(chi33[2]==0){
			    Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK6];
			    V=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5];
			    diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			    diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			    else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			    break;
			  }
			  else{
			    Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK6];
			    Wb1=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK3];
			    diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			    diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			    else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			    break;
			  }
			case 3 :
			  if(chi33[2]==0){
			    Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5];
			    V=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK6];
			    diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			    diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			    else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			    break;
			  }
			  else{
			    Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5];
			    Wb1=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK3];
			    diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			    diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			    else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			    break;
			  }
			case 4 :
			  if(chi33[2]==0){
			    Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK6];
			    V=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK5];
			    diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			    diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			    else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			    break;
			  }
			  else{
			    Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK6];
			    Wb1=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK3];
			    diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			    diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			    else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			    break;
			  }
			case 5 :
			  if(chi33[2]==0){
			    Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK6];
			    V=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK4];
			    diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			    diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			    else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			    break;
			  }
			  else{
			    Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK6];
			    Wb1=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK3];
			    diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			    diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			    else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			    break;
			  }
			case 6 :
			  if(chi33[2]==0){
			    Wb2=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5];
			    V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK6];
			    diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			    diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			    else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			    break;
			  }
			  else{
			    Wb2=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK4];
			    Wb1=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK3];
			    diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			    diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			    else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			    break;
			  }
			case 7 :
			  if(chi33[2]==0){
			    Wb2=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK6];
			    V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK5];
			    diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			    diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			    else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			    break;
			  }
			  else{
			    Wb2=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK4];
			    Wb1=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK3];
			    diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			    diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			    else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			    break;
			  }
			case 8 :
			  if(chi33[2]==0){
			    Wb2=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK6];
			    V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4];
			    diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			    diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			    else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			    break;
			  }
			  else{
			    Wb2=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK5];
			    Wb1=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK3];
			    diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			    diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			    else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			    break;
			  }
			case 9 :
			  if(chi33[2]==0){
			    Wb2=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK6];
			    V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK1];
			    diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			    diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			    else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			    break;
			  }
			  else{
			    Wb2=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK5];
			    Wb1=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK3];
			    diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			    diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			    else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			    break;
			  }
			}
		    }
		    else{
		      if(chi44[0]<chi && chi44[0]<chi22[0] && chi44[0]<chi33[0] && chi44[0]<chi11[0] && chi44[0]<chi55[0] && chi44[0]<chi66[0]) {
			chi=chi44[0];
			int temp= chi44[1];
			Wb1= (*t_jetWJetP4)[0]+(*t_jetAK4P4)[*itAK4];
			switch(temp)
			  {
			  case 0 :
			    if(chi44[2]==0){
			      Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK1];
			      V=(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK6];
			      diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			      diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			      if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			      else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			      break;
			    }
			    else{
			      Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK1];
			      Wb1=(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK4];
			      diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			      diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			      if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			      else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			      break;
			    }
			  case 1:
			    if(chi44[2]==0){
			      Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK5];
			      V=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK6];
			      diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			      diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			      if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			      else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			      break;
			    }
			    else{
			      Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK5];
			      Wb1=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK4];
			      diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			      diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			      if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			      else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			      break;
			    }
			  case 2:
			    if(chi44[2]==0){
			      Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK6];
			      V=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK5];
			      diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			      diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			      if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			      else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			      break;
			    }
			    else{
			      Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK6];
			      Wb1=(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK4];
			      diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			      diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			      if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			      else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			      break;
			    }
			  case 3 :
			    if(chi44[2]==0){
			      Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK5];
			      V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK6];
			      diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			      diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			      if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			      else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			      break;
			    }
			    else{
			      Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK5];
			      Wb1=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK4];
			      diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			      diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			      if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			      else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			      break;
			    }
			  case 4 :
			    if(chi44[2]==0){
			      Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK6];
			      V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK5];
			      diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			      diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			      if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			      else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			      break;
			    }
			    else{
			      Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK6];
			      Wb1=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK4];
			      diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			      diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			      if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			      else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			      break;
			    }
			  case 5 :
			    if(chi44[2]==0){
			      Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK6];
			      V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK1];
			      diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			      diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			      if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			      else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			      break;
			    }
			    else{
			      Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK6];
			      Wb1=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK4];
			      diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			      diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			      if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			      else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			      break;
			    }
			  case 6 :
			    if(chi44[2]==0){
			      Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK5];
			      V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK6];
			      diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			      diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			      if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			      else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			      break;
			    }
			    else{
			      Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK1];
			      Wb1=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK4];
			      diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			      diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			      if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			      else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			      break;
			    }
			  case 7 :
			    if(chi44[2]==0){
			      Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK6];
			      V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK5];
			      diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			      diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			      if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			      else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			      break;
			    }
			    else{
			      Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK1];
			      Wb1=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK4];
			      diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			      diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			      if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			      else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			      break;
			    }
			  case 8 :
			    if(chi44[2]==0){
			      Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK6];
			      V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK1];
			      diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			      diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			      if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			      else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			      break;
			    }
			    else{
			      Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK5];
			      Wb1=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK4];
			      diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			      diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			      if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			      else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			      break;
			    }
			  case 9 :
			    if(chi44[2]==0){
			      Wb2=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK6];
			      V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3];
			      diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			      diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			      if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			      else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			      break;
			    }
			    else{
			      Wb2=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK5];
			      Wb1=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4];
			      diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			      diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			      if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			      else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			      break;
			    }
			  }
		      }
		      else{
			if(chi55[0]<chi && chi55[0]<chi22[0] && chi55[0]<chi33[0] && chi55[0]<chi44[0] && chi55[0]<chi11[0] && chi55[0]<chi66[0]) {
			  chi=chi55[0];
			  Wb1= (*t_jetWJetP4)[0]+(*t_jetAK4P4)[*itAK5];
			  int temp=chi55[1];
			  switch(temp)
			    {
			    case 0 :
			      if(chi55[2]==0){
				Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4];
				V=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK6];
				diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				break;
			      }
			      else{
				Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4];
				Wb1=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK5];
				diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				break;
			      }
			    case 1:
			      if(chi55[2]==0){
				Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK1];
				V=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK6];
				diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				break;
			      }
			      else{
				Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK1];
				Wb1=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK5];
				diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				break;
			      }
			    case 2:
			      if(chi55[2]==0){
				Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK6];
				V=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK1];
				diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				break;
			      }
			      else{
				Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK6];
				Wb1=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK5];
				diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				break;
			      }
			    case 3 :
			      if(chi55[2]==0){
				Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK1];
				V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK6];
				diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				break;
			      }
			      else{
				Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK1];
				Wb1=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK5];
				diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				break;
			      } 
			    case 4 :
			      if(chi55[2]==0){
				Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK6];
				V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK1];
				diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				break;
			      }
			      else{
				Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK6];
				Wb1=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK5];
				diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				break;
			      }
			    case 5 :
			      if(chi55[2]==0){
				Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK6];
				V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4];
				diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				break;
			      }
			      else{
				Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK6];
				Wb1=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5];
				diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				break;
			      }
			    case 6 :
			      if(chi55[2]==0){
				Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK1];
				V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK6];
				diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				break;
			      }
			      else{
				Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK4];
				Wb1=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK5];
				diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				break;
			      }
			    case 7 :
			      if(chi55[2]==0){
				Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK6];
				V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK1];
				diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				break;
			      }
			      else{
				Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK4];
				Wb1=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK5];
				diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				break;
			      }
			    case 8 :
			      if(chi55[2]==0){
				Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK6];
				V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4];
				diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				break;
			      }
			      else{
				Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK1];
				Wb1=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5];
				diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				break;
			      }
			    case 9 :
			      if(chi55[2]==0){
				Wb2=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK6];
				V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3];
				diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				break;
			      }
			      else{
				Wb2=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK6]+(*t_jetAK4P4)[*itAK1];
				Wb1=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK5];
				diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				break;
			      }
			    }
			}
			else{
			  if(chi66[0]<chi && chi66[0]<chi22[0] && chi66[0]<chi33[0] && chi66[0]<chi44[0] && chi66[0]<chi55[0] && chi66[0]<chi11[0]) {
			    chi=chi66[0];
			    Wb1= (*t_jetWJetP4)[0]+(*t_jetAK4P4)[*itAK6];
			    int temp=chi66[1];
			    switch(temp)
			      {
			      case 0 :
				if(chi66[2]==0){
				  Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4];
				  V=(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK1];
				  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				  break;
				}
				else{
				  Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4];
				  Wb1=(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK6];
				  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				  break;
				}
			      case 1:
				if(chi66[2]==0){
				  Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK5];
				  V=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK1];
				  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				  break;
				}
				else{
				  Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK5];
				  Wb1=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK6];
				  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				  break;
				}
			      case 2:
				if(chi66[2]==0){
				  Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK1];
				  V=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5];
				  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				  break;
				}
				else{
				  Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK1];
				  Wb1=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK6];
				  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				  break;
				}
			      case 3 :
				if(chi66[2]==0){
				  Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5];
				  V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK1];
				  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				  break;
				}
				else{
				  Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5];
				  Wb1=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK6];
				  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				  break;
				}
			      case 4 :
				if(chi66[2]==0){
				  Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK1];
				  V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK5];
				  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				  break;
				}
				else{
				  Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK1];
				  Wb1=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK6];
				  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				  break;
				}
			      case 5 :
				if(chi66[2]==0){
				  Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK1];
				  V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4];
				  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				  break;
				}
				else{
				  Wb2=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK1];
				  Wb1=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK6];
				  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				  break;
				}
			      case 6 :
				if(chi66[2]==0){
				  Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5];
				  V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK1];
				  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				  break;
				}
				else{
				  Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK4];
				  Wb1=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK6];
				  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				  break;
				}
			      case 7 :
				if(chi66[2]==0){
				  Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK1];
				  V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK5];
				  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				  break;
				}
				else{
				  Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK4];
				  Wb1=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK6];
				  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				  break;
				}
			      case 8 :
				if(chi66[2]==0){
				  Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK1];
				  V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4];
				  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				  break;
				}
				else{
				  Wb2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK5];
				  Wb1=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK6];
				  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				  break; 
				}
			      case 9 :
				if(chi66[2]==0){
				  Wb2=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK1];
				  V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3];
				  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				  break;
				}
				else{
				  Wb2=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK5];
				  Wb1=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK6];
				  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
				  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
				  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
				  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
				  break;
				}
			      }
			  }
			}
		      }
		    }
		  }
		}
	      }
	      ctr1++;
	      if(ctr1>cleanJetAK4.size()-6){break;}
	    }
	    ctr2++;
	    if(ctr2>cleanJetAK4.size()-6){break;}
	  }
	  ctr3++;
	  if(ctr3>cleanJetAK4.size()-6){break;} 
	}
	ctr4++;
	if(ctr4>cleanJetAK4.size()-6){break;}
      }
      ctr5++;
      if(ctr5>cleanJetAK4.size()-6){break;}
    }
    if(chi>240. && chi!=10000)chi=239.;

    //if(TPfl.M()<50.0 || TPfh.M()<50.)cout<<chi<<" meow "<<endl;
    chi_output.push_back(chi);
    chi_output.push_back(TPfl.M());
    chi_output.push_back(TPfh.M());
    return chi_output;
  }

  else if(!NrecoTop && NrecoV>=2 && cleanJetAK4.size()>=4 ){
    //cout<<"3\n";
    for(unsigned int itW=0 ; ; ++itW) {
      for (unsigned int itW1 = itW+1 ;itW1<t_jetWJetP4->size(); ++itW1){
	for (vector<int>::iterator itAK1 = cleanJetAK4.begin() ; ; ++itAK1){
	  for (vector<int>::iterator itAK2 = itAK1+1 ; ; ++itAK2){
	    for (vector<int>::iterator itAK3 = itAK2+1 ; ; ++itAK3){
	      for (vector<int>::iterator itAK4 = itAK3+1 ; itAK4 != cleanJetAK4.end(); ++itAK4){

		if((*t_jetAK4CSV)[*itAK1]>0.5426 && (*t_jetAK4CSV)[*itAK2]>0.5426){
		  Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK1];
		  Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK2];
		  V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4];
			       
		  chi_V2[0]=pow((V.M()-88.5)/(9.5),2.)+pow((Wb2.M()-184)/(18.1),2.)+pow((Wb1.M()-184)/(18.1),2.);//+pow(((*t_jetWJetP4)[itW].M()-80.84)/(9.4),2.)+pow(((*t_jetWJetP4)[itW1].M()-80.84)/(9.4),2.);
		  if(V.M()>65 && V.M()<105){
		    Wb1=V+(*t_jetAK4P4)[*itAK1];
		    chi_1=pow(((*t_jetWJetP4)[itW].M()-88.5)/(9.5),2.)+pow((Wb2.M()-184)/(18.1),2.)+pow((Wb1.M()-171.3)/(17.7),2.);
		    if(chi_1<chi_V2[0]){chi_V2[0]=chi_1; type[0]=1;}
		    Wb2=V+(*t_jetAK4P4)[*itAK2];
		    Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK1];
		    chi_2=pow(((*t_jetWJetP4)[itW1].M()-88.5)/(9.5),2.)+pow((Wb1.M()-184)/(18.1),2.)+pow((Wb2.M()-171.3)/(17.7),2.);
		    if(chi_2<chi_V2[0]){chi_V2[0]=chi_2; type[0]=2;}
		  }
		  Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK2];
		  Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK1];
		  V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4];
		  chi_V2[1]=pow((V.M()-88.5)/(9.5),2.)+pow((Wb2.M()-184)/(18.1),2.)+pow((Wb1.M()-184)/(18.1),2.);//+pow((v1.M()-88.5)/(9.5),2.0);//+pow(((*t_jetWJetP4)[itW].M()-80.84)/(9.4),2.)+pow(((*t_jetWJetP4)[itW1].M()-80.84)/(9.4),2.); 
		}

		if((*t_jetAK4CSV)[*itAK1]>0.5426 && (*t_jetAK4CSV)[*itAK3]>0.5426){
		  Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK1];
		  Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK3];
		  V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4];
		  chi_V2[2]=pow((V.M()-88.5)/(9.5),2.)+pow((Wb2.M()-184)/(18.1),2.)+pow((Wb1.M()-184)/(18.1),2.);//+pow((v1.M()-88.5)/(9.5),2.0);//+pow(((*t_jetWJetP4)[itW].M()-80.84)/(9.4),2.);//+pow(((*t_jetWJetP4)[itW1].M()-80.84)/(9.4),2.);
			       
		  if(V.M()>65 && V.M()<105){
		    Wb1=V+(*t_jetAK4P4)[*itAK1];
		    chi_1=pow(((*t_jetWJetP4)[itW].M()-88.5)/(9.5),2.)+pow((Wb2.M()-184)/(18.1),2.)+pow((Wb1.M()-171.3)/(17.7),2.);
		    if(chi_1<chi_V2[2]){chi_V2[2]=chi_1; type[1]=1;}
			     
		    chi_2=pow(((*t_jetWJetP4)[itW1].M()-88.5)/(9.5),2.)+pow((Wb1.M()-184)/(18.1),2.)+pow((Wb2.M()-171.3)/(17.7),2.);
		    if(chi_2<chi_V2[2]){chi_V2[2]=chi_2; type[1]=2;}
		  }

		  Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK3];
		  Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK1];
		  V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4];
		  chi_V2[3]=pow((V.M()-88.5)/(9.5),2.)+pow((Wb2.M()-184)/(18.1),2.)+pow((Wb1.M()-184)/(18.1),2.);//+pow((v1.M()-88.5)/(9.5),2.0);//+pow(((*t_jetWJetP4)[itW].M()-80.84)/(9.4),2.);//+pow(((*t_jetWJetP4)[itW1].M()-80.84)/(9.4),2.);
		}

		if((*t_jetAK4CSV)[*itAK1]>0.5426 && (*t_jetAK4CSV)[*itAK4]>0.5426){
		  Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK1];
		  Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK4];
		  V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3];
		  chi_V2[4]=pow((V.M()-88.5)/(9.5),2.)+pow((Wb2.M()-184)/(18.1),2.)+pow((Wb1.M()-184)/(18.1),2.);//+pow((v1.M()-88.5)/(9.5),2.0);//+pow(((*t_jetWJetP4)[itW].M()-80.84)/(9.4),2.);//+pow(((*t_jetWJetP4)[itW1].M()-80.84)/(9.4),2.);

		  if(V.M()>65 && V.M()<105){
				 
		    Wb1=V+(*t_jetAK4P4)[*itAK1];
		    chi_1=pow(((*t_jetWJetP4)[itW].M()-88.5)/(9.5),2.)+pow((Wb2.M()-184)/(18.1),2.)+pow((Wb1.M()-171.3)/(17.7),2.);
		    if(chi_1<chi_V2[4]){chi_V2[4]=chi_1; type[2]=1;}
			     
		    Wb2=V+(*t_jetAK4P4)[*itAK4];
		    Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK1];
		    chi_2=pow(((*t_jetWJetP4)[itW1].M()-88.5)/(9.5),2.)+pow((Wb1.M()-184)/(18.1),2.)+pow((Wb2.M()-171.3)/(17.7),2.);
		    if(chi_2<chi_V2[4]){chi_V2[4]=chi_2; type[2]=2;}
		  }

		  Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK4];
		  Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK1];
		  V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK2];
		  chi_V2[5]=pow((V.M()-88.5)/(9.5),2.)+pow((Wb2.M()-184)/(18.1),2.)+pow((Wb1.M()-184)/(18.1),2.);//+pow((v1.M()-88.5)/(9.5),2.0);//+pow(((*t_jetWJetP4)[itW].M()-80.84)/(9.4),2.);//+pow(((*t_jetWJetP4)[itW1].M()-80.84)/(9.4),2.);
		}
		if((*t_jetAK4CSV)[*itAK2]>0.5426 && (*t_jetAK4CSV)[*itAK3]>0.5426){
		  Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK2];
		  Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK3];
		  V=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK4];
		  chi_V2[6]=pow((V.M()-88.5)/(9.5),2.)+pow((Wb2.M()-184)/(18.1),2.)+pow((Wb1.M()-184)/(18.1),2.);//+pow((v1.M()-88.5)/(9.5),2.0);//+pow(((*t_jetWJetP4)[itW].M()-80.84)/(9.4),2.);//+pow(((*t_jetWJetP4)[itW1].M()-80.84)/(9.4),2.);
		  if(V.M()>65 && V.M()<105){
		    Wb1=V+(*t_jetAK4P4)[*itAK2];
		    chi_1=pow(((*t_jetWJetP4)[itW].M()-88.5)/(9.5),2.)+pow((Wb2.M()-184)/(18.1),2.)+pow((Wb1.M()-171.3)/(17.7),2.);
		    if(chi_1<chi_V2[6]){chi_V2[6]=chi_1; type[3]=1;}
			     
		    Wb2=V+(*t_jetAK4P4)[*itAK3];
		    Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK2];
		    chi_2=pow(((*t_jetWJetP4)[itW1].M()-88.5)/(9.5),2.)+pow((Wb1.M()-184)/(18.1),2.)+pow((Wb2.M()-171.3)/(17.7),2.);
		    if(chi_2<chi_V2[6]){chi_V2[6]=chi_2; type[3]=2;}
		  }

		  Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK3];
		  Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK2];
		  V=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK4];
		  chi_V2[7]=pow((V.M()-88.5)/(9.5),2.)+pow((Wb2.M()-184)/(18.1),2.)+pow((Wb1.M()-184)/(18.1),2.);//+pow((v1.M()-88.5)/(9.5),2.0);//+pow(((*t_jetWJetP4)[itW].M()-80.84)/(9.4),2.);//+pow(((*t_jetWJetP4)[itW1].M()-80.84)/(9.4),2.);
		}

		if((*t_jetAK4CSV)[*itAK2]>0.5426 && (*t_jetAK4CSV)[*itAK4]>0.5426){
		  Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK2];
		  Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK4];
		  V=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK3];
		  chi_V2[8]=pow((V.M()-88.5)/(9.5),2.)+pow((Wb2.M()-184)/(18.1),2.)+pow((Wb1.M()-184)/(18.1),2.);//+pow((v1.M()-88.5)/(9.5),2.0);//+pow(((*t_jetWJetP4)[itW].M()-80.84)/(9.4),2.);//+pow(((*t_jetWJetP4)[itW1].M()-80.84)/(9.4),2.);
		  if(V.M()>65 && V.M()<105){
		    Wb1=V+(*t_jetAK4P4)[*itAK2];
				 
		    chi_1=pow(((*t_jetWJetP4)[itW].M()-88.5)/(9.5),2.)+pow((Wb2.M()-184)/(18.1),2.)+pow((Wb1.M()-171.3)/(17.7),2.);
				 
		    if(chi_1<chi_V2[8]){chi_V2[8]=chi_1; type[4]=1;}
			     
		    Wb2=V+(*t_jetAK4P4)[*itAK4];
		    Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK2];
		    chi_2=pow(((*t_jetWJetP4)[itW1].M()-88.5)/(9.5),2.)+pow((Wb1.M()-184)/(18.1),2.)+pow((Wb2.M()-171.3)/(17.7),2.);
		    if(chi_2<chi_V2[8]){chi_V2[8]=chi_2; type[4]=2;}
		  }


		  Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK4];
		  Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK2];
		  V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK1];
		  chi_V2[9]=pow((V.M()-88.5)/(9.5),2.)+pow((Wb2.M()-184)/(18.1),2.)+pow((Wb1.M()-184)/(18.1),2.);//+pow((v1.M()-88.5)/(9.5),2.0);//+pow(((*t_jetWJetP4)[itW].M()-80.84)/(9.4),2.);//+pow(((*t_jetWJetP4)[itW1].M()-80.84)/(9.4),2.);
		}

		if((*t_jetAK4CSV)[*itAK3]>0.5426 && (*t_jetAK4CSV)[*itAK4]>0.5426){
		  Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK3];
		  Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK4];
		  V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK1];
		  chi_V2[10]=pow((V.M()-88.5)/(9.5),2.)+pow((Wb2.M()-184)/(18.1),2.)+pow((Wb1.M()-184)/(18.1),2.);//+pow((v1.M()-88.5)/(9.5),2.0);//+pow(((*t_jetWJetP4)[itW].M()-80.84)/(9.4),2.);//+pow(((*t_jetWJetP4)[itW1].M()-80.84)/(9.4),2.);
			       
		  if(V.M()>65 && V.M()<105){
		    Wb1=V+(*t_jetAK4P4)[*itAK3];

		    chi_1=pow(((*t_jetWJetP4)[itW].M()-88.5)/(9.5),2.)+pow((Wb2.M()-184)/(18.1),2.)+pow((Wb1.M()-171.3)/(17.7),2.);
		    if(chi_1<chi_V2[10]){chi_V2[10]=chi_1; type[5]=1;}
			     
		    Wb2=V+(*t_jetAK4P4)[*itAK4];
		    Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK3];
		    chi_2=pow(((*t_jetWJetP4)[itW1].M()-88.5)/(9.5),2.)+pow((Wb1.M()-184)/(18.1),2.)+pow((Wb2.M()-171.3)/(17.7),2.);
		    if(chi_2<chi_V2[10]){chi_V2[10]=chi_2; type[10]=2;}
		  }

			       
		  Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK4];
		  Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK3];
		  V=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK2];
		  chi_V2[11]=pow((V.M()-88.5)/(9.5),2.)+pow((Wb2.M()-184)/(18.1),2.)+pow((Wb1.M()-184)/(18.1),2.);//+pow((v1.M()-88.5)/(9.5),2.0);//+pow(((*t_jetWJetP4)[itW].M()-80.84)/(9.4),2.);//+pow(((*t_jetWJetP4)[itW1].M()-80.84)/(9.4),2.);
			       
		}
		for(int i=0;i<12;i++)chi_idx[i]=i;
		for(int i=1;i<12;++i)
		  {
		    for(int j=0;j<(12-i);++j)
		      if(chi_V2[j]>chi_V2[j+1])
			{
			  temp=chi_V2[j];
			  chi_V2[j]=chi_V2[j+1];
			  chi_V2[j+1]=temp;
				       
			  temp=chi_idx[j];
			  chi_idx[j]=chi_idx[j+1];
			  chi_idx[j+1]=temp;
			}
		  }
			     
		if(chi_V2[0]<chi){
		  chi=chi_V2[0];
		  switch(chi_idx[0])
		    {
		    case 0 :
		      if(type[0]==1){
			//Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK1];
			Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK2];
			V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4];
			Wb1=V+(*t_jetAK4P4)[*itAK1];			     
			diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[itW]+Wb2).M());
			diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[itW]+Wb1).M());
			if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[itW]+Wb2;}
			else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[itW]+Wb1;}
			break;
		      }
		      else{ 
			if(type[0]==2){
			  V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4];
			  Wb2=V+(*t_jetAK4P4)[*itAK2];
			  Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK1];
			  diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[itW1]+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[itW1]+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[itW1]+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[itW1]+Wb1;}
			  break;
			}
			else{
			  Wb1=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK2];
			  Wb2=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK1];
			  V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4];
			  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			  break;
			}
		      }
		    case 1:
		      Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK2];
		      Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK1];
		      V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4];
		      diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
		      diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
		      if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
		      else {TPfl=v1+Wb2;TPfh=V+Wb1;}
		      break;
		    case 2:
		      if(type[1]==1){
			//Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK1];                                                          
			Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK3];
			V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4];
			Wb1=V+(*t_jetAK4P4)[*itAK1];
			diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[itW]+Wb2).M());
			diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[itW]+Wb1).M());
			if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[itW]+Wb2;}
			else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[itW]+Wb1;}
			break;
		      }
		      else{ 
			if(type[1]==2){
			  V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4];
			  Wb2=V+(*t_jetAK4P4)[*itAK3];
			  Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK1];
			  diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[itW1]+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[itW1]+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[itW1]+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[itW1]+Wb1;}
			  break;
			}
			else{
			  Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK1];
			  Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK3];
			  V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4];
			  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			  break;
			}
		      }
		    case 3 :
		      Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK3];
		      Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK1];
		      V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4];
		      diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
		      diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
		      if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
		      else {TPfl=v1+Wb2;TPfh=V+Wb1;}
		      break;
		    case 4 :
		      if(type[2]==1){
			//Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK1];                                                          
			Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK4];
			V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3];
			Wb1=V+(*t_jetAK4P4)[*itAK1];
			diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[itW]+Wb2).M());
			diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[itW]+Wb1).M());
			if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[itW]+Wb2;}
			else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[itW]+Wb1;}
			break;
		      }
		      else{
			if(type[2]==2){
			  V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3];
			  Wb2=V+(*t_jetAK4P4)[*itAK4];
			  Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK1];
			  diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[itW1]+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[itW1]+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[itW1]+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[itW1]+Wb1;}
			  break;
			}
			else{
			  Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK1];
			  Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK4];
			  V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3];
			  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			  break;
			}
		      }
		    case 5 :
		      Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK4];
		      Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK1];
		      V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK2];
		      diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
		      diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
		      if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
		      else {TPfl=v1+Wb2;TPfh=V+Wb1;}
		      break;
		    case 6 :
		      if(type[3]==1){
			//Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK1];                                                          
			Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK3];
			V=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK4];
			Wb1=V+(*t_jetAK4P4)[*itAK2];
			diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[itW]+Wb2).M());
			diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[itW]+Wb1).M());
			if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[itW]+Wb2;}
			else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[itW]+Wb1;}
			break;
		      }
		      else{
			if(type[3]==2){
			  V=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK4];
			  Wb2=V+(*t_jetAK4P4)[*itAK3];
			  Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK2];
			  diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[itW1]+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[itW1]+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[itW1]+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[itW1]+Wb1;}
			  break;
			}
			else{
			  Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK2];
			  Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK3];
			  V=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK4];
			  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			  break;
			}
		      }
		    case 7 :
		      Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK3];
		      Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK2];
		      V=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK4];
		      diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
		      diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
		      if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
		      else {TPfl=v1+Wb2;TPfh=V+Wb1;}
		      break;
		    case 8 :
		      if(type[4]==1){
			//Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK1];                                                          
			Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK4];
			V=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK3];
			Wb1=V+(*t_jetAK4P4)[*itAK2];
			diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[itW]+Wb2).M());
			diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[itW]+Wb1).M());
			if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[itW]+Wb2;}
			else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[itW]+Wb1;}
			break;
		      }
		      else{
			if(type[4]==2){
			  V=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK3];
			  Wb2=V+(*t_jetAK4P4)[*itAK4];
			  Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK2];
			  diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[itW1]+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[itW1]+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[itW1]+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[itW1]+Wb1;}
			  break;
			}
			else{
			  Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK2];
			  Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK4];
			  V=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK3];
			  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			  break;
			}
		      }
		    case 9 :
		      Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK4];
		      Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK2];
		      V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK1];
		      diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
		      diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
		      if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
		      else {TPfl=v1+Wb2;TPfh=V+Wb1;}
		      break;
		    case 10 :
		      if(type[5]==1){
			//Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK1];                                                          
			Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK4];
			V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK1];
			Wb1=V+(*t_jetAK4P4)[*itAK3];
			diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[itW]+Wb2).M());
			diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[itW]+Wb1).M());
			if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[itW]+Wb2;}
			else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[itW]+Wb1;}
			break;
		      }
		      else{
			if(type[5]==2){
			  V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK1];
			  Wb2=V+(*t_jetAK4P4)[*itAK4];
			  Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK3];
			  diff1=fabs((v1+Wb1).M()-((*t_jetWJetP4)[itW1]+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-((*t_jetWJetP4)[itW1]+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=(*t_jetWJetP4)[itW1]+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=(*t_jetWJetP4)[itW1]+Wb1;}
			  break;
			}
			else{
			  Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK3];
			  Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK4];
			  V=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK1];
			  diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
			  diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
			  if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
			  else {TPfl=v1+Wb2;TPfh=V+Wb1;}
			  break;
			}
		      }
		    case 11 :
		      Wb1=(*t_jetWJetP4)[itW]+(*t_jetAK4P4)[*itAK4];
		      Wb2=(*t_jetWJetP4)[itW1]+(*t_jetAK4P4)[*itAK3];
		      V=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK2];
		      diff1=fabs((v1+Wb1).M()-(V+Wb2).M());
		      diff2=fabs((v1+Wb2).M()-(V+Wb1).M());
		      if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+Wb2;}
		      else {TPfl=v1+Wb2;TPfh=V+Wb1;}
		      break;
		    }
		}
	      }
	      ctr1++;
	      if(ctr1>cleanJetAK4.size()-4){break;}
	    }
	    ctr2++;
	    if(ctr2>cleanJetAK4.size()-4){break;}
	  }
	  ctr3++;
	  if(ctr3>cleanJetAK4.size()-4){break;} 
	}
      }
      ctr4++;
      if(ctr4>=t_jetWJetP4->size()-1){break;}
    }
    //if(chi>240)cout<<"alert: "<<chi<<" "<<TPfl.M()<<" "<<TPfh.M()<<endl;
    if(chi>240. && chi!=10000)chi=239.;
    chi_output.push_back(chi);
    chi_output.push_back(TPfl.M());
    chi_output.push_back(TPfh.M());
    return chi_output;
  }
  else if(NrecoTop>=1 && NrecoV==1 && cleanJetAK4.size()>=3){
    for(unsigned int itT=0; itT <t_jetTopJetP4->size(); ++itT) {
      //cout<<"=======================================\n";
      //cout<<"size "<<cleanJetAK4.size()<<endl;
		     
      for(vector<int>::iterator itAK1 = cleanJetAK4.begin() ; ; ++itAK1) {
	//cout<<"itak1 : "<<*itAK1<<endl;
	for (vector<int>::iterator itAK2 = itAK1+1; ; ++itAK2){
	  //cout<<"\n itak2 "<<*itAK2<<endl;
	  for (vector<int>::iterator itAK3 = itAK2+1; itAK3 != cleanJetAK4.end(); ++itAK3){  		
	    //cout<<"itaK3 "<<*itAK3<<" ";
	    chi_1=10000;
	    chi_2=10000;
	    chi_3=10000;
	    Wb1=(*t_jetWJetP4)[0]+(*t_jetAK4P4)[*itAK1];
	    Wb2=(*t_jetWJetP4)[0]+(*t_jetAK4P4)[*itAK2];
	    Wb=(*t_jetWJetP4)[0]+(*t_jetAK4P4)[*itAK3];
	    V1=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3];
	    V2=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK1];
	    V=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK2];
			 
	    //cout<<"ooh\n";
	    if((*t_jetAK4CSV)[*itAK1]){
	      chi_1=pow((Wb1.M()-184)/(18.1),2.)+pow(((*t_jetTopJet_SoftDropMass)[itT]-179.7)/(18.4),2.)+pow((V1.M()-88.5)/(9.5),2.);  
	      Wbc=V1+(*t_jetAK4P4)[*itAK1];
	      if(V1.M()>65. && V1.M()<105.)chi_V2[0]=pow(((*t_jetTopJet_SoftDropMass)[itT]-179.7)/(18.4),2.)+pow(((*t_jetWJet_MassPruned)[0]-88.5)/(9.5),2.)+pow((Wbc.M()-171.3)/(17.7),2.);
	    }
	    if((*t_jetAK4CSV)[*itAK2]){
	      chi_2=pow((Wb2.M()-184)/(18.1),2.)+pow(((*t_jetTopJet_SoftDropMass)[itT]-179.7)/(18.4),2.)+pow((V2.M()-88.5)/(9.5),2.);
	      Wbc=V2+(*t_jetAK4P4)[*itAK2];
	      if(V2.M()>65. && V2.M()<105.)chi_V2[1]=pow(((*t_jetTopJet_SoftDropMass)[itT]-179.7)/(18.4),2.)+pow(((*t_jetWJet_MassPruned)[0]-88.5)/(9.5),2.)+pow((Wbc.M()-171.3)/(17.7),2.);

	    }
	    if((*t_jetAK4CSV)[*itAK3]){
	      chi_3=pow((Wb.M()-184)/(18.1),2.)+pow(((*t_jetTopJet_SoftDropMass)[itT]-179.7)/(18.4),2.)+pow((V.M()-88.5)/(9.5),2.);
	      Wbc=V+(*t_jetAK4P4)[*itAK3];
	      if(V.M()>65. && V.M()<105.)chi_V2[2]=pow(((*t_jetTopJet_SoftDropMass)[itT]-179.7)/(18.4),2.)+pow(((*t_jetWJet_MassPruned)[0]-88.5)/(9.5),2.)+pow((Wbc.M()-172.)/(17.6),2.);
	    }
	    //cout<<chi_1<<" "<<chi_2<<" "<<chi_3<<endl;
	    if(chi_1<chi && chi_1<chi_2 && chi_1<chi_3 && chi_1<chi_V2[0] && chi_1<chi_V2[1] && chi_1<chi_V2[2]){
	      chi=chi_1;
	      diff1=fabs((v1+Wb1).M()-(V1+(*t_jetTopJetP4)[itT]).M());
	      diff2=fabs((v1 + (*t_jetTopJetP4)[itT]).M()-(V1+Wb1).M());
	      if(diff1<diff2){TPfl=v1+Wb1;TPfh=V1+(*t_jetTopJetP4)[itT];}
	      else {TPfl=v1+(*t_jetTopJetP4)[itT];TPfh=V1+Wb1;}
	    }
	    else if(chi_2<chi && chi_2<chi_1 && chi_2<chi_3 && chi_2<chi_V2[0] && chi_2<chi_V2[1] && chi_2<chi_V2[2]){
	      chi=chi_2;
	      diff1=fabs((v1+Wb2).M()-(V2+(*t_jetTopJetP4)[itT]).M());
	      diff2=fabs((v1 + (*t_jetTopJetP4)[itT]).M()-(V2+Wb2).M());
	      if(diff1<diff2){TPfl=v1+Wb2;TPfh=V2+(*t_jetTopJetP4)[itT];}
	      else {TPfl=v1+(*t_jetTopJetP4)[itT];TPfh=V2+Wb2;}
	    }
	    else if(chi_3<chi && chi_3<chi_1 && chi_3<chi_2 && chi_2<chi_V2[0] && chi_2<chi_V2[1] && chi_2<chi_V2[2]){
	      chi=chi_3;
	      diff1=fabs((v1+Wb2).M()-(V+(*t_jetTopJetP4)[itT]).M());
	      diff2=fabs((v1 + (*t_jetTopJetP4)[itT]).M()-(V+Wb2).M());
	      if(diff1<diff2){TPfl=v1+Wb;TPfh=V+(*t_jetTopJetP4)[itT];}
	      else {TPfl=v1+(*t_jetTopJetP4)[itT];TPfh=V+Wb;}
	    }
	    else if(chi_V2[0]<chi && chi_V2[0]<chi_1 && chi_V2[0]<chi_2 && chi_V2[0]<chi_3 && chi_V2[0]<chi_V2[1] && chi_V2[0]<chi_V2[2]){
	      chi=chi_V2[0];
	      Wbc=V1+(*t_jetAK4P4)[*itAK1];
	      diff1=fabs((v1+Wbc).M()-((*t_jetWJetP4)[0]+(*t_jetTopJetP4)[itT]).M());
	      diff2=fabs((v1 + (*t_jetTopJetP4)[itT]).M()-((*t_jetWJetP4)[0]+Wbc).M());
	      if(diff1<diff2){TPfl=v1+Wbc;TPfh=V+(*t_jetTopJetP4)[itT];}
	      else {TPfl=v1+(*t_jetTopJetP4)[itT];TPfh=V+Wb;}
	    }
	    else if(chi_V2[1]<chi && chi_V2[1]<chi_1 && chi_V2[1]<chi_2 && chi_V2[1]<chi_3 && chi_V2[1]<chi_V2[0] && chi_V2[1]<chi_V2[2]){
	      chi=chi_V2[1];
	      Wbc=V2+(*t_jetAK4P4)[*itAK2];
	      diff1=fabs((v1+Wbc).M()-((*t_jetWJetP4)[0]+(*t_jetTopJetP4)[itT]).M());
	      diff2=fabs((v1 + (*t_jetTopJetP4)[itT]).M()-((*t_jetWJetP4)[0]+Wbc).M());
	      if(diff1<diff2){TPfl=v1+Wbc;TPfh=(*t_jetWJetP4)[0]+(*t_jetTopJetP4)[itT];}
	      else {TPfl=v1+(*t_jetTopJetP4)[itT];TPfh=(*t_jetWJetP4)[0]+Wbc;}
	    }
	    else if(chi_V2[2]<chi && chi_V2[2]<chi_1 && chi_V2[2]<chi_2 && chi_V2[2]<chi_3 && chi_V2[2]<chi_V2[0] && chi_V2[2]<chi_V2[1]){
	      chi=chi_V2[2];
	      Wbc=V+(*t_jetAK4P4)[*itAK3];
	      diff1=fabs((v1+Wbc).M()-((*t_jetWJetP4)[0]+(*t_jetTopJetP4)[itT]).M());
	      diff2=fabs((v1 + (*t_jetTopJetP4)[itT]).M()-((*t_jetWJetP4)[0]+Wbc).M());
	      if(diff1<diff2){TPfl=v1+Wbc;TPfh=(*t_jetWJetP4)[0]+(*t_jetTopJetP4)[itT];}
	      else {TPfl=v1+(*t_jetTopJetP4)[itT];TPfh=(*t_jetWJetP4)[0]+Wbc;}
	    }
	  }
	  ctr1++;
	  if(ctr1>cleanJetAK4.size()-3){break;}
	}
	ctr2++;
	if(ctr2>cleanJetAK4.size()-3){break;}
      }
    }
    //cout<<" Jets size: "<<cleanJetAK4.size()<<" ctr1: "<<ctr1<<" ctr2: "<<ctr2<<endl;
    if(chi>240. && chi!=10000)chi=239.;
    //if(TPfl.M()<50.0 || TPfh.M()<50.)cout<<chi<<" 2 "<<endl;
    chi_output.push_back(chi);
    chi_output.push_back(TPfl.M());
    chi_output.push_back(TPfh.M());
    return chi_output;
  }
  else if(NrecoTop==1 && !NrecoV && cleanJetAK4.size()>=5){		 
    for (vector<int>::iterator itAK1 = cleanJetAK4.begin() ; ; ++itAK1){
      for (vector<int>::iterator itAK2 = itAK1+1 ; ; ++itAK2){
	for (vector<int>::iterator itAK3 = itAK2+1 ; ; ++itAK3){
	  for (vector<int>::iterator itAK4 = itAK3+1 ; ; ++itAK4){
	    for (vector<int>::iterator itAK5 = itAK4+1 ;itAK5!= cleanJetAK4.end() ; ++itAK5){	   

	      //cout<<cleanJetAK4.size()<<endl;
	      //cout<<*itAK1<<" "<<*itAK2<<" "<<*itAK3<<" "<<*itAK4<<" "<<*itAK5<<endl;
	      if((*t_jetAK4CSV)[*itAK3]>0.5426 || (*t_jetAK4CSV)[*itAK4]>0.5426 || (*t_jetAK4CSV)[*itAK5]>0.5426 )chi_V2[0]=chi_cal(*itAK1,*itAK2,*itAK3,*itAK4,*itAK5);   
	      if((*t_jetAK4CSV)[*itAK2]>0.5426 || (*t_jetAK4CSV)[*itAK4]>0.5426 || (*t_jetAK4CSV)[*itAK5]>0.5426 )chi_V2[1]=chi_cal(*itAK1,*itAK3,*itAK2,*itAK4,*itAK5);
	      if((*t_jetAK4CSV)[*itAK3]>0.5426 || (*t_jetAK4CSV)[*itAK2]>0.5426 || (*t_jetAK4CSV)[*itAK5]>0.5426 )chi_V2[2]=chi_cal(*itAK1,*itAK4,*itAK3,*itAK2,*itAK5);
	      if((*t_jetAK4CSV)[*itAK3]>0.5426 || (*t_jetAK4CSV)[*itAK4]>0.5426 || (*t_jetAK4CSV)[*itAK2]>0.5426 )chi_V2[3]=chi_cal(*itAK1,*itAK5,*itAK3,*itAK4,*itAK2);
	      if((*t_jetAK4CSV)[*itAK1]>0.5426 || (*t_jetAK4CSV)[*itAK4]>0.5426 || (*t_jetAK4CSV)[*itAK5]>0.5426 )chi_V2[4]=chi_cal(*itAK3,*itAK2,*itAK1,*itAK4,*itAK5);
	      if((*t_jetAK4CSV)[*itAK1]>0.5426 || (*t_jetAK4CSV)[*itAK3]>0.5426 || (*t_jetAK4CSV)[*itAK5]>0.5426 )chi_V2[5]=chi_cal(*itAK4,*itAK2,*itAK1,*itAK3,*itAK5);
	      if((*t_jetAK4CSV)[*itAK1]>0.5426 || (*t_jetAK4CSV)[*itAK4]>0.5426 || (*t_jetAK4CSV)[*itAK3]>0.5426 )chi_V2[6]=chi_cal( *itAK5,*itAK2,*itAK1,*itAK4,*itAK3);
	      if((*t_jetAK4CSV)[*itAK1]>0.5426 || (*t_jetAK4CSV)[*itAK2]>0.5426 || (*t_jetAK4CSV)[*itAK5]>0.5426 )chi_V2[7]=chi_cal( *itAK3,*itAK4,*itAK1,*itAK2,*itAK5);
	      if((*t_jetAK4CSV)[*itAK1]>0.5426 || (*t_jetAK4CSV)[*itAK2]>0.5426 || (*t_jetAK4CSV)[*itAK4]>0.5426 )chi_V2[8]=chi_cal(*itAK3,*itAK5,*itAK1,*itAK2,*itAK4);
	      if((*t_jetAK4CSV)[*itAK1]>0.5426 || (*t_jetAK4CSV)[*itAK2]>0.5426 || (*t_jetAK4CSV)[*itAK3]>0.5426 )chi_V2[9]=chi_cal(*itAK5,*itAK4,*itAK1,*itAK2,*itAK3);
	      //for(int i=0;i<10;i++)cout<<chi_V2[i]<<" ";
	      //cout<<"\n";
	      for(int i=0;i<10;i++)chi_idx[i]=i;
	      //for(int i=0;i<10;i++)cout<<chi_idx[i]<<" ";
	      // cout<<"\n";
	      for(int i=1;i<10;++i)
		{
		  for(int j=0;j<(10-i);++j)
		    if(chi_V2[j]>chi_V2[j+1])
		      {
			temp=chi_V2[j];
			chi_V2[j]=chi_V2[j+1];
			chi_V2[j+1]=temp;
				     
			temp=chi_idx[j];
			chi_idx[j]=chi_idx[j+1];
			chi_idx[j+1]=temp;
		      }
		}
	      //for(int i=0;i<10;i++)cout<<chi_V2[i]<<" ";
	      //cout<<"\n"<<chi_idx[0]<<endl;
	      //cout<<"=================\n";
	      //cout<<V2[0]<<" "<<chi_idx[0]<<endl;
	      if(chi_V2[0]<chi){
		chi=chi_V2[0];
		switch(chi_idx[0])
		  {
		  case 0 :
		    Wb1=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5];
		    V=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK2];
		    diff1=fabs((v1+Wb1).M()-(V+(*t_jetTopJetP4)[0]).M());
		    diff2=fabs((v1+(*t_jetTopJetP4)[0]).M()-(V+Wb1).M());
		    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+(*t_jetTopJetP4)[0];}
		    else {TPfl=v1+(*t_jetTopJetP4)[0];TPfh=V+Wb1;}
		    break;
		  case 1:
		    Wb1=(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5];
		    V=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK3];
		    diff1=fabs((v1+Wb1).M()-(V+(*t_jetTopJetP4)[0]).M());
		    diff2=fabs((v1+(*t_jetTopJetP4)[0]).M()-(V+Wb1).M());
		    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+(*t_jetTopJetP4)[0];}
		    else {TPfl=v1+(*t_jetTopJetP4)[0];TPfh=V+Wb1;}
		    break;
		  case 2:
		    Wb1=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK5];
		    V=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK4];
		    diff1=fabs((v1+Wb1).M()-(V+(*t_jetTopJetP4)[0]).M());
		    diff2=fabs((v1+(*t_jetTopJetP4)[0]).M()-(V+Wb1).M());
		    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+(*t_jetTopJetP4)[0];}
		    else {TPfl=v1+(*t_jetTopJetP4)[0];TPfh=V+Wb1;}
		    break;
		  case 3 :
		    Wb1=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK2];
		    V=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK5];
		    diff1=fabs((v1+Wb1).M()-(V+(*t_jetTopJetP4)[0]).M());
		    diff2=fabs((v1+(*t_jetTopJetP4)[0]).M()-(V+Wb1).M());
		    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+(*t_jetTopJetP4)[0];}
		    else {TPfl=v1+(*t_jetTopJetP4)[0];TPfh=V+Wb1;}
		    break;
		  case 4 :
		    Wb1=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK5];
		    V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK2];
		    diff1=fabs((v1+Wb1).M()-(V+(*t_jetTopJetP4)[0]).M());
		    diff2=fabs((v1+(*t_jetTopJetP4)[0]).M()-(V+Wb1).M());
		    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+(*t_jetTopJetP4)[0];}
		    else {TPfl=v1+(*t_jetTopJetP4)[0];TPfh=V+Wb1;}
		    break;
		  case 5 :
		    Wb1=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK5];
		    V=(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK2];
		    diff1=fabs((v1+Wb1).M()-(V+(*t_jetTopJetP4)[0]).M());
		    diff2=fabs((v1+(*t_jetTopJetP4)[0]).M()-(V+Wb1).M());
		    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+(*t_jetTopJetP4)[0];}
		    else {TPfl=v1+(*t_jetTopJetP4)[0];TPfh=V+Wb1;}
		    break;
		  case 6 :
		    Wb1=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4]+(*t_jetAK4P4)[*itAK1];
		    V=(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK2];
		    diff1=fabs((v1+Wb1).M()-(V+(*t_jetTopJetP4)[0]).M());
		    diff2=fabs((v1+(*t_jetTopJetP4)[0]).M()-(V+Wb1).M());
		    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+(*t_jetTopJetP4)[0];}
		    else {TPfl=v1+(*t_jetTopJetP4)[0];TPfh=V+Wb1;}
		    break;
		  case 7 :
		    Wb1=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK5];
		    V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK4];
		    diff1=fabs((v1+Wb1).M()-(V+(*t_jetTopJetP4)[0]).M());
		    diff2=fabs((v1+(*t_jetTopJetP4)[0]).M()-(V+Wb1).M());
		    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+(*t_jetTopJetP4)[0];}
		    else {TPfl=v1+(*t_jetTopJetP4)[0];TPfh=V+Wb1;}
		    break;
		  case 8 :
		    Wb1=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK4];
		    V=(*t_jetAK4P4)[*itAK3]+(*t_jetAK4P4)[*itAK5];
		    diff1=fabs((v1+Wb1).M()-(V+(*t_jetTopJetP4)[0]).M());
		    diff2=fabs((v1+(*t_jetTopJetP4)[0]).M()-(V+Wb1).M());
		    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+(*t_jetTopJetP4)[0];}
		    else {TPfl=v1+(*t_jetTopJetP4)[0];TPfh=V+Wb1;}
		    break;
		  case 9 :
		    Wb1=(*t_jetAK4P4)[*itAK1]+(*t_jetAK4P4)[*itAK2]+(*t_jetAK4P4)[*itAK3];
		    V=(*t_jetAK4P4)[*itAK5]+(*t_jetAK4P4)[*itAK4];
		    diff1=fabs((v1+Wb1).M()-(V+(*t_jetTopJetP4)[0]).M());
		    diff2=fabs((v1+(*t_jetTopJetP4)[0]).M()-(V+Wb1).M());
		    if(diff1<diff2){TPfl=v1+Wb1;TPfh=V+(*t_jetTopJetP4)[0];}
		    else {TPfl=v1+(*t_jetTopJetP4)[0];TPfh=V+Wb1;}
		    break;
		  }
	      } 
	    }
	    ctr1++;
	    if(ctr1>cleanJetAK4.size()-5){break;}
	  }
	  ctr2++;
	  if(ctr2>cleanJetAK4.size()-5){break;}
	}
	ctr3++;
	if(ctr3>cleanJetAK4.size()-5){break;}
      }
      ctr4++;
      if(ctr4>cleanJetAK4.size()-5){break;}
    }
    if(chi>240. && chi!=10000)chi=239.;
		 
    //h_Chi2->Fill(chi);
    //h_Chi2_0t2V->Fill(chi);
    //if(TPfl.M()<50.0 || TPfh.M()<50.)cout<<chi<<" 2 "<<endl;
    chi_output.push_back(chi);
    chi_output.push_back(TPfl.M());
    chi_output.push_back(TPfh.M());
    return chi_output;	       
  }
  chi_output.push_back(10000);
  chi_output.push_back(10000);
  chi_output.push_back(10000);
  return chi_output; 
}

#endif // #ifdef TPrimeAna_cxx

/*
Inherited class analysis example:
Shows how to define an analysis class. Inheritance occurs through mother class RazorTuplizer
This example shows the structure of the inheritance but runs all the methods from the mother class.
*/

#include "RazorTuplizer.h"

class RazorAna : public RazorTuplizer{
public:
  //analyzer constructor and destructor
  explicit RazorAna(const edm::ParameterSet&);
  ~RazorAna();
  
  void resetBranches();
  void setBranches();
  
  void enableEventInfoBranches();
  void enablePileUpBranches();
  void enableMuonBranches();
  void enableElectronBranches();
  void enableTauBranches();
  void enablePhotonBranches();
  void enableJetBranches();
  void enableJetAK8Branches();
  void enableMetBranches();
  void enableRazorBranches();
  void enableGenParticles();

  //Re-defining select objects and fill tree branches 
  bool fillEventInfo(const edm::Event& iEvent);
  bool fillPileUp();//Fill summary PU info
  bool fillMuons();//Fills muon 4-momentum only. PT > 5GeV
  bool fillElectrons();//Fills Ele 4-momentum only. PT > 5GeV
  bool fillTaus();//Fills Tau 4-momentum only. PT > 20GeV
  bool fillPhotons();//Fills photon 4-momentum only. PT > 20GeV && ISO < 0.3
  bool fillJets();//Fills AK4 Jet 4-momentum, CSV, and CISV. PT > 20GeV
  bool fillJetsAK8();//Fills AK5 Jet 4-momentum.
  bool fillMet();//Fills MET(mag, phi)
  bool fillRazor();//Fills MR and RSQ
  bool fillGenParticles();

  bool isAncestor(const reco::Candidate*, const reco::Candidate*);

protected:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&);
  
  //PU
  int nBunchXing;
  int BunchXing[99];
  int nPU[99];
  float nPUmean[99];

  //Mu
  int muonCharge[99];//muon charge
  bool muonIsLoose[99];
  bool muonIsTight[99];
  float muon_d0[99];//transverse impact paramenter
  float muon_dZ[99];//impact parameter
  float muon_ip3d[99];//3d impact paramenter
  float muon_ip3dSignificance[99];//3d impact paramenter/error
  unsigned int muonType[99];//muonTypeBit: global, tracker, standalone 
  float muon_sumChargedHadronPt[99];//pfISO dr04
  float muon_sumChargedParticlePt[99];//pfISO dr04
  float muon_sumNeutralHadronEt[99];//pfISO dr04
  float muon_sumPhotonEt[99];//pfISO dr04
  
  //Ele
  float eleCharge[99];
  float eleE_SC[99];
  //float SC_ElePt[99]; 
  float eleEta_SC[99];
  float elePhi_SC[99];
  float eleSigmaIetaIeta[99];
  float eleFull5x5SigmaIetaIeta[99];
  float eleR9[99];
  float ele_dEta[99];
  float ele_dPhi[99];
  float ele_HoverE[99];
  float ele_d0[99];
  float ele_dZ[99];
  float ele_sumChargedHadronPt[99];
  float ele_sumNeutralHadronEt[99];
  float ele_sumPhotonEt[99];
  int ele_MissHits[99];
  int ele_ConvRejec[99];
  float ele_OneOverEminusOneOverP[99];
  float ele_RegressionE[99];
  float ele_CombineP4[99];
  
  //Taus

  //Photons
  float phoSigmaIetaIeta[99];
  float phoFull5x5SigmaIetaIeta[99];
  float phoR9[99];
  float pho_HoverE[99];
  float pho_sumChargedHadronPt[99];
  float pho_sumNeutralHadronEt[99];
  float pho_sumPhotonEt[99];
  int pho_isConversion[99];
  float pho_RegressionE[99];
  float pho_IDMVA[99];
  
  //Jets
  float jetMass[99];
  float jetJetArea[99];
  float jetPileupE[99];  
  
  //Event Info
  float pvX;
  float pvY;
  float pvZ;

  //MET
  float sumMET;
  float UncMETdpx;
  float UncMETdpy;
  float UncMETdSumEt;

  //Gen Info
  unsigned int nGenParticle;
  int motherIndex[99];
  int motherId[99];
  int gParticleId[99];
  float gParticleE[99];
  float gParticlePt[99];
  float gParticleEta[99];
  float gParticlePhi[99];
  float gParticleVx[99];
  float gParticleVy[99];
  float gParticleVz[99];

};

//define this as a plug-in
DEFINE_FWK_MODULE(RazorAna);

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
  
  void enableMuonBranches();
  void enableElectronBranches();
  void enableTauBranches();
  void enablePhotonBranches();
  void enableJetBranches();
  void enableJetAK8Branches();
  void enableMetBranches();
  void enableRazorBranches();

  //Re-defining select objects and fill tree branches 
  bool fillMuons();//Fills muon 4-momentum only. PT > 5GeV
  bool fillElectrons();//Fills Ele 4-momentum only. PT > 5GeV
  bool fillTaus();//Fills Tau 4-momentum only. PT > 20GeV
  bool fillPhotons();//Fills photon 4-momentum only. PT > 20GeV && ISO < 0.3
  bool fillJets();//Fills AK4 Jet 4-momentum, CSV, and CISV. PT > 20GeV
  bool fillJetsAK8();//Fills AK5 Jet 4-momentum.
  bool fillMet();//Fills MET(mag, phi)
  bool fillRazor();//Fills MR and RSQ

protected:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&);
  
  //Mu
  float muonIsLoose[99];
  float muonIsTight[99];
  
  //Ele
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
  float ele_MissHits[99];
  float ele_ConvRejec[99];
  float ele_OneOverEminusOneOverP[99];
  float ele_RegressionE[99];
  float ele_CombineP4[99];
};

//define this as a plug-in
DEFINE_FWK_MODULE(RazorAna);

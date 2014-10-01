/*
Example of inherited class from class RazorTuplizer
Defines just constructor and destructor
*/

#include "RazorInhExa.hh"

//Defines Constructor
RazorAna::RazorAna(const edm::ParameterSet& iConfig) : RazorTuplizer(iConfig){};
//Defines Destructor
RazorAna::~RazorAna(){};
//Define custom analyzer
void RazorAna::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;
  
  //initialize
  resetBranches();
  loadEvent(iEvent); //loads objects and resets tree branches

  bool isGoodEvent =
    fillEventInfo(iEvent)
    && fillMuons()
    && fillElectrons()
    && fillTaus()
    && fillPhotons()
    && fillJets()
    && fillJetsAK8()
    && fillMet()
    && fillRazor();

  //fill the tree if the event wasn't rejected
  if(isGoodEvent) outputTree->Fill();
};

void RazorAna::resetBranches(){
  RazorTuplizer::resetBranches();
  //Re-set newly defined variables;
  for(int j = 0; j < 99; j++){
    //Mu
    muonIsLoose[j] = -99.0;
    muonIsTight[j] = -99.0;
    
    //Ele
    eleE_SC[j] = -99.0;
    //elePt_SC[j] = -99.0;
    eleEta_SC[j] = -99.0;
    elePhi_SC[j] = -99.0;
    eleSigmaIetaIeta[j] = -99.0;
    eleFull5x5SigmaIetaIeta[j] = -99.0;
    eleR9[j] = -99;
    ele_dEta[j] = -99;
    ele_dPhi[j] = -99;
    ele_HoverE[j] = -99;
    ele_d0[j] = -99;
    ele_dZ[j] = -99;
    ele_sumChargedHadronPt[j] = -99.0;
    ele_sumNeutralHadronEt[j] = -99.0;
    ele_sumPhotonEt[j] = -99.0;
    ele_MissHits[j] = -99;
    ele_ConvRejec[j] = -99;
    ele_OneOverEminusOneOverP[j] = -99;
    ele_RegressionE[j] = -99;
    ele_CombineP4[j] = -99;
  }
};

void RazorAna::setBranches(){
  enableEventInfoBranches();
  enableMuonBranches();
  enableElectronBranches();
  enableTauBranches();
  enablePhotonBranches();
  enableJetBranches();
  enableJetAK8Branches();
  enableMetBranches();
  enableRazorBranches();
};

void RazorAna::enableMuonBranches(){
  RazorTuplizer::enableMuonBranches();
  outputTree->Branch("muonIsLoose", muonIsLoose,"muonIsLoose[nMuons]/F");
  outputTree->Branch("muonIsTight", muonIsTight,"muonIsTight[nMuons]/F");
};

void RazorAna::enableElectronBranches(){
  RazorTuplizer::enableElectronBranches();
  outputTree->Branch("eleCharge", eleCharge, "eleCharge[nElectrons]/F");
  outputTree->Branch("EleE_SC", eleE_SC,"eleE_SC[nElectrons]/F");
  //outputTree->Branch("SC_ElePt", SC_ElePt,"SC_ElePt[nElectrons]/F");
  outputTree->Branch("eleEta_SC", eleEta_SC,"eleEta_SC[nElectrons]/F");
  outputTree->Branch("elePhi_SC", elePhi_SC,"elePhi_SC[nElectrons]/F");
  outputTree->Branch("eleSigmaIetaIeta", eleSigmaIetaIeta, "eleSigmaIetaIeta[nElectrons]/F");
  outputTree->Branch("eleFull5x5SigmaIetaIeta", eleFull5x5SigmaIetaIeta, "eleFull5x5SigmaIetaIeta[nElectrons]/F");
  outputTree->Branch("eleR9", eleR9, "eleR9[nElectrons]/F");
  outputTree->Branch("ele_dEta", ele_dEta, "ele_dEta[nElectrons]/F");
  outputTree->Branch("ele_dPhi", ele_dPhi, "ele_dPhi[nElectrons]/F");
  outputTree->Branch("ele_HoverE", ele_HoverE, "ele_HoverE[nElectrons]/F");
  outputTree->Branch("ele_d0", ele_d0, "ele_d0[nElectrons]/F");
  outputTree->Branch("ele_dZ", ele_dZ, "ele_dZ[nElectrons]/F");
  outputTree->Branch("ele_sumChargedHadronPt", ele_sumChargedHadronPt, "ele_sumChargedHadronPt[nElectrons]/F");
  outputTree->Branch("ele_sumNeutralHadronEt", ele_sumNeutralHadronEt, "ele_sumNeutralHadronEt[nElectrons]/F");
  outputTree->Branch("ele_sumPhotonEt", ele_sumPhotonEt, "ele_sumPhotonEt[nElectrons]/F");
  outputTree->Branch("ele_MissHits", ele_MissHits, "ele_MissHits[nElectrons]/F");
  outputTree->Branch("ele_ConvRejec", ele_ConvRejec, "ele_ConvRejec[nElectrons]/I");
  outputTree->Branch("ele_OneOverEminusOneOverP", ele_OneOverEminusOneOverP, "ele_OneOverEminusOneOverP[nElectrons]/F");
  outputTree->Branch("ele_RegressionE", ele_RegressionE, "ele_RegressionE[nElectrons]/F");
  outputTree->Branch("ele_TrackRegressionE", ele_TrackRegressionE, "ele_TrackRegressionE[nElectrons]/F");
  //outputTree->Branch("", , "[nElectrons]/F");
}

void RazorAna::enableTauBranches(){
  RazorTuplizer::enableTauBranches();
};

void RazorAna::enablePhotonBranches(){
  RazorTuplizer::enablePhotonBranches();
};

void RazorAna::enableJetBranches(){
  RazorTuplizer::enableJetBranches();
};

void RazorAna::enableJetAK8Branches(){
  RazorTuplizer::enableJetAK8Branches();
};

void RazorAna::enableMetBranches(){
  RazorTuplizer::enableMetBranches();
}

void RazorAna::enableRazorBranches(){
  RazorTuplizer::enableRazorBranches();
}

/*
Re-defining Fill methods (will not use mother class methods)
*/
bool RazorAna::fillMuons(){
  //PV required for Tight working point
  const reco::Vertex &PV = vertices->front();
  for(const pat::Muon &mu : *muons){
    if(mu.pt() < 5) continue;
    muonE[nMuons] = mu.energy();
    muonPt[nMuons] = mu.pt();
    muonEta[nMuons] = mu.eta();
    muonPhi[nMuons] = mu.phi();
    muonIsLoose[nMuons] = mu.isLooseMuon();
    muonIsTight[nMuons] = mu.isTightMuon(PV);
    nMuons++;
  }

  return true;
};

bool RazorAna::fillElectrons(){
  for(const pat::Electron &ele : *electrons){
    if(ele.pt() < 5) continue;
    eleE[nElectrons] = ele.energy();
    elePt[nElectrons] = ele.pt();
    eleEta[nElectrons] = ele.eta();
    elePhi[nElectrons] = ele.phi();
    eleCharge[nElectrons] = ele.gsfTrack()->charge();
    eleE_SC[nElectrons] = ele.superCluster()->energy();
    //SC_ElePt[nElectrons] = ele.superCluster()->pt();
    eleEta_SC[nElectrons] = ele.superCluster()->eta();
    elePhi_SC[nElectrons] = ele.superCluster()->phi();
    eleSigmaIetaIeta[nElectrons] = ele.sigmaIetaIeta();
    eleFull5x5SigmaIetaIeta[nElectrons] = ele.full5x5_sigmaIetaIeta();
    eleR9[nElectrons] = ele.r9();
    ele_dEta[nElectrons] = ele.deltaEtaSuperClusterTrackAtVtx();
    ele_dPhi[nElectrons] = ele.deltaPhiSuperClusterTrackAtVtx();
    ele_HoverE[nElectrons] = ele.hcalOverEcalBc();
    ele_sumChargedHadronPt[nElectrons] = ele.pfIsolationVariables().sumChargedHadronPt;
    ele_sumNeutralHadronEt[nElectrons] = ele.pfIsolationVariables().sumNeutralHadronEt;
    ele_sumPhotonEt[nElectrons] = ele.pfIsolationVariables().sumPhotonEt;
    ele_MissHits[nElectrons] = ele.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits();
    ele_ConvRejec[nElectrons] = ele.convFlags();
    ele_OneOverEminusOneOverP[nElectrons] = 1./ele.correctedEcalEnergy()  -  1./ele.trackMomentumAtVtx().R();
    ele_RegressionE[nElectrons] = ele.ecalRegressionEnergy();
    ele_TrackRegressionE[nElectrons] = ele.ecalTrackRegressionEnergy();
    //ele_CombineP4[nElectrons] = ele.p4(2).E();
    nElectrons++;
  }
  
  return true;
};

bool RazorAna::fillTaus(){
  for (const pat::Tau &tau : *taus) {
    if (tau.pt() < 20) continue;
    tauE[nTaus] = tau.energy();
    tauPt[nTaus] = tau.pt();
    tauEta[nTaus] = tau.eta();
    tauPhi[nTaus] = tau.phi();
    nTaus++;
  }

  return true;
};

bool RazorAna::fillPhotons(){
  for (const pat::Photon &pho : *photons) {
    if (pho.pt() < 20 or pho.chargedHadronIso()/pho.pt() > 0.3) continue;
    phoE[nPhotons] = pho.energy();
    phoPt[nPhotons] = pho.pt();
    phoEta[nPhotons] = pho.eta();
    phoPhi[nPhotons] = pho.phi();
    nPhotons++;
  }

  return true;
};

bool RazorAna::fillJets(){
  for (const pat::Jet &j : *jets) {
    if (j.pt() < 20) continue;
    jetE[nJets] = j.energy();
    jetPt[nJets] = j.pt();
    jetEta[nJets] = j.eta();
    jetPhi[nJets] = j.phi();
    jetCSV[nJets] = j.bDiscriminator("combinedSecondaryVertexBJetTags");
    jetCISV[nJets] = j.bDiscriminator("combinedInclusiveSecondaryVertexBJetTags");
    nJets++;
  }

  return true;
};

bool RazorAna::fillJetsAK8(){
  for (const pat::Jet &j : *jetsAK8) {
    fatJetE[nFatJets] = j.energy();
    fatJetPt[nFatJets] = j.pt();
    fatJetEta[nFatJets] = j.eta();
    fatJetPhi[nFatJets] = j.phi();
    nFatJets++;
  }

  return true;
};

bool RazorAna::fillMet(){
  const pat::MET &Met = mets->front();
  metPt = Met.pt();
  metPhi = Met.phi();

  return true;
};

bool RazorAna::fillRazor(){
  //get good jets for razor calculation                                  
  vector<TLorentzVector> goodJets;
  for (const pat::Jet &j : *jets) {
    if (j.pt() < 20) continue;
    //add to goodJets vector                                          
    TLorentzVector newJet(j.px(), j.py(), j.pz(), j.energy());
    goodJets.push_back(newJet);
  }
  //MET                                                                                                        
  const pat::MET &Met = mets->front();
  
  //compute the razor variables using the selected jets and MET
  if(goodJets.size() > 1){
    vector<TLorentzVector> hemispheres = getHemispheres(goodJets);
    TLorentzVector pfMet(Met.px(), Met.py(), 0.0, 0.0);
    MR = computeMR(hemispheres[0], hemispheres[1]);
    RSQ = computeR2(hemispheres[0], hemispheres[1], pfMet);
  }

  return true;
};

//------ Method called once each job just before starting event loop ------//                                                               
void RazorAna::beginJob(){
  setBranches();
};

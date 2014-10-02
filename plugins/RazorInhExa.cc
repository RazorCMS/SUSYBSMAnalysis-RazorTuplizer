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
    muonCharge[j] = -99.0;
    muonIsLoose[j] = -99.0;
    muonIsTight[j] = -99.0;
    muon_d0[j] = -99.0;
    muon_dZ[j] = -99.0;
    muon_ip3d[j] = -99.0;
    muon_ip3dSignificance[j] = -99.0;
    muonType[j] = 0;
    muon_sumChargedHadronPt[j] = -99.0;
    muon_sumChargedParticlePt[j] = -99.0;
    muon_sumNeutralHadronEt[j] = -99.0;
    muon_sumPhotonEt[j] = -99.0;
    
    //Ele
    eleE_SC[j] = -99.0;
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

    //Taus
    
    //Photons
    phoSigmaIetaIeta[j] = -99.0;
    phoFull5x5SigmaIetaIeta[j] = -99.0;
    phoR9[j] = -99.0;
    pho_HoverE[j] = -99.0;
    pho_sumChargedHadronPt[j] = -99.0;
    pho_sumNeutralHadronEt[j] = -99.0;
    pho_sumPhotonEt[j] = -99.0;
    pho_isConversion[j] = -99;
    pho_RegressionE[j] = -99.0;
    pho_pfMVA[j] = -99.0;
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
  outputTree->Branch("muonCharge", muonCharge, "muonCharge[nMuons]/F");
  outputTree->Branch("muonIsLoose", muonIsLoose,"muonIsLoose[nMuons]/F");
  outputTree->Branch("muonIsTight", muonIsTight,"muonIsTight[nMuons]/F");
  outputTree->Branch("muon_d0", muon_d0, "muon_d0[nMuons]/F");
  outputTree->Branch("muon_dZ", muon_dZ, "muon_dZ[nMuons]/F");
  outputTree->Branch("muon_ip3d", muon_ip3d, "muon_ip3d[nMuons]/F");
  outputTree->Branch("muon_ip3dSignificance", muon_ip3dSignificance, "muon_ip3dSignificance[nMuons]/F");
  outputTree->Branch("muonType", muonType, "muonType[nMuons]/s");
  outputTree->Branch("muon_sumChargedHadronPt", muon_sumChargedHadronPt, "muon_sumChargedHadronPt[nMuons]/F");
  outputTree->Branch("muon_sumChargedParticlePt", muon_sumChargedParticlePt, "muon_sumChargedParticlePt[nMuons]/F");
  outputTree->Branch("muon_sumNeutralHadronEt", muon_sumNeutralHadronEt, "muon_sumNeutralHadronEt[nMuons]/F");
  outputTree->Branch("muon_sumPhotonEt", muon_sumPhotonEt, "muon_sumPhotonEt[nMuons]/F");
  //outputTree->Branch("", , "[nMuons]/F");
  //outputTree->Branch("", , "[nMuons]/F");
  //outputTree->Branch("", , "[nMuons]/F");
  //outputTree->Branch("", , "[nMuons]/F");
  //outputTree->Branch("", , "[nMuons]/F");
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
  outputTree->Branch("phoSigmaIetaIeta", phoSigmaIetaIeta, "phoSigmaIetaIeta[nPhotons]/F");
  outputTree->Branch("phoFull5x5SigmaIetaIeta", phoFull5x5SigmaIetaIeta, "phoFull5x5SigmaIetaIeta[nPhotons]/F");
  outputTree->Branch("phoR9", phoR9, "phoR9[nPhotons]/F");
  outputTree->Branch("pho_HoverE", pho_HoverE, "pho_HoverE[nPhotons]/F");
  outputTree->Branch("pho_sumChargedHadronPt", pho_sumChargedHadronPt, "pho_sumChargedHadronPt[nPhotons]/F");
  outputTree->Branch("pho_sumNeutralHadronEt", pho_sumNeutralHadronEt, "pho_sumNeutralHadronEt[nPhotons]/F");
  outputTree->Branch("pho_sumPhotonEt", pho_sumPhotonEt, "pho_sumPhotonEt[nPhotons]/F");
  outputTree->Branch("pho_isConversion", pho_isConversion, "pho_isConversion[nPhotons]/I");
  outputTree->Branch("pho_RegressionE", pho_RegressionE, "pho_RegressionE[nPhotons]/F");
  outputTree->Branch("pho_pfMVA", pho_pfMVA, "pho_pfMVA[nPhotons]/F");
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
    muonCharge[nMuons] = mu.charge();
    muonIsLoose[nMuons] = mu.isLooseMuon();
    muonIsTight[nMuons] = mu.isTightMuon(PV);
    muon_d0[nMuons] = -mu.muonBestTrack()->dxy(PV.position());
    muon_dZ[nMuons] = mu.muonBestTrack()->dz(PV.position());
    muon_ip3d[nMuons] = mu.dB(pat::Muon::PV3D);
    muon_ip3dSignificance[nMuons] = mu.dB(pat::Muon::PV3D)/mu.edB(pat::Muon::PV3D);
    muonType[nMuons] = mu.isMuon() || mu.isGlobalMuon();
    muon_sumChargedHadronPt[nMuons] = mu.pfIsolationR04().sumChargedHadronPt;
    muon_sumChargedParticlePt[nMuons] = mu.pfIsolationR04().sumChargedParticlePt;
    muon_sumNeutralHadronEt[nMuons] = mu.pfIsolationR04().sumNeutralHadronEt;
    muon_sumPhotonEt[nMuons] =  mu.pfIsolationR04().sumPhotonEt;
    nMuons++;
  }

  return true;
};

bool RazorAna::fillElectrons(){
  const reco::Vertex &PV = vertices->front();
  for(const pat::Electron &ele : *electrons){
    if(ele.pt() < 5) continue;
    eleE[nElectrons] = ele.energy();
    elePt[nElectrons] = ele.pt();
    eleEta[nElectrons] = ele.eta();
    elePhi[nElectrons] = ele.phi();
    //eleCharge[nElectrons] = ele.gsfTrack()->charge();
    eleCharge[nElectrons] = ele.charge();
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
    ele_d0[nElectrons] = -ele.gsfTrack().get()->dxy(PV.position());
    ele_dZ[nElectrons] = ele.gsfTrack().get()->dz(PV.position());
    ele_sumChargedHadronPt[nElectrons] = ele.pfIsolationVariables().sumChargedHadronPt;
    ele_sumNeutralHadronEt[nElectrons] = ele.pfIsolationVariables().sumNeutralHadronEt;
    ele_sumPhotonEt[nElectrons] = ele.pfIsolationVariables().sumPhotonEt;
    ele_MissHits[nElectrons] = ele.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits();
    ele_ConvRejec[nElectrons] = ele.convFlags();
    ele_OneOverEminusOneOverP[nElectrons] = 1./ele.correctedEcalEnergy()  -  1./ele.trackMomentumAtVtx().R();
    ele_RegressionE[nElectrons] = ele.ecalRegressionEnergy();
    ele_TrackRegressionE[nElectrons] = ele.ecalTrackRegressionEnergy();
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
    phoSigmaIetaIeta[nPhotons] = pho.sigmaIetaIeta();
    phoFull5x5SigmaIetaIeta[nPhotons] = pho.full5x5_sigmaIetaIeta();
    phoR9[nPhotons] = pho.r9();
    pho_sumChargedHadronPt[nPhotons] = pho.chargedHadronIso();
    pho_sumNeutralHadronEt[nPhotons] = pho.neutralHadronIso();
    pho_sumPhotonEt[nPhotons] = pho.photonIso();
    pho_isConversion[nPhotons] = pho.hasConversionTracks();
    pho_RegressionE[nPhotons] = pho.getCorrectedEnergy(reco::Photon::P4type::regression1);
    pho_pfMVA[nPhotons] = pho.pfMVA();
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

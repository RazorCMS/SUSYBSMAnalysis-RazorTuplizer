/*
Example of inherited class from class RazorTuplizer
Defines just constructor and destructor
*/

#include "RazorInhExa.hh"

//Defines Constructor
RazorAna::RazorAna(const edm::ParameterSet& iConfig) : RazorTuplizer(iConfig){

  std::vector<std::string> myTrigWeights;
  myTrigWeights.push_back(edm::FileInPath("EgammaAnalysis/ElectronTools/data/CSA14/TrigIDMVA_25ns_EB_BDT.weights.xml").fullPath().c_str());
  myTrigWeights.push_back(edm::FileInPath("EgammaAnalysis/ElectronTools/data/CSA14/TrigIDMVA_25ns_EE_BDT.weights.xml").fullPath().c_str());

  myMVATrig = new EGammaMvaEleEstimatorCSA14();
  myMVATrig->initialize("BDT",
			EGammaMvaEleEstimatorCSA14::kTrig,
			true,
			myTrigWeights);

  std::vector<std::string> myNonTrigWeights;
  myNonTrigWeights.push_back(edm::FileInPath("EgammaAnalysis/ElectronTools/data/CSA14/EIDmva_EB_5_25ns_BDT.weights.xml").fullPath().c_str());
  myNonTrigWeights.push_back(edm::FileInPath("EgammaAnalysis/ElectronTools/data/CSA14/EIDmva_EE_5_25ns_BDT.weights.xml").fullPath().c_str());
  myNonTrigWeights.push_back(edm::FileInPath("EgammaAnalysis/ElectronTools/data/CSA14/EIDmva_EB_10_25ns_BDT.weights.xml").fullPath().c_str());
  myNonTrigWeights.push_back(edm::FileInPath("EgammaAnalysis/ElectronTools/data/CSA14/EIDmva_EE_10_25ns_BDT.weights.xml").fullPath().c_str());
  
  myMVANonTrig = new EGammaMvaEleEstimatorCSA14();
  myMVANonTrig->initialize("BDT",
			EGammaMvaEleEstimatorCSA14::kNonTrig,
			true,
			myNonTrigWeights);

};

//Defines Destructor
RazorAna::~RazorAna(){};
//Define custom analyzer
void RazorAna::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;
  
  //initialize
  resetBranches();
  loadEvent(iEvent); //loads objects and resets tree branches

  NEvents->Fill(0); //increment event count

  bool isGoodEvent =
    fillEventInfo(iEvent)
    && fillPileUp()
    && fillMuons()
    && fillElectrons()
    && fillTaus()
    && fillPhotons()
    && fillJets()
    && fillJetsAK8()
    && fillMet()
    && fillRazor()
    && fillMC()
    && fillGenParticles();

  //fill the tree if the event wasn't rejected
  if(isGoodEvent) RazorEvents->Fill();
};

void RazorAna::resetBranches(){
  RazorTuplizer::resetBranches();
  //Re-set newly defined variables;
  //PU
  nBunchXing = 0;
  for(int j = 0; j < 99; j++){
    //PU
    BunchXing[j] = -99;
    nPU[j] = -99;
    nPUmean[j] = -99.0;
    
    //Mu
    muonCharge[j] = -99;
    muonIsLoose[j] = false;
    muonIsTight[j] = false;
    muon_d0[j] = -99.0;
    muon_dZ[j] = -99.0;
    muon_ip3d[j] = -99.0;
    muon_ip3dSignificance[j] = -99.0;
    muonType[j] = 0;
    muonQuality[j] = 0;
    muon_relIso04DBetaCorr[j] = -99.0;
    
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
    ele_relIsoDBetaCorr[j] = -99.0;
    ele_MissHits[j] = -99;
    ele_PassConvVeto[j] = false;
    ele_OneOverEminusOneOverP[j] = -99.0;
    ele_IDMVATrig[j] = -99.0;
    ele_IDMVANonTrig[j] = -99.0;
    ele_RegressionE[j] = -99.0;
    ele_CombineP4[j] = -99.0;

    //Taus
    tau_IsLoose[j] = false;
    tau_IsMedium[j] = false;
    tau_IsTight[j] = false;
    tau_passEleVetoLoose[j] = false;
    tau_passEleVetoMedium[j] = false;
    tau_passEleVetoTight[j] = false;
    tau_passMuVetoLoose[j] = false;
    tau_passMuVetoMedium[j] = false;
    tau_passMuVetoTight[j] = false;
    tau_ID[j] = 0;
    tau_combinedIsoDeltaBetaCorr3Hits[j] = -99.0;
    tau_eleVetoMVA[j] = -99.0;
    tau_eleVetoCategory[j] = -1;
    tau_muonVetoMVA[j] = -99.0;
    tau_isoMVAnewDMwLT[j] = -99.0;
    tau_isoMVAnewDMwoLT[j] = -99.0;
    tau_leadCandPt[j] = -99.0;
    tau_leadCandID[j] = 0;
    tau_leadChargedHadrCandPt[j] = -99.0;
    tau_leadChargedHadrCandID[j] = 0;

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
    pho_IDMVA[j] = -99.0;

    //Jets
    jetMass[j] =  -99.0;
    jetJetArea[j] = -99.0;
    jetPileupE[j] = -99.0;
    jetPileupId[j] = -99.0;

    //Event Info
    pvX = -99.0;
    pvY = -99.0;
    pvZ = -99.0;

    //MET
    sumMET = -99.0;
    UncMETdpx = -99.0;
    UncMETdpy = -99.0;
    UncMETdSumEt = -99.0;
    
    //GenInfo
    nGenParticle = 0;
    gParticleMotherId[j] = -99999;
    gParticleMotherIndex[j] = -99999;
    gParticleId[j] = -99999;
    gParticleStatus[j] = -99999;
    gParticleE[j] = -99999.0;
    gParticlePt[j] = -99999.0;
    gParticleEta[j] = -99999.0;
    gParticlePhi[j] = -99999.0;
  }
};

void RazorAna::setBranches(){
  enableEventInfoBranches();
  enablePileUpBranches();
  enableMuonBranches();
  enableElectronBranches();
  enableTauBranches();
  enablePhotonBranches();
  enableJetBranches();
  enableJetAK8Branches();
  enableMetBranches();
  enableRazorBranches();
  enableMCBranches();
  enableGenParticles();
  
};


void RazorAna::enableEventInfoBranches(){
  RazorTuplizer::enableEventInfoBranches();
  RazorEvents->Branch("pvX", &pvX, "pvX/F");
  RazorEvents->Branch("pvY", &pvY, "pvY/F");
  RazorEvents->Branch("pvZ", &pvZ, "pvZ/F");
};

void RazorAna::enablePileUpBranches(){
  RazorEvents->Branch("nBunchXing", &nBunchXing, "nBunchXing/I");
  RazorEvents->Branch("BunchXing", BunchXing, "BunchXing[nBunchXing]/I");
  RazorEvents->Branch("nPU", nPU, "nPU[nBunchXing]/I");
  RazorEvents->Branch("nPUmean", nPUmean, "nPUmean[nBunchXing]/F");
};

void RazorAna::enableMuonBranches(){
  RazorTuplizer::enableMuonBranches();
  RazorEvents->Branch("muonCharge", muonCharge, "muonCharge[nMuons]/I");
  RazorEvents->Branch("muonIsLoose", muonIsLoose,"muonIsLoose[nMuons]/O");
  RazorEvents->Branch("muonIsTight", muonIsTight,"muonIsTight[nMuons]/O");
  RazorEvents->Branch("muon_d0", muon_d0, "muon_d0[nMuons]/F");
  RazorEvents->Branch("muon_dZ", muon_dZ, "muon_dZ[nMuons]/F");
  RazorEvents->Branch("muon_ip3d", muon_ip3d, "muon_ip3d[nMuons]/F");
  RazorEvents->Branch("muon_ip3dSignificance", muon_ip3dSignificance, "muon_ip3dSignificance[nMuons]/F");
  RazorEvents->Branch("muonType", muonType, "muonType[nMuons]/s");
  RazorEvents->Branch("muonQuality", muonQuality, "muonQuality[nMuons]/i");
  RazorEvents->Branch("muon_relIso04DBetaCorr", muon_relIso04DBetaCorr, "muon_relIso04DBetaCorr[nMuons]/F");
};

void RazorAna::enableElectronBranches(){
  RazorTuplizer::enableElectronBranches();
  RazorEvents->Branch("eleCharge", eleCharge, "eleCharge[nElectrons]/F");
  //RazorEvents->Branch("EleE_SC", eleE_SC,"eleE_SC[nElectrons]/F");
  RazorEvents->Branch("eleEta_SC", eleEta_SC,"eleEta_SC[nElectrons]/F");
  //RazorEvents->Branch("elePhi_SC", elePhi_SC,"elePhi_SC[nElectrons]/F");
  RazorEvents->Branch("eleSigmaIetaIeta", eleSigmaIetaIeta, "eleSigmaIetaIeta[nElectrons]/F");
  RazorEvents->Branch("eleFull5x5SigmaIetaIeta", eleFull5x5SigmaIetaIeta, "eleFull5x5SigmaIetaIeta[nElectrons]/F");
  RazorEvents->Branch("eleR9", eleR9, "eleR9[nElectrons]/F");
  RazorEvents->Branch("ele_dEta", ele_dEta, "ele_dEta[nElectrons]/F");
  RazorEvents->Branch("ele_dPhi", ele_dPhi, "ele_dPhi[nElectrons]/F");
  RazorEvents->Branch("ele_HoverE", ele_HoverE, "ele_HoverE[nElectrons]/F");
  RazorEvents->Branch("ele_d0", ele_d0, "ele_d0[nElectrons]/F");
  RazorEvents->Branch("ele_dZ", ele_dZ, "ele_dZ[nElectrons]/F");
  RazorEvents->Branch("ele_relIsoDBetaCorr", ele_relIsoDBetaCorr, "ele_relIsoDBetaCorr[nElectrons]/F");
  RazorEvents->Branch("ele_MissHits", ele_MissHits, "ele_MissHits[nElectrons]/I");
  RazorEvents->Branch("ele_PassConvVeto", ele_PassConvVeto, "ele_PassConvVeto[nElectrons]/O");
  RazorEvents->Branch("ele_OneOverEminusOneOverP", ele_OneOverEminusOneOverP, "ele_OneOverEminusOneOverP[nElectrons]/F");
  RazorEvents->Branch("ele_IDMVATrig", ele_IDMVATrig, "ele_IDMVATrig[nElectrons]/F");
  RazorEvents->Branch("ele_IDMVANonTrig", ele_IDMVANonTrig, "ele_IDMVANonTrig[nElectrons]/F");
  RazorEvents->Branch("ele_RegressionE", ele_RegressionE, "ele_RegressionE[nElectrons]/F");
  RazorEvents->Branch("ele_CombineP4", ele_CombineP4, "ele_CombineP4[nElectrons]/F");
};

void RazorAna::enableTauBranches(){
  RazorTuplizer::enableTauBranches();
  RazorEvents->Branch("tau_IsLoose", tau_IsLoose, "tau_IsLoose[nTaus]/O");
  RazorEvents->Branch("tau_IsMedium", tau_IsMedium, "tau_IsMedium[nTaus]/O");
  RazorEvents->Branch("tau_IsTight", tau_IsTight, "tau_IsTight[nTaus]/O");
  RazorEvents->Branch("tau_passEleVetoLoose", tau_passEleVetoLoose, "tau_passEleVetoLoose[nTaus]/O");
  RazorEvents->Branch("tau_passEleVetoMedium", tau_passEleVetoMedium, "tau_passEleVetoMedium[nTaus]/O");
  RazorEvents->Branch("tau_passEleVetoTight", tau_passEleVetoTight, "tau_passEleVetoTight[nTaus]/O");
  RazorEvents->Branch("tau_passMuVetoLoose", tau_passMuVetoLoose, "tau_passMuVetoLoose[nTaus]/O");
  RazorEvents->Branch("tau_passMuVetoMedium", tau_passMuVetoMedium, "tau_passMuVetoMedium[nTaus]/O");
  RazorEvents->Branch("tau_passMuVetoTight", tau_passMuVetoTight, "tau_passMuVetoTight[nTaus]/O");
  RazorEvents->Branch("tau_ID", tau_ID, "tau_ID[nTaus]/i");
  RazorEvents->Branch("tau_combinedIsoDeltaBetaCorr3Hits", tau_combinedIsoDeltaBetaCorr3Hits, "tau_combinedIsoDeltaBetaCorr3Hits[nTaus]/F");
  RazorEvents->Branch("tau_eleVetoMVA", tau_eleVetoMVA, "tau_eleVetoMVA[nTaus]/F");
  RazorEvents->Branch("tau_eleVetoCategory", tau_eleVetoCategory, "tau_eleVetoCategory[nTaus]/I");
  RazorEvents->Branch("tau_muonVetoMVA", tau_muonVetoMVA, "tau_muonVetoMVA[nTaus]/F");
  RazorEvents->Branch("tau_isoMVAnewDMwLT", tau_isoMVAnewDMwLT, "tau_isoMVAnewDMwLT[nTaus]/F");
  RazorEvents->Branch("tau_isoMVAnewDMwoLT", tau_isoMVAnewDMwoLT, "tau_isoMVAnewDMwoLT[nTaus]/F");
  RazorEvents->Branch("tau_leadCandPt", tau_leadCandPt, "tau_leadCandPt[nTaus]/F");
  RazorEvents->Branch("tau_leadCandID", tau_leadCandID, "tau_leadCandID[nTaus]/I");
  RazorEvents->Branch("tau_leadChargedHadrCandPt", tau_leadChargedHadrCandPt, "tau_leadChargedHadrCandPt[nTaus]/F");
  RazorEvents->Branch("tau_leadChargedHadrCandID", tau_leadChargedHadrCandID, "tau_leadChargedHadrCandID[nTaus]/I"); 
};

void RazorAna::enablePhotonBranches(){
  RazorTuplizer::enablePhotonBranches();
  RazorEvents->Branch("phoSigmaIetaIeta", phoSigmaIetaIeta, "phoSigmaIetaIeta[nPhotons]/F");
  RazorEvents->Branch("phoFull5x5SigmaIetaIeta", phoFull5x5SigmaIetaIeta, "phoFull5x5SigmaIetaIeta[nPhotons]/F");
  RazorEvents->Branch("phoR9", phoR9, "phoR9[nPhotons]/F");
  RazorEvents->Branch("pho_HoverE", pho_HoverE, "pho_HoverE[nPhotons]/F");
  RazorEvents->Branch("pho_sumChargedHadronPt", pho_sumChargedHadronPt, "pho_sumChargedHadronPt[nPhotons]/F");
  RazorEvents->Branch("pho_sumNeutralHadronEt", pho_sumNeutralHadronEt, "pho_sumNeutralHadronEt[nPhotons]/F");
  RazorEvents->Branch("pho_sumPhotonEt", pho_sumPhotonEt, "pho_sumPhotonEt[nPhotons]/F");
  RazorEvents->Branch("pho_isConversion", pho_isConversion, "pho_isConversion[nPhotons]/I");
  RazorEvents->Branch("pho_RegressionE", pho_RegressionE, "pho_RegressionE[nPhotons]/F");
  RazorEvents->Branch("pho_IDMVA", pho_IDMVA, "pho_IDMVA[nPhotons]/F");
};

void RazorAna::enableJetBranches(){
  RazorTuplizer::enableJetBranches();
  RazorEvents->Branch("jetMass", jetMass, "jetMass[nJets]/F");
  RazorEvents->Branch("jetJetArea", jetJetArea, "jetJetArea[nJets]/F");
  RazorEvents->Branch("jetPileupE", jetPileupE, "jetPileupE[nJets]/F");
  RazorEvents->Branch("jetPileupId", jetPileupId, "jetPileupId[nJets]/F");
};

void RazorAna::enableJetAK8Branches(){
  RazorTuplizer::enableJetAK8Branches();
};

void RazorAna::enableMetBranches(){
  RazorTuplizer::enableMetBranches();
  RazorEvents->Branch("sumMET", &sumMET, "sumMET/F");
}

void RazorAna::enableRazorBranches(){
  RazorTuplizer::enableRazorBranches();
}
void RazorAna::enableGenParticles(){
  RazorEvents->Branch("nGenParticle", &nGenParticle, "nGenParticle/s");
  RazorEvents->Branch("gParticleMotherId", gParticleMotherId, "gParticleMotherId[nGenParticle]/I");
  RazorEvents->Branch("gParticleMotherIndex", gParticleMotherIndex, "gParticleMotherIndex[nGenParticle]/I");
  RazorEvents->Branch("gParticleId", gParticleId, "gParticleId[nGenParticle]/I");
  RazorEvents->Branch("gParticleStatus", gParticleStatus, "gParticleStatus[nGenParticle]/I");
  RazorEvents->Branch("gParticleE", gParticleE, "gParticleE[nGenParticle]/F");
  RazorEvents->Branch("gParticlePt", gParticlePt, "gParticlePt[nGenParticle]/F");
  RazorEvents->Branch("gParticleEta", gParticleEta, "gParticleEta[nGenParticle]/F");
  RazorEvents->Branch("gParticlePhi", gParticlePhi, "gParticlePhi[nGenParticle]/F");
}

/*
Re-defining Fill methods (will not use mother class methods)
*/

bool RazorAna::fillEventInfo(const edm::Event& iEvent){
  //store basic event info
  runNum = iEvent.id().run();
  lumiNum = iEvent.luminosityBlock();
  eventNum = iEvent.id().event();

  //select the primary vertex, if any
  if (vertices->empty()) return false; // skip the event if no PV found
  const reco::Vertex &PV = vertices->front();
  pvX = PV.x();
  pvY = PV.y();
  pvZ= PV.z();
  
  nPV = 0;
  //Check for good vertices
  for(unsigned int i = 0; i < vertices->size(); i++){
    if(vertices->at(i).isValid() && !vertices->at(i).isFake())nPV++;
  }
  if(nPV == 0)return false;

  return true;  
};

bool RazorAna::fillPileUp(){
  for(const PileupSummaryInfo &pu : *puInfo){
    BunchXing[nBunchXing] = pu.getBunchCrossing();
    nPU[nBunchXing] = pu.getPU_NumInteractions();
    nPUmean[nBunchXing] = pu.getTrueNumInteractions();
    nBunchXing++;
    //std::cout << "BC: " << pu.getBunchCrossing() << std::endl;
  }
  return true;
};

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
    muonType[nMuons] = mu.isMuon() + mu.isGlobalMuon() + mu.isTrackerMuon() + mu.isStandAloneMuon()
      + mu.isCaloMuon() + mu.isPFMuon() + mu.isRPCMuon();
    muonQuality[nMuons] = 
      muon::isGoodMuon(mu,muon::All)
    + muon::isGoodMuon(mu,muon::AllGlobalMuons)
    + muon::isGoodMuon(mu,muon::AllStandAloneMuons)
    + muon::isGoodMuon(mu,muon::AllTrackerMuons)
    + muon::isGoodMuon(mu,muon::TrackerMuonArbitrated)
    + muon::isGoodMuon(mu,muon::AllArbitrated)      
    + muon::isGoodMuon(mu,muon::GlobalMuonPromptTight)      
    + muon::isGoodMuon(mu,muon::TMLastStationLoose)      
    + muon::isGoodMuon(mu,muon::TMLastStationTight)      
    + muon::isGoodMuon(mu,muon::TM2DCompatibilityLoose)      
    + muon::isGoodMuon(mu,muon::TM2DCompatibilityTight)      
    + muon::isGoodMuon(mu,muon::TMOneStationLoose)      
    + muon::isGoodMuon(mu,muon::TMOneStationTight)      
    + muon::isGoodMuon(mu,muon::TMLastStationOptimizedLowPtLoose)      
    + muon::isGoodMuon(mu,muon::TMLastStationOptimizedLowPtTight)      
    + muon::isGoodMuon(mu,muon::GMTkChiCompatibility)      
    + muon::isGoodMuon(mu,muon::GMStaChiCompatibility)      
    + muon::isGoodMuon(mu,muon::GMTkKinkTight)      
    + muon::isGoodMuon(mu,muon::TMLastStationAngLoose)      
    + muon::isGoodMuon(mu,muon::TMLastStationAngTight)      
    + muon::isGoodMuon(mu,muon::TMOneStationAngLoose)      
    + muon::isGoodMuon(mu,muon::TMOneStationAngTight)      
    + muon::isGoodMuon(mu,muon::TMLastStationOptimizedBarrelLowPtLoose)      
    + muon::isGoodMuon(mu,muon::TMLastStationOptimizedBarrelLowPtTight)
    + muon::isGoodMuon(mu,muon::RPCMuLoose);       
    muon_relIso04DBetaCorr[nMuons] = ( mu.pfIsolationR04().sumChargedHadronPt + fmax(0, mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt - 0.5* mu.pfIsolationR04().sumPUPt) ) / mu.pt();
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
    eleCharge[nElectrons] = ele.charge();
    eleE_SC[nElectrons] = ele.superCluster()->energy();
    eleEta_SC[nElectrons] = ele.superCluster()->eta();
    elePhi_SC[nElectrons] = ele.superCluster()->phi();
    eleSigmaIetaIeta[nElectrons] = ele.sigmaIetaIeta();
    eleFull5x5SigmaIetaIeta[nElectrons] = ele.full5x5_sigmaIetaIeta();
    eleR9[nElectrons] = ele.r9();
    ele_dEta[nElectrons] = ele.deltaEtaSuperClusterTrackAtVtx();
    ele_dPhi[nElectrons] = ele.deltaPhiSuperClusterTrackAtVtx();
    ele_HoverE[nElectrons] = ele.hcalOverEcal();
    ele_d0[nElectrons] = -ele.gsfTrack().get()->dxy(PV.position());
    ele_dZ[nElectrons] = ele.gsfTrack().get()->dz(PV.position());
    ele_relIsoDBetaCorr[nElectrons] = ( ele.pfIsolationVariables().sumChargedHadronPt + fmax(0,ele.pfIsolationVariables().sumNeutralHadronEt +  ele.pfIsolationVariables().sumPhotonEt - 0.5*ele.pfIsolationVariables().sumPUPt) ) / ele.pt();
    ele_MissHits[nElectrons] = ele.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits();

    //Conversion Veto
    ele_PassConvVeto[nElectrons] = false;
    if( beamSpot.isValid() && conversions.isValid() ) {
      ele_PassConvVeto[nElectrons] = !ConversionTools::hasMatchedConversion(ele,conversions,
									    beamSpot->position());
    } else {
      cout << "\n\nERROR!!! conversions not found!!!\n";
    }
  
    // 1/E - 1/P
    if( ele.ecalEnergy() == 0 ){
      ele_OneOverEminusOneOverP[nElectrons] = 1e30;
    } else if( !std::isfinite(ele.ecalEnergy())){
      ele_OneOverEminusOneOverP[nElectrons] = 1e30;
    } else {
    ele_OneOverEminusOneOverP[nElectrons] = 1./ele.ecalEnergy()  -  ele.eSuperClusterOverP()/ele.ecalEnergy();    
    }

    //ID MVA
    ele_IDMVATrig[nElectrons] = myMVATrig->mvaValue(ele,false);
    ele_IDMVANonTrig[nElectrons] = myMVANonTrig->mvaValue(ele,false);

    ele_RegressionE[nElectrons] = ele.ecalRegressionEnergy();
    ele_CombineP4[nElectrons] = ele.ecalTrackRegressionEnergy();
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
    
    tau_IsLoose[nTaus] = bool(tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"));
    tau_IsMedium[nTaus] = bool(tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"));
    tau_IsTight[nTaus] = bool(tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"));
    tau_passEleVetoLoose[nTaus] = bool(tau.tauID("againstElectronLooseMVA5"));
    tau_passEleVetoMedium[nTaus] = bool(tau.tauID("againstElectronMediumMVA5"));
    tau_passEleVetoTight[nTaus] = bool(tau.tauID("againstElectronTightMVA5"));
    tau_passMuVetoLoose[nTaus] = bool(tau.tauID("againstMuonLooseMVA"));
    tau_passMuVetoMedium[nTaus] = bool(tau.tauID("againstMuonMediumMVA"));
    tau_passMuVetoTight[nTaus] = bool(tau.tauID("againstMuonTightMVA") );  
  
    tau_combinedIsoDeltaBetaCorr3Hits[nTaus] = tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
    tau_eleVetoMVA[nTaus] = tau.tauID("againstElectronMVA5raw") ;
    tau_eleVetoCategory[nTaus] = tau.tauID("againstElectronMVA5category");
    tau_muonVetoMVA[nTaus] = tau.tauID("againstMuonMVAraw");
    tau_isoMVAnewDMwLT[nTaus] = tau.tauID("byIsolationMVA3newDMwLTraw");
    tau_isoMVAnewDMwoLT[nTaus] = tau.tauID("byIsolationMVA3newDMwoLTraw") ; 

    tau_ID[nTaus] = 
      bool(tau.tauID("decayModeFinding")) +
      bool(tau.tauID("decayModeFindingNewDMs")) +
      bool(tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")) +
      bool(tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits")) +
      bool(tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits")) +
      bool(tau.tauID("againstElectronVLooseMVA5")) +
      bool(tau.tauID("againstElectronLooseMVA5")) +
      bool(tau.tauID("againstElectronMediumMVA5")) +
      bool(tau.tauID("againstElectronTightMVA5")) +
      bool(tau.tauID("againstElectronVTightMVA5")) +
      bool(tau.tauID("againstMuonLoose3")) +
      bool(tau.tauID("againstMuonTight3")) +
      bool(tau.tauID("againstMuonLooseMVA")) +
      bool(tau.tauID("againstMuonMediumMVA")) +
      bool(tau.tauID("againstMuonTightMVA")) +
      bool(tau.tauID("byVLooseIsolationMVA3newDMwLT")) +
      bool(tau.tauID("byLooseIsolationMVA3newDMwLT")) +
      bool(tau.tauID("byMediumIsolationMVA3newDMwLT")) +
      bool(tau.tauID("byTightIsolationMVA3newDMwLT")) +
      bool(tau.tauID("byVTightIsolationMVA3newDMwLT")) +
      bool(tau.tauID("byVVTightIsolationMVA3newDMwLT"));

    tau_leadCandPt[nTaus] = 0;
    tau_leadCandID[nTaus] = 0;
    tau_leadChargedHadrCandPt[nTaus] = 0;
    tau_leadChargedHadrCandID[nTaus] = 0;
    if (tau.leadCand().isNonnull()) {
      tau_leadCandPt[nTaus] = tau.leadCand()->pt();
      tau_leadCandID[nTaus] = tau.leadCand()->pdgId();
    }
    if (tau.leadChargedHadrCand().isNonnull()) { 
      tau_leadChargedHadrCandPt[nTaus] = tau.leadChargedHadrCand()->pt();
      tau_leadChargedHadrCandID[nTaus] = tau.leadChargedHadrCand()->pdgId();
    }
      
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
    pho_HoverE[nPhotons] = pho.hadronicOverEm();
    pho_sumChargedHadronPt[nPhotons] = pho.chargedHadronIso();
    pho_sumNeutralHadronEt[nPhotons] = pho.neutralHadronIso();
    pho_sumPhotonEt[nPhotons] = pho.photonIso();
    pho_isConversion[nPhotons] = pho.hasConversionTracks();
    pho_RegressionE[nPhotons] = pho.getCorrectedEnergy(reco::Photon::P4type::regression1);
    pho_IDMVA[nPhotons] = pho.pfMVA();
    /*
    const reco::Candidate* genPhoton = pho.genPhoton();
    if(genPhoton != NULL)std::cout << "======>gen PT: " << genPhoton->pt() <<
      " recoPT: " << pho.pt() << std::endl;
    */
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
    jetMass[nJets] = j.mass();
    jetJetArea[nJets] = j.jetArea();
    jetPileupE[nJets] = j.pileup();
    jetPileupId[nJets] = j.userFloat("pileupJetId:fullDiscriminant");
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
  sumMET = Met.sumEt();
  genMetPt = Met.genMET()->pt();
  genMetPhi = Met.genMET()->phi();
  return true;
};

bool RazorAna::fillRazor(){
  //get good jets for razor calculation                                  
  vector<TLorentzVector> goodJets;
  for (const pat::Jet &j : *jets) {
    if (j.pt() < 40) continue;
    if (fabs(j.eta()) > 3.0) continue;
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

bool RazorAna::fillGenParticles(){
  std::vector<const reco::Candidate*> prunedV;//Allows easier comparison for mother finding
  //Fills selected gen particles
  for(size_t i=0; i<prunedGenParticles->size();i++){
    if(
       (abs((*prunedGenParticles)[i].pdgId()) >= 1 && abs((*prunedGenParticles)[i].pdgId()) <= 6 
    	&& ( (*prunedGenParticles)[i].status() < 30 	     
	     )
	)
       || (abs((*prunedGenParticles)[i].pdgId()) >= 11 && abs((*prunedGenParticles)[i].pdgId()) <= 16)
       || (abs((*prunedGenParticles)[i].pdgId()) == 21 
    	   && (*prunedGenParticles)[i].status() < 30
    	   )
       || (abs((*prunedGenParticles)[i].pdgId()) >= 22 && abs((*prunedGenParticles)[i].pdgId()) <= 25
    	   && ( (*prunedGenParticles)[i].status() < 30
		)
	   )
       || (abs((*prunedGenParticles)[i].pdgId()) >= 32 && abs((*prunedGenParticles)[i].pdgId()) <= 42)
       || (abs((*prunedGenParticles)[i].pdgId()) >= 1000001 && abs((*prunedGenParticles)[i].pdgId()) <= 1000039)
       ){
      prunedV.push_back(&(*prunedGenParticles)[i]);
    }
    
    //cout << i << " : " << (*prunedGenParticles)[i].pdgId() << " " << (*prunedGenParticles)[i].status() << " " << (*prunedGenParticles)[i].pt() << "\n";
    //if (prunedV.size()<99) prunedV.push_back(&(*prunedGenParticles)[i]); //keep all pruned particles
  }

  //Total number of gen particles
  nGenParticle = prunedV.size();
  //Look for mother particle and Fill gen variables
  for(unsigned int i = 0; i < prunedV.size(); i++){
    gParticleId[i] = prunedV[i]->pdgId();
    gParticleStatus[i] = prunedV[i]->status();
    gParticleE[i] = prunedV[i]->energy();
    gParticlePt[i] = prunedV[i]->pt();
    gParticleEta[i] = prunedV[i]->eta();
    gParticlePhi[i] = prunedV[i]->phi();
    gParticleMotherId[i] = 0;
    gParticleMotherIndex[i] = -1;
    if(prunedV[i]->numberOfMothers() > 0){
      
      //find the ID of the first mother that has a different ID than the particle itself
      const reco::Candidate* firstMotherWithDifferentID = findFirstMotherWithDifferentID(prunedV[i]);
      if (firstMotherWithDifferentID) {
	gParticleMotherId[i] = firstMotherWithDifferentID->pdgId();
      }

      //find the mother and keep going up the mother chain if the ID's are the same
      const reco::Candidate* originalMotherWithSameID = findOriginalMotherWithSameID(prunedV[i]);
      for(unsigned int j = 0; j < prunedV.size(); j++){	
	if(prunedV[j] == originalMotherWithSameID){
	  gParticleMotherIndex[i] = j;
	  break;
	}
      }
    } else {
      gParticleMotherIndex[i] = -1;
    }
  }
  return true;
};


const reco::Candidate* RazorAna::findFirstMotherWithDifferentID(const reco::Candidate *particle){

  if( particle == 0 ){
    printf("ERROR! null candidate pointer, this should never happen\n");
    return 0;
  }

  // Is this the first parent with a different ID? If yes, return, otherwise
  // go deeper into recursion
  if (particle->numberOfMothers() > 0 && particle->pdgId() != 0) {
    if (particle->pdgId() == particle->mother(0)->pdgId()) {
      return findFirstMotherWithDifferentID(particle->mother(0));
    } else {
      return particle->mother(0);
    }
  }

  return 0;
}

const reco::Candidate* RazorAna::findOriginalMotherWithSameID(const reco::Candidate *particle){

  if( particle == 0 ){
    printf("ERROR! null candidate pointer, this should never happen\n");
    return 0;
  }

  // Is there another parent with the same ID? If yes, go deeper into recursion
  if (particle->numberOfMothers() > 0 && particle->pdgId() != 0) {
    if (particle->mother(0)->numberOfMothers() == 0 || 
	(particle->mother(0)->numberOfMothers() > 0 && particle->mother(0)->mother(0)->pdgId() != particle->mother(0)->pdgId())
	) {
      return particle->mother(0);
    } else {
      return findOriginalMotherWithSameID(particle->mother(0));
    }
  }
  return 0;
}


bool RazorAna::isAncestor(const reco::Candidate* ancestor, const reco::Candidate* particle){
  //particle is already the ancestor
  if(ancestor == particle ) return true;
  //otherwise loop on mothers, if any and return true if the ancestor is found
  for(size_t i=0;i< particle->numberOfMothers();i++){
    if(isAncestor(ancestor,particle->mother(i))) return true;
  }
  //if we did not return yet, then particle and ancestor are not relatives
  return false;
};

//------ Method called once each job just before starting event loop ------//
void RazorAna::beginJob(){
  setBranches();
};

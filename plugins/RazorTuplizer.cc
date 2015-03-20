// -*- C++ -*-
// Class:      RazorTuplizer
/*
  Description: Base class for miniAOD analysis with CRAB
*/
//         Author:  Caltech razor team
//         Created:  Thu, 17 Jul 2014 15:00:06 GMT

#include "RazorTuplizer.h"
//------ Constructors and destructor ------//
RazorTuplizer::RazorTuplizer(const edm::ParameterSet& iConfig): 
  //get inputs from config file
  enableTriggerInfo_(iConfig.getParameter<bool> ("enableTriggerInfo")),
  triggerPathNamesFile_(iConfig.getParameter<string> ("triggerPathNamesFile")),
  verticesToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  muonsToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronsToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  tausToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  photonsToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
  jetsToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  jetsAK8Token_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetsAK8"))),
  packedPFCandsToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packedPfCands"))),
  prunedGenParticlesToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("prunedGenParticles"))),
  packedGenParticlesToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packedGenParticles"))),
  genJetsToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets"))),
  triggerBitsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"))),
  triggerObjectsToken_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
  triggerPrescalesToken_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("triggerPrescales"))),     
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
  metFilterBitsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilterBits"))),
  lheInfoToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheInfo"))),
  genInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genInfo"))),
  puInfoToken_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("puInfo"))),
  hcalNoiseInfoToken_(consumes<HcalNoiseSummary>(iConfig.getParameter<edm::InputTag>("hcalNoiseInfo"))),
  secondaryVerticesToken_(consumes<vector<reco::VertexCompositePtrCandidate> >(iConfig.getParameter<edm::InputTag>("secondaryVertices"))),
  rhoAllToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoAll"))),
  rhoFastjetAllToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetAll"))),
  rhoFastjetAllCaloToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetAllCalo"))),
  rhoFastjetCentralCaloToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetCentralCalo"))),
  rhoFastjetCentralChargedPileUpToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetCentralChargedPileUp"))),
  rhoFastjetCentralNeutralToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetCentralNeutral"))),
  beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
  ebRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("ebRecHits"))),
  eeRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("eeRecHits"))),
  esRecHitsToken_(consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > >(iConfig.getParameter<edm::InputTag>("esRecHits"))),
  ebeeClustersToken_(consumes<vector<reco::CaloCluster> >(iConfig.getParameter<edm::InputTag>("ebeeClusters"))),
  esClustersToken_(consumes<vector<reco::CaloCluster> >(iConfig.getParameter<edm::InputTag>("esClusters"))),
  conversionsToken_(consumes<vector<reco::Conversion> >(iConfig.getParameter<edm::InputTag>("conversions"))),
  singleLegConversionsToken_(consumes<vector<reco::Conversion> >(iConfig.getParameter<edm::InputTag>("singleLegConversions"))),
  gedGsfElectronCoresToken_(consumes<vector<reco::GsfElectronCore> >(iConfig.getParameter<edm::InputTag>("gedGsfElectronCores"))),
  gedPhotonCoresToken_(consumes<vector<reco::PhotonCore> >(iConfig.getParameter<edm::InputTag>("gedPhotonCores"))),
  superClustersToken_(consumes<vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("superClusters"))),
  lostTracksToken_(consumes<vector<pat::PackedCandidate> >(iConfig.getParameter<edm::InputTag>("lostTracks")))
{
  //declare the TFileService for output
  edm::Service<TFileService> fs;
  
  //set up output tree
  RazorEvents = fs->make<TTree>("RazorEvents", "selected miniAOD information");
  NEvents = fs->make<TH1F>("NEvents",";;NEvents;",1,-0.5,0.5);

  //set up electron MVA ID
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

  //set up photon MVA ID
  std::vector<std::string> myPhotonMVAWeights;
  myPhotonMVAWeights.push_back(edm::FileInPath("SUSYBSMAnalysis/RazorTuplizer/data/PhotonIDMVA_PHYS14_EB.weights.xml").fullPath().c_str());
  myPhotonMVAWeights.push_back(edm::FileInPath("SUSYBSMAnalysis/RazorTuplizer/data/PhotonIDMVA_PHYS14_EE.weights.xml").fullPath().c_str());
  std::vector<std::string> myPhotonMVAMethodNames;
  myPhotonMVAMethodNames.push_back("BDTG photons barrel");
  myPhotonMVAMethodNames.push_back("BDTG photons endcap");

  myPhotonMVA = new EGammaMvaPhotonEstimator();
  myPhotonMVA->initialize(myPhotonMVAMethodNames,myPhotonMVAWeights,
			  EGammaMvaPhotonEstimator::kPhotonMVATypeDefault);



  //Read in HLT Trigger Path List from config file
  for (int i = 0; i<NTriggersMAX; ++i) triggerPathNames[i] = "";
  ifstream myfile (edm::FileInPath(triggerPathNamesFile_.c_str()).fullPath().c_str()) ;
  if (myfile.is_open()) {
    string line;
    int index;
    string hltpathname;

    while(myfile>>index>>hltpathname) {
      // loop only eneterd if both c_tmp AND gamma_tmp
      // can be retrieved from the file.
      if (index < NTriggersMAX) {
	triggerPathNames[index] = hltpathname;
      }    
    }    
    myfile.close();
  } else {
    cout << "ERROR!!! Could not open trigger path name file : " << edm::FileInPath(triggerPathNamesFile_.c_str()).fullPath().c_str() << "\n";
  }
  
  if(enableTriggerInfo_) {
    cout << "\n";
    cout << "****************** Trigger Paths Defined For Razor Ntuple ******************\n";    
    for (int i = 0; i<NTriggersMAX; ++i) {
      if (triggerPathNames[i] != "") cout << "Trigger " << i << " " << triggerPathNames[i] << "\n";
    }
    cout << "****************************************************************************\n";    
    cout << "\n";
  }

}

RazorTuplizer::~RazorTuplizer()
{
}

//------ Enable the desired set of branches ------//
void RazorTuplizer::setBranches(){
  enableEventInfoBranches();
  enablePileUpBranches();
  enableMuonBranches();
  enableElectronBranches();
  enableTauBranches();
  enableIsoPFCandidateBranches();
  enablePhotonBranches();
  enableJetBranches();
  enableJetAK8Branches();
  enableMetBranches();
  enableRazorBranches();
  if (enableTriggerInfo_) enableTriggerBranches();
  enableMCBranches();
  enableGenParticleBranches();
}

void RazorTuplizer::enableEventInfoBranches(){
  RazorEvents->Branch("nPV", &nPV, "nPV/I");
  RazorEvents->Branch("runNum", &runNum, "runNum/I");
  RazorEvents->Branch("lumiNum", &lumiNum, "lumiNum/I");
  RazorEvents->Branch("eventNum", &eventNum, "eventNum/I");
  RazorEvents->Branch("pvX", &pvX, "pvX/F");
  RazorEvents->Branch("pvY", &pvY, "pvY/F");
  RazorEvents->Branch("pvZ", &pvZ, "pvZ/F");
  RazorEvents->Branch("fixedGridRhoAll", &fixedGridRhoAll, "fixedGridRhoAll/F");
  RazorEvents->Branch("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll, "fixedGridRhoFastjetAll/F");
  RazorEvents->Branch("fixedGridRhoFastjetAllCalo", &fixedGridRhoFastjetAllCalo, "fixedGridRhoFastjetAllCalo/F");
  RazorEvents->Branch("fixedGridRhoFastjetCentralCalo", &fixedGridRhoFastjetCentralCalo, "fixedGridRhoFastjetCentralCalo/F");
  RazorEvents->Branch("fixedGridRhoFastjetCentralChargedPileUp", &fixedGridRhoFastjetCentralChargedPileUp, "fixedGridRhoFastjetCentralChargedPileUp/F");
  RazorEvents->Branch("fixedGridRhoFastjetCentralNeutral", &fixedGridRhoFastjetCentralNeutral, "fixedGridRhoFastjetCentralNeutral/F");
}

void RazorTuplizer::enablePileUpBranches(){
  RazorEvents->Branch("nBunchXing", &nBunchXing, "nBunchXing/I");
  RazorEvents->Branch("BunchXing", BunchXing, "BunchXing[nBunchXing]/I");
  RazorEvents->Branch("nPU", nPU, "nPU[nBunchXing]/I");
  RazorEvents->Branch("nPUmean", nPUmean, "nPUmean[nBunchXing]/F");
};

void RazorTuplizer::enableMuonBranches(){
  RazorEvents->Branch("nMuons", &nMuons,"nMuons/I");
  RazorEvents->Branch("muonE", muonE,"muonE[nMuons]/F");
  RazorEvents->Branch("muonPt", muonPt,"muonPt[nMuons]/F");
  RazorEvents->Branch("muonEta", muonEta,"muonEta[nMuons]/F");
  RazorEvents->Branch("muonPhi", muonPhi,"muonPhi[nMuons]/F");
  RazorEvents->Branch("muonCharge", muonCharge, "muonCharge[nMuons]/I");
  RazorEvents->Branch("muonIsLoose", muonIsLoose,"muonIsLoose[nMuons]/O");
  RazorEvents->Branch("muonIsTight", muonIsTight,"muonIsTight[nMuons]/O");
  RazorEvents->Branch("muon_d0", muon_d0, "muon_d0[nMuons]/F");
  RazorEvents->Branch("muon_dZ", muon_dZ, "muon_dZ[nMuons]/F");
  RazorEvents->Branch("muon_ip3d", muon_ip3d, "muon_ip3d[nMuons]/F");
  RazorEvents->Branch("muon_ip3dSignificance", muon_ip3dSignificance, "muon_ip3dSignificance[nMuons]/F");
  RazorEvents->Branch("muonType", muonType, "muonType[nMuons]/i");
  RazorEvents->Branch("muonQuality", muonQuality, "muonQuality[nMuons]/i");
  RazorEvents->Branch("muon_pileupIso", muon_pileupIso, "muon_pileupIso[nMuons]/F");
  RazorEvents->Branch("muon_chargedIso", muon_chargedIso, "muon_chargedIso[nMuons]/F");
  RazorEvents->Branch("muon_photonIso", muon_photonIso, "muon_photonIso[nMuons]/F");
  RazorEvents->Branch("muon_neutralHadIso", muon_neutralHadIso, "muon_neutralHadIso[nMuons]/F");
  RazorEvents->Branch("muon_ptrel", muon_ptrel, "muon_ptrel[nMuons]/F");
  RazorEvents->Branch("muon_miniiso", muon_miniiso, "muon_miniiso[nMuons]/F");
}

void RazorTuplizer::enableElectronBranches(){
  RazorEvents->Branch("nElectrons", &nElectrons,"nElectrons/I");
  RazorEvents->Branch("eleE", eleE,"eleE[nElectrons]/F");
  RazorEvents->Branch("elePt", elePt,"elePt[nElectrons]/F");
  RazorEvents->Branch("eleEta", eleEta,"eleEta[nElectrons]/F");
  RazorEvents->Branch("elePhi", elePhi,"elePhi[nElectrons]/F");
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
  RazorEvents->Branch("ele_pileupIso", ele_pileupIso, "ele_pileupIso[nElectrons]/F");
  RazorEvents->Branch("ele_chargedIso", ele_chargedIso, "ele_chargedIso[nElectrons]/F");
  RazorEvents->Branch("ele_photonIso", ele_photonIso, "ele_photonIso[nElectrons]/F");
  RazorEvents->Branch("ele_neutralHadIso", ele_neutralHadIso, "ele_neutralHadIso[nElectrons]/F");
  RazorEvents->Branch("ele_MissHits", ele_MissHits, "ele_MissHits[nElectrons]/I");
  RazorEvents->Branch("ele_PassConvVeto", ele_PassConvVeto, "ele_PassConvVeto[nElectrons]/O");
  RazorEvents->Branch("ele_OneOverEminusOneOverP", ele_OneOverEminusOneOverP, "ele_OneOverEminusOneOverP[nElectrons]/F");
  RazorEvents->Branch("ele_IDMVATrig", ele_IDMVATrig, "ele_IDMVATrig[nElectrons]/F");
  RazorEvents->Branch("ele_IDMVANonTrig", ele_IDMVANonTrig, "ele_IDMVANonTrig[nElectrons]/F");
  RazorEvents->Branch("ele_RegressionE", ele_RegressionE, "ele_RegressionE[nElectrons]/F");
  RazorEvents->Branch("ele_CombineP4", ele_CombineP4, "ele_CombineP4[nElectrons]/F");
  RazorEvents->Branch("ele_ptrel", ele_ptrel, "ele_ptrel[nElectrons]/F");
  RazorEvents->Branch("ele_miniiso", ele_miniiso, "ele_miniiso[nElectrons]/F");
}

void RazorTuplizer::enableTauBranches(){
  RazorEvents->Branch("nTaus", &nTaus,"nTaus/I");
  RazorEvents->Branch("tauE", tauE,"tauE[nTaus]/F");
  RazorEvents->Branch("tauPt", tauPt,"tauPt[nTaus]/F");
  RazorEvents->Branch("tauEta", tauEta,"tauEta[nTaus]/F");
  RazorEvents->Branch("tauPhi", tauPhi,"tauPhi[nTaus]/F");
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
}

void RazorTuplizer::enableIsoPFCandidateBranches(){
  RazorEvents->Branch("nIsoPFCandidates", &nIsoPFCandidates, "nIsoPFCandidates/i");
  RazorEvents->Branch("isoPFCandidatePt", isoPFCandidatePt, "isoPFCandidatePt[nIsoPFCandidates]/F");
  RazorEvents->Branch("isoPFCandidateEta", isoPFCandidateEta, "isoPFCandidateEta[nIsoPFCandidates]/F");
  RazorEvents->Branch("isoPFCandidatePhi", isoPFCandidatePhi, "isoPFCandidatePhi[nIsoPFCandidates]/F");
  RazorEvents->Branch("isoPFCandidateIso04", isoPFCandidateIso04, "isoPFCandidateIso04[nIsoPFCandidates]/F");
  RazorEvents->Branch("isoPFCandidateD0", isoPFCandidateD0, "isoPFCandidateD0[nIsoPFCandidates]/F");
  RazorEvents->Branch("isoPFCandidatePdgId", isoPFCandidatePdgId, "isoPFCandidatePdgId[nIsoPFCandidates]/I");  
}

void RazorTuplizer::enablePhotonBranches(){
  RazorEvents->Branch("nPhotons", &nPhotons,"nPhotons/I");
  RazorEvents->Branch("phoE", phoE,"phoE[nPhotons]/F");
  RazorEvents->Branch("phoPt", phoPt,"phoPt[nPhotons]/F");
  RazorEvents->Branch("phoEta", phoEta,"phoEta[nPhotons]/F");
  RazorEvents->Branch("phoPhi", phoPhi,"phoPhi[nPhotons]/F");
  RazorEvents->Branch("phoSigmaIetaIeta", phoSigmaIetaIeta, "phoSigmaIetaIeta[nPhotons]/F");
  RazorEvents->Branch("phoFull5x5SigmaIetaIeta", phoFull5x5SigmaIetaIeta, "phoFull5x5SigmaIetaIeta[nPhotons]/F");
  RazorEvents->Branch("phoR9", phoR9, "phoR9[nPhotons]/F");
  RazorEvents->Branch("pho_HoverE", pho_HoverE, "pho_HoverE[nPhotons]/F");
  RazorEvents->Branch("pho_sumChargedHadronPt", pho_sumChargedHadronPt, "pho_sumChargedHadronPt[nPhotons]/F");
  RazorEvents->Branch("pho_sumNeutralHadronEt", pho_sumNeutralHadronEt, "pho_sumNeutralHadronEt[nPhotons]/F");
  RazorEvents->Branch("pho_sumPhotonEt", pho_sumPhotonEt, "pho_sumPhotonEt[nPhotons]/F");
  RazorEvents->Branch("pho_isConversion", pho_isConversion, "pho_isConversion[nPhotons]/O");
  RazorEvents->Branch("pho_passEleVeto", pho_passEleVeto, "pho_passEleVeto[nPhotons]/O");
  RazorEvents->Branch("pho_RegressionE", pho_RegressionE, "pho_RegressionE[nPhotons]/F");
  RazorEvents->Branch("pho_RegressionEUncertainty", pho_RegressionEUncertainty, "pho_RegressionEUncertainty[nPhotons]/F");
  RazorEvents->Branch("pho_IDMVA", pho_IDMVA, "pho_IDMVA[nPhotons]/F");
  RazorEvents->Branch("pho_superClusterEta", pho_superClusterEta, "pho_superClusterEta[nPhotons]/F");
  RazorEvents->Branch("pho_superClusterPhi", pho_superClusterPhi, "pho_superClusterPhi[nPhotons]/F");
  RazorEvents->Branch("pho_hasPixelSeed", pho_hasPixelSeed, "pho_hasPixelSeed[nPhotons]/O");
}

void RazorTuplizer::enableJetBranches(){
  RazorEvents->Branch("nJets", &nJets,"nJets/I");
  RazorEvents->Branch("jetE", jetE,"jetE[nJets]/F");
  RazorEvents->Branch("jetPt", jetPt,"jetPt[nJets]/F");
  RazorEvents->Branch("jetEta", jetEta,"jetEta[nJets]/F");
  RazorEvents->Branch("jetPhi", jetPhi,"jetPhi[nJets]/F");
  RazorEvents->Branch("jetCSV", jetCSV,"jetCSV[nJets]/F");
  RazorEvents->Branch("jetCISV", jetCISV,"jetCISV[nJets]/F");
  RazorEvents->Branch("jetMass", jetMass, "jetMass[nJets]/F");
  RazorEvents->Branch("jetJetArea", jetJetArea, "jetJetArea[nJets]/F");
  RazorEvents->Branch("jetPileupE", jetPileupE, "jetPileupE[nJets]/F");
  RazorEvents->Branch("jetPileupId", jetPileupId, "jetPileupId[nJets]/F");
  RazorEvents->Branch("jetPartonFlavor", jetPartonFlavor, "jetPartonFlavor[nJets]/I");
  RazorEvents->Branch("jetHadronFlavor", jetHadronFlavor, "jetHadronFlavor[nJets]/I");
}

void RazorTuplizer::enableJetAK8Branches(){
  RazorEvents->Branch("nFatJets", &nFatJets,"nFatJets/i");
  RazorEvents->Branch("fatJetE", fatJetE,"fatJetE[nFatJets]/F");
  RazorEvents->Branch("fatJetPt", fatJetPt,"fatJetPt[nFatJets]/F");
  RazorEvents->Branch("fatJetEta", fatJetEta,"fatJetEta[nFatJets]/F");
  RazorEvents->Branch("fatJetPhi", fatJetPhi,"fatJetPhi[nFatJets]/F");
  RazorEvents->Branch("fatJetPrunedM", fatJetPrunedM,"fatJetPrunedM[nFatJets]/F");
  RazorEvents->Branch("fatJetTrimmedM", fatJetTrimmedM,"fatJetTrimmedM[nFatJets]/F");
  RazorEvents->Branch("fatJetFilteredM", fatJetFilteredM,"fatJetFilteredM[nFatJets]/F");
  RazorEvents->Branch("fatJetTau1", fatJetTau1,"fatJetTau1[nFatJets]/F");
  RazorEvents->Branch("fatJetTau2", fatJetTau2,"fatJetTau2[nFatJets]/F");
  RazorEvents->Branch("fatJetTau3", fatJetTau3,"fatJetTau3[nFatJets]/F");
}

void RazorTuplizer::enableMetBranches(){
  RazorEvents->Branch("metPt", &metPt, "metPt/F");
  RazorEvents->Branch("metPhi", &metPhi, "metPhi/F");
  RazorEvents->Branch("sumMET", &sumMET, "sumMET/F");
  RazorEvents->Branch("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, "Flag_HBHENoiseFilter/O");
  RazorEvents->Branch("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, "Flag_CSCTightHaloFilter/O");
  RazorEvents->Branch("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, "Flag_hcalLaserEventFilter/O");
  RazorEvents->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, "Flag_EcalDeadCellTriggerPrimitiveFilter/O");
  RazorEvents->Branch("Flag_goodVertices", &Flag_goodVertices, "Flag_goodVertices/O");
  RazorEvents->Branch("Flag_trackingFailureFilter", &Flag_trackingFailureFilter, "Flag_trackingFailureFilter/O");
  RazorEvents->Branch("Flag_eeBadScFilter", &Flag_eeBadScFilter, "Flag_eeBadScFilter/O");
  RazorEvents->Branch("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, "Flag_ecalLaserCorrFilter/O");
  RazorEvents->Branch("Flag_trkPOGFilters", &Flag_trkPOGFilters, "Flag_trkPOGFilters/O");
  RazorEvents->Branch("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, "Flag_trkPOG_manystripclus53X/O");
  RazorEvents->Branch("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, "Flag_trkPOG_toomanystripclus53X/O");
  RazorEvents->Branch("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, "Flag_trkPOG_logErrorTooManyClusters/O");
  RazorEvents->Branch("Flag_METFilters", &Flag_METFilters, "Flag_METFilters/O");  
}

void RazorTuplizer::enableRazorBranches(){
  RazorEvents->Branch("MR", &MR, "MR/F");
  RazorEvents->Branch("RSQ", &RSQ, "RSQ/F");
  RazorEvents->Branch("MR_AK8", &MR_AK8, "MR_AK8/F");
  RazorEvents->Branch("RSQ_AK8", &RSQ_AK8, "RSQ_AK8/F");
}

void RazorTuplizer::enableTriggerBranches(){
  nameHLT = new std::vector<std::string>; nameHLT->clear();
  //RazorEvents->Branch("nameHLT", "std::vector<std::string>",&nameHLT);
  RazorEvents->Branch("HLTDecision", &triggerDecision, ("HLTDecision[" + std::to_string(NTriggersMAX) +  "]/O").c_str());
}

void RazorTuplizer::enableMCBranches(){
  RazorEvents->Branch("nGenJets", &nGenJets, "nGenJets/I");
  RazorEvents->Branch("genJetE", genJetE, "genJetE[nGenJets]/F");
  RazorEvents->Branch("genJetPt", genJetPt, "genJetPt[nGenJets]/F");
  RazorEvents->Branch("genJetEta", genJetEta, "genJetEta[nGenJets]/F");
  RazorEvents->Branch("genJetPhi", genJetPhi, "genJetPhi[nGenJets]/F");
  RazorEvents->Branch("genMetPt", &genMetPt, "genMetPt/F");
  RazorEvents->Branch("genMetPhi", &genMetPhi, "genMetPhi/F");
  RazorEvents->Branch("genVertexZ", &genVertexZ, "genVertexZ/F");
  RazorEvents->Branch("genWeight", &genWeight, "genWeight/F");
  RazorEvents->Branch("genSignalProcessID", &genSignalProcessID, "genSignalProcessID/i");
  RazorEvents->Branch("genQScale", &genQScale, "genQScale/F");
  RazorEvents->Branch("genAlphaQCD", &genAlphaQCD, "genAlphaQCD/F");
  RazorEvents->Branch("genAlphaQED", &genAlphaQED, "genAlphaQED/F");
}

void RazorTuplizer::enableGenParticleBranches(){
  RazorEvents->Branch("nGenParticle", &nGenParticle, "nGenParticle/I");
  RazorEvents->Branch("gParticleMotherId", gParticleMotherId, "gParticleMotherId[nGenParticle]/I");
  RazorEvents->Branch("gParticleMotherIndex", gParticleMotherIndex, "gParticleMotherIndex[nGenParticle]/I");
  RazorEvents->Branch("gParticleId", gParticleId, "gParticleId[nGenParticle]/I");
  RazorEvents->Branch("gParticleStatus", gParticleStatus, "gParticleStatus[nGenParticle]/I");
  RazorEvents->Branch("gParticleE", gParticleE, "gParticleE[nGenParticle]/F");
  RazorEvents->Branch("gParticlePt", gParticlePt, "gParticlePt[nGenParticle]/F");
  RazorEvents->Branch("gParticleEta", gParticleEta, "gParticleEta[nGenParticle]/F");
  RazorEvents->Branch("gParticlePhi", gParticlePhi, "gParticlePhi[nGenParticle]/F");
}

//------ Load the miniAOD objects and reset tree variables for each event ------//
void RazorTuplizer::loadEvent(const edm::Event& iEvent){
  //load all miniAOD objects for the current event
  iEvent.getByToken(triggerBitsToken_, triggerBits);
  iEvent.getByToken(triggerObjectsToken_, triggerObjects);
  iEvent.getByToken(triggerPrescalesToken_, triggerPrescales);
  iEvent.getByToken(metFilterBitsToken_, metFilterBits);
  iEvent.getByToken(verticesToken_, vertices);
  iEvent.getByToken(packedPFCandsToken_, packedPFCands);
  iEvent.getByToken(muonsToken_, muons);
  iEvent.getByToken(electronsToken_, electrons);
  iEvent.getByToken(photonsToken_, photons);
  iEvent.getByToken(tausToken_, taus);
  iEvent.getByToken(jetsToken_, jets);
  iEvent.getByToken(jetsAK8Token_, jetsAK8);
  iEvent.getByToken(metToken_, mets);
  iEvent.getByToken(prunedGenParticlesToken_,prunedGenParticles);
  iEvent.getByToken(packedGenParticlesToken_,packedGenParticles);
  iEvent.getByToken(genJetsToken_,genJets);
  iEvent.getByToken(lheInfoToken_, lheInfo);
  iEvent.getByToken(genInfoToken_,genInfo);
  iEvent.getByToken(puInfoToken_,puInfo);
  iEvent.getByToken(hcalNoiseInfoToken_,hcalNoiseInfo);
  iEvent.getByToken(secondaryVerticesToken_,secondaryVertices);
  iEvent.getByToken(rhoAllToken_,rhoAll);
  iEvent.getByToken(rhoFastjetAllToken_,rhoFastjetAll);
  iEvent.getByToken(rhoFastjetAllCaloToken_,rhoFastjetAllCalo);
  iEvent.getByToken(rhoFastjetCentralCaloToken_,rhoFastjetCentralCalo);
  iEvent.getByToken(rhoFastjetCentralChargedPileUpToken_,rhoFastjetCentralChargedPileUp);
  iEvent.getByToken(rhoFastjetCentralNeutralToken_,rhoFastjetCentralNeutral);
  iEvent.getByToken(beamSpotToken_,beamSpot);
  iEvent.getByToken(ebRecHitsToken_,ebRecHits);
  iEvent.getByToken(eeRecHitsToken_,eeRecHits);
  iEvent.getByToken(esRecHitsToken_,esRecHits);
  iEvent.getByToken(ebeeClustersToken_,ebeeClusters);
  iEvent.getByToken(esClustersToken_,esClusters);
  iEvent.getByToken(conversionsToken_,conversions);
  iEvent.getByToken(singleLegConversionsToken_,singleLegConversions);
  iEvent.getByToken(gedGsfElectronCoresToken_,gedGsfElectronCores);
  iEvent.getByToken(gedPhotonCoresToken_, gedPhotonCores);
  iEvent.getByToken(superClustersToken_,superClusters);
  iEvent.getByToken(lostTracksToken_,lostTracks);
  
}   

//called by the loadEvent() method
void RazorTuplizer::resetBranches(){
    //reset tree variables
    nBunchXing = 0;
    nMuons = 0;
    nElectrons = 0;
    nTaus = 0;
    nPhotons = 0;
    nJets = 0;
    nFatJets = 0;
    nGenJets = 0;
    nGenParticle = 0;

    //nameHLT->clear();

    for(int i = 0; i < 99; i++){
        //PU
        BunchXing[i] = -99;
        nPU[i] = -99;
        nPUmean[i] = -99.0;

        //Muon
        muonE[i] = 0.0;
        muonPt[i] = 0.0;
        muonEta[i] = 0.0;
        muonPhi[i] = 0.0;
        muonCharge[i] = -99;
        muonIsLoose[i] = false;
        muonIsTight[i] = false;
        muon_d0[i] = -99.0;
        muon_dZ[i] = -99.0;
        muon_ip3d[i] = -99.0;
        muon_ip3dSignificance[i] = -99.0;
        muonType[i] = 0;
        muonQuality[i] = 0;
        muon_pileupIso[i] = -99.0;
        muon_chargedIso[i] = -99.0;
        muon_photonIso[i] = -99.0;
        muon_neutralHadIso[i] = -99.0;
	muon_ptrel[i] = -99.0;
	muon_miniiso[i] = -99.0;

        //Electron
        eleE[i] = 0.0;
        elePt[i] = 0.0;
        eleEta[i] = 0.0;
        elePhi[i] = 0.0;
        eleE_SC[i] = -99.0;
        eleEta_SC[i] = -99.0;
        elePhi_SC[i] = -99.0;
        eleSigmaIetaIeta[i] = -99.0;
        eleFull5x5SigmaIetaIeta[i] = -99.0;
        eleR9[i] = -99;
        ele_dEta[i] = -99;
        ele_dPhi[i] = -99;
        ele_HoverE[i] = -99;
        ele_d0[i] = -99;
        ele_dZ[i] = -99;
	ele_pileupIso[i] = -99.0;
        ele_chargedIso[i] = -99.0;
        ele_photonIso[i] = -99.0;
        ele_neutralHadIso[i] = -99.0;
	ele_MissHits[i] = -99;
        ele_PassConvVeto[i] = false;
        ele_OneOverEminusOneOverP[i] = -99.0;
        ele_IDMVATrig[i] = -99.0;
        ele_IDMVANonTrig[i] = -99.0;
        ele_RegressionE[i] = -99.0;
        ele_CombineP4[i] = -99.0;
	ele_ptrel[i] = -99.0;
	ele_miniiso[i] = -99.0;

        //Tau
        tauE[i] = 0.0;
        tauPt[i] = 0.0;
        tauEta[i] = 0.0;
        tauPhi[i] = 0.0;
        tau_IsLoose[i] = false;
        tau_IsMedium[i] = false;
        tau_IsTight[i] = false;
        tau_passEleVetoLoose[i] = false;
        tau_passEleVetoMedium[i] = false;
        tau_passEleVetoTight[i] = false;
        tau_passMuVetoLoose[i] = false;
        tau_passMuVetoMedium[i] = false;
        tau_passMuVetoTight[i] = false;
        tau_ID[i] = 0;
        tau_combinedIsoDeltaBetaCorr3Hits[i] = -99.0;
        tau_eleVetoMVA[i] = -99.0;
        tau_eleVetoCategory[i] = -1;
        tau_muonVetoMVA[i] = -99.0;
        tau_isoMVAnewDMwLT[i] = -99.0;
        tau_isoMVAnewDMwoLT[i] = -99.0;
        tau_leadCandPt[i] = -99.0;
        tau_leadCandID[i] = 0;
        tau_leadChargedHadrCandPt[i] = -99.0;
        tau_leadChargedHadrCandID[i] = 0;

        //IsoPFCandidates
        nIsoPFCandidates = 0;
        isoPFCandidatePt[i] = -99.0;
        isoPFCandidateEta[i] = -99.0;
        isoPFCandidatePhi[i] = -99.0;
        isoPFCandidateIso04[i] = -99.0;
        isoPFCandidateD0[i] = -99.0;
        isoPFCandidatePdgId[i] = 0;

        //Photon
        phoE[i] = 0.0;
        phoPt[i] = 0.0;
        phoEta[i] = 0.0;
        phoPhi[i] = 0.0;
        phoSigmaIetaIeta[i] = -99.0;
        phoFull5x5SigmaIetaIeta[i] = -99.0;
        phoR9[i] = -99.0;
        pho_HoverE[i] = -99.0;
        pho_sumChargedHadronPt[i] = -99.0;
        pho_sumNeutralHadronEt[i] = -99.0;
        pho_sumPhotonEt[i] = -99.0;
        pho_isConversion[i] = false;
        pho_passEleVeto[i] = false;    
        pho_RegressionE[i] = -99.0;
        pho_RegressionEUncertainty[i] = -99.0;
        pho_IDMVA[i] = -99.0;
        pho_superClusterEta[i] = -99.0;
        pho_superClusterPhi[i] = -99.0;
        pho_hasPixelSeed[i] = false;

        //Jet
        jetE[i] = 0.0;
        jetPt[i] = 0.0;
        jetEta[i] = 0.0;
        jetPhi[i] = 0.0;
        jetCSV[i] = 0.0;
        jetCISV[i] = 0.0;
        jetMass[i] =  -99.0;
        jetJetArea[i] = -99.0;
        jetPileupE[i] = -99.0;
        jetPileupId[i] = -99.0;
        jetPartonFlavor[i] = 0;
        jetHadronFlavor[i] = 0;

        //AK8 Jet
        fatJetE[i] = 0.0;
        fatJetPt[i] = 0.0;
        fatJetEta[i] = 0.0;
        fatJetPhi[i] = 0.0;
        fatJetPrunedM[i] = 0.0;
        fatJetTrimmedM[i] = 0.0;
        fatJetFilteredM[i] = 0.0;
        fatJetTau1[i] = 0.0;
        fatJetTau2[i] = 0.0;
        fatJetTau3[i] = 0.0;

        genJetE[i] = 0.0;
        genJetPt[i] = 0.0;
        genJetEta[i] = 0.0;
        genJetPhi[i] = 0.0;
    }

    for(int i = 0; i < 500; i++){
        //Gen Particle
        gParticleMotherId[i] = -99999;
        gParticleMotherIndex[i] = -99999;
        gParticleId[i] = -99999;
        gParticleStatus[i] = -99999;
        gParticleE[i] = -99999.0;
        gParticlePt[i] = -99999.0;
        gParticleEta[i] = -99999.0;
        gParticlePhi[i] = -99999.0;

    }

    //MET
    metPt = -999;
    metPhi = -999;
    sumMET = -99.0;
    UncMETdpx = -99.0;
    UncMETdpy = -99.0;
    UncMETdSumEt = -99.0;

    genMetPt = -999;
    genMetPhi = -999;

    MR = -999;
    RSQ = -999;

    //Event
    nPV = -1;
    eventNum = 0;
    lumiNum = 0;
    runNum = 0;
    pvX = -99.0;
    pvY = -99.0;
    pvZ = -99.0;
    fixedGridRhoAll = -99.0;
    fixedGridRhoFastjetAll = -99.0;
    fixedGridRhoFastjetAllCalo = -99.0;
    fixedGridRhoFastjetCentralCalo = -99.0;
    fixedGridRhoFastjetCentralChargedPileUp = -99.0;
    fixedGridRhoFastjetCentralNeutral = -99.0;

}

//------ Methods to fill tree variables ------//

bool RazorTuplizer::fillEventInfo(const edm::Event& iEvent){
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

  //get rho
  fixedGridRhoAll = *rhoAll;
  fixedGridRhoFastjetAll = *rhoFastjetAll;
  fixedGridRhoFastjetAllCalo = *rhoFastjetAllCalo;
  fixedGridRhoFastjetCentralCalo = *rhoFastjetCentralCalo;
  fixedGridRhoFastjetCentralChargedPileUp = *rhoFastjetCentralChargedPileUp;
  fixedGridRhoFastjetCentralNeutral = *rhoFastjetCentralNeutral;

  return true;
}

bool RazorTuplizer::fillPileUp(){
  for(const PileupSummaryInfo &pu : *puInfo){
    BunchXing[nBunchXing] = pu.getBunchCrossing();
    nPU[nBunchXing] = pu.getPU_NumInteractions();
    nPUmean[nBunchXing] = pu.getTrueNumInteractions();
    nBunchXing++;
    //std::cout << "BC: " << pu.getBunchCrossing() << std::endl;
  }
  return true;
};

bool RazorTuplizer::fillMuons(){
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
    muon_pileupIso[nMuons] = mu.pfIsolationR04().sumPUPt;
    muon_chargedIso[nMuons] = mu.pfIsolationR04().sumChargedHadronPt;
    muon_photonIso[nMuons] = mu.pfIsolationR04().sumPhotonEt;
    muon_neutralHadIso[nMuons] = mu.pfIsolationR04().sumNeutralHadronEt;
    muon_ptrel[nMuons] = getLeptonPtRel( jets, &mu );
    muon_miniiso[nMuons] = getPFMiniIsolation(packedPFCands, dynamic_cast<const reco::Candidate *>(&mu), 0.05, 0.2, 10., false, false);

    nMuons++;
  }

  return true;
};

bool RazorTuplizer::fillElectrons(){
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
    ele_pileupIso[nElectrons] = ele.pfIsolationVariables().sumPUPt;
    ele_chargedIso[nElectrons] = ele.pfIsolationVariables().sumChargedHadronPt;
    ele_photonIso[nElectrons] = ele.pfIsolationVariables().sumPhotonEt;
    ele_neutralHadIso[nElectrons] = ele.pfIsolationVariables().sumNeutralHadronEt;
    ele_MissHits[nElectrons] = ele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);

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
    ele_ptrel[nElectrons] = getLeptonPtRel( jets, &ele );
    ele_miniiso[nElectrons] = getPFMiniIsolation(packedPFCands, dynamic_cast<const reco::Candidate *>(&ele), 0.05, 0.2, 10., false, false);

    nElectrons++;
  }
  
  return true;
};

bool RazorTuplizer::fillTaus(){
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

bool RazorTuplizer::fillIsoPFCandidates(){
  for (const pat::PackedCandidate &candidate : *packedPFCands) {

    if (candidate.charge() != 0 && candidate.pt() > 5 && candidate.fromPV() == 3 ) {
      double tmpIsoPFNoPU = 0;
      double tmpIsoPFPU = 0;
      for (const pat::PackedCandidate &isoCandidate : *packedPFCands) {	
	if ( (candidate.pdgId() != 1 && candidate.pdgId() != 2)
	    && deltaR(candidate.eta(), candidate.phi(), isoCandidate.eta(), isoCandidate.phi()) < 0.4
	     && !(candidate.eta() == isoCandidate.eta() && candidate.phi() == isoCandidate.phi())
	    ) {
	  if (candidate.fromPV() == 2 || candidate.fromPV() == 3) {
	    tmpIsoPFNoPU += isoCandidate.pt();
	  } else if (candidate.fromPV() == 0) {
	    tmpIsoPFPU += isoCandidate.pt();
	  }
	}
      }

      if ( 
	  (candidate.pt() > 50 ) ||
	  (candidate.pt() > 20 && (tmpIsoPFNoPU - 0.5*tmpIsoPFPU)/candidate.pt() < 3.0) ||
	  (candidate.pt() <= 20 && tmpIsoPFNoPU - 0.5*tmpIsoPFPU < 25)
	   ) {

	isoPFCandidatePt[nIsoPFCandidates] = candidate.pt();
	isoPFCandidateEta[nIsoPFCandidates] = candidate.eta();
	isoPFCandidatePhi[nIsoPFCandidates] = candidate.phi();
	isoPFCandidateIso04[nIsoPFCandidates] = max(0.0, tmpIsoPFNoPU - 0.5*tmpIsoPFPU) ;
	isoPFCandidateD0[nIsoPFCandidates] = candidate.dxy();
	isoPFCandidatePdgId[nIsoPFCandidates] = candidate.pdgId();
	

	nIsoPFCandidates++;

	//For Debugging
	// cout << "\n";
	// cout << candidate.charge() << " " << candidate.pdgId() << " " << candidate.pt() << " " << candidate.eta() << " " << candidate.phi() 
	//      << " | " << candidate.fromPV() << " " << candidate.dz() ;
	// cout << " | " << tmpIsoPFNoPU << " " <<tmpIsoPFPU <<  " " <<  tmpIsoPFNoPU - 0.5*tmpIsoPFPU << " " << (tmpIsoPFNoPU - 0.5*tmpIsoPFPU)/candidate.pt()  ;
	// cout << "\n";
	for (const pat::PackedCandidate &isoCandidate : *packedPFCands) {	
	  if ( (candidate.pdgId() != 1 && candidate.pdgId() != 2)
	       && deltaR(candidate.eta(), candidate.phi(), isoCandidate.eta(), isoCandidate.phi()) < 0.4
	       && !(candidate.eta() == isoCandidate.eta() && candidate.phi() == isoCandidate.phi())
	       ) {
	    //cout << "isoCandidate " << isoCandidate.pdgId() << " " << isoCandidate.pt() << " " << isoCandidate.eta() << " " << isoCandidate.phi() << " | " << deltaR(candidate.eta(), candidate.phi(), isoCandidate.eta(), isoCandidate.phi()) << " " << candidate.fromPV() << "\n";	  
	  }
	}

      } // if candidate passes isolation
    }
  }
  return true;
}

bool RazorTuplizer::fillPhotons(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  noZS::EcalClusterLazyTools *lazyToolnoZS = new noZS::EcalClusterLazyTools(iEvent, iSetup, ebRecHitsToken_, eeRecHitsToken_);
    
  for (const pat::Photon &pho : *photons) {
    if (pho.pt() < 20) continue;

    std::vector<float> vCov = lazyToolnoZS->localCovariances( *(pho.superCluster()->seed()) );

    phoE[nPhotons] = pho.energy();
    phoPt[nPhotons] = pho.pt();
    phoEta[nPhotons] = pho.eta();
    phoPhi[nPhotons] = pho.phi();
    phoSigmaIetaIeta[nPhotons] = pho.see();
    phoFull5x5SigmaIetaIeta[nPhotons] = (isnan(vCov[0]) ? 0. : sqrt(vCov[0]));

    //phoR9[nPhotons] = pho.r9();
    //Use the noZS version of this according to Emanuele
    phoR9[nPhotons] = lazyToolnoZS->e3x3( *(pho.superCluster()->seed()) ) / pho.superCluster()->rawEnergy();

    pho_HoverE[nPhotons] = pho.hadTowOverEm();
    pho_isConversion[nPhotons] = pho.hasConversionTracks();
    pho_passEleVeto[nPhotons] = !hasMatchedPromptElectron(pho.superCluster(),electrons, 
									   conversions, beamSpot->position());

    //Don't use default miniAOD quantities for now
    // pho_sumChargedHadronPt[nPhotons] = pho.chargedHadronIso();
    // pho_sumNeutralHadronEt[nPhotons] = pho.neutralHadronIso();
    // pho_sumPhotonEt[nPhotons] = pho.photonIso();

    //**********************************************************
    //Compute PF isolation
    //absolute uncorrected isolations with footprint removal
    //**********************************************************
    const float coneSizeDR = 0.3;
    const float dxyMax = 0.1;
    const float dzMax = 0.2;
    float chargedIsoSum = 0;
    float neutralHadronIsoSum = 0;
    float photonIsoSum = 0;

    // First, find photon direction with respect to the good PV
    const reco::Vertex &pv = vertices->front();
    math::XYZVector photon_directionWrtVtx(pho.superCluster()->x() - pv.x(),
					   pho.superCluster()->y() - pv.y(),
					   pho.superCluster()->z() - pv.z());
    // Loop over all PF candidates
    for (const pat::PackedCandidate &candidate : *packedPFCands) {

     // Check if this candidate is within the isolation cone
      float dR=deltaR(photon_directionWrtVtx.Eta(),photon_directionWrtVtx.Phi(),
		      candidate.eta(), candidate.phi());
      if( dR > coneSizeDR ) continue;

      // Check if this candidate is not in the footprint
      bool inFootprint = false;      
      for (auto itr : pho.associatedPackedPFCandidates()) {	
	if ( &(*itr) == &candidate) {
	  inFootprint = true;
	}
      }     
      if( inFootprint ) continue;


      // Find candidate type
      reco::PFCandidate::ParticleType thisCandidateType = reco::PFCandidate::X;

      // the neutral hadrons and charged hadrons can be of pdgId types
      // only 130 (K0L) and +-211 (pi+-) in packed candidates
      const int pdgId = candidate.pdgId();
      if( pdgId == 22 )
	thisCandidateType = reco::PFCandidate::gamma;
      else if( abs(pdgId) == 130) // PDG ID for K0L
	thisCandidateType = reco::PFCandidate::h0;
      else if( abs(pdgId) == 211) // PDG ID for pi+-
	thisCandidateType = reco::PFCandidate::h;
      

      // Increment the appropriate isolation sum
      if( thisCandidateType == reco::PFCandidate::h ){
	// for charged hadrons, additionally check consistency
	// with the PV
	float dxy = -999, dz = -999;
	dz = candidate.pseudoTrack().dz(pv.position());
	dxy =candidate.pseudoTrack().dxy(pv.position());
	if (fabs(dz) > dzMax) continue;
	if(fabs(dxy) > dxyMax) continue;
	// The candidate is eligible, increment the isolaiton
	chargedIsoSum += candidate.pt();
      }
      if( thisCandidateType == reco::PFCandidate::h0 )
	neutralHadronIsoSum += candidate.pt();
      if( thisCandidateType == reco::PFCandidate::gamma )
	photonIsoSum += candidate.pt();
    }
    pho_sumChargedHadronPt[nPhotons] = chargedIsoSum;
    pho_sumNeutralHadronEt[nPhotons] = neutralHadronIsoSum;
    pho_sumPhotonEt[nPhotons] = photonIsoSum;
    

    //*****************************************************************
    //Computer Worst Isolation Looping over all vertices
    //*****************************************************************
    const double ptMin = 0.0;
    const float dRvetoBarrel = 0.0;
    const float dRvetoEndcap = 0.0;    
    float dRveto = 0;
    if (pho.isEB()) dRveto = dRvetoBarrel;
    else dRveto = dRvetoEndcap;
    
    float worstIsolation = 999;
    std::vector<float> allIsolations;
    for(unsigned int ivtx=0; ivtx<vertices->size(); ++ivtx) {
    
      // Shift the photon according to the vertex
      reco::VertexRef vtx(vertices, ivtx);
      math::XYZVector photon_directionWrtVtx(pho.superCluster()->x() - vtx->x(),
					     pho.superCluster()->y() - vtx->y(),
					     pho.superCluster()->z() - vtx->z());
    
      float sum = 0;
      // Loop over all PF candidates
      for (const pat::PackedCandidate &candidate : *packedPFCands) {
		
	//require that PFCandidate is a charged hadron
	const int pdgId = candidate.pdgId();
	if( abs(pdgId) != 211) continue;

	if (candidate.pt() < ptMin)
	  continue;
      
	float dxy = -999, dz = -999;
	dz = candidate.pseudoTrack().dz(pv.position());
	dxy =candidate.pseudoTrack().dxy(pv.position());
	if( fabs(dxy) > dxyMax) continue;     
	if ( fabs(dz) > dzMax) continue;
	
	float dR = deltaR(photon_directionWrtVtx.Eta(), photon_directionWrtVtx.Phi(), 
			  candidate.eta(),      candidate.phi());
	if(dR > coneSizeDR || dR < dRveto) continue;
	
	sum += candidate.pt();
      }
      
      allIsolations.push_back(sum);
    }

    if( allIsolations.size()>0 )
      worstIsolation = * std::max_element( allIsolations.begin(), allIsolations.end() );
    

    //*****************************************************************
    //Photon ID MVA variable
    //*****************************************************************
    pho_IDMVA[nPhotons] = myPhotonMVA->mvaValue( pho,  *rhoAll, photonIsoSum, chargedIsoSum, worstIsolation,
						 lazyToolnoZS, false);
				       
    pho_RegressionE[nPhotons] = pho.getCorrectedEnergy(reco::Photon::P4type::regression1);
    pho_RegressionEUncertainty[nPhotons] = pho.getCorrectedEnergyError(reco::Photon::P4type::regression1);
    
    pho_superClusterEta[nPhotons] = pho.superCluster()->eta();
    pho_superClusterPhi[nPhotons] = pho.superCluster()->phi();
    pho_hasPixelSeed[nPhotons] = pho.hasPixelSeed();

    /*
    const reco::Candidate* genPhoton = pho.genPhoton();
    if(genPhoton != NULL)std::cout << "======>gen PT: " << genPhoton->pt() <<
      " recoPT: " << pho.pt() << std::endl;
    */
    nPhotons++;
  }

  delete lazyToolnoZS;
  return true;
};

bool RazorTuplizer::fillJets(){
  for (const pat::Jet &j : *jets) {
    if (j.pt() < 20) continue;
    jetE[nJets] = j.correctedP4(0).E();
    jetPt[nJets] = j.correctedP4(0).Pt();
    jetEta[nJets] = j.correctedP4(0).Eta();
    jetPhi[nJets] = j.correctedP4(0).Phi();
    jetCSV[nJets] = j.bDiscriminator("pfCombinedSecondaryVertexBJetTags");
    jetCISV[nJets] = j.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags");
    jetMass[nJets] = j.correctedP4(0).M();
    jetJetArea[nJets] = j.jetArea();
    jetPileupE[nJets] = j.pileup();
    jetPileupId[nJets] = j.userFloat("pileupJetId:fullDiscriminant");
    jetPartonFlavor[nJets] = j.partonFlavour();
    jetHadronFlavor[nJets] = j.hadronFlavour();
    nJets++;
  }

  return true;
};

bool RazorTuplizer::fillJetsAK8(){
  for (const pat::Jet &j : *jetsAK8) {
    fatJetE[nFatJets] = j.correctedP4(0).E();
    fatJetPt[nFatJets] = j.correctedP4(0).Pt();
    fatJetEta[nFatJets] = j.correctedP4(0).Eta();
    fatJetPhi[nFatJets] = j.correctedP4(0).Phi();
    fatJetPrunedM[nFatJets] = (float) j.userFloat("ak8PFJetsCHSPrunedLinks");
    fatJetTrimmedM[nFatJets] = (float) j.userFloat("ak8PFJetsCHSTrimmedLinks");
    fatJetFilteredM[nFatJets] = (float) j.userFloat("ak8PFJetsCHSFilteredLinks");
    fatJetTau1[nFatJets] =  (float) j.userFloat("NjettinessAK8:tau1");
    fatJetTau2[nFatJets] =  (float) j.userFloat("NjettinessAK8:tau2");
    fatJetTau3[nFatJets] =  (float) j.userFloat("NjettinessAK8:tau3");
    nFatJets++;
  }

  return true;
};

bool RazorTuplizer::fillMet(const edm::Event& iEvent){
  const pat::MET &Met = mets->front();
  metPt = Met.pt();
  metPhi = Met.phi();
  sumMET = Met.sumEt();

  //MET filters
  const edm::TriggerNames &metNames = iEvent.triggerNames(*metFilterBits);
  for(unsigned int i = 0, n = metFilterBits->size(); i < n; ++i){
    if(strcmp(metNames.triggerName(i).c_str(), "Flag_trackingFailureFilter") == 0)
      Flag_trackingFailureFilter = metFilterBits->accept(i);
    else if(strcmp(metNames.triggerName(i).c_str(), "Flag_goodVertices") == 0)
      Flag_goodVertices = metFilterBits->accept(i);
    else if(strcmp(metNames.triggerName(i).c_str(), "Flag_CSCTightHaloFilter") == 0)
      Flag_CSCTightHaloFilter = metFilterBits->accept(i);
    else if(strcmp(metNames.triggerName(i).c_str(), "Flag_trkPOGFilters") == 0)
      Flag_trkPOGFilters = metFilterBits->accept(i);
    else if(strcmp(metNames.triggerName(i).c_str(), "Flag_trkPOG_logErrorTooManyClusters") == 0)
      Flag_trkPOG_logErrorTooManyClusters = metFilterBits->accept(i);
    else if(strcmp(metNames.triggerName(i).c_str(), "Flag_EcalDeadCellTriggerPrimitiveFilter") == 0)
      Flag_EcalDeadCellTriggerPrimitiveFilter = metFilterBits->accept(i);
    else if(strcmp(metNames.triggerName(i).c_str(), "Flag_ecalLaserCorrFilter") == 0)
      Flag_ecalLaserCorrFilter = metFilterBits->accept(i);
    else if(strcmp(metNames.triggerName(i).c_str(), "Flag_trkPOG_manystripclus53X") == 0)
      Flag_trkPOG_manystripclus53X = metFilterBits->accept(i);
    else if(strcmp(metNames.triggerName(i).c_str(), "Flag_eeBadScFilter") == 0)
      Flag_eeBadScFilter = metFilterBits->accept(i);
    else if(strcmp(metNames.triggerName(i).c_str(), "Flag_METFilters") == 0)
      Flag_METFilters = metFilterBits->accept(i);
    else if(strcmp(metNames.triggerName(i).c_str(), "Flag_HBHENoiseFilter") == 0)
      Flag_HBHENoiseFilter = metFilterBits->accept(i);
    else if(strcmp(metNames.triggerName(i).c_str(), "Flag_trkPOG_toomanystripclus53X") == 0)
      Flag_trkPOG_toomanystripclus53X = metFilterBits->accept(i);
    else if(strcmp(metNames.triggerName(i).c_str(), "Flag_hcalLaserEventFilter") == 0)
      Flag_hcalLaserEventFilter = metFilterBits->accept(i);
  }
  
  return true;
};

bool RazorTuplizer::fillMC(){
    for(const reco::GenJet &j : *genJets){
        genJetE[nGenJets] = j.energy();
        genJetPt[nGenJets] = j.pt();
        genJetEta[nGenJets] = j.eta();
        genJetPhi[nGenJets] = j.phi();
        nGenJets++;
    }

    const pat::MET &Met = mets->front();
    genMetPt = Met.genMET()->pt();
    genMetPhi = Met.genMET()->phi();

    bool foundGenVertex = false;
    for(size_t i=0; i<prunedGenParticles->size();i++){
      if (!foundGenVertex) {
	for (unsigned int j=0; j<(*prunedGenParticles)[i].numberOfDaughters(); ++j) {
	  const reco::Candidate *dau = (*prunedGenParticles)[i].daughter(j);
	  if (dau) {
	    genVertexX = dau->vx();
	    genVertexY = dau->vy();
	    genVertexZ = dau->vz();
	    foundGenVertex = true;
	    break;
	  }
	}
      }
    }
 
    genWeight = genInfo->weight();
    genSignalProcessID = genInfo->signalProcessID();
    genQScale = genInfo->qScale();
    genAlphaQCD = genInfo->alphaQCD();
    genAlphaQED = genInfo->alphaQED();

    return true;
}

bool RazorTuplizer::fillRazor(){ 
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

  //compute the razor variables using the selected jets and MET
  if(goodJets.size() > 1){
    vector<TLorentzVector> hemispheres = getHemispheres(goodJets);
    TLorentzVector pfMet(Met.px(), Met.py(), 0.0, 0.0);
    MR = computeMR(hemispheres[0], hemispheres[1]);
    RSQ = computeR2(hemispheres[0], hemispheres[1], pfMet);
  }
                                
  vector<TLorentzVector> goodJets_AK8;
    for (const pat::Jet &j : *jetsAK8) {
      if (j.pt() < 40) continue;
      if (fabs(j.eta()) > 3.0) continue;
      TLorentzVector newJet(j.px(), j.py(), j.pz(), j.energy());
      goodJets_AK8.push_back(newJet);
    }

  //compute the razor variables using the selected jets and MET
  if(goodJets.size() > 1){
    vector<TLorentzVector> hemispheres = getHemispheres(goodJets);
    TLorentzVector pfMet(Met.px(), Met.py(), 0.0, 0.0);
    MR_AK8 = computeMR(hemispheres[0], hemispheres[1]);
    RSQ_AK8 = computeR2(hemispheres[0], hemispheres[1], pfMet);
  }

  
  return true;
}

bool RazorTuplizer::fillGenParticles(){
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

bool RazorTuplizer::fillTrigger(const edm::Event& iEvent){
  //fill trigger information

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  //std::cout << "\n === TRIGGER PATHS === " << std::endl;

  // //********************************************************************
  // //Option to save all HLT path names in the ntuple per event
  // //Expensive option in terms of ntuple size
  // //********************************************************************
  // for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
  //   string hltPathNameReq = "HLT_";   
  //   if (triggerBits->accept(i)) 
  //     if ((names.triggerName(i)).find(hltPathNameReq) != string::npos) 
  // 	nameHLT->push_back(names.triggerName(i));

  //   std::cout << "Trigger " << names.triggerName(i) << 
  //     ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
  //     ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)") 
  //   	      << std::endl;
  // }

  //********************************************************************
  // Save trigger decisions in array of booleans
  //********************************************************************
  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
    string hltPathNameReq = "HLT_";   
    
    if ((names.triggerName(i)).find(hltPathNameReq) == string::npos) continue;
    if ((names.triggerName(i)).find_last_of("_") == string::npos) continue;
    int lastUnderscorePos = (names.triggerName(i)).find_last_of("_");
    string hltPathNameWithoutVersionNumber = (names.triggerName(i)).substr(0,lastUnderscorePos);
    
    for (unsigned int j = 0; j < NTriggersMAX; ++j) {
      if (triggerPathNames[j] == "") continue;
      if (hltPathNameWithoutVersionNumber == triggerPathNames[j]) {
	triggerDecision[j] = triggerBits->accept(i);
      }
    }    
  }
    
  return true;
}


//------ Method called for each event ------//

void RazorTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;
  
  //initialize
  resetBranches();
  loadEvent(iEvent); //loads objects and resets tree branches
  
  NEvents->Fill(0); //increment event count

  //filler methods should fill relevant tree variables and return false if the event should be rejected
  bool isGoodEvent = 
    fillEventInfo(iEvent)
    && fillPileUp()
    && fillMuons() 
    && fillElectrons()
    && fillTaus()
    && fillIsoPFCandidates()
    && fillPhotons(iEvent,iSetup)
    && fillJets()
    && fillJetsAK8()
    && fillMet(iEvent)
    && fillRazor()
    && fillMC()
    && fillGenParticles();
  //NOTE: if any of the above functions return false, the event will be rejected immediately with no further processing

  if (enableTriggerInfo_) isGoodEvent = (isGoodEvent && fillTrigger(iEvent));
  
  //fill the tree if the event wasn't rejected
  if(isGoodEvent) RazorEvents->Fill();
}

//------ Method called once each job just before starting event loop ------//
void RazorTuplizer::beginJob(){
  setBranches();
}

//------ Method called once each job just after ending the event loop ------//
void RazorTuplizer::endJob(){
}

//define this as a plug-in
DEFINE_FWK_MODULE(RazorTuplizer);

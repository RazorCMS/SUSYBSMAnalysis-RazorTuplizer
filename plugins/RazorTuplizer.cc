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
  isData_(iConfig.getParameter<bool> ("isData")),
  useGen_(iConfig.getParameter<bool> ("useGen")),  
  isFastsim_(iConfig.getParameter<bool> ("isFastsim")),  
  enableTriggerInfo_(iConfig.getParameter<bool> ("enableTriggerInfo")),
  triggerPathNamesFile_(iConfig.getParameter<string> ("triggerPathNamesFile")),
  eleHLTFilterNamesFile_(iConfig.getParameter<string> ("eleHLTFilterNamesFile")),
  muonHLTFilterNamesFile_(iConfig.getParameter<string> ("muonHLTFilterNamesFile")),
  photonHLTFilterNamesFile_(iConfig.getParameter<string> ("photonHLTFilterNamesFile")),
  verticesToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  tracksTag_(consumes<edm::View<reco::Track> >(iConfig.getParameter<edm::InputTag>("tracks"))),
  trackTimeTag_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("trackTime"))),
  trackTimeResoTag_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("trackTimeReso"))),
  muonsToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronsToken_(consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  tausToken_(consumes<reco::PFTauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  photonsToken_(consumes<reco::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
  jetsToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  jetsPuppiToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jetsPuppi"))),
  jetsAK8Token_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("jetsAK8"))),
  PFCandsToken_(consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
  PFClustersToken_(consumes<reco::PFClusterCollection>(iConfig.getParameter<edm::InputTag>("pfClusters"))),
  //genParticlesToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"))),
  //genParticlesToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"))),
  genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
  genJetsToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets"))),
  triggerBitsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"))),
  hepMCToken_(consumes<edm::HepMCProduct>(iConfig.getParameter<edm::InputTag>("hepMC"))),
  //triggerObjectsToken_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
  //triggerPrescalesToken_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("triggerPrescales"))),     
  metToken_(consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
//  metNoHFToken_(consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("metsNoHF"))),
  metPuppiToken_(consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("metsPuppi"))),
  metFilterBitsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilterBits"))),
  //hbheNoiseFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("hbheNoiseFilter"))),
  //hbheTightNoiseFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("hbheTightNoiseFilter"))),
  //hbheIsoNoiseFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("hbheIsoNoiseFilter"))),
  //badChargedCandidateFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("BadChargedCandidateFilter"))),
  //badMuonFilterToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("BadMuonFilter"))),
//  lheRunInfoTag_(iConfig.getParameter<edm::InputTag>("lheInfo")),
//  lheRunInfoToken_(consumes<LHERunInfoProduct,edm::InRun>(lheRunInfoTag_)),
//  lheInfoToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheInfo"))),
  genInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genInfo"))),
  genLumiHeaderToken_(consumes<GenLumiInfoHeader,edm::InLumi>(edm::InputTag("generator",""))),
  puInfoToken_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("puInfo"))),
  //hcalNoiseInfoToken_(consumes<HcalNoiseSummary>(iConfig.getParameter<edm::InputTag>("hcalNoiseInfo"))),
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
  gedPhotonCoresToken_(consumes<vector<reco::PhotonCore> >(iConfig.getParameter<edm::InputTag>("gedPhotonCores")))
  //superClustersToken_(consumes<vector<reco::SuperCluster> >(iConfig.getParameter<edm::InputTag>("superClusters"))),
//  lostTracksToken_(consumes<vector<reco::PFCandidate> >(iConfig.getParameter<edm::InputTag>("lostTracks")))
{
  //declare the TFileService for output
  edm::Service<TFileService> fs;
  
  //set up output tree
  RazorEvents = fs->make<TTree>("RazorEvents", "selected miniAOD information");
  NEvents = fs->make<TH1F>("NEvents",";;NEvents;",1,-0.5,0.5);
  if (useGen_) {
    sumWeights = fs->make<TH1D>("sumWeights",";;sumWeights;",1,-0.5,0.5);
    sumScaleWeights = fs->make<TH1D>("sumScaleWeights",";;sumScaleWeights;",9,-0.5,8.5);
    sumPdfWeights = fs->make<TH1D>("sumPdfWeights",";;sumPdfWeights;",100,-0.5,99.5);
    sumAlphasWeights = fs->make<TH1D>("sumAlphasWeights",";;sumAlphasWeights;",2,-0.5,1.5);
    
    sumWeights->Sumw2();
    sumScaleWeights->Sumw2();
    sumPdfWeights->Sumw2();
    sumAlphasWeights->Sumw2();
  }
  else {
    sumWeights = 0;
    sumScaleWeights = 0;
    sumPdfWeights = 0;
    sumAlphasWeights = 0;
  }

  //set up electron MVA ID
  std::vector<std::string> myTrigWeights;
  myTrigWeights.push_back(edm::FileInPath("SUSYBSMAnalysis/RazorTuplizer/data/TrigIDMVA_25ns_EB_BDT.weights.xml").fullPath().c_str());
  myTrigWeights.push_back(edm::FileInPath("SUSYBSMAnalysis/RazorTuplizer/data/TrigIDMVA_25ns_EE_BDT.weights.xml").fullPath().c_str());

  myMVATrig = new EGammaMvaEleEstimatorCSA14();
  myMVATrig->initialize("BDT",
			EGammaMvaEleEstimatorCSA14::kTrig,
			true,
			myTrigWeights);

  std::vector<std::string> myNonTrigWeights;
  myNonTrigWeights.push_back(edm::FileInPath("SUSYBSMAnalysis/RazorTuplizer/data/EIDmva_EB1_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml").fullPath().c_str());
  myNonTrigWeights.push_back(edm::FileInPath("SUSYBSMAnalysis/RazorTuplizer/data/EIDmva_EB2_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml").fullPath().c_str());
  myNonTrigWeights.push_back(edm::FileInPath("SUSYBSMAnalysis/RazorTuplizer/data/EIDmva_EE_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml").fullPath().c_str());
  myNonTrigWeights.push_back(edm::FileInPath("SUSYBSMAnalysis/RazorTuplizer/data/EIDmva_EB1_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml").fullPath().c_str());
  myNonTrigWeights.push_back(edm::FileInPath("SUSYBSMAnalysis/RazorTuplizer/data/EIDmva_EB2_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml").fullPath().c_str());
  myNonTrigWeights.push_back(edm::FileInPath("SUSYBSMAnalysis/RazorTuplizer/data/EIDmva_EE_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml").fullPath().c_str());
  
  myMVANonTrig = new ElectronMVAEstimatorRun2NonTrig();
  myMVANonTrig->initialize("BDTG method",
			ElectronMVAEstimatorRun2NonTrig::kPHYS14,
			true,
			myNonTrigWeights);

  //set up photon MVA ID
  std::vector<std::string> myPhotonMVAWeights;
  myPhotonMVAWeights.push_back(edm::FileInPath("SUSYBSMAnalysis/RazorTuplizer/data/PhotonIDMVA_Spring15_50ns_v0_EB.weights.xml").fullPath().c_str());
  myPhotonMVAWeights.push_back(edm::FileInPath("SUSYBSMAnalysis/RazorTuplizer/data/PhotonIDMVA_Spring15_50ns_v0_EE.weights.xml").fullPath().c_str());
  std::vector<std::string> myPhotonMVAMethodNames;
  myPhotonMVAMethodNames.push_back("BDTG photons barrel");
  myPhotonMVAMethodNames.push_back("BDTG photons endcap");

  myPhotonMVA = new EGammaMvaPhotonEstimator();
  myPhotonMVA->initialize(myPhotonMVAMethodNames,myPhotonMVAWeights,
			  EGammaMvaPhotonEstimator::kPhotonMVATypeDefault);


  //*****************************************************************************************
  //Read in HLT Trigger Path List from config file
  //*****************************************************************************************
  for (int i = 0; i<NTriggersMAX; ++i) triggerPathNames[i] = "";
  ifstream myfile (edm::FileInPath(triggerPathNamesFile_.c_str()).fullPath().c_str()) ;
  if (myfile.is_open()) {
    string line;
    int index;
    string hltpathname;

    while(myfile>>index>>hltpathname) {
      
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

  //*****************************************************************************************
  //Read in Electron HLT Filters List from config file
  //*****************************************************************************************
  for (int i = 0; i<MAX_ElectronHLTFilters; ++i) eleHLTFilterNames[i] = "";
  ifstream myEleHLTFilterFile (edm::FileInPath(eleHLTFilterNamesFile_.c_str()).fullPath().c_str()) ;
  if (myEleHLTFilterFile.is_open()) {
    char tmp[1024];
    string line;
    int index;
    string hltfiltername;

    while(myEleHLTFilterFile>>line) {
      
      if ( line.empty() || line.substr(0,1) == "#") {
	myEleHLTFilterFile.getline(tmp,1024);
	continue;
      }

      index = atoi(line.c_str());
      myEleHLTFilterFile >> hltfiltername;
      
      if (index < MAX_ElectronHLTFilters) {
	eleHLTFilterNames[index] = hltfiltername;
      }    
    }    
    myEleHLTFilterFile.close();
  } else {
    cout << "ERROR!!! Could not open trigger path name file : " << edm::FileInPath(eleHLTFilterNamesFile_.c_str()).fullPath().c_str() << "\n";
  }
  
  if(enableTriggerInfo_) {
    cout << "\n";
    cout << "****************** Electron HLT Filters Defined For Razor Ntuple ******************\n";    
    for (int i = 0; i<MAX_ElectronHLTFilters; ++i) {
      if (eleHLTFilterNames[i] != "") cout << "Ele HLT Filters " << i << " " << eleHLTFilterNames[i] << "\n";
    }
    cout << "****************************************************************************\n";    
    cout << "\n";
  }


  //*****************************************************************************************
  //Read in Muon HLT Filters List from config file
  //*****************************************************************************************
  for (int i = 0; i<MAX_MuonHLTFilters; ++i) muonHLTFilterNames[i] = "";
  ifstream myMuonHLTFilterFile (edm::FileInPath(muonHLTFilterNamesFile_.c_str()).fullPath().c_str()) ;
  if (myMuonHLTFilterFile.is_open()) {
    char tmp[1024];
    string line;
    int index;
    string hltfiltername;

    while(myMuonHLTFilterFile>>line) {
      
      if ( line.empty() || line.substr(0,1) == "#") {
	myMuonHLTFilterFile.getline(tmp,1024);
	continue;
      }

      index = atoi(line.c_str());
      myMuonHLTFilterFile >> hltfiltername;
      
      if (index < MAX_MuonHLTFilters) {
	muonHLTFilterNames[index] = hltfiltername;
      }    
    }    
    myMuonHLTFilterFile.close();
  } else {
    cout << "ERROR!!! Could not open trigger path name file : " << edm::FileInPath(muonHLTFilterNamesFile_.c_str()).fullPath().c_str() << "\n";
  }
  
  if(enableTriggerInfo_) {
    cout << "\n";
    cout << "****************** Muon HLT Filters Defined For Razor Ntuple ******************\n";    
    for (int i = 0; i<MAX_MuonHLTFilters; ++i) {
      if (muonHLTFilterNames[i] != "") cout << "Muon HLT Filters " << i << " " << muonHLTFilterNames[i] << "\n";
    }
    cout << "****************************************************************************\n";    
    cout << "\n";
  }



  //*****************************************************************************************
  //Read in Photon HLT Filters List from config file
  //*****************************************************************************************
  for (int i = 0; i<MAX_PhotonHLTFilters; ++i) photonHLTFilterNames[i] = "";
  ifstream myPhotonHLTFilterFile (edm::FileInPath(photonHLTFilterNamesFile_.c_str()).fullPath().c_str()) ;
  if (myPhotonHLTFilterFile.is_open()) {
    char tmp[1024];
    string line;
    int index;
    string hltfiltername;

    while(myPhotonHLTFilterFile>>line) {
      
      if ( line.empty() || line.substr(0,1) == "#") {
	myPhotonHLTFilterFile.getline(tmp,1024);
	continue;
      }

      index = atoi(line.c_str());
      myPhotonHLTFilterFile >> hltfiltername;
      
      if (index < MAX_PhotonHLTFilters) {
	photonHLTFilterNames[index] = hltfiltername;
      }    
    }    
    myPhotonHLTFilterFile.close();
  } else {
    cout << "ERROR!!! Could not open trigger path name file : " << edm::FileInPath(photonHLTFilterNamesFile_.c_str()).fullPath().c_str() << "\n";
  }
  
  if(enableTriggerInfo_) {
    cout << "\n";
    cout << "****************** Photon HLT Filters Defined For Razor Ntuple ******************\n";    
    for (int i = 0; i<MAX_PhotonHLTFilters; ++i) {
      if (photonHLTFilterNames[i] != "") cout << "Photon HLT Filters " << i << " " << photonHLTFilterNames[i] << "\n";
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
  enablePVAllBranches();
  enablePileUpBranches();
  enableMuonBranches();
  enableElectronBranches();
  enableTauBranches();
  enableIsoPFCandidateBranches();
  enablePhotonBranches();
  enableJetBranches();
  enableJetAK8Branches();
  enableMetBranches();
  
  if (enableTriggerInfo_) enableTriggerBranches();
  enableMCBranches();
  enableGenParticleBranches();
}

void RazorTuplizer::enableEventInfoBranches(){
  RazorEvents->Branch("isData", &isData, "isData/O");
  RazorEvents->Branch("nPV", &nPV, "nPV/I");
  RazorEvents->Branch("runNum", &runNum, "runNum/i");
  RazorEvents->Branch("lumiNum", &lumiNum, "lumiNum/i");
  RazorEvents->Branch("eventNum", &eventNum, "eventNum/i");
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

void RazorTuplizer::enablePVAllBranches() {
  RazorEvents->Branch("nPVAll", &nPVAll,"nPVAll/I");
  RazorEvents->Branch("pvAllX", pvAllX,"pvAllX[nPVAll]/F");
  RazorEvents->Branch("pvAllXError", pvAllXError,"pvAllXError[nPVAll]/F");
  RazorEvents->Branch("pvAllY", pvAllY,"pvAllY[nPVAll]/F");
  RazorEvents->Branch("pvAllYError", pvAllYError,"pvAllYError[nPVAll]/F");
  RazorEvents->Branch("pvAllZ", pvAllZ,"pvAllZ[nPVAll]/F");
  RazorEvents->Branch("pvAllZError", pvAllZError,"pvAllZError[nPVAll]/F");
  RazorEvents->Branch("pvAllT", pvAllT,"pvAllT[nPVAll]/F");
  RazorEvents->Branch("pvAllTError", pvAllTError,"pvAllTError[nPVAll]/F");
  RazorEvents->Branch("pvAllLogSumPtSq", pvAllLogSumPtSq,"pvAllLogSumPtSq[nPVAll]/F");
  RazorEvents->Branch("pvAllSumPx", pvAllSumPx,"pvAllSumPx[nPVAll]/F");
  RazorEvents->Branch("pvAllSumPy", pvAllSumPy,"pvAllSumPy[nPVAll]/F");

  RazorEvents->Branch("pvAllLogSumPtSq_dt", pvAllLogSumPtSq_dt,"pvAllLogSumPtSq_dt[nPVAll]/F");
  RazorEvents->Branch("pvAllSumPx_dt", pvAllSumPx_dt,"pvAllSumPx_dt[nPVAll]/F");
  RazorEvents->Branch("pvAllSumPy_dt", pvAllSumPy_dt,"pvAllSumPy_dt[nPVAll]/F");

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
  RazorEvents->Branch("muonIsMedium", muonIsMedium,"muonIsMedium[nMuons]/O");
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
  RazorEvents->Branch("muon_chargedMiniIso", muon_chargedMiniIso, "muon_chargedMiniIso[nMuons]/F");
  RazorEvents->Branch("muon_photonAndNeutralHadronMiniIso", muon_photonAndNeutralHadronMiniIso, "muon_photonAndNeutralHadronMiniIso[nMuons]/F");
  RazorEvents->Branch("muon_chargedPileupMiniIso", muon_chargedPileupMiniIso, "muon_chargedPileupMiniIso[nMuons]/F");
  RazorEvents->Branch("muon_activityMiniIsoAnnulus", muon_activityMiniIsoAnnulus, "muon_activityMiniIsoAnnulus[nMuons]/F");
  RazorEvents->Branch("muon_passSingleMuTagFilter", muon_passSingleMuTagFilter, "muon_passSingleMuTagFilter[nMuons]/O");
  RazorEvents->Branch("muon_passHLTFilter", &muon_passHLTFilter, Form("muon_passHLTFilter[nMuons][%d]/O",MAX_MuonHLTFilters));
  RazorEvents->Branch("muon_validFractionTrackerHits", muon_validFractionTrackerHits, "muon_validFractionTrackerHits[nMuons]/F");
  RazorEvents->Branch("muon_isGlobal", muon_isGlobal,"muon_isGlobal[nMuons]/O");
  RazorEvents->Branch("muon_normChi2", muon_normChi2,"muon_normChi2[nMuons]/F");
  RazorEvents->Branch("muon_chi2LocalPosition", muon_chi2LocalPosition,"muon_chi2LocalPosition[nMuons]/F");
  RazorEvents->Branch("muon_kinkFinder", muon_kinkFinder,"muon_kinkFinder[nMuons]/F");
  RazorEvents->Branch("muon_segmentCompatability", muon_segmentCompatability,"muon_segmentCompatability[nMuons]/F");
  RazorEvents->Branch("muonIsICHEPMedium", muonIsICHEPMedium,"muonIsICHEPMedium[nMuons]/O");
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
  RazorEvents->Branch("ele_ip3d", ele_ip3d, "ele_ip3d[nElectrons]/F");
  RazorEvents->Branch("ele_ip3dSignificance", ele_ip3dSignificance, "ele_ip3dSignificance[nElectrons]/F");
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
  RazorEvents->Branch("ele_chargedMiniIso", ele_chargedMiniIso, "ele_chargedMiniIso[nElectrons]/F");
  RazorEvents->Branch("ele_photonAndNeutralHadronMiniIso", ele_photonAndNeutralHadronMiniIso, "ele_photonAndNeutralHadronMiniIso[nElectrons]/F");
  RazorEvents->Branch("ele_chargedPileupMiniIso", ele_chargedPileupMiniIso, "ele_chargedPileupMiniIso[nElectrons]/F");
  RazorEvents->Branch("ele_activityMiniIsoAnnulus", ele_activityMiniIsoAnnulus, "ele_activityMiniIsoAnnulus[nElectrons]/F");
  RazorEvents->Branch("ele_passSingleEleTagFilter", ele_passSingleEleTagFilter, "ele_passSingleEleTagFilter[nElectrons]/O");
  RazorEvents->Branch("ele_passTPOneTagFilter", ele_passTPOneTagFilter, "ele_passTPOneTagFilter[nElectrons]/O");
  RazorEvents->Branch("ele_passTPTwoTagFilter", ele_passTPTwoTagFilter, "ele_passTPTwoTagFilter[nElectrons]/O");
  RazorEvents->Branch("ele_passTPOneProbeFilter", ele_passTPOneProbeFilter, "ele_passTPOneProbeFilter[nElectrons]/O");
  RazorEvents->Branch("ele_passTPTwoProbeFilter", ele_passTPTwoProbeFilter, "ele_passTPTwoProbeFilter[nElectrons]/O");
  RazorEvents->Branch("ele_passHLTFilter", &ele_passHLTFilter, Form("ele_passHLTFilter[nElectrons][%d]/O",MAX_ElectronHLTFilters));
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
  RazorEvents->Branch("tau_chargedIsoPtSum", tau_chargedIsoPtSum, "tau_chargedIsoPtSum[nTaus]/F");
  RazorEvents->Branch("tau_neutralIsoPtSum", tau_neutralIsoPtSum, "tau_neutralIsoPtSum[nTaus]/F");
  RazorEvents->Branch("tau_puCorrPtSum", tau_puCorrPtSum, "tau_puCorrPtSum[nTaus]/F");
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
  RazorEvents->Branch("pho_sumChargedHadronPtAllVertices", &pho_sumChargedHadronPtAllVertices,Form("pho_sumChargedHadronPtAllVertices[nPhotons][%d]/F",MAX_NPV));
  RazorEvents->Branch("pho_sumChargedHadronPt", &pho_sumChargedHadronPt, "pho_sumChargedHadronPt[nPhotons]/F");
  RazorEvents->Branch("pho_sumNeutralHadronEt", pho_sumNeutralHadronEt, "pho_sumNeutralHadronEt[nPhotons]/F");
  RazorEvents->Branch("pho_sumPhotonEt", pho_sumPhotonEt, "pho_sumPhotonEt[nPhotons]/F");
  RazorEvents->Branch("pho_sumWorstVertexChargedHadronPt", pho_sumWorstVertexChargedHadronPt, "pho_sumWorstVertexChargedHadronPt[nPhotons]/F");
  RazorEvents->Branch("pho_pfIsoChargedHadronIso", pho_pfIsoChargedHadronIso, "pho_pfIsoChargedHadronIso[nPhotons]/F");
  RazorEvents->Branch("pho_pfIsoChargedHadronIsoWrongVtx", pho_pfIsoChargedHadronIsoWrongVtx, "pho_pfIsoChargedHadronIsoWrongVtx[nPhotons]/F");
  RazorEvents->Branch("pho_pfIsoNeutralHadronIso", pho_pfIsoNeutralHadronIso, "pho_pfIsoNeutralHadronIso[nPhotons]/F");
  RazorEvents->Branch("pho_pfIsoPhotonIso", pho_pfIsoPhotonIso, "pho_pfIsoPhotonIso[nPhotons]/F");
  RazorEvents->Branch("pho_pfIsoModFrixione", pho_pfIsoModFrixione, "pho_pfIsoModFrixione[nPhotons]/F");
  RazorEvents->Branch("pho_pfIsoSumPUPt", pho_pfIsoSumPUPt, "pho_pfIsoSumPUPt[nPhotons]/F");
  RazorEvents->Branch("pho_isConversion", pho_isConversion, "pho_isConversion[nPhotons]/O");
  RazorEvents->Branch("pho_passEleVeto", pho_passEleVeto, "pho_passEleVeto[nPhotons]/O");
  RazorEvents->Branch("pho_RegressionE", pho_RegressionE, "pho_RegressionE[nPhotons]/F");
  RazorEvents->Branch("pho_RegressionEUncertainty", pho_RegressionEUncertainty, "pho_RegressionEUncertainty[nPhotons]/F");
  RazorEvents->Branch("pho_IDMVA", pho_IDMVA, "pho_IDMVA[nPhotons]/F");
  RazorEvents->Branch("pho_superClusterEnergy", pho_superClusterEnergy, "pho_superClusterEnergy[nPhotons]/F");
  RazorEvents->Branch("pho_superClusterRawEnergy", pho_superClusterRawEnergy, "pho_superClusterRawEnergy[nPhotons]/F");
  RazorEvents->Branch("pho_superClusterEta", pho_superClusterEta, "pho_superClusterEta[nPhotons]/F");
  RazorEvents->Branch("pho_superClusterPhi", pho_superClusterPhi, "pho_superClusterPhi[nPhotons]/F");
  RazorEvents->Branch("pho_superClusterX", pho_superClusterX, "pho_superClusterX[nPhotons]/F");
  RazorEvents->Branch("pho_superClusterY", pho_superClusterY, "pho_superClusterY[nPhotons]/F");
  RazorEvents->Branch("pho_superClusterZ", pho_superClusterZ, "pho_superClusterZ[nPhotons]/F");

  RazorEvents->Branch("pho_superClusterSeedX", pho_superClusterSeedX, "pho_superClusterSeedX[nPhotons]/F");
  RazorEvents->Branch("pho_superClusterSeedY", pho_superClusterSeedY, "pho_superClusterSeedY[nPhotons]/F");
  RazorEvents->Branch("pho_superClusterSeedZ", pho_superClusterSeedZ, "pho_superClusterSeedZ[nPhotons]/F");
  RazorEvents->Branch("pho_superClusterSeedT", pho_superClusterSeedT, "pho_superClusterSeedT[nPhotons]/F");
  
//  RazorEvents->Branch("pho_superClusterSeedXError", pho_superClusterSeedXError, "pho_superClusterSeedXError[nPhotons]/F");
//  RazorEvents->Branch("pho_superClusterSeedYError", pho_superClusterSeedYError, "pho_superClusterSeedYError[nPhotons]/F");
//  RazorEvents->Branch("pho_superClusterSeedZError", pho_superClusterSeedZError, "pho_superClusterSeedZError[nPhotons]/F");
//  RazorEvents->Branch("pho_superClusterSeedTError", pho_superClusterSeedT, "pho_superClusterSeedT[nPhotons]/F");
 

  RazorEvents->Branch("pho_hasPixelSeed", pho_hasPixelSeed, "pho_hasPixelSeed[nPhotons]/O");
  RazorEvents->Branch("pho_passHLTFilter", &pho_passHLTFilter, Form("pho_passHLTFilter[nPhotons][%d]/O",MAX_PhotonHLTFilters));
  RazorEvents->Branch("pho_convType", pho_convType, "pho_convType[nPhotons]/I");
  RazorEvents->Branch("pho_convTrkZ", pho_convTrkZ, "pho_convTrkZ[nPhotons]/F");
  RazorEvents->Branch("pho_convTrkClusZ", pho_convTrkClusZ, "pho_convTrkClusZ[nPhotons]/F");
  RazorEvents->Branch("pho_vtxSumPx", &pho_vtxSumPx,Form("pho_vtxSumPx[nPhotons][%d]/F",MAX_NPV));
  RazorEvents->Branch("pho_vtxSumPy", &pho_vtxSumPy,Form("pho_vtxSumPy[nPhotons][%d]/F",MAX_NPV));



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
  RazorEvents->Branch("jetPileupIdFlag", jetPileupIdFlag, "jetPileupIdFlag[nJets]/I");
  RazorEvents->Branch("jetPassIDLoose", jetPassIDLoose, "jetPassIDLoose[nJets]/O");
  RazorEvents->Branch("jetPassIDTight", jetPassIDTight, "jetPassIDTight[nJets]/O");
  RazorEvents->Branch("jetPassMuFrac", jetPassMuFrac, "jetPassMuFrac[nJets]/O");
  RazorEvents->Branch("jetPassEleFrac", jetPassEleFrac, "jetPassEleFrac[nJets]/O");
  RazorEvents->Branch("jetPartonFlavor", jetPartonFlavor, "jetPartonFlavor[nJets]/I");
  RazorEvents->Branch("jetHadronFlavor", jetHadronFlavor, "jetHadronFlavor[nJets]/I");
  RazorEvents->Branch("jetChargedEMEnergyFraction", jetChargedEMEnergyFraction, "jetChargedEMEnergyFraction[nJets]/F"); 
  RazorEvents->Branch("jetNeutralEMEnergyFraction", jetNeutralEMEnergyFraction, "jetNeutralEMEnergyFraction[nJets]/F"); 
  RazorEvents->Branch("jetChargedHadronEnergyFraction", jetChargedHadronEnergyFraction, "jetChargedHadronEnergyFraction[nJets]/F"); 
  RazorEvents->Branch("jetNeutralHadronEnergyFraction", jetNeutralHadronEnergyFraction, "jetNeutralHadronEnergyFraction[nJets]/F"); 
  RazorEvents->Branch("jetMuonEnergyFraction", jetMuonEnergyFraction, "jetMuonEnergyFraction[nJets]/F"); 
  RazorEvents->Branch("jetHOEnergyFraction", jetHOEnergyFraction, "jetHOEnergyFraction[nJets]/F");
  RazorEvents->Branch("jetHFHadronEnergyFraction", jetHFHadronEnergyFraction, "jetHFHadronEnergyFraction[nJets]/F");
  RazorEvents->Branch("jetHFEMEnergyFraction",jetHFEMEnergyFraction, "jetHFEMEnergyFraction[nJets]/F");
  RazorEvents->Branch("jetAllMuonPt", jetAllMuonPt,"jetAllMuonPt[nJets]/F");
  RazorEvents->Branch("jetAllMuonEta", jetAllMuonEta,"jetAllMuonEta[nJets]/F");
  RazorEvents->Branch("jetAllMuonPhi", jetAllMuonPhi,"jetAllMuonPhi[nJets]/F");
  RazorEvents->Branch("jetAllMuonM", jetAllMuonM,"jetAllMuonM[nJets]/F");
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
  RazorEvents->Branch("metType0Pt", &metType0Pt, "metType0Pt/F");
  RazorEvents->Branch("metType0Phi", &metType0Phi, "metType0Phi/F");
  RazorEvents->Branch("metType1Pt", &metType1Pt, "metType1Pt/F");
  RazorEvents->Branch("metType1Phi", &metType1Phi, "metType1Phi/F");
  RazorEvents->Branch("metType0Plus1Pt", &metType0Plus1Pt, "metType0Plus1Pt/F");
  RazorEvents->Branch("metType0Plus1Phi", &metType0Plus1Phi, "metType0Plus1Phi/F");
  RazorEvents->Branch("metNoHFPt", &metNoHFPt, "metNoHFPt/F");
  RazorEvents->Branch("metNoHFPhi", &metNoHFPhi, "metNoHFPhi/F");
  RazorEvents->Branch("metPuppiPt", &metPuppiPt, "metPuppiPt/F");
  RazorEvents->Branch("metPuppiPhi", &metPuppiPhi, "metPuppiPhi/F");
  RazorEvents->Branch("metCaloPt", &metCaloPt, "metCaloPt/F");
  RazorEvents->Branch("metCaloPhi", &metCaloPhi, "metCaloPhi/F");

  RazorEvents->Branch("metType1PtJetResUp", &metType1PtJetResUp, "metType1PtJetResUp/F");
  RazorEvents->Branch("metType1PtJetResDown", &metType1PtJetResDown, "metType1PtJetResDown/F");
  RazorEvents->Branch("metType1PtJetEnUp", &metType1PtJetEnUp, "metType1PtJetEnUp/F");
  RazorEvents->Branch("metType1PtJetEnDown", &metType1PtJetEnDown, "metType1PtJetEnDown/F");
  RazorEvents->Branch("metType1PtMuonEnUp", &metType1PtMuonEnUp, "metType1PtMuonEnUp/F");
  RazorEvents->Branch("metType1PtMuonEnDown", &metType1PtMuonEnDown, "metType1PtMuonEnDown/F");
  RazorEvents->Branch("metType1PtElectronEnUp", &metType1PtElectronEnUp, "metType1PtElectronEnUp/F");
  RazorEvents->Branch("metType1PtElectronEnDown", &metType1PtElectronEnDown, "metType1PtElectronEnDown/F");
  RazorEvents->Branch("metType1PtTauEnUp", &metType1PtTauEnUp, "metType1PtTauEnUp/F");
  RazorEvents->Branch("metType1PtTauEnDown", &metType1PtTauEnDown, "metType1PtTauEnDown/F");
  RazorEvents->Branch("metType1PtUnclusteredEnUp", &metType1PtUnclusteredEnUp, "metType1PtUnclusteredEnUp/F");
  RazorEvents->Branch("metType1PtUnclusteredEnDown", &metType1PtUnclusteredEnDown, "metType1PtUnclusteredEnDown/F");
  RazorEvents->Branch("metType1PtPhotonEnUp", &metType1PtPhotonEnUp, "metType1PtPhotonEnUp/F");
  RazorEvents->Branch("metType1PtPhotonEnDown", &metType1PtPhotonEnDown, "metType1PtPhotonEnDown/F");

  RazorEvents->Branch("metType1PhiJetResUp", &metType1PhiJetResUp, "metType1PhiJetResUp/F");
  RazorEvents->Branch("metType1PhiJetResDown", &metType1PhiJetResDown, "metType1PhiJetResDown/F");
  RazorEvents->Branch("metType1PhiJetEnUp", &metType1PhiJetEnUp, "metType1PhiJetEnUp/F");
  RazorEvents->Branch("metType1PhiJetEnDown", &metType1PhiJetEnDown, "metType1PhiJetEnDown/F");
  RazorEvents->Branch("metType1PhiMuonEnUp", &metType1PhiMuonEnUp, "metType1PhiMuonEnUp/F");
  RazorEvents->Branch("metType1PhiMuonEnDown", &metType1PhiMuonEnDown, "metType1PhiMuonEnDown/F");
  RazorEvents->Branch("metType1PhiElectronEnUp", &metType1PhiElectronEnUp, "metType1PhiElectronEnUp/F");
  RazorEvents->Branch("metType1PhiElectronEnDown", &metType1PhiElectronEnDown, "metType1PhiElectronEnDown/F");
  RazorEvents->Branch("metType1PhiTauEnUp", &metType1PhiTauEnUp, "metType1PhiTauEnUp/F");
  RazorEvents->Branch("metType1PhiTauEnDown", &metType1PhiTauEnDown, "metType1PhiTauEnDown/F");
  RazorEvents->Branch("metType1PhiUnclusteredEnUp", &metType1PhiUnclusteredEnUp, "metType1PhiUnclusteredEnUp/F");
  RazorEvents->Branch("metType1PhiUnclusteredEnDown", &metType1PhiUnclusteredEnDown, "metType1PhiUnclusteredEnDown/F");
  RazorEvents->Branch("metType1PhiPhotonEnUp", &metType1PhiPhotonEnUp, "metType1PhiPhotonEnUp/F");
  RazorEvents->Branch("metType1PhiPhotonEnDown", &metType1PhiPhotonEnDown, "metType1PhiPhotonEnDown/F");

  RazorEvents->Branch("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, "Flag_HBHENoiseFilter/O");
  RazorEvents->Branch("Flag_HBHETightNoiseFilter", &Flag_HBHETightNoiseFilter, "Flag_HBHETightNoiseFilter/O");
  RazorEvents->Branch("Flag_HBHEIsoNoiseFilter", &Flag_HBHEIsoNoiseFilter, "Flag_HBHEIsoNoiseFilter/O");
  //RazorEvents->Branch("Flag_badChargedCandidateFilter", &Flag_badChargedCandidateFilter, "Flag_badChargedCandidateFilter/O");
  //RazorEvents->Branch("Flag_badMuonFilter", &Flag_badMuonFilter, "Flag_badMuonFilter/O");
  RazorEvents->Branch("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, "Flag_CSCTightHaloFilter/O");
  RazorEvents->Branch("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, "Flag_hcalLaserEventFilter/O");
  RazorEvents->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, "Flag_EcalDeadCellTriggerPrimitiveFilter/O");
  RazorEvents->Branch("Flag_EcalDeadCellBoundaryEnergyFilter", &Flag_EcalDeadCellBoundaryEnergyFilter, "Flag_EcalDeadCellBoundaryEnergyFilter/O");
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
  RazorEvents->Branch("HLTDecision", &triggerDecision, ("HLTDecision[" + std::to_string(NTriggersMAX) +  "]/O").c_str());
  RazorEvents->Branch("HLTPrescale", &triggerHLTPrescale, ("HLTPrescale[" + std::to_string(NTriggersMAX) +  "]/I").c_str());
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
  RazorEvents->Branch("genVertexX", &genVertexX, "genVertexX/F");
  RazorEvents->Branch("genVertexY", &genVertexY, "genVertexY/F");
  RazorEvents->Branch("genVertexT", &genVertexT, "genVertexT/F");
  RazorEvents->Branch("genWeight", &genWeight, "genWeight/F");
  RazorEvents->Branch("genSignalProcessID", &genSignalProcessID, "genSignalProcessID/i");
  RazorEvents->Branch("genQScale", &genQScale, "genQScale/F");
  RazorEvents->Branch("genAlphaQCD", &genAlphaQCD, "genAlphaQCD/F");
  RazorEvents->Branch("genAlphaQED", &genAlphaQED, "genAlphaQED/F");
  scaleWeights = new std::vector<float>; scaleWeights->clear();
  pdfWeights = new std::vector<float>; pdfWeights->clear();
  alphasWeights = new std::vector<float>; alphasWeights->clear();
  if (isFastsim_) {
    RazorEvents->Branch("lheComments", "std::string",&lheComments);
  }
  RazorEvents->Branch("scaleWeights", "std::vector<float>",&scaleWeights);
  RazorEvents->Branch("pdfWeights", "std::vector<float>",&pdfWeights);
  RazorEvents->Branch("alphasWeights", "std::vector<float>",&alphasWeights);
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
  iEvent.getByToken(hepMCToken_, hepMC);
//  iEvent.getByToken(triggerObjectsToken_, triggerObjects);
//  iEvent.getByToken(triggerPrescalesToken_, triggerPrescales);
  iEvent.getByToken(metFilterBitsToken_, metFilterBits);
  iEvent.getByToken(verticesToken_, vertices);
  iEvent.getByToken(tracksTag_,tracks);
  iEvent.getByToken(trackTimeTag_,times);
  iEvent.getByToken(trackTimeResoTag_,timeResos);
  iEvent.getByToken(PFCandsToken_, pfCands);
  iEvent.getByToken(PFClustersToken_, pfClusters);
  iEvent.getByToken(muonsToken_, muons);
  iEvent.getByToken(electronsToken_, electrons);
  iEvent.getByToken(photonsToken_, photons);
  iEvent.getByToken(tausToken_, taus);
  iEvent.getByToken(jetsToken_, jets);
  iEvent.getByToken(jetsPuppiToken_, jetsPuppi);
  iEvent.getByToken(jetsAK8Token_, jetsAK8);
  iEvent.getByToken(metToken_, mets);
  //iEvent.getByToken(metNoHFToken_, metsNoHF);
  iEvent.getByToken(metPuppiToken_, metsPuppi);
//  iEvent.getByToken(hcalNoiseInfoToken_,hcalNoiseInfo);
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
//  iEvent.getByToken(superClustersToken_,superClusters);
//  iEvent.getByToken(lostTracksToken_,lostTracks);
//  iEvent.getByToken(hbheNoiseFilterToken_, hbheNoiseFilter);
//  iEvent.getByToken(hbheTightNoiseFilterToken_, hbheTightNoiseFilter);
//  iEvent.getByToken(hbheIsoNoiseFilterToken_, hbheIsoNoiseFilter);
  //iEvent.getByToken(badChargedCandidateFilterToken_, badChargedCandidateFilter);
  //iEvent.getByToken(badMuonFilterToken_, badMuonFilter);

  if (useGen_) {
//    iEvent.getByToken(genParticlesToken_,genParticles);
    iEvent.getByToken(genParticlesToken_,genParticles);
    iEvent.getByToken(genJetsToken_,genJets);

    //for Spring16 fastsim, this has been changed and removed
//    if (!isFastsim_) iEvent.getByToken(lheInfoToken_, lheInfo);     
    
    iEvent.getByToken(genInfoToken_,genInfo);
    iEvent.getByToken(puInfoToken_,puInfo);
  }
  

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

    for(int i = 0; i < NTriggersMAX; i++){
      triggerDecision[i] = false;
      triggerHLTPrescale[i] = 0;
    }

    for (int i=0; i < MAX_NPV; ++i) {
      //pvAll
      pvAllX[i] = 0.;
      pvAllY[i] = 0.;
      pvAllZ[i] = 0.;
      pvAllT[i] = 0.;
      pvAllXError[i] = 0.;
      pvAllYError[i] = 0.;
      pvAllZError[i] = 0.;
      pvAllTError[i] = 0.;

      pvAllLogSumPtSq[i] = 0.;
      pvAllSumPx[i] = 0.;
      pvAllSumPy[i] = 0.;
      
      pvAllLogSumPtSq_dt[i] = 0.;
      pvAllSumPx_dt[i] = 0.;
      pvAllSumPy_dt[i] = 0.;

    }
    
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
        muonIsMedium[i] = false;
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
	muon_chargedMiniIso[i] = -99.0;
	muon_photonAndNeutralHadronMiniIso[i] = -99.0;
	muon_chargedPileupMiniIso[i] = -99.0;
	muon_activityMiniIsoAnnulus[i] = -99.0;
	muon_passSingleMuTagFilter[i] = false;
	for (int q=0;q<MAX_MuonHLTFilters;q++) muon_passHLTFilter[i][q] = false;
	muon_validFractionTrackerHits[i] = -99.0;
        muon_isGlobal[i] = false;
        muon_normChi2[i] = -99.0;
        muon_chi2LocalPosition[i] = -99.0;
        muon_kinkFinder[i] = -99.0;
        muon_segmentCompatability[i] = -99.0;
        muonIsICHEPMedium[i] = false;

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
	ele_ip3d[i] = -99;
	ele_ip3dSignificance[i] = -99;
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
	ele_chargedMiniIso[i] = -99.0;
	ele_photonAndNeutralHadronMiniIso[i] = -99.0;
	ele_chargedPileupMiniIso[i] = -99.0;
	ele_activityMiniIsoAnnulus[i] = -99.0;
	ele_passSingleEleTagFilter[i] = false;
	ele_passTPOneTagFilter[i] = false;
	ele_passTPTwoTagFilter[i] = false;
	ele_passTPOneProbeFilter[i] = false;
	ele_passTPTwoProbeFilter[i] = false;
	for (int q=0;q<MAX_ElectronHLTFilters;q++) ele_passHLTFilter[i][q] = false;

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
	tau_chargedIsoPtSum[i] = -99.0;
	tau_neutralIsoPtSum[i] = -99.0;
	tau_puCorrPtSum[i] = -99.0;
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
	pho_sumWorstVertexChargedHadronPt[i] = -99.0;
	pho_pfIsoChargedHadronIso[i] = -99.0;
	pho_pfIsoChargedHadronIsoWrongVtx[i] = -99.0;
	pho_pfIsoNeutralHadronIso[i] = -99.0;
	pho_pfIsoPhotonIso[i] = -99.0;
	pho_pfIsoModFrixione[i] = -99.0;
	pho_pfIsoSumPUPt[i] = -99.0;       
	pho_isConversion[i] = false;
        pho_passEleVeto[i] = false;    
        pho_RegressionE[i] = -99.0;
        pho_RegressionEUncertainty[i] = -99.0;
        pho_IDMVA[i] = -99.0;
	pho_superClusterEnergy[i] = -99.0;
	pho_superClusterRawEnergy[i] = -99.0;
        pho_superClusterEta[i]    = -99.0;
        pho_superClusterPhi[i]    = -99.0;
	pho_superClusterX[i]      = -99.0;
	pho_superClusterY[i]      = -99.0;
	pho_superClusterZ[i]      = -99.0;
	pho_superClusterSeedX[i]      = -99.0;
	pho_superClusterSeedY[i]      = -99.0;
	pho_superClusterSeedZ[i]      = -99.0;
	pho_superClusterSeedT[i]      = -99.0;
	pho_superClusterSeedXError[i]      = -99.0;
	pho_superClusterSeedYError[i]      = -99.0;
	pho_superClusterSeedZError[i]      = -99.0;
	pho_superClusterSeedTError[i]      = -99.0;



        pho_hasPixelSeed[i] = false;
	for (int q=0;q<MAX_PhotonHLTFilters;q++) pho_passHLTFilter[i][q] = false;
        pho_convType[i] = -99;
        pho_convTrkZ[i] = -99.;
        pho_convTrkClusZ[i] = -99.;
        
        for (int ipv=0; ipv < MAX_NPV; ++ipv) {
          pho_sumChargedHadronPtAllVertices[i][ipv] = -99.0;
          pho_vtxSumPx[i][ipv] = 0.;
          pho_vtxSumPy[i][ipv] = 0.;
          

	}

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
	jetPileupIdFlag[i] = -1;
	jetPassIDLoose[i] = false;
	jetPassIDTight[i] = false;
	jetPassMuFrac[i] = false;
	jetPassEleFrac[i] = false;
        jetPartonFlavor[i] = 0;
        jetHadronFlavor[i] = 0;
	jetChargedEMEnergyFraction[i] = -99.0;
	jetNeutralEMEnergyFraction[i] = -99.0;
	jetChargedHadronEnergyFraction[i] = -99.0;
	jetNeutralHadronEnergyFraction[i] = -99.0;
	jetMuonEnergyFraction[i] = -99.0;
	jetHOEnergyFraction[i] = -99.0;
	jetHFHadronEnergyFraction[i] = -99.0;
	jetHFEMEnergyFraction[i] = -99.0;
        jetAllMuonPt[i] = 0.0;
        jetAllMuonEta[i] = 0.0;
        jetAllMuonPhi[i] = 0.0;
        jetAllMuonM[i] = 0.0;
	
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
    metType0Pt = -99.0;
    metType0Phi = -99.0;
    metType1Pt = -99.0;
    metType1Phi = -99.0;
    metType0Plus1Pt = -99.0;
    metType0Plus1Phi = -99.0;
    metPtRecomputed = -99.0;
    metPhiRecomputed = -99.0;
    metNoHFPt = -99.0;
    metNoHFPhi = -99.0;
    metPuppiPt = -99.0;
    metPuppiPhi = -99.0;
    metCaloPt = -999;
    metCaloPhi = -999;
    Flag_HBHENoiseFilter = false;
    Flag_HBHETightNoiseFilter = false;
    Flag_HBHEIsoNoiseFilter = false;
    //Flag_badChargedCandidateFilter = false;
    //Flag_badMuonFilter = false;
    Flag_CSCTightHaloFilter = false;
    Flag_hcalLaserEventFilter = false;
    Flag_EcalDeadCellTriggerPrimitiveFilter = false;
    Flag_EcalDeadCellBoundaryEnergyFilter = false;
    Flag_goodVertices = false;
    Flag_trackingFailureFilter = false;
    Flag_eeBadScFilter = false;
    Flag_ecalLaserCorrFilter = false;
    Flag_trkPOGFilters = false;  
    Flag_trkPOG_manystripclus53X = false;
    Flag_trkPOG_toomanystripclus53X = false;
    Flag_trkPOG_logErrorTooManyClusters = false;
    Flag_METFilters = false;

    metType1PtJetResUp=-999.;
    metType1PtJetResDown=-999.;
    metType1PtJetEnUp=-999.;
    metType1PtJetEnDown=-999.;
    metType1PtMuonEnUp=-999.;
    metType1PtMuonEnDown=-999.;
    metType1PtElectronEnUp=-999.;
    metType1PtElectronEnDown=-999.;
    metType1PtTauEnUp=-999.;
    metType1PtTauEnDown=-999.;
    metType1PtUnclusteredEnUp=-999.;
    metType1PtUnclusteredEnDown=-999.;
    metType1PtPhotonEnUp=-999.;
    metType1PtPhotonEnDown=-999.;
    metType1PtMETUncertaintySize=-999.;
    metType1PtJetResUpSmear=-999.;
    metType1PtJetResDownSmear=-999.;
    metType1PtMETFullUncertaintySize=-999.;
  
    metType1PhiJetResUp=-999.;
    metType1PhiJetResDown=-999.;
    metType1PhiJetEnUp=-999.;
    metType1PhiJetEnDown=-999.;
    metType1PhiMuonEnUp=-999.;
    metType1PhiMuonEnDown=-999.;
    metType1PhiElectronEnUp=-999.;
    metType1PhiElectronEnDown=-999.;
    metType1PhiTauEnUp=-999.;
    metType1PhiTauEnDown=-999.;
    metType1PhiUnclusteredEnUp=-999.;
    metType1PhiUnclusteredEnDown=-999.;
    metType1PhiPhotonEnUp=-999.;
    metType1PhiPhotonEnDown=-999.;
    metType1PhiMETUncertaintySize=-999.;
    metType1PhiJetResUpSmear=-999.;
    metType1PhiJetResDownSmear=-999.;
    metType1PhiMETFullUncertaintySize=-999.;

    genMetPt = -999;
    genMetPhi = -999;
    genVertexZ = -999;
    genVertexX = -999;
    genVertexY = -999;
    genVertexT = -999;
    genWeight = 1;
    genSignalProcessID = -999;
    genQScale = -999;
    genAlphaQCD = -999;
    genAlphaQED = -999;
    scaleWeights->clear();
    pdfWeights->clear();
    alphasWeights->clear();

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
  isData = isData_;
  runNum = iEvent.id().run();
  lumiNum = iEvent.luminosityBlock();
  eventNum = iEvent.id().event();
  
  //select the primary vertex, if any
  nPV = 0;
  myPV = &(vertices->front());
  bool foundPV = false;
  for(unsigned int i = 0; i < vertices->size(); i++){
    if(vertices->at(i).isValid() && !vertices->at(i).isFake()){
      if (!foundPV) {
	myPV = &(vertices->at(i));
	foundPV = true;
      }
      nPV++;
    }
  }
  
  pvX = myPV->x();
  pvY = myPV->y();
  pvZ = myPV->z();

  //get rho
  fixedGridRhoAll = *rhoAll;
  fixedGridRhoFastjetAll = *rhoFastjetAll;
  fixedGridRhoFastjetAllCalo = *rhoFastjetAllCalo;
  fixedGridRhoFastjetCentralCalo = *rhoFastjetCentralCalo;
  fixedGridRhoFastjetCentralChargedPileUp = *rhoFastjetCentralChargedPileUp;
  fixedGridRhoFastjetCentralNeutral = *rhoFastjetCentralNeutral;

  return true;
}

bool RazorTuplizer::fillPVAll() {
  
  nPVAll = std::min(int(vertices->size()),int(MAX_NPV));
  
  for (int ipv = 0; ipv < nPVAll; ++ipv) {
    const reco::Vertex &vtx = vertices->at(ipv);
    pvAllX[ipv] = vtx.x();
    pvAllY[ipv] = vtx.y();
    pvAllZ[ipv] = vtx.z();
    pvAllT[ipv] = vtx.t();
    
    pvAllXError[ipv] = vtx.xError();
    pvAllYError[ipv] = vtx.yError();
    pvAllZError[ipv] = vtx.zError();
    pvAllTError[ipv] = vtx.tError();


  }
  
  double pvAllSumPtSqD[MAX_NPV];
  double pvAllSumPxD[MAX_NPV];
  double pvAllSumPyD[MAX_NPV];
 
  double pvAllSumPtSqD_dt[MAX_NPV];
  double pvAllSumPxD_dt[MAX_NPV];
  double pvAllSumPyD_dt[MAX_NPV];
  
  for (int ipv=0; ipv<nPVAll; ++ipv) {
    pvAllSumPtSqD[ipv] = 0.;
    pvAllSumPxD[ipv] = 0.;
    pvAllSumPyD[ipv] = 0.;

    pvAllSumPtSqD_dt[ipv] = 0.;
    pvAllSumPxD_dt[ipv] = 0.;
    pvAllSumPyD_dt[ipv] = 0.;
 
  }
  
  //for (const reco::PFCandidate &pfcand : *pfCands) {
  for (const reco::PFCandidate &pfcand : *pfCands) {
    if (pfcand.charge()==0) continue;
    double mindz = std::numeric_limits<double>::max();
    int ipvmin = -1;
    for (int ipv = 0; ipv < nPVAll; ++ipv) {
      const reco::Vertex &vtx = vertices->at(ipv);
      //double dz = std::abs(pfcand.dz(vtx.position()));
      double dz = std::abs(pfcand.vz()-vtx.z());
      if (dz<mindz) {
        mindz = dz;
        ipvmin = ipv;
      }
    }
        
    if (mindz<0.2 && ipvmin>=0 && ipvmin<MAX_NPV) {
      pvAllSumPtSqD[ipvmin] += pfcand.pt()*pfcand.pt();
      pvAllSumPxD[ipvmin] += pfcand.px();
      pvAllSumPyD[ipvmin] += pfcand.py();
    }
  //deltat cut between the track and the primary vertex: |deltat|<60ps
  //for pv: vtx.time()
  //for tracks: https://github.com/lgray/cmssw/blob/muon_iso_timing/TimingAnalysis/MuonIsoTiming/plugins/MuonIsoTiming.cc#L278
  //
  
  bool passDeltaTcut = false;

  const reco::Vertex &vtx_ipvmin = vertices->at(ipvmin);
  
  //auto ref = pfcand.trackRef();
 
  //for( unsigned i = 0; i < tracks->size(); ++i ) {
  //auto ref = tracks->refAt(i);
   //std::cout<<"PFcand here......"<<std::endl;
//    if(times->find(pfcand.trackRef())!=times->end())
    if(times->contains(pfcand.trackRef().id()))
     {
     const float time = (*times)[pfcand.trackRef()];
     const float timeReso = (*timeResos)[pfcand.trackRef()] != 0.f ? (*timeResos)[pfcand.trackRef()] : 0.170f;	
     if(timeReso < 0.02)
	{
	   if(std::abs(time-vtx_ipvmin.t())>0.06)
	   {
		passDeltaTcut = true;
	   } 
	}
      }
  //  }
    if (passDeltaTcut) {
      pvAllSumPtSqD_dt[ipvmin] += pfcand.pt()*pfcand.pt();
      pvAllSumPxD_dt[ipvmin] += pfcand.px();
      pvAllSumPyD_dt[ipvmin] += pfcand.py();
    }
    
 
  }
  
  for (int ipv=0; ipv<nPVAll; ++ipv) {
    pvAllLogSumPtSq[ipv] = log(pvAllSumPtSqD[ipv]);
    pvAllSumPx[ipv] = pvAllSumPxD[ipv];
    pvAllSumPy[ipv] = pvAllSumPyD[ipv];    
  
    pvAllLogSumPtSq_dt[ipv] = log(pvAllSumPtSqD_dt[ipv]);
    pvAllSumPx_dt[ipv] = pvAllSumPxD_dt[ipv];
    pvAllSumPy_dt[ipv] = pvAllSumPyD_dt[ipv];    

  }
  
  
  
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
  for(const reco::Muon &mu : *muons){
    if(mu.pt() < 5) continue;
    muonE[nMuons] = mu.energy();
    muonPt[nMuons] = mu.pt();
    muonEta[nMuons] = mu.eta();
    muonPhi[nMuons] = mu.phi();
    muonCharge[nMuons] = mu.charge();
//    muonIsLoose[nMuons] = mu.isLooseMuon();
//    muonIsMedium[nMuons] = mu.isMediumMuon();
//    muonIsTight[nMuons] = mu.isTightMuon(*myPV);
    muon_d0[nMuons] = -mu.muonBestTrack()->dxy(myPV->position());
    muon_dZ[nMuons] = mu.muonBestTrack()->dz(myPV->position());
//    muon_ip3d[nMuons] = mu.dB(reco::Muon::PV3D);
//    muon_ip3dSignificance[nMuons] = mu.dB(reco::Muon::PV3D)/mu.edB(reco::Muon::PV3D);
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
    //muon_ptrel[nMuons] = getLeptonPtRel( jets, &mu );
//    tuple<double,double,double> PFMiniIso = getPFMiniIsolation(pfCands, dynamic_cast<const reco::Candidate *>(&mu), 0.05, 0.2, 10., false, false);
//    muon_chargedMiniIso[nMuons] = std::get<0>(PFMiniIso);
//    muon_photonAndNeutralHadronMiniIso[nMuons] = std::get<1>(PFMiniIso);
//    muon_chargedPileupMiniIso[nMuons] = std::get<2>(PFMiniIso);
//    muon_activityMiniIsoAnnulus[nMuons] = ActivityPFMiniIsolationAnnulus( pfCands, dynamic_cast<const reco::Candidate *>(&mu), 0.4, 0.05, 0.2, 10.);
    muon_validFractionTrackerHits[nMuons] = (mu.innerTrack().isNonnull() ? mu.track()->validFraction() : -99.0);
    muon_isGlobal[nMuons] = muon::isGoodMuon(mu,muon::AllGlobalMuons);
    muon_normChi2[nMuons] = ( muon::isGoodMuon(mu,muon::AllGlobalMuons) ? mu.globalTrack()->normalizedChi2() : -99.0);
    muon_chi2LocalPosition[nMuons] = mu.combinedQuality().chi2LocalPosition;
    muon_kinkFinder[nMuons] = mu.combinedQuality().trkKink;
    muon_segmentCompatability[nMuons] = muon::segmentCompatibility(mu);

    bool isGoodGlobal = mu.isGlobalMuon() && mu.globalTrack()->normalizedChi2() < 3 && mu.combinedQuality().chi2LocalPosition < 12 && mu.combinedQuality().trkKink < 20;
    muonIsICHEPMedium[nMuons] = muon::isLooseMuon(mu) && muon_validFractionTrackerHits[nMuons] > 0.49 && muon::segmentCompatibility(mu) > (isGoodGlobal ? 0.303 : 0.451);

    //*************************************************
    //Trigger Object Matching
    //*************************************************
    bool passTagMuonFilter = false;
/*
    for (pat::TriggerObjectStandAlone trigObject : *triggerObjects) {

      if (deltaR(trigObject.eta(), trigObject.phi(),mu.eta(),mu.phi()) > 0.3) continue;

      //check single muon filters
      if ( trigObject.hasFilterLabel("hltL3fL1sMu25L1f0Tkf27QL3trkIsoFiltered0p09") ||
	   trigObject.hasFilterLabel("hltL3fL1sMu20Eta2p1L1f0Tkf24QL3trkIsoFiltered0p09") ||
	   trigObject.hasFilterLabel("hltL3fL1sMu16Eta2p1L1f0Tkf20QL3trkIsoFiltered0p09") ||
	   trigObject.hasFilterLabel("hltL3fL1sMu16L1f0Tkf20QL3trkIsoFiltered0p09") ||
	   trigObject.hasFilterLabel("hltL3crIsoL1sMu25L1f0L2f10QL3f27QL3trkIsoFiltered0p09") ||
	   trigObject.hasFilterLabel("hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09") ||
	   trigObject.hasFilterLabel("hltL3crIsoL1sMu16Eta2p1L1f0L2f10QL3f20QL3trkIsoFiltered0p09") ||
	   trigObject.hasFilterLabel("hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09")
	   ) passTagMuonFilter = true;

      //check all filters
      for ( int q=0; q<MAX_MuonHLTFilters;q++) {
	if (trigObject.hasFilterLabel(muonHLTFilterNames[q].c_str())) muon_passHLTFilter[nMuons][q] = true;
      }

    }
*/

    muon_passSingleMuTagFilter[nMuons] = passTagMuonFilter;

    nMuons++;
  }

  return true;
};

bool RazorTuplizer::fillElectrons(){
  for(const reco::GsfElectron &ele : *electrons){
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
    ele_d0[nElectrons] = -ele.gsfTrack().get()->dxy(myPV->position());
    ele_dZ[nElectrons] = ele.gsfTrack().get()->dz(myPV->position());
    //ele_ip3d[nElectrons] = ele.dB(reco::GsfElectron::PV3D);
    //ele_ip3dSignificance[nElectrons] = ele.dB(reco::GsfElectron::PV3D)/ele.edB(reco::GsfElectron::PV3D);   
    ele_pileupIso[nElectrons] = ele.pfIsolationVariables().sumPUPt;
    ele_chargedIso[nElectrons] = ele.pfIsolationVariables().sumChargedHadronPt;
    ele_photonIso[nElectrons] = ele.pfIsolationVariables().sumPhotonEt;
    ele_neutralHadIso[nElectrons] = ele.pfIsolationVariables().sumNeutralHadronEt;
    ele_MissHits[nElectrons] = ele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);

    //*************************************************
    //Conversion Veto
    //*************************************************
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

    //*************************************************
    //ID MVA
    //*************************************************
    ele_IDMVATrig[nElectrons]    = myMVATrig->mvaValue(ele, false);
    ele_IDMVANonTrig[nElectrons] = myMVANonTrig->mvaValue(ele,conversions, beamSpot->position(),false);

    //ele_RegressionE[nElectrons] = ele.ecalRegressionEnergy();
    //ele_CombineP4[nElectrons]   = ele.ecalTrackRegressionEnergy();

    //ele_ptrel[nElectrons]   = getLeptonPtRel( jets, &ele );
//    tuple<double,double,double> PFMiniIso = getPFMiniIsolation(pfCands, dynamic_cast<const reco::Candidate *>(&ele), 0.05, 0.2, 10., false, false);
//    ele_chargedMiniIso[nElectrons] = std::get<0>(PFMiniIso);
//    ele_photonAndNeutralHadronMiniIso[nElectrons] = std::get<1>(PFMiniIso);
//    ele_chargedPileupMiniIso[nElectrons] = std::get<2>(PFMiniIso);
//    ele_activityMiniIsoAnnulus[nElectrons] = ActivityPFMiniIsolationAnnulus( pfCands, dynamic_cast<const reco::Candidate *>(&ele), 0.4, 0.05, 0.2, 10.);

    //*************************************************
    //Trigger Object Matching
    //*************************************************
    bool passSingleEleTagFilter = false;
    bool passTPOneTagFilter= false;
    bool passTPTwoTagFilter= false;
    bool passTPOneProbeFilter= false;
    bool passTPTwoProbeFilter= false;

/*
    for (pat::TriggerObjectStandAlone trigObject : *triggerObjects) {
      if (deltaR(trigObject.eta(), trigObject.phi(),ele.eta(),ele.phi()) > 0.3) continue;

      //check Single ele filters
      if (trigObject.hasFilterLabel("hltEle23WPLooseGsfTrackIsoFilter")  ||
	  trigObject.hasFilterLabel("hltEle27WPLooseGsfTrackIsoFilter")  ||
	  trigObject.hasFilterLabel("hltEle27WPTightGsfTrackIsoFilter")  ||
	  trigObject.hasFilterLabel("hltEle32WPLooseGsfTrackIsoFilter")  ||
	  trigObject.hasFilterLabel("hltEle32WPTightGsfTrackIsoFilter")
	  ) {
	passSingleEleTagFilter = true;
      }
      
      //check Tag and Probe Filters
      if (trigObject.hasFilterLabel("hltEle25WP60Ele8TrackIsoFilter")) passTPOneTagFilter = true;
      if (trigObject.hasFilterLabel("hltEle25WP60SC4TrackIsoFilter")) passTPTwoTagFilter = true;
      if (trigObject.hasFilterLabel("hltEle25WP60Ele8Mass55Filter")) passTPOneProbeFilter = true;
      if (trigObject.hasFilterLabel("hltEle25WP60SC4Mass55Filter")) passTPTwoProbeFilter = true;

      //check all filters
      for ( int q=0; q<MAX_ElectronHLTFilters;q++) {
	if (trigObject.hasFilterLabel(eleHLTFilterNames[q].c_str())) ele_passHLTFilter[nElectrons][q] = true;
      }

    }
  */

    ele_passSingleEleTagFilter[nElectrons] = passSingleEleTagFilter;
    ele_passTPOneTagFilter[nElectrons] = passTPOneTagFilter;
    ele_passTPTwoTagFilter[nElectrons] = passTPTwoTagFilter;
    ele_passTPOneProbeFilter[nElectrons] = passTPOneProbeFilter;
    ele_passTPTwoProbeFilter[nElectrons] = passTPTwoProbeFilter;


    nElectrons++;
  }
  
  return true;
};

bool RazorTuplizer::fillTaus(){
  for (const reco::PFTau &tau : *taus) {
    if (tau.pt() < 20) continue;
    tauE[nTaus] = tau.energy();
    tauPt[nTaus] = tau.pt();
    tauEta[nTaus] = tau.eta();
    tauPhi[nTaus] = tau.phi();
    /*
    tau_IsLoose[nTaus] = bool(tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"));
    tau_IsMedium[nTaus] = bool(tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"));
    tau_IsTight[nTaus] = bool(tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"));
    tau_passEleVetoLoose[nTaus] = bool(tau.tauID("againstElectronLooseMVA6"));
    tau_passEleVetoMedium[nTaus] = bool(tau.tauID("againstElectronMediumMVA6"));
    tau_passEleVetoTight[nTaus] = bool(tau.tauID("againstElectronTightMVA6"));
    tau_passMuVetoLoose[nTaus] = bool(tau.tauID("againstMuonLoose3"));
    //tau_passMuVetoMedium[nTaus] = bool(tau.tauID("")); //doesn't exist anymore in miniAOD 2015 v2
    tau_passMuVetoTight[nTaus] = bool(tau.tauID("againstMuonTight3") );  
    tau_combinedIsoDeltaBetaCorr3Hits[nTaus] = tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
    tau_chargedIsoPtSum[nTaus] = tau.tauID("chargedIsoPtSum");
    tau_neutralIsoPtSum[nTaus] = tau.tauID("neutralIsoPtSum");
    tau_puCorrPtSum[nTaus] = tau.tauID("puCorrPtSum");
    tau_eleVetoMVA[nTaus] = tau.tauID("againstElectronMVA6Raw") ;
    tau_eleVetoCategory[nTaus] = tau.tauID("againstElectronMVA6category");
    //tau_muonVetoMVA[nTaus] = tau.tauID("againstMuonMVAraw"); //doesn't exist anymore in miniAOD 2015 v2
    tau_isoMVAnewDMwLT[nTaus] = tau.tauID("byIsolationMVArun2v1DBnewDMwLTraw");
    //tau_isoMVAnewDMwoLT[nTaus] = tau.tauID("byIsolationMVA3newDMwoLTraw") ; //doesn't exist anymore in miniAOD 2015 v2 

    tau_ID[nTaus] = 
      bool(tau.tauID("decayModeFinding")) +
      bool(tau.tauID("decayModeFindingNewDMs")) +
      bool(tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")) +
      bool(tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits")) +
      bool(tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits")) +
      bool(tau.tauID("againstElectronVLooseMVA6")) +
      bool(tau.tauID("againstElectronLooseMVA6")) +
      bool(tau.tauID("againstElectronMediumMVA6")) +
      bool(tau.tauID("againstElectronTightMVA6")) +
      bool(tau.tauID("againstElectronVTightMVA6")) +
      bool(tau.tauID("againstMuonLoose3")) +
      bool(tau.tauID("againstMuonTight3")) +
      bool(tau.tauID("byVLooseIsolationMVArun2v1DBnewDMwLT")) +
      bool(tau.tauID("byLooseIsolationMVArun2v1DBnewDMwLT")) +
      bool(tau.tauID("byMediumIsolationMVArun2v1DBnewDMwLT")) +
      bool(tau.tauID("byTightIsolationMVArun2v1DBnewDMwLT")) +
      bool(tau.tauID("byVTightIsolationMVArun2v1DBnewDMwLT")) +
      bool(tau.tauID("byVVTightIsolationMVArun2v1DBnewDMwLT"));
*/
    tau_leadCandPt[nTaus] = 0;
    tau_leadCandID[nTaus] = 0;
    tau_leadChargedHadrCandPt[nTaus] = 0;
    tau_leadChargedHadrCandID[nTaus] = 0;
  
  if (tau.leadPFCand().isNonnull()) {
      tau_leadCandPt[nTaus] = tau.leadPFCand()->pt();
      tau_leadCandID[nTaus] = tau.leadPFCand()->pdgId();
    }


    if (tau.leadPFChargedHadrCand().isNonnull()) { 
      tau_leadChargedHadrCandPt[nTaus] = tau.leadPFChargedHadrCand()->pt();
      tau_leadChargedHadrCandID[nTaus] = tau.leadPFChargedHadrCand()->pdgId();
    }

      
    nTaus++;
  }

  return true;
};

bool RazorTuplizer::fillIsoPFCandidates(){

//  for (const reco::PFCandidate &candidate : *pfCands) {

/*
    if (candidate.charge() != 0 && candidate.pt() > 5 && candidate.fromPV() == 3 ) {
      double tmpIsoPFNoPU = 0;
      double tmpIsoPFPU = 0;
      for (const reco::PFCandidate &isoCandidate : *pfCands) {	
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
	//isoPFCandidateD0[nIsoPFCandidates] = candidate.dxy();
	isoPFCandidatePdgId[nIsoPFCandidates] = candidate.pdgId();
	

	nIsoPFCandidates++;

	//For Debugging
	// cout << "\n";
	// cout << candidate.charge() << " " << candidate.pdgId() << " " << candidate.pt() << " " << candidate.eta() << " " << candidate.phi() 
	//      << " | " << candidate.fromPV() << " " << candidate.dz() ;
	// cout << " | " << tmpIsoPFNoPU << " " <<tmpIsoPFPU <<  " " <<  tmpIsoPFNoPU - 0.5*tmpIsoPFPU << " " << (tmpIsoPFNoPU - 0.5*tmpIsoPFPU)/candidate.pt()  ;
	// cout << "\n";
	for (const reco::PFCandidate &isoCandidate : *pfCands) {	
	  if ( (candidate.pdgId() != 1 && candidate.pdgId() != 2)
	       && deltaR(candidate.eta(), candidate.phi(), isoCandidate.eta(), isoCandidate.phi()) < 0.4
	       && !(candidate.eta() == isoCandidate.eta() && candidate.phi() == isoCandidate.phi())
	       ) {
	    //cout << "isoCandidate " << isoCandidate.pdgId() << " " << isoCandidate.pt() << " " << isoCandidate.eta() << " " << isoCandidate.phi() << " | " << deltaR(candidate.eta(), candidate.phi(), isoCandidate.eta(), isoCandidate.phi()) << " " << candidate.fromPV() << "\n";	  
	  }
	}

      } // if candidate passes isolation
    }*/
 // }

  return true;
}

bool RazorTuplizer::fillPhotons(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  noZS::EcalClusterLazyTools *lazyToolnoZS = new noZS::EcalClusterLazyTools(iEvent, iSetup, ebRecHitsToken_, eeRecHitsToken_);
    
  for (const reco::Photon &pho : *photons) {
    //if (pho.pt() < 20) continue;

    std::vector<float> vCov = lazyToolnoZS->localCovariances( *(pho.superCluster()->seed()) );

    //-------------------------------------------------
    //default photon 4-mometum already vertex corrected
    //-------------------------------------------------
    //phoE[nPhotons] = pho.getCorrectedEnergy(reco::Photon::P4type::ecal_standard);
    phoE[nPhotons]   = pho.energy();
    phoPt[nPhotons]  = pho.pt();
    phoEta[nPhotons] = pho.eta(); //correct this for the vertex
    phoPhi[nPhotons] = pho.phi(); //correct this for the vertex

    /*std::cout << "phoE: " << pho.energy() << " phoCorr En:" << pho.getCorrectedEnergy(reco::Photon::P4type::regression2) << " un: " 
	      << pho.getCorrectedEnergyError(reco::Photon::P4type::regression2) << " " 
	      << pho.getCorrectedEnergyError( pho.getCandidateP4type() ) << std::endl;
    */

    phoSigmaIetaIeta[nPhotons] = pho.sigmaIetaIeta();
    phoFull5x5SigmaIetaIeta[nPhotons] = pho.full5x5_sigmaIetaIeta();    

    //phoR9[nPhotons] = pho.r9();
    //Use the noZS version of this according to Emanuele
    phoR9[nPhotons] = pho.full5x5_r9();

    pho_HoverE[nPhotons] = pho.hadTowOverEm();
    pho_isConversion[nPhotons] = pho.hasConversionTracks();
    pho_passEleVeto[nPhotons] = !hasMatchedPromptElectron(pho.superCluster(),electrons, 
									   conversions, beamSpot->position());

    //**********************************************************
    // Fill default miniAOD isolation quantities
    //**********************************************************
    pho_pfIsoChargedHadronIso[nPhotons] = pho.chargedHadronIso();
    pho_pfIsoChargedHadronIsoWrongVtx[nPhotons] = pho.chargedHadronIsoWrongVtx();
    pho_pfIsoNeutralHadronIso[nPhotons] = pho.neutralHadronIso();
    pho_pfIsoPhotonIso[nPhotons] = pho.photonIso();
    pho_pfIsoModFrixione[nPhotons] = pho.getPflowIsolationVariables().modFrixione;
    pho_pfIsoSumPUPt[nPhotons] = pho.sumPUPt();
    
    //**********************************************************
    //Compute PF isolation
    //absolute uncorrected isolations with footprint removal
    //**********************************************************
    const float coneSizeDR = 0.3;
    const float dxyMax = 0.1;
    const float dzMax = 0.2;
    float chargedIsoSumAllVertices[MAX_NPV];
    for (int q=0;q<MAX_NPV;++q) chargedIsoSumAllVertices[q] = 0.0;
    float chargedIsoSum = 0;
    float neutralHadronIsoSum = 0;
    float photonIsoSum = 0;

    // First, find photon direction with respect to the good PV
    math::XYZVector photon_directionWrtVtx(pho.superCluster()->x() - myPV->x(),
					   pho.superCluster()->y() - myPV->y(),
					   pho.superCluster()->z() - myPV->z());
    // Loop over all PF candidates
    for (const reco::PFCandidate &candidate : *pfCands) {

     // Check if this candidate is within the isolation cone
      float dR=deltaR(photon_directionWrtVtx.Eta(),photon_directionWrtVtx.Phi(),
		      candidate.eta(), candidate.phi());
      if( dR > coneSizeDR ) continue;

      // Check if this candidate is not in the footprint
      /*
      bool inFootprint = false;      
      for (auto itr : pho.associatedPackedPFCandidates()) {	
	if ( &(*itr) == &candidate) {
	  inFootprint = true;
	}
      }     
      if( inFootprint ) continue;

	*/
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

	//For the primary vertex
	dz = candidate.trackRef()->dz(myPV->position());
	dxy =candidate.trackRef()->dxy(myPV->position());
	if (fabs(dz) <= dzMax && fabs(dxy) <= dxyMax) {
	  chargedIsoSum += candidate.pt();
	}

	//loop over all vertices
	for(int q = 0; q < nPVAll; q++){
	  if(!(vertices->at(q).isValid() && !vertices->at(q).isFake())) continue;

	  dz = candidate.trackRef()->dz(vertices->at(q).position());
	  dxy =candidate.trackRef()->dxy(vertices->at(q).position());
	  if (fabs(dz) > dzMax) continue;
	  if(fabs(dxy) > dxyMax) continue;
	  // The candidate is eligible, increment the isolation
	  chargedIsoSumAllVertices[q] += candidate.pt();
	}
      }
      if( thisCandidateType == reco::PFCandidate::h0 )
	neutralHadronIsoSum += candidate.pt();
      if( thisCandidateType == reco::PFCandidate::gamma )
	photonIsoSum += candidate.pt();
    }

    //fill the proper variables
    for(int q = 0; q < nPVAll; q++) {
      pho_sumChargedHadronPtAllVertices[nPhotons][q] = chargedIsoSumAllVertices[q];
    }
    pho_sumChargedHadronPt[nPhotons] = chargedIsoSum;
    pho_sumNeutralHadronEt[nPhotons] = neutralHadronIsoSum;
    pho_sumPhotonEt[nPhotons] = photonIsoSum;
    

    //*****************************************************************
    //Compute Worst Isolation Looping over all vertices
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
      for (const reco::PFCandidate &candidate : *pfCands) {
		
	//require that PFCandidate is a charged hadron
	const int pdgId = candidate.pdgId();
	if( abs(pdgId) != 211) continue;

	if (candidate.pt() < ptMin)
	  continue;
      
	float dxy = -999, dz = -999;
	dz = candidate.trackRef()->dz(myPV->position());
	dxy =candidate.trackRef()->dxy(myPV->position());
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
    
    pho_sumWorstVertexChargedHadronPt[nPhotons] = worstIsolation;

    //*****************************************************************
    //Photon ID MVA variable
    //*****************************************************************
    pho_IDMVA[nPhotons] = myPhotonMVA->mvaValue( pho,  *rhoAll, photonIsoSum, chargedIsoSum, worstIsolation,
						 lazyToolnoZS, false);
				       
    //pho_RegressionE[nPhotons] = pho.getCorrectedEnergy(reco::Photon::P4type::regression1);
    //pho_RegressionEUncertainty[nPhotons] = pho.getCorrectedEnergyError(reco::Photon::P4type::regression1);
    
    //---------------------
    //Use Latest Regression
    //---------------------
    pho_RegressionE[nPhotons]            = pho.getCorrectedEnergy( pho.getCandidateP4type() );
    pho_RegressionEUncertainty[nPhotons] = pho.getCorrectedEnergyError( pho.getCandidateP4type() );
    
    //default photon 4-momentum is already corrected.
    //compute photon corrected 4-mometum 
    /*
      TVector3 phoPos( pho.superCluster()->x(), pho.superCluster()->y(), pho.superCluster()->z() );
      TVector3 vtxPos( pvX, pvY, pvZ );
      TLorentzVector phoP4 = photonP4FromVtx( vtxPos, phoPos, pho_RegressionE[nPhotons] );
      std::cout << "etaDefault: " << phoEta[nPhotons] << " CP: " << phoP4.Eta() << " phiDefault: " 
      << phoPhi[nPhotons] << " CP: " << phoP4.Phi() << std::endl;
      phoEta[nPhotons] = phoP4.Eta();
      phoPhi[nPhotons] = phoP4.Phi();
    */
    
    //-----------------------
    // super cluster position
    //-----------------------  
    pho_superClusterEnergy[nPhotons] = pho.superCluster()->energy();
    pho_superClusterRawEnergy[nPhotons] = pho.superCluster()->rawEnergy();
    pho_superClusterEta[nPhotons]    = pho.superCluster()->eta();
    pho_superClusterPhi[nPhotons]    = pho.superCluster()->phi();
    pho_superClusterX[nPhotons]      = pho.superCluster()->x();
    pho_superClusterY[nPhotons]      = pho.superCluster()->y();
    pho_superClusterZ[nPhotons]      = pho.superCluster()->z();
    pho_hasPixelSeed[nPhotons]       = pho.hasPixelSeed();

    pho_superClusterSeedX[nPhotons]      = pho.superCluster()->seed()->x();
    pho_superClusterSeedY[nPhotons]      = pho.superCluster()->seed()->y();
    pho_superClusterSeedZ[nPhotons]      = pho.superCluster()->seed()->z();
  
    for (const reco::PFCluster &pfcluster : *pfClusters) {
	if(pfcluster.seed() == pho.superCluster()->seed()->seed())
	{
	pho_superClusterSeedT[nPhotons] = pfcluster.time();
	}
    } 

    //Detector DetId_this = pho.superCluster()->seed()->seed();
 
    //*************************************************
    //Trigger Object Matching
    //*************************************************
    /*
    for (pat::TriggerObjectStandAlone trigObject : *triggerObjects) {
      if (deltaR(trigObject.eta(), trigObject.phi(),pho.eta(),pho.phi()) > 0.3) continue;
    
      //check all filters
      for ( int q=0; q<MAX_PhotonHLTFilters;q++) {
	if (trigObject.hasFilterLabel(photonHLTFilterNames[q].c_str())) pho_passHLTFilter[nPhotons][q] = true;
      }
    }
    */
    //conversion matching for beamspot pointing
    const reco::Conversion *convmatch = 0;
    double drmin = std::numeric_limits<double>::max();
    //double leg conversions
    for (const reco::Conversion &conv : *conversions) {
      if (conv.refittedPairMomentum().rho()<10.) continue;
      if (!conv.conversionVertex().isValid()) continue;
      if (TMath::Prob(conv.conversionVertex().chi2(),  conv.conversionVertex().ndof())<1e-6) continue;
      
      math::XYZVector mom(conv.refittedPairMomentum());      
      math::XYZPoint scpos(pho.superCluster()->position());
      math::XYZPoint cvtx(conv.conversionVertex().position());
      math::XYZVector cscvector = scpos - cvtx;
      
      double dr = reco::deltaR(mom,cscvector);
      
      if (dr<drmin && dr<0.1) {
        drmin = dr;
        convmatch = &conv;
      }
    }
    if (!convmatch) {
      drmin = std::numeric_limits<double>::max();
      //single leg conversions
      for (const reco::Conversion &conv : *singleLegConversions) {      
        math::XYZVector mom(conv.tracksPin()[0]);      
        math::XYZPoint scpos(pho.superCluster()->position());
        math::XYZPoint cvtx(conv.conversionVertex().position());
        math::XYZVector cscvector = scpos - cvtx;
        
        double dr = reco::deltaR(mom,cscvector);
        
        if (dr<drmin && dr<0.1) {
          drmin = dr;
          convmatch = &conv;
        }
      }
    }
    
    //matched conversion, compute conversion type
    //and extrapolation to beamline
    //*FIXME* Both of these additional two requirements are inconsistent and make the conversion
    //selection depend on poorly defined criteria, but we keep them for sync purposes
    //if (convmatch && pho.hasConversionTracks() && conversions->size()>0) {
    if (convmatch){// && pho.hasConversionTracks() && conversions->size()>0) {
      int ntracks = convmatch->nTracks();
      
      math::XYZVector mom(ntracks==2 ? convmatch->refittedPairMomentum() : convmatch->tracksPin()[0]);      
      math::XYZPoint scpos(pho.superCluster()->position());
      math::XYZPoint cvtx(convmatch->conversionVertex().position());
      math::XYZVector cscvector = scpos - cvtx;

      double z = cvtx.z();
      double rho = cvtx.rho();

      int legtype = ntracks==2 ? 0 : 1;
      int dettype = pho.isEB() ? 0 : 1;
      int postype =0;
      
      if (pho.isEB()) {
        if (rho<15.) {
          postype = 0;
        }
        else if (rho>=15. && rho<60.) {
          postype = 1;
        }
        else {
          postype = 2;
        }
      }
      else {
        if (std::abs(z) < 50.) {
          postype = 0;
        }
        else if (std::abs(z) >= 50. && std::abs(z) < 100.) {
          postype = 1;
        }
        else {
          postype = 2;
        }
      }
      
      pho_convType[nPhotons] = legtype + 2*dettype + 4*postype;
      pho_convTrkZ[nPhotons] = cvtx.z() - ((cvtx.x()-beamSpot->x0())*mom.x()+(cvtx.y()-beamSpot->y0())*mom.y())/mom.rho() * mom.z()/mom.rho();
      pho_convTrkClusZ[nPhotons] = cvtx.z() - ((cvtx.x()-beamSpot->x0())*cscvector.x()+(cvtx.y()-beamSpot->y0())*cscvector.y())/cscvector.rho() * cscvector.z()/cscvector.rho();
    }

    nPhotons++;
  }
  
  double pho_vtxSumPxD[OBJECTARRAYSIZE][MAX_NPV];
  double pho_vtxSumPyD[OBJECTARRAYSIZE][MAX_NPV];
  

  for (int ipho = 0; ipho<nPhotons; ++ipho) {
    for (int ipv = 0; ipv<nPVAll; ++ipv) {
      pho_vtxSumPxD[ipho][ipv] = 0.;
      pho_vtxSumPyD[ipho][ipv] = 0.;

	}
  }
  
  
  //fill information on tracks to exclude around photons for vertex selection purposes
  for (const reco::PFCandidate &pfcand : *pfCands) {
    if (pfcand.charge()==0) continue;
    double mindz = std::numeric_limits<double>::max();
    int ipvmin = -1;
    for (int ipv = 0; ipv < nPVAll; ++ipv) {
      const reco::Vertex &vtx = vertices->at(ipv);
      //double dz = std::abs(pfcand.dz(vtx.position()));
      double dz = std::abs(pfcand.vz()-vtx.z());
      if (dz<mindz) {
        mindz = dz;
        ipvmin = ipv;
      }
    }
    
    if (mindz<0.2 && ipvmin>=0 && ipvmin<MAX_NPV) {
      const reco::Vertex &vtx = vertices->at(ipvmin);
      for (int ipho = 0; ipho < nPhotons; ++ipho) {
        const reco::Photon &pho = photons->at(ipho);
        math::XYZVector phodir(pho.superCluster()->x()-vtx.x(),pho.superCluster()->y()-vtx.y(),pho.superCluster()->z()-vtx.z());
        double dr = reco::deltaR(phodir, pfcand);
        if (dr<0.05) {
          pho_vtxSumPxD[ipho][ipvmin] += pfcand.px();
          pho_vtxSumPyD[ipho][ipvmin] += pfcand.py();
        }
	//add addition dt cut here:
	//
      }
    }
  }
  
  for (int ipho = 0; ipho<nPhotons; ++ipho) {
    for (int ipv = 0; ipv<nPVAll; ++ipv) {
      pho_vtxSumPx[ipho][ipv] = pho_vtxSumPxD[ipho][ipv];
      pho_vtxSumPy[ipho][ipv] = pho_vtxSumPyD[ipho][ipv];

    }
  }
  
  
  delete lazyToolnoZS;
  return true;
};

bool RazorTuplizer::fillJets(){
  for (const reco::PFJet &j : *jets) {
    if (j.pt() < 10) continue;

    jetE[nJets] = j.energy();
    jetPt[nJets] = j.pt();
    jetEta[nJets] = j.eta();
    jetPhi[nJets] = j.phi();

    jetMass[nJets] = j.mass();
    
   // jetE[nJets] = j.correctedP4(0).E();
   // jetPt[nJets] = j.correctedP4(0).Pt();
   // jetEta[nJets] = j.correctedP4(0).Eta();
   // jetPhi[nJets] = j.correctedP4(0).Phi();
   
    //jetCSV[nJets] = j.bDiscriminator("pfCombinedSecondaryVertexBJetTags");
    //jetCISV[nJets] = j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    //jetMass[nJets] = j.correctedP4(0).M();
    jetJetArea[nJets] = j.jetArea();
    jetPileupE[nJets] = j.pileup();
    //jetPileupId[nJets] = j.userFloat("pileupJetId:fullDiscriminant");
    jetPileupIdFlag[nJets] = 0;
    jetPassIDLoose[nJets] = passJetID(&j, 0);
    jetPassIDTight[nJets] = passJetID(&j, 1);
    jetPassMuFrac[nJets]  = ( j.muonEnergyFraction() < 0.80 );
    jetPassEleFrac[nJets]  = ( j.electronEnergyFraction() < 0.90 );
 /*
    if (useGen_) {
      jetPartonFlavor[nJets] = j.partonFlavour();
      jetHadronFlavor[nJets] = j.hadronFlavour();
    } else 
    */
    {
      jetPartonFlavor[nJets] = -999;
      jetHadronFlavor[nJets] = -999;
    }
 

    //extra jet information (may be only needed for debugging)
    jetChargedEMEnergyFraction[nJets] = j.chargedEmEnergyFraction();
    jetNeutralEMEnergyFraction[nJets] = j.neutralEmEnergyFraction();
    jetChargedHadronEnergyFraction[nJets] = j.chargedHadronEnergyFraction();
    jetNeutralHadronEnergyFraction[nJets] = j.neutralHadronEnergyFraction();
    jetMuonEnergyFraction[nJets] =  j.muonEnergyFraction();
    jetHOEnergyFraction[nJets] =  j.hoEnergyFraction();
    jetHFHadronEnergyFraction[nJets] =  j.HFHadronEnergyFraction();
    jetHFEMEnergyFraction[nJets] =  j.HFEMEnergyFraction();
	
    //save muon vector for JEC's
    reco::Candidate::LorentzVector AllMuonP4; AllMuonP4.SetPxPyPzE(0,0,0,0);
    const std::vector<reco::CandidatePtr> & cands = j.daughterPtrVector();
    for ( std::vector<reco::CandidatePtr>::const_iterator cand = cands.begin();
	  cand != cands.end(); ++cand ) {
      const reco::PFCandidate *pfcand = dynamic_cast<const reco::PFCandidate *>(cand->get());
      const reco::Candidate *mu = (pfcand != 0 ? ( pfcand->muonRef().isNonnull() ? pfcand->muonRef().get() : 0) : cand->get());
      if ( mu != 0 && (mu->isGlobalMuon() || mu->isStandAloneMuon()) ) {
	reco::Candidate::LorentzVector muonP4 = (*cand)->p4();
	AllMuonP4 = AllMuonP4 + muonP4;
      }
    }
    jetAllMuonPt[nJets] = AllMuonP4.Pt();
    jetAllMuonEta[nJets] = AllMuonP4.Eta();
    jetAllMuonPhi[nJets] = AllMuonP4.Phi();
    jetAllMuonM[nJets] = AllMuonP4.M();
   
    nJets++;
  }

  return true;
};

bool RazorTuplizer::fillJetsAK8(){
  for (const reco::PFJet &j : *jetsAK8) {
    fatJetE[nFatJets] = j.energy();
    fatJetPt[nFatJets] = j.pt();
    fatJetEta[nFatJets] = j.eta();
    fatJetPhi[nFatJets] = j.phi();

   // fatJetE[nFatJets] = j.correctedP4(0).E();
   // fatJetPt[nFatJets] = j.correctedP4(0).Pt();
   // fatJetEta[nFatJets] = j.correctedP4(0).Eta();
   // fatJetPhi[nFatJets] = j.correctedP4(0).Phi();
    //fatJetPrunedM[nFatJets] = (float) j.userFloat("ak8PFJetsCHSPrunedLinks");
    //fatJetTrimmedM[nFatJets] = (float) j.userFloat("ak8PFJetsCHSTrimmedLinks");
    //fatJetFilteredM[nFatJets] = (float) j.userFloat("ak8PFJetsCHSFilteredLinks");
    //fatJetPrunedM[nFatJets] = (float) j.userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSPrunedMass");                                                     
    //fatJetTrimmedM[nFatJets] = (float) j.userFloat("ak8PFJetsCHSTrimmedMass");
    //fatJetFilteredM[nFatJets] = (float) j.userFloat("ak8PFJetsCHSFilteredMass");  
    //fatJetTau1[nFatJets] =  (float) j.userFloat("NjettinessAK8Puppi:tau1");
    //fatJetTau2[nFatJets] =  (float) j.userFloat("NjettinessAK8Puppi:tau2");
    //fatJetTau3[nFatJets] =  (float) j.userFloat("NjettinessAK8Puppi:tau3");
    nFatJets++;
  }

  return true;
};

bool RazorTuplizer::fillMet(const edm::Event& iEvent){
  const reco::PFMET &Met = mets->front();
  
  //metPt = Met.uncorPt();
  //metPhi = Met.uncorPhi();
  sumMET = Met.sumEt();
  metType0Pt = 0;
  metType0Phi = 0;
  metType1Pt = Met.pt();
  metType1Phi = Met.phi();
  metType0Plus1Pt = 0;
  metType0Plus1Phi = 0;
  //metCaloPt = Met.caloMETPt();
  //metCaloPhi = Met.caloMETPhi();
/*
  if(!isData_)
    {
      metType1PtJetResUp           = Met.shiftedPt(reco::PFMET::METUncertainty::JetResUp, reco::PFMET::METCorrectionLevel::Type1);
      metType1PtJetResDown         = Met.shiftedPt(reco::PFMET::METUncertainty::JetResDown, reco::PFMET::METCorrectionLevel::Type1);
      metType1PtJetEnUp            = Met.shiftedPt(reco::PFMET::METUncertainty::JetEnUp, reco::PFMET::METCorrectionLevel::Type1);
      metType1PtJetEnDown          = Met.shiftedPt(reco::PFMET::METUncertainty::JetEnDown, reco::PFMET::METCorrectionLevel::Type1);
      metType1PtMuonEnUp           = Met.shiftedPt(reco::PFMET::METUncertainty::MuonEnUp, reco::PFMET::METCorrectionLevel::Type1);
      metType1PtMuonEnDown         = Met.shiftedPt(reco::PFMET::METUncertainty::MuonEnDown, reco::PFMET::METCorrectionLevel::Type1);
      metType1PtElectronEnUp       = Met.shiftedPt(reco::PFMET::METUncertainty::ElectronEnUp, reco::PFMET::METCorrectionLevel::Type1);
      metType1PtElectronEnDown     = Met.shiftedPt(reco::PFMET::METUncertainty::ElectronEnDown, reco::PFMET::METCorrectionLevel::Type1);
      metType1PtTauEnUp	           = Met.shiftedPt(reco::PFMET::METUncertainty::TauEnUp, reco::PFMET::METCorrectionLevel::Type1);
      metType1PtTauEnDown          = Met.shiftedPt(reco::PFMET::METUncertainty::TauEnDown, reco::PFMET::METCorrectionLevel::Type1);
      metType1PtUnclusteredEnUp    = Met.shiftedPt(reco::PFMET::METUncertainty::UnclusteredEnUp, reco::PFMET::METCorrectionLevel::Type1);
      metType1PtUnclusteredEnDown  = Met.shiftedPt(reco::PFMET::METUncertainty::UnclusteredEnDown, reco::PFMET::METCorrectionLevel::Type1);
      metType1PtPhotonEnUp         = Met.shiftedPt(reco::PFMET::METUncertainty::PhotonEnUp, reco::PFMET::METCorrectionLevel::Type1);
      metType1PtPhotonEnDown       = Met.shiftedPt(reco::PFMET::METUncertainty::PhotonEnDown, reco::PFMET::METCorrectionLevel::Type1);
      // metType1PtMETUncertaintySize = Met.shiftedPt(reco::PFMET::METUncertainty::METUncertaintySize, reco::PFMET::METCorrectionLevel::Type1);
      // metType1PtJetResUpSmear     = Met.shiftedPt(reco::PFMET::METUncertainty::JetResUpSmear, reco::PFMET::METCorrectionLevel::Type1);	      
      // metType1PtJetResDownSmear   = Met.shiftedPt(reco::PFMET::METUncertainty::JetResDownSmear, reco::PFMET::METCorrectionLevel::Type1);	      
      // metType1PtMETFullUncertaintySize = Met.shiftedPt(reco::PFMET::METUncertainty::METFullUncertaintySize, reco::PFMET::METCorrectionLevel::Type1);	     
      
      metType1PhiJetResUp          = Met.shiftedPhi(reco::PFMET::METUncertainty::JetResUp, reco::PFMET::METCorrectionLevel::Type1);
      metType1PhiJetResDown        = Met.shiftedPhi(reco::PFMET::METUncertainty::JetResDown, reco::PFMET::METCorrectionLevel::Type1);
      metType1PhiJetEnUp           = Met.shiftedPhi(reco::PFMET::METUncertainty::JetEnUp, reco::PFMET::METCorrectionLevel::Type1);	      
      metType1PhiJetEnDown         = Met.shiftedPhi(reco::PFMET::METUncertainty::JetEnDown, reco::PFMET::METCorrectionLevel::Type1);	      
      metType1PhiMuonEnUp          = Met.shiftedPhi(reco::PFMET::METUncertainty::MuonEnUp, reco::PFMET::METCorrectionLevel::Type1);	      
      metType1PhiMuonEnDown        = Met.shiftedPhi(reco::PFMET::METUncertainty::MuonEnDown, reco::PFMET::METCorrectionLevel::Type1);	 
      metType1PhiElectronEnUp      = Met.shiftedPhi(reco::PFMET::METUncertainty::ElectronEnUp, reco::PFMET::METCorrectionLevel::Type1);	      
      metType1PhiElectronEnDown    = Met.shiftedPhi(reco::PFMET::METUncertainty::ElectronEnDown, reco::PFMET::METCorrectionLevel::Type1);	      
      metType1PhiTauEnUp           = Met.shiftedPhi(reco::PFMET::METUncertainty::TauEnUp, reco::PFMET::METCorrectionLevel::Type1);	      
      metType1PhiTauEnDown         = Met.shiftedPhi(reco::PFMET::METUncertainty::TauEnDown, reco::PFMET::METCorrectionLevel::Type1);	      
      metType1PhiUnclusteredEnUp   = Met.shiftedPhi(reco::PFMET::METUncertainty::UnclusteredEnUp, reco::PFMET::METCorrectionLevel::Type1);	      
      metType1PhiUnclusteredEnDown = Met.shiftedPhi(reco::PFMET::METUncertainty::UnclusteredEnDown, reco::PFMET::METCorrectionLevel::Type1);	      
      metType1PhiPhotonEnUp        = Met.shiftedPhi(reco::PFMET::METUncertainty::PhotonEnUp, reco::PFMET::METCorrectionLevel::Type1);	      
      metType1PhiPhotonEnDown      = Met.shiftedPhi(reco::PFMET::METUncertainty::PhotonEnDown, reco::PFMET::METCorrectionLevel::Type1);	      
      // metType1PhiMETUncertaintySize = Met.shiftedPhi(reco::PFMET::METUncertainty::METUncertaintySize, reco::PFMET::METCorrectionLevel::Type1);	      
      // metType1PhiJetResUpSmear     = Met.shiftedPhi(reco::PFMET::METUncertainty::JetResUpSmear, reco::PFMET::METCorrectionLevel::Type1);	      
      // metType1PhiJetResDownSmear   = Met.shiftedPhi(reco::PFMET::METUncertainty::JetResDownSmear, reco::PFMET::METCorrectionLevel::Type1);	      
      // metType1PhiMETFullUncertaintySize = Met.shiftedPhi(reco::PFMET::METUncertainty::METFullUncertaintySize, reco::PFMET::METCorrectionLevel::Type1);	     
    }
 */
 
  const reco::PFMET &MetPuppi = metsPuppi->front();
  //const reco::PFMET &MetNoHF = metsNoHF->front();
  metPuppiPt = MetPuppi.pt();
  metPuppiPhi = MetPuppi.phi();
  //metNoHFPt = MetNoHF.pt();
  //metNoHFPhi = MetNoHF.phi();
  
  //MET filters
  if (!isFastsim_) {
    const edm::TriggerNames &metNames = iEvent.triggerNames(*metFilterBits);
    
    //*******************************************************************************
    //For Debug printout
    //*******************************************************************************
    // for (unsigned int i = 0, n = metFilterBits->size(); i < n; ++i) {
    // 	std::cout << "MET Filter " << metNames.triggerName(i).c_str() << "\n";
    // }
    
    for(unsigned int i = 0, n = metFilterBits->size(); i < n; ++i){
      if(strcmp(metNames.triggerName(i).c_str(), "Flag_trackingFailureFilter") == 0)
	Flag_trackingFailureFilter = metFilterBits->accept(i);
      else if(strcmp(metNames.triggerName(i).c_str(), "Flag_goodVertices") == 0)
	Flag_goodVertices = metFilterBits->accept(i);
      else if(strcmp(metNames.triggerName(i).c_str(), "Flag_globalTightHalo2016Filter") == 0)
	Flag_CSCTightHaloFilter = metFilterBits->accept(i);
      else if(strcmp(metNames.triggerName(i).c_str(), "Flag_trkPOGFilters") == 0)
	Flag_trkPOGFilters = metFilterBits->accept(i);
      else if(strcmp(metNames.triggerName(i).c_str(), "Flag_trkPOG_logErrorTooManyClusters") == 0)
	Flag_trkPOG_logErrorTooManyClusters = metFilterBits->accept(i);
      else if(strcmp(metNames.triggerName(i).c_str(), "Flag_EcalDeadCellTriggerPrimitiveFilter") == 0)
	Flag_EcalDeadCellTriggerPrimitiveFilter = metFilterBits->accept(i);
      else if(strcmp(metNames.triggerName(i).c_str(), "Flag_EcalDeadCellBoundaryEnergyFilter") == 0)
	Flag_EcalDeadCellBoundaryEnergyFilter = metFilterBits->accept(i);
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
       else if(strcmp(metNames.triggerName(i).c_str(), "Flag_HBHENoiseIsoFilter") == 0)
         Flag_HBHEIsoNoiseFilter = metFilterBits->accept(i);
      else if(strcmp(metNames.triggerName(i).c_str(), "Flag_trkPOG_toomanystripclus53X") == 0)
	Flag_trkPOG_toomanystripclus53X = metFilterBits->accept(i);
      else if(strcmp(metNames.triggerName(i).c_str(), "Flag_hcalLaserEventFilter") == 0)
	Flag_hcalLaserEventFilter = metFilterBits->accept(i);     
    } //loop over met filters

    //use custom hbhefilter, because miniAOD filters are problematic.
    //Flag_HBHENoiseFilter = *hbheNoiseFilter;
    //Flag_HBHETightNoiseFilter = *hbheTightNoiseFilter;
    //Flag_HBHEIsoNoiseFilter = *hbheIsoNoiseFilter;
    //Flag_badChargedCandidateFilter = *badChargedCandidateFilter;
    //Flag_badMuonFilter = *badMuonFilter;
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

    //const reco::PFMET &Met = mets->front();
    //genMetPt = Met.genMET()->pt();
    //genMetPhi = Met.genMET()->phi();

    bool foundGenVertex = false;
    for(size_t i=0; i<genParticles->size();i++){
      if (!foundGenVertex) {
	for (unsigned int j=0; j<(*genParticles)[i].numberOfDaughters(); ++j) {
	  const reco::Candidate *dau = (*genParticles)[i].daughter(j);
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
     
    //save genVertex time info
    const HepMC::GenEvent *this_mc = hepMC->GetEvent();

    auto origin_vtx = (*this_mc->vertices_begin())->position();
    //xyz0Ptr->SetXYZ(origin.x() * mmToCm, origin.y() * mmToCm, origin.z() * mmToCm);
    //*t0Ptr = origin.t() * mmToNs;
    genVertexT = origin_vtx.t();
 
    genWeight = genInfo->weight();
    genSignalProcessID = genInfo->signalProcessID();
    genQScale = genInfo->qScale();
    genAlphaQCD = genInfo->alphaQCD();
    genAlphaQED = genInfo->alphaQED();
    

    if (isFastsim_) {

      //get lhe weights for systematic uncertainties:      
      double nomlheweight = genInfo->weights()[0];
      
      //fill scale variation weights
      if (genInfo->weights().size()>=10) {  
	for (unsigned int iwgt=1; iwgt<10; ++iwgt) {
	  //normalize to 
	  double wgtval = genInfo->weights()[iwgt]*genWeight/genInfo->weights()[1];
	  scaleWeights->push_back(wgtval);
	}
      }
   
      //fill pdf variation weights
      if (firstPdfWeight>=0 && lastPdfWeight>=0 && lastPdfWeight<int(genInfo->weights().size()) && (lastPdfWeight-firstPdfWeight+1)==100) {
        
	//fill pdf variation weights after converting with mc2hessian transformation
	std::array<double, 100> inpdfweights;
	for (int iwgt=firstPdfWeight, ipdf=0; iwgt<=lastPdfWeight; ++iwgt, ++ipdf) {
	  inpdfweights[ipdf] = genInfo->weights()[iwgt]/genInfo->weights()[firstPdfWeight-1];
	}
        
	std::array<double, 60> outpdfweights;
	pdfweightshelper.DoMC2Hessian(inpdfweights.data(),outpdfweights.data());
        
	for (unsigned int iwgt=0; iwgt<60; ++iwgt) {
	  double wgtval = outpdfweights[iwgt]*genWeight;
	  pdfWeights->push_back(wgtval);
	}       
              
	//fill alpha_s variation weights
	if (firstAlphasWeight>=0 && lastAlphasWeight>=0 && lastAlphasWeight<int(genInfo->weights().size())) {
	  for (int iwgt = firstAlphasWeight; iwgt<=lastAlphasWeight; ++iwgt) {
	    double wgtval = genInfo->weights()[iwgt]*genWeight/nomlheweight;
	    alphasWeights->push_back(wgtval);
	  }
	}        
        
      }
    } else {

/*
      if (lheInfo.isValid() && lheInfo->weights().size()>0) {
	
	double nomlheweight = lheInfo->weights()[0].wgt;
	
	//fill scale variation weights
	if (lheInfo->weights().size()>=9) {  
	  for (unsigned int iwgt=0; iwgt<9; ++iwgt) {
	    double wgtval = lheInfo->weights()[iwgt].wgt*genWeight/nomlheweight;
	    scaleWeights->push_back(wgtval);
	  }
	}
	
	//fill pdf variation weights
	if (firstPdfWeight>=0 && lastPdfWeight>=0 && lastPdfWeight<int(lheInfo->weights().size()) && (lastPdfWeight-firstPdfWeight+1)==100) {
	  
	  //fill pdf variation weights after converting with mc2hessian transformation
	  std::array<double, 100> inpdfweights;
	  for (int iwgt=firstPdfWeight, ipdf=0; iwgt<=lastPdfWeight; ++iwgt, ++ipdf) {
	    inpdfweights[ipdf] = lheInfo->weights()[iwgt].wgt/nomlheweight;
	  }
	  
	  std::array<double, 60> outpdfweights;
	  pdfweightshelper.DoMC2Hessian(inpdfweights.data(),outpdfweights.data());
	  
	  for (unsigned int iwgt=0; iwgt<60; ++iwgt) {
	    double wgtval = outpdfweights[iwgt]*genWeight;
	    pdfWeights->push_back(wgtval);
	  }       
	  
	  //fill alpha_s variation weights
	  if (firstAlphasWeight>=0 && lastAlphasWeight>=0 && lastAlphasWeight<int(lheInfo->weights().size())) {
	    for (int iwgt = firstAlphasWeight; iwgt<=lastAlphasWeight; ++iwgt) {
	      double wgtval = lheInfo->weights()[iwgt].wgt*genWeight/nomlheweight;
	      alphasWeights->push_back(wgtval);
	    }
	  }      	
	}
      }
	*/  
    }
   
    //fill sum of weights histograms
    sumWeights->Fill(0.,genWeight);
    
    for (unsigned int iwgt=0; iwgt<scaleWeights->size(); ++iwgt) {
      sumScaleWeights->Fill(double(iwgt),(*scaleWeights)[iwgt]);
    }
    for (unsigned int iwgt=0; iwgt<pdfWeights->size(); ++iwgt) {
      sumPdfWeights->Fill(double(iwgt),(*pdfWeights)[iwgt]);     
    }
    for (unsigned int iwgt=0; iwgt<alphasWeights->size(); ++iwgt) {
      sumAlphasWeights->Fill(double(iwgt),(*alphasWeights)[iwgt]);
    }        

    return true;
}

bool RazorTuplizer::fillRazor(){ 
  //get good jets for razor calculation
  vector<TLorentzVector> goodJets;
  for (const reco::PFJet &j : *jets) {
    if (j.pt() < 40) continue;
    if (fabs(j.eta()) > 3.0) continue;
    //add to goodJets vector
    TLorentzVector newJet(j.px(), j.py(), j.pz(), j.energy());
    goodJets.push_back(newJet);
  }
  //MET
  const reco::PFMET &Met = mets->front();
  
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
    for (const reco::PFJet &j : *jetsAK8) {
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
  for(size_t i=0; i<genParticles->size();i++){
    if(
       (abs((*genParticles)[i].pdgId()) >= 1 && abs((*genParticles)[i].pdgId()) <= 6 
    	&& ( (*genParticles)[i].status() < 30 	     
	     )
	)
       || (abs((*genParticles)[i].pdgId()) >= 11 && abs((*genParticles)[i].pdgId()) <= 16)
       || (abs((*genParticles)[i].pdgId()) == 21 
    	   && (*genParticles)[i].status() < 30
    	   )
       || (abs((*genParticles)[i].pdgId()) >= 22 && abs((*genParticles)[i].pdgId()) <= 25
    	   && ( (*genParticles)[i].status() < 30
		)
	   )
       || (abs((*genParticles)[i].pdgId()) >= 32 && abs((*genParticles)[i].pdgId()) <= 42)
       || (abs((*genParticles)[i].pdgId()) >= 1000001 && abs((*genParticles)[i].pdgId()) <= 1000039)
       ){
      prunedV.push_back(&(*genParticles)[i]);
    }    
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
  // std::cout << "\n === TRIGGER PATHS === " << std::endl;

  //********************************************************************
  //Option to save all HLT path names in the ntuple per event
  //Expensive option in terms of ntuple size
  //********************************************************************
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
	//triggerHLTPrescale[j] = triggerPrescales->getPrescaleForIndex(i);
      }
    }
  }

  // //********************************************************************
  // // Print Trigger Objects
  // //********************************************************************  
  // for (pat::TriggerObjectStandAlone trigObject : *triggerObjects) {
  //   cout << "triggerObj: " << trigObject.pt() << " " << trigObject.eta() << " " << trigObject.phi() << "\n";
  //   for(int j=0; j<int(trigObject.filterLabels().size());j++) {
  //     //trigObject.unpackPathNames(names);
  //     //cout << "filter: " << (trigObject.pathNames())[j] << " " << (trigObject.filterLabels())[j] << "\n";
  //     cout << "filter: " << (trigObject.filterLabels())[j] << "\n";
  //   }    
  // }


  return true;
}

//------ Method called for each run ------//

void RazorTuplizer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
  
  //read LHE header if present and determine which weights to read for pdf and alphas uncertainties based on the 
  //central pdf set used
  //This is semi-hardcoded for now to work with current centrally produced samples
  //covering SUSY signal samples, LO madgraph, NLO madgraph_aMC@NLO, and NLO powheg
  //generated with nnpdf30
  //More robust selection will require some parsing of <initrwgt> block
  
  if (useGen_) {
  
    if (isFastsim_) {

      //hardcode this for now. All SUSY signals use 5 flavor scheme, LO madgraph
      //to do it properly, we would need to parse the LHE header
      int pdfidx = 263000; 
            
      if (pdfidx == 263000) {
	//NNPDF30_lo_as_0130 (nf5) for LO madgraph samples and SUSY signals
	pdfweightshelper.Init(100,60,edm::FileInPath("SUSYBSMAnalysis/RazorTuplizer/data/NNPDF30_lo_as_0130_hessian_60.csv"));
	firstPdfWeight = 11;
	lastPdfWeight = 110;
	firstAlphasWeight = -1;
	lastAlphasWeight = -1;      
      }
      else if (pdfidx == 263400) {
	//NNPdf30_lo_as_0130_nf4 for LO madgraph samples
	pdfweightshelper.Init(100,60,edm::FileInPath("SUSYBSMAnalysis/RazorTuplizer/data/NNPDF30_lo_as_0130_nf_4_hessian_60.csv"));
	firstPdfWeight = 112;
	lastPdfWeight = 211;
	firstAlphasWeight = -1;
	lastAlphasWeight = -1;            
      }
      else if (pdfidx == 260000 || pdfidx == -1) {
	//NNPdf30_nlo_as_0118 (nf5) for NLO powheg samples
	//(work around apparent bug in current powheg samples by catching "-1" as well)
	pdfweightshelper.Init(100,60,edm::FileInPath("SUSYBSMAnalysis/RazorTuplizer/data/NNPDF30_nlo_as_0118_hessian_60.csv"));
	firstPdfWeight = 10;
	lastPdfWeight = 109;
	firstAlphasWeight = 110;
	lastAlphasWeight = 111; 
      }
      else if (pdfidx == 292200) {
	//NNPdf30_nlo_as_0118 (nf5) with built-in alphas variations for NLO aMC@NLO samples
	pdfweightshelper.Init(100,60,edm::FileInPath("SUSYBSMAnalysis/RazorTuplizer/data/NNPDF30_nlo_as_0118_hessian_60.csv"));
	firstPdfWeight = 10;
	lastPdfWeight = 109;
	firstAlphasWeight = 110;
	lastAlphasWeight = 111; 
      }   
      else if (pdfidx == 292000) {
	//NNPdf30_nlo_as_0118_nf4 with built-in alphas variations for NLO aMC@NLO samples
	pdfweightshelper.Init(100,60,edm::FileInPath("SUSYBSMAnalysis/RazorTuplizer/data/NNPDF30_nlo_as_0118_nf_4_hessian_60.csv"));
	firstPdfWeight = 10;
	lastPdfWeight = 109;
	firstAlphasWeight = 110;
	lastAlphasWeight = 111;
      }
      else {
	firstPdfWeight = -1;
	lastPdfWeight = -1;
	firstAlphasWeight = -1;
	lastAlphasWeight = -1;
      }    
    } else {
//      edm::Handle<LHERunInfoProduct> lheRunInfo;
//      iRun.getByLabel(lheRunInfoTag_,lheRunInfo);
   /* 
      if (lheRunInfo.isValid()) {
	int pdfidx = lheRunInfo->heprup().PDFSUP.first;
      
	if (pdfidx == 263000) {
	  //NNPDF30_lo_as_0130 (nf5) for LO madgraph samples and SUSY signals
	  pdfweightshelper.Init(100,60,edm::FileInPath("SUSYBSMAnalysis/RazorTuplizer/data/NNPDF30_lo_as_0130_hessian_60.csv"));
	  firstPdfWeight = 10;
	  lastPdfWeight = 109;
	  firstAlphasWeight = -1;
	  lastAlphasWeight = -1;      
	}
	else if (pdfidx == 263400) {
	  //NNPdf30_lo_as_0130_nf4 for LO madgraph samples
	  pdfweightshelper.Init(100,60,edm::FileInPath("SUSYBSMAnalysis/RazorTuplizer/data/NNPDF30_lo_as_0130_nf_4_hessian_60.csv"));
	  firstPdfWeight = 111;
	  lastPdfWeight = 210;
	  firstAlphasWeight = -1;
	  lastAlphasWeight = -1;            
	}
	else if (pdfidx == 260000 || pdfidx == -1) {
	  //NNPdf30_nlo_as_0118 (nf5) for NLO powheg samples
	  //(work around apparent bug in current powheg samples by catching "-1" as well)
	  pdfweightshelper.Init(100,60,edm::FileInPath("SUSYBSMAnalysis/RazorTuplizer/data/NNPDF30_nlo_as_0118_hessian_60.csv"));
	  firstPdfWeight = 9;
	  lastPdfWeight = 108;
	  firstAlphasWeight = 109;
	  lastAlphasWeight = 110; 
	}
	else if (pdfidx == 292200) {
	  //NNPdf30_nlo_as_0118 (nf5) with built-in alphas variations for NLO aMC@NLO samples
	  pdfweightshelper.Init(100,60,edm::FileInPath("SUSYBSMAnalysis/RazorTuplizer/data/NNPDF30_nlo_as_0118_hessian_60.csv"));
	  firstPdfWeight = 9;
	  lastPdfWeight = 108;
	  firstAlphasWeight = 109;
	  lastAlphasWeight = 110; 
	}   
	else if (pdfidx == 292000) {
	  //NNPdf30_nlo_as_0118_nf4 with built-in alphas variations for NLO aMC@NLO samples
	  pdfweightshelper.Init(100,60,edm::FileInPath("SUSYBSMAnalysis/RazorTuplizer/data/NNPDF30_nlo_as_0118_nf_4_hessian_60.csv"));
	  firstPdfWeight = 9;
	  lastPdfWeight = 108;
	  firstAlphasWeight = 109;
	  lastAlphasWeight = 110;
	}
	else {
	  firstPdfWeight = -1;
	  lastPdfWeight = -1;
	  firstAlphasWeight = -1;
	  lastAlphasWeight = -1;
	}    
      }*/
    }    
  }
  else {
    firstPdfWeight = -1;
    lastPdfWeight = -1;
    firstAlphasWeight = -1;
    lastAlphasWeight = -1;
  }      
  
}


//------ Method called for each lumi block ------//
void RazorTuplizer::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {

  if (useGen_) {  
    iLumi.getByToken(genLumiHeaderToken_,genLumiHeader);
  }
  
  //fill lhe comment lines with SUSY model parameter information
  lheComments = "";
  if (genLumiHeader.isValid() && isFastsim_) {
    lheComments = genLumiHeader->configDescription();    
  }    

  //Below: To check what the weights mean
  //for (uint i=0; i<genLumiHeader->weightNames().size(); ++i) cout << "Weights: " << genLumiHeader->weightNames()[i] << "\n";

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
    && fillPVAll()
    && fillMuons() 
    && fillElectrons()
    && fillTaus()
    && fillIsoPFCandidates()
    && fillPhotons(iEvent,iSetup)
    && fillJets()
    && fillJetsAK8()
    && fillMet(iEvent);
  //NOTE: if any of the above functions return false, the event will be rejected immediately with no further processing
  
  bool isGoodMCEvent = true;
  if (useGen_) {
    isGoodMCEvent = fillMC()
      && fillPileUp()
      && fillGenParticles();
  }
  
  isGoodEvent = isGoodEvent&&isGoodMCEvent;
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

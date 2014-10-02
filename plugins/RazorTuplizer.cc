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
  outputTree = fs->make<TTree>("outputTree", "selected miniAOD information");
  //setBranches();
}

RazorTuplizer::~RazorTuplizer()
{
}

//------ Enable the desired set of branches ------//
void RazorTuplizer::setBranches(){
  enableEventInfoBranches();
  enableMuonBranches();
  enableElectronBranches();
  enableTauBranches();
  enablePhotonBranches();
  enableJetBranches();
  enableJetAK8Branches();
  enableMetBranches();
  enableRazorBranches();
  enableMCBranches();
}

void RazorTuplizer::enableEventInfoBranches(){
  outputTree->Branch("nPV", &nPV, "nPV/I");
  outputTree->Branch("runNum", &runNum, "runNum/I");
  outputTree->Branch("lumiNum", &lumiNum, "lumiNum/I");
  outputTree->Branch("eventNum", &eventNum, "eventNum/I");
}

void RazorTuplizer::enableMuonBranches(){
  outputTree->Branch("nMuons", &nMuons,"nMuons/I");
  outputTree->Branch("muonE", muonE,"muonE[nMuons]/F");
  outputTree->Branch("muonPt", muonPt,"muonPt[nMuons]/F");
  outputTree->Branch("muonEta", muonEta,"muonEta[nMuons]/F");
  outputTree->Branch("muonPhi", muonPhi,"muonPhi[nMuons]/F");
}

void RazorTuplizer::enableElectronBranches(){
  outputTree->Branch("nElectrons", &nElectrons,"nElectrons/I");
  outputTree->Branch("eleE", eleE,"eleE[nElectrons]/F");
  outputTree->Branch("elePt", elePt,"elePt[nElectrons]/F");
  outputTree->Branch("eleEta", eleEta,"eleEta[nElectrons]/F");
  outputTree->Branch("elePhi", elePhi,"elePhi[nElectrons]/F");
}

void RazorTuplizer::enableTauBranches(){
  outputTree->Branch("nTaus", &nTaus,"nTaus/I");
  outputTree->Branch("tauE", tauE,"tauE[nTaus]/F");
  outputTree->Branch("tauPt", tauPt,"tauPt[nTaus]/F");
  outputTree->Branch("tauEta", tauEta,"tauEta[nTaus]/F");
  outputTree->Branch("tauPhi", tauPhi,"tauPhi[nTaus]/F");
}

void RazorTuplizer::enablePhotonBranches(){
  outputTree->Branch("nPhotons", &nPhotons,"nPhotons/I");
  outputTree->Branch("phoE", phoE,"phoE[nPhotons]/F");
  outputTree->Branch("phoPt", phoPt,"phoPt[nPhotons]/F");
  outputTree->Branch("phoEta", phoEta,"phoEta[nPhotons]/F");
  outputTree->Branch("phoPhi", phoPhi,"phoPhi[nPhotons]/F");
}

void RazorTuplizer::enableJetBranches(){
  outputTree->Branch("nJets", &nJets,"nJets/I");
  outputTree->Branch("jetE", jetE,"jetE[nJets]/F");
  outputTree->Branch("jetPt", jetPt,"jetPt[nJets]/F");
  outputTree->Branch("jetEta", jetEta,"jetEta[nJets]/F");
  outputTree->Branch("jetPhi", jetPhi,"jetPhi[nJets]/F");
  outputTree->Branch("jetCSV", jetCSV,"jetCSV[nJets]/F");
  outputTree->Branch("jetCISV", jetCISV,"jetCISV[nJets]/F");
}

void RazorTuplizer::enableJetAK8Branches(){
  outputTree->Branch("nFatJets", &nFatJets,"nFatJets/I");
  outputTree->Branch("fatJetE", fatJetE,"fatJetE[nFatJets]/F");
  outputTree->Branch("fatJetPt", fatJetPt,"fatJetPt[nFatJets]/F");
  outputTree->Branch("fatJetEta", fatJetEta,"fatJetEta[nFatJets]/F");
  outputTree->Branch("fatJetPhi", fatJetPhi,"fatJetPhi[nFatJets]/F");
}

void RazorTuplizer::enableMetBranches(){
  outputTree->Branch("metPt", &metPt, "metPt/D");
  outputTree->Branch("metPhi", &metPhi, "metPhi/D");
}

void RazorTuplizer::enableRazorBranches(){
  outputTree->Branch("MR", &MR, "MR/D");
  outputTree->Branch("RSQ", &RSQ, "RSQ/D");
}

void RazorTuplizer::enableMCBranches(){
  outputTree->Branch("genMetPt", &genMetPt, "genMetPt/D");
  outputTree->Branch("genMetPhi", &genMetPhi, "genMetPhi/D");
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
  nMuons = 0;
  nElectrons = 0;
  nTaus = 0;
  nPhotons = 0;
  nJets = 0;
  nFatJets = 0;
  
  for(int i = 0; i < 99; i++){
    muonE[i] = 0.0;
    muonPt[i] = 0.0;
    muonEta[i] = 0.0;
    muonPhi[i] = 0.0;
    
    eleE[i] = 0.0;
    elePt[i] = 0.0;
    eleEta[i] = 0.0;
    elePhi[i] = 0.0;
    
    tauE[i] = 0.0;
    tauPt[i] = 0.0;
    tauEta[i] = 0.0;
    tauPhi[i] = 0.0;
    
    phoE[i] = 0.0;
    phoPt[i] = 0.0;
    phoEta[i] = 0.0;
    phoPhi[i] = 0.0;
    
    jetE[i] = 0.0;
    jetPt[i] = 0.0;
    jetEta[i] = 0.0;
    jetPhi[i] = 0.0;
    jetCSV[i] = 0.0;
    jetCISV[i] = 0.0;
    
    fatJetE[i] = 0.0;
    fatJetPt[i] = 0.0;
    fatJetEta[i] = 0.0;
    fatJetPhi[i] = 0.0;
  }

    metPt = -999;
    metPhi = -999;

    genMetPt = -999;
    genMetPhi = -999;
    
    MR = -999;
    RSQ = -999;
    
    nPV = -1;
    eventNum = 0;
    lumiNum = 0;
    runNum = 0;
}

//------ Methods to fill tree variables ------//

bool RazorTuplizer::fillEventInfo(const edm::Event& iEvent){
  //store basic event info
  runNum = iEvent.id().run();
  lumiNum = iEvent.luminosityBlock();
  eventNum = iEvent.id().event();
  
  //select the primary vertex, if any
  if (vertices->empty()) return false; // skip the event if no PV found
  //const reco::Vertex &PV = vertices->front();
  nPV = vertices->size();
  
  return true;
}

bool RazorTuplizer::fillMuons(){
  for(const pat::Muon &mu : *muons){
    if(mu.pt() < 5 || !mu.isLooseMuon()) continue;
    muonE[nMuons] = mu.energy();
    muonPt[nMuons] = mu.pt();
    muonEta[nMuons] = mu.eta();
    muonPhi[nMuons] = mu.phi();
    nMuons++;
  }
  
  return true;
}

bool RazorTuplizer::fillElectrons(){
  for(const pat::Electron &ele : *electrons){
    if(ele.pt() < 5) continue;
    eleE[nElectrons] = ele.energy();
    elePt[nElectrons] = ele.pt();
    eleEta[nElectrons] = ele.eta();
    elePhi[nElectrons] = ele.phi();
    nElectrons++;
  }   
  
  return true;
}

bool RazorTuplizer::fillTaus(){
  for (const pat::Tau &tau : *taus) {
    if (tau.pt() < 20) continue;
    tauE[nTaus] = tau.energy();
    tauPt[nTaus] = tau.pt();
    tauEta[nTaus] = tau.eta();
    tauPhi[nTaus] = tau.phi();
    nTaus++;
  }
  
  return true;
}

bool RazorTuplizer::fillPhotons(){
  for (const pat::Photon &pho : *photons) {
    if (pho.pt() < 20 or pho.chargedHadronIso()/pho.pt() > 0.3) continue;
    phoE[nPhotons] = pho.energy();
    phoPt[nPhotons] = pho.pt();
    phoEta[nPhotons] = pho.eta();
    phoPhi[nPhotons] = pho.phi();
    nPhotons++;
  }
  
  return true;
}

bool RazorTuplizer::fillJets(){
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
}

bool RazorTuplizer::fillJetsAK8(){
  for (const pat::Jet &j : *jetsAK8) {
    fatJetE[nFatJets] = j.energy();
    fatJetPt[nFatJets] = j.pt();
    fatJetEta[nFatJets] = j.eta();
    fatJetPhi[nFatJets] = j.phi();
    nFatJets++;
  }
  
  return true;
}

bool RazorTuplizer::fillMet(){
  const pat::MET &Met = mets->front();
  metPt = Met.pt();
  metPhi = Met.phi();
  
  return true;
}

bool RazorTuplizer::fillMC(){
    const pat::MET &Met = mets->front();
    genMetPt = Met.genMET()->pt();
    genMetPhi = Met.genMET()->phi();

    return true;
}

bool RazorTuplizer::fillRazor(){ 
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
}

//------ Method called for each event ------//

void RazorTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;
  
  //initialize
  resetBranches();
  loadEvent(iEvent); //loads objects and resets tree branches
  
  //filler methods should fill relevant tree variables and return false if the event should be rejected
  bool isGoodEvent = 
    fillEventInfo(iEvent)
    && fillMuons() 
    && fillElectrons()
    && fillTaus()
    && fillPhotons()
    && fillJets()
    && fillJetsAK8()
    && fillMet()
    && fillRazor()
    && fillMC();
  
  //fill the tree if the event wasn't rejected
  if(isGoodEvent) outputTree->Fill();
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

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
  if (enableTriggerInfo_) enableTriggerBranches();
  enableMCBranches();
}

void RazorTuplizer::enableEventInfoBranches(){
  RazorEvents->Branch("nPV", &nPV, "nPV/I");
  RazorEvents->Branch("runNum", &runNum, "runNum/I");
  RazorEvents->Branch("lumiNum", &lumiNum, "lumiNum/I");
  RazorEvents->Branch("eventNum", &eventNum, "eventNum/I");
}

void RazorTuplizer::enableMuonBranches(){
  RazorEvents->Branch("nMuons", &nMuons,"nMuons/I");
  RazorEvents->Branch("muonE", muonE,"muonE[nMuons]/F");
  RazorEvents->Branch("muonPt", muonPt,"muonPt[nMuons]/F");
  RazorEvents->Branch("muonEta", muonEta,"muonEta[nMuons]/F");
  RazorEvents->Branch("muonPhi", muonPhi,"muonPhi[nMuons]/F");
}

void RazorTuplizer::enableElectronBranches(){
  RazorEvents->Branch("nElectrons", &nElectrons,"nElectrons/I");
  RazorEvents->Branch("eleE", eleE,"eleE[nElectrons]/F");
  RazorEvents->Branch("elePt", elePt,"elePt[nElectrons]/F");
  RazorEvents->Branch("eleEta", eleEta,"eleEta[nElectrons]/F");
  RazorEvents->Branch("elePhi", elePhi,"elePhi[nElectrons]/F");
}

void RazorTuplizer::enableTauBranches(){
  RazorEvents->Branch("nTaus", &nTaus,"nTaus/I");
  RazorEvents->Branch("tauE", tauE,"tauE[nTaus]/F");
  RazorEvents->Branch("tauPt", tauPt,"tauPt[nTaus]/F");
  RazorEvents->Branch("tauEta", tauEta,"tauEta[nTaus]/F");
  RazorEvents->Branch("tauPhi", tauPhi,"tauPhi[nTaus]/F");
}

void RazorTuplizer::enablePhotonBranches(){
  RazorEvents->Branch("nPhotons", &nPhotons,"nPhotons/I");
  RazorEvents->Branch("phoE", phoE,"phoE[nPhotons]/F");
  RazorEvents->Branch("phoPt", phoPt,"phoPt[nPhotons]/F");
  RazorEvents->Branch("phoEta", phoEta,"phoEta[nPhotons]/F");
  RazorEvents->Branch("phoPhi", phoPhi,"phoPhi[nPhotons]/F");
}

void RazorTuplizer::enableJetBranches(){
  RazorEvents->Branch("nJets", &nJets,"nJets/I");
  RazorEvents->Branch("jetE", jetE,"jetE[nJets]/F");
  RazorEvents->Branch("jetPt", jetPt,"jetPt[nJets]/F");
  RazorEvents->Branch("jetEta", jetEta,"jetEta[nJets]/F");
  RazorEvents->Branch("jetPhi", jetPhi,"jetPhi[nJets]/F");
  RazorEvents->Branch("jetCSV", jetCSV,"jetCSV[nJets]/F");
  RazorEvents->Branch("jetCISV", jetCISV,"jetCISV[nJets]/F");
}

void RazorTuplizer::enableJetAK8Branches(){
  RazorEvents->Branch("nFatJets", &nFatJets,"nFatJets/I");
  RazorEvents->Branch("fatJetE", fatJetE,"fatJetE[nFatJets]/F");
  RazorEvents->Branch("fatJetPt", fatJetPt,"fatJetPt[nFatJets]/F");
  RazorEvents->Branch("fatJetEta", fatJetEta,"fatJetEta[nFatJets]/F");
  RazorEvents->Branch("fatJetPhi", fatJetPhi,"fatJetPhi[nFatJets]/F");
}

void RazorTuplizer::enableMetBranches(){
  RazorEvents->Branch("metPt", &metPt, "metPt/F");
  RazorEvents->Branch("metPhi", &metPhi, "metPhi/F");
}

void RazorTuplizer::enableRazorBranches(){
  RazorEvents->Branch("MR", &MR, "MR/F");
  RazorEvents->Branch("RSQ", &RSQ, "RSQ/F");
}

void RazorTuplizer::enableTriggerBranches(){
  nameHLT = new std::vector<std::string>; nameHLT->clear();
  RazorEvents->Branch("nameHLT", "std::vector<std::string>",&nameHLT);
}

void RazorTuplizer::enableMCBranches(){
  RazorEvents->Branch("nGenJets", &nGenJets, "nGenJets/I");
  RazorEvents->Branch("genJetE", genJetE, "genJetE[nGenJets]/F");
  RazorEvents->Branch("genJetPt", genJetPt, "genJetPt[nGenJets]/F");
  RazorEvents->Branch("genJetEta", genJetEta, "genJetEta[nGenJets]/F");
  RazorEvents->Branch("genJetPhi", genJetPhi, "genJetPhi[nGenJets]/F");

  RazorEvents->Branch("genMetPt", &genMetPt, "genMetPt/F");
  RazorEvents->Branch("genMetPhi", &genMetPhi, "genMetPhi/F");
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
  nGenJets = 0;
  
  //nameHLT->clear();
  
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

    genJetE[i] = 0.0;
    genJetPt[i] = 0.0;
    genJetEta[i] = 0.0;
    genJetPhi[i] = 0.0;
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
  
  nPV = 0;
  //Check for good vertices
  for(unsigned int i = 0; i < vertices->size(); i++){
    if(vertices->at(i).isValid() && !vertices->at(i).isFake())nPV++;
  }
  if(nPV == 0)return false;
  
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

bool RazorTuplizer::fillTrigger(const edm::Event& iEvent){
  //fill trigger information

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  //  std::cout << "\n === TRIGGER PATHS === " << std::endl;

  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
    string hltPathNameReq = "HLT_";
    if (triggerBits->accept(i)) 
      if ((names.triggerName(i)).find(hltPathNameReq) != string::npos) 
	nameHLT->push_back(names.triggerName(i));
    //std::cout << "Trigger " << names.triggerName(i) << 
    //  ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
    //  ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)") 
    //	      << std::endl;
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
    && fillMuons() 
    && fillElectrons()
    && fillTaus()
    && fillPhotons()
    && fillJets()
    && fillJetsAK8()
    && fillMet()
    && fillRazor()
    && fillMC();

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

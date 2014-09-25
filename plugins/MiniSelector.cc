// -*- C++ -*-
// Class:      MiniSelector
/*
Description: Base class for miniAOD analysis with CRAB
*/
// Original Author:  Dustin James Anderson
//         Created:  Thu, 17 Jul 2014 15:00:06 GMT

#include "MiniSelector.h"
// constructors and destructor
MiniSelector::MiniSelector(const edm::ParameterSet& iConfig): 
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
}

MiniSelector::~MiniSelector()
{
}

// member functions

void MiniSelector::loadEvent(const edm::Event& iEvent){
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

//------ Method called for each event ------//

void MiniSelector::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    using namespace edm;

    //initialize
    loadEvent(iEvent);

}

//------ Method called once each job just before starting event loop ------//
void MiniSelector::beginJob(){
}

//------ Method called once each job just after ending the event loop ------//
void MiniSelector::endJob(){
}

//------ HELPER FUNCTIONS ------//

vector<TLorentzVector> MiniSelector::getHemispheres(vector<TLorentzVector> jets){
    int nJets = jets.size();
    vector<TLorentzVector> possibleHem1s; //holds possible hemisphere combinations
    vector<TLorentzVector> possibleHem2s;

    //stolen from https://github.com/pierinim/BSMatLHC/blob/master/BSMApp/src/CMS/CMSHemisphere.cc
    int nComb = pow(2, nJets);

    //step 1: store all possible partitions of the input jets
    int j_count;
    for(int i = 1; i < nComb-1; i++){ //note we omit the trivial hemisphere combinations (0 and nComb-1)
        TLorentzVector j_temp1, j_temp2;
        int itemp = i;
        j_count = nComb/2;
        int count = 0;
        while(j_count > 0){ //decompose i into binary '1's and '0's ; put the '1' jets into j_temp1 and the '0' jets into j_temp2
            if(itemp/j_count == 1){
                j_temp1 += jets[count];
            } else {
                j_temp2 += jets[count];
            }
            itemp -= j_count*(itemp/j_count); //note this is always (0 or 1)*j_count
            j_count /= 2;
            count++;
        }
        possibleHem1s.push_back(j_temp1);
        possibleHem2s.push_back(j_temp2);
    }
 
    //step 2: choose the partition that minimizes m1^2 + m2^2
    double mMin = -1;
    TLorentzVector myHem1;
    TLorentzVector myHem2;
    for(size_t i=0; i < possibleHem1s.size(); i++){
        double mTemp = possibleHem1s[i].M2() + possibleHem2s[i].M2();
        if(mMin < 0 || mTemp < mMin){
            mMin = mTemp;
            myHem1 = possibleHem1s[i];
            myHem2 = possibleHem2s[i];
        }
    }

    //return the hemispheres in decreasing order of pt
    vector<TLorentzVector> hemsOut;
    if(myHem1.Pt() > myHem2.Pt()){
        hemsOut.push_back(myHem1);
        hemsOut.push_back(myHem2);
    } else {
        hemsOut.push_back(myHem2);
        hemsOut.push_back(myHem1);
    }

    return hemsOut;
}

double MiniSelector::computeMR(TLorentzVector hem1, TLorentzVector hem2){
    return sqrt(pow(hem1.P() + hem2.P(), 2) - pow(hem1.Pz() + hem2.Pz(), 2));
}

double MiniSelector::computeR2(TLorentzVector hem1, TLorentzVector hem2, TLorentzVector pfMet){
    double mR = computeMR(hem1, hem2);
    double term1 = pfMet.Pt()/2*(hem1.Pt() + hem2.Pt());
    double term2 = pfMet.Px()/2*(hem1.Px() + hem2.Px()) + pfMet.Py()/2*(hem1.Py() + hem2.Py()); //dot product of MET with (p1T + p2T)
    double mTR = sqrt(term1 - term2);
    return (mTR / mR) * (mTR / mR);
}

bool MiniSelector::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle){
    //particle is already the ancestor
    if(ancestor == particle ) return true;

    //otherwise loop on mothers, if any and return true if the ancestor is found
    for(size_t i=0;i< particle->numberOfMothers();i++)
    {
        if(isAncestor(ancestor,particle->mother(i))) return true;
    }
    //if we did not return yet, then particle and ancestor are not relatives
    return false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniSelector);

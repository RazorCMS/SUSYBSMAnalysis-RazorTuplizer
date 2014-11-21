// -*- C++ -*-
// Class:      RazorTuplizer
/*
Description: Base class for miniAOD analysis with CRAB
*/
//         Author:  Caltech razor team
//         Created:  Thu, 17 Jul 2014 15:00:06 GMT

#ifndef RAZORTUPLIZER_H
#define RAZORTUPLIZER_H

// system include files
#include <memory>
#include <string>
#include <vector>

using namespace std;

// CMSSW framework includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//CMSSW package includes
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "EgammaAnalysis/ElectronTools/interface/EGammaMvaEleEstimatorCSA14.h"

//ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

//------ Class declaration ------//

class RazorTuplizer : public edm::EDAnalyzer {
public:
  //analyzer constructor and destructor
  explicit RazorTuplizer(const edm::ParameterSet&);
  ~RazorTuplizer();
  
  void loadEvent(const edm::Event& iEvent); //call at the beginning of each event to get input handles from the python config
  virtual void resetBranches();
  
  //enable desired output variables
  virtual void enableEventInfoBranches();
  virtual void setBranches();
  virtual void enableMuonBranches();
  virtual void enableElectronBranches();
  virtual void enableTauBranches();
  virtual void enablePhotonBranches();
  virtual void enableJetBranches();
  virtual void enableJetAK8Branches();
  virtual void enableMetBranches();
  virtual void enableRazorBranches();
  virtual void enableTriggerBranches();
  virtual void enableMCBranches();
  
  //select objects and fill tree branches
  virtual bool fillEventInfo(const edm::Event& iEvent);
  virtual bool fillMuons();//Fills looseID muon 4-momentum only. PT > 5GeV
  virtual bool fillElectrons();//Fills Ele 4-momentum only. PT > 5GeV
  virtual bool fillTaus();//Fills Tau 4-momentum only. PT > 20GeV
  virtual bool fillPhotons();//Fills photon 4-momentum only. PT > 20GeV && ISO < 0.3
  virtual bool fillJets();//Fills AK5 Jet 4-momentum, CSV, and CISV. PT > 20GeV 
  virtual bool fillJetsAK8();//Fills AK8 Jet 4-momentum.
  virtual bool fillMet();//Fills MET(mag, phi)
  virtual bool fillRazor();//Fills MR and RSQ
  virtual bool fillTrigger(const edm::Event& iEvent);//Fills trigger information
  virtual bool fillMC();
  
  //------ HELPER FUNCTIONS ------//
  
  //splits jets into two hemisperes for razor variable calculation
  //(minimizes sum of mass^2's of hemispheres)
  vector<TLorentzVector> getHemispheres(vector<TLorentzVector> jets);
  //compute M_R using two hemispheres
  float computeMR(TLorentzVector hem1, TLorentzVector hem2);
  //compute R^2 using two hemispheres and MET vector
  float computeR2(TLorentzVector hem1, TLorentzVector hem2, TLorentzVector pfMet);
  //returns true if particle 1 is an ancestor of particle 2, false otherwise
  //(takes two members of prunedGenParticles)
  bool isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle);

protected:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  //----- Member data ------//
  
  //EDM tokens for each miniAOD input object
  edm::EDGetTokenT<reco::VertexCollection> verticesToken_;
  edm::EDGetTokenT<pat::MuonCollection> muonsToken_;
  edm::EDGetTokenT<pat::ElectronCollection> electronsToken_;
  edm::EDGetTokenT<pat::TauCollection> tausToken_;
  edm::EDGetTokenT<pat::PhotonCollection> photonsToken_;
  edm::EDGetTokenT<pat::JetCollection> jetsToken_;
  edm::EDGetTokenT<pat::JetCollection> jetsAK8Token_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> packedPFCandsToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenParticlesToken_;
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenParticlesToken_;
  edm::EDGetTokenT<reco::GenJetCollection> genJetsToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBitsToken_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectsToken_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescalesToken_;
  edm::EDGetTokenT<pat::METCollection> metToken_;
  edm::EDGetTokenT<edm::TriggerResults> metFilterBitsToken_;
  edm::EDGetTokenT<LHEEventProduct> lheInfoToken_;
  edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puInfoToken_;
  edm::EDGetTokenT<HcalNoiseSummary> hcalNoiseInfoToken_;
  edm::EDGetTokenT<vector<reco::VertexCompositePtrCandidate> > secondaryVerticesToken_;
  edm::EDGetTokenT<double> rhoAllToken_;
  edm::EDGetTokenT<double> rhoFastjetAllToken_;
  edm::EDGetTokenT<double> rhoFastjetAllCaloToken_;
  edm::EDGetTokenT<double> rhoFastjetCentralCaloToken_;
  edm::EDGetTokenT<double> rhoFastjetCentralChargedPileUpToken_;
  edm::EDGetTokenT<double> rhoFastjetCentralNeutralToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > ebRecHitsToken_;
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > eeRecHitsToken_;
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > esRecHitsToken_;
  edm::EDGetTokenT<vector<reco::CaloCluster> > ebeeClustersToken_;
  edm::EDGetTokenT<vector<reco::CaloCluster> > esClustersToken_;
  edm::EDGetTokenT<vector<reco::Conversion> > conversionsToken_;
  edm::EDGetTokenT<vector<reco::Conversion> > singleLegConversionsToken_;
  edm::EDGetTokenT<vector<reco::GsfElectronCore> > gedGsfElectronCoresToken_;
  edm::EDGetTokenT<vector<reco::PhotonCore> > gedPhotonCoresToken_;
  edm::EDGetTokenT<vector<reco::SuperCluster> > superClustersToken_;
  edm::EDGetTokenT<vector<pat::PackedCandidate> > lostTracksToken_;
  
  //EDM handles for each miniAOD input object
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  edm::Handle<edm::TriggerResults> metFilterBits;
  edm::Handle<reco::VertexCollection> vertices;
  edm::Handle<pat::PackedCandidateCollection> packedPFCands;
  edm::Handle<pat::MuonCollection> muons;
  edm::Handle<pat::ElectronCollection> electrons;
  edm::Handle<pat::PhotonCollection> photons;
  edm::Handle<pat::TauCollection> taus;
  edm::Handle<pat::JetCollection> jets;
  edm::Handle<pat::JetCollection> jetsAK8;
  edm::Handle<pat::METCollection> mets;
  edm::Handle<edm::View<reco::GenParticle> > prunedGenParticles;
  edm::Handle<edm::View<pat::PackedGenParticle> > packedGenParticles;
  edm::Handle<reco::GenJetCollection> genJets;
  edm::Handle<LHEEventProduct> lheInfo;
  edm::Handle<GenEventInfoProduct> genInfo;
  edm::Handle<std::vector<PileupSummaryInfo> > puInfo;
  edm::Handle<HcalNoiseSummary> hcalNoiseInfo;
  edm::Handle<vector<reco::VertexCompositePtrCandidate> > secondaryVertices;
  edm::Handle<double> rhoAll;
  edm::Handle<double> rhoFastjetAll;
  edm::Handle<double> rhoFastjetAllCalo;
  edm::Handle<double> rhoFastjetCentralCalo;
  edm::Handle<double> rhoFastjetCentralChargedPileUp;
  edm::Handle<double> rhoFastjetCentralNeutral;
  edm::Handle<reco::BeamSpot> beamSpot;
  edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > ebRecHits;
  edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > eeRecHits;
  edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > esRecHits;
  edm::Handle<vector<reco::CaloCluster> > ebeeClusters;
  edm::Handle<vector<reco::CaloCluster> > esClusters;
  edm::Handle<vector<reco::Conversion> > conversions;
  edm::Handle<vector<reco::Conversion>> singleLegConversions;
  edm::Handle<vector<reco::GsfElectronCore> > gedGsfElectronCores;
  edm::Handle<vector<reco::PhotonCore> > gedPhotonCores;
  edm::Handle<vector<reco::SuperCluster> > superClusters;
  edm::Handle<vector<pat::PackedCandidate> > lostTracks;
  
  //output tree
  TTree *RazorEvents;
  TH1F *NEvents;

  //------ Variables for tree ------//

  //Muons
  int nMuons;
  float muonE[99];
  float muonPt[99];
  float muonEta[99];
  float muonPhi[99];
  
  //Electrons
  int nElectrons;
  float eleE[99];
  float elePt[99];
  float eleEta[99];
  float elePhi[99];
  
  //Taus
  int nTaus;
  float tauE[99];
  float tauPt[99];
  float tauEta[99];
  float tauPhi[99];
  
  //Photons
  int nPhotons;
  float phoE[99];
  float phoPt[99];
  float phoEta[99];
  float phoPhi[99];
  
  //AK4 Jets
  int nJets;
  float jetE[99];
  float jetPt[99];
  float jetEta[99];
  float jetPhi[99];
  float jetCSV[99];
  float jetCISV[99];
  
  //AK8 Jets
  int nFatJets;
  float fatJetE[99];
  float fatJetPt[99];
  float fatJetEta[99];
  float fatJetPhi[99];
  
  //MET 
  float metPt;
  float metPhi;

  //MC
  int nGenJets;
  float genJetE[99];
  float genJetPt[99];
  float genJetEta[99];
  float genJetPhi[99];

  float genMetPt;
  float genMetPhi;
  
  //razor variables
  float MR, RSQ;
  float MR_AK8, RSQ_AK8;
  
  //event info
  int nPV;
  int runNum;
  int lumiNum;
  int eventNum;

  //trigger info
  vector<string>  *nameHLT;

};

#endif

import FWCore.ParameterSet.Config as cms

#------ Setup ------#

#initialize the process
process = cms.Process("razorAna")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("Configuration.EventContent.EventContent_cff")

#load input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'root://xrootd-cms.infn.it//store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/003E832C-8AFC-E311-B7AA-002590596490.root' #PU 20, 25 ns CSA14 TTJets sample
        #'root://eoscms//eos/cms/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU40bx25_POSTLS170_V7-v2/00000/00800BE3-E826-E411-AD01-20CF3019DEE9.root' #PU 40, 25ns
        'root://xrootd.unl.edu://store/mc/Phys14DR/BlackHole_BH1_MD2000_MBH3000_n4_13TeV_TuneZ2star-blackmax/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/2E32617C-D56F-E411-9F23-003048D3CD92.root' #PHYS14 sample
    )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

#TFileService for output 
process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string("razorNtupleAna.root"),
                                   closeFileFast = cms.untracked.bool(True)
)

#load run conditions
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

#------ Declare the correct global tag ------#

#global tag for CSA14 25ns 20 PU scenario
process.GlobalTag.globaltag = 'PHYS14_25_V1::All'
#global tag for CSA14 50ns 40 PU scenario
#process.GlobalTag.globaltag = 'PHYS14_50_V1::All'

#------ If we add any inputs beyond standard miniAOD event content, import them here ------#

#------ Analyzer ------#

#list input collections
process.ntuples = cms.EDAnalyzer('RazorAna', 

    enableTriggerInfo = cms.bool(False),

    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    taus = cms.InputTag("slimmedTaus"),
    photons = cms.InputTag("slimmedPhotons"),
    jets = cms.InputTag("slimmedJets"),
    jetsAK8 = cms.InputTag("slimmedJetsAK8"),
    mets = cms.InputTag("slimmedMETs"),
    packedPfCands = cms.InputTag("packedPFCandidates"),

    packedGenParticles = cms.InputTag("packedGenParticles"),
    prunedGenParticles = cms.InputTag("prunedGenParticles"),
    genJets = cms.InputTag("slimmedGenJets", "", "PAT"),

    triggerBits = cms.InputTag("TriggerResults","","HLT"),
    triggerPrescales = cms.InputTag("patTrigger"),
    triggerObjects = cms.InputTag("selectedPatTrigger"),
    metFilterBits = cms.InputTag("TriggerResults", "", "PAT"),

    lheInfo = cms.InputTag("externalLHEProducer", "", "LHE"),
    genInfo = cms.InputTag("generator", "", "SIM"),
    #puInfo = cms.InputTag("mixData", "", "HLT"), #uncomment for samples with pre-mixed pileup
    puInfo = cms.InputTag("addPileupInfo", "", "HLT"), #uncomment if no pre-mixing
    hcalNoiseInfo = cms.InputTag("hcalnoise", "", "RECO"),

    secondaryVertices = cms.InputTag("slimmedSecondaryVertices", "", "PAT"),

    rhoAll = cms.InputTag("fixedGridRhoAll", "", "RECO"),
    rhoFastjetAll = cms.InputTag("fixedGridRhoFastjetAll", "", "RECO"),
    rhoFastjetAllCalo = cms.InputTag("fixedGridRhoFastjetAllCalo", "", "RECO"),
    rhoFastjetCentralCalo = cms.InputTag("fixedGridRhoFastjetCentralCalo", "", "RECO"),
    rhoFastjetCentralChargedPileUp = cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp", "", "RECO"),
    rhoFastjetCentralNeutral = cms.InputTag("fixedGridRhoFastjetCentralNeutral", "", "RECO"),

    beamSpot = cms.InputTag("offlineBeamSpot", "", "RECO"),

    ebRecHits = cms.InputTag("reducedEgamma", "reducedEBRecHits", "PAT"),
    eeRecHits = cms.InputTag("reducedEgamma", "reducedEERecHits", "PAT"),
    esRecHits = cms.InputTag("reducedEgamma", "reducedESRecHits", "PAT"),
    ebeeClusters = cms.InputTag("reducedEgamma", "reducedEBEEClusters", "PAT"),
    esClusters = cms.InputTag("reducedEgamma", "reducedESClusters", "PAT"),
    conversions = cms.InputTag("reducedEgamma", "reducedConversions", "PAT"),
    singleLegConversions = cms.InputTag("reducedEgamma", "reducedSingleLegConversions", "PAT"),
    gedGsfElectronCores = cms.InputTag("reducedEgamma", "reducedGedGsfElectronCores", "PAT"),
    gedPhotonCores = cms.InputTag("reducedEgamma", "reducedGedPhotonCores", "PAT"),
    superClusters = cms.InputTag("reducedEgamma", "reducedSuperClusters", "PAT"),

    lostTracks = cms.InputTag("lostTracks", "", "PAT")
)

#run
process.p = cms.Path(
        process.ntuples)

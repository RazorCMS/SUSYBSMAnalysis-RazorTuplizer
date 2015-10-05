import FWCore.ParameterSet.Config as cms

#------ Setup ------#

#initialize the process
process = cms.Process("razorTuplizer")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("Configuration.EventContent.EventContent_cff")

#load input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'/store/data/Run2015D/SingleMuon/MINIAOD/PromptReco-v3/000/256/677/00000/CEAE1A74-3A5F-E511-821F-02163E013938.root'
        '/store/data/Run2015D/HTMHT/MINIAOD/PromptReco-v3/000/256/936/00000/5CDEE11A-2F62-E511-B317-02163E014526.root'        
    )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

#TFileService for output 
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("razorNtuple.root"),
    closeFileFast = cms.untracked.bool(True)
)

#load run conditions
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")

#------ Declare the correct global tag ------#

#Global Tag for Run2015D
process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v2'

#------ If we add any inputs beyond standard miniAOD event content, import them here ------#

process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False) 
process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")

#------ Analyzer ------#

#list input collections
process.ntuples = cms.EDAnalyzer('RazorTuplizer', 
    isData = cms.bool(True),    
    useGen = cms.bool(False),
    isFastsim = cms.bool(False),
    enableTriggerInfo = cms.bool(True),                                 
    triggerPathNamesFile = cms.string("SUSYBSMAnalysis/RazorTuplizer/data/RazorHLTPathnames.dat"),
    eleHLTFilterNamesFile = cms.string("SUSYBSMAnalysis/RazorTuplizer/data/RazorElectronHLTFilterNames.dat"),
    muonHLTFilterNamesFile = cms.string("SUSYBSMAnalysis/RazorTuplizer/data/RazorMuonHLTFilterNames.dat"),
    photonHLTFilterNamesFile = cms.string("SUSYBSMAnalysis/RazorTuplizer/data/RazorPhotonHLTFilterNames.dat"),

    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    taus = cms.InputTag("slimmedTaus"),
    photons = cms.InputTag("slimmedPhotons"),
    jets = cms.InputTag("slimmedJets"),
    jetsPuppi = cms.InputTag("slimmedJetsPuppi"),
    jetsAK8 = cms.InputTag("slimmedJetsAK8"),
    mets = cms.InputTag("slimmedMETs"),
    metsNoHF = cms.InputTag("slimmedMETsNoHF"),
    metsPuppi = cms.InputTag("slimmedMETsPuppi"),
    packedPfCands = cms.InputTag("packedPFCandidates"),

    packedGenParticles = cms.InputTag("packedGenParticles"),
    prunedGenParticles = cms.InputTag("prunedGenParticles"),
    genJets = cms.InputTag("slimmedGenJets", "", "RECO"),

    triggerBits = cms.InputTag("TriggerResults","","HLT"),
    triggerPrescales = cms.InputTag("patTrigger"),
    triggerObjects = cms.InputTag("selectedPatTrigger"),
    metFilterBits = cms.InputTag("TriggerResults", "", "RECO"),
    hbheNoiseFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"),
    hbheTightNoiseFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun2Tight"),
    hbheIsoNoiseFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHEIsoNoiseFilterResult"),

    lheInfo = cms.InputTag("externalLHEProducer", "", "LHE"),
    genInfo = cms.InputTag("generator", "", "SIM"),
    puInfo = cms.InputTag("addPileupInfo", "", "HLT"), #uncomment if no pre-mixing
    #puInfo = cms.InputTag("mixData", "", "HLT"), #uncomment for samples with pre-mixed pileup
    hcalNoiseInfo = cms.InputTag("hcalnoise", "", "RECO"),

    secondaryVertices = cms.InputTag("slimmedSecondaryVertices", "", "RECO"),

    rhoAll = cms.InputTag("fixedGridRhoAll", "", "RECO"),
    rhoFastjetAll = cms.InputTag("fixedGridRhoFastjetAll", "", "RECO"),
    rhoFastjetAllCalo = cms.InputTag("fixedGridRhoFastjetAllCalo", "", "RECO"),
    rhoFastjetCentralCalo = cms.InputTag("fixedGridRhoFastjetCentralCalo", "", "RECO"),
    rhoFastjetCentralChargedPileUp = cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp", "", "RECO"),
    rhoFastjetCentralNeutral = cms.InputTag("fixedGridRhoFastjetCentralNeutral", "", "RECO"),

    beamSpot = cms.InputTag("offlineBeamSpot", "", "RECO"),

    ebRecHits = cms.InputTag("reducedEgamma", "reducedEBRecHits", "RECO"),
    eeRecHits = cms.InputTag("reducedEgamma", "reducedEERecHits", "RECO"),
    esRecHits = cms.InputTag("reducedEgamma", "reducedESRecHits", "RECO"),
    ebeeClusters = cms.InputTag("reducedEgamma", "reducedEBEEClusters", "RECO"),
    esClusters = cms.InputTag("reducedEgamma", "reducedESClusters", "RECO"),
    conversions = cms.InputTag("reducedEgamma", "reducedConversions", "RECO"),
    singleLegConversions = cms.InputTag("reducedEgamma", "reducedSingleLegConversions", "RECO"),
    gedGsfElectronCores = cms.InputTag("reducedEgamma", "reducedGedGsfElectronCores", "RECO"),
    gedPhotonCores = cms.InputTag("reducedEgamma", "reducedGedPhotonCores", "RECO"),
    superClusters = cms.InputTag("reducedEgamma", "reducedSuperClusters", "RECO"),

    lostTracks = cms.InputTag("lostTracks", "", "RECO")
)

#run
process.p = cms.Path( process.HBHENoiseFilterResultProducer*
                      process.ntuples)

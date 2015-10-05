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
        '/store/mc/RunIISpring15FSPremix/SMS-T1bbbb_mGluino-1150_mLSP-400to975-1100to1125_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/MCRUN2_74_V9-v1/50000/320B2CE3-BE5D-E511-8C9E-B083FED76C6C.root'
    )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

#TFileService for output 
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("razorNtuple.root"),
    closeFileFast = cms.untracked.bool(True)
)

#load run conditions
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

#------ Declare the correct global tag ------#

#global tag for CSA14 25ns 20 PU (asymptotic alignment and calibration) scenario
#process.GlobalTag.globaltag = 'PLS170_V7AN1::All'
#global tag for CSA14 50ns 40 PU (more pessimistic alignment and calibration) scenario
#process.GlobalTag.globaltag = 'PLS170_V6AN1::All'
#global tag for PHYS14 asymptotic 25ns scenario
process.GlobalTag.globaltag = 'MCRUN2_74_V9::All'

#------ If we add any inputs beyond standard miniAOD event content, import them here ------#

#process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
#process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)

#------ Analyzer ------#

#list input collections
process.ntuples = cms.EDAnalyzer('RazorTuplizer', 
    isData = cms.bool(False),    
    useGen = cms.bool(True),
    isFastsim = cms.bool(True),
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
    genJets = cms.InputTag("slimmedGenJets", "", "PAT"),

    triggerBits = cms.InputTag("TriggerResults","","HLT"),
    triggerPrescales = cms.InputTag("patTrigger"),
    triggerObjects = cms.InputTag("selectedPatTrigger"),
    metFilterBits = cms.InputTag("TriggerResults", "", "PAT"),
    hbheNoiseFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"),
    hbheTightNoiseFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun2Tight"),
    hbheIsoNoiseFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHEIsoNoiseFilterResult"),

    lheInfo = cms.InputTag("source", "", "LHEFile"),
    genInfo = cms.InputTag("generator", "", "HLT"),
    puInfo = cms.InputTag("addPileupInfo", "", "HLT"), #uncomment if no pre-mixing
    #puInfo = cms.InputTag("mixData", "", "HLT"), #uncomment for samples with pre-mixed pileup
    hcalNoiseInfo = cms.InputTag("hcalnoise", "", "HLT"),

    secondaryVertices = cms.InputTag("slimmedSecondaryVertices", "", "PAT"),

    rhoAll = cms.InputTag("fixedGridRhoAll", "", "HLT"),
    rhoFastjetAll = cms.InputTag("fixedGridRhoFastjetAll", "", "HLT"),
    rhoFastjetAllCalo = cms.InputTag("fixedGridRhoFastjetAllCalo", "", "HLT"),
    rhoFastjetCentralCalo = cms.InputTag("fixedGridRhoFastjetCentralCalo", "", "HLT"),
    rhoFastjetCentralChargedPileUp = cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp", "", "HLT"),
    rhoFastjetCentralNeutral = cms.InputTag("fixedGridRhoFastjetCentralNeutral", "", "HLT"),

    beamSpot = cms.InputTag("offlineBeamSpot", "", "HLT"),

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
process.p = cms.Path( #process.HBHENoiseFilterResultProducer*
                      process.ntuples)

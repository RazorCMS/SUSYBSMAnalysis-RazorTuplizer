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
	#PU0, timing
	'/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/069E588E-4AA1-E611-AEF1-0025905B856E.root'
	#'/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/MINIAODSIM/81X_upgrade2023_realistic_v3_2023D3Timing13TeV-v1/10000/C8A548B5-7199-E611-8D6D-0CC47A4D7658.root',
	#'/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/MINIAODSIM/81X_upgrade2023_realistic_v3_2023D3Timing13TeV-v1/10000/D4C10DB4-7199-E611-B0CB-0CC47A78A45A.root'	
	#PU140, timing
	#'/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/MINIAODSIM/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/687D8419-FCA1-E611-8ACA-0025905A6118.root',
	#'/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/MINIAODSIM/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/B0190804-1DA2-E611-80C9-0CC47A4D7616.root',
	#'/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/MINIAODSIM/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/E2E812E5-8EA2-E611-B504-0CC47A4D76D6.root'
	#PU200, timing
 	#'/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/MINIAODSIM/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/EA9A2A23-C1A6-E611-82F1-0CC47A7C34D0.root',
	#'/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/MINIAODSIM/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/EEF3D529-C1A6-E611-91C8-0CC47A4D7698.root'	
   )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
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


#process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_v3'
process.GlobalTag.globaltag = '81X_upgrade2023_realistic_v3'

#------ If we add any inputs beyond standard miniAOD event content, import them here ------#

#process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
#process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
#process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False) 
#process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")

#process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
#process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
#process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
#process.BadChargedCandidateFilter.taggingMode = cms.bool(True)

#process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
#process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
#process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
#process.BadPFMuonFilter.taggingMode = cms.bool(True)

#------ Analyzer ------#

#list input collections
process.ntuples = cms.EDAnalyzer('RazorTuplizer', 
    isData = cms.bool(False),    
    useGen = cms.bool(True),
    isFastsim = cms.bool(False),
    enableTriggerInfo = cms.bool(True),                                 
    triggerPathNamesFile = cms.string("SUSYBSMAnalysis/RazorTuplizer/data/RazorHLTPathnames.dat"),
    eleHLTFilterNamesFile = cms.string("SUSYBSMAnalysis/RazorTuplizer/data/RazorElectronHLTFilterNames.dat"),
    muonHLTFilterNamesFile = cms.string("SUSYBSMAnalysis/RazorTuplizer/data/RazorMuonHLTFilterNames.dat"),
    photonHLTFilterNamesFile = cms.string("SUSYBSMAnalysis/RazorTuplizer/data/RazorPhotonHLTFilterNames.dat"),

    vertices = cms.InputTag("offlinePrimaryVertices4D"), # for timing case
    #vertices = cms.InputTag("offlinePrimaryVerticesWithBS"),  # for non-timing case
    
    muons = cms.InputTag("muons"),
    electrons = cms.InputTag("gedGsfElectrons"),
    taus = cms.InputTag("hpsPFTauProducer"),
    photons = cms.InputTag("gedPhotons"),
    jets = cms.InputTag("ak4PFJetsCHS"),
    jetsPuppi = cms.InputTag("ak4PFJets"),
    jetsAK8 = cms.InputTag("ak8PFJetsCHS"),
    mets = cms.InputTag("pfMet"),
    #metsNoHF = cms.InputTag("pfMet30"),
    metsPuppi = cms.InputTag("pfMet"),
    pfCands = cms.InputTag("particleFlow","","RECO"),

    #packedPfCands = cms.InputTag("packedPFCandidates"),

    genParticles = cms.InputTag("genParticles"),

    #packedGenParticles = cms.InputTag("packedGenParticles"),
    #prunedGenParticles = cms.InputTag("prunedGenParticles"),
    genJets = cms.InputTag("ak4GenJets"),

    triggerBits = cms.InputTag("TriggerResults","","HLT"),
    hepMC = cms.InputTag("generatorSmeared", "", "SIM"),
    
    #triggerPrescales = cms.InputTag("patTrigger"),
    
    #triggerObjects = cms.InputTag("selectedPatTrigger"),

    metFilterBits = cms.InputTag("TriggerResults", "", "RECO"),

    #hbheNoiseFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"),
    #hbheTightNoiseFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun2Tight"),
    #hbheIsoNoiseFilter = cms.InputTag("HBHENoiseFilterResultProducer","HBHEIsoNoiseFilterResult"),
    
    #BadChargedCandidateFilter = cms.InputTag("BadChargedCandidateFilter",""),
    #BadMuonFilter = cms.InputTag("BadPFMuonFilter",""),

    #lheInfo = cms.InputTag("externalLHEProducer", "", ""),
    genInfo = cms.InputTag("generator", "", "SIM"),
   
    tracks = cms.InputTag("generalTracks", "", "RECO"), 
    trackTime = cms.InputTag("trackTimeValueMapProducer","generalTracksConfigurableFlatResolutionModel"),
    trackTimeReso = cms.InputTag("trackTimeValueMapProducer","generalTracksConfigurableFlatResolutionModelResolution"),

    puInfo = cms.InputTag("slimmedAddPileupInfo", "", "RECO"), #uncomment if no pre-mixing
    #puInfo = cms.InputTag("mixData", "", "HLT"), #uncomment for samples with pre-mixed pileup
    #hcalNoiseInfo = cms.InputTag("hcalnoise", "", "RECO"),

    secondaryVertices = cms.InputTag("inclusiveSecondaryVertices", "", "RECO"),

    rhoAll = cms.InputTag("fixedGridRhoAll", "", "RECO"),
    
    rhoFastjetAll = cms.InputTag("fixedGridRhoFastjetAll", "", "RECO"),
    rhoFastjetAllCalo = cms.InputTag("fixedGridRhoFastjetAllCalo", "", "RECO"),
    rhoFastjetCentralCalo = cms.InputTag("fixedGridRhoFastjetCentralCalo", "", "RECO"),
    rhoFastjetCentralChargedPileUp = cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp", "", "RECO"),
    rhoFastjetCentralNeutral = cms.InputTag("fixedGridRhoFastjetCentralNeutral", "", "RECO"),

    beamSpot = cms.InputTag("offlineBeamSpot", "", "RECO"),
    pfClusters = cms.InputTag("particleFlowClusterECAL","","RECO"),
    #ebRecHits = cms.InputTag("reducedEgamma", "reducedEBRecHits", "RECO"),
    ebRecHits = cms.InputTag("ecalRecHit", "EcalRecHitsEB", "RECO"),
    eeRecHits = cms.InputTag("ecalRecHit", "EcalRecHitsEE", "RECO"),
    esRecHits = cms.InputTag("ecalRecHit", "EcalRecHitsES", "RECO"),
    #ebeeClusters = cms.InputTag("reducedEgamma", "reducedEBEEClusters", "RECO"),
    ebeeClusters = cms.InputTag("particleFlowEGamma", "EBEEClusters", "RECO"),
    esClusters = cms.InputTag("particleFlowEGamma", "ESClusters", "RECO"),
    #conversions = cms.InputTag("reducedEgamma", "reducedConversions", "RECO"),
    conversions = cms.InputTag("allConversions", "", "RECO"),
    
    #singleLegConversions = cms.InputTag("reducedEgamma", "reducedSingleLegConversions", "RECO"),
    singleLegConversions = cms.InputTag("particleFlowEGamma", "", "RECO"),
    
    gedGsfElectronCores = cms.InputTag("gedGsfElectronCores", "", "RECO"),
    gedPhotonCores = cms.InputTag("gedPhotonCore", "", "RECO"),
    #superClusters = cms.InputTag("reducedEgamma", "reducedSuperClusters", "RECO"),

    #lostTracks = cms.InputTag("lostTracks", "", "RECO")
)

#run
process.p = cms.Path( #process.HBHENoiseFilterResultProducer*
                      #process.BadChargedCandidateFilter*
                      #process.BadPFMuonFilter*
                      process.ntuples)

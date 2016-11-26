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
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/069E588E-4AA1-E611-AEF1-0025905B856E.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/1001BF85-48A1-E611-82B2-0CC47A4D76AC.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/386984EA-42A1-E611-B533-0025905B8598.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/3A248F3C-47A1-E611-89B6-0CC47A78A42C.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/3AEBA995-48A1-E611-BB55-0CC47A4C8E82.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/3EFE1EDE-41A1-E611-B9F7-0025905B8598.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/44F7E9CD-23A1-E611-B5BF-0025905A608A.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/4A540BF1-44A1-E611-B25F-0CC47A78A33E.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/5036E3CD-43A1-E611-B8E6-0CC47A7C351E.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/50B6402D-3EA1-E611-AC09-0025905A48FC.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/5663D6A6-21A1-E611-8464-0CC47A7C3444.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/5848AFC6-3CA1-E611-B7E1-0CC47A7C345C.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/62020C77-43A1-E611-961E-0CC47A4C8F30.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/685FA7A7-13A2-E611-BF6D-0025905A6104.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/74737119-4AA1-E611-93BB-0025905A6110.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/7A9BC26A-3DA1-E611-B910-0CC47A7452DA.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/7EA6BFC7-42A1-E611-83A5-0CC47A78A3B4.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/84D88B7D-3DA1-E611-9F71-0CC47A7C345C.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/8C92DCC4-42A1-E611-89FB-0CC47A7C34D0.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/9023B38D-31A1-E611-8DA5-0CC47A4D76CC.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/929F0BF2-50A1-E611-880A-0CC47A7452DA.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/96C21A8A-37A1-E611-8DF5-0CC47A4D76D6.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/A6AE5CC3-47A1-E611-B2B9-0CC47A78A360.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/AADB9AAD-3EA1-E611-AC46-0CC47A4C8EE2.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/ACE49A3B-41A1-E611-B165-0CC47A4C8EBA.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/B6051312-3BA1-E611-B1A5-0CC47A7C3628.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/B6D5F61D-3FA1-E611-B39C-0CC47A4C8F26.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/CA98E6CE-35A1-E611-BAA9-0CC47A4D768C.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/CEFAE993-48A1-E611-AB9D-0025905A6094.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/D014FF52-44A1-E611-9B94-0CC47A4C8E38.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/D46483C2-47A1-E611-8610-0CC47A78A45A.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/D64BFE1F-3BA1-E611-A119-0CC47A4D7658.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/D6D5D749-4EA1-E611-BD4F-0CC47A4D760A.root',
	'/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/DC45E1D2-3BA1-E611-AC00-0025905A60E4.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/EC8245CE-43A1-E611-ABBF-0CC47A78A340.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/ECB79E95-3EA1-E611-8B01-0CC47A4D76A0.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/F2C17FE6-29A1-E611-9977-0CC47A78A2F6.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/FAD00866-41A1-E611-AD2B-0025905B858C.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU140r1-v1/10000/FE34A90F-49A1-E611-842A-0CC47A4D76B2.root'
  )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

#TFileService for output 
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("razorNtuple_PU140_NoTiming.root"),
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

    #vertices = cms.InputTag("offlinePrimaryVertices4D"), # for timing case
    vertices = cms.InputTag("offlinePrimaryVerticesWithBS"),  # for non-timing case
    
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

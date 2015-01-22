from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola__Spring14miniaod-PU40bx25_POSTLS170_V7-v2__11Nov2014__V3' #change the request name for each new task
config.General.workArea = 'crab'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'python/razorTuplizer.py'
config.JobType.outputFiles = ['razorNtuple.root']

config.section_("Data")
#25ns CSA14 scenario (20 PU, processed with CMSSW_7_0_6_patch1)
#config.Data.inputDataset = '/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM'
#config.Data.inputDataset = '/SMS-T1tttt_2J_mGl-1500_mLSP-100_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM'

#config.Data.inputDataset = '/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM'

#config.Data.inputDataset = '/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU20bx25_POSTLS170_V5-v2/MINIAODSIM'
#config.Data.inputDataset = '/WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM'
#config.Data.inputDataset = '/WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM'
#config.Data.inputDataset = '/WJetsToLNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM'
#
#config.Data.inputDataset = '/DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-50_HT-200to400_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-50_HT-400to600_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-50_HT-600toInf_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM'
#
#config.Data.inputDataset = '/ZJetsToNuNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM'
#config.Data.inputDataset = '/ZJetsToNuNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM'
#config.Data.inputDataset = '/ZJetsToNuNu_HT-400to600_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU20bx25_POSTLS170_V5-v1/MINIAODSIM'
#config.Data.inputDataset = '/ZJetsToNuNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU20bx25_POSTLS170_V5-v2/MINIAODSIM'

#50ns CSA14 scenario (40 PU, processed with CMSSW_7_0_6_patch1)
#config.Data.inputDataset = '/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/StoreResults-Spring14dr_PU_S14_POSTLS170_V6AN1_miniAOD706p1_814812ec83fce2f620905d2bb30e9100-v2/USER'
#config.Data.inputDataset = '/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/xuchen-Spring14dr-PU_S14_POSTLS170_V6AN1-miniAOD706p1_t2-814812ec83fce2f620905d2bb30e9100/USER'
#config.Data.inputDataset = '/SMS-T1tttt_2J_mGl-1500_mLSP-100_Tune4C_13TeV-madgraph-tauola/duanders-Spring14dr-PU_S14_POSTLS170_V6AN1-miniAOD706p1-814812ec83fce2f620905d2bb30e9100/USER'
#config.Data.inputDataset = '/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/duanders-Spring14dr-PU_S14_POSTLS170_V6AN1-miniAOD706p1-v2-814812ec83fce2f620905d2bb30e9100/USER'

#CSA14 scenario 3 (40 PU, 25 ns)
config.Data.inputDataset = '/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU40bx25_POSTLS170_V7-v2/MINIAODSIM'

config.Data.dbsUrl = 'global' #change this according to the DBS instance (usually 'global') of the target dataset
#config.Data.dbsUrl = 'phys03' 
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 3
config.Data.publication = False
#config.Data.publishDbsUrl = 'phys03' #enable for publishing
#config.Data.publishDataName = 'razorNtuple'
#config.Data.ignoreLocality = False #disable AAA
config.Data.ignoreLocality = True #enable AAA

config.section_("Site")
config.Site.storageSite = 'T2_US_Caltech'

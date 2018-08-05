import FWCore.ParameterSet.Config as cms

#------ Setup ------#

#initialize the process
process = cms.Process("razorTuplizer")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("Configuration.EventContent.EventContent_cff")

#load input files
readFiles = cms.untracked.vstring()
readFiles.extend( [
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/02527FB1-5CAD-E711-979A-24BE05C63681.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/04EB2BA6-66AD-E711-9003-24BE05C4D8C1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/0A13405A-5CAD-E711-A5F0-E0071B7A48A0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/0C06CE0E-5CAD-E711-A171-E0071B73B6E0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/0C69EDEE-62AD-E711-97DA-5065F3816251.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/0E1C6227-60AD-E711-A218-24BE05C626C1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/0E4450B3-5EAD-E711-B3DA-E0071B73B6D0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/0EE0318F-5FAD-E711-BC9C-24BE05C63721.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/10040861-5FAD-E711-A141-E0071B73C640.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/1014411E-63AD-E711-B78F-24BE05C4D821.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/104E907A-64AD-E711-8030-5065F381C1D1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/10666440-5DAD-E711-A561-24BE05C666B1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/1425AC18-61AD-E711-8CEE-E0071B7A3540.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/1A3AEEA1-66AD-E711-8761-24BE05C3CCE1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/1C136F44-5EAD-E711-A008-E0071B6C9DA0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/1C3E7769-5BAD-E711-9729-E0071B73B6E0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/1C5A8B25-65AD-E711-8EB2-E0071B7A3540.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/1CDF6498-61AD-E711-AE0D-E0071B7AF7C0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/208646A9-5FAD-E711-B3BA-24BE05CE2D41.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/20C66166-5BAD-E711-8D13-E0071B7AC760.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/20E5A37D-5BAD-E711-95FD-E0071B73B6F0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/20F8BD87-60AD-E711-AC3F-E0071B6CAD10.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/22A3B70E-60AD-E711-944B-E0071B7B2320.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/2654F07B-64AD-E711-81EF-5065F3818291.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/26765C8D-60AD-E711-8B36-E0071B7A48A0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/2678D709-67AD-E711-B79E-24BE05C3CBE1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/267A28F5-5FAD-E711-BF86-E0071B7AC710.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/26DE3AAA-66AD-E711-AA39-24BE05C38C91.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/285C74F7-5AAD-E711-BCD6-5065F3818261.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/2A25B55D-5CAD-E711-924A-24BE05C626B1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/2A4FC186-61AD-E711-BDF6-24BE05C616E1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/2AA2A529-5FAD-E711-903A-5065F3812261.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/2C51B100-67AD-E711-AB07-E0071B6C9DD0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/2EFBF57F-60AD-E711-B7B4-E0071B7A08F0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/30EA7399-5FAD-E711-B93C-24BE05C4D8F1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/3216AB20-5FAD-E711-8AA4-E0071B7AF7C0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/32AA99A9-67AD-E711-92E1-E0071B73B6E0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/342B9069-5BAD-E711-814B-E0071B73B6E0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/3488CF5E-5BAD-E711-9318-E0071B73C620.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/3632D58F-5FAD-E711-A079-24BE05C3B8B1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/3A57A924-5FAD-E711-A9FC-E0071B6C9DD0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/3E22D72C-5FAD-E711-B3BC-E0071B7AC760.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/3E275D7C-61AD-E711-AA4A-24BE05C3B8B1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/3E48CD1F-5FAD-E711-9CC2-E0071B7A9810.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/4024E6B2-5CAD-E711-9E8B-24BE05CEDC81.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/441AD8B2-5CAD-E711-8092-24BE05CEDC81.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/44B7B91E-63AD-E711-83C6-5065F37D4131.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/466D6DB0-65AD-E711-AD02-24BE05CEDC71.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/480A6EBA-5CAD-E711-B7D0-5065F38102F1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/48D69BAC-66AD-E711-8B40-24BE05C63651.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/4C5A1436-5FAD-E711-9697-E0071B741D70.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/4E76CD89-61AD-E711-870E-E0071B740D80.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/504EF3F8-5BAD-E711-9577-E0071B73C610.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/52AAA9DD-62AD-E711-857E-24BE05C63721.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/540EF4A2-60AD-E711-8CFA-E0071B6CAD20.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/562CEC4E-65AD-E711-830C-24BE05C666B1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/58C94C1E-5FAD-E711-8AC9-24BE05C3CCE1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/5A78B31F-5FAD-E711-8286-24BE05C6C7E1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/5A934542-5DAD-E711-9B33-E0071B74AC10.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/5C27F1F7-5BAD-E711-9048-24BE05C6C7E1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/5EC5E28A-5FAD-E711-929F-24BE05C4D8F1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/5EFA16E9-60AD-E711-8B7F-5065F37D4131.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/60960EEA-60AD-E711-A3AD-E0071B749C80.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/60A020CB-63AD-E711-B652-E0071B73C620.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/60F6A734-5EAD-E711-8F63-E0071B7AC750.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/6683A505-5CAD-E711-979C-E0071B73C610.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/6883E82C-5EAD-E711-B1AF-E0071B7AC700.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/68AB9CAC-66AD-E711-9571-24BE05C63651.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/68B7FE26-60AD-E711-8682-24BE05C626C1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/6C5CEB24-5EAD-E711-861E-E0071B7A9810.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/6CAB848A-60AD-E711-9763-24BE05CECB71.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/703BBA25-60AD-E711-A33F-24BE05CEBD61.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/70642FD4-63AD-E711-8D0C-E0071B7AF7C0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/7269769C-66AD-E711-9A8A-24BE05C666B1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/72992B41-5FAD-E711-8CDA-E0071B7A4550.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/72E01C4C-5DAD-E711-AD7F-E0071B740D80.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/74C4E530-65AD-E711-BD8C-24BE05BDBE21.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/76409D55-5BAD-E711-8A87-5065F37D4131.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/78025AD2-6CAD-E711-9BDF-E0071B7A4550.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/8006141F-5FAD-E711-99E1-24BE05C666B1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/8009D28E-65AD-E711-9129-E0071B7AC770.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/80A79CAC-66AD-E711-AD23-24BE05C63651.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/82776A54-5BAD-E711-A881-24BE05C626B1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/8413200A-5CAD-E711-952D-E0071B73B6E0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/8609E229-5FAD-E711-A2D1-E0071B7A5650.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/861DBBC1-5CAD-E711-AAAA-5065F381A2E1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/862CCC0E-63AD-E711-9B0E-24BE05BDBE21.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/867548D7-5CAD-E711-A916-E0071B6C9DC0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/86768DD9-5AAD-E711-B060-24BE05C4C891.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/8686E025-60AD-E711-B5BE-24BE05CEBD61.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/8696AEBB-5FAD-E711-854E-E0071B7AC770.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/88C89F7C-60AD-E711-B203-24BE05C4D851.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/8A7ACFAB-61AD-E711-8CA7-E0071B7A4550.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/8AE75033-5FAD-E711-AC6D-24BE05C45BF1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/8C05A972-66AD-E711-A03C-24BE05C4D821.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/8C5FE903-5EAD-E711-A42C-24BE05CEDC71.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/8CCCA12D-5FAD-E711-89D8-24BE05C49891.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/8EEA8412-63AD-E711-AFB8-24BE05CEEDE1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/90031697-5FAD-E711-B9EB-E0071B74AC00.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/900B1A02-66AD-E711-91AF-24BE05CEFDF1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/92DA2BA6-66AD-E711-89D6-24BE05C4D8C1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/941A2BB2-5AAD-E711-B3FB-5065F381E251.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/9456DA25-60AD-E711-9C1B-24BE05CEBD61.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/94A23B33-5FAD-E711-B65F-24BE05C3DB21.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/9835A0C5-6CAD-E711-964C-E0071B73C620.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/987F001F-5FAD-E711-8FA9-E0071B7A45D0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/98B4A734-5EAD-E711-AFE5-E0071B7AC750.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/9A8EC6AE-5CAD-E711-83C8-24BE05C46B01.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/9C4BCF1F-60AD-E711-AA58-24BE05C626C1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/9CB54DE6-60AD-E711-97B6-E0071B74AC00.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/9E80FFB6-5CAD-E711-8F4C-E0071B73B6D0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/9EF3F60E-5EAD-E711-B9AC-E0071B7A58B0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/A22CD4C0-5AAD-E711-B224-E0071B7A45D0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/A27EF11D-66AD-E711-AF8F-24BE05CEED91.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/A47F5038-62AD-E711-9AB7-24BE05C44BB1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/A6753E80-61AD-E711-90E6-E0071B74AC00.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/A6DBBE81-60AD-E711-BE3C-24BE05CEADD1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/A6FAE225-60AD-E711-9E22-24BE05CEBD61.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/A88AA64B-5DAD-E711-A788-24BE05CEFDF1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/B00AA83D-64AD-E711-BDC4-24BE05C68681.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/B017D11F-60AD-E711-B0AF-24BE05C626C1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/B030AFC2-66AD-E711-BF40-24BE05C38C91.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/B06029A9-5DAD-E711-9B76-E0071B7A9800.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/B0A387AE-5CAD-E711-8AC6-5065F38102F1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/B0A3F01D-66AD-E711-8DB4-24BE05CEED91.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/B206E540-67AD-E711-A523-24BE05C3DB21.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/B286C2A2-61AD-E711-A280-E0071B749CA0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/B2DF043F-5CAD-E711-AB54-E0071B741D70.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/B474E65D-5BAD-E711-85A3-E0071B7A3540.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/B4CBCD2F-5FAD-E711-9DD3-24BE05C6E591.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/B64CA12E-5FAD-E711-B833-E0071B73C600.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/B678CE42-5BAD-E711-A523-24BE05C63721.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/B89EF81E-5FAD-E711-831F-5065F3816291.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/BA6D4764-5EAD-E711-A687-E0071B6C9DC0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/BAAF2C41-5FAD-E711-B84A-E0071B7A4550.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/BC30299C-61AD-E711-BAB8-E0071B7B2320.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/C2A1BF9E-5BAD-E711-A059-E0071B7AC710.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/C44A7C3F-5DAD-E711-BFCD-5065F38122D1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/C4BBC947-5FAD-E711-A80B-E0071B7B2320.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/C6217739-5DAD-E711-A015-24BE05BDBE21.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/C6995527-5FAD-E711-965F-E0071B749C80.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/C8273F62-64AD-E711-9FDC-E0071B7A8570.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/C870DB5C-61AD-E711-9D7B-E0071B7AC710.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/C8D63E0C-5CAD-E711-BADA-E0071B73B6F0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/CA3903B6-5BAD-E711-B520-E0071B73C640.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/CA4CDA90-5FAD-E711-B469-E0071B7AC770.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/D212454B-62AD-E711-88E4-24BE05C488E1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/D2EFE07A-64AD-E711-A519-24BE05C46B01.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/D66E94B2-5AAD-E711-A1D9-E0071B7AE500.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/D805D90F-5CAD-E711-880D-E0071B73B6E0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/D84C0647-65AD-E711-BDE1-24BE05C4D851.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/D8F72B01-62AD-E711-9A72-E0071B73B6D0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/DA01579A-5FAD-E711-AD44-E0071B6C9DD0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/DA838637-66AD-E711-8319-24BE05CEED91.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/DA9F7E3E-5FAD-E711-8C82-5065F38152E1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/DC20086D-5BAD-E711-B590-24BE05BDBE21.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/DC9D0E27-5FAD-E711-BC04-24BE05C668E1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/DE9F7EEB-62AD-E711-B561-24BE05C616C1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/E043D73D-5DAD-E711-B9B6-24BE05C3CCE1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/E24FA7B2-5AAD-E711-A2AC-E0071B7AE500.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/E2F70E2B-5FAD-E711-8888-24BE05BD0F81.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/E4DD053D-5FAD-E711-B8F5-24BE05C4D8F1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/E4FB1C9B-5DAD-E711-9580-E0071B73C630.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/E85C527F-61AD-E711-BC1B-E0071B73B6D0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/E87D887D-61AD-E711-989D-E0071B7AC700.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/EC03A1C1-5CAD-E711-9D58-24BE05C4D821.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/EE0588AE-5CAD-E711-99AE-5065F38102F1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/EE1D5C80-62AD-E711-99BC-E0071B7A48A0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/F012B040-66AD-E711-A6A2-E0071B7A5650.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/F206FA89-5FAD-E711-8ED0-24BE05C4D8F1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/F250EE42-5EAD-E711-98EA-E0071B6C9DA0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/F457982A-5FAD-E711-8188-24BE05C4C891.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/F621BB90-64AD-E711-9427-E0071B6C9DD0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/F821A224-65AD-E711-923B-E0071B740D80.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/FC4FDE37-66AD-E711-B0FE-24BE05CEED91.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/FC542E2B-5FAD-E711-894A-24BE05C46B01.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/FCBB6F81-65AD-E711-A247-E0071B74AC10.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/FE2B8FBC-5CAD-E711-B4BA-24BE05CEDC81.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/FE3DE6A0-5FAD-E711-8D08-5065F38172A1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU140EA1000-v1/10000/FEF5C448-5DAD-E711-AD0A-24BE05C4D851.root' ] );

process.source = cms.Source("PoolSource",
    fileNames =  readFiles  )
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
process.GlobalTag.globaltag = '93X_upgrade2023_realistic_v2'

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
    #eeTimeDigi=cms.InputTag("mix", "EETimeDigi", "HLT"),
    #ebTimeDigi=cms.InputTag("mix", "EBTimeDigi", "HLT"),
    #superClusters = cms.InputTag("reducedEgamma", "reducedSuperClusters", "RECO"),

    #lostTracks = cms.InputTag("lostTracks", "", "RECO")
)

#run
process.p = cms.Path( #process.HBHENoiseFilterResultProducer*
                      #process.BadChargedCandidateFilter*
                      #process.BadPFMuonFilter*
                      process.ntuples)

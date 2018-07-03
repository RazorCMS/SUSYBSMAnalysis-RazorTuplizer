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
#	'/store/relval/CMSSW_9_0_0_pre4/RelValZMM_14/GEN-SIM-RECO/PU25ns_90X_upgrade2023_realistic_v3_D4TPU200c2-v1/10000/021BDAF8-0EF0-E611-81D0-0CC47A7C35F4.root'
#        '/store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/00D62D92-62AA-E711-994A-5065F38142E1.root'
#        '/store/relval/CMSSW_9_3_2/RelValH125GGgluonfusion_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17PU200EA1000-v1/10000/00932E58-6DAD-E711-BD20-5065F381B271.root'

       'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/00D62D92-62AA-E711-994A-5065F38142E1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/02110E2F-53AA-E711-BC98-24BE05BD4F81.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/029547F8-4EAA-E711-A90C-4C79BA180B95.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/02B84BB5-4CAA-E711-BD96-5065F381E201.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/04C0368A-50AA-E711-9BA5-4C79BA181331.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/04FAD1F7-53AA-E711-9CBD-E0071B73B6E0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/060ADDD3-50AA-E711-A8EC-24BE05C6D5A1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/0A9FD9D8-4CAA-E711-89EE-4C79BA320DCF.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/0C6514FB-4FAA-E711-8B0F-4C79BA18099D.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/0EF2FFFB-51AA-E711-9D4B-E0071B7A5870.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/1011F9E8-59AA-E711-A5C3-4C79BA18105F.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/1099F5B3-4BAA-E711-B5FF-E0071B7A5680.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/10B2F8E8-59AA-E711-BC98-4C79BA18105F.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/18485C7D-50AA-E711-B7F5-4C79BA320475.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/18ACFFFB-51AA-E711-BB42-E0071B7A5870.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/18E78535-4DAA-E711-A2CD-24BE05CEEDB1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/1A01E1C4-50AA-E711-B054-24BE05CEECD1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/1A1CA88B-52AA-E711-8B3D-E0071B7A68B0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/1A5664EC-55AA-E711-B36F-E0071B7A3830.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/1A6B8C18-51AA-E711-AB03-4C79BA18183D.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/1A8D28B7-4EAA-E711-9610-24BE05C33C81.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/1AC42705-52AA-E711-A35C-E0071B7A48A0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/1C3324FB-4DAA-E711-846C-4C79BA181847.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/1C7BF0D6-53AA-E711-9870-24BE05C3CCE1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/1E03760A-5CAA-E711-87E6-4C79BA181443.root',
#'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/1EAF3D5A-52AA-E711-9-E0071B74BDE0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/6600388F-52AA-E711-8EDB-E0071B7A28E0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/661700FC-51AA-E711-AF7E-E0071B7A5870.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/66860B2B-53AA-E711-B0A6-24BE05C61601.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/669F729D-51AA-E711-B0F3-4C79BA180A73.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/66E227B7-4EAA-E711-92F2-24BE05C33C81.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/6A0EA08D-52AA-E711-B95C-E0071B7A75A0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/6A4C22B8-4EAA-E711-827B-E0071B7A7640.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/6AAE3D4D-50AA-E711-9D87-24BE05C6E591.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/6C46DC66-53AA-E711-97CC-4C79BA1811A1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/6C733DE0-50AA-E711-BC8E-5065F38182E1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/6E4EC9C0-4DAA-E711-9ADD-4C79BA1810B9.root',
        'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/6E665304-4FAA-E711-8291-4C79BA320D8F.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/6E87C7C3-4FAA-E711-A7F8-5065F381C1D1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/6EF68DC5-4FAA-E711-B895-24BE05CE4D91.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/7085D7C6-50AA-E711-8E45-24BE05C61601.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/70A0CE39-4EAA-E711-B7A2-24BE05C6E7C1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/70DDA205-58AA-E711-9CAA-4C79BA18128D.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/7210C7C3-4FAA-E711-BCB1-5065F381C1D1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/741C89C5-50AA-E711-8D86-24BE05C6E7C1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/745484C5-4FAA-E711-A26A-24BE05CE4D91.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/74E3AAD4-53AA-E711-9627-24BE05BD0F81.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/76C7C7CD-52AA-E711-9C2F-4C79BA181211.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/781900FC-51AA-E711-86E5-E0071B7A5870.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/7842368A-50AA-E711-9983-4C79BA181331.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/7869D940-52AA-E711-9EAA-4C79BA3201DF.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/7890BB65-51AA-E711-8FD4-24BE05CE1E01.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/7A2AFB3D-52AA-E711-9A3D-4C79BA260295.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/7AC5C565-51AA-E711-9D5F-24BE05CE1E01.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/7E6E01D1-52AA-E711-A634-4C79BA320DCF.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/7E8A428E-51AA-E711-A763-4C79BA181357.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/8095770D-5EAA-E711-9AE5-4C79BA180929.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/82B251C5-52AA-E711-9F45-4C79BA181327.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/82D76A60-51AA-E711-AFBC-24BE05C6E7E1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/840DF93E-52AA-E711-8300-4C79BA180A01.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/84F56EC5-4FAA-E711-90AD-24BE05CE4D91.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/885B7804-52AA-E711-AE7E-E0071B7A48A0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/9043ED8F-52AA-E711-84A0-E0071B7AA680.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/90521BD9-53AA-E711-8635-24BE05CEFB41.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/90EF69F9-4CAA-E711-BC78-4C79BA1810B9.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/929D77C5-4FAA-E711-A818-24BE05CE4D91.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/985A1C57-51AA-E711-B017-24BE05C3EC61.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/98DCA0F3-53AA-E711-A7E9-E0071B7A9800.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/9A4F9316-51AA-E711-91ED-4C79BA180D3F.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/9A680B88-53AA-E711-AA1E-E0071B73C600.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/9AE28535-4DAA-E711-86B1-24BE05CEEDB1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/9E3722C4-55AA-E711-AAD8-E0071B7A3830.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/9E40674D-65AA-E711-8354-E0071B7B2380.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/A061E655-4CAA-E711-B3D7-E0071B695B81.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/A29046D4-50AA-E711-B06B-5065F381E251.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/A6B6643A-52AA-E711-A7A6-4C79BA181187.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/A6E9C565-51AA-E711-94D9-24BE05CE1E01.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/A848D940-52AA-E711-9A0E-4C79BA3201DF.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/A8CDD4FE-4FAA-E711-A58B-4C79BA180A9F.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/AA24A455-53AA-E711-9BEF-E0071B7A5650.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/AAAD21A0-52AA-E711-AAA3-4C79BA260295.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/ACA19316-51AA-E711-936F-4C79BA180D3F.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/AE4ECFBE-54AA-E711-9D81-E0071B7AC770.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/AE671CB9-4EAA-E711-B92A-E0071B7A7640.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/AEAF30AF-4CAA-E711-93F3-5065F38122A1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/AECE7905-52AA-E711-B20D-24BE05C6C741.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/B00C21A0-52AA-E711-9F5B-4C79BA260295.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/B07B47B7-5AAA-E711-AAB3-E0071B7AC750.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/B0A87234-4EAA-E711-8FF0-E0071B695B81.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/B41BB975-54AA-E711-A3EB-24BE05C6D5A1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/B6AC7A05-52AA-E711-918E-24BE05C6C741.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/B8E08C18-51AA-E711-9055-4C79BA18183D.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/BEA012FB-4FAA-E711-8605-4C79BA18099D.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/BEA147F8-4EAA-E711-B654-4C79BA180B95.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/C054678E-52AA-E711-B63E-E0071B7A25E0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/C06CFCD8-50AA-E711-8F86-5065F37D4131.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/C0BC649E-4BAA-E711-80BD-5065F3816201.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/C0C0A58E-52AA-E711-9DD9-E0071B7A25E0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/C0F74B8E-52AA-E711-9D57-E0071B7A75A0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/C22AFABB-55AA-E711-A7F0-E0071B740DA0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/C24E9300-4FAA-E711-8D83-4C79BA1811A1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/C2CBAB8E-4EAA-E711-92BC-E0071B7AC770.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/C4927935-4DAA-E711-80B1-24BE05CEEDB1.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/C4E2EC8F-52AA-E711-B00B-E0071B7AA680.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_9_3_2/RelValZMM_14/GEN-SIM-RECO/PU25ns_93X_upgrade2023_realistic_v2_2023D17TimingPU200-v1/10000/C6425934-5065F3816251.root'
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
    #superClusters = cms.InputTag("reducedEgamma", "reducedSuperClusters", "RECO"),

    #lostTracks = cms.InputTag("lostTracks", "", "RECO")
)

#run
process.p = cms.Path( #process.HBHENoiseFilterResultProducer*
                      #process.BadChargedCandidateFilter*
                      #process.BadPFMuonFilter*
                      process.ntuples)

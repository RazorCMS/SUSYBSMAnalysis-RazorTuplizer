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
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/0048030D-28A6-E611-8F01-0CC47A7C3638.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/00C58A95-22A6-E611-A502-0025905B85B6.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/02184C3C-24A6-E611-8452-0025905B855A.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/02D97C84-2BA6-E611-B1BF-0CC47A7C356A.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/04093083-24A6-E611-B0D1-0CC47A7C3444.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/04D4D179-1BA6-E611-A2F4-0CC47A4D767E.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/080CFEFD-38A6-E611-9FAD-0CC47A78A3EE.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/083EC81F-1BA6-E611-A502-0025905B85DC.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/083ED33E-23A6-E611-9109-0CC47A7C361E.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/0A5C2C1A-33A6-E611-8F12-0CC47A78A42E.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/0AA2F5D5-21A6-E611-A55A-0CC47A7C353E.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/0EE88101-25A6-E611-817B-0CC47A7C34C4.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/0EEAF622-28A6-E611-A6C7-0025905A608E.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/10C6592F-21A6-E611-B17E-0CC47A7C3604.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/125D9A30-2DA6-E611-8096-0025905A612A.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/14340F59-25A6-E611-99D9-0CC47A78A3F4.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/148A741C-1FA6-E611-AE26-003048FFCBB8.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/14B207D3-20A6-E611-8D9D-0CC47A7C35D2.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/18B1AAA6-29A6-E611-A220-0025905B857E.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/18B40FB2-1FA6-E611-B06C-0CC47A7C354C.root',   
	'/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/1A5705C1-25A6-E611-8A57-0CC47A7C34EE.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/1E6E7475-28A6-E611-BBC5-0CC47A4D76AA.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/1ED7C424-29A6-E611-A479-0CC47A78A468.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/1EE2CC70-1FA6-E611-88F3-0025905B8580.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/20BCE7D1-21A6-E611-85BE-0025905B85BA.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/222CA77E-1BA6-E611-8BC5-0CC47A4C8F06.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/22C43512-25A6-E611-9320-0025905B85F6.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/2448A4FE-25A6-E611-932E-0CC47A7C3428.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/2A552ED1-1FA6-E611-9787-0025905B85FC.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/2A933163-1AA6-E611-9BC1-0025905A6090.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/2C33DD89-28A6-E611-9AEC-0025905A609A.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/2CA903B4-25A6-E611-9836-0CC47A78A458.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/2EA25B1A-26A6-E611-B7ED-0CC47A7C357A.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/2ED0214A-25A6-E611-A19A-0CC47A7C3450.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/2EDF20EF-2BA6-E611-9574-0CC47A78A41C.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/32C6B30E-1BA6-E611-B5C6-0CC47A4D7690.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/34135B81-20A6-E611-B554-0CC47A4D7604.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/34B5590C-31A6-E611-912D-0CC47A7C35A4.root',
	'/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/34C81FFE-1AA6-E611-B153-0CC47A78A3EE.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/34F49668-26A6-E611-A771-0CC47A7C3430.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/36014B41-21A6-E611-9F3E-0CC47A7C3434.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/36824C2A-1FA6-E611-8127-0025905A48D0.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/36863B59-25A6-E611-AE50-0CC47A4D7654.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/386E0790-29A6-E611-8C57-0CC47A7C345C.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/3EDFB619-24A6-E611-AF87-0CC47A7C361E.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/402609DA-22A6-E611-86B1-0CC47A78A3B4.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/40E5DA19-20A6-E611-89DB-0CC47A7C3604.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/423D2412-1BA6-E611-BB1F-0CC47A4C8E82.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/426AC659-1DA6-E611-9180-0CC47A7C3472.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/44C30F41-19A6-E611-A4BC-0025905B85CC.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/44D00B53-53A6-E611-946F-0CC47A7C353E.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/463AB921-1BA6-E611-8B9F-0025905A6090.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/46A43D35-19A6-E611-9D78-0CC47A4C8E82.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/46D87B6E-1DA6-E611-9202-0CC47A7C357E.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/480108A4-1CA6-E611-8F60-0CC47A78A3F8.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/482A964E-1AA6-E611-B11E-0CC47A78A440.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/485C1AFC-24A6-E611-84F4-0CC47A78A414.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/48698ED7-22A6-E611-82D7-0CC47A7C3424.root',
	'/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/489A9C6E-42A6-E611-911A-0CC47A7C3612.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/48CB9D8F-1CA6-E611-ADEE-0CC47A7C3472.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/4AEC797A-28A6-E611-B0A4-0CC47A4C8F12.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/523F5A13-1DA6-E611-90BA-0CC47A7C3410.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/52BF0480-26A6-E611-80ED-0CC47A4D7626.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/5463FC54-1DA6-E611-8780-0CC47A7C3434.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/54C43C63-1DA6-E611-8F27-0CC47A4D7640.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/5691E321-1BA6-E611-9453-0025905B858C.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/5C6AC299-2BA6-E611-9C35-0CC47A7C3612.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/5E16CC97-20A6-E611-8935-0025905A6094.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/5E3EFEDE-1FA6-E611-AE31-0025905A48EC.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/623CB463-1FA6-E611-8077-0CC47A4C8E26.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/62B7565A-1FA6-E611-AB94-0CC47A7C357E.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/640E4237-2BA6-E611-ACEE-0CC47A7C34D0.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/645E9387-2FA6-E611-9413-0CC47A4D7628.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/649E1804-25A6-E611-A092-0CC47A7C34A6.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/668C07C1-83A6-E611-A54B-0CC47A78A45A.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/6698CC49-25A6-E611-8D90-0CC47A7C357E.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/66B59EC0-23A6-E611-BD74-0025905A612C.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/6882FA1B-20A6-E611-849E-0CC47A7C351E.root',
	'/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/6893B9A5-25A6-E611-A6C7-0CC47A7C340C.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/6A32E072-1DA6-E611-9C5F-0CC47A7C3636.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/6C4137FE-24A6-E611-8370-0CC47A78A456.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/6E4F7E37-21A6-E611-85C4-0CC47A7C345C.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/6EC09022-26A6-E611-A36E-0025905A6134.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/7041BEA8-2BA6-E611-9276-0CC47A4C8E34.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/72B72E18-24A6-E611-A578-0CC47A7C3444.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/748A791A-1BA6-E611-9C7A-0025905A60CA.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/7C01A8AF-1EA6-E611-A904-0025905B858C.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/7C443E47-21A6-E611-AA08-0CC47A7C3432.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/80A5DB6E-1FA6-E611-8E8C-0025905A611C.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/820206DF-2AA6-E611-9E9D-0CC47A7C347A.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/822A4820-1BA6-E611-A289-0025905B85DC.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/84913A7D-26A6-E611-9F43-0CC47A7C3410.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/86522294-1EA6-E611-B0FB-0CC47A4C8E46.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/88122193-22A6-E611-B05D-0025905B8590.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/88A9AFFA-26A6-E611-B307-0025905B8560.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/8A0280CA-21A6-E611-9371-0025905B858C.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/8A219F6F-2AA6-E611-BEE5-0CC47A78A3EE.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/8A60FDF2-22A6-E611-BB2F-0025905A60F8.root',
	'/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/90212264-25A6-E611-8A12-0025905B8586.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/92762468-2CA6-E611-905B-0CC47A7C3420.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/94C1253B-24A6-E611-ADCA-0025905B85A2.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/964316C0-25A6-E611-B242-0CC47A7C34E6.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/984908CC-1FA6-E611-82C0-0025905A497A.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/985AE086-26A6-E611-9E2E-0CC47A4C8E86.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/987BAB68-26A6-E611-BA4D-0CC47A7C353E.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/98C1DB9F-19A6-E611-B015-0CC47A7C3404.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/9EB1B4AA-21A6-E611-9D75-0CC47A7C3628.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/A2A3B01B-28A6-E611-883C-0CC47A7C354C.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/A43A7C54-1CA6-E611-B4A5-0025905A48E4.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/A446C4A8-19A6-E611-AF2E-0025905A612A.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/A48C6839-30A6-E611-8B54-0CC47A7C3420.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/A6FB67EE-19A6-E611-AC42-0CC47A7C356A.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/A6FBBF37-23A6-E611-9EBC-0CC47A7C346E.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/AA9F50B0-1CA6-E611-9709-0CC47A4D7658.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/AADBDE88-24A6-E611-A976-0CC47A7C3450.root',
	'/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/AC4CC2AE-25A6-E611-8B62-0CC47A7C3444.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/AC8BD113-22A6-E611-ADBB-0CC47A7C35E0.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/AE56B86D-28A6-E611-9D78-0CC47A7C3636.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/AE62E8C7-23A6-E611-9719-0025905B85EC.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/B07A6EB6-21A6-E611-ABF8-0CC47A7C3412.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/B4A8C568-2CA6-E611-938A-0025905B859A.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/B637AEE8-19A6-E611-A20F-0CC47A4D768C.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/B8245379-1DA6-E611-AB2C-0025905B8586.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/B854CA9D-1AA6-E611-B940-0CC47A7C353E.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/B858042F-1EA6-E611-B560-0CC47A7C340C.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/BA30D59C-23A6-E611-9567-0CC47A78A468.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/BA334DD3-26A6-E611-9555-0CC47A7C3430.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/BA48998F-22A6-E611-A013-0025905B85DE.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/BA7EC456-27A6-E611-84DB-0CC47A4D75EC.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/BC0567B8-1CA6-E611-A949-0025905B8598.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/BC46C36B-26A6-E611-9E2F-0CC47A7C340C.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/C0884220-1BA6-E611-A908-0025905B8566.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/C0C7FA76-22A6-E611-81F5-0CC47A4D76AC.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/C40F1320-24A6-E611-8870-0CC47A7C3572.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/C47CDEFB-1EA6-E611-87EE-0CC47A4C8E8A.root',
	'/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/C63C2A85-28A6-E611-982E-0CC47A4C8F30.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/C68F64A5-25A6-E611-86CD-0CC47A7C3430.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/C8886D11-1FA6-E611-8FCB-0025905A6104.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/C88C6718-24A6-E611-A50F-0CC47A7C35A8.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/C89BDD38-24A6-E611-B7F7-0025905A6104.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/CC1E12B5-23A6-E611-9F53-0CC47A78A408.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/D2AA082E-22A6-E611-BE82-0025905A60D0.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/D2D0108D-26A6-E611-8706-0CC47A4D76B2.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/D6125624-1DA6-E611-8BB7-0025905A6090.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/D636D053-1DA6-E611-8685-0CC47A7C3430.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/D64A2CD4-22A6-E611-B4A4-0CC47A4D7654.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/D67D7F25-19A6-E611-ABC2-0CC47A7C357A.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/D688EB74-1BA6-E611-878D-0CC47A78A426.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/D8CB4021-20A6-E611-93D4-0CC47A7C34A6.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/DA790E79-22A6-E611-B1E2-0CC47A7C354A.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/DC51E334-2BA6-E611-B5AF-0CC47A7C3412.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/DC69E463-1FA6-E611-BDAA-0CC47A78A4BA.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/E0C0F811-20A6-E611-B442-0CC47A7C34A6.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/E297542B-1DA6-E611-B69C-003048FFD76C.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/E29FF027-1DA6-E611-943F-0025905B8592.root',
	'/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/E494AE41-1CA6-E611-8E4D-0CC47A4D764A.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/E6AE2E13-45A6-E611-9091-0025905B85EC.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/E6B050C4-1AA6-E611-A7C1-0025905A6060.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/E8EF2902-1DA6-E611-A103-0CC47A7C35F8.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/EA5ADC69-1DA6-E611-8F5F-0CC47A7C35B2.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/ECAB8154-21A6-E611-9BDA-0CC47A4D7638.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/EE3D6750-2EA6-E611-B523-0CC47A7C35F8.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/EE4F3E0C-20A6-E611-A4B6-0CC47A7C3604.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/EE4F715E-31A6-E611-A9CF-0CC47A7C3412.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/EE5F2B84-28A6-E611-9FE7-0025905A60E4.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/EEA11377-1FA6-E611-BB79-0025905B8564.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/EECE6EEE-24A6-E611-8C51-0CC47A7C3428.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/EEF46BA7-25A6-E611-ABA6-0CC47A7C35C8.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/F07E6640-24A6-E611-ADCA-003048FFCBB2.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/F0A7DD40-29A6-E611-9AA0-0CC47A7C3422.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/F0C9F6F3-26A6-E611-A406-0CC47A78A42C.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/F43C3D16-1FA6-E611-B35F-0025905B85F6.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/F4D69E9C-29A6-E611-85C0-0CC47A7C3424.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/F6DBED0E-1BA6-E611-AE7E-0CC47A4D7644.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/F81DF137-23A6-E611-A107-0CC47A4D76D0.root',
	'/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/F84FAA09-20A6-E611-91E6-0CC47A7C3408.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/FA548CF3-4CA6-E611-829E-0CC47A7C3628.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/FC0E981D-24A6-E611-A1F3-0CC47A7C345C.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/FE00BF66-1AA6-E611-A93B-0025905B8612.root',
       '/store/relval/CMSSW_8_1_0_pre15/RelValH125GGgluonfusion_13/GEN-SIM-RECO/PU25ns_81X_upgrade2023_realistic_v3_2023D3Timing13TeVPU200r4-v1/10000/FEC38F92-2DA6-E611-89A1-0025905B859A.root'
)
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

#TFileService for output 
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("razorNtuple_PU200_NoTiming.root"),
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

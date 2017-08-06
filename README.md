SUSYBSMAnalysis-RazorTuplizer
=============================

Razor ntuplizer for running over LHC Run 2 miniAOD compatible with CMSSW_8_0_X

-----------------------------------
Instructions for compiling in CMSSW
-----------------------------------

    cmsrel CMSSW_8_0_26_patch1
    cd CMSSW_8_0_26_patch1/src
    cmsenv
    git cms-merge-topic cms-met:METRecipe_80X_part2 -u
    git cms-merge-topic rafaellopesdesa:Regression80XEgammaAnalysis_v2
    git cms-merge-topic ikrav:egm_id_80X_v2
    git clone git@github.com:cms-jet/JetToolbox.git JMEAnalysis/JetToolbox -b jetToolbox_80X_V3
    git clone git@github.com:RazorCMS/SUSYBSMAnalysis-RazorTuplizer SUSYBSMAnalysis/RazorTuplizer
    cp SUSYBSMAnalysis/RazorTuplizer/data/Spring16_GeneralPurpose_V1/ RecoEgamma/ElectronIdentification/data/ -r
    cp SUSYBSMAnalysis/RazorTuplizer/data/Spring16_HZZ_V1/ RecoEgamma/ElectronIdentification/data/ -r
    scram b

---------------------    
Running the ntuplizer
---------------------

    cmsRun python/razorTuplizer_MC_reHLT.py

    
Before running, check python/razorTuplizer.py to make sure that the correct global tag is defined. (process.GlobalTag.globaltag = ...)

To run using CRAB3:

    source /cvmfs/cms.cern.ch/crab3/crab.sh
    crab submit -c crabConfigRazorTuplizer.py

---------------------------------------
Including new bad track MET filter
---------------------------------------

As of 26/06/2016 (80X miniAOD v2), need to first run

git cms-merge-topic -u cms-met:CMSSW_8_0_X-METFilterUpdate

to get new bad charged candidate and bad muon MET filters.

---------------------------------------
Including new bad track MET filter
---------------------------------------

As of 23/11/2016 (80X miniAOD v2), to pick up the Spring16 electron MVA,
need to run

git cms-merge-topic ikrav:egm_id_80X_v2
cp SUSYBSMAnalysis/RazorTuplizer/data/Spring16_GeneralPurpose_V1/ RecoEgamma/ElectronIdentification/data/ -r
cp SUSYBSMAnalysis/RazorTuplizer/data/Spring16_HZZ_V1/ RecoEgamma/ElectronIdentification/data/ -r


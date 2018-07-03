SUSYBSMAnalysis-RazorTuplizer
=============================

Razor ntuplizer for running over LHC Run 2 GEN-SIM-RECO compatible with CMSSW_9_3_2

-----------------------------------
Instructions for compiling in CMSSW
-----------------------------------

    cmsrel CMSSW_9_3_2
    cd CMSSW_9_3_2/src
    cmsenv
    git clone -b CMSSW_9_3_2 git@github.com:RazorCMS/SUSYBSMAnalysis-RazorTuplizer SUSYBSMAnalysis/RazorTuplizer
    scram b -j20

---------------------    
Running the ntuplizer
---------------------

	cmsRun python/razorTuplizer_MC_25ns_MiniAODV2.py

    
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


--------------------------------------------
Running the ntuplizer for photon corrections
--------------------------------------------

	cmsRun python/razorTuplizer_MC_25ns_MiniAODV2_PhoCorr.py

SUSYBSMAnalysis-RazorTuplizer
=============================

Razor ntuplizer for running over LHC Run 2 miniAOD

Instructions for compiling in CMSSW
--------------

    cmsrel CMSSW_7_4_7
    cd CMSSW_7_4_7/src
    cmsenv
    git clone git@github.com:RazorCMS/SUSYBSMAnalysis-RazorTuplizer SUSYBSMAnalysis/RazorTuplizer
    scram b
    
Running the ntuplizer
--------------

    cmsRun python/razorTuplizer.py
    
Before running, check python/razorTuplizer.py to make sure that the correct global tag is defined. (process.GlobalTag.globaltag = ...)

To run using CRAB3:

    source /cvmfs/cms.cern.ch/crab3/crab.sh
    crab submit -c crabConfigRazorTuplizer.py

SUSYBSMAnalysis-RazorTuplizer
=============================

Razor ntuplizer for running over LHC Run 2 miniAOD

Instructions for compiling in CMSSW
--------------

    cmsrel CMSSW_7_2_0_pre6
    cd CMSSW_7_2_0_pre6/src
    git clone git@github.com:RazorCMS/SUSYBSMAnalysis-RazorTuplizer SUSYBSMAnalysis/RazorTuplizer
    scram b

Running the ntuplizer
--------------

    cmsRun python/razorTuplizer.py

To run using CRAB3:

    source /cvmfs/cms.cern.ch/crab3/crab.sh
    crab submit -c crabConfigRazorTuplizer.py

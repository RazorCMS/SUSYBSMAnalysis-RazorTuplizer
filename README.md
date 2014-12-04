SUSYBSMAnalysis-RazorTuplizer
=============================

Razor ntuplizer for running over LHC Run 2 miniAOD

Instructions for compiling in CMSSW
--------------

    cmsrel CMSSW_7_2_0
    cd CMSSW_7_2_0/src
    cmsenv
    git cms-merge-topic sixie:CMSSWTagsForRazorNtupler_V1.3
    git clone git@github.com:RazorCMS/SUSYBSMAnalysis-RazorTuplizer SUSYBSMAnalysis/RazorTuplizer
    scram b
    
For tags corresponding to V1.2 or earlier, you must use CMSSW_7_0_6_patch1. For tags corresponding to V1.3 or later, you can use CMSSW 72X.

You must run the cms-merge-topic command before cloning the RazorTuplizer code, because it requires an empty CMSSW src directory to work.

Running the BASE ntuplizer
--------------

    cmsRun python/razorTuplizer.py
    
Before running, check python/razorTuplizer.py to make sure that the correct global tag is defined. (process.GlobalTag.globaltag = ...)

To run using CRAB3:

    source /cvmfs/cms.cern.ch/crab3/crab.sh
    crab submit -c crabConfigRazorTuplizer.py

Running custom ntuplizer Example
--------------

    cmsRun python/razorAnaTuple.py

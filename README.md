B1;2cSUSYBSMAnalysis-RazorTuplizer
=============================

Razor ntuplizer for running over LHC Run 2 miniAOD

Instructions for compiling in CMSSW
--------------

    cmsrel CMSSW_7_2_0_pre6
    cd CMSSW_7_2_0_pre6/src
    git clone git@github.com:RazorCMS/SUSYBSMAnalysis-RazorTuplizer SUSYBSMAnalysis/RazorTuplizer
    scram b
    
Use CMSSW_7_0_6_patch1 instead of CMSSW_7_2_0_pre6 for now, because 72X does not appear to be compatible with the 70X miniAOD samples.

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

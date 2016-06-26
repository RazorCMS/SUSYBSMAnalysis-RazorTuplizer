SUSYBSMAnalysis-RazorTuplizer
=============================

Razor ntuplizer for running over LHC Run 2 miniAOD compatible with CMSSW_7_6_X

---------------------------------------
!!!! Including new bad track MET filter
 ---------------------------------------

As of 26/06/2016 (80X miniAOD v2), need to first run

git cms-merge-topic -u cms-met:CMSSW_8_0_X-METFilterUpdate

to get new bad charged candidate and bad muon MET filters.

-----------------------------------
Instructions for compiling in CMSSW
-----------------------------------

    cmsrel CMSSW_7_6_3
    cd CMSSW_7_6_3/src
    cmsenv
    git clone git@github.com:RazorCMS/SUSYBSMAnalysis-RazorTuplizer SUSYBSMAnalysis/RazorTuplizer
    scram b

---------------------    
Running the ntuplizer
---------------------

	cmsRun python/razorTuplizer_MC_25ns_MiniAODV2.py

    
Before running, check python/razorTuplizer.py to make sure that the correct global tag is defined. (process.GlobalTag.globaltag = ...)

To run using CRAB3:

    source /cvmfs/cms.cern.ch/crab3/crab.sh
    crab submit -c crabConfigRazorTuplizer.py



--------------------------------------------------------
In order to use the Energy smearing and scale correction
--------------------------------------------------------

https://twiki.cern.ch/twiki/bin/viewauth/CMS/EGMSmearer


	cmsrel CMSSW_7_6_3_patch2
	cd CMSSW_7_6_3_patch2/src
	cmsenv
	git cms-merge-topic -u matteosan1:smearer_76X
	git clone git@github.com:RazorCMS/SUSYBSMAnalysis-RazorTuplizer SUSYBSMAnalysis/RazorTuplizer
	scram b

--------------------------------------------
Running the ntuplizer for photon corrections
--------------------------------------------

	cmsRun python/razorTuplizer_MC_25ns_MiniAODV2_PhoCorr.py

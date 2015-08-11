#include <TFile.h>
#include "SUSYBSMAnalysis/RazorTuplizer/interface/EGammaMvaPhotonEstimator.h"
#include <cmath>
#include <vector>


//--------------------------------------------------------------------------------------------------
EGammaMvaPhotonEstimator::EGammaMvaPhotonEstimator() :
fisInitialized(kFALSE),
fMVAType(kPhotonMVATypeDefault),
fNMVABins(0)

{
  // Constructor.  
}

//--------------------------------------------------------------------------------------------------
EGammaMvaPhotonEstimator::~EGammaMvaPhotonEstimator()
{
  for (unsigned int i=0;i<fTMVAReader.size(); ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];
  }
}

//--------------------------------------------------------------------------------------------------
void EGammaMvaPhotonEstimator::initialize( std::vector<std::string> methodName,
					   std::vector<std::string> weightsfiles,
					   EGammaMvaPhotonEstimator::MVAType type				       
					   ) {

  //clean up first
  for (unsigned int i=0;i<fTMVAReader.size(); ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];    
  }
  fTMVAReader.clear();
  fMethodnames.clear();

  //initialize
  fisInitialized = kTRUE;
  fMVAType = type;

  //Define expected number of bins
  UInt_t ExpectedNBins = 0;
  if (type == kPhotonMVATypeDefault) {
    ExpectedNBins = 2;
  }

  fNMVABins = ExpectedNBins;
  
  //Check number of weight files given
  if (fNMVABins != weightsfiles.size() ) {
    std::cout << "Error: Expected Number of bins = " << fNMVABins << " does not equal to weightsfiles.size() = " 
              << weightsfiles.size() << std::endl; 
 
   #ifndef STANDALONE
    assert(fNMVABins == weightsfiles.size());
   #endif 
  }
  //Check number of weight files given
  if (methodName.size() != weightsfiles.size() ) {
    std::cout << "Error: Size mismatch between methodName vector and weightsfiles vector \n";
   #ifndef STANDALONE
    assert(methodName.size() == weightsfiles.size());
   #endif 
  }

  //Loop over all bins
  for (unsigned int i=0;i<fNMVABins; ++i) {

    TMVA::Reader *tmpTMVAReader = new TMVA::Reader( "!Color:!Silent:Error" );  
    tmpTMVAReader->SetVerbose(kTRUE);
  
    if (type == kPhotonMVATypeDefault) {
            
      // Add all the vars, we take the string with variable name from the weights file (the Expression field)
      tmpTMVAReader->AddVariable("recoPhi"   , &fMVAVar_Phi);
      tmpTMVAReader->AddVariable("r9"        , &fMVAVar_R9);
      tmpTMVAReader->AddVariable("sieieFull5x5", &fMVAVar_Sieie);
      tmpTMVAReader->AddVariable("sieipFull5x5", &fMVAVar_Sieip);
      tmpTMVAReader->AddVariable("e1x3Full5x5/e5x5Full5x5"        , &fMVAVar_E1x3OverE5x5);
      tmpTMVAReader->AddVariable("e2x2Full5x5/e5x5Full5x5"        , &fMVAVar_E2x2OverE5x5);
      tmpTMVAReader->AddVariable("e2x5Full5x5/e5x5Full5x5"        , &fMVAVar_E2x5MaxOverE5x5);
      tmpTMVAReader->AddVariable("recoSCEta" , &fMVAVar_Sceta);
      tmpTMVAReader->AddVariable("rawE"      , &fMVAVar_RawE);
      tmpTMVAReader->AddVariable("scEtaWidth", &fMVAVar_ScEtaWidth);
      tmpTMVAReader->AddVariable("scPhiWidth", &fMVAVar_ScPhiWidth);
      if (i == 1) {
	tmpTMVAReader->AddVariable("esEn/rawE" , &fMVAVar_ESEnOverRawE);
	tmpTMVAReader->AddVariable("esRR"      , &fMVAVar_ESEffSigmaRR);   
      }
      tmpTMVAReader->AddVariable("rho"       , &fMVAVar_Rho);
      tmpTMVAReader->AddVariable("phoIsoRaw" , &fMVAVar_PhotonIsoRaw);
      tmpTMVAReader->AddVariable("chIsoRaw"  , &fMVAVar_ChargedIsoRaw);
      tmpTMVAReader->AddVariable("chWorstRaw", &fMVAVar_WorstVertexChargedIsoRaw);
      // Add spectators
      tmpTMVAReader->AddSpectator("recoPt" , &fMVAVar_pt);
      tmpTMVAReader->AddSpectator("recoEta", &fMVAVar_eta);
          
    } 

    tmpTMVAReader->BookMVA(methodName[i] , weightsfiles[i]);
    std::cout << "MVABin " << i << " : MethodName = " << methodName[i]
              << " , type == " << type << " , "
              << "Load weights file : " << weightsfiles[i] 
              << std::endl;
    fTMVAReader.push_back(tmpTMVAReader);
    fMethodnames.push_back(methodName[i]);
  }
  std::cout << "Photon ID MVA Initialization Completed\n";

}


//--------------------------------------------------------------------------------------------------
UInt_t EGammaMvaPhotonEstimator::GetMVABin( double sceta) const {
  
    //Default is to return the first bin
    unsigned int bin = 0;

    if (fMVAType == EGammaMvaPhotonEstimator::kPhotonMVATypeDefault ) {
      bin = 0;
      if (sceta < 1.479 ) bin = 0;
      else bin = 1;
    } else {
      cout << "Error: Photon MVA Type provided = " << fMVAType << " is not implemented\n";
    }

    return bin;
}





//--------------------------------------------------------------------------------------------------

Double_t EGammaMvaPhotonEstimator::mvaValue(const pat::Photon& pho, 				
					    double varRho,
					    double varPhotonIsoRaw,
					    double varChargedIsoRaw,
					    double varWorstVertexChargedIsoRaw,
					    noZS::EcalClusterLazyTools *myEcalClusterLazyTool,
					    bool printDebug) {
  
  if (!fisInitialized) { 
    std::cout << "Error: EGammaMvaPhotonEstimator not properly initialized.\n"; 
    return -9999;
  }

  if ( fMVAType != EGammaMvaPhotonEstimator::kPhotonMVATypeDefault) {
    std::cout << "Error: This method should be called for MVA Type kPhotonMVATypeDefault only" << endl;
    return -9999;
  }
  
  //Compute cluster shape variables
  float see = -999;
  std::vector<float> vCov = myEcalClusterLazyTool->localCovariances( *(pho.superCluster()->seed()) );
  see = (isnan(vCov[0]) ? 0. : sqrt(vCov[0]));
  float sep = vCov[1];
  double E5x5 = myEcalClusterLazyTool-> e5x5(*(pho.superCluster()->seed()));
  

  fMVAVar_Phi = pho.phi();
  fMVAVar_R9 = pho.r9();
  fMVAVar_Sieie = see;
  fMVAVar_Sieip = sep;
  if (E5x5 > 0) {
    fMVAVar_E1x3OverE5x5 = myEcalClusterLazyTool-> e1x3(*(pho.superCluster()->seed())) / E5x5;
    fMVAVar_E2x2OverE5x5 = myEcalClusterLazyTool-> e2x2(*(pho.superCluster()->seed())) / E5x5;
    fMVAVar_E2x5MaxOverE5x5 = myEcalClusterLazyTool-> e2x5Max(*(pho.superCluster()->seed())) / E5x5;
  } else {
    fMVAVar_E1x3OverE5x5 = 0;
    fMVAVar_E2x2OverE5x5 = 0;
    fMVAVar_E2x5MaxOverE5x5 = 0;
  }
  fMVAVar_Sceta = pho.superCluster()->eta();
  fMVAVar_RawE = pho.superCluster()->rawEnergy();
  fMVAVar_ScEtaWidth = pho.superCluster()->etaWidth();
  fMVAVar_ScPhiWidth = pho.superCluster()->phiWidth();
  fMVAVar_Rho = varRho ;
  fMVAVar_PhotonIsoRaw = varPhotonIsoRaw;
  fMVAVar_ChargedIsoRaw = varChargedIsoRaw;
  fMVAVar_WorstVertexChargedIsoRaw = varWorstVertexChargedIsoRaw;
  fMVAVar_ESEnOverRawE = pho.superCluster()->preshowerEnergy() / pho.superCluster()->rawEnergy();
  fMVAVar_ESEffSigmaRR = myEcalClusterLazyTool->eseffsirir( *(pho.superCluster()) );

  fMVAVar_pt = pho.pt();
  fMVAVar_eta = pho.eta();
 
  // evaluate  
  Double_t mva = -9999;  
  mva = fTMVAReader[GetMVABin(fMVAVar_Sceta)]->EvaluateMVA(fMethodnames[GetMVABin(fMVAVar_Sceta)]);

  if(printDebug) {
    cout << " *** Photon MVA Evaluator : fMVAType " << fMVAType 
	 << " Bin = " << GetMVABin(fMVAVar_Sceta) << "\n";
    cout << "phi " << fMVAVar_Phi
	 << " r9 " << fMVAVar_R9
	 << " sieie " << fMVAVar_Sieie    
	 << " sieip " << fMVAVar_Sieip 
	 << " E1x3OverE5x5 " << fMVAVar_E1x3OverE5x5
	 << " E2x2OverE5x5 " << fMVAVar_E2x2OverE5x5
	 << " E2x5MaxOverE5x5 " << fMVAVar_E2x5MaxOverE5x5
	 << " Sceta " << fMVAVar_Sceta
	 << " RawE " << fMVAVar_RawE
	 << " ScEtaWidth " << fMVAVar_ScEtaWidth
	 << " ScPhiWidth " << fMVAVar_ScPhiWidth
	 << " Rho " << fMVAVar_Rho
	 << " PhotonIsoRaw " << fMVAVar_PhotonIsoRaw
	 << " ChargedIsoRaw "  << fMVAVar_ChargedIsoRaw
	 << " WorstVertexChargedIsoRaw " << fMVAVar_WorstVertexChargedIsoRaw
	 << " ESEnOverRawE " << fMVAVar_ESEnOverRawE
	 << " ESEffSigmaRR " << fMVAVar_ESEffSigmaRR    
	 << " pt " << fMVAVar_pt
	 << " eta " << fMVAVar_eta << "\n";
    cout << " ### MVA " << mva << "\n";   
  }

  return mva;
}





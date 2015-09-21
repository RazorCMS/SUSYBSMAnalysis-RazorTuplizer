#include <TFile.h>
#include "SUSYBSMAnalysis/RazorTuplizer/interface/ElectronMVAEstimatorRun2NonTrig.h"
#include <cmath>
#include <vector>


//--------------------------------------------------------------------------------------------------
ElectronMVAEstimatorRun2NonTrig::ElectronMVAEstimatorRun2NonTrig() :
fMethodname("BDTG method"),
fisInitialized(kFALSE),
fMVAType(kPHYS14)
{
  // Constructor.  
}

//--------------------------------------------------------------------------------------------------
ElectronMVAEstimatorRun2NonTrig::~ElectronMVAEstimatorRun2NonTrig()
{
  for (unsigned int i=0;i<fTMVAReader.size(); ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];
  }
}

//--------------------------------------------------------------------------------------------------
void ElectronMVAEstimatorRun2NonTrig::initialize( std::string methodName,
                                       std::string weightsfile,
                                       ElectronMVAEstimatorRun2NonTrig::MVAType type)
{
  
  std::vector<std::string> tempWeightFileVector;
  tempWeightFileVector.push_back(weightsfile);
  initialize(methodName,type,kFALSE,tempWeightFileVector);
}


//--------------------------------------------------------------------------------------------------
void ElectronMVAEstimatorRun2NonTrig::initialize( std::string methodName,
                                       ElectronMVAEstimatorRun2NonTrig::MVAType type,
                                       Bool_t useBinnedVersion,
				       std::vector<std::string> weightsfiles
  ) {

  //clean up first
  for (unsigned int i=0;i<fTMVAReader.size(); ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];
  }
  fTMVAReader.clear();

  //initialize
  fisInitialized = kTRUE;
  fMVAType = type;
  fMethodname = methodName;
 
  
  //Check number of weight files given
  if ( nCategories != weightsfiles.size() ) {
    std::cout << "Error: Expected Number of bins = " << nCategories << " does not equal to weightsfiles.size() = " 
              << weightsfiles.size() << std::endl; 
 
    assert(nCategories == weightsfiles.size());

  }

  //Loop over all bins
  for (unsigned int i=0;i<nCategories; ++i) {
  
    TMVA::Reader *tmpTMVAReader = new TMVA::Reader( "!Color:!Silent:Error" );  
    tmpTMVAReader->SetVerbose(kTRUE);

  
    // Pure ECAL -> shower shapes
    tmpTMVAReader->AddVariable("ele_oldsigmaietaieta", &fMVAVar_see);
    tmpTMVAReader->AddVariable("ele_oldsigmaiphiiphi", &fMVAVar_spp);
    tmpTMVAReader->AddVariable("ele_oldcircularity",   &fMVAVar_OneMinusE1x5E5x5);
    tmpTMVAReader->AddVariable("ele_oldr9",            &fMVAVar_R9);
    tmpTMVAReader->AddVariable("ele_scletawidth",      &fMVAVar_etawidth);
    tmpTMVAReader->AddVariable("ele_sclphiwidth",      &fMVAVar_phiwidth);
    tmpTMVAReader->AddVariable("ele_he",               &fMVAVar_HoE);
    // Endcap only variables
    if( i == CAT_EE_PT5to10 || i == CAT_EE_PT10plus ) 
      tmpTMVAReader->AddVariable("ele_psEoverEraw",    &fMVAVar_PreShowerOverRaw);
    
    //Pure tracking variables
    tmpTMVAReader->AddVariable("ele_kfhits",           &fMVAVar_kfhits);
    tmpTMVAReader->AddVariable("ele_kfchi2",           &fMVAVar_kfchi2);
    tmpTMVAReader->AddVariable("ele_gsfchi2",          &fMVAVar_gsfchi2);
    
    // Energy matching
    tmpTMVAReader->AddVariable("ele_fbrem",           &fMVAVar_fbrem);

    tmpTMVAReader->AddVariable("ele_gsfhits",         &fMVAVar_gsfhits);
    tmpTMVAReader->AddVariable("ele_expected_inner_hits",             &fMVAVar_expectedMissingInnerHits);
    tmpTMVAReader->AddVariable("ele_conversionVertexFitProbability",  &fMVAVar_convVtxFitProbability);

    tmpTMVAReader->AddVariable("ele_ep",              &fMVAVar_EoP);
    tmpTMVAReader->AddVariable("ele_eelepout",        &fMVAVar_eleEoPout);
    tmpTMVAReader->AddVariable("ele_IoEmIop",         &fMVAVar_IoEmIoP);
  
    // Geometrical matchings
    tmpTMVAReader->AddVariable("ele_deltaetain",      &fMVAVar_deta);
    tmpTMVAReader->AddVariable("ele_deltaphiin",      &fMVAVar_dphi);
    tmpTMVAReader->AddVariable("ele_deltaetaseed",    &fMVAVar_detacalo);
    
    // Spectator variables  
    tmpTMVAReader->AddSpectator("ele_pT",             &fMVAVar_pt);
    tmpTMVAReader->AddSpectator("ele_isbarrel",       &fMVAVar_isBarrel);
    tmpTMVAReader->AddSpectator("ele_isendcap",       &fMVAVar_isEndcap);
    tmpTMVAReader->AddSpectator("scl_eta",            &fMVAVar_SCeta);

    tmpTMVAReader->AddSpectator("ele_eClass",                 &fMVAVar_eClass);
    tmpTMVAReader->AddSpectator("ele_pfRelIso",               &fMVAVar_pfRelIso);
    tmpTMVAReader->AddSpectator("ele_expected_inner_hits",    &fMVAVar_expectedInnerHits);
    tmpTMVAReader->AddSpectator("ele_vtxconv",                &fMVAVar_vtxconv);
    tmpTMVAReader->AddSpectator("mc_event_weight",            &fMVAVar_mcEventWeight);
    tmpTMVAReader->AddSpectator("mc_ele_CBmatching_category", &fMVAVar_mcCBmatchingCategory);

    tmpTMVAReader->BookMVA(fMethodname , weightsfiles[i]);
    std::cout << "MVABin " << i << " : MethodName = " << fMethodname 
              << " , type == " << type << " , "
              << "Load weights file : " << weightsfiles[i] 
              << std::endl;
    fTMVAReader.push_back(tmpTMVAReader);
  }
  std::cout << "Electron ID MVA Completed\n";

}


//--------------------------------------------------------------------------------------------------
UInt_t ElectronMVAEstimatorRun2NonTrig::GetMVABin( double eta, double pt) const {
  
    //Default is to return the first bin
    unsigned int bin = 0;
    
    if (fMVAType == ElectronMVAEstimatorRun2NonTrig::kPHYS14 ) {
      bin = 0;
      if (pt < 10 && fabs(eta) < 0.8) bin = CAT_EB1_PT5to10;
      if (pt < 10 && fabs(eta) >= 0.8 && fabs(eta) < 1.479) bin = CAT_EB2_PT5to10;
      if (pt < 10 && fabs(eta) >= 1.479) bin = CAT_EE_PT5to10;
      if (pt >= 10 && fabs(eta) < 0.8) bin = CAT_EB1_PT10plus;
      if (pt >= 10 && fabs(eta) >= 0.8 && fabs(eta) < 1.479) bin = CAT_EB2_PT10plus;
      if (pt >= 10 && fabs(eta) >= 1.479) bin = CAT_EE_PT10plus;
    }

    return bin;
}



//--------------------------------------------------------------------------------------------------
//MVA Value
Double_t ElectronMVAEstimatorRun2NonTrig::mvaValue(const reco::GsfElectron& ele, 
						   const reco::Vertex& vertex, 
						   const TransientTrackBuilder& transientTrackBuilder,					
						   noZS::EcalClusterLazyTools myEcalCluster,
						   edm::Handle<vector<reco::Conversion> > conversions,
						   const math::XYZPoint beamspotPosition,
						   bool printDebug) {
  
  if (!fisInitialized) { 
    std::cout << "Error: ElectronMVAEstimatorRun2NonTrig not properly initialized.\n"; 
    return -9999;
  }
  
  bool validKF= false; 
  reco::TrackRef myTrackRef = ele.closestCtfTrackRef();
  validKF = (myTrackRef.isAvailable());
  validKF = (myTrackRef.isNonnull());  

  // Pure ECAL -> shower shapes
  std::vector<float> vCov = myEcalCluster.localCovariances(*(ele.superCluster()->seed())) ;
  fMVAVar_see = ele.full5x5_sigmaIetaIeta();
  fMVAVar_spp = ele.full5x5_sigmaIphiIphi();
  fMVAVar_OneMinusE1x5E5x5 = (ele.full5x5_e5x5()) !=0. ? 1.0 - ele.full5x5_e1x5()/ele.full5x5_e5x5() : -1.0 ;
  fMVAVar_R9              =  ele.full5x5_r9();
  fMVAVar_etawidth        =  ele.superCluster()->etaWidth();
  fMVAVar_phiwidth        =  ele.superCluster()->phiWidth();
  fMVAVar_HoE             =  ele.hadronicOverEm();
  //Endcap only variables
  fMVAVar_PreShowerOverRaw=  ele.superCluster()->preshowerEnergy() / ele.superCluster()->rawEnergy();

  // Pure tracking variables
  fMVAVar_kfhits          =  (validKF) ? myTrackRef->hitPattern().trackerLayersWithMeasurement() : -1. ; 
  fMVAVar_kfchi2          =  (validKF) ? myTrackRef->normalizedChi2() : 0 ;
  fMVAVar_gsfchi2         =  ele.gsfTrack()->normalizedChi2();  

  // Energy matching
  fMVAVar_fbrem           =  ele.fbrem();
  fMVAVar_gsfhits         =  ele.gsfTrack()->hitPattern().trackerLayersWithMeasurement();
  fMVAVar_expectedMissingInnerHits = ele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);

  reco::ConversionRef conv_ref = ConversionTools::matchedConversion(ele,conversions, beamspotPosition);
  double vertexFitProbability = -1.; 
  if(!conv_ref.isNull()) {
    const reco::Vertex &vtx = conv_ref.get()->conversionVertex(); if (vtx.isValid()) {
      vertexFitProbability = TMath::Prob( vtx.chi2(), vtx.ndof());
    } 
  }
  fMVAVar_convVtxFitProbability = vertexFitProbability;

  fMVAVar_EoP             =  ele.eSuperClusterOverP();
  fMVAVar_eleEoPout       =  ele.eEleClusterOverPout();
  float pAtVertex = ele.trackMomentumAtVtx().R();
  fMVAVar_IoEmIoP         =  (1.0/ele.ecalEnergy()) - (1.0 / pAtVertex);

  // Geometrical matchings
  fMVAVar_deta            =  ele.deltaEtaSuperClusterTrackAtVtx();
  fMVAVar_dphi            =  ele.deltaPhiSuperClusterTrackAtVtx();
  fMVAVar_detacalo        =  ele.deltaEtaSeedClusterTrackAtCalo();

  // Spectators
  fMVAVar_pt              =  ele.pt();                          
  constexpr float ebeeSplit = 1.479;
  fMVAVar_isBarrel        =  (ele.superCluster()->eta()<ebeeSplit);
  fMVAVar_isEndcap        =  (ele.superCluster()->eta()>=ebeeSplit);
  fMVAVar_SCeta           =  ele.superCluster()->eta();

  // The spectator variables below were examined for training, but
  // are not necessary for evaluating the discriminator, so they are
  // given dummy values (the specator variables above are also unimportant).
  // They are introduced only to match the definition of the discriminator 
  // in the weights file.
  constexpr unsigned nines = 999;
  fMVAVar_eClass               = nines;
  fMVAVar_pfRelIso             = nines;
  fMVAVar_expectedInnerHits    = nines;
  fMVAVar_vtxconv              = nines;
  fMVAVar_mcEventWeight        = nines;
  fMVAVar_mcCBmatchingCategory = nines;

  // evaluate
  bindVariables();
  Double_t mva = -9999;  
  mva = fTMVAReader[GetMVABin(fMVAVar_SCeta,fMVAVar_pt)]->EvaluateMVA(fMethodname);

  if(printDebug) {
    cout << " *** Inside the class fMethodname " << fMethodname << " fMVAType " << fMVAType << endl;
    cout << " bin " << GetMVABin(fMVAVar_SCeta,fMVAVar_pt) 
	 << " fbrem " <<  fMVAVar_fbrem  
      	 << " kfchi2 " << fMVAVar_kfchi2  
	 << " mykfhits " << fMVAVar_kfhits  
	 << " gsfchi2 " << fMVAVar_gsfchi2  
	 << " deta " <<  fMVAVar_deta  
	 << " dphi " << fMVAVar_dphi  
      	 << " detacalo " << fMVAVar_detacalo  
	 << " see " << fMVAVar_see  
	 << " spp " << fMVAVar_spp  
	 << " etawidth " << fMVAVar_etawidth  
	 << " phiwidth " << fMVAVar_phiwidth  
	 << " OneMinusE1x5E5x5 " << fMVAVar_OneMinusE1x5E5x5  
	 << " R9 " << fMVAVar_R9  
	 << " HoE " << fMVAVar_HoE  
	 << " EoP " << fMVAVar_EoP  
	 << " IoEmIoP " << fMVAVar_IoEmIoP  
	 << " eleEoPout " << fMVAVar_eleEoPout  
	 << " eta " << fMVAVar_SCeta  
	 << " pt " << fMVAVar_pt << endl;
    cout << " ### MVA " << mva << endl;
  }

  return mva;
}



Double_t ElectronMVAEstimatorRun2NonTrig::mvaValue(const pat::Electron& ele,
						   edm::Handle<vector<reco::Conversion> > conversions,
						   const math::XYZPoint beamspotPosition,
						   bool printDebug) {
    
    if (!fisInitialized) {
        std::cout << "Error: ElectronMVAEstimatorRun2NonTrig not properly initialized.\n";
        return -9999;
    }
    
    bool validKF= false;
    reco::TrackRef myTrackRef = ele.closestCtfTrackRef();
    validKF = (myTrackRef.isAvailable());
    validKF = (myTrackRef.isNonnull());
    
    // Pure ECAL -> shower shapes
    fMVAVar_see = ele.full5x5_sigmaIetaIeta(); 
    fMVAVar_spp = ele.full5x5_sigmaIphiIphi();  
    fMVAVar_OneMinusE1x5E5x5        =  (ele.full5x5_e5x5()) !=0. ? 1.0-(ele.full5x5_e1x5()/ele.full5x5_e5x5()) : -1.0 ;
    fMVAVar_R9              =  ele.full5x5_r9();
    fMVAVar_etawidth        =  ele.superCluster()->etaWidth();
    fMVAVar_phiwidth        =  ele.superCluster()->phiWidth();
    fMVAVar_HoE             =  ele.hadronicOverEm();
    //Endcap only variables
    fMVAVar_PreShowerOverRaw=  ele.superCluster()->preshowerEnergy() / ele.superCluster()->rawEnergy();

    // Pure tracking variables
    fMVAVar_kfhits          =  (validKF) ? myTrackRef->hitPattern().trackerLayersWithMeasurement() : -1. ;
    fMVAVar_kfchi2          =  (validKF) ? myTrackRef->normalizedChi2() : 0 ;
    fMVAVar_gsfchi2         =  ele.gsfTrack()->normalizedChi2();
    
    // Energy matching
    fMVAVar_fbrem           =  ele.fbrem();
    fMVAVar_gsfhits         =  ele.gsfTrack()->found();
    fMVAVar_expectedMissingInnerHits = ele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);

    reco::ConversionRef conv_ref = ConversionTools::matchedConversion(ele,conversions, beamspotPosition);
    double vertexFitProbability = -1.; 
    if(!conv_ref.isNull()) {
      const reco::Vertex &vtx = conv_ref.get()->conversionVertex(); if (vtx.isValid()) {
	vertexFitProbability = TMath::Prob( vtx.chi2(), vtx.ndof());
      } 
    }
    fMVAVar_convVtxFitProbability = vertexFitProbability;

    fMVAVar_EoP             =  ele.eSuperClusterOverP();
    fMVAVar_IoEmIoP         =  (1.0/ele.ecalEnergy()) - (1.0 / ele.p());  
    fMVAVar_eleEoPout       =  ele.eEleClusterOverPout();
       
    // Geometrical matchings
    fMVAVar_deta            =  ele.deltaEtaSuperClusterTrackAtVtx();
    fMVAVar_dphi            =  ele.deltaPhiSuperClusterTrackAtVtx();
    fMVAVar_detacalo        =  ele.deltaEtaSeedClusterTrackAtCalo();
               
    // Spectators
    fMVAVar_pt              =  ele.pt();
    fMVAVar_isBarrel        =  (ele.superCluster()->eta()<1.479);
    fMVAVar_isEndcap        =  (ele.superCluster()->eta()>=1.479);
     fMVAVar_SCeta             =  ele.superCluster()->eta();
           
     // The spectator variables below were examined for training, but
     // are not necessary for evaluating the discriminator, so they are
     // given dummy values (the specator variables above are also unimportant).
     // They are introduced only to match the definition of the discriminator 
     // in the weights file.
     fMVAVar_eClass               = 999;
     fMVAVar_pfRelIso             = 999;
     fMVAVar_expectedInnerHits    = 999;
     fMVAVar_vtxconv              = 999;
     fMVAVar_mcEventWeight        = 999;
     fMVAVar_mcCBmatchingCategory = 999;
     
     // evaluate
     bindVariables();
     Double_t mva = -9999;
     mva = fTMVAReader[GetMVABin(fMVAVar_SCeta,fMVAVar_pt)]->EvaluateMVA(fMethodname);
     
     if(printDebug) {
       cout << "Ele : " << ele.pt() << " " << ele.eta() << " " << ele.phi() << "\n";
       cout << " *** Inside the class fMethodname " << fMethodname << " fMVAType " << fMVAType << endl;
       cout << " fbrem " <<  fMVAVar_fbrem
	    << " kfchi2 " << fMVAVar_kfchi2
	    << " mykfhits " << fMVAVar_kfhits
	    << " gsfchi2 " << fMVAVar_gsfchi2
	    << " deta " <<  fMVAVar_deta
	    << " dphi " << fMVAVar_dphi
	    << " detacalo " << fMVAVar_detacalo
	    << " see " << fMVAVar_see
	    << " spp " << fMVAVar_spp
	    << " etawidth " << fMVAVar_etawidth  
	    << " phiwidth " << fMVAVar_phiwidth  
	    << " OneMinusE1x5E5x5 " << fMVAVar_OneMinusE1x5E5x5  
	    << " R9 " << fMVAVar_R9  
	    << " HoE " << fMVAVar_HoE  
	    << " EoP " << fMVAVar_EoP  
	    << " IoEmIoP " << fMVAVar_IoEmIoP  
	    << " eleEoPout " << fMVAVar_eleEoPout  
	    << " eta " << fMVAVar_SCeta
	    << " pt " << fMVAVar_pt << endl;
       cout << " ### MVA " << mva << endl;
     }
     
     return mva;
}



void ElectronMVAEstimatorRun2NonTrig::bindVariables() {

  // this binding is needed for variables that sometime diverge. 


  if(fMVAVar_fbrem < -1.)
    fMVAVar_fbrem = -1.;	
  
  fMVAVar_deta = fabs(fMVAVar_deta);
  if(fMVAVar_deta > 0.06)
    fMVAVar_deta = 0.06;
  
  fMVAVar_dphi = fabs(fMVAVar_dphi);
  if(fMVAVar_dphi > 0.6)
    fMVAVar_dphi = 0.6;
 
  if(fMVAVar_EoP > 20.)
    fMVAVar_EoP = 20.;
  
  if(fMVAVar_eleEoPout > 20.)
    fMVAVar_eleEoPout = 20.;
  
  fMVAVar_detacalo = fabs(fMVAVar_detacalo);
  if(fMVAVar_detacalo > 0.2)
    fMVAVar_detacalo = 0.2;
  
  if(fMVAVar_OneMinusE1x5E5x5 < -1.)
    fMVAVar_OneMinusE1x5E5x5 = -1;
  
  if(fMVAVar_OneMinusE1x5E5x5 > 2.)
    fMVAVar_OneMinusE1x5E5x5 = 2.; 
  
  if(fMVAVar_R9 > 5)
    fMVAVar_R9 = 5;
  
  if(fMVAVar_gsfchi2 > 200.)
    fMVAVar_gsfchi2 = 200;
  
  
  if(fMVAVar_kfchi2 > 10.)
    fMVAVar_kfchi2 = 10.;
    
  return;
}









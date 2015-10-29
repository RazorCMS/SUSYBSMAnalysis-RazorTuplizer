#include "RazorTuplizer.h"

//------ Auxiliary tools for RazorTuplizer class ------//

vector<TLorentzVector> RazorTuplizer::getHemispheres(vector<TLorentzVector> jets){
  int nJets = jets.size();
  vector<TLorentzVector> possibleHem1s; //holds possible hemisphere combinations
  vector<TLorentzVector> possibleHem2s;

  //stolen from https://github.com/pierinim/BSMatLHC/blob/master/BSMApp/src/CMS/CMSHemisphere.cc
  int nComb = pow(2, nJets);
  
  //step 1: store all possible partitions of the input jets
  int j_count;
  for(int i = 1; i < nComb-1; i++){ //note we omit the trivial hemisphere combinations (0 and nComb-1)
    TLorentzVector j_temp1, j_temp2;
    int itemp = i;
    j_count = nComb/2;
    int count = 0;
    while(j_count > 0){ //decompose i into binary '1's and '0's ; put the '1' jets into j_temp1 and the '0' jets into j_temp2
      if(itemp/j_count == 1){
	j_temp1 += jets[count];
      } else {
	j_temp2 += jets[count];
      }
      itemp -= j_count*(itemp/j_count); //note this is always (0 or 1)*j_count
      j_count /= 2;
      count++;
    }
    possibleHem1s.push_back(j_temp1);
    possibleHem2s.push_back(j_temp2);
  }
  
  //step 2: choose the partition that minimizes m1^2 + m2^2
  double mMin = -1;
  TLorentzVector myHem1;
  TLorentzVector myHem2;
  for(size_t i=0; i < possibleHem1s.size(); i++){
    double mTemp = possibleHem1s[i].M2() + possibleHem2s[i].M2();
    if(mMin < 0 || mTemp < mMin){
      mMin = mTemp;
      myHem1 = possibleHem1s[i];
      myHem2 = possibleHem2s[i];
    }
  }
  
  //return the hemispheres in decreasing order of pt
  vector<TLorentzVector> hemsOut;
  if(myHem1.Pt() > myHem2.Pt()){
    hemsOut.push_back(myHem1);
    hemsOut.push_back(myHem2);
  } else {
    hemsOut.push_back(myHem2);
    hemsOut.push_back(myHem1);
  }
  
  return hemsOut;
}

float RazorTuplizer::computeMR(TLorentzVector hem1, TLorentzVector hem2){
  return sqrt(pow(hem1.P() + hem2.P(), 2) - pow(hem1.Pz() + hem2.Pz(), 2));
}

float RazorTuplizer::computeR2(TLorentzVector hem1, TLorentzVector hem2, TLorentzVector pfMet){
  double mR = computeMR(hem1, hem2);
  double term1 = pfMet.Pt()/2*(hem1.Pt() + hem2.Pt());
  double term2 = pfMet.Px()/2*(hem1.Px() + hem2.Px()) + pfMet.Py()/2*(hem1.Py() + hem2.Py()); //dot product of MET with (p1T + p2T)
  double mTR = sqrt(term1 - term2);
  return (mTR / mR) * (mTR / mR);
}

bool RazorTuplizer::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle){
  //particle is already the ancestor
  if(ancestor == particle ) return true;
  
  //otherwise loop on mothers, if any and return true if the ancestor is found
  for(size_t i=0;i< particle->numberOfMothers();i++){
    if(isAncestor(ancestor,particle->mother(i))) return true;
  }
  //if we did not return yet, then particle and ancestor are not relatives
  return false;
}

const reco::Candidate* RazorTuplizer::findFirstMotherWithDifferentID(const reco::Candidate *particle){

  if( particle == 0 ){
    printf("ERROR! null candidate pointer, this should never happen\n");
    return 0;
  }

  // Is this the first parent with a different ID? If yes, return, otherwise
  // go deeper into recursion
  if (particle->numberOfMothers() > 0 && particle->pdgId() != 0) {
    if (particle->pdgId() == particle->mother(0)->pdgId()
	&& particle->mother(0)->status() != 11  // prevent infinite loop for sherpa documentation gluons
	) {
      return findFirstMotherWithDifferentID(particle->mother(0));
    } else {
      return particle->mother(0);
    }
  }

  return 0;
}

const reco::Candidate* RazorTuplizer::findOriginalMotherWithSameID(const reco::Candidate *particle){

  if( particle == 0 ){
    printf("ERROR! null candidate pointer, this should never happen\n");
    return 0;
  }

  // Is there another parent with the same ID? If yes, go deeper into recursion
  if (particle->numberOfMothers() > 0 && particle->pdgId() != 0) {    
    if (particle->mother(0)->numberOfMothers() == 0 || 
	particle->mother(0)->status() == 11 ||  // prevent infinite loop for sherpa documentation gluons
	(particle->mother(0)->numberOfMothers() > 0 && particle->mother(0)->mother(0)->pdgId() != particle->mother(0)->pdgId())
	) {
      return particle->mother(0);
    } else {
      return findOriginalMotherWithSameID(particle->mother(0));
    }
  }

  return 0;
}

//A copy from RecoEgamma/EgammaTools/src/ConversionTools.cc
//temporary solution because I'm not sure how to convert a PAT electron handle
//to a GSF electron handle
bool RazorTuplizer::hasMatchedPromptElectron(const reco::SuperClusterRef &sc, const edm::Handle<std::vector<pat::Electron> > &eleCol,
					const edm::Handle<reco::ConversionCollection> &convCol, const math::XYZPoint &beamspot, 
					float lxyMin, float probMin, unsigned int nHitsBeforeVtxMax) {

   //check if a given SuperCluster matches to at least one GsfElectron having zero expected inner hits
   //and not matching any conversion in the collection passing the quality cuts
 
   if (sc.isNull()) return false;
   
   for (std::vector<pat::Electron>::const_iterator it = eleCol->begin(); it!=eleCol->end(); ++it) {
     //match electron to supercluster
     if (it->superCluster()!=sc) continue;
 
     //check expected inner hits
     if (it->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 0) continue;
 
     //check if electron is matching to a conversion
     if (ConversionTools::hasMatchedConversion(*it,convCol,beamspot,lxyMin,probMin,nHitsBeforeVtxMax)) continue;
        
     return true;
   }
   
   return false;
  
}


tuple<double,double,double> RazorTuplizer::getPFMiniIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
					 const reco::Candidate* ptcl,
					 double r_iso_min, double r_iso_max, double kt_scale,
					 bool use_pfweight, bool charged_only) {
  
  //if (ptcl->pt()<5.) return 99999.;
  double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.);
  if(ptcl->isElectron()) {
    if (fabs(ptcl->eta())>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}
  } else if(ptcl->isMuon()) {
    deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;
  } else {
    //deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01; // maybe use muon cones??
  }
  double iso_nh(0.); double iso_ch(0.);
  double iso_ph(0.); double iso_pu(0.);
  double ptThresh(0.5);
  if(ptcl->isElectron()) ptThresh = 0;
  double r_iso = max(r_iso_min,min(r_iso_max, kt_scale/ptcl->pt()));
  for (const pat::PackedCandidate &pfc : *pfcands) {
    if (abs(pfc.pdgId())<7) continue;
    double dr = deltaR(pfc, *ptcl);
    if (dr > r_iso) continue;
    ////////////////// NEUTRALS /////////////////////////
    if (pfc.charge()==0){
      if (pfc.pt()>ptThresh) {
	double wpf(1.);
	if (use_pfweight){
	  double wpv(0.), wpu(0.);
	  for (const pat::PackedCandidate &jpfc : *pfcands) {
	    double jdr = deltaR(pfc, jpfc);
	    if (pfc.charge()!=0 || jdr<0.00001) continue;
	    double jpt = jpfc.pt();
	    if (pfc.fromPV()>1) wpv *= jpt/jdr;
	    else wpu *= jpt/jdr;
	  }
	  wpv = log(wpv);
	  wpu = log(wpu);
	  wpf = wpv/(wpv+wpu);
	}
	/////////// PHOTONS ////////////
	if (abs(pfc.pdgId())==22) {
	  if(dr < deadcone_ph) continue;
	  iso_ph += wpf*pfc.pt();
	  /////////// NEUTRAL HADRONS ////////////
	} else if (abs(pfc.pdgId())==130) {
	  if(dr < deadcone_nh) continue;
	  iso_nh += wpf*pfc.pt();
	}
      }
      ////////////////// CHARGED from PV /////////////////////////
    } else if (pfc.fromPV()>1){
      if (abs(pfc.pdgId())==211) {
	if(dr < deadcone_ch) continue;
	iso_ch += pfc.pt();
      }
      ////////////////// CHARGED from PU /////////////////////////
    } else {
      if (pfc.pt()>ptThresh){
	if(dr < deadcone_pu) continue;
	iso_pu += pfc.pt();
      }
    }
  }
  double iso(0.);
  if (charged_only){
    iso = iso_ch;
  } else {
    iso = iso_ph + iso_nh;
    if (!use_pfweight) iso -= 0.5*iso_pu;
    if (iso>0) iso += iso_ch;
    else iso = iso_ch;
  }
  iso = iso/ptcl->pt();

  //return pair of numbers for charged mini-iso and photon+neutralHadron mini-iso
  tuple<double,double,double> result;
  std::get<0>(result) =  iso_ch;
  std::get<1>(result) =  iso_ph + iso_nh;
  std::get<2>(result) =  iso_pu;

  return result;
}

double RazorTuplizer::ActivityPFMiniIsolationAnnulus(edm::Handle<pat::PackedCandidateCollection> pfcands,
						     const reco::Candidate* ptcl,
						     double dROuterSize,
						     double r_iso_min, double r_iso_max, double kt_scale) {  
  double iso_nh(0.); double iso_ch(0.);
  double iso_ph(0.); double iso_pu(0.);
  double ptThresh(0.5);
  if(ptcl->isElectron()) ptThresh = 0;
  double r_iso = max(r_iso_min,min(r_iso_max, kt_scale/ptcl->pt()));
  for (const pat::PackedCandidate &pfc : *pfcands) {
    if (abs(pfc.pdgId())<7) continue;
    double dr = deltaR(pfc, *ptcl);

    //select the annulus between mini-isolation cone size (r_iso) and dROuterSize
    if (!(dr < dROuterSize && dr >= r_iso)) continue;

    ////////////////// NEUTRALS /////////////////////////
    if (pfc.charge()==0){
      if (pfc.pt()>ptThresh) {

	/////////// PHOTONS ////////////
	if (abs(pfc.pdgId())==22) {
	  iso_ph += pfc.pt();
	  /////////// NEUTRAL HADRONS ////////////
	} else if (abs(pfc.pdgId())==130) {
	  iso_nh += pfc.pt();
	}
      }
      ////////////////// CHARGED from PV /////////////////////////
    } else if (pfc.fromPV()>1){
      if (abs(pfc.pdgId())==211) {
	iso_ch += pfc.pt();
      }
      ////////////////// CHARGED from PU /////////////////////////
    } else {
      if (pfc.pt()>ptThresh){
	iso_pu += pfc.pt();
      }
    }
  }
  double activity(0.);
  activity = iso_ch + fmax( 0.0, iso_ph + iso_nh - 0.5*iso_pu);

  return activity;
}

//**************************************************************
//Compute ptRel for leptons
//1) find closest jet
//2) subtract lepton from jet
//3) project lepton momentum perpendicular to closest jet
//**************************************************************
double RazorTuplizer::getLeptonPtRel(edm::Handle<pat::JetCollection> jets, const reco::Candidate* lepton) {

    const pat::Jet *closestJet = 0;
    double minDR = 9999;
    for (const pat::Jet &j : *jets) {
      if (j.pt() < 20) continue;
      double tmpDR = deltaR(j.eta(),j.phi(),lepton->eta(),lepton->phi());
      if (tmpDR < minDR) {
	minDR = tmpDR;
	closestJet = &j;
      }
    }

    //if no jet was found nearby, return some large default value
    if (!closestJet) return 9999;

    TLorentzVector closestJetFourVector(closestJet->px(),closestJet->py(),closestJet->pz(),closestJet->energy());    
    for (unsigned int i = 0, n = closestJet->numberOfSourceCandidatePtrs(); i < n; ++i) {
      
      const pat::PackedCandidate &candidate = dynamic_cast<const pat::PackedCandidate &>(*(closestJet->sourceCandidatePtr(i)));
      bool isPartOfLepton = false;

      if (lepton->isMuon()) {
	// muon candidate pointers to the PF candidate is null in miniAOD. 
	// we will match by relative pt difference and deltaR. thresholds at 0.1% and 0.001 in DR were tuned by eye
	if (abs(candidate.pdgId()) == 13 
	    && fabs(candidate.pt() - lepton->pt()) / lepton->pt() < 0.001
	    && deltaR(candidate.eta() , candidate.phi(), lepton->eta() , lepton->phi()) < 0.001
	    ) isPartOfLepton = true;
      }
      if (lepton->isElectron()) {
	for (auto itr : ((pat::Electron*)lepton)->associatedPackedPFCandidates()) {
	  if ( &(*itr) == &candidate) {
	    isPartOfLepton = true;
	    break;	  
	  }	
	}
      }
      //if the PF candidate is part of the muon, subtract its momentum from the jet momentum
      if (isPartOfLepton) {
	closestJetFourVector.SetPxPyPzE( closestJetFourVector.Px() - candidate.px(), 
					 closestJetFourVector.Py() - candidate.py(),
					 closestJetFourVector.Pz() - candidate.pz(),
					 closestJetFourVector.E() - candidate.energy());
      }
    }
    TLorentzVector lepFourVector(lepton->px(),lepton->py(),lepton->pz(),lepton->energy());    
    return lepFourVector.Perp(closestJetFourVector.Vect());
}

TLorentzVector RazorTuplizer::photonP4FromVtx( TVector3 vtx, TVector3 phoPos, double E )
{
  TVector3 p_hat;//Corrected photon direction
  p_hat = phoPos - vtx;
  TVector3 phoP3 = p_hat.Unit()*E;
  TLorentzVector phoP4;
  phoP4.SetVectM( phoP3, .0 );
  return phoP4;
};

bool RazorTuplizer::passJetID( const pat::Jet *jet, int cutLevel) {
  bool result = false;

  double NHF = jet->neutralHadronEnergyFraction();
  double NEMF = jet->neutralEmEnergyFraction();
  int NumConst = jet->chargedMultiplicity() + jet->neutralMultiplicity() ;
  double CHF = jet->chargedHadronEnergyFraction();
  double MUF = jet->muonEnergyFraction();
  double CEMF = jet->chargedEmEnergyFraction();
  int NumNeutralParticles =jet->neutralMultiplicity();
  int CHM = jet->chargedMultiplicity();

  //Loose
  if (cutLevel == 0) {
    if ( fabs(jet->eta()) <= 2.4) {
      if ( NHF  < 0.99 && NEMF < 0.99 && NumConst > 1 
	   && CHF > 0 && CHM > 0 && CEMF < 0.99 ) result = true;	   
    } else if( fabs(jet->eta()) <= 3.0)  {
      if ( NHF  < 0.99 && NEMF < 0.99 && NumConst > 1 ) result = true;	  
    } else {
      if ( NEMF < 0.90 && NumNeutralParticles > 10 ) result = true;	  
    }
  } 

  //Tight
  else if (cutLevel == 1) {
    if ( fabs(jet->eta()) <= 2.4) {
      if ( NHF  < 0.90 && NEMF < 0.90 && NumConst > 1 
	   && CHF > 0 && CHM > 0 && CEMF < 0.99 ) result = true;	   
    } else if( fabs(jet->eta()) <= 3.0)  {
      if ( NHF  < 0.90 && NEMF < 0.90 && NumConst > 1 ) result = true;	  
    } else {
      if ( NEMF < 0.90 && NumNeutralParticles > 10 ) result = true;	  
    }
  }

  //Tight Lep Veto
  else if (cutLevel == 2) {
    if ( fabs(jet->eta()) <= 2.4) {
      if ( NHF  < 0.90 && NEMF < 0.90 && NumConst > 1 
	   && CHF > 0 && CHM > 0 && CEMF < 0.99 && MUF < 0.8 ) result = true;	   
    } else if( fabs(jet->eta()) <= 3.0)  {
      if ( NHF  < 0.90 && NEMF < 0.90 && NumConst > 1 ) result = true;	  
    } else {
      if ( NEMF < 0.90 && NumNeutralParticles > 10 ) result = true;	  
    }
  }

  return result;
}

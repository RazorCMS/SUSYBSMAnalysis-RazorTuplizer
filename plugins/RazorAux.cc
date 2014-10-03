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



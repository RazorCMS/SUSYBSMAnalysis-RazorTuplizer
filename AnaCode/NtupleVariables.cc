//////////////////////////////////////////////////////////
// 
// Original Author : Irene Dutta
//                   IISER Pune
// Date Created    : Wed May 29, 2016
//////////////////////////////////////////////////////////
#define NtupleVariables_cxx

#include "NtupleVariables.h"
#include <TH2.h>
#include <TStyle.h>
#include<iostream>

double NtupleVariables::DeltaPhi(double phi1, double phi2) {
  double result = phi1 - phi2;
  while (result > M_PI)    result -= 2 * M_PI;
  while (result <= -M_PI)  result += 2 * M_PI;
  return result;
}

double NtupleVariables::DeltaR(double eta1, double phi1, double eta2, double phi2) {
  double deta = eta1 - eta2;
  double dphi = DeltaPhi(phi1, phi2);
  return std::sqrt(deta*deta + dphi*dphi);
}

#include "EvtGenModels/EvtGResPseudoTensor.hh"

namespace Gamapola{
  GResPseudoTensor::GResPseudoTensor(): GInterfaceForMathFunctions(),
  fCG_Kstr(-100),
  fCG_Rho(-100),
  fCoeffDelta(-100)
  {}
  GResPseudoTensor::~GResPseudoTensor()
  {
//     std::cout << "GCurrents destructor calling. . ." << std::endl;
  }
  
  void GResPseudoTensor::GInitializeResPseudoTensor(const int& charge, const std::string& resName)
  {
      fCG_Kstr = (charge==0)?sqrt(2.)/3.:-2./3.;
      fCG_Rho = (charge==0)?1./sqrt(3.):-1./sqrt(6.);
      fCoeffDelta = (charge==0)?1.:0.;
      if(resName == "K2_1600")
      {
          Dict1D<double> database;
          database["mass"] = kMK2_1600;
          database["width"] = kGammaK2_1600;
          database["left"] = std::pow((kMBmeson * kMBmeson - kMK2_1600 * kMK2_1600)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270),1.5);
          database["right"] = std::pow((kMBmeson * kMBmeson - kMK2_1600 * kMK2_1600)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270),1.5);
          fRes[resName] = database;
      }
      else if(resName == "K2_1770")
      {
          Dict1D<double> database;
          database["mass"] = kMK2_1770;
          database["width"] = kGammaK2_1770;
          database["left"] = std::pow((kMBmeson * kMBmeson - kMK2_1770 * kMK2_1770)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270),1.5);
          database["right"] = std::pow((kMBmeson * kMBmeson - kMK2_1770 * kMK2_1770)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270),1.5);
          fRes[resName] = database;
      }
  }
  
  std::complex<double> GResPseudoTensor::GKKstr(double s, double sij)
  {
      return std::complex<double>(0,1/std::sqrt(s/2.)) * fCG_Kstr * GBWKstr(sij,0);
  }
  std::complex<double> GResPseudoTensor::GKRho(double s, double sij)
  {
      return std::complex<double>(0,1/std::sqrt(s/2.)) * fCG_Rho * GBWRho(sij);
  }
  
  std::complex<double> GResPseudoTensor::GCKstr1(double s, double sKpi1, double fV, double hV)
  {
      return 2 * ( ( fV + hV * std::sqrt(s) * (fEnergy1 + fEnergy3) ) * (kMPion1*kMPion1 - kMKaon*kMKaon)/sKpi1 -
      hV * std::sqrt(s) * (fEnergy1 - fEnergy3) );
  }
  std::complex<double> GResPseudoTensor::GCKstr2(double s, double sKpi2, double fV, double hV)
  {
      return 2 * ( ( fV + hV * std::sqrt(s) * (fEnergy2 + fEnergy3) ) * (kMPion2*kMPion2 - kMKaon*kMKaon)/sKpi2 -
      hV * std::sqrt(s) * (fEnergy2 - fEnergy3) );
  }
  std::complex<double> GResPseudoTensor::GCRho(double s, double sPi1Pi2, double fV, double hV)
  {
//       std::cout << ( fV + hV * std::sqrt(s) * (fEnergy1 + fEnergy2) ) << "  " << fV << "  " << hV << '\n';
      return 2 * ( ( fV + hV * std::sqrt(s) * (fEnergy1 + fEnergy2) ) * (kMPion1*kMPion1 - kMPion2*kMPion2)/sPi1Pi2 -
      hV * std::sqrt(s) * (fEnergy1 - fEnergy2) );
  }
  
  std::complex<double> GResPseudoTensor::GC6( Dict1D<double>& pars, double s, double sKpi1, double sKpi2)
  {
      auto&& kPi1Pi2 = GKRho(s, fSij);
      auto&& kKPi1 = GKKstr(s, sKpi1);
      auto&& kKPi2 = GKKstr(s, sKpi2);
      
      auto&& cPi1Pi2 = GCRho(s, fSij, pars["fV_Rho"], pars["hV_Rho"]);
      auto&& cKpi1 = GCKstr1(s, sKpi1, pars["fV_K*"], pars["hV_K*"]);
      auto&& cKpi2 = GCKstr2(s, sKpi2, pars["fV_K*"], pars["hV_K*"]);

//       std::cout << cPi1Pi2 << "  " << kPi1Pi2 << '\n';
//       std::cout << kMPion1 << "  " << kMPion2 << '\n';
      auto&& coeffRhoPi1Pi2 = kPi1Pi2 * ( cPi1Pi2 * (fE0Pi + fE0Pj) - 2. * pars["fV_rho"] * fE0Pi );
      auto&& coeffKstrKpi1 = fCoeffDelta * kKPi1 * ( pars["fV_K*"] * fE0Pi + 2. * cKpi1 * (fE0Pi + fE0Pj) );
      auto&& coeffKstrKpi2 = kKPi2 * pars["fV_K*"] * fE0Pj;
      
      return coeffRhoPi1Pi2 + coeffKstrKpi1 + coeffKstrKpi2;
  }
  
  std::complex<double> GResPseudoTensor::GC7( Dict1D<double>& pars, double s, double sKpi1, double sKpi2)
  {
      fArg = 3.;
      auto&& kPi1Pi2 = GKRho(s, fSij);
      auto&& kKPi1 = GKKstr(s, sKpi1);
      auto&& kKPi2 = GKKstr(s, sKpi2);
      
      auto&& cPi1Pi2 = GCRho(s, fSij, pars["fV_Rho"], pars["hV_Rho"]);
      auto&& cKpi1 = GCKstr1(s, sKpi1, pars["fV_K*"], pars["hV_K*"]);
      auto&& cKpi2 = GCKstr2(s, sKpi2, pars["fV_K*"], pars["hV_K*"]);

//       std::cout << pars["fV_rho"] << "  "<< pars["hV_rho"] << '\n';
      auto&& coeffRhoPi1Pi2 = fCG_Rho * kPi1Pi2 * ( cPi1Pi2 * (fE0Pi + fE0Pj) + 2 * pars["fV_rho"] * fE0Pj );
      auto&& coeffKstrKpi1 = 2 * fCoeffDelta * kKPi1 * cKpi1 * pars["fV_K*"];
      auto&& coeffKstrKpi2 = kKPi2 * ( 2. * pars["fV_K*"] * ( fE0Pi + fE0Pj ) + cKpi2 * fE0Pj );
      
      return coeffRhoPi1Pi2 + coeffKstrKpi1 + coeffKstrKpi2;
  }
  
  std::vector<std::complex<double>> GResPseudoTensor::GMVec( Dict1D<double>& pars, double* x, int nRes, const std::string& resName)
  {
      auto s = x[4];
      auto sKpi1 = x[2];
      auto spi1pi2 = x[3];
      auto sKpi2 = GSij(s, sKpi1, spi1pi2);
      
      auto&& mass = fRes[resName]["mass"];
      auto&& width = fRes[resName]["width"];
      auto&& leftCoeff = fRes[resName]["left"];
      auto&& rightCoeff = fRes[resName]["right"];
    
      auto&& c6 = GC6(pars, s, sKpi1, sKpi2) * GBWK1bw(s, mass, width) * (kMBmeson*kMBmeson-s);
      auto&& c7 = GC7(pars, s, sKpi1, sKpi2) * GBWK1bw(s, mass, width) * (kMBmeson*kMBmeson-s);
      
//       std::cout << resName << "  " << mass << "  " << width << '\n';
      for(auto&& iAxis = 0; iAxis < 3; ++iAxis)
      {
          auto&& mVecCoord = c6 * fMomiV[iAxis] - c7 * fMomjV[iAxis];
          fJVec3L[nRes][iAxis] = mVecCoord * leftCoeff;
          fJVec3R[nRes][iAxis] = mVecCoord * rightCoeff;
      }
      
      return std::vector<std::complex<double>>{c6 * fMomiV[0] - c7 * fMomjV[0], 
                                               c6 * fMomiV[1] - c7 * fMomjV[1], 
                                               c6 * fMomiV[2] - c7 * fMomjV[2]};
  }
  
}

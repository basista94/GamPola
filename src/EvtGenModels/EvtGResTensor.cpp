#include "EvtGenModels/EvtGResTensor.hh"

namespace Gamapola{
  GResTensor::GResTensor(): GInterfaceForMathFunctions(),
  fCG_Kstr(-100),
  fCG_Rho(-100),
  fCoeffDelta(-100)
  {}
  GResTensor::~GResTensor()
  {
//     std::cout << "GCurrents destructor calling. . ." << std::endl;
  }
  
  void GResTensor::GInitializeResTensor(const int& charge, const std::string& resName)
  {
      fCG_Kstr = (charge==0)?sqrt(2.)/3.:-2./3.;
      fCG_Rho = (charge==0)?1./sqrt(3.):-1./sqrt(6.);
      fCoeffDelta = (charge==0)?1.:0.;
      if(resName == "K2_1430")
      {
          Dict1D<double> database;
          database["mass"] = kMK2_1430;
          database["width"] = kGammaK2_1430;
          database["left"] = std::pow((kMBmeson * kMBmeson - kMK2_1430 * kMK2_1430)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270),1.5);
          database["right"] = -std::pow((kMBmeson * kMBmeson - kMK2_1430 * kMK2_1430)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270),1.5);
          fRes[resName] = database;
      }
  }
  
  std::complex<double> GResTensor::GC4Kst(double* p, double s, double sij, double sik, int nRes)
  {
    std::complex<double> coeff = std::complex<double>(0, sqrt(2.)) * fgKstarKpi * sqrt(s) * fFF[nRes][0];
    std::complex<double> c4Kstr = fCG_Kstr * GBWKstr(sik,0) * coeff;
//     std::cout << "C4 gg5Kstr " << fFF[nRes][0] << "  " << fCCouplings[nRes][0] << std::endl;
//     std::cout << "c4Kst " << fCCouplings[nRes][0] * c4Kstr << '\n';
    cC[0][nRes][0][0][0] = -c4Kstr;
    
    return fCCouplings[nRes][0] * c4Kstr;
  }
  std::complex<double> GResTensor::GC5Kst(double* p, double s, double sij, double sik, int nRes)
  {
    std::complex<double> coeff = std::complex<double>(0, sqrt(2.)) * fgKstarKpi * sqrt(s) * fFF[nRes][0];
    std::complex<double> c5Kstr = fCG_Kstr * GBWKstr(sij,0) * coeff * fCoeffDelta;
//     std::cout << "C5 gg5Kstr " << fFF[nRes][0] << "  " << fCCouplings[nRes][0] << std::endl;
//     std::cout << "c5Kstr " << fCCouplings[nRes][0] * c5Kstr << '\n';
    cC[0][nRes][0][1][0] = c5Kstr;
    
    return fCCouplings[nRes][0] * c5Kstr;
  }
  std::complex<double> GResTensor::GC4Rho(double* p, double s, double sij, double sik, int nRes)
  {
    std::complex<double> coeff = std::complex<double>(0, sqrt(2.)) * fgRhoPiPi * sqrt(s) * fFF[nRes][1];
    std::complex<double> c4Rho = fCG_Rho * GBWRho(fSij) * coeff;
    cC[0][nRes][1][0][0] = c4Rho;
    cC[0][nRes][1][1][0] = c4Rho;
//     std::cout << "C4 gg5Rho " << fFF[nRes][1] << "  " << fCCouplings[nRes][1] << std::endl;
//     std::cout << "c4Rho " << fCCouplings[nRes][1] * c4Rho << '\n';
    return fCCouplings[nRes][1] * c4Rho;
  }  
  std::complex<double> GResTensor::GC4(double* p, double s, double sij, double sik, int nRes)
  {
    return (GC4Rho(p, s, sij, sik, nRes) - GC4Kst(p, s, sij, sik, nRes));
  }
  std::complex<double> GResTensor::GC5(double* p, double s, double sij, double sik, int nRes)
  {
    return (GC5Kst(p,s, sij, sik, nRes) + GC4Rho(p, s, sij, sik, nRes));
  }
  std::complex<double>* GResTensor::GKVec(double* p, double* x, int nRes, const std::string& resName)
  {
    
    auto s = x[4];
    auto sKpi1 = x[2];
    auto spi1pi2 = x[3];
    auto sKpi2 = GSij(s, sKpi1, spi1pi2);
    auto&& mass = fRes[resName]["mass"];
    auto&& width = fRes[resName]["width"];
    auto&& leftCoeff = fRes[resName]["left"];
    auto&& rightCoeff = fRes[resName]["right"];
//     std::cout << resName << "  " << mass << "  " << width << '\n';
      
//     std::cout << "****************************1430****************************************" << '\n';
    for(size_t iAxis = 0; iAxis < 3; iAxis++)
    {
      fJVec3[nRes][iAxis] = fDelta2[nRes] * (GC4(p,s, sKpi1, sKpi2, nRes) * fK2KinC4combo[iAxis] +
                                                     GC5(p,s, sKpi1, sKpi2, nRes) * fK2KinC5combo[iAxis] ) * GBWK1bw(s,mass,width) * (kMBmeson*kMBmeson-s)/(2*std::sqrt(2)*s);
      fJVec3L[nRes][iAxis] = fJVec3[nRes][iAxis] * leftCoeff;
      fJVec3R[nRes][iAxis] = fJVec3[nRes][iAxis] * rightCoeff;
      
//       std::cout << "delta2: " << fDelta2[nRes] << " current:  " << (GC4(p,s, sKpi1, sKpi2, nRes) * fK2KinC4combo[iAxis] +
//                                                      GC5(p,s, sKpi1, sKpi2, nRes) * fK2KinC5combo[iAxis] ) << " bw: " << GBWK1bw(s,mass,width) << " mb^2-s: " << (kMBmeson*kMBmeson-s)/(2*std::sqrt(2)*s) << '\n';
    }
    for(int jC = 0; jC < 2; jC++)
        {
          cRL[0][nRes][jC][0][0][0] = kPhaseSpaceFactor * cC[0][nRes][jC][0][0] * 
          (fEpsilonR[0] * fK2KinC4combo[0] + fEpsilonR[1] * fK2KinC4combo[1] + fEpsilonR[2] * fK2KinC4combo[2]) * fRight[nRes] / std::sqrt(2.) * GBWK1bw(s,mass,width) * (kMBmeson*kMBmeson-s)/sqrt(s)/(2*std::sqrt(2)*s);
        
          cRL[0][nRes][jC][1][0][0] = kPhaseSpaceFactor * cC[0][nRes][jC][1][0] * 
          (fEpsilonR[0] * fK2KinC5combo[0] + fEpsilonR[1] * fK2KinC5combo[1] + fEpsilonR[2] * fK2KinC5combo[2]) * fRight[nRes] / std::sqrt(2.) * GBWK1bw(s,mass,width) * (kMBmeson*kMBmeson-s)/sqrt(s)/(2*std::sqrt(2)*s);
          
          cRL[0][nRes][jC][0][1][0] = kPhaseSpaceFactor * cC[0][nRes][jC][0][0] * 
          (fEpsilonL[0] * fK2KinC4combo[0] + fEpsilonL[1] * fK2KinC4combo[1] + fEpsilonL[2] * fK2KinC4combo[2]) * fLeft[nRes] / std::sqrt(2.) * GBWK1bw(s,mass,width) * (kMBmeson*kMBmeson-s)/sqrt(s)/(2*std::sqrt(2)*s);
        
          cRL[0][nRes][jC][1][1][0] = kPhaseSpaceFactor * cC[0][nRes][jC][1][0] * 
          (fEpsilonL[0] * fK2KinC5combo[0] + fEpsilonL[1] * fK2KinC5combo[1] + fEpsilonL[2] * fK2KinC5combo[2]) * fLeft[nRes] / std::sqrt(2.) * GBWK1bw(s,mass,width) * (kMBmeson*kMBmeson-s)/sqrt(s)/(2*std::sqrt(2)*s);
//           std::cout << cRL[0][nRes][jC][0][0][0] << "   "  << 0 << "  " << 0 << "  " << jC << "  " << nRes << std::endl;
//           std::cout << cRL[0][nRes][jC][1][0][0] << "   "  << 1 << "  " << 0 << "  " << jC << "  " << nRes << std::endl;
//           std::cout << cRL[0][nRes][jC][0][1][0] << "   "  << 0 << "  " << 1 << "  " << jC << "  " << nRes << std::endl;
//           std::cout << cRL[0][nRes][jC][1][1][0] << "   "  << 1 << "  " << 1 << "  " << jC << "  " << nRes << std::endl;
        }
//        for(size_t iRL = 0; iRL < 2; iRL++) // {RL}
//           std::cout << fDelta2[nRes] * (cRL[0][nRes][0][0][iRL][0] + cRL[0][nRes][1][0][iRL][0] + 
//             cRL[0][nRes][0][1][iRL][0] + cRL[0][nRes][1][1][iRL][0]) << "  " << iRL << "  " << nRes << std::endl;   
//    std::cout << fDelta2[nRes] * (cRL[0][nRes][0][0][0][0] + cRL[0][nRes][1][0][0][0] + cRL[0][nRes][0][1][0][0] + cRL[0][nRes][1][1][0][0]) << 
//    "  " << fJVec3[nRes][charge][0] * fEpsilonR[0] + fJVec3[nRes][charge][1] * fEpsilonR[1] + fJVec3[nRes][charge][2] * fEpsilonR[2] << std::endl;
    return fJVec3[nRes];
  }
}

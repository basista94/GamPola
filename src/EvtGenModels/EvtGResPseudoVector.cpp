#include "EvtGenModels/EvtGResPseudoVector.hh"

namespace Gamapola{
  GResPseudoVector::GResPseudoVector(): GInterfaceForMathFunctions(),
  fCG_Kstr(-100),
  fCG_Rho(-100),
  fCoeffDelta(-100)
  {}

  GResPseudoVector::~GResPseudoVector()
  {
//     std::cout << "GCurrents destructor calling. . ." << std::endl;
  }
  
  void GResPseudoVector::GInitializeResPseudoVector(const int& charge, const std::string& resName)
  {
      fCG_Kstr = (charge==0)?sqrt(2.)/3.:-2./3.;
      fCG_Rho = (charge==0)?1./sqrt(3.):-1./sqrt(6.);
      fCoeffDelta = (charge==0)?1.:0.;
      if(resName == "K*_1410")
      {
          Dict1D<double> database;
          database["mass"] = kMKst_1410;
          database["width"] = kGammaKst_1410;
          database["left"] = std::pow((kMBmeson * kMBmeson - kMKst_1410 * kMKst_1410)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270),1.5);
          database["right"] = -std::pow((kMBmeson * kMBmeson - kMKst_1410 * kMKst_1410)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270),1.5);
          fRes[resName] = database;
      }
      else if(resName == "K*_1680")
      {
          Dict1D<double> database;
          database["mass"] = kMKst_1680;
          database["width"] = kGammaKst_1680;
          database["left"] = std::pow((kMBmeson * kMBmeson - kMKst_1680 * kMKst_1680)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270),1.5);
          database["right"] = -std::pow((kMBmeson * kMBmeson - kMKst_1680 * kMKst_1680)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270),1.5);
          fRes[resName] = database;
      }
  }
  
  std::complex<double> GResPseudoVector::GC3Kst(double* p, double s, double sij, double sik, int nRes)
  {
    std::complex<double> coeff = std::complex<double>(0,2.) * sqrt(s) * fFF[nRes][0];
    std::complex<double> c3Kstr = fCG_Kstr * fgKstarKpi * coeff * (GBWKstr(sik,0) + fCoeffDelta * GBWKstr(sij,0));
//     if(c3Kstr != c3Kstr)
//         std::cout << coeff << "  " << GBWKstr(sik,0) << "  " << fCoeffDelta << "  " << GBWKstr(sij,0) << '\n';
    cC[0][nRes][0][0][0] = c3Kstr;
//     std::cout << "C3 gg3Kstr " << fFF[nRes][0] << "   "<< fCCouplings[nRes][0] << std::endl;
//     std::cout << std::complex<double>(0,-2.) * sqrt(s) * fFF[nRes][0] * fgKstarKpi * GBWKstr(sij,0) * fCG_Kstr << std::endl;
//     std::cout << fCCouplings[nRes][0] * c3Kstr * fPiPjVec[0] * fEpsilonL[0] << std::endl;
//     std::cout << fCG_Kstr << "  " << fgKstarKpi << "   " << GBWKstr(sij,0) << "   " << fFF[nRes][0] << std::endl;
//     std::cout << fCCouplings[nRes][0] * c3Kstr << '\n';
    return fCCouplings[nRes][0] * c3Kstr;
  }
  std::complex<double> GResPseudoVector::GC3Rho(double* p, double s, double sij, double sik, int nRes)
  {
    std::complex<double> coeff = std::complex<double>(0,2.) * sqrt(s) * fFF[nRes][1];
    std::complex<double> c3Rho = fCG_Rho * fgRhoPiPi * GBWRho(fSij) * coeff;
    cC[0][nRes][1][0][0] = c3Rho;
//     std::cout << "C3 gg3Rho " << fFF[nRes][1] << "  " << fCCouplings[nRes][1] << std::endl;
//     std::cout << fCCouplings[nRes][1] * c3Rho * fPiPjVec[0] * fEpsilonL[0] << std::endl;
//     std::cout << fCCouplings[nRes][1] * c3Rho << '\n';
    return fCCouplings[nRes][1] * c3Rho;
  }
  std::complex<double> GResPseudoVector::GC3(double* p, double s, double sij, double sik, int nRes)
  {
    std::complex<double> currentC3 = (GC3Kst(p,s, sij, sik, nRes) + GC3Rho(p,s, sij, sik, nRes));
    return currentC3;
  }
  std::complex<double>* GResPseudoVector::GLVec(double* p, double* x, int nRes, const std::string& resName)
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
      
   for(size_t iAxis = 0; iAxis < 3; iAxis++)
   {
      fJVec3[nRes][iAxis] = fDelta2[nRes] * GC3(p,s, sKpi1, sKpi2, nRes) * fPiPjVec[iAxis] * GBWK1bw(s, mass, width);
      fJVec3L[nRes][iAxis] = fJVec3[nRes][iAxis] * leftCoeff;
      fJVec3R[nRes][iAxis] = fJVec3[nRes][iAxis] * rightCoeff;
//       if(fJVec3[nRes][charge][iAxis] != fJVec3[nRes][charge][iAxis])
//         std::cout << fDelta2[nRes] << "  " << fc3[nRes-2][charge] << "  " << fPiPjVec[iAxis] << '\n';
//       std::cout << fJVec3[nRes][charge][iAxis] << '\n';
   }
   for(int jC = 0; jC < 2; jC++)
   {
          
    cRL[0][nRes][jC][0][0][0] = kPhaseSpaceFactor * cC[0][nRes][jC][0][0] * (fEpsilonR[0] * fPiPjVec[0] + fEpsilonR[1] * fPiPjVec[1]
    + fEpsilonR[2] * fPiPjVec[2] ) * fRight[nRes] / std::sqrt(2.) * GBWK1bw(s, mass, width)/sqrt(s);
        
    cRL[0][nRes][jC][0][1][0] = kPhaseSpaceFactor * cC[0][nRes][jC][0][0] * (fEpsilonL[0] * fPiPjVec[0] + fEpsilonL[1] * fPiPjVec[1]
    + fEpsilonL[2] * fPiPjVec[2]) * fLeft[nRes] / std::sqrt(2.) * GBWK1bw(s, mass, width)/sqrt(s);
//     std::cout << cRL[0][nRes][jC][0][0][0] << "   "  << 0 << "  " << 0 << "  " << jC << "  " << nRes << std::endl;
//     std::cout << cRL[0][nRes][jC][0][1][0] << "   "  << 0 << "  " << 1 << "  " << jC << "  " << nRes << std::endl;
   }
//    std::cout << (cRL[0][nRes][0][0][1][0] + cRL[0][nRes][1][0][1][0]) * fDelta2[nRes] <<
//    "   " << 
//    fJVec3[nRes][charge][0] * fEpsilonL[0] + fJVec3[nRes][charge][1] * fEpsilonL[1] + fJVec3[nRes][charge][2] * fEpsilonL[2]  << std::endl;

//    for(size_t iCurr = 0; iCurr < 2; iCurr++) // C1,2
//         for(size_t iRL = 0; iRL < 2; iRL++) // {RL}
//           std::cout << fDelta2[nRes] * (cRL[0][nRes][0][0][iRL][0] + cRL[0][nRes][1][0][iRL][0]) << "  " << iRL << "  " << nRes << std::endl;
    return fJVec3[nRes];
  }
}

#include "EvtGenModels/EvtGResVector.hh"

namespace Gamapola{
  GResVector::GResVector(): GInterfaceForMathFunctions(),
  fCG_Kstr(-100),
  fCG_Rho(-100),
  fCoeffDelta(-100)
  {
  }
  GResVector::~GResVector()
  {
//     std::cout << "GCurrents destructor calling. . ." << std::endl;
  }
  
  void GResVector::GInitializeResVector(const int& charge, const std::string& resName)
  {
      fCG_Kstr = (charge==0)?sqrt(2.)/3.:-2./3.;
      fCG_Rho = (charge==0)?1./sqrt(3.):-1./sqrt(6.);
      fCoeffDelta = (charge==0)?1.:0.;
      if(resName == "K1_1270")
      {
          Dict1D<double> database;
          database["mass"] = kMK1_1270;
          database["width"] = kGammaK1_1270;
          database["left"] = 1.;
          database["right"] = 1.;
          fRes[resName] = database;
          std::cout << "1270" << '\n';
      }
      else if(resName == "K1_1400")
      {
          Dict1D<double> database;
          database["mass"] = kMK1_1400;
          database["width"] = kGammaK1_1400;
          database["left"] = (kMBmeson * kMBmeson - kMK1_1400 * kMK1_1400)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270);
          database["right"] = (kMBmeson * kMBmeson - kMK1_1400 * kMK1_1400)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270);
          fRes[resName] = database;
          std::cout << "1400" << '\n';
      }
  }
  
  std::complex<double> GResVector::GC1Kst(double* p, double s, double sij, double sik, int nRes)
  {
    std::complex<double> expPhaseSKstar(cos(p[3]), sin(p[3]) );
    std::complex<double> expPhaseSRho(cos(p[4]), sin(p[4]) );
    std::complex<double> expPhaseDRho(cos(p[5]), sin(p[5]) );
    std::complex<double> ff2(p[6], p[10]);
    double sSmallBrace = (1 + ( kMPion1 * kMPion1 - kMKaon * kMKaon ) / sik );
    double sBigBrace = ( 2 * fPiPjVec4 - sSmallBrace * ( sqrt(s) * fEnergy1 - kMPion1 * kMPion1 ) );
    std::complex<double> c1Kstr = fCG_Kstr * fgKstarKpi * 
    ( (fA2[nRes][0] * sSmallBrace + fB2[nRes][0] * sBigBrace) * GBWKstr(sik,1) - 2. * fA[nRes][0] * GBWKstr(sij,0) * fCoeffDelta );
//     if(c1Kstr != c1Kstr)
//     std::cout << fCG_Kstr[nRes] << "  " << fgKstarKpi << "  " << fA2[nRes][0] << "  " << sSmallBrace << "  " << fB2[nRes][0] << "  " <<
//     sBigBrace << "  " << GBWKstr(sik,1) << "  " << 2. * fA[nRes][0] << "  " << GBWKstr(sij,0) << "  " << fCoeffDelta[nRes] << '\n';
    
    for(int iphSKstr = 0; iphSKstr < 2; iphSKstr++)
    {
      cC[0][nRes][0][0][iphSKstr] = fCG_Kstr * fgKstarKpi * 
      ( (cF[0][nRes][0][0][1][iphSKstr] * sSmallBrace + cF[0][nRes][0][1][1][iphSKstr] * sBigBrace) * GBWKstr(sik,1)
      - 2. * cF[0][nRes][0][0][0][iphSKstr] * GBWKstr(sij,0) * fCoeffDelta );
    
      cC[1][nRes][0][0][iphSKstr] = fCG_Kstr * fgKstarKpi * 
      ( (cF[1][nRes][0][0][1][iphSKstr] * sSmallBrace + cF[1][nRes][0][1][1][iphSKstr] * sBigBrace) * GBWKstr(sik,1)
      - 2. * cF[1][nRes][0][0][0][iphSKstr] * GBWKstr(sij,0) * fCoeffDelta );
    }
//     std::cout << "C1 of K*" << c1Kstr << "  " << nRes << std::endl;
//     std::complex<double> c1Kstr = (cC[0][nRes][0][0][0] * sin(p[2]) + cC[1][nRes][0][0][0] * cos(p[2])) + 
//     (cC[0][nRes][0][0][1] * sin(p[2]) + cC[1][nRes][0][0][1] * cos(p[2])) * expPhaseSKstar;
//     std::cout << c1Kstr << '\n';
    return c1Kstr;
  }
  std::complex<double> GResVector::GC2Kst(double* p, double s, double sij, double sik, int nRes)
  {
    std::complex<double> expPhaseSKstar(cos(p[3]), sin(p[3]) );
    std::complex<double> expPhaseSRho(cos(p[4]), sin(p[4]) );
    std::complex<double> expPhaseDRho(cos(p[5]), sin(p[5]) );
    std::complex<double> ff2(p[6], p[10]);
    double sSmallBrace = (1 + ( kMPion2 * kMPion2 - kMKaon * kMKaon ) / sij );
    double sBigBrace = ( 2 * fPiPjVec4 - sSmallBrace * ( sqrt(s) * fEnergy2 - kMPion2 * kMPion2 ) );
    std::complex<double> c2Kstr = fCG_Kstr * fgKstarKpi * 
    ( (fA[nRes][0] * sSmallBrace + fB[nRes][0] * sBigBrace) * GBWKstr(sij,0) * fCoeffDelta
    - 2. * fA2[nRes][0] * GBWKstr(sik,1) );
//     if(c2Kstr != c2Kstr)
//       std::cout << fCG_Kstr[nRes] << "  " << fgKstarKpi << "  " << fA[nRes][0] << "  " << sSmallBrace << "  "  <<
//       fB[nRes][0] << "  " << sBigBrace << "  " << GBWKstr(sij,0) << "  " <<  fCoeffDelta[nRes] << '\n';
//     std::cout << c2Kstr << std::endl;
    for(int iphSKstr = 0; iphSKstr < 2; iphSKstr++)
    {
      cC[0][nRes][0][1][iphSKstr] = fCG_Kstr * fgKstarKpi * 
      ( (cF[0][nRes][0][0][0][iphSKstr] * sSmallBrace + cF[0][nRes][0][1][0][iphSKstr] * sBigBrace) * GBWKstr(sij,0) * fCoeffDelta
      - 2. * cF[0][nRes][0][0][1][iphSKstr] * GBWKstr(sik,1) );
      
      cC[1][nRes][0][1][iphSKstr] = fCG_Kstr * fgKstarKpi * 
      ( (cF[1][nRes][0][0][0][iphSKstr] * sSmallBrace + cF[1][nRes][0][1][0][iphSKstr] * sBigBrace) * GBWKstr(sij,0) * fCoeffDelta
      - 2. * cF[1][nRes][0][0][1][iphSKstr] * GBWKstr(sik,1) );
    }
//     std::cout << "C2 of K*" << c2Kstr << "  " << nRes << std::endl;
//     std::complex<double> c2Kstr = (cC[0][nRes][0][1][0] * sin(p[2]) + cC[1][nRes][0][1][0] * cos(p[2])) + 
//     (cC[0][nRes][0][1][1] * sin(p[2]) + cC[1][nRes][0][1][1] * cos(p[2])) * expPhaseSKstar;
//     std::cout << c2Kstr << '\n';
    return c2Kstr;
  }
  std::complex<double> GResVector::GC1Rho(double* p, double s, double sij, double sik, int nRes)
  {
    std::complex<double> expPhaseSKstar(cos(p[3]), sin(p[3]) );
    std::complex<double> expPhaseSRho(cos(p[4]), sin(p[4]) );
    std::complex<double> expPhaseDRho(cos(p[5]), sin(p[5]) );
    std::complex<double> ff2(p[6], p[10]);
    double sVerySmallBrace = sqrt( s ) * ( fEnergy1 - fEnergy2 );
    std::complex<double> c1Rho = fCG_Rho * fgRhoPiPi * ( fA[nRes][1] - fB[nRes][1] * sVerySmallBrace) * GBWRho(fSij);
//     if(c1Rho != c1Rho)
//       std::cout << fCG_Rho[nRes] << "  " << fgRhoPiPi << "  " << fA[nRes][1] << "  "<< fB[nRes][1] << "  " << sVerySmallBrace << "  " << GBWRho(fSij) << '\n';
    for(int iphRho = 0; iphRho < 2; iphRho++)
    {
      cC[0][nRes][1][0][iphRho] = fCG_Rho * fgRhoPiPi * ( cF[0][nRes][1][0][0][iphRho] - 
      cF[0][nRes][1][1][0][iphRho] * sVerySmallBrace) * GBWRho(fSij);
      
      cC[1][nRes][1][0][iphRho] = fCG_Rho * fgRhoPiPi * ( cF[1][nRes][1][0][0][iphRho] - 
      cF[1][nRes][1][1][0][iphRho] * sVerySmallBrace) * GBWRho(fSij);
    }
//     std::cout << "C1 of Rho" << c1Rho << "  " << nRes << std::endl;
//     std::complex<double> c1Rho = (cC[0][nRes][1][0][0] * sin(p[2]) + cC[1][nRes][1][0][0] * cos(p[2])) * expPhaseDRho + 
//     (cC[0][nRes][1][0][1] * sin(p[2]) + cC[1][nRes][1][0][1] * cos(p[2])) * expPhaseDRho * expPhaseSRho;
//     std::cout << c1Rho << '\n';
    return c1Rho;
  }
  std::complex<double> GResVector::GC2Rho(double* p, double s, double sij, double sik, int nRes)
  {
   std::complex<double> expPhaseSKstar(cos(p[3]), sin(p[3]) );
    std::complex<double> expPhaseSRho(cos(p[4]), sin(p[4]) );
    std::complex<double> expPhaseDRho(cos(p[5]), sin(p[5]) );
    std::complex<double> ff2(p[6], p[10]);
    double sVerySmallBrace = sqrt( s ) * ( fEnergy1 - fEnergy2 );
    std::complex<double> c2Rho = fCG_Rho * fgRhoPiPi * ( fA[nRes][1] + fB[nRes][1] * sVerySmallBrace) * GBWRho(fSij);
//     if(isnan(c2Rho.real()))
//       std::cout << "nan" << '\n';
//     std::cout << c2Rho << std::endl;
    for(int iphRho = 0; iphRho < 2; iphRho++)
    {
      cC[0][nRes][1][1][iphRho] = fCG_Rho * fgRhoPiPi * ( cF[0][nRes][1][0][0][iphRho] + 
      cF[0][nRes][1][1][0][iphRho] * sVerySmallBrace) * GBWRho(fSij);
      
      cC[1][nRes][1][1][iphRho] = fCG_Rho * fgRhoPiPi * ( cF[1][nRes][1][0][0][iphRho] + 
      cF[1][nRes][1][1][0][iphRho] * sVerySmallBrace) * GBWRho(fSij);
    }
//     std::complex<double> c2Rho = (cC[0][nRes][1][1][0] * sin(p[2]) + cC[1][nRes][1][1][0] * cos(p[2])) * expPhaseDRho + 
//     (cC[0][nRes][1][1][1] * sin(p[2]) + cC[1][nRes][1][1][1] * cos(p[2])) * expPhaseDRho * expPhaseSRho;
//     std::cout << "C2 of Rho" << c2Rho << "  " << nRes << std::endl;
//     std::cout << c2Rho << '\n';
    return c2Rho;
  }
  std::complex<double> GResVector::GC1Kappa(double* p, double s, double sij, double sik, int nRes)
  {
    std::complex<double> ff2(p[6], p[10]);
    std::complex<double> c1Kappa = -fCG_Kstr * 2 * fgK1KappaPi[nRes] * GBWKappa(sik);
    for(int iph = 0; iph < 2; iph++)
    {
      cC[0][nRes][2][0][iph] = c1Kappa;
      cC[1][nRes][2][0][iph] = c1Kappa;
    }
//     if(c1Kappa != c1Kappa)
//      std::cout << " C1Kappa " << c1Kappa << "  " << nRes << "  " << fCoeffDelta[nRes] << "   " << fCG_Kstr[nRes] << std::endl;
//     std::cout << ff2 * c1Kappa << '\n';
    return ff2 * c1Kappa;
  }
  std::complex<double> GResVector::GC2Kappa(double* p, double s, double sij, double sik, int nRes)
  {
    std::complex<double> ff2(p[6], p[10]);
    std::complex<double> c2Kappa = -fCG_Kstr * 2 * fgK1KappaPi[nRes] * GBWKappa(sij) * fCoeffDelta;
    for(int iph = 0; iph < 2; iph++)
    {
      cC[0][nRes][2][1][iph] = c2Kappa;
      cC[1][nRes][2][1][iph] = c2Kappa;
    }
//     std::cout << ff2 * c2Kappa << '\n';
//     if(c2Kappa != c2Kappa)
//       std::cout << " C2Kappa " << c2Kappa << "  " << nRes << "  " << fCoeffDelta[nRes] << "   " << fCG_Kstr[nRes] << std::endl;
    return ff2 * c2Kappa;
  }
  
  std::complex<double> GResVector::GC1(double* p, double s, double sij, double sik, int nRes)
  {
   std::complex<double> expPhaseSKstar(cos(p[3]), sin(p[3]) );
    std::complex<double> expPhaseSRho(cos(p[4]), sin(p[4]) );
    std::complex<double> expPhaseDRho(cos(p[5]), sin(p[5]) );
    std::complex<double> ff2(p[6], p[10]);
    std::complex<double> currentC1 = ( GC1Kst(p,s, sij, sik, nRes) + GC1Rho(p,s, sij, sik, nRes) + 
    GC1Kappa(p,s, sij, sik, nRes) );
    return currentC1;
  }
  std::complex<double> GResVector::GC2(double* p, double s, double sij, double sik, int nRes)
  {
    std::complex<double> expPhaseSKstar(cos(p[3]), sin(p[3]) );
    std::complex<double> expPhaseSRho(cos(p[4]), sin(p[4]) );
    std::complex<double> expPhaseDRho(cos(p[5]), sin(p[5]) );
    std::complex<double> ff2(p[6], p[10]);
    std::complex<double> currentC2 = ( GC2Kst(p,s, sij, sik, nRes) + GC2Rho(p,s, sij, sik, nRes) +
    GC2Kappa(p,s, sij, sik, nRes) );
    return currentC2;
  }
  std::complex<double>* GResVector::GJVec(double* p, double* x, int nRes, const std::string& resName)
  {
    auto s = x[4];
    auto sKpi1 = x[2];
    auto spi1pi2 = x[3];
    auto sKpi2 = GSij(s, sKpi1, spi1pi2);
    
    std::complex<double> expPhaseSKstar(cos(p[3]), sin(p[3]) );
    std::complex<double> expPhaseSRho(cos(p[4]), sin(p[4]) );
    std::complex<double> expPhaseDRho(cos(p[5]), sin(p[5]) );
    std::complex<double> ff2(p[6], p[10]);
    
    auto&& mass = fRes[resName]["mass"];
    auto&& width = fRes[resName]["width"];
    auto&& leftCoeff = fRes[resName]["left"];
    auto&& rightCoeff = fRes[resName]["right"];
//     std::cout << resName << "  " << mass << "  " << width << '\n';
//     std::cout << "*************************1+************************************" << '\n';
    for(size_t iAxis = 0; iAxis < 3; iAxis++)
    {
      fJVec3[nRes][iAxis] = fDelta2[nRes] * (GC1(p,s, sKpi1, sKpi2, nRes) * fMomiV[iAxis] - GC2(p,s, sKpi1, sKpi2, nRes) * fMomjV[iAxis]) * GBWK1bw(s,mass,width); 
      fJVec3L[nRes][iAxis] = fJVec3[nRes][iAxis] * leftCoeff;
      fJVec3R[nRes][iAxis] = fJVec3[nRes][iAxis] * rightCoeff;
      
//       std::cout << "delta2: " << fDelta2[nRes] << " current:  " << (GC1(p,s, sKpi1, sKpi2, nRes) * fMomiV[iAxis] - GC2(p,s, sKpi1, sKpi2, nRes) * fMomjV[iAxis]) << " bw: " << GBWK1bw(s,mass,width) << '\n';
    }
//     std::cout << "In vector: " << fJVec3L[nRes][0] << " " << fJVec3R[nRes][0] << " " << 0 << '\n';
    for(int iTrig = 0; iTrig < 2; iTrig++)
      for(int jC = 0; jC < 3; jC++)
        for(int iPh = 0; iPh < 2; iPh++)
        {
          
          cRL[iTrig][nRes][jC][0][0][iPh] = kPhaseSpaceFactor * (fEpsilonR[0] * fMomiV[0] * cC[iTrig][nRes][jC][0][iPh]
          + fEpsilonR[1] * fMomiV[1] * cC[iTrig][nRes][jC][0][iPh]
          + fEpsilonR[2] * fMomiV[2] * cC[iTrig][nRes][jC][0][iPh]) * fRight[nRes] / std::sqrt(2.) * GBWK1bw(s,mass,width)/sqrt(s);
        
          cRL[iTrig][nRes][jC][1][0][iPh] = kPhaseSpaceFactor * (-fEpsilonR[0] * fMomjV[0] * cC[iTrig][nRes][jC][1][iPh]
          - fEpsilonR[1] * fMomjV[1] * cC[iTrig][nRes][jC][1][iPh]
          - fEpsilonR[2] * fMomjV[2] * cC[iTrig][nRes][jC][1][iPh]) * fRight[nRes] / std::sqrt(2.) * GBWK1bw(s,mass,width)/sqrt(s);
          
          cRL[iTrig][nRes][jC][0][1][iPh] = kPhaseSpaceFactor * (fEpsilonL[0] * fMomiV[0] * cC[iTrig][nRes][jC][0][iPh]
          + fEpsilonL[1] * fMomiV[1] * cC[iTrig][nRes][jC][0][iPh]
          + fEpsilonL[2] * fMomiV[2] * cC[iTrig][nRes][jC][0][iPh]) * fLeft[nRes] / std::sqrt(2.) * GBWK1bw(s,mass,width)/sqrt(s);
          
          cRL[iTrig][nRes][jC][1][1][iPh] = kPhaseSpaceFactor * (-fEpsilonL[0] * fMomjV[0] * cC[iTrig][nRes][jC][1][iPh]
          - fEpsilonL[1] * fMomjV[1] * cC[iTrig][nRes][jC][1][iPh]
          - fEpsilonL[2] * fMomjV[2] * cC[iTrig][nRes][jC][1][iPh]) * fLeft[nRes] / std::sqrt(2.) * GBWK1bw(s,mass,width)/sqrt(s);
        }
    
    return fJVec3[nRes];
  }
}

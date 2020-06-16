#include "EvtGenModels/EvtGDecayRate.hh"

namespace Gamapola{
  GDecayRate::~GDecayRate()
  {
//     std::cout << "GCurrents destructor calling. . ." << std::endl;
  }
  double GDecayRate::GJ2(Dict2D<double>& pars2d, double* p, double* x)
  { 
    auto&& lambda = p[7];
    auto&& s = x[4];
    std::complex<double> matrixElementLX(0.,0.);
    std::complex<double> matrixElementLY(0.,0.);
    std::complex<double> matrixElementLZ(0.,0.);
    
    std::complex<double> matrixElementRX(0.,0.);
    std::complex<double> matrixElementRY(0.,0.);
    std::complex<double> matrixElementRZ(0.,0.);
    
//     std::cout << pars2d["K2_1660"]["fV_rho"] << '\n';
//     std::cout << "start?" << '\n';
    GJVec(p, x, 0, "K1_1270");
    GJVec(p, x, 1, "K1_1400");
    GLVec(p, x, 2, "K*_1410");
    GLVec(p, x, 3, "K*_1680");
    GKVec(p, x, 4, "K2_1430");
    GMVec(pars2d["K2_1660"], x, 5, "K2_1600");
    GNVec(pars2d["K2_1770"], x, 6, "K2_1770");
//     std::cout << "end?" << '\n';
    for(auto nRes = 0; nRes < 5; nRes++)
    { 
      matrixElementLX += fJVec3L[nRes][0];
      matrixElementLY += fJVec3L[nRes][1];
      matrixElementLZ += fJVec3L[nRes][2];

      matrixElementRX += fJVec3R[nRes][0];
      matrixElementRY += fJVec3R[nRes][1];
      matrixElementRZ += fJVec3R[nRes][2];
//       std::cout << fRight[nRes] * fJVec3[nRes][0] - fJVec3R[nRes][0] << "  " << nRes << '\n';
//       std::cout << nRes << '\n';
//       std::cout << "X: " << fJVec3L[nRes][0] << " " << fJVec3R[nRes][0] << " " << 0 << " " << nRes << '\n';
//       std::cout << "Y: " << fJVec3L[nRes][1] << " " << fJVec3R[nRes][1] << " " << 1 << " " << nRes << '\n';
//       std::cout << "Z: " << fJVec3L[nRes][2] << " " << fJVec3R[nRes][2] << " " << 2 << " " << nRes << '\n';
    }
    
    std::complex<double> scprL = kPhaseSpaceFactor*(matrixElementLX * fEpsilonL[0] + matrixElementLY * fEpsilonL[1] + matrixElementLZ * fEpsilonL[2]) / sqrt(s);
    std::complex<double> scprR = kPhaseSpaceFactor*(matrixElementRX * fEpsilonR[0] + matrixElementRY * fEpsilonR[1] + matrixElementRZ * fEpsilonR[2]) / sqrt(s);
//     if(scprL != scprL)
//     {
//       std::cout << matrixElementLX << "  " << matrixElementLY << "  " << matrixElementLZ << '\n';
//       std::cout << matrixElementRX << "  " << matrixElementRY << "  " << matrixElementRZ << '\n';
//     }
    //     std::cout << scprR << "  " << scprL << std::endl;
//     std::cout << kPhaseSpaceFactor << '\n';
    
    double pdf = -( ( 1 - lambda ) / 2. * std::real( scprL * std::conj( scprL ) ) + 
                         ( 1 + lambda ) / 2. * std::real( scprR * std::conj( scprR ) ) ); // could be...
    fOmega = ( -std::real( scprL * std::conj( scprL ) ) + std::real( scprR * std::conj( scprR ) ) ) / 
             ( std::real( scprR * std::conj( scprR ) ) + std::real( scprL * std::conj( scprL ) ) );
     
//     if(pdf != pdf)
//         std::cout << lambda << "  " << scprL << "   " << scprR << '\n';
//     std::cout << " Independent pdf value: " << pdf << std::endl;
    return pdf;
//     std::cout << matrixElementX << "  " << matrixElementY << "  " << matrixElementZ << std::endl;
//       return -1e-5*std::real(matrixElementX * std::conj(matrixElementX) + 
//              matrixElementY * std::conj(matrixElementY) + 
//              matrixElementZ * std::conj(matrixElementZ)) / s;
  }
}

#include "EvtGenModels/EvtGQuarkPairCreationModel.hh"

namespace Gamapola{
  GQuarkPairCreationModel::GQuarkPairCreationModel(): GInterfaceForMathFunctions(),
  fArg(0)
  {
    std::cout << "GQuarkPairCreationModel constructor calling. . ." << std::endl;
  }
  GQuarkPairCreationModel::~GQuarkPairCreationModel()
  {
//     std::cout << "GQuarkPairCreationModel destructor calling. . ." << std::endl;
  }
  void GQuarkPairCreationModel::GQPCMConstanstsCalculation(int nRes)
  {
    double ra = kRK1;
    double rb = kRKst;
    double rc = kRpi;
    double r2 = (ra*ra+rb*rb+rc*rc);
    double coeff = -4. / sqrt( sqrt( TMath::Pi() * TMath::Pi() * TMath::Pi() * TMath::Pi() * TMath::Pi() ) / 3. );
    fcoeff[nRes][0] = coeff * sqrt( ra * ra * ra * ra * ra ) * sqrt( rb * rc * rb * rc * rb * rc ) / sqrt( r2 * r2 * r2 * r2 * r2 );
    fbracePart[nRes][0] = ( r2 + ra * ra ) * ( rb * rb + rc * rc ) / ( 4 * r2 );
    fexpPart[nRes][0] = ra * ra * (rb * rb + rc * rc) / (8. * r2);
    
    fcoeff[nRes][1] = fcoeff[nRes][0];
    fbracePart[nRes][1] = fbracePart[nRes][0];
    fexpPart[nRes][1] = fexpPart[nRes][0];
    
    fcoeffAdd[nRes] = -sqrt(2.) / ( 3. * TMath::Pi() * sqrt( sqrt( TMath::Pi() ) ) ) * ra * ra * rb * rc * rb * sqrt( ra ) * sqrt( rb * rc  ) / 
    ( r2 * r2 * sqrt( r2 ) );
    
    fcoeff[nRes][2] = sqrt(2.) * fcoeffAdd[nRes] * (3 * ra * ra - 2 * ( rb * rb + rc * rc ) ) / r2;
    
    fbracePart[nRes][2] = ra * ra * ( rb * rb + rc * rc ) * ( 2 * ra * ra + rb * rb + rc * rc ) / 
    ( 4 * r2 * (3 * ra * ra - 2 * ( rb * rb + rc * rc ) ) );
    fexpPart[nRes][2] = fexpPart[nRes][0];
    fsqrt2 = sqrt(2.);
//     std::cout << r2 << "  " << fcoeffAdd[nRes] << "   " << fbracePart[nRes][2] << "   " << fexpPart[nRes][2] << std::endl;
  }
  double GQuarkPairCreationModel::GI0( int nRes, int nFnc)
  {    
    double mom2 = fMomentaA[nFnc] * fMomentaA[nFnc];
    double brace = 1 - mom2 * fbracePart[nRes][nFnc];
    double expo = exp( -mom2 * fexpPart[nRes][nFnc] );
//     std::cout << fMomentaA[nFnc] << std::endl;
    return fcoeff[nRes][nFnc] * brace * expo;
  }
  double GQuarkPairCreationModel::GI1(int nRes, int nFnc)
  {
    double mom2 = fMomentaA[nFnc] * fMomentaA[nFnc];
    double expo = exp( -mom2 * fexpPart[nRes][nFnc] );
    return -fcoeff[nRes][nFnc] * expo;
  }
  double GQuarkPairCreationModel::GMS(double* p, int nRes, int nFunc, int nSwap)
  {
    double gammaQPC = p[0];
    double f2 = p[1];
    double mom2 = fMomentaA[nFunc] * fMomentaA[nFunc];
    double sMS = sqrt( 1.5 ) * ( 2 * fI1[0][nFunc] - fI0[0][nFunc] ) / 18.;
    if(nRes == 0 && nFunc == 1)
      sExp[nSwap][nRes][nFunc] = -mom2;
    else
    {
      sExp[nSwap][nRes][nFunc] = fMomentaConst[nRes][nFunc] * fMomentaConst[nRes][nFunc] - mom2;
    }
//     sExp[nSwap][nRes][nFunc] = fMomentaConst[nRes][nFunc] * fMomentaConst[nRes][nFunc] - mom2;
//     std::cout << mom2 << "  " << fMomentaConst[nRes][nFunc] << '\n';
    double MS = gammaQPC * sMS * exp( sExp[nSwap][nRes][nFunc] * f2 );
//     std::cout << mom2 << "  " << std::endl;
    return MS;
  }
  double GQuarkPairCreationModel::GMD(double* p, int nRes, int nFunc, int nSwap)
  {
    double gammaQPC = p[0];
    double f2 = p[1];
    double mom2 = fMomentaA[nFunc] * fMomentaA[nFunc];
    double sMD = sqrt( 1.5 ) * ( fI1[0][nFunc] + fI0[0][nFunc] ) / 18.;
    double MD = gammaQPC * sMD * exp( sExp[nSwap][nRes][nFunc] * f2 );
//     std::cout << gammaQPC << " " << sMD << " " << sExp[nSwap][nRes][nFunc] << " " << f2 << "  " << nRes << " " <<'\n';
    return MD;
  }
//   Non-relativistic S and D matrix elements ************************************************************
  std::complex<double> GQuarkPairCreationModel::GMSKstarNR1270(double *p, int nSwap)
  {
    double thetaK1 = p[2];
    cM[0][0][0][0][nSwap] = fMS[0][0] * fsqrt2;
    cM[1][0][0][0][nSwap] = -fMS[0][0];
//     std::cout << fMS[0][0] * ( fsqrt2 * sin( thetaK1 ) - cos( thetaK1 ) ) << std::endl;
    return std::complex<double>(p[29], p[30]) * (cM[0][0][0][0][nSwap] * fSinThK1 + cM[1][0][0][0][nSwap] * fCosThK1);
//     return fMS[0][0] * ( fsqrt2 * sin( thetaK1 ) - cos( thetaK1 ) )/* * fExponentialM[0][0]*/;
  }
  std::complex<double> GQuarkPairCreationModel::GMSKstarNR1400(double *p, int nSwap)
  {
    double thetaK1 = p[2];
    cM[0][1][0][0][nSwap] = fMS[0][0];
    cM[1][1][0][0][nSwap] = fMS[0][0] * fsqrt2;
//     std::cout << fMS[0][0]  << std::endl;
    return std::complex<double>(p[33], p[34]) * (cM[0][1][0][0][nSwap] * fSinThK1 + cM[1][1][0][0][nSwap] * fCosThK1);
//     return fMS[0][0] * ( fsqrt2 * cos( thetaK1 ) + sin( thetaK1 ) )/* * fExponentialM[1][0]*/;
  }
  std::complex<double> GQuarkPairCreationModel::GMDKstarNR1270(double *p, int nSwap)
  {
    double thetaK1 = p[2];
    cM[0][0][0][1][nSwap] = -fMD[0][0];
    cM[1][0][0][1][nSwap] = -fsqrt2 * fMD[0][0];
    return std::complex<double>(p[31], p[32]) * (cM[0][0][0][1][nSwap] * fSinThK1 + cM[1][0][0][1][nSwap] * fCosThK1);
//     return fMD[0][0] * ( -sin( thetaK1 ) - fsqrt2 * cos( thetaK1 ) )/* * fExponentialM[0][0]*/;
  }
  std::complex<double> GQuarkPairCreationModel::GMDKstarNR1400(double *p, int nSwap)
  {
    double thetaK1 = p[2];
    cM[0][1][0][1][nSwap] = fMD[0][0] * fsqrt2;
    cM[1][1][0][1][nSwap] = -fMD[0][0];
    return std::complex<double>(p[35], p[36]) * (cM[0][1][0][1][nSwap] * fSinThK1 + cM[1][1][0][1][nSwap] * fCosThK1);
//     return fMD[0][0] * ( -cos( thetaK1 ) + fsqrt2 * sin( thetaK1 ) )/* * fExponentialM[1][0]*/;
  }
  std::complex<double> GQuarkPairCreationModel::GMSRhoNR1270(double *p, int nSwap)
  {
    double thetaK1 = p[2];
    cM[0][0][1][0][nSwap] = fMS[0][1] * fsqrt2;
    cM[1][0][1][0][nSwap] = fMS[0][1];
//     std::cout << fMS[0][1] << std::endl;
    return std::complex<double>(p[37], p[38]) * (cM[0][0][1][0][nSwap] * fSinThK1 + cM[1][0][1][0][nSwap] * fCosThK1);
//     return fMS[0][1] * ( fsqrt2 * sin( thetaK1 ) + cos( thetaK1 ) )/* * fExponentialM[0][1]*/;
  }
  std::complex<double> GQuarkPairCreationModel::GMSRhoNR1400(double *p, int nSwap)
  {
    double thetaK1 = p[2];
    cM[0][1][1][0][nSwap] = -fMS[0][1];
    cM[1][1][1][0][nSwap] = fMS[0][1] * fsqrt2;
//     std::cout << fExponentialM[1][1] << std::endl;
    return std::complex<double>(p[41], p[42]) * (cM[0][1][1][0][nSwap] * fSinThK1 + cM[1][1][1][0][nSwap] * fCosThK1);
//     return fMS[0][1] * ( fsqrt2 * cos( thetaK1 ) - sin( thetaK1 ) )/* * fExponentialM[1][1]*/;
  }
  std::complex<double> GQuarkPairCreationModel::GMDRhoNR1270(double *p, int nSwap)
  {
    double thetaK1 = p[2];
    cM[0][0][1][1][nSwap] = -fMD[0][1];
    cM[1][0][1][1][nSwap] = fMD[0][1] * fsqrt2;
    return std::complex<double>(p[39], p[40]) * (cM[0][0][1][1][nSwap] * fSinThK1 + cM[1][0][1][1][nSwap] * fCosThK1);
//     return fMD[0][1] * ( -sin( thetaK1 ) + fsqrt2 * cos( thetaK1 ) ) /** fExponentialM[0][1]*/;
  }
  std::complex<double> GQuarkPairCreationModel::GMDRhoNR1400(double *p, int nSwap)
  {
    double thetaK1 = p[2];
    cM[0][1][1][1][nSwap] = -fMD[0][1] * fsqrt2;
    cM[1][1][1][1][nSwap] = -fMD[0][1];
    return std::complex<double>(p[43], p[44]) * (cM[0][1][1][1][nSwap] * fSinThK1 + cM[1][1][1][1][nSwap] * fCosThK1);
//     return fMD[0][1] * ( -cos( thetaK1 ) - fsqrt2 * sin( thetaK1 ) )/* * fExponentialM[1][1]*/;
  }
//   ***************************Relativistic S and D matrix elements******************************************************************
  std::complex<double> GQuarkPairCreationModel::GMSKstar1270(double *p, int nSwap)
  {
    cM[0][0][0][0][nSwap] *= fEnergyCoeff[0][0];
    cM[1][0][0][0][nSwap] *= fEnergyCoeff[0][0];
//     std::cout << fEnergyCoeff[0][0] * fMK0starNR[0][0] << std::endl;
    return fEnergyCoeff[0][0] * fMK0starNR[0][0];
  }
  std::complex<double> GQuarkPairCreationModel::GMSKstar1400(double *p, int nSwap)
  {
    cM[0][1][0][0][nSwap] *= fEnergyCoeff[0][0];
    cM[1][1][0][0][nSwap] *= fEnergyCoeff[0][0];
//     std::cout << sMS[nSwap][1][0] << std::endl;
    return fEnergyCoeff[0][0] * fMK0starNR[1][0];
  }
  std::complex<double> GQuarkPairCreationModel::GMDKstar1270(double *p, int nSwap)
  {
    cM[0][0][0][1][nSwap] *= fEnergyCoeff[0][0];
    cM[1][0][0][1][nSwap] *= fEnergyCoeff[0][0];
    return fEnergyCoeff[0][0] * fMK0starNR[0][1];
  }
  std::complex<double> GQuarkPairCreationModel::GMDKstar1400(double *p, int nSwap)
  {
    cM[0][1][0][1][nSwap] *= fEnergyCoeff[0][0];
    cM[1][1][0][1][nSwap] *= fEnergyCoeff[0][0];
    return fEnergyCoeff[0][0] * fMK0starNR[1][1];
  }
  std::complex<double> GQuarkPairCreationModel::GMSRho1270(double *p, int nSwap)
  {
    cM[0][0][1][0][nSwap] *= fEnergyCoeff[0][1];
    cM[1][0][1][0][nSwap] *= fEnergyCoeff[0][1];
//     std::cout << "GMSRho1270  " << cM[0][0][1][0][nSwap] << "  " << cM[1][0][1][0][nSwap] << "  " << fEnergyCoeff[0][1] * fMRho0NR[0][0] << '\n';
    return fEnergyCoeff[0][1] * fMRho0NR[0][0];
  }
  std::complex<double> GQuarkPairCreationModel::GMSRho1400(double *p, int nSwap)
  {
    cM[0][1][1][0][nSwap] *= fEnergyCoeff[0][1];
    cM[1][1][1][0][nSwap] *= fEnergyCoeff[0][1];
//     std::cout << "GMSRho1400  " << cM[0][1][1][0][nSwap] << "  " << cM[1][1][1][0][nSwap] << "  " << fEnergyCoeff[0][1] * fMRho0NR[1][0] << '\n';
    return fEnergyCoeff[0][1] * fMRho0NR[1][0];
  }
  std::complex<double> GQuarkPairCreationModel::GMDRho1270(double *p, int nSwap)
  {
    cM[0][0][1][1][nSwap] *= fEnergyCoeff[0][1]; // 0 when comparing with Emi
    cM[1][0][1][1][nSwap] *= fEnergyCoeff[0][1]; // 0 when comparing with Emi
    return fEnergyCoeff[0][1] * fMRho0NR[0][1];  // 0 when comparing with Emi
  }
  std::complex<double> GQuarkPairCreationModel::GMDRho1400(double *p, int nSwap)
  {
    cM[0][1][1][1][nSwap] *= fEnergyCoeff[0][1];  // 0 when comparing with Emi
    cM[1][1][1][1][nSwap] *= fEnergyCoeff[0][1];  // 0 when comparing with Emi
    return fEnergyCoeff[0][1] * fMRho0NR[1][1];  // 0 when comparing with Emi
  }
  
    double GQuarkPairCreationModel::GMPBNR(double *p, int nSwap)
  {
    double gammaQPC = p[0];
    double f2 = p[1];
    double mom2 = fMomentaA[2] * fMomentaA[2];
    double expo = exp( -mom2 * fexpPart[0][2] ); 
    double expo2 = exp( -mom2 * f2);
    fMPBNR[nSwap] = gammaQPC * fcoeff[0][2] * expo * expo2 * 
    (1 - fbracePart[0][2] * mom2) * fMomentaA[2];
    return fMPBNR[nSwap];
  }
  double GQuarkPairCreationModel::GMPANR(double *p, int nSwap)
  {
    double gammaQPC = p[0];
    double f2 = p[1];
    double mom2 = fMomentaA[2] * fMomentaA[2];
    double expo = exp( -mom2 * fexpPart[0][2] ); 
    double expo2 = exp( -mom2 * f2);
    fMPANR[nSwap] = gammaQPC * fcoeff[0][2] * expo * expo2 * fMomentaA[2];
    return fMPANR[nSwap];
  }
}

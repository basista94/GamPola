#ifndef GINTERFACEFORMATHFUNCTIONS_HH
#define GINTERFACEFORMATHFUNCTIONS_HH

#include <iostream>
#include <cmath>
#include <TMath.h>
#include <complex>
#include <cstring>
#include "Math/WrappedMultiTF1.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "EvtGenModels/EvtGPhysicsPDG.hh"
// #include "EvtGenModels/EvtGGrammSmidt.hh"
#include "EvtGenModels/EvtGVector4D.hh"
#include "EvtGenModels/EvtGVector3D.hh"

namespace Gamapola
{
  class GInterfaceForMathFunctions{
  public:
    GInterfaceForMathFunctions();
    virtual ~GInterfaceForMathFunctions();
    virtual void GSetResonances(const int& charge) = 0;
    virtual double GOmegaPDF1(double* x, double* p) = 0;
//*************************************QPCM***************************
    virtual double GProcessingComputationOfPDF(double* x, double* p) = 0;
//     virtual double GMomentaA(double s, double sij, double massMes) = 0;
//     virtual double GPiPjVec4(double s, double sij, double sjk) = 0;
//     virtual void GPiPjVec3() = 0;
    virtual std::complex<double> GScPrd(const double* pi, const std::complex<double>* pj) = 0;
//     virtual double GEnergyA(double s, double sij, double massMes) = 0;
    virtual double GSij(double s, double sij, double sik) = 0;
//     virtual double GPiPj() = 0;
    virtual void GCalculateKinematics(double* p, const int& charge) = 0;
    virtual double GI0(int nRes, int nFnc) = 0;
    virtual double GI1(int nRes, int nFnc) = 0;
    virtual double GMS(double* p, int nRes, int nFnc, int nSwap) = 0;
    virtual double GMD(double* p, int nRes, int nFnc, int nSwap) = 0;
//     virtual double GExponentialM(double *p, int nRes, int nFnc) = 0;

//     Non-relativistic S and D waves
    virtual std::complex<double> GMSKstarNR1270(double *p, int nSwap) = 0;
    virtual std::complex<double> GMDKstarNR1270(double *p, int nSwap) = 0;
    virtual std::complex<double> GMSRhoNR1270(double *p, int nSwap) = 0;
    virtual std::complex<double> GMDRhoNR1270(double *p, int nSwap) = 0;
    virtual std::complex<double> GMSKstarNR1400(double *p, int nSwap) = 0;
    virtual std::complex<double> GMDKstarNR1400(double *p, int nSwap) = 0;
    virtual std::complex<double> GMSRhoNR1400(double *p, int nSwap) = 0;
    virtual std::complex<double> GMDRhoNR1400(double *p, int nSwap) = 0;
//     Relativistic S and D waves
    virtual std::complex<double> GMSKstar1270(double *p, int nSwap) = 0;
    virtual std::complex<double> GMDKstar1270(double *p, int nSwap) = 0;
    virtual std::complex<double> GMSRho1270(double *p, int nSwap) = 0;
    virtual std::complex<double> GMDRho1270(double *p, int nSwap) = 0;
    virtual std::complex<double> GMSKstar1400(double *p, int nSwap) = 0;
    virtual std::complex<double> GMDKstar1400(double *p, int nSwap) = 0;
    virtual std::complex<double> GMSRho1400(double *p, int nSwap) = 0;
    virtual std::complex<double> GMDRho1400(double *p, int nSwap) = 0;
    
    virtual double GMPBNR(double *p, int nSwap) = 0;
    virtual double GMPANR(double *p, int nSwap) = 0;
    
    virtual std::complex<double> GAKstar(double* x,double* p, int nRes, int nSwap) = 0;
    virtual std::complex<double> GBKstar(double* x,double* p, int nRes, int nSwap) = 0;
    virtual std::complex<double> GARho(double* x,double* p, int nRes, int nSwap) = 0;
    virtual std::complex<double> GBRho(double* x,double* p, int nRes, int nSwap) = 0;
    
    virtual void GCouplingConstantsCalculation() = 0;
//     Breit-Wigner
    virtual std::complex<double> GBWKstr(double sij, int nSwap) = 0;
    virtual std::complex<double> GBWRho(double sij) = 0;
    virtual std::complex<double> GBWKappa(double sij) = 0;
    virtual std::complex<double> GBWK1(double s, int nRes) = 0;
    virtual std::complex<double> GBWK1bw(const double& s, const double& mass, const double& width) = 0;
    virtual double GWidth(const double& s, const double& sbc, int nRes) = 0;
    
    virtual void GInitializeResVector(const int& charge, const std::string& resName) = 0;
    virtual void GInitializeResPseudoVector(const int& charge, const std::string& resName) = 0;
    virtual void GInitializeResTensor(const int& charge, const std::string& resName) = 0;
    virtual void GInitializeResPseudoTensor(const int& charge, const std::string& resName) = 0;
    virtual void GInitializeResPseudoTensor2Tensor(const int& charge, const std::string& resName) = 0;
    
    virtual std::complex<double> GC1(double* p, double s, double sij, double sjk, int nRes) = 0;
    virtual std::complex<double> GC2(double* p, double s, double sij, double sjk, int nRes) = 0;
    virtual std::complex<double> GC3(double* p, double s, double sij, double sjk, int nRes) = 0;
    virtual std::complex<double> GC4(double* p, double s, double sij, double sjk, int nRes) = 0;
    virtual std::complex<double> GC5(double* p, double s, double sij, double sjk, int nRes) = 0;
    
    
    virtual double GJ2(Dict2D<double>& pars2d, double* p, double* x) = 0;
    virtual std::complex<double>* GJVec(double* p, double* x, int nRes, const std::string& resName) = 0;
    virtual std::complex<double>* GLVec(double* p, double* x, int nRes, const std::string& resName) = 0;
    virtual std::complex<double>* GKVec(double* p, double* x, int nRes, const std::string& resName) = 0;
    virtual std::vector<std::complex<double>> GMVec( Dict1D<double>& pars, double* x, int nRes, const std::string& resName) = 0;
    virtual void GNVec( Dict1D<double>& pars, double* x, int nRes, const std::string& resName) = 0;
    virtual std::vector<std::complex<double> > GGetKinematicalCoefficients() = 0;
  protected:
    static const int fMaxNRes = 7;
    static const int fMaxNDecays = 3;
    static const int fMaxNCharges = 3;
    static const int fMaxNWaves = 3;
    int fArg;
    double fMasses[4];
    const std::complex<double>* fNormalizationIntegrals;
    int fNIntegrals;
//     std::vector<double (Gamapola::GMathFunctions::*) (double* x, double* p)> fFunctionPointers; 
//     For I0, I1:
    int fIndex[fMaxNRes];
    std::vector<int> fIndexFor2Decays;
    std::vector<std::vector<std::vector<double> > > fConstParsForI01; 
    std::vector<std::vector<std::vector<double> > > fConstParsForMomA;
    std::vector<std::vector<double> >  fMomentaConst;
    double fMomentaA[fMaxNDecays];
    double fEnergyA[fMaxNDecays];
    double fEnergyCoeff[fMaxNRes][fMaxNDecays];
    std::vector<std::complex<double> > fKinematicalCoefficients;
    double fResMasses[fMaxNRes];
    double fResWidths[fMaxNRes];
    double fResBrs[fMaxNRes][fMaxNDecays];
    double fI0[fMaxNRes][fMaxNDecays];
    double fI1[fMaxNRes][fMaxNDecays];
    double fMS[fMaxNRes][fMaxNDecays];
    double fMD[fMaxNRes][fMaxNDecays];
    std::vector<double> fMPBNR;
    std::vector<double> fMPANR;
    
    std::complex<double> fMK0starNR[fMaxNRes][fMaxNDecays];
    std::complex<double> fMRho0NR[fMaxNRes][fMaxNDecays];
    std::complex<double> fMK0star[fMaxNRes][fMaxNDecays];
    std::complex<double> fMRho0[fMaxNRes][fMaxNDecays];
    
    std::complex<double> fA[fMaxNRes][fMaxNDecays];
    std::complex<double> fB[fMaxNRes][fMaxNDecays];
    std::complex<double> fA2[fMaxNRes][fMaxNDecays];
    std::complex<double> fB2[fMaxNRes][fMaxNDecays];
    
    std::complex<double> fc4Rho;
    std::complex<double> fc4Kstr;
    std::complex<double> fc5Kstr;
    
    std::complex<double> fc1[2][fMaxNCharges];
    std::complex<double> fc2[2][fMaxNCharges];
    std::complex<double> fc3[2][fMaxNCharges];
    std::complex<double> fc4[2][fMaxNCharges];
    std::complex<double> fc5[2][fMaxNCharges];
    
    std::complex<double> fJVec3[fMaxNRes][3];
    std::complex<double> fJVec3L[fMaxNRes][3];
    std::complex<double> fJVec3R[fMaxNRes][3];
    
    std::complex<double> (Gamapola::GInterfaceForMathFunctions::*fFunctionPointersKstrNR[fMaxNRes][fMaxNWaves])(double* p, int nSwap); 
    std::complex<double> (Gamapola::GInterfaceForMathFunctions::*fFunctionPointersKstr[fMaxNRes][fMaxNWaves]) (double* p, int nSwap);
    std::complex<double> (Gamapola::GInterfaceForMathFunctions::*fFunctionPointersRhoNR[fMaxNRes][fMaxNWaves]) (double* p, int nSwap);
    std::complex<double> (Gamapola::GInterfaceForMathFunctions::*fFunctionPointersRho[fMaxNRes][fMaxNWaves]) (double* p, int nSwap);
    
//     Coupling GCouplingConstantsCalculation
    double fgKstarKpi;
    double fgRhoPiPi;
    double fgK1KappaPi[2];
    const char** fResNames; 
    double fCGKstr[fMaxNRes];
    double fCGRho[fMaxNRes];
    double fDelta[fMaxNRes];
    static constexpr double fK = 1.;
    static constexpr double fRho = 1.;
    std::complex<double> fDelta2[fMaxNRes];
    std::complex<double> fDeltaKappa[fMaxNRes];
    std::complex<double> fCCouplings[fMaxNRes][fMaxNDecays];
    double fFF[fMaxNRes][fMaxNDecays];
    double fK1[2];
//     Model parameters
    int fCharge[fMaxNRes];
    double fOmega;
//     Some kinematical variables
    double fEnergy1;
    double fEnergy2;
    double fEnergy3;
    double fEnergy4;
    
    double fPiPjVec4;
    double fPiPjVec[3];
    double fSij;
    
    double fMomi;
    double fMomj;
    double fMomk;
    double fMoml;
    
    double fMomij;
    double fPhii;
    double fPhij;
    double fPhik;
    
    double fMomiV[3];
    double fMomjV[3];
    double fMomkV[3];
    double fMomlV[3];
    
    GVector4D fqK;
    GVector4D fqPi1;
    GVector4D fqPi2;
    
    std::complex<double> fE0PiPj;
    std::complex<double> fE0Pi;
    std::complex<double> fE0Pj;
    std::complex<double> fK2KinC4combo[3];
    std::complex<double> fK2KinC5combo[3];
//     Polarization vectors
    std::complex<double> fEpsilonL[3];
    std::complex<double> fEpsilonR[3];
    std::complex<double> fEpsilon0[3];
    double fSinThK1;
    double fCosThK1;
//     Some intermediate shit for acceleration of symbolic calculation
    double sExp[2][fMaxNRes][fMaxNDecays];
    double cM[2][2][2][2][2]; // {sin,cos}; {1270,1400}; {K*,Ro}; {S,D}; {0,1}(nSwap)
    std::complex<double> cF[2][2][2][2][2][2]; // {sin,cos}; {1270,1400}; {K*,Ro}; {f,h}; {0,1}(nSwap); {iPhS,iPhD}
    std::complex<double> cC[2][fMaxNRes][3][2][2]; // {sin,cos}; {1270,1400}; {K*,Ro, Kappa}; {C1,C2}; {iPhS,iPhD}
    std::complex<double> cRL[2][fMaxNRes][3][2][2][2]; // {sin,cos}; {1270,1400}; {K*,Ro, Kappa}; {C1,C2}; {R,L}; {iPhS,iPhD}
    
    double fcoeff[fMaxNRes][fMaxNDecays];
    double fbracePart[fMaxNRes][fMaxNDecays];
    double fexpPart[fMaxNRes][fMaxNDecays];
    double fcoeffAdd[2];
    double fbracePartAdd[2];
    double fexpPartAdd[2];
    double fLeft[fMaxNRes];
    double fRight[fMaxNRes];
    double fHadronicFF[fMaxNRes][fMaxNDecays];
    int fNRes;
//     GGrammSmidt* fGS;
    GVector4D fq_gamma;
    GVector4D fq_K;
    GVector4D fq_pi1;
    GVector4D fq_pi2;
    
    double kMKaon;// = 0.49367700;
    double kMPion1;// = 0.13957018;
    double kMPion2;// = 0.13957018;//0.134976;
  };
}
#endif

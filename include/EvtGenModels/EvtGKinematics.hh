#ifndef GKINEMATICS_HH
#define GKINEMATICS_HH

#include "EvtGenModels/EvtGInterfaceForMathFunctions.hh"
#include "EvtGenModels/EvtGSimpleRandomEngine.hh"
#include "EvtGenModels/EvtGVector4D.hh"
#include "EvtGenModels/EvtGVector3D.hh"

namespace Gamapola
{
  class GKinematics: virtual public GInterfaceForMathFunctions{
  public:
    GKinematics();
    virtual ~GKinematics();
    static double GSjkMin(const double& sij, const double& mi, const double& mj, const double& mk, const double& M);
    static double GSjkMax(const double& sij, const double& mi, const double& mj, const double& mk, const double& M);
    static double GMomentaA(double s, double sij, double massMes);
    std::complex<double> GScPrd(const double* pi, const std::complex<double>* pj);
    static double GEnergyA(double s, double sij, double massMes);
    double GSij(double s, double sij, double sik);
    void GCalculateKinematics(double* p, const int& charge);
//     double GExponentialM(double *p, int nRes, int nFnc);
//     double GCosThetaij(const double& cosThstar, const double& sij, const double& mi, 
//                        const double& mj, const double& mk, const double& M);
    void GKinematicsConstanstsCalculation(const double& mRes, const int& nRes);
    void GPhaseSpaceTo4Vectors();
    static std::vector<double> G4VectorsToPhaseSpace(const GVector4D& qK, const GVector4D& qpi1, const GVector4D& qpi2, const GVector4D& qGamma, const int& charge);
    
    std::vector<double> GGetPhaseSpace();
    const GVector4D& GGet4VecGamma() const { return fq_gamma; }
    const GVector4D& GGet4VecK() const { return fq_K; }
    const GVector4D& GGet4VecPi1() const { return fq_pi1; }
    const GVector4D& GGet4VecPi2() const { return fq_pi2; }
  private:
    double GPi1Pi2Vec4(double s, double sij, double sjk);
    void GPi1Pi2Vec3();
    double GPiPj();
    int fArg;
    GSimpleRandomEngine* fRandomEngine;
    double f_Sij;
    double fCosTheta;
  };
}
#endif

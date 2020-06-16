#ifndef GQUARKPAIRCREATIONMODEL_HH
#define GQUARKPAIRCREATIONMODEL_HH

#include "EvtGenModels/EvtGInterfaceForMathFunctions.hh"

namespace Gamapola
{
  class GQuarkPairCreationModel: virtual public GInterfaceForMathFunctions{
  public:
    GQuarkPairCreationModel();
    virtual ~GQuarkPairCreationModel();
    double GI0(int nRes, int nFnc);
    double GI1(int nRes, int nFnc);
    double GMS(double* p, int nRes, int nFnc, int nSwap);
    double GMD(double* p, int nRes, int nFnc, int nSwap);
    //     Non-relativistic S and D waves
    std::complex<double> GMSKstarNR1270(double *p, int nSwap);
    std::complex<double> GMDKstarNR1270(double *p, int nSwap);
    std::complex<double> GMSRhoNR1270(double *p, int nSwap);
    std::complex<double> GMDRhoNR1270(double *p, int nSwap);
    std::complex<double> GMSKstarNR1400(double *p, int nSwap);
    std::complex<double> GMDKstarNR1400(double *p, int nSwap);
    std::complex<double> GMSRhoNR1400(double *p, int nSwap);
    std::complex<double> GMDRhoNR1400(double *p, int nSwap);
    double GMPBNR(double *p, int nSwap);
    double GMPANR(double *p, int nSwap);
//     Relativistic S and D waves
    std::complex<double> GMSKstar1270(double *p, int nSwap);
    std::complex<double> GMDKstar1270(double *p, int nSwap);
    std::complex<double> GMSRho1270(double *p, int nSwap);
    std::complex<double> GMDRho1270(double *p, int nSwap);
    std::complex<double> GMSKstar1400(double *p, int nSwap);
    std::complex<double> GMDKstar1400(double *p, int nSwap);
    std::complex<double> GMSRho1400(double *p, int nSwap);
    std::complex<double> GMDRho1400(double *p, int nSwap);
    void GQPCMConstanstsCalculation(int nRes);
  private:
    int fArg;
    double fsqrt2;
    double fsqrt15;
  };
}
#endif

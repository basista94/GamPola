#ifndef GFORMFACTORS_HH
#define GFORMFACTORS_HH

#include "EvtGenModels/EvtGInterfaceForMathFunctions.hh"

namespace Gamapola
{
  class GFormFactors: virtual public GInterfaceForMathFunctions{
  public:
    GFormFactors();
    virtual ~GFormFactors();
    std::complex<double> GAKstar(double* x,double* p, int nRes, int nSwap);
    std::complex<double> GBKstar(double* x,double* p, int nRes, int nSwap);
    std::complex<double> GARho(double* x,double* p, int nRes, int nSwap);
    std::complex<double> GBRho(double* x,double* p, int nRes, int nSwap);
  private:
    int fArg;
    double fsqrt05;
  };
}
#endif

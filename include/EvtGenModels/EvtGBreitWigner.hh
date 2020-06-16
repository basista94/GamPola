#ifndef GBREITWIGNER_HH
#define GBREITWIGNER_HH

#include "EvtGenModels/EvtGInterfaceForMathFunctions.hh"

namespace Gamapola
{
  class GBreitWigner: virtual public GInterfaceForMathFunctions{
  public:
    GBreitWigner();
    virtual ~GBreitWigner();
    std::complex<double> GBWKstr(double sij, int nSwap);
    std::complex<double> GBWRho(double sij);
    std::complex<double> GBWK1(double s, int nRes);
    std::complex<double> GBWKappa(double sij);
    std::complex<double> GBWK1bw(const double& s, const double& mass, const double& width);
    double GWidth(const double& s, const double& sbc, int nRes);
  private:
    int fArg;
    double GBlattWeisskopfFactor(const double& q, const double& q0) ;
    double GWidthBarrierKstr(const double& sbc);
  };
}
#endif

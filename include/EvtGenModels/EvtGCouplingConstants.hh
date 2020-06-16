#ifndef GCOUPLINGCONSTANTS_HH
#define GCOUPLINGCONSTANTS_HH

#include "EvtGenModels/EvtGInterfaceForMathFunctions.hh"

namespace Gamapola
{
  class GCouplingConstants: virtual public GInterfaceForMathFunctions{
  public:
    GCouplingConstants();
    virtual ~GCouplingConstants();
    void GCouplingConstantsCalculation();
  private:
    int fArg;
    double fs;
  };
}
#endif

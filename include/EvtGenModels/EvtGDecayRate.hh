#ifndef GDECAYRATE_HH
#define GDECAYRATE_HH
#include "EvtGenModels/EvtGInterfaceForMathFunctions.hh"

namespace Gamapola
{
    class GDecayRate: virtual public GInterfaceForMathFunctions
    {
    public:
        virtual ~GDecayRate();
        double GJ2(Dict2D<double>& pars2d, double* p, double* x);
    };
}
#endif

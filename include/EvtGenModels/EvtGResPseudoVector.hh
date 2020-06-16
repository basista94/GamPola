#ifndef GPSEUDOVECTOR_HH
#define GPSEUDOVECTOR_HH
#include "EvtGenModels/EvtGInterfaceForMathFunctions.hh"

namespace Gamapola
{
    class GResPseudoVector: virtual public GInterfaceForMathFunctions
    {
    public:
        GResPseudoVector();
        virtual ~GResPseudoVector();
        
        void GInitializeResPseudoVector(const int& charge, const std::string& resName);
        std::complex<double>* GLVec(double* p, double* x, int nRes, const std::string& resName);
        std::complex<double> GC3(double* p, double s, double sij, double sjk, int nRes);
    private:
        std::complex<double> GC3Kst(double* p, double s, double sij, double sjk, int nRes);
        std::complex<double> GC3Rho(double* p, double s, double sij, double sjk, int nRes);
        
        double fCG_Kstr;
        double fCG_Rho;
        double fCoeffDelta;
        Dict2D<double> fRes;
    };
}
#endif

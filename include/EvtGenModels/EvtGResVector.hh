#ifndef GRESVECTOR_HH
#define GRESVECTOR_HH

#include "EvtGenModels/EvtGInterfaceForMathFunctions.hh"

namespace Gamapola
{
    class GResVector: virtual public GInterfaceForMathFunctions
    {
    public:
        GResVector();
        virtual ~GResVector();
        void GInitializeResVector(const int& charge, const std::string& resName);
        std::complex<double>* GJVec(double* p, double* x, int nRes, const std::string& resName);
    private:
        std::complex<double> GC1(double* p, double s, double sij, double sjk, int nRes);
        std::complex<double> GC2(double* p, double s, double sij, double sjk, int nRes);
        std::complex<double> GC1Kst(double* p, double s, double sij, double sjk, int nRes);
        std::complex<double> GC2Kst(double* p, double s, double sij, double sjk, int nRes);
        std::complex<double> GC1Rho(double* p, double s, double sij, double sjk, int nRes);
        std::complex<double> GC2Rho(double* p, double s, double sij, double sjk, int nRes);
        std::complex<double> GC1Kappa(double* p, double s, double sij, double sjk, int nRes);
        std::complex<double> GC2Kappa(double* p, double s, double sij, double sjk, int nRes);
        
        double fCG_Kstr;
        double fCG_Rho;
        double fCoeffDelta;
        Dict2D<double> fRes;
    };
}

#endif

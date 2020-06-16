#ifndef GTENSORRES_HH
#define GTENSORRES_HH

#include "EvtGenModels/EvtGInterfaceForMathFunctions.hh"

namespace Gamapola
{
    class GResTensor: virtual public GInterfaceForMathFunctions
    {
    public:
        GResTensor();
        virtual ~GResTensor();
        
        void GInitializeResTensor(const int& charge, const std::string& resName);
        std::complex<double>* GKVec(double* p, double* x, int nRes, const std::string& resName);
    private:
        std::complex<double> GC4(double* p, double s, double sij, double sjk, int nRes);
        std::complex<double> GC5(double* p, double s, double sij, double sjk, int nRes);
        std::complex<double> GC4Kst(double* p, double s, double sij, double sjk, int nRes);
        std::complex<double> GC5Kst(double* p, double s, double sij, double sjk, int nRes);
        std::complex<double> GC4Rho(double* p, double s, double sij, double sjk, int nRes);
        
        double fCG_Kstr;
        double fCG_Rho;
        double fCoeffDelta;
        Dict2D<double> fRes;
    };
}

#endif

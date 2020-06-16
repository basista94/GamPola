#ifndef GPSEUDOTENSORRES_HH
#define GPSEUDOTENSORRES_HH

#include "EvtGenModels/EvtGInterfaceForMathFunctions.hh"
#include <unordered_map>
namespace Gamapola
{
    class GResPseudoTensor: virtual public GInterfaceForMathFunctions
    {
    public:
        GResPseudoTensor();
        virtual ~GResPseudoTensor();
        
        void GInitializeResPseudoTensor(const int& charge, const std::string& resName);
        std::vector<std::complex<double>> GMVec( Dict1D<double>& pars, double* x, int nRes, const std::string& resName);
    private:
        
        void GCaching(Dict1D<double>& pars, double s, double sij, double sjk);
        std::complex<double> GC6( Dict1D<double>& pars, double s, double sij, double sjk);
        std::complex<double> GC7( Dict1D<double>& pars, double s, double sij, double sjk);
        
        std::complex<double> GKKstr(double s, double sij);
        std::complex<double> GKRho(double s, double sij);
        
        std::complex<double> GCKstr1(double s, double sij, double fV, double hV);
        std::complex<double> GCKstr2(double s, double sij, double fV, double hV);
        std::complex<double> GCRho(double s, double sij, double fV, double hV);
        
        double fCG_Kstr;
        double fCG_Rho;
        double fCoeffDelta;
        double fArg;
        Dict2D<double> fRes;
    };
}

#endif

#ifndef GPSEUDOTENSORRES2TENSOR_HH
#define GPSEUDOTENSORRES2TENSOR_HH

#include "EvtGenModels/EvtGInterfaceForMathFunctions.hh"
#include <unordered_map>
namespace Gamapola
{
    class GResPseudoTensor2Tensor: virtual public GInterfaceForMathFunctions
    {
    public:
        GResPseudoTensor2Tensor();
        virtual ~GResPseudoTensor2Tensor();
        
        void GInitializeResPseudoTensor2Tensor(const int& charge, const std::string& resName);
        void GNVec( Dict1D<double>& pars, double* x, int nRes, const std::string& resName);
    private:
        
        std::vector<std::complex<double>> GN1Vec( Dict1D<double>& pars, double* x, int nRes, const std::string& resName);
        std::vector<std::complex<double>> GN2Vec( Dict1D<double>& pars, double* x, int nRes, const std::string& resName);
        
        double fCG_Kstr;
        double fCG_Rho;
        double fCoeffDelta;
        double fArg;
        Dict2D<double> fRes;
    };
}

#endif

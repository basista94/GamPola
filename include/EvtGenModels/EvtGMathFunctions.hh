#ifndef GMATHFUNCTIONS_HH
#define GMATHFUNCTIONS_HH

#include "EvtGenModels/EvtGInterfaceForMathFunctions.hh"
#include "EvtGenModels/EvtGKinematics.hh"
#include "EvtGenModels/EvtGBreitWigner.hh"
#include "EvtGenModels/EvtGCouplingConstants.hh"
#include "EvtGenModels/EvtGQuarkPairCreationModel.hh"
#include "EvtGenModels/EvtGFormFactors.hh"
// #include "EvtGenModels/EvtGCurrents.hh"
#include "EvtGenModels/EvtGResPseudoVector.hh"
#include "EvtGenModels/EvtGResTensor.hh"
#include "EvtGenModels/EvtGResVector.hh"
#include "EvtGenModels/EvtGResPseudoTensor.hh"
#include "EvtGenModels/EvtGResPseudoTensor2Tensor.hh"
#include "EvtGenModels/EvtGDecayRate.hh"

namespace Gamapola
{
  class GMathFunctions: public GKinematics, 
                        public GBreitWigner, 
                        public GCouplingConstants,
                        public GQuarkPairCreationModel,
                        public GFormFactors,
                        public GResPseudoVector,
                        public GResTensor,
                        public GResVector,
                        public GResPseudoTensor,
                        public GResPseudoTensor2Tensor,
                        public GDecayRate {
  public:
    GMathFunctions();
    virtual ~GMathFunctions();
    double GOmegaPDF1(double* x, double* p);
    double GProcessingComputationOfPDF(double* x, double* p);
    void GSetResonances(const int& charge);
    const std::vector<double (Gamapola::GMathFunctions::*) (double* x, double* p)>& GGetFunctionPointers() const
    {
      return fFunctionPointers;
    }
    std::vector<std::complex<double> > GGetKinematicalCoefficients();
  private:
    std::vector<double (Gamapola::GMathFunctions::*) (double* x, double* p)> fFunctionPointers; 
    void GSetupParameters(double* p);
    void GSetupDefaultParameters(double* p);
    Dict2D<double> fPars2D;
  };
}
#endif

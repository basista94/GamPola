#ifndef GINTERFACEFORMINIMIZATION_H
#define GINTERFACEFORMINIMIZATION_H

#include <iostream>
#include <cstring>
#include <sstream>
#include <fstream>
#include <string>
#include <cstdlib>
// #include <ionamip>
#include "EvtGenModels/EvtGMathFunctions.hh"
#include "GSymbolicExpressionsAnalyzer.h"
#include "Math/WrappedMultiTF1.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "TMinuitMinimizer.h"

namespace Gamapola
{
  class GInterfaceForMinimization{
  public:
    GInterfaceForMinimization();
    virtual ~GInterfaceForMinimization();
    virtual void GSetFunction(double (Gamapola::GMathFunctions::*memFunction)(double* x, double* p));
    virtual void GSetFunctionParameters(const int& nPars, double* pars);
    virtual void GSetFunctionVariables(const int& nVars, double* initialVarValues, 
                               double* lowVarsLimits, double* upVarsLimits, 
                               double* varsSteps, const std::vector<std::string>& kindsOfVars, const std::vector<std::string>& namesOfVars, bool randomize = false);
    virtual void GSetNameOfMinimizer(const char* nameOfMinimizer);
    virtual void GSetNRuns(const int& nRuns);
    virtual void GGenerateEvents(const double& low = std::numeric_limits<float>::min(), const double& up = std::numeric_limits<float>::max()) = 0;
//     void GMinimizeWithMinuit();
    void GMinimizeWithMonteCarlo();
    
    void GSetPDFExpression(const ex& pdfExpression);
//     int GInitializeSymbolicComputation(const int& nPars, const char** resNames, ex* modelPars, double* valsM, int* charges);
//     void GWriteNormalizationIntegrals(const char* filename, const char* reachivationMode, const int& nPars, const char** resNames, ex* modelPars, double* valsM, int* charges);
//     void GReadNormalizationIntegrals(const char* filename, const int& nPars, const char** resNames, ex* modelPars, double* valsM, int* charges);
    void GSetEventsNumber(const int& nEvents);
    void GSetDecayMode(const int& charge);
    const std::vector<std::vector<double> >& GGetErrorMatrix() const
    {
      return fErrorMatrix;
    }
    
    const double& GGetMinFuncValue() const
    {
      return fMinFunctionValue;
    }
    const double& GGetMaxFuncValue() const
    {
      return fMaxFunctionValue;
    }
    double* GGetFitParameters()
    {
      return fFitParameters;
    }
    TF1* GGetFuncPointer() const
    {
      return fFuncPointer;
    }
    const double* GGetCovarianceMatrix()
    {
      return fCovarianceMatrix;
    }
    const int& GGetStatus() const
    {
      return fStatus;
    }
    const double& GGetEdm()
    {
      return fEdm;
    }
    void GPrintErrorMatrix();
    double GPartialDerivative(double (Gamapola::GMathFunctions::*memFunction)(double* x, double* p), int nVar);
  protected:
    int fArg;
    double (Gamapola::GMathFunctions::*fMemFunction)(double* x, double* p);
    std::shared_ptr<GMathFunctions> fMf;
    GSymbolicExpressionsAnalyzer* fSea;
    GFUNCP_CUBA fPointersOnFunctionsOfKinPars;
    GFUNCP_CUBA fPointersOnFunctionsOfModelPars;
    const char* fMemFunctionName;
    const char* fMemFunctionClassName;
    const char* fNameOfMinimizer;
    const char* fReachivationMode;
    
    int fNParameters;
    std::vector<double> fParsValues;
    
    int fNVariables;
    std::vector<double> fInitialVarValues;
    std::vector<double> fLowVarsLimits;
    std::vector<double> fUpVarsLimits; 
    std::vector<double> fVarsSteps;
    std::vector<std::string> fKindsOfVars;
    std::vector<std::string> fNamesOfVars;
    
    double fMinFunctionValue;
    double fMaxFunctionValue;
    double* fFitParameters;
    double* fCovarianceMatrix;
    TF1* fFuncPointer;
    std::vector<std::vector<double> > fErrorMatrix;
    std::vector<std::string> fparNames;
    int fStatus;
    //     some generated values
    double fIntegralOfPDF;
    std::complex<double>* fNormalizationIntegrals;
    int fNIntegrals;
    bool fReadWrite;
    const char* fFileWithIntegrals;
    int fNRuns;
    double fEdm;
    bool fRandomize;
    int fNEvents;
    int fCharge;
  };
}
#endif

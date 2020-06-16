#ifndef GFITTINGPARAMETERS_H
#define GFITTINGPARAMETERS_H

#include <iostream>
#include "EvtGenModels/EvtGMathFunctions.hh"
#include "Math/WrappedMultiTF1.h"
#include "Minuit2/Minuit2Minimizer.h"

namespace Gamapola
{
  class GFittingParameters{
  public:
    GFittingParameters();
    virtual ~GFittingParameters();
    void GSetFunction(double (Gamapola::GMathFunctions::*memFunction)(double* x, double* p),
                      const char* memFunctionName, const char* memFunctionClassName);
    void GSetFunctionParameters(const int& nPars, double* pars);
    void GSetFunctionVariables(const int& nVars, double* initialVarValues, 
                               double* lowVarsLimits, double* upVarsLimits, 
                               double* varsSteps, const char* kindsOfVars[], const char* namesOfVars[]);
    void GMinimize();
    const double& GGetMinFuncValue() const
    {
      return fMinFunctionValue;
    }
    const double* GGetFitParameters() const
    {
      return fFitParameters;
    }
    TF1* GGetFuncPointer() const
    {
      return fFuncPointer;
    }
    
  private:
    int fArg;
    double (Gamapola::GMathFunctions::*fMemFunction)(double* x, double* p);
    const char* fMemFunctionName;
    const char* fMemFunctionClassName;
    
    int fNParameters;
    double* fParsValues;
    
    int fNVariables;
    double* fInitialVarValues;
    double* fLowVarsLimits;
    double* fUpVarsLimits; 
    double* fVarsSteps;
    const char** fKindsOfVars;
    const char** fNamesOfVars;
    
    double fMinFunctionValue;
    const double* fFitParameters;
    TF1* fFuncPointer;
  };
}
#endif

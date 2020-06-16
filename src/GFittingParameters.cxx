#include "GFittingParameters.h"

namespace Gamapola{
  GFittingParameters::GFittingParameters():
  fArg(0)
  {
    std::cout << "GFittingParameters constructor calling. . ." << std::endl;
  }
  GFittingParameters::~GFittingParameters()
  {
    fMemFunctionName = NULL;
    fMemFunctionClassName = NULL;
    fParsValues = NULL;
    fInitialVarValues = NULL;
    fLowVarsLimits = NULL;
    fUpVarsLimits = NULL;
    fVarsSteps = NULL;
    fFitParameters = NULL;
    fFuncPointer = NULL;
    
    delete fMemFunctionName;
    delete fMemFunctionClassName;
    delete fParsValues;
    delete fInitialVarValues;
    delete fLowVarsLimits;
    delete fUpVarsLimits;
    delete fVarsSteps;
    delete fFitParameters;
    delete fFuncPointer;
//     std::cout << "GFittingParameters destructor calling. . ." << std::endl;
  }
  void GFittingParameters::GSetFunction(double (Gamapola::GMathFunctions::*memFunction)(double* x, double* p),
                      const char* memFunctionName, const char* memFunctionClassName)
  {
    fMemFunction = memFunction;
    fMemFunctionName = memFunctionName;
    fMemFunctionClassName = memFunctionClassName;
    std::cout << fMemFunctionName << "  " << fMemFunctionClassName << std::endl;
  }
  void GFittingParameters::GSetFunctionParameters(const int& nPars, double* pars)
  {
    fNParameters = nPars;
    fParsValues = pars;
//     for(int i = 0; i < fNParameters; i++)
//       std::cout << i << "   " << fParsValues[i] << std::endl;
  }
  
  void GFittingParameters::GSetFunctionVariables(const int& nVars, double* initialVarValues, 
                               double* lowVarsLimits, double* upVarsLimits, 
                               double* varsSteps, const char* kindsOfVars[], const char* namesOfVars[])
  {
    fInitialVarValues = initialVarValues;
    fLowVarsLimits = lowVarsLimits;
    fUpVarsLimits = upVarsLimits;
    fVarsSteps = varsSteps;
    fKindsOfVars = kindsOfVars;
    fNamesOfVars = namesOfVars;
    fNVariables = nVars;
    std::cout << fInitialVarValues[0] << "  " << fLowVarsLimits[0] << "  " << fUpVarsLimits[0] <<
    "  " << fVarsSteps[0] << "  " << fKindsOfVars[0] << "  " << fNamesOfVars[0] << std::endl;
  }
  
  void GFittingParameters::GMinimize()
  {
    GMathFunctions * mf = new GMathFunctions();
    TF1* pdf = new TF1(fMemFunctionName, mf, fMemFunction, 0, 0, fNParameters, fMemFunctionClassName, fMemFunctionName);
    pdf->SetParameters(fParsValues);
    ROOT::Minuit2::Minuit2Minimizer globalMin(ROOT::Minuit2::kMigrad);
    ROOT::Math::WrappedMultiTF1 g1(*pdf, fNVariables);
    globalMin.SetFunction(g1);
    globalMin.SetMaxIterations(1000);
    globalMin.SetTolerance(0.01);
    for(int i = 0; i < fNVariables; i++)
      globalMin.SetVariable(i, fNamesOfVars[i], fInitialVarValues[i],fVarsSteps[i]);
    globalMin.Minimize();
    double* varsArr = new double(fNVariables);
    for(int i = 0; i < fNVariables; i++)
      varsArr[i] = globalMin.X()[i];
    fMinFunctionValue = globalMin.MinValue();
    fFitParameters = varsArr;
    fFuncPointer = pdf;
  }
}
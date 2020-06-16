#ifndef GSYMBOLICEXPRESSIONSANALYZER_H
#define GSYMBOLICEXPRESSIONSANALYZER_H

#include <iostream>
// #include "GInterfaceForSymbolicComputations.h"
#include "GSymbolicMathFunctions.h"
#include "GExCompiler.h"
#include <fstream>

namespace Gamapola
{
  class GSymbolicExpressionsAnalyzer: 
//   virtual public GInterfaceForSymbolicComputations,
  virtual public GSymbolicMathFunctions{
  public:
    GSymbolicExpressionsAnalyzer(const std::string& fileName);
    virtual ~GSymbolicExpressionsAnalyzer();
    void GSetSymbolicKinematicalVariables(ex* symKinVars, double* valsK, int* flagsK, const int& nVars);
    void GSetSymbolicModelParameters(const int& nPars, const char** resNames, ex* modelPars, double* valsM, int* charges);
    void GGiNACAnalyzing();
    void GTesting();
    void GSetReadWriteStatus(bool readwrite);
    void GSetRearchivateStatus(const char* reachivationMode);
    const GFUNCP_CUBA& GGetPointersOnFunctionsOfKinPars() const
    {
      return fnormIntsCalc;
    }
    const GFUNCP_CUBA& GGetPointersOnFunctionsOfModelPars() const
    {
      return fnormIntsCalc2;
    }
    void GArchivateExpressions(ex pdfSymb, const char* templName, const char* archName, const char* archName2);
    std::vector<ex> GUnArchivateExpressions(lst syms, const char* templName, const char* archName);
    const int& GGetNumberOfNormalizationIntegrals() const
    {
      return fNumberOfNormalizationIntegrals;
    }
  private:
//     ex* fSymbKinVars;
//     ex* fSymbModelPars;
//     lst fListOfKinVars;
//     lst fListOfModPars;
    int fArg;
    bool fReadWriteStatus;
    const char* fRearchivationMode;
    std::string fFileName;
//     GFUNCP_CUBA fnormIntsCalc;
//     GFUNCP_CUBA fnormIntsCalc2;
//     ex fpdfSymb;
//     lst fCh;
//     std::vector<ex> fExpressions;
//     std::vector<std::complex<double> > fModelParsVals;
//     int fNumberOfNormalizationIntegrals;
  };
}
#endif

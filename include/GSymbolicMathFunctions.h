#ifndef GSYMBOLICMATHFUNCTIONS_H
#define GSYMBOLICMATHFUNCTIONS_H

#include <iostream>
#include <cmath>
#include <TMath.h>
#include <complex>
#include <cstring>
#include "Math/WrappedMultiTF1.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "EvtGenModels/EvtGPhysicsPDG.hh"
#include "ginac/ginac.h"
#include "EvtGenModels/EvtGMathFunctions.hh"
#include "GExCompiler.h"

using namespace GiNaC;
namespace Gamapola
{
  class GSymbolicMathFunctions: virtual public GMathFunctions{
  public:
    GSymbolicMathFunctions();
    virtual ~GSymbolicMathFunctions();
    void GSymbQPMC(ex* p);
    void GSymbFormFactors(ex* p);
    void GSymbCurrents(ex* p);
  protected:
    ex fcomp[2][2][2][2][2][2];
    ex fc[2][2][2][3][5][2];
    ex fph[2][2][3][2][2];
    ex fsubRes[2][2][3][5];
    ex fij[2][2][5];
    ex fres[2][5];
    ex frl[2];
    ex* fSymbKinVars;
    ex* fSymbModelPars;
    lst fListOfKinVars;
    lst fListOfModPars;
    int fArg;
    GFUNCP_CUBA fnormIntsCalc;
    GFUNCP_CUBA fnormIntsCalc2;
    ex fpdfSymb;
    lst fCh;
    lst fkk;
    lst foriginal;
    lst fsyms;
    lst fmod;
    std::vector<ex> fExpressions;
    std::vector<std::complex<double> > fModelParsVals;
    std::vector<std::complex<double> > fKinVals;
    int fNumberOfNormalizationIntegrals;
    int fNModelPars;
  };
}
#endif

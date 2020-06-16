#ifndef GEVENTSANALYZER_H
#define GEVENTSANALYZER_H

#include <iostream>
#include "GInterfaceForMinimization.h"

namespace Gamapola
{
  class GEventsAnalyzer: virtual public GInterfaceForMinimization{
  public:
    GEventsAnalyzer(const int& charge);
    ~GEventsAnalyzer();
    void GAnalyzeEvents();
    void GGenerateEvents(const double& low = std::numeric_limits<float>::min(), const double& up = std::numeric_limits<float>::max());
    
    int GInitializeSymbolicComputation(const int& nPars, const char** resNames, ex* modelPars, double* valsM, int* charges);
    void GWriteNormalizationIntegrals(const char* filename, const char* reachivationMode, const int& nPars, const char** resNames, ex* modelPars, double* valsM, int* charges);
    void GReadNormalizationIntegrals(const char* filename, const int& nPars, const char** resNames, ex* modelPars, double* valsM, int* charges);
    double GLogLikelihoodPDF1(double* x, double* p);
    double GPDF(double* x, double* p);
    double GDenominator(double* x);
    double GMatrixElementSquareNormalized(double* x, double* p);
    const std::vector<double>& GGetFitParameters() const { return fFitParameters; }
    const std::vector<double>& GGetFitParametersErrors() const { return fFitParametersErrors; }
    std::vector<double> GGenerateModelPars();
    std::vector<std::complex<double>> EvaluateModelCoeffs(double* modelPars);
    std::vector<std::complex<double>> EvaluateKinematicsCoeffs(double*x, double*p);
  private:
    void GMinimizeWithMinuit();
    int fArg;
    std::shared_ptr<GMathFunctions> fMf;
    std::vector<double> fFitParameters;
    std::vector<double> fFitParametersErrors;
  };
}
#endif

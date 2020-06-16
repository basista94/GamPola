#include "GInterfaceForMinimization.h"

namespace Gamapola{
  GInterfaceForMinimization::GInterfaceForMinimization():
  fArg(0),
  fMemFunction(nullptr),
  fMf(new GMathFunctions()),
  fSea(nullptr),
  fPointersOnFunctionsOfKinPars(nullptr),
  fPointersOnFunctionsOfModelPars(nullptr),
  fMemFunctionName(nullptr),
  fMemFunctionClassName(nullptr),
  fNameOfMinimizer(nullptr),
  fReachivationMode(""),
  
  fNParameters(0),
  fParsValues(0),
    
  fNVariables(0),
  fInitialVarValues(0),
  fLowVarsLimits(0),
  fUpVarsLimits(0),
  fVarsSteps(0),
  fKindsOfVars(0),
  fNamesOfVars(0),
    
  fMinFunctionValue(0),
  fMaxFunctionValue(0),
  fFitParameters(0),
  fCovarianceMatrix(0),
  fFuncPointer(0), 
  fErrorMatrix(0),
  fparNames(0),
  fStatus(0),
  fIntegralOfPDF(0),
  fNormalizationIntegrals(nullptr),
  fNIntegrals(0),
  fReadWrite(false),
  fFileWithIntegrals(nullptr),
  fNRuns(1e6),
  fEdm(0),
  fRandomize(0),
  fNEvents(0)
  {
    std::cout << "GInterfaceForMinimization constructor calling. . ." << std::endl;
  }
  GInterfaceForMinimization::~GInterfaceForMinimization()
  {
    delete [] fFitParameters;
//     if(fFuncPointer)
//       delete fFuncPointer;
    delete [] fCovarianceMatrix;
//     std::cout << "GInterfaceForMinimization destructor calling. . ." << std::endl;
  }
  
  void GInterfaceForMinimization::GSetEventsNumber(const int& nEvents)
  {
    fNEvents = nEvents;
  }
  void GInterfaceForMinimization::GSetFunction(double (Gamapola::GMathFunctions::*memFunction)(double* x, double* p))
  {
    fMemFunction = memFunction;
  }
  void GInterfaceForMinimization::GSetFunctionParameters(const int& nPars, double* pars)
  {
    fNParameters = nPars;    
    fParsValues.insert(fParsValues.end(), &pars[0], &pars[nPars]);
  }
  
  void GInterfaceForMinimization::GSetFunctionVariables(const int& nVars, double* initialVarValues, 
                               double* lowVarsLimits, double* upVarsLimits, 
                               double* varsSteps, const std::vector<std::string>& kindsOfVars, const std::vector<std::string>& namesOfVars, bool randomize)
  {
    fInitialVarValues.insert(fInitialVarValues.end(), &initialVarValues[0], &initialVarValues[nVars]);
    fLowVarsLimits.insert(fLowVarsLimits.end(), &lowVarsLimits[0], &lowVarsLimits[nVars]);
    fUpVarsLimits.insert(fUpVarsLimits.end(), &upVarsLimits[0], &upVarsLimits[nVars]);
    fVarsSteps.insert(fVarsSteps.end(), &varsSteps[0], &varsSteps[nVars]);
    fKindsOfVars = kindsOfVars;
    fNamesOfVars = namesOfVars;
    fNVariables = nVars;
    fFitParameters = new double[fNVariables];
    fCovarianceMatrix = new double[fNVariables*fNVariables];
    fRandomize = randomize;
    std::cout << this << "  " << fLowVarsLimits[0] << '\n';
    
//     for(auto i = 0; i < fNVariables; ++i)
//         std::cout << fNamesOfVars[i] << "  " << fKindsOfVars[i] << "  " << kindsOfVars[i] << '\n';
  }
  
  
  void GInterfaceForMinimization::GSetNameOfMinimizer(const char* nameOfMinimizer)
  {
    fNameOfMinimizer = nameOfMinimizer;
  }
  
  void GInterfaceForMinimization::GSetDecayMode(const int& charge)
  {
    fMf->GSetResonances(charge);
    fCharge = charge;
  }
  
  void GInterfaceForMinimization::GSetNRuns(const int& nRuns)
  {
    fNRuns = nRuns;
  }
  void GInterfaceForMinimization::GPrintErrorMatrix()
  {
    std::ofstream ofs;
    ofs.open("CovarianceMatrix.txt", std::ios::out);
//     ofs << "____________________________________________________Error Matrix_________________________________________________" << std::endl;
    std::cout << "EEEE" << '\n';
    std::vector<double> errorVector;
    for(int i = 0; i < fNVariables; ++i)
    {
      for(int j = 0; j < fNVariables; ++j)
        if(*(fCovarianceMatrix + i*fNVariables + j) != 0)
          errorVector.push_back(*(fCovarianceMatrix + i*fNVariables + j));
    }
    for(size_t i = 0; i < fErrorMatrix.size(); i++)
    {
      ofs << std::setw(10) << fparNames[i] << "   ";
      for(size_t j = 0; j < fErrorMatrix.size(); ++j)
      {
//           std::cout << i << " " << fErrorMatrix.size() << " " << j << "  " << errorVector.size() << '\n';
        fErrorMatrix[i][j] = errorVector[ i * fErrorMatrix.size() + j];
        if(i != j)
          fErrorMatrix[i][j] = fErrorMatrix[i][j] / 
          sqrt(errorVector[ i * fErrorMatrix.size() + i] * errorVector[ j * fErrorMatrix.size() + j]);
        else
          fErrorMatrix[i][j] = sqrt(errorVector[ i * fErrorMatrix.size() + j]);
        ofs << fErrorMatrix[i][j] << "   ";
      }
      ofs << std::endl;
    }
     ofs << "_________________________________________________________________________________________________________________" << std::endl;
     ofs.close();
  }
//   only applied for p.d.f > 0, for -2Log(L) won't work...
  void GInterfaceForMinimization::GMinimizeWithMonteCarlo() // non universal method!!!
  {
//     fFuncPointer = new TF1(fMemFunctionName, fMf, fMemFunction, 0, 0, fNParameters, fMemFunctionClassName, fMemFunctionName);
    double maxFunctionValue = -1e300;
    double minFunctionValue = 1e300;
    ROOT::Math::ParamFunctor functor = ROOT::Math::ParamFunctor(fMf,fMemFunction);    
    double pdfVal = 0.; 
    for(int i = 0; i < fNRuns; ++i)
    {
      
      auto&& inVars = fMf->GGetPhaseSpace();
//       std::cout << "Before" << '\n';
      pdfVal = functor(&inVars[0], &fParsValues[0]);
//       pdfVal = (*fMf)->GProcessingComputationOfPDF(fInitialVarValues, fParsValues);
      if(maxFunctionValue < pdfVal)
      {
        maxFunctionValue = pdfVal;
        std::cout << "MAX: " <<  maxFunctionValue << std::endl;
//         for(int iVars = 0; iVars < fNVariables; iVars++)
//           varsMax[iVars] = inVars[iVars];
      }
      if(minFunctionValue > pdfVal)
      {
        minFunctionValue = pdfVal;
        std::cout << "MIN: " << minFunctionValue << "  " << i << std::endl;
        for(int iVars = 0; iVars < fNVariables; iVars++)
        {
//           varsMin[iVars] = inVars[iVars];
          fFitParameters[iVars] = inVars[iVars];
        }
      }
    }
    fMinFunctionValue = minFunctionValue;
    fMaxFunctionValue = maxFunctionValue;
  }
}

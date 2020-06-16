#include "GEventsAnalyzer.h"

namespace Gamapola{
  GEventsAnalyzer::GEventsAnalyzer(const int& charge):GInterfaceForMinimization(),
  fArg(0),
  fMf(new GMathFunctions())
  {
    fMf->GSetResonances(charge);
    std::cout << "GEventsAnalyzer constructor calling. . ." << std::endl;
  }
  GEventsAnalyzer::~GEventsAnalyzer()
  {
//     std::cout << "GEventsAnalyzer destructor calling. . ." << std::endl;
  }
  double GEventsAnalyzer::GMatrixElementSquareNormalized(double* x, double* p)
  {      
    return -fMf->GProcessingComputationOfPDF(p,x);
  }
  double GEventsAnalyzer::GDenominator(double* x)
  {
      std::complex<double> sumC(0.,0.);
      std::complex<double> modCoeffs[fNIntegrals];
      std::vector<std::complex<double>> pC(26);
      for(int i = 0; i < 3; i++)
        pC[i] = std::complex<double>(x[i],0);
    
      pC[3] = std::complex<double>(cos(x[3]), sin(x[3]));
      pC[4] = std::complex<double>(cos(x[4]), sin(x[4]));
      pC[5] = std::complex<double>(cos(x[5]), sin(x[5]));
      pC[6] = std::complex<double>(x[6], x[10]);
      pC[7] = std::complex<double>(x[7],0.);
      pC[8] = std::complex<double>(x[8],x[9]);
      pC[9] = std::complex<double>(x[11],x[12]);
      pC[10] = std::complex<double>(x[13],x[14]);
      pC[11] = std::complex<double>(x[15],x[16]);
      
      pC[12] = std::complex<double>(x[17],x[18]);
      pC[13] = std::complex<double>(x[19],x[20]);
      pC[14] = std::complex<double>(x[21],x[22]);
      pC[15] = std::complex<double>(x[23],x[24]);
      pC[16] = std::complex<double>(x[25],x[26]);
      pC[17] = std::complex<double>(x[27],x[28]);
      
      pC[18] = std::complex<double>(x[29],x[30]);
      pC[19] = std::complex<double>(x[31],x[32]);
      pC[20] = std::complex<double>(x[33],x[34]);
      pC[21] = std::complex<double>(x[35],x[36]);
      pC[22] = std::complex<double>(x[37],x[38]);
      pC[23] = std::complex<double>(x[39],x[40]);
      pC[24] = std::complex<double>(x[41],x[42]);
      pC[25] = std::complex<double>(x[43],x[44]);
//     if(fModelCoeffs)
      fPointersOnFunctionsOfModelPars(&pC[0], &modCoeffs[0]);
//     else 
//       throw("AA");
      for(int i = 0 ; i < fNIntegrals; i++)
      {
        sumC += modCoeffs[i] * fNormalizationIntegrals[i];
//       std::cout << modCoeffs[i] << "  " << fNormalizationIntegrals[i] << std::endl;
      }
      return std::real(sumC);
  }
  double GEventsAnalyzer::GPDF(double* x, double* p)
  {
      std::complex<double> sumC(0.,0.);
      std::complex<double> modCoeffs[fNIntegrals];
      std::vector<std::complex<double>> pC(26);
      for(int i = 0; i < 3; i++)
        pC[i] = std::complex<double>(x[i],0);
    
      pC[3] = std::complex<double>(cos(x[3]), sin(x[3]));
      pC[4] = std::complex<double>(cos(x[4]), sin(x[4]));
      pC[5] = std::complex<double>(cos(x[5]), sin(x[5]));
      pC[6] = std::complex<double>(x[6], x[10]);
      pC[7] = std::complex<double>(x[7],0.);
      pC[8] = std::complex<double>(x[8],x[9]);
      pC[9] = std::complex<double>(x[11],x[12]);
      pC[10] = std::complex<double>(x[13],x[14]);
      pC[11] = std::complex<double>(x[15],x[16]);
      
      pC[12] = std::complex<double>(x[17],x[18]);
      pC[13] = std::complex<double>(x[19],x[20]);
      pC[14] = std::complex<double>(x[21],x[22]);
      pC[15] = std::complex<double>(x[23],x[24]);
      pC[16] = std::complex<double>(x[25],x[26]);
      pC[17] = std::complex<double>(x[27],x[28]);
      
      pC[18] = std::complex<double>(x[29],x[30]);
      pC[19] = std::complex<double>(x[31],x[32]);
      pC[20] = std::complex<double>(x[33],x[34]);
      pC[21] = std::complex<double>(x[35],x[36]);
      pC[22] = std::complex<double>(x[37],x[38]);
      pC[23] = std::complex<double>(x[39],x[40]);
      pC[24] = std::complex<double>(x[41],x[42]);
      pC[25] = std::complex<double>(x[43],x[44]);
//     if(fModelCoeffs)
      fPointersOnFunctionsOfModelPars(&pC[0], &modCoeffs[0]);
//     else 
//       throw("AA");
      for(int i = 0 ; i < fNIntegrals; i++)
      {
        sumC += modCoeffs[i] * fNormalizationIntegrals[i];
//       std::cout << modCoeffs[i] << "  " << fNormalizationIntegrals[i] << std::endl;
      }
      return GMatrixElementSquareNormalized(x,p) / std::real(sumC);
  }
  double GEventsAnalyzer::GLogLikelihoodPDF1(double* x, double* p)
  {
    std::complex<double> sumC(0.,0.);
    std::complex<double> modCoeffs[fNIntegrals];
    std::vector<std::complex<double>> pC(26);
    for(int i = 0; i < 3; i++)
      pC[i] = std::complex<double>(x[i],0);
    
    pC[3] = std::complex<double>(cos(x[3]), sin(x[3]));
    pC[4] = std::complex<double>(cos(x[4]), sin(x[4]));
    pC[5] = std::complex<double>(cos(x[5]), sin(x[5]));
    pC[6] = std::complex<double>(x[6], x[10]);
    pC[7] = std::complex<double>(x[7],0.);
    pC[8] = std::complex<double>(x[8],x[9]);
    pC[9] = std::complex<double>(x[11],x[12]);
    pC[10] = std::complex<double>(x[13],x[14]);
    pC[11] = std::complex<double>(x[15],x[16]);
    
    pC[12] = std::complex<double>(x[17],x[18]);
    pC[13] = std::complex<double>(x[19],x[20]);
    pC[14] = std::complex<double>(x[21],x[22]);
    pC[15] = std::complex<double>(x[23],x[24]);
    pC[16] = std::complex<double>(x[25],x[26]);
    pC[17] = std::complex<double>(x[27],x[28]);
    
    pC[18] = std::complex<double>(x[29],x[30]);
    pC[19] = std::complex<double>(x[31],x[32]);
    pC[20] = std::complex<double>(x[33],x[34]);
    pC[21] = std::complex<double>(x[35],x[36]);
    pC[22] = std::complex<double>(x[37],x[38]);
    pC[23] = std::complex<double>(x[39],x[40]);
    pC[24] = std::complex<double>(x[41],x[42]);
    pC[25] = std::complex<double>(x[43],x[44]);
//     if(fModelCoeffs)
      fPointersOnFunctionsOfModelPars(&pC[0], &modCoeffs[0]);
//     else 
//       throw("AA");
    for(int i = 0 ; i < fNIntegrals; i++)
    {
      sumC += modCoeffs[i] * fNormalizationIntegrals[i];
//       std::cout << modCoeffs[i] << "  " << fNormalizationIntegrals[i] << std::endl;
    }
    
    double loglik = 0;
    int nEvents = (int)p[0];
    int nPars = (int)p[1];
    for(int i = 0; i < nEvents; i++)
    {
      double pdf = GMatrixElementSquareNormalized(x,p+2+i*nPars);
//       std::cout << pdf << '\n';
//       std::cout << "SHOOO " << pdf << '\n';
//       std::cout << fabs(pdf) << '\n';
//       if(fabs(pdf)<1e-13 || pdf > 5)
//       {
//         std::cout << pdf << " " << *(p+2+i*nPars) << "  " << *(p+2+i*nPars+1) << "  " << *(p+2+i*nPars+2) << "  " <<
//         *(p+2+i*nPars+3) << "  " << *(p+2+i*nPars+4) << "  " << *(p+2+i*nPars+5) <<'\n';
//         pdf = 1.;
//       }
//       pdf = (pdf == 0.)?1:pdf;
      loglik += -2*std::log(pdf);
//       std::cout << loglik << '\n';
    }
//     std::cout << loglik << "  " << loglik + 2*nEvents*TMath::Log(std::real(sumC)) << '\n';
//     std::cout << " In loglik: " << std::endl;
//     std::cout << x[21] << "  " << x[22] << "   " << fCCouplings[3][0] << std::endl;
//     std::cout << x[23] << "  " << x[24] << "   " << fCCouplings[3][1] << std::endl;
//     std::cout << x[25] << "  " << x[26] << "   " << fCCouplings[4][0] << std::endl;
//     std::cout << x[27] << "  " << x[28] << "   " << fCCouplings[4][1] << std::endl;
//     std::cout << loglik + 2*nEvents*TMath::Log(std::real(sumC)) <<
//     pC[14] << "  " << pC[15] << "  " << pC[16] << "  "  << pC[17] << "  " << pC[3] << "  " << pC[4] << "  " << pC[5] << std::endl;
//     std::cout << loglik << "   " << 2*nEvents*TMath::Log(std::real(sumC)) << "  " << pC[7] << " denominator: " << sumC << '\n';
    
    return loglik + 2*nEvents*std::log(std::real(sumC));
  }
  void GEventsAnalyzer::GAnalyzeEvents()
  {
    std::complex<double> modelParsComplex[26];
    std::complex<double> modelCoeffs[fNIntegrals];    
    for(size_t iPar = 0; iPar < 3; iPar++)
      modelParsComplex[iPar] = std::complex<double>(fParsValues[iPar],0);
    
    modelParsComplex[3] = std::complex<double>(cos(fParsValues[3]),sin(fParsValues[3]));
    modelParsComplex[4] = std::complex<double>(cos(fParsValues[4]),sin(fParsValues[4]));
    modelParsComplex[5] = std::complex<double>(cos(fParsValues[5]),sin(fParsValues[5]));
    modelParsComplex[6] = std::complex<double>(fParsValues[6],fParsValues[10]);
    modelParsComplex[7] = std::complex<double>(fParsValues[7], 0);
    modelParsComplex[8] = std::complex<double>(fParsValues[8],fParsValues[9]);
    modelParsComplex[9] = std::complex<double>(fParsValues[11],fParsValues[12]);
    modelParsComplex[10] = std::complex<double>(fParsValues[13],fParsValues[14]);
    modelParsComplex[11] = std::complex<double>(fParsValues[15],fParsValues[16]);
    modelParsComplex[12] = std::complex<double>(fParsValues[17],fParsValues[18]);
    modelParsComplex[13] = std::complex<double>(fParsValues[19],fParsValues[20]);
    modelParsComplex[14] = std::complex<double>(fParsValues[21],fParsValues[22]);
    modelParsComplex[15] = std::complex<double>(fParsValues[23],fParsValues[24]);
    modelParsComplex[16] = std::complex<double>(fParsValues[25],fParsValues[26]);
    modelParsComplex[17] = std::complex<double>(fParsValues[27],fParsValues[28]);
    
    modelParsComplex[18] = std::complex<double>(fParsValues[29],fParsValues[30]);
    modelParsComplex[19] = std::complex<double>(fParsValues[31],fParsValues[32]);
    modelParsComplex[20] = std::complex<double>(fParsValues[33],fParsValues[34]);
    modelParsComplex[21] = std::complex<double>(fParsValues[35],fParsValues[36]);
    modelParsComplex[22] = std::complex<double>(fParsValues[37],fParsValues[38]);
    modelParsComplex[23] = std::complex<double>(fParsValues[39],fParsValues[40]);
    modelParsComplex[24] = std::complex<double>(fParsValues[41],fParsValues[42]);
    modelParsComplex[25] = std::complex<double>(fParsValues[43],fParsValues[44]);
    
    
    std::cout << " Number od norm integrals: " << fNIntegrals << std::endl;
    if(fPointersOnFunctionsOfModelPars)
      fPointersOnFunctionsOfModelPars(&modelParsComplex[0], &modelCoeffs[0]);

    if(!fNameOfMinimizer)
      GMinimizeWithMinuit();
  }
  
  std::vector<std::complex<double>> GEventsAnalyzer::EvaluateKinematicsCoeffs(double*x, double*p)
  {
      fMf->GProcessingComputationOfPDF(x,p);
      std::vector<std::complex<double>> kinCoeffs;
      kinCoeffs.resize(fNIntegrals);
      if(fPointersOnFunctionsOfKinPars)
      {
        fPointersOnFunctionsOfKinPars(&fMf->GGetKinematicalCoefficients()[0], &kinCoeffs[0]);
      }
      return kinCoeffs;
  }
  std::vector<std::complex<double>> GEventsAnalyzer::EvaluateModelCoeffs(double* modelPars)
  {
      std::vector<std::complex<double>> modelParsComplex;
      modelParsComplex.resize(26);
      
      std::vector<std::complex<double>> modelCoeffs;
      modelCoeffs.resize(fNIntegrals);
      
      for(size_t iPar = 0; iPar < 3; iPar++)
          modelParsComplex[iPar] = std::complex<double>(modelPars[iPar],0);
      
      modelParsComplex[3] = std::complex<double>(cos(modelPars[3]),sin(modelPars[3]));
      modelParsComplex[4] = std::complex<double>(cos(modelPars[4]),sin(modelPars[4]));
      modelParsComplex[5] = std::complex<double>(cos(modelPars[5]),sin(modelPars[5]));
      modelParsComplex[6] = std::complex<double>(modelPars[6],modelPars[10]);
      modelParsComplex[7] = std::complex<double>(modelPars[7], 0);
      modelParsComplex[8] = std::complex<double>(modelPars[8],modelPars[9]);
      modelParsComplex[9] = std::complex<double>(modelPars[11],modelPars[12]);
      modelParsComplex[10] = std::complex<double>(modelPars[13],modelPars[14]);
      modelParsComplex[11] = std::complex<double>(modelPars[15],modelPars[16]);
      modelParsComplex[12] = std::complex<double>(modelPars[17],modelPars[18]);
      modelParsComplex[13] = std::complex<double>(modelPars[19],modelPars[20]);
      modelParsComplex[14] = std::complex<double>(modelPars[21],modelPars[22]);
      modelParsComplex[15] = std::complex<double>(modelPars[23],modelPars[24]);
      modelParsComplex[16] = std::complex<double>(modelPars[25],modelPars[26]);
      modelParsComplex[17] = std::complex<double>(modelPars[27],modelPars[28]);
      modelParsComplex[18] = std::complex<double>(modelPars[29],modelPars[30]);
      modelParsComplex[19] = std::complex<double>(modelPars[31],modelPars[32]);
      modelParsComplex[20] = std::complex<double>(modelPars[33],modelPars[34]);
      modelParsComplex[21] = std::complex<double>(modelPars[35],modelPars[36]);
      modelParsComplex[22] = std::complex<double>(modelPars[37],modelPars[38]);
      modelParsComplex[23] = std::complex<double>(modelPars[39],modelPars[40]);
      modelParsComplex[24] = std::complex<double>(modelPars[41],modelPars[42]);
      modelParsComplex[25] = std::complex<double>(modelPars[43],modelPars[44]);
      
      fPointersOnFunctionsOfModelPars(&modelParsComplex[0], &modelCoeffs[0]);
      
      return modelCoeffs;
  }
  void GEventsAnalyzer::GGenerateEvents(const double& low , const double& up)
  {    
    srand (time(0));
    int nEvents = 0, nAllEvents = 0;
    double ui;
    GMinimizeWithMonteCarlo();
    std::cout << -fMinFunctionValue << std::endl;
    std::cout.precision(16);
    double maxPDFValue = -fMinFunctionValue;
    fNormalizationIntegrals = new std::complex<double>[fNIntegrals];
    if(!fReadWrite)
      for(int i = 0; i < fNIntegrals; i++)
        fNormalizationIntegrals[i] = 0.;
    double pdf;
    std::complex<double> modelParsComplex[26];
    std::complex<double> kinCoeffs[fNIntegrals];
    std::complex<double> modelCoeffs[fNIntegrals];    
    for(size_t iPar = 0; iPar < 3; iPar++)
      modelParsComplex[iPar] = std::complex<double>(fParsValues[iPar],0);
    
    modelParsComplex[3] = std::complex<double>(cos(fParsValues[3]),sin(fParsValues[3]));
    modelParsComplex[4] = std::complex<double>(cos(fParsValues[4]),sin(fParsValues[4]));
    modelParsComplex[5] = std::complex<double>(cos(fParsValues[5]),sin(fParsValues[5]));
    modelParsComplex[6] = std::complex<double>(fParsValues[6],fParsValues[10]);
    modelParsComplex[7] = std::complex<double>(fParsValues[7], 0);
    modelParsComplex[8] = std::complex<double>(fParsValues[8],fParsValues[9]);
    modelParsComplex[9] = std::complex<double>(fParsValues[11],fParsValues[12]);
    modelParsComplex[10] = std::complex<double>(fParsValues[13],fParsValues[14]);
    modelParsComplex[11] = std::complex<double>(fParsValues[15],fParsValues[16]);
    modelParsComplex[12] = std::complex<double>(fParsValues[17],fParsValues[18]);
    modelParsComplex[13] = std::complex<double>(fParsValues[19],fParsValues[20]);
    modelParsComplex[14] = std::complex<double>(fParsValues[21],fParsValues[22]);
    modelParsComplex[15] = std::complex<double>(fParsValues[23],fParsValues[24]);
    modelParsComplex[16] = std::complex<double>(fParsValues[25],fParsValues[26]);
    modelParsComplex[17] = std::complex<double>(fParsValues[27],fParsValues[28]);
    
    modelParsComplex[18] = std::complex<double>(fParsValues[29],fParsValues[30]);
    modelParsComplex[19] = std::complex<double>(fParsValues[31],fParsValues[32]);
    modelParsComplex[20] = std::complex<double>(fParsValues[33],fParsValues[34]);
    modelParsComplex[21] = std::complex<double>(fParsValues[35],fParsValues[36]);
    modelParsComplex[22] = std::complex<double>(fParsValues[37],fParsValues[38]);
    modelParsComplex[23] = std::complex<double>(fParsValues[39],fParsValues[40]);
    modelParsComplex[24] = std::complex<double>(fParsValues[41],fParsValues[42]);
    modelParsComplex[25] = std::complex<double>(fParsValues[43],fParsValues[44]);
    
    
    std::cout << " Number od norm integrals: " << fNIntegrals << std::endl;
    if(fPointersOnFunctionsOfModelPars)
    {
        std::cout << &modelParsComplex[0] << '\n';
        std::cout << &modelCoeffs[0] << '\n';
        std::cout << "dd" << '\n';
      fPointersOnFunctionsOfModelPars(&modelParsComplex[0], &modelCoeffs[0]);
    }
    clock_t tStart = clock();
    while(nEvents < fNEvents)
    {
      auto&& inVars = fMf->GGetPhaseSpace();
      ui = double(rand())/RAND_MAX;
      nAllEvents++;
      auto&& pdf = -fMf->GProcessingComputationOfPDF(&inVars[0], &fParsValues[0]);
      if(!fReadWrite && (std::sqrt(inVars[4]) >= 0 ) && (std::sqrt(inVars[4]) <= up))
      {
        fPointersOnFunctionsOfKinPars(&fMf->GGetKinematicalCoefficients()[0], kinCoeffs );
//         std::complex<double> pdfC(0.,0.);
        for(int i = 0; i < fNIntegrals; i++)
        {
//          pdfC += kinCoeffs[i] * modelCoeffs[i];
//           std::cout << kinCoeffs[i] << "  " << modelCoeffs[i] << std::endl;
          fNormalizationIntegrals[i] += kinCoeffs[i];
        }
//        std::cout << std::real(pdfC) << "   "/*<< (*fMf)->GGetKinematicalCoefficients()[0]*/ << "  " << pdf << std::endl;
      }
      if(ui*maxPDFValue < pdf)
      {
        nEvents++;
        if(nEvents%1000 == 0 )
                  std::cout << nEvents << std::endl;
      }
    }
    std::cout << "Total number of events: " << nAllEvents << std::endl;
    if(!fReadWrite)
    {
      std::ofstream ost;
      ost.open((std::string(fFileWithIntegrals)+"_dir/"+std::string(fFileWithIntegrals)).c_str());
      for(int iInts = 0; iInts < fNIntegrals; ++iInts)
        ost << std::real(fNormalizationIntegrals[iInts]) << "  " << std::imag(fNormalizationIntegrals[iInts]) << std::endl;
      ost.close();
    }
    std::cout << "Time taken: " << (double)(clock() - tStart)/CLOCKS_PER_SEC << "s " << std::endl;
  }
  int GEventsAnalyzer::GInitializeSymbolicComputation(const int& nPars, const char** resNames, ex* modelPars, double* valsM, int* charges)
  {
//     std::cout << "aloooooo *****************************************nrejfnejhebfjef" << '\n';
//     fSea = seaPointer;
    GSymbolicExpressionsAnalyzer* sea = new GSymbolicExpressionsAnalyzer(fFileWithIntegrals);
    sea->GSetSymbolicModelParameters(nPars, resNames, modelPars, valsM, charges);
    sea->GSetReadWriteStatus(fReadWrite);
    sea->GSetRearchivateStatus(fReachivationMode);
    sea->GGiNACAnalyzing();
    auto&& nIntegrals = sea->GGetNumberOfNormalizationIntegrals();
    fPointersOnFunctionsOfKinPars = sea->GGetPointersOnFunctionsOfKinPars();
    fPointersOnFunctionsOfModelPars = sea->GGetPointersOnFunctionsOfModelPars();
    
    delete sea;

    
    return nIntegrals;
  }
  void GEventsAnalyzer::GWriteNormalizationIntegrals(const char* filename, const char* rearchivationMode,
      const int& nPars, const char** resNames, ex* modelPars, double* valsM, int* charges)
  {
    fReadWrite = false;
    fFileWithIntegrals = filename;
    fReachivationMode = rearchivationMode;
    std::cout << rearchivationMode << '\n';
    fNIntegrals = GInitializeSymbolicComputation(nPars, resNames, modelPars, valsM, charges);
  }
  
  void GEventsAnalyzer::GReadNormalizationIntegrals(const char* filename,
      const int& nPars, const char** resNames, ex* modelPars, double* valsM, int* charges)
  {
    fReadWrite = true;
    fFileWithIntegrals = filename;
    
    int integralsNumber = 0;
    std::cout.precision(16);
    std::vector<std::complex<double> > partialIntegrals;
    std::cout << this << '\n';
      std::ifstream ist((std::string(fFileWithIntegrals)+"_dir/"+std::string(fFileWithIntegrals)).c_str());
      std::string line;
      std::cout << fFileWithIntegrals << std::endl;
      
      while (std::getline(ist, line))
      {
        double integRe = 0.0;
        double integIm = 0.0;
        std::istringstream iss(line);
        iss >> integRe >> integIm;
//         std::cout << integRe << "  " << integIm << "  " << integralsNumber << std::endl;
        partialIntegrals.push_back(std::complex<double>(integRe, integIm));
        integralsNumber++;
      }
//       ist.close();
      fNIntegrals = integralsNumber;
      fNormalizationIntegrals = new std::complex<double>[fNIntegrals];
      
      for(int i = 0; i < fNIntegrals; i++)
        fNormalizationIntegrals[i] = partialIntegrals[i];
    GInitializeSymbolicComputation(nPars, resNames, modelPars, valsM, charges);
  }
  std::vector<double> GEventsAnalyzer::GGenerateModelPars()
  {
      std::vector<double> model(fNVariables, 0.0);
      model = fInitialVarValues;
      for(int i = 0; i < fNVariables; ++i)
        if(std::strcmp(fKindsOfVars[i].c_str(), "fixed") != 0 && fRandomize)
        {
            model[i] = fLowVarsLimits[i] + (fUpVarsLimits[i] - fLowVarsLimits[i]) * (double)rand()/RAND_MAX;
        }
      return model;
  }
  void GEventsAnalyzer::GMinimizeWithMinuit()
  {
      fparNames.clear();
      fFitParameters.clear();
      fFitParametersErrors.clear();
    std::cout << "Pognali" << '\n';
    std::cout << fMf << '\n';
//     std::cout << *fMf << '\n';
    fFuncPointer = new TF1("GLogLikelihoodPDF1", this, &Gamapola::GEventsAnalyzer::GLogLikelihoodPDF1, 0, 0, fNParameters, "GEventsAnalyzer", "GLogLikelihoodPDF1");
    fFuncPointer->SetParameters(&fParsValues[0]);
    
    TMinuitMinimizer globalMin(ROOT::Minuit::kMigrad);
    ROOT::Math::WrappedMultiTF1 g1(*fFuncPointer, fNVariables);
    globalMin.SetFunction(g1);
    unsigned int numberOfFreedVariables = 0;
    std::cout << "minuit" << '\n';
    for(int i = 0; i < fNVariables; ++i)
    {
      if(std::strcmp(fKindsOfVars[i].c_str(), "fixed") != 0 && fRandomize)
      {
          fInitialVarValues[i] = fLowVarsLimits[i] + (fUpVarsLimits[i] - fLowVarsLimits[i]) * (double)rand()/RAND_MAX;
      }
//       std::cout << fNamesOfVars[i] << "   !" << fKindsOfVars[i] << "!  " << fInitialVarValues[i] << std::endl;
      if(std::strcmp(fKindsOfVars[i].c_str(), "limited") == 0)
      {
        globalMin.SetLimitedVariable(i, fNamesOfVars[i].c_str(), fInitialVarValues[i],fVarsSteps[i], fLowVarsLimits[i], fUpVarsLimits[i]);
        ++numberOfFreedVariables;
        fparNames.push_back(fNamesOfVars[i]);
      }
      else if(std::strcmp(fKindsOfVars[i].c_str(), "fixed") == 0)
        globalMin.SetFixedVariable(i, fNamesOfVars[i].c_str(), fInitialVarValues[i]);
      else
      {
        globalMin.SetVariable(i, fNamesOfVars[i].c_str(), fInitialVarValues[i],fVarsSteps[i]);
        fparNames.push_back(fNamesOfVars[i]);
        numberOfFreedVariables++;
      }
    }
    fErrorMatrix.resize(numberOfFreedVariables);
    for(size_t i = 0; i < numberOfFreedVariables; ++i)
      fErrorMatrix[i].resize(numberOfFreedVariables, 0.);
    globalMin.Minimize();
    globalMin.GetCovMatrix(fCovarianceMatrix);
    std::cout << std::scientific;
    GPrintErrorMatrix();
    globalMin.PrintResults();
    for(int i = 0; i < fNVariables; ++i)
    {
      fFitParameters.push_back(globalMin.X()[i]);
      fFitParametersErrors.push_back(globalMin.Errors()[i]);
//       std::cout << fFitParameters[i] << "  " << globalMin.Errors()[i] << '\n';
    }
    fMinFunctionValue = globalMin.MinValue();
    fStatus = globalMin.Status();
    fEdm = globalMin.Edm();
    if(fFuncPointer)
      delete fFuncPointer;
  }
}

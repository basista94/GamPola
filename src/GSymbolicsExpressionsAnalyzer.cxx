#include "GSymbolicExpressionsAnalyzer.h"
#include <dlfcn.h>
#include <sys/stat.h> 
#include <sys/types.h> 
#include <stdlib.h>

namespace Gamapola{
//   bool cl_inhibit_floating_point_underflow = true;
  GSymbolicExpressionsAnalyzer::GSymbolicExpressionsAnalyzer(const std::string& fileName): 
//   GInterfaceForSymbolicComputations(),
  GSymbolicMathFunctions(),
  fArg(0),
  fFileName(fileName)
  {
    std::cout << "GSymbolicExpressionsAnalyzer constructor calling. . ." << std::endl;
  }
  GSymbolicExpressionsAnalyzer::~GSymbolicExpressionsAnalyzer()
  {
    std::cout << "GSymbolicExpressionsAnalyzer destructor calling. . ." << std::endl;
  }
  void GSymbolicExpressionsAnalyzer::GSetSymbolicKinematicalVariables(ex* symKinVars, double* vals, 
                                                                      int* flags, const int& nVars)
  {
    fSymbKinVars = symKinVars;
    for(int i = 0; i < nVars; i++)
    {
      if(flags[i] == 1)
        fListOfKinVars.append(fSymbKinVars[i] == vals[i]);
    }
  }
  void GSymbolicExpressionsAnalyzer::GSetSymbolicModelParameters(const int& nPars, const char** resNames, ex* modelPars, double* valsM, int* charges)
  {
    GSetResonances(charges[0]);
//     double inVars3[] = {-0.263587,  std::atan(1),1.217044, 0.353116, 0., 2.89};
//     inVars3[4] = inVars3[5] - inVars3[2] - inVars3[3] + kMKaon * kMKaon + kMPion * kMPion + kMPion * kMPion;
    auto&& inVars = GGetPhaseSpace();
    double inVars3[] = {-0.4869676164756378,  0.9661186931133823,  0.5927130966324399, 0.4618010422370061,  0.7757106193798375,  1.547548107629819};
    std::cout << "ddffssf" << GProcessingComputationOfPDF(&inVars[0], valsM) << '\n';
    
    fNModelPars = nPars;
    fSymbModelPars = modelPars;
    
    std::cout << "******************************" << '\n';
    fModelParsVals.push_back(std::complex<double>(valsM[0],0.)); 
    fModelParsVals.push_back(std::complex<double>(valsM[1],0.)); 
    fModelParsVals.push_back(std::complex<double>(valsM[2],0.)); 
    fModelParsVals.push_back(std::complex<double>(cos(valsM[3]),sin(valsM[3]))); 
    fModelParsVals.push_back(std::complex<double>(cos(valsM[4]),sin(valsM[4]))); 
    fModelParsVals.push_back(std::complex<double>(cos(valsM[5]),sin(valsM[5]))); 
    fModelParsVals.push_back(std::complex<double>(valsM[6],valsM[10])); 
    fModelParsVals.push_back(std::complex<double>(valsM[7],0.)); 
    fModelParsVals.push_back(std::complex<double>(valsM[8],valsM[9])); 
    fModelParsVals.push_back(std::complex<double>(valsM[11],valsM[12]));
    fModelParsVals.push_back(std::complex<double>(valsM[13],valsM[14]));
    fModelParsVals.push_back(std::complex<double>(valsM[15],valsM[16]));
    fModelParsVals.push_back(std::complex<double>(valsM[17],valsM[18]));
    fModelParsVals.push_back(std::complex<double>(valsM[19],valsM[20]));
    fModelParsVals.push_back(std::complex<double>(valsM[21],valsM[22]));
    fModelParsVals.push_back(std::complex<double>(valsM[23],valsM[24]));
    fModelParsVals.push_back(std::complex<double>(valsM[25],valsM[26]));
    fModelParsVals.push_back(std::complex<double>(valsM[27],valsM[28]));

    
    fModelParsVals.push_back(std::complex<double>(valsM[29],valsM[30]));
    fModelParsVals.push_back(std::complex<double>(valsM[31],valsM[32]));
    fModelParsVals.push_back(std::complex<double>(valsM[33],valsM[34]));
    fModelParsVals.push_back(std::complex<double>(valsM[35],valsM[36]));
    fModelParsVals.push_back(std::complex<double>(valsM[37],valsM[38]));
    fModelParsVals.push_back(std::complex<double>(valsM[39],valsM[40]));
    fModelParsVals.push_back(std::complex<double>(valsM[41],valsM[42]));
    fModelParsVals.push_back(std::complex<double>(valsM[43],valsM[44]));

    fListOfModPars.append(fSymbModelPars[0] == valsM[0]);
    fListOfModPars.append(fSymbModelPars[1] == valsM[1]);
    fListOfModPars.append(fSymbModelPars[2] == valsM[2]);
    fListOfModPars.append(fSymbModelPars[3] == cos(valsM[3]) + I * sin(valsM[3]));
    fListOfModPars.append(fSymbModelPars[4] == cos(valsM[4]) + I * sin(valsM[4]));
    fListOfModPars.append(fSymbModelPars[5] == cos(valsM[5]) + I * sin(valsM[5]));
    fListOfModPars.append(fSymbModelPars[6] == valsM[6] + I*valsM[10]);
    fListOfModPars.append(fSymbModelPars[7] == valsM[7]);
    fListOfModPars.append(fSymbModelPars[8] == valsM[8]+I*valsM[9]);
    fListOfModPars.append(fSymbModelPars[9] == valsM[11]+I*valsM[12]);
    fListOfModPars.append(fSymbModelPars[10] == valsM[13]+I*valsM[14]);
    fListOfModPars.append(fSymbModelPars[11] == valsM[15]+I*valsM[16]);
    fListOfModPars.append(fSymbModelPars[12] == valsM[17]+I*valsM[18]);
    fListOfModPars.append(fSymbModelPars[13] == valsM[19]+I*valsM[20]);
    fListOfModPars.append(fSymbModelPars[14] == valsM[21]+I*valsM[22]);
    fListOfModPars.append(fSymbModelPars[15] == valsM[23]+I*valsM[24]);
    fListOfModPars.append(fSymbModelPars[16] == valsM[25]+I*valsM[26]);
    fListOfModPars.append(fSymbModelPars[17] == valsM[27]+I*valsM[28]);

    // *******************************extra parameters*************************************
    fListOfModPars.append(fSymbModelPars[18] == valsM[29]+I*valsM[30]);
    fListOfModPars.append(fSymbModelPars[19] == valsM[31]+I*valsM[32]);
    fListOfModPars.append(fSymbModelPars[20] == valsM[33]+I*valsM[34]);
    fListOfModPars.append(fSymbModelPars[21] == valsM[35]+I*valsM[36]);
    fListOfModPars.append(fSymbModelPars[22] == valsM[37]+I*valsM[38]);
    fListOfModPars.append(fSymbModelPars[23] == valsM[39]+I*valsM[40]);
    fListOfModPars.append(fSymbModelPars[24] == valsM[41]+I*valsM[42]);
    fListOfModPars.append(fSymbModelPars[25] == valsM[43]+I*valsM[44]);
  }
  void GSymbolicExpressionsAnalyzer::GTesting()
  {
      //Checking*********************************************************************************************************************
  double m12 = 0.692;
  double m13 = 0.992;
  double m2 = 0.140, m3 = 0.140, M = 1.272;
  double s12Min = m12*m12;
  double s23Min = (m2 + m3) * (m2 + m3);
  double s13Min = m13*m13;
  double sMax = (M)*(M);
  double inVars3[] = {-0.263587,  std::atan(1),1.217044, 0.353116, 0., 2.89};
  inVars3[4] = inVars3[5] - inVars3[2] - inVars3[3] + kMKaon * kMKaon + kMPion1 * kMPion1 + kMPion2 * kMPion2;
    
  double modelPars3[45] = {
          4.,       4.,    TMath::Pi()/3.,  TMath::Pi()/3., TMath::Pi()/3., 0.,      0.54,      0.789,     0.54,      0.98,          0.123,
    // gammaQPC(0), f2(1),     thetaK1(2),    phiDK*(3),      phiSRho(4),     phiDRho(5),       gg(6),   lambda(7),  ff(8) ,    ffIm(9),   ggIm(10)
         0.45,     -0.65,      0.98,       -0.32,         0.62,        0.91,
//      ff2(11)  ff2Im(12)  ff3(13)   ff3Im(14)    ff4(15)  ff4Im(16)
       1.0,          0.,             9.0,           -2.7,      
//     gg3Kstr(17)  gg3ImKstr(18)   gg3Rho(19)   gg3ImRho(20)  
       1.0,          0.,           1.4,            -0.93,
//     gg4Kstr(21)  gg4ImKstr(22)    gg4Rho(23)   gg4ImRho(24)
       1.0,          0.,              2.34,           1.7,
//   gg5Kstr(25)  gg5ImKstr(26)   gg5Rho(27)       gg5ImRho(28)
         1.0,          0.0,              1.0,           0.0,
         1.0,          0.0,              1.0,           0.0,
         1.0,          0.0,              1.0,           0.0,
         1.0,          0.0,              1.0,           0.0
  }; 
  const int nRes = 5;
  const char* resNames[nRes] = {"K1_1400", "K*_1410","K*_1680", "K2_1430","K1_1270"};
  int charges[nRes] = {0,0,0,0,0};
  ex p[26] =
  {
    symbol("gQPMC"), symbol("f2"),  symbol("thK1"), symbol("phDK"), 
    symbol("phSR"),  symbol("phDR"), symbol("gg"),       symbol("lambda"),
    symbol("ff"),    symbol("ff2"),  symbol("ff3"),      symbol("ff4"),
    symbol("gg3Kstr"), symbol("gg3Rho"), 
    symbol("gg4Kstr"), symbol("gg4Rho"),
    symbol("gg5Kstr"), symbol("gg5Rho"),
    
    symbol("kstr1270S"),    symbol("kstr1270D"),  symbol("kstr1400S"),      symbol("kstr1400D"),
    symbol("rho1270S"),    symbol("rho1270D"),  symbol("rho1400S"),      symbol("rho1400D")
  }; 
  GSetResonances(charges[0]);
  std::cout << "Numerical pdf: " << GProcessingComputationOfPDF(inVars3, modelPars3) << std::endl;
//   GSetSymbolicModelParameters(p, modelPars3, 18);
//   GSymbolicMathFunctions::GSymbCurrents(fSymbModelPars);
//   GGiNACAnalyzing();
  }
  void GSymbolicExpressionsAnalyzer::GGiNACAnalyzing()
  {     
//     std::cout << fListOfModPars << '\n';
    lst coeffsList, coeffsList2;
    std::cout.precision(16);
    GExCompiler* ec = new GExCompiler();
    std::cout << "*********************************" << fRearchivationMode << '\n';
    auto&& dirName = std::string(fFileName+"_dir");
    auto&& kinGar = dirName+"/"+fFileName+"kinematicsCoeffs.gar";
    auto&& couplingsGar = dirName+"/"+fFileName+"couplingCoeffs.gar";
        
    if (std::strcmp(fRearchivationMode, "rearchivate") == 0)
    {
        struct stat st;
        if(stat(dirName.c_str(), &st) == 0)
        {
            std::cout << dirName << '\n';
            system(std::string("rm -rf "+dirName).c_str());
        }
        if (mkdir(dirName.c_str(), 0777) == -1) 
            std::cerr << "Error :  " << strerror(errno) << '\n'; 
        else
            std::cout << "Directory created"; 
        
        for(auto i = 0; i < 26; ++i)
            std::cout << fSymbModelPars[i] << '\n';
        GSymbolicMathFunctions::GSymbCurrents(fSymbModelPars);
//         std::cout << fpdfSymb << '\n';
        GArchivateExpressions(fpdfSymb, "foo", kinGar.c_str(), couplingsGar.c_str());
    }
    if(!fReadWriteStatus)
    {
      GSymbolicMathFunctions::GSymbCurrents(fSymbModelPars);
      std::vector<ex> myCoeffs = GUnArchivateExpressions(fsyms, "foo", kinGar.c_str());
      std::vector<ex> basisSymb = GUnArchivateExpressions(fmod, "foo", couplingsGar.c_str());
      fNumberOfNormalizationIntegrals = int(myCoeffs.size());
      double sum = 0;
      for(size_t i = 0; i < basisSymb.size(); i++)
      {
//       std::cout << myCoeffs[i].subs(fCh) << "  " << basisSymb[i].subs(fListOfModPars) << std::endl;
        sum += ex_to<numeric>(myCoeffs[i].subs(fCh) * basisSymb[i].subs(fListOfModPars)).to_double();
      }
      std::cout << sum << "  "  << fpdfSymb.subs(fListOfModPars).subs(fCh) << std::endl;
      for (size_t i = 0; i < basisSymb.size(); ++i)
      {
        coeffsList.append(myCoeffs[i]);
        coeffsList2.append(basisSymb[i]);
      }
      std::vector<std::complex<double> > aaa;
      std::vector<std::complex<double> > bbb;
      aaa.resize(basisSymb.size(),0);
      bbb.resize(basisSymb.size(),0);
      std::cout << "Compilation phase ..." << std::endl;
      auto&& kinCxx = dirName+"/"+fFileName+"kinematicsCoeffs.cxx";
      auto&& couplingsCxx = dirName+"/"+fFileName+"couplingCoeffs.cxx";
      ec->GCompileEx(coeffsList,fsyms,fnormIntsCalc, kinCxx.c_str());
      ec->GCompileEx(coeffsList2,fmod,fnormIntsCalc2, couplingsCxx.c_str());
//       std::cout << "Done! " << std::endl;     
//       fnormIntsCalc(&fKinVals[0], &aaa[0]);
//       fnormIntsCalc2(&fModelParsVals[0], &bbb[0]);
//       std::complex<double> sum2 = 0.;
//       for(size_t i = 0; i < basisSymb.size(); i++)
//       {
//         sum2 += aaa[i] * bbb[i];
//       } 
//       std::cout << "Recovered pdf: " << sum2 << std::endl;
    }
    else
    {
        auto&& couplingsCxxSo = dirName+"/"+fFileName+"couplingCoeffs.cxx.so";
        auto&& kinematicsCxxSo = dirName+"/"+fFileName+"kinematicsCoeffs.cxx.so";
        fnormIntsCalc2 = (GFUNCP_CUBA)ec->GLinkSoFile(couplingsCxxSo.c_str(), 1);
        fnormIntsCalc = (GFUNCP_CUBA)ec->GLinkSoFile(kinematicsCxxSo.c_str(), 1);
    }
    
    delete ec;
  }
  void GSymbolicExpressionsAnalyzer::GArchivateExpressions(ex pdfSymb, const char* templName, const char* archName, const char* archName2)
  {
    std::cout << "Extraction phase ... " << std::endl;
//     std::cout << "  " << "pdfSymb" << pdfSymb << '\n';
    ex pdfEval = pdfSymb.subs(fCh);
    int nOp = 0;
    std::vector<ex> basisSymb;
    for (size_t i = 0; i != pdfEval.nops(); ++i)
    {
      ex inter = 1;
      for(size_t j = 0; j != pdfEval.op(i).nops(); ++j)
      {
        if(!is_a<numeric>(pdfEval.op(i).op(j)))
          inter *= pdfEval.op(i).op(j);
        nOp++;
      }
      basisSymb.push_back(inter);
//       std::cout << inter << '\n';
    }  
    std::vector<ex> vecE;
    std::vector<ex> myCoeffs;
    myCoeffs.resize(basisSymb.size(), 0);
    for(size_t i = 0; i !=pdfSymb.nops(); i++)
      vecE.push_back(pdfSymb.op(i));
    
    std::vector<ex> semiKin;
    std::vector<ex> semiModel;
    ex semK, semM;
    for(size_t i = 0; i < vecE.size(); i++)
    {
      ex iKin = 1;
      ex iModel = 1;
      ex iNumb = 1;
      semK = vecE[i].subs(fListOfModPars);
      semM = vecE[i].subs(fCh);
//       std::cout << semM << "  " << semK << std::endl;
      for(size_t j = 0; j != semK.nops(); ++j)
        if(!is_a<numeric>(semK.op(j)))
          iKin *= semK.op(j);
      for(size_t j = 0; j != semM.nops(); ++j)
        if(!is_a<numeric>(semM.op(j)))
          iModel *= semM.op(j); 
      for(size_t j = 0; j != vecE[i].nops(); ++j)
        if(is_a<numeric>(vecE[i].op(j)))
          iNumb = vecE[i].op(j);
      for(size_t j = 0; j < basisSymb.size(); j++)
        if(iModel.is_equal(basisSymb[j]))
          myCoeffs[j] += iKin * iNumb;
    }
    double sum = 0;
    for(size_t i = 0; i < myCoeffs.size(); i++)
    {
//       std::cout << myCoeffs[i].subs(fCh) << "  " << basisSymb[i].subs(fListOfModPars) << std::endl;
      sum += ex_to<numeric>(myCoeffs[i].subs(fCh) * basisSymb[i].subs(fListOfModPars)).to_double();
    }
    std::cout << sum << "  " << pdfSymb.subs(fCh).subs(fListOfModPars) << std::endl;
    std::cout << "Done! " << std::endl; 
    std::cout << "Archiving ..." << std::endl;
    archive a, a2;
    for(size_t i = 0; i < myCoeffs.size(); i++)
    {
      std::stringstream ss;//create a stringstream
      ss << i;//add number to the stream
      std::string s(templName);
      a.archive_ex(myCoeffs[i], (s+ss.str()).c_str());
      a2.archive_ex(basisSymb[i],(s+ss.str()).c_str());
//       std::cout << basisSymb[i] << '\n';
    }
    std::ofstream out(archName);
    out << a;
    out.close();
    std::ofstream out2(archName2);
    out2 << a2;
    out2.close();
    std::cout << "Done!" << std::endl;
  }
  std::vector<ex> GSymbolicExpressionsAnalyzer::GUnArchivateExpressions(lst syms, const char* templName, const char* archName)
  {
    std::cout << "Unarchiving..." << std::endl;
    std::vector<ex> expressions;
    archive a;
    std::ifstream in(archName);
    in >> a;
    ex exi = symbol("exi");
    int i = 0;
    int exe = 0;
    while(exe == 0)
    {
      try
      {
        std::string s(templName);
        std::stringstream ss;//create a stringstream
        ss << i;//add number to the stream
        exi = a.unarchive_ex(syms, (s+ss.str()).c_str());
        expressions.push_back(exi);
      }
      catch(...)
      {
        std::cout << "In archive " << archName << " were(was) " << i << " expression(s) " << std::endl;
        exe = i;
      }
      i++;
    }
    std::cout << "Done!" << std::endl;
    return expressions;
  }
  void GSymbolicExpressionsAnalyzer::GSetReadWriteStatus(bool readwrite)
  {
    fReadWriteStatus = readwrite;
  }
  
  void GSymbolicExpressionsAnalyzer::GSetRearchivateStatus(const char* reachivationMode)
  {
    fRearchivationMode = reachivationMode;
  }
}

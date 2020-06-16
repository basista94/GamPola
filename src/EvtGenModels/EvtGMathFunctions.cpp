#include "EvtGenModels/EvtGMathFunctions.hh"
#include <unordered_map>

namespace Gamapola{
  GMathFunctions::GMathFunctions():
  GKinematics(), 
  GBreitWigner(), 
  GCouplingConstants(), 
  GQuarkPairCreationModel(),
  GFormFactors(),
  GResVector(),
  GResTensor(),
  GResPseudoVector(),
  GResPseudoTensor(),
  GResPseudoTensor2Tensor(),
  GDecayRate()
  {
//     fGS = new GGrammSmidt();
    std::cout << "GMathFunctions constructor calling. . ." << std::endl;
  }
  GMathFunctions::~GMathFunctions()
  {
//     delete fGS;
//     std::cout << "GMathFunctions destructor calling. . ." << std::endl;
  }
  void GMathFunctions::GSetResonances(const int& charge)
  {
      const int nRes = 5;
      const char* resNames[] = {"K1_1400", "K*_1410","K*_1680", "K2_1430","K1_1270", "K2_1600", "K2_1770"};
      kMKaon = 0.49367700;
      kMPion1 = (charge == 1)?0.13957018:0.134976;
      kMPion2 = 0.13957018;
    fEpsilonL[0] = std::complex<double>(1., 0.)/sqrt(2.);
    fEpsilonL[1] = std::complex<double>(0., -1.)/sqrt(2.);
    fEpsilonL[2] = std::complex<double>(0., 0.)/sqrt(2.);
    
    fEpsilonR[0] = std::complex<double>(-1., 0.)/sqrt(2.);
    fEpsilonR[1] = std::complex<double>(0., -1.)/sqrt(2.);
    fEpsilonR[2] = std::complex<double>(0., 0.)/sqrt(2.);
    
    fEpsilon0[0] = std::complex<double>(0., 0.);
    fEpsilon0[1] = std::complex<double>(0., 0.);
    fEpsilon0[2] = std::complex<double>(1., 0.);
    fConstParsForI01.resize(10);
    fConstParsForMomA.resize(10);        
    fMomentaConst.resize(fConstParsForI01.size());
    std::vector<double> pars;
    std::vector<std::vector<double> > parsForRes;
//     fResNames = resNames;
    fNRes = nRes;
    
    GInitializeResVector(charge, "K1_1270");
    GInitializeResVector(charge, "K1_1400");
    GInitializeResTensor(charge, "K2_1430");
    GInitializeResPseudoVector(charge, "K*_1410");
    GInitializeResPseudoVector(charge, "K*_1680");
    GInitializeResPseudoTensor(charge, "K2_1600");
    GInitializeResPseudoTensor2Tensor(charge, "K2_1770");
    
    for(int i = 0; i < nRes; i++)
    {
      if(std::strcmp(resNames[i],"K2_1600") == 0 || std::strcmp(resNames[i],"K2_1770") == 0)
          continue;
      if(std::strcmp(resNames[i],"K1_1270") == 0 ||
         std::strcmp(resNames[i],"K1_1400") == 0 || 
         std::strcmp(resNames[i],"K*_1410") == 0 ||
         std::strcmp(resNames[i],"K*_1680") == 0 ||
         std::strcmp(resNames[i],"K2_1430") == 0 ||
         std::strcmp(resNames[i],"K2_1600") == 0 ||
         std::strcmp(resNames[i],"K2_1770") == 0)
      {
        double mK1;
        double gammaK1;
        double br1, br2;
        int iRes;
        
        if(std::strcmp(resNames[i],"K1_1270") == 0)
          iRes = 0;
        else if(std::strcmp(resNames[i],"K1_1400") == 0)
          iRes = 1;
        else if(std::strcmp(resNames[i],"K*_1410") == 0)
          iRes = 2;
        else if(std::strcmp(resNames[i],"K*_1680") == 0)
          iRes = 3;
        else if(std::strcmp(resNames[i],"K2_1430") == 0)
          iRes = 4;
        else if(std::strcmp(resNames[i],"K2_1600") == 0)
          iRes = 5;
        else if(std::strcmp(resNames[i],"K2_1770") == 0)
          iRes = 6;
        
        fCharge[iRes] = charge;
//         fCGKstr[iRes] = (fCharge[iRes]==0)?sqrt(2.)/3.:-2./3.;
//         fCGRho[iRes] = (fCharge[iRes]==0)?1./sqrt(3.):-1./sqrt(6.);
//         fDelta[iRes] = (fCharge[iRes]==0)?1.:0.;
        double coeff = 1.*sqrt((kMBmeson * kMBmeson - kMK1_1400 * kMK1_1400) * (kMBmeson * kMBmeson - kMK1_1400 * kMK1_1400)
        * (kMBmeson * kMBmeson - kMK1_1400 * kMK1_1400) / ( (kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270) * 
          (kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270) * (kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270) ) );
        fK1[1] = coeff;
        fK1[0] = 1.;
        fFF[2][0] = kFF1Kst_1410;
        fFF[2][1] = kFF2Kst_1410;
        fFF[3][0] = kFF1Kst_1680;
        fFF[3][1] = kFF2Kst_1680;
        fFF[4][0] = kFF1K2_1430;
        fFF[4][1] = kFF2K2_1430;
        
//         fLeft[0] = 1;
//         fLeft[1] = (kMBmeson * kMBmeson - kMK1_1400 * kMK1_1400)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270);
//         fLeft[2] = (kMBmeson * kMBmeson - kMKst_1410 * kMKst_1410)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270);
//         fLeft[3] = (kMBmeson * kMBmeson - kMKst_1680 * kMKst_1680)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270);
//         fLeft[4] = 1;
        
        fLeft[0] = 1;
        fLeft[1] = (kMBmeson * kMBmeson - kMK1_1400 * kMK1_1400)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270);
        fLeft[2] = std::pow((kMBmeson * kMBmeson - kMKst_1410 * kMKst_1410)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270),1.5);
        fLeft[3] = std::pow((kMBmeson * kMBmeson - kMKst_1680 * kMKst_1680)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270),1.5);
        fLeft[4] = std::pow((kMBmeson * kMBmeson - kMK2_1430 * kMK2_1430)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270),1.5);
        fLeft[5] = std::pow((kMBmeson * kMBmeson - kMK2_1600 * kMK2_1600)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270),1.5);
//         std::cout << fLeft[4] << '\n';
        fRight[0] =  1;
        fRight[1] =  (kMBmeson * kMBmeson - kMK1_1400 * kMK1_1400)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270);
        fRight[2] = -std::pow((kMBmeson * kMBmeson - kMKst_1410 * kMKst_1410)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270),1.5);
        fRight[3] = -std::pow((kMBmeson * kMBmeson - kMKst_1680 * kMKst_1680)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270),1.5);
        fRight[4] = -std::pow((kMBmeson * kMBmeson - kMK2_1430 * kMK2_1430)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270),1.5);
        fRight[5] = std::pow((kMBmeson * kMBmeson - kMK2_1600 * kMK2_1600)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270),1.5);
    
//         fRight[0] =  1;
//         fRight[1] =  (kMBmeson * kMBmeson - kMK1_1400 *  kMK1_1400)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270);
//         fRight[2] = -(kMBmeson * kMBmeson - kMKst_1410 *  kMKst_1410)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270);
//         fRight[3] = -(kMBmeson * kMBmeson - kMKst_1680 * kMKst_1680)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270);
//         fRight[4] = -1;
        
        fDelta2[iRes] = std::complex<double>(1.,0.);
//         fIndex.push_back(iRes);
        fIndex[iRes] = iRes;
        std::cout << "Number of resonances: " << nRes << "  " << fCharge[iRes] << std::endl;
        fIndexFor2Decays.push_back(iRes);
        GCouplingConstantsCalculation();
        switch(iRes)
        {
          case 0:
            mK1 = kMK1_1270;
            gammaK1 = kGammaK1_1270;
            break;
          case 1:
            mK1 = kMK1_1400;
            gammaK1 = kGammaK1_1400;
            break;
          case 2:
            mK1 = kMKst_1410;
            gammaK1 = kGammaKst_1410;
            break;
          case 3:
            mK1 = kMKst_1680;
            gammaK1 = kGammaKst_1680;
            break;
          case 4:
            mK1 = kMK2_1430;
            gammaK1 = kGammaK2_1430;
            break;
          case 5:
            mK1 = kMK2_1600;
            gammaK1 = kGammaK2_1600;
            break;
          case 6:
            mK1 = kMK2_1770;
            gammaK1 = kGammaK2_1770;
            break;
        }          
          pars.push_back(kMPion1); pars.push_back(kMK0star_892); pars.push_back(mK1); parsForRes.push_back(pars);
          pars.clear();
          pars.push_back(kMKaon); pars.push_back(kMRho0_775); pars.push_back(mK1); parsForRes.push_back(pars);
          pars.clear();  
          pars.push_back(kMPion1); pars.push_back(kMK0star_1430); pars.push_back(mK1);
          parsForRes.push_back(pars);
          
          fConstParsForMomA[iRes] = parsForRes;
          pars.clear();
          parsForRes.clear();

          fResMasses[iRes] = mK1;
          fResWidths[iRes] = gammaK1;
//           fResBrs[iRes] = 
          fMomentaConst[iRes].resize(fConstParsForMomA[iRes].size());

        if(iRes == 0)
        {
          fFunctionPointersKstrNR[0][0] = &GInterfaceForMathFunctions::GMSKstarNR1270;
          fFunctionPointersKstrNR[0][1] = &GInterfaceForMathFunctions::GMDKstarNR1270;
          fFunctionPointersRhoNR[0][0] = &GInterfaceForMathFunctions::GMSRhoNR1270;
          fFunctionPointersRhoNR[0][1] = &GInterfaceForMathFunctions::GMDRhoNR1270;
          fFunctionPointersKstr[0][0] = &GInterfaceForMathFunctions::GMSKstar1270;
          fFunctionPointersKstr[0][1] = &GInterfaceForMathFunctions::GMDKstar1270;
          fFunctionPointersRho[0][0] = &GInterfaceForMathFunctions::GMSRho1270;
          fFunctionPointersRho[0][1] = &GInterfaceForMathFunctions::GMDRho1270;
        }
        else if(iRes == 1)
        {
          fFunctionPointersKstrNR[1][0] = &GInterfaceForMathFunctions::GMSKstarNR1400;
          fFunctionPointersKstrNR[1][1] = &GInterfaceForMathFunctions::GMDKstarNR1400;  
          fFunctionPointersRhoNR[1][0] = &GInterfaceForMathFunctions::GMSRhoNR1400;
          fFunctionPointersRhoNR[1][1] = &GInterfaceForMathFunctions::GMDRhoNR1400;      
          fFunctionPointersKstr[1][0] = &GInterfaceForMathFunctions::GMSKstar1400;
          fFunctionPointersKstr[1][1] = &GInterfaceForMathFunctions::GMDKstar1400;    
          fFunctionPointersRho[1][0] = &GInterfaceForMathFunctions::GMSRho1400;
          fFunctionPointersRho[1][1] = &GInterfaceForMathFunctions::GMDRho1400;
        }
        
        GKinematicsConstanstsCalculation(mK1, iRes);
        GQPCMConstanstsCalculation(iRes);
      }
    }
  }
  double GMathFunctions::GOmegaPDF1(double* x, double* p)
  {
    return fOmega;
  }
  void GMathFunctions::GSetupDefaultParameters(double* p)
  {
      p[2] = 1.042; p[3] = 0.; p[4] = 0.; p[5] = 0.; 
      
      p[29] = 1.; p[30] = 0.; p[31] = 1.; p[32] = 0.; p[37] = 1.; p[38] = 0.; p[39] = 1.; p[40] = 0.;
      
      p[33] = 1.; p[34] = 0.; p[35] = 1.; p[36] = 0.; p[41] = 1.; p[42] = 0.; p[43] = 1.; p[44] = 0.;
      p[8] = 0.47; p[9] = 0.;
      
      p[11] = 0.78; p[12] = 0.; p[17] = 1.; p[18] = 0.; p[19] = 1.; p[20] = 0.;
      
      p[13] = 1.24; p[14] = 0.; p[21] = 1.; p[22] = 0.; p[23] = 1.; p[24] = 0.;
      
      p[15] = 0.44; p[16] = 0.; p[25] = 1.; p[26] = 0.; p[27] = 1.; p[28] = 0;
  }
  void GMathFunctions::GSetupParameters(double* p)
  {
//       Dict1D<double> parsFF;
//       parsFF["cosThK1"] = cos(p[2]);
//       parsFF["sinThK1"] = sin(p[2]);
//       parsFF["expPhaseSKstar_Re"] = cos(p[3]);
//       parsFF["expPhaseSKstar_Im"] = sin(p[3]);
//       parsFF["expPhaseSRho_Re"] = cos(p[4]);
//       parsFF["expPhaseSRho_Im"] = sin(p[4]);
//       parsFF["expPhaseDRho_Re"] = cos(p[5]);
//       parsFF["expPhaseDRho_Im"] = sin(p[5]);
//       
//       Dict1D<double> parsK_1270;
//       parsK_1270["K*_S_Re"] = p[29];
//       parsK_1270["K*_S_Im"] = p[30];
//       parsK_1270["K*_D_Re"] = p[31];
//       parsK_1270["K*_D_Im"] = p[32];
//       parsK_1270["rho_S_Re"] = p[37];
//       parsK_1270["rho_S_Im"] = p[38];
//       parsK_1270["rho_D_Re"] = p[39];
//       parsK_1270["rho_D_Im"] = p[40];
//       
//       Dict1D<double> parsK_1400;
//       parsK_1400["K*_S_Re"] = p[33];
//       parsK_1400["K*_S_Im"] = p[34];
//       parsK_1400["K*_D_Re"] = p[35];
//       parsK_1400["K*_D_Im"] = p[36];
//       parsK_1400["rho_S_Re"] = p[41];
//       parsK_1400["rho_S_Im"] = p[42];
//       parsK_1400["rho_D_Re"] = p[43];
//       parsK_1400["rho_D_Im"] = p[44];
//       parsK_1400["ff_Re"] = p[8];
//       parsK_1400["ff_Im"] = p[9];
//       
//       Dict1D<double> parsK_1410;
//       parsK_1410["ff_Re"] = p[11];
//       parsK_1410["ff_Im"] = p[12];
//       parsK_1410["K*_Re"] = p[17];
//       parsK_1410["K*_Im"] = p[18];
//       parsK_1410["rho_Re"] = p[19];
//       parsK_1410["rho_Im"] = p[20];
//       
//       Dict1D<double> parsK_1680;
//       parsK_1680["ff_Re"] = p[13];
//       parsK_1680["ff_Im"] = p[14];
//       parsK_1680["K*_Re"] = p[21];
//       parsK_1680["K*_Im"] = p[22];
//       parsK_1680["rho_Re"] = p[23];
//       parsK_1680["rho_Im"] = p[24];
//       
//       Dict1D<double> parsK2_1430;
//       parsK2_1430["ff_Re"] = p[15];
//       parsK2_1430["ff_Im"] = p[16];
//       parsK2_1430["K*_Re"] = p[25];
//       parsK2_1430["K*_Im"] = p[26];
//       parsK2_1430["rho_Re"] = p[27];
//       parsK2_1430["rho_Im"] = p[28];
            
      Dict1D<double> parsK_1600;
      parsK_1600["fV_K*"] = 1;
      parsK_1600["hV_K*"] = 0.2;
      parsK_1600["fV_Rho"] = 2;
      parsK_1600["hV_Rho"] = 4;
      
      Dict1D<double> parsK_1770;
      parsK_1770["fT_Re"] = 1;
      parsK_1770["fT_Im"] = 3.;
      parsK_1770["fV_K*"] = 1;
      parsK_1770["hV_K*"] = 0.2;
      parsK_1770["fV_Rho"] = 2;
      parsK_1770["hV_Rho"] = 4.0;
    
//       fPars2D["FF"] = parsFF;
//       fPars2D["K1_1270"] = parsK_1270;
//       fPars2D["K1_1400"] = parsK_1400;
//       fPars2D["K*_1410"] = parsK_1410;
//       fPars2D["K*_1680"] = parsK_1680;
//       fPars2D["K2_1430"] = parsK2_1430;
      fPars2D["K2_1660"] = parsK_1600;
      fPars2D["K2_1770"] = parsK_1770;
  }
  double GMathFunctions::GProcessingComputationOfPDF(double* x, double* p)
  {
    auto s = x[4];
    auto spi1pi2 = x[3];
    auto sKpi1 = x[2];
    
    GSetupParameters(p);
    double pdf;
    fDelta2[1] = std::complex<double>(p[8],p[9]);
    fDelta2[2] = std::complex<double>(p[11],p[12]);
    fDelta2[3] = std::complex<double>(p[13],p[14]);
    fDelta2[4] = std::complex<double>(p[15],p[16]);
    fCCouplings[2][0] = std::complex<double>(p[17],p[18]);
    fCCouplings[2][1] = std::complex<double>(p[19],p[20]);
    fCCouplings[3][0] = std::complex<double>(p[21],p[22]);
    fCCouplings[3][1] = std::complex<double>(p[23],p[24]);
    fCCouplings[4][0] = std::complex<double>(p[25],p[26]);
    fCCouplings[4][1] = std::complex<double>(p[27],p[28]);    
    fSinThK1 = sin(p[2]);
    fCosThK1 = cos(p[2]);
    auto sKpi2 = GSij(s, spi1pi2, sKpi1);
    fSij = GSij(s, sKpi2, sKpi1);
//     std::cout << " *********************************************** " << '\n';
//     std::cout << fSij << "  " << spi1pi2 << '\n';
//     std::cout << s << "  " << x[4] << '\n';
//     std::cout << spi1pi2 << "  " << x[3] << '\n';
//     std::cout << sKpi1 << "  " << x[2] << '\n';
//     std::cout << s - sKpi2 - sKpi1 + kMKaon * kMKaon + kMPion1 * kMPion1 + kMPion2 * kMPion2 << '\n';
    
//     std::cout << x[0] << "  " << x[1] << "  " << sKpi1 << "  " << spi1pi2 << "   " << sKpi2 << "  " << s << '\n';
    
    for(size_t iRes = 0; iRes < 2; iRes++)
    {
//       int iRes = fIndex[ires];
      fMomentaA[0] = GKinematics::GMomentaA(s, sKpi1, fConstParsForMomA[iRes][0][0]);
      fEnergyA[0] = GKinematics::GEnergyA(s, sKpi1, fConstParsForMomA[iRes][0][0]);
      fMomentaA[1] = GKinematics::GMomentaA(s,spi1pi2, fConstParsForMomA[iRes][1][0]);
      fEnergyA[1] = GKinematics::GEnergyA(s, spi1pi2, fConstParsForMomA[iRes][1][0]);
      fEnergyCoeff[0][0] = 8 * sqrt ( TMath::Pi() * TMath::Pi() * TMath::Pi() * 
        GKinematics::GEnergyA(s, sKpi1, fConstParsForMomA[iRes][0][0]) * 
        GKinematics::GEnergyA(s, fConstParsForMomA[iRes][0][0] * fConstParsForMomA[iRes][0][0], sqrt(sKpi1)) * sqrt( s ) );

      fEnergyCoeff[0][1] = 8 * sqrt ( TMath::Pi() * TMath::Pi() * TMath::Pi() * 
        GKinematics::GEnergyA(s, spi1pi2, fConstParsForMomA[iRes][1][0]) * 
        GKinematics::GEnergyA(s, fConstParsForMomA[iRes][1][0] * fConstParsForMomA[iRes][1][0], sqrt(spi1pi2)) * sqrt( s ) );
      for(size_t i = 0; i < 2/*fMomentaConst[fIndex[ires]].size()-1*/; i++)
      { 
        fI0[0][i] = GI0(iRes, i);
        fI1[0][i] = GI1(iRes, i);
        fMS[0][i] = GMS(p,iRes,i,0);
        fMD[0][i] = GMD(p,iRes,i,0);
//         std::cout <<"S,D waves: " << fMS[0][i] << "   " << fMD[0][i] << std::endl;
      }
    for(int iWave = 0; iWave < 2; iWave++)
    {
//       Non-relativistic matrix elements
      fMK0starNR[iRes][iWave] = (this->*fFunctionPointersKstrNR[iRes][iWave])(p,0);
      fMRho0NR[iRes][iWave] = (this->*fFunctionPointersRhoNR[iRes][iWave])(p,0);   
//       Relativistic matrix elements
      fMK0star[iRes][iWave] = (this->*fFunctionPointersKstr[iRes][iWave])(p,0);
      fMRho0[iRes][iWave] = (this->*fFunctionPointersRho[iRes][iWave])(p,0);  // why zero?
//       std::cout << "MK0star: " << fMK0star[iRes][iWave] << " for wave " << iWave << " res " << iRes << std::endl;
//       std::cout << "MRho0: " << fMRho0[iRes][iWave] << " for wave " << iWave << " res " << iRes << std::endl;
    }
      fA[iRes][0] = GAKstar(x,p,iRes, 0);
      fA[iRes][1] = GARho(x,p,iRes, 0);
      fB[iRes][0] = GBKstar(x,p,iRes,0);
      fB[iRes][1] = GBRho(x,p,iRes,0);
//       std::cout << "A for K*: " << fA[iRes][0] << " B for K*: " << fB[iRes][0] << std::endl;
//       std::cout << "A for Rho: " << fA[iRes][1] << " B for Rho: " << fB[iRes][1] << std::endl;
// // //  *********************************************************************************************************************  
    double sij = sKpi1;
    fMomentaA[0] = GKinematics::GMomentaA(s, sKpi2, fConstParsForMomA[iRes][0][0]);
    fEnergyA[0] = GKinematics::GEnergyA(s, sKpi2, fConstParsForMomA[iRes][0][0]);
    fEnergyCoeff[0][0] = 8 * sqrt ( TMath::Pi() * TMath::Pi() * TMath::Pi() * 
    GKinematics::GEnergyA(s, sKpi2, fConstParsForMomA[iRes][0][0]) * 
    GKinematics::GEnergyA(s, fConstParsForMomA[iRes][0][0] * fConstParsForMomA[iRes][0][0], sqrt(sKpi2)) * sqrt( s ) );
//     std::cout << std::endl;
//     std::cout << "After swaping arguments for Res: " << iRes << std::endl;
      for(size_t i = 0; i < 2/*fMomentaConst[fIndex[ires]].size()-1*/; i++)
      {
        fI0[0][i] = GI0(iRes, i);
        fI1[0][i] = GI1(iRes, i);
        fMS[0][i] = GMS(p,iRes,i,1);
        fMD[0][i] = GMD(p,iRes,i,1);
//         std::cout <<"S,D waves: " << fMS[0][i] << "   " << fMD[0][i] << std::endl;
      }
    for(int iWave = 0; iWave < 2; iWave++)
    {
//       Non-relativistic matrix elements
      fMK0starNR[iRes][iWave] = (this->*fFunctionPointersKstrNR[iRes][iWave])(p,1);
      fMRho0NR[iRes][iWave] = (this->*fFunctionPointersRhoNR[iRes][iWave])(p,1);   
//       Relativistic matrix elements
      fMK0star[iRes][iWave] = (this->*fFunctionPointersKstr[iRes][iWave])(p,1);
      fMRho0[iRes][iWave] = (this->*fFunctionPointersRho[iRes][iWave])(p,1);  // why zero?
//       std::cout << "MK0star: " << fMK0star[iRes][iWave] << " for wave " << iWave << " res " << iRes << std::endl;
//       std::cout << "MRho0: " << fMRho0[iRes][iWave] << " for wave " << iWave << " res " << iRes << std::endl;
    }
      x[2] = sKpi2;
      fA2[iRes][0] = GAKstar(x,p,iRes, 1);
      fB2[iRes][0] = GBKstar(x,p,iRes,1);
      x[2] = sij;
//     std::cout << "A2 for K*: " << fA2[iRes][0] << " B2 for K*: " << fB2[iRes][0] << std::endl; 
    }
    GCalculateKinematics(x,fCharge[0]);
    Dict2D<double> pars2d;
//     std::cout << fPars2D["K2_1660"]["fV_rho"] << '\n';
    pdf = GJ2(fPars2D,p,x);
//     std::cout << pdf << '\n';
    if(pdf != pdf)
    {
        std::cout << pdf << '\n';
        pdf = 0.;
    }
        return pdf;
        
  }
  std::vector<std::complex<double> > GMathFunctions::GGetKinematicalCoefficients() 
  {
    std::vector<std::complex<double> > kinCoeffs;
    kinCoeffs.resize(240, 0);
    int iii = 0;
//     std::cout << "**********************************************" << std::endl;
    for(size_t iTrig = 0; iTrig < 2; iTrig++) // {sin, cos}
      for(size_t iCurr = 0; iCurr < 2; iCurr++) // C1,2
        for(size_t iRL = 0; iRL < 2; iRL++) // {RL}
          for(size_t iMod = 0; iMod < 3; iMod++) // {K*ro}
            for(size_t iRes = 0; iRes < 5; iRes++) // {K1}
              for(size_t iPh = 0; iPh < 2; iPh++) // {phase S,D}
              {
                if(isnan(std::abs(cRL[iTrig][iRes][iMod][iCurr][iRL][iPh])))
                {
                  std::fill(kinCoeffs.begin(), kinCoeffs.end(), 0);
                  break;
                }
                kinCoeffs[iii] = cRL[iTrig][iRes][iMod][iCurr][iRL][iPh]; // cRL[iTrig][nRes][jC][0][0][iPh]
//                 std::cout << kinCoeffs[iii] << "  " << iTrig << "  " << iRes << "  " << iMod << "  " << iCurr
//                 <<"  " << iRL << "  " << iPh << std::endl;
                iii++;
              }     
    return kinCoeffs;
  }
}

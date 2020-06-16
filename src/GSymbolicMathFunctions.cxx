#include "GSymbolicMathFunctions.h"

namespace Gamapola{
  GSymbolicMathFunctions::GSymbolicMathFunctions(): GMathFunctions(),
  fArg(0),
  fNumberOfNormalizationIntegrals(0)
  {
    for(size_t iTrigo = 0; iTrigo < 2; iTrigo++) // {sin, cos}
      for(size_t iCurr = 0; iCurr < 2; iCurr++)
        for(size_t iRL = 0; iRL < 2; iRL++)
          for(size_t iMod = 0; iMod < 3; iMod++)
            for(size_t iRes = 0; iRes < 5; iRes++)
              for(size_t iPh = 0; iPh < 2; iPh++)
              {
                std::stringstream ss;//create a stringstream
                ss << iTrigo << iCurr << iRL << iMod << iRes << iPh;//add number to the stream
                fc[iTrigo][iCurr][iRL][iMod][iRes][iPh] = symbol(("c"+ss.str()).c_str());
//                 std::cout << c[i1][i2][i3][i4][i5][i6] << std::endl;
              }
  }
  GSymbolicMathFunctions::~GSymbolicMathFunctions()
  {
    std::cout << "GSymbolicMathFunctions destructor calling. . ." << std::endl;
  }
  void GSymbolicMathFunctions::GSymbQPMC(ex *p)
  {
    lst original, ch, syms, kk;
    for(size_t iTrigo = 0; iTrigo < 2; iTrigo++) // {sin, cos}
      for(size_t iCurr = 0; iCurr < 2; iCurr++)
        for(size_t iRL = 0; iRL < 2; iRL++)
          for(size_t iMod = 0; iMod < 3; iMod++)
            for(size_t iRes = 0; iRes < 5; iRes++)
              for(size_t iPh = 0; iPh < 2; iPh++)
              {
                double rI = (double)rand()/RAND_MAX;
                double iI = (double)rand()/RAND_MAX;
                syms.append(fc[iTrigo][iCurr][iRL][iMod][iRes][iPh]);
                original.append(fc[iTrigo][iCurr][iRL][iMod][iRes][iPh]);
                fKinVals.push_back(cRL[iTrigo][iRes][iMod][iCurr][iRL][iPh]);
                ch.append(fc[iTrigo][iCurr][iRL][iMod][iRes][iPh] == std::real(cRL[iTrigo][iRes][iMod][iCurr][iRL][iPh])
                +I*std::imag(cRL[iTrigo][iRes][iMod][iCurr][iRL][iPh]));
                
                kk.append(fc[iTrigo][iCurr][iRL][iMod][iRes][iPh] == std::real(cRL[iTrigo][iRes][iMod][iCurr][iRL][iPh])
                +I*std::imag(cRL[iTrigo][iRes][iMod][iCurr][iRL][iPh]));
//                 std::cout << cRL[iTrigo][iRes][iMod][iCurr][iRL][iPh] << " " << iTrigo << "  " <<  std::endl;
              }
    fCh = ch;
    fkk = kk;
    fsyms = syms;
    foriginal = original;
    auto extra = 0;
    
    for(size_t iCurr = 0; iCurr < 2; iCurr++)
      for(size_t iRL = 0; iRL < 2; iRL++)
      {
        extra = 0;
        for(size_t iMod = 0; iMod < 2; iMod++)
          for(size_t iRes = 0; iRes < 2; iRes++)
            for(size_t iPh = 0; iPh < 2; iPh++)
            {
              fph[iCurr][iRL][iMod][iRes][iPh] = fSymbModelPars[18+extra] * (fc[0][iCurr][iRL][iMod][iRes][iPh] * sin(2[fSymbModelPars]) + 
              fc[1][iCurr][iRL][iMod][iRes][iPh] * cos(2[fSymbModelPars])); // {sin}, {C1}, {Rl}, {K*ro}, {k1}, {phases}
//               std::cout << " Curr1: " << fph[iCurr][iRL][iMod][iRes][iPh].subs(fkk).subs(fListOfModPars) << std::endl;
              ++extra;
            }
      }
    // 0 --- K*, 1270, s
    // 1 --- K*, 1270, d
    // 2 --- K*, 1400, s
    // 3 --- K*, 1400, d
    // 4 --- rho, 1270, s
    // 5 --- rho, 1270, d
    // 6 --- rho, 1400, s
    // 7 --- rho, 1400, d
            
    for(size_t iCurr = 0; iCurr < 2; iCurr++)
      for(size_t iRL = 0; iRL < 2; iRL++)
        for(size_t iRes = 0; iRes < 2; iRes++)
          for(size_t iPh = 0; iPh < 2; iPh++)
          {
            fph[iCurr][iRL][2][iRes][iPh] = fc[0][iCurr][iRL][2][iRes][iPh]; // {sin}, {C1}, {Rl}, {K*ro}, {k1}, {phases}
//             std::cout << " Curr2: " << fph[iCurr][iRL][2][iRes][iPh].subs(fkk).subs(fListOfModPars) << std::endl;
          }        
  }
  void GSymbolicMathFunctions::GSymbFormFactors(ex* p)
  {
    GSymbQPMC(p);
    for(size_t iCurr = 0; iCurr < 2; iCurr++)
      for(size_t iRL = 0; iRL < 2; iRL++)
        for(size_t iRes = 0; iRes < 2; iRes++)
          {
            fsubRes[iCurr][iRL][0][iRes] = fph[iCurr][iRL][0][iRes][0] + fSymbModelPars[3] * fph[iCurr][iRL][0][iRes][1];
            
            fsubRes[iCurr][iRL][1][iRes] = fph[iCurr][iRL][1][iRes][0] * fSymbModelPars[4] + 
            fSymbModelPars[4] * fSymbModelPars[5] * fph[iCurr][iRL][1][iRes][1];
            
            fsubRes[iCurr][iRL][2][iRes] = fph[iCurr][iRL][2][iRes][0]; // {C1}, {Rl}, {K*ro}, {k1}, {phases}
            
//             std::cout << fsubRes[iCurr][iRL][0][iRes].subs(fkk).subs(fListOfModPars)
//               << "  " << fsubRes[iCurr][iRL][1][iRes].subs(fkk).subs(fListOfModPars)
//               << "  " << fsubRes[iCurr][iRL][2][iRes].subs(fkk).subs(fListOfModPars) << std::endl;
          }
    for(size_t iCurr = 0; iCurr < 2; iCurr++)
      for(size_t iRL = 0; iRL < 2; iRL++)
        for(size_t iMod = 0; iMod < 3; iMod++)
          for(size_t iRes = 2; iRes < 5; iRes++)
          { 
            fsubRes[iCurr][iRL][iMod][iRes] = fc[0][iCurr][iRL][iMod][iRes][0]; // {C1}, {Rl}, {K*ro}, {k1}, {phases}
//             std::cout << fsubRes[iCurr][iRL][iMod][iRes].subs(fkk).subs(fListOfModPars) << "  " <<
//             iCurr << "  " << iRL << "  " << iMod << "  " << iRes << std::endl;
//             std::cout << subRes[i1][i2][i3][0].subs(kk).subs(fListOfModPars)
//               << "  " << subRes[i1][i2][i3][1].subs(kk).subs(fListOfModPars)
//               << "  " << subRes[i1][i2][i3][2].subs(kk).subs(fListOfModPars) << std::endl;
          }
  }
  void GSymbolicMathFunctions::GSymbCurrents( ex* p)
  {
    GSymbFormFactors(p);
    lst mod;
    for(size_t iCurr = 0; iCurr < 2; iCurr++)
      for(size_t iRL = 0; iRL < 2; iRL++)
        {
          fij[iCurr][iRL][0] = fsubRes[iCurr][iRL][0][0] + fsubRes[iCurr][iRL][1][0] + fSymbModelPars[6] * fsubRes[iCurr][iRL][2][0];
          fij[iCurr][iRL][1] = fsubRes[iCurr][iRL][0][1] + fsubRes[iCurr][iRL][1][1];
          fij[iCurr][iRL][2] = fSymbModelPars[12] * fsubRes[iCurr][iRL][0][2] + fSymbModelPars[13] * fsubRes[iCurr][iRL][1][2];
          fij[iCurr][iRL][3] = fSymbModelPars[14] * fsubRes[iCurr][iRL][0][3] + fSymbModelPars[15] * fsubRes[iCurr][iRL][1][3];
          fij[iCurr][iRL][4] = fSymbModelPars[16] * fsubRes[iCurr][iRL][0][4] + fSymbModelPars[17] * fsubRes[iCurr][iRL][1][4];
//           std::cout << fij[iCurr][iRL][2].subs(fkk).subs(fListOfModPars) << "  " << iCurr << "  " << iRL << "  " << 2 <<
//           std::endl << fij[iCurr][iRL][3].subs(fkk).subs(fListOfModPars) << "  " << iCurr << "  " << iRL << "  " << 3 <<
//           std::endl << fij[iCurr][iRL][4].subs(fkk).subs(fListOfModPars) << "  " << iCurr << "  " << iRL << "  " << 4 << std::endl;
          // {C1}, {Rl}, {K*ro}, {k1}
        }
    for(size_t iRL = 0; iRL < 2; iRL++)
    {
      for(size_t iRes = 0; iRes < 2; iRes++)
      {
        fres[iRL][iRes] = fij[0][iRL][iRes] + fij[1][iRL][iRes]; // {C1}, {Rl}, {k1}
//         std::cout << fres[iRL][iRes].subs(fkk).subs(fListOfModPars) << "  " << iRL << "  " << iRes << std::endl;
        fres[iRL][iRes+2] = fij[0][iRL][iRes+2];
//         std::cout << fres[iRL][iRes+2].subs(fkk).subs(fListOfModPars) << "  " << iRL << "  " << iRes+2 << std::endl;
//         std::cout << fij[1][iRL][iRes] << "  " << iRes << std::endl;
      }
      fres[iRL][4] = fij[0][iRL][4] + fij[1][iRL][4]; // {C1}, {Rl}, {k1}
//       std::cout << fres[iRL][4].subs(fkk).subs(fListOfModPars) << "  " << iRL << "  " << 4 << std::endl;
    }
    for(size_t iRL = 0; iRL < 2; iRL++)
    {
      frl[iRL] = fres[iRL][0] + fSymbModelPars[8] * fres[iRL][1] + fSymbModelPars[9] * fres[iRL][2]+ 
      fSymbModelPars[10] * fres[iRL][3] + fSymbModelPars[11] * fres[iRL][4]; //  {Rl}, {k1}
//       std::cout << fSymbModelPars[9].subs(fListOfModPars) << std::endl;
//       std::cout <<  (/*fSymbModelPars[9] * */fres[iRL][2]).subs(fkk).subs(fListOfModPars)
//       << "   " <<  (/*fSymbModelPars[10] **/ fres[iRL][3]).subs(fkk).subs(fListOfModPars)
//       << "   " <<  (/*fSymbModelPars[11] * */fres[iRL][4]).subs(fkk).subs(fListOfModPars) << std::endl;
// //       std::cout << (fres[iRL][0] + fSymbModelPars[8] * fres[iRL][1]).subs(fkk).subs(fListOfModPars) << std::endl;
    }       
    fpdfSymb = (frl[0] * conjugate(frl[0]) + frl[1] * conjugate(frl[1]) - 
    fSymbModelPars[7] * (frl[1] * conjugate(frl[1]) - frl[0] * conjugate(frl[0]))).expand();  
    std::cout << fkk[0] << std::endl;
    std::cout << "Symbolic: " << fpdfSymb.subs(fkk).subs(fListOfModPars) <<"  " << fNModelPars << std::endl;
    for(size_t i = 0; i < fNModelPars; i++)
      mod.append(fSymbModelPars[i]);
    fmod = mod;
  }
}

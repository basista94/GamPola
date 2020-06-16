#include "EvtGenModels/EvtGResPseudoTensor2Tensor.hh"

namespace Gamapola{
    GResPseudoTensor2Tensor::GResPseudoTensor2Tensor(): GInterfaceForMathFunctions(),
    fCG_Kstr(-100),
    fCG_Rho(-100),
    fCoeffDelta(-100)
    {}
    GResPseudoTensor2Tensor::~GResPseudoTensor2Tensor()
    {
    }
    void GResPseudoTensor2Tensor::GInitializeResPseudoTensor2Tensor(const int& charge, const std::string& resName)
    {
      fCG_Kstr = (charge==0)?sqrt(2.)/3.:-2./3.;
      fCG_Rho = (charge==0)?1./sqrt(3.):-1./sqrt(6.);
      fCoeffDelta = (charge==0)?1.:0.;
      
      GInitializeResPseudoTensor(charge, resName);
      if(resName == "K2_1770")
      {
          Dict1D<double> database;
          database["mass"] = kMK2_1770;
          database["width"] = kGammaK2_1770;
          database["left"] = std::pow((kMBmeson * kMBmeson - kMK2_1770 * kMK2_1770)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270),1.5);
          database["right"] = std::pow((kMBmeson * kMBmeson - kMK2_1770 * kMK2_1770)/(kMBmeson * kMBmeson - kMK1_1270 * kMK1_1270),1.5);
          fRes[resName] = database;
      }
    }
    
    
    
    void GResPseudoTensor2Tensor::GNVec( Dict1D<double>& pars, double* x, int nRes, const std::string& resName)
    {
        for(auto&& iAxis = 0; iAxis < 3; ++iAxis)
        {
            fJVec3L[nRes][iAxis] = 0.;
            fJVec3R[nRes][iAxis] = 0.;
            fJVec3[nRes][iAxis] = 0.;
        }
        GMVec(pars, x, nRes, resName);
        
        auto&& pK2_1 = fqK + fqPi1;
        auto&& pKpK2_1 = fqK * pK2_1;
        auto&& sK2_1 = pK2_1 * pK2_1;
        auto&& e0pK = fEpsilon0[1] * fqK.get(1) + fEpsilon0[2] * fqK.get(2) + fEpsilon0[3] * fqK.get(3);
        auto&& e0pK2_1 = fEpsilon0[1] * pK2_1.get(1) + fEpsilon0[2] * pK2_1.get(2) + fEpsilon0[3] * pK2_1.get(3);
        
        auto&& pK2_2 = fqK + fqPi2;
        auto&& pKpK2_2 = fqK * pK2_2;
        auto&& sK2_2 = pK2_2 * pK2_2;
        auto&& e0pK2_2 = fEpsilon0[1] * pK2_2.get(1) + fEpsilon0[2] * pK2_2.get(2) + fEpsilon0[3] * pK2_2.get(3);
        
        auto&& mass = fRes[resName]["mass"];
        auto&& width = fRes[resName]["width"];
        auto&& leftCoeff = fRes[resName]["left"];
        auto&& rightCoeff = fRes[resName]["right"];
        auto&& s = (fqK + fqPi1 + fqPi2)*(fqK + fqPi1 + fqPi2);
        
        auto&& bw = std::complex<double>(pars["fT_Re"], pars["fT_Im"]) * GBWK1bw(s, mass, width) * (kMBmeson*kMBmeson-s);
//         std::cout << "here" << '\n';
//         for(auto&& iAxis = 0; iAxis < 3; ++iAxis)
//             std::cout << fJVec3L[nRes][iAxis] << "   " << fJVec3R[nRes][iAxis] << "  " << iAxis << '\n';
        
        for(auto&& iAxis = 0; iAxis < 3; ++iAxis)
        {
//             std::cout << iAxis << '\n';
            auto&& nKpi1 = fCG_Kstr * fCoeffDelta * ( e0pK * ( pK2_1.get(iAxis) * pKpK2_1 - sK2_1 * fqK.get(iAxis) ) + 
            e0pK2_1 * (fqK.get(iAxis) * (kMKaon * kMKaon * sK2_1 + 2 * pKpK2_1 * pKpK2_1)) / (3. * sK2_1) - e0pK2_1 * pKpK2_1 * fqK.get(iAxis) ) * GBWK1bw(sK2_1, kMK2_1430, kGammaK2_1430);
            
            auto&& nKpi2 = fCG_Kstr * ( e0pK * ( pK2_2.get(iAxis) * pKpK2_2 - sK2_2 * fqK.get(iAxis) ) + 
            e0pK2_2 * (fqK.get(iAxis) * (kMKaon * kMKaon * sK2_2 + 2 * pKpK2_2 * pKpK2_2)) / (3. * sK2_2) - e0pK2_2 * pKpK2_2 * fqK.get(iAxis) ) * GBWK1bw(sK2_2, kMK2_1430, kGammaK2_1430);
            
            fJVec3[nRes][iAxis] += (nKpi1+nKpi2) * bw;
        
            fJVec3L[nRes][iAxis] += fJVec3[nRes][iAxis] * leftCoeff;
            fJVec3R[nRes][iAxis] += fJVec3[nRes][iAxis] * rightCoeff;
        }
    }
}

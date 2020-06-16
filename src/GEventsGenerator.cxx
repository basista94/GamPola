#include "GEventsGenerator.h"
#include "TTree.h"
#include "TFile.h"

namespace Gamapola{
  GEventsGenerator::GEventsGenerator():GInterfaceForMinimization(),
  fArg(0),
//   
  fLowCosThetaLimit(0),
  fUpCosThetaLimit(0),
  fLowPhiLimit(0),
  fUpPhiLimit(0),
//     generated kinematic variables
  fCosThetaPDF(0),
  fPhiPDF(0),
  fPDF1(0),
  fOmega(0),
  fOmega2(0),
//     generated flat kinematic variables
  fCosThetaPDFFlat(0),
  fPhiPDFFlat(0),
  fPDF1Flat(0),
  fOmegaFlat(0),
  fOmega2Flat(0),
  fKaon_4V(new std::vector<GVector4D>),
  fPion1_4V(new std::vector<GVector4D>),
  fPion2_4V(new std::vector<GVector4D>),
  fPhoton_4V(new std::vector<GVector4D>)
  {
    std::cout << "GEventsGenerator constructor calling. . ." << std::endl;
  }
  GEventsGenerator::~GEventsGenerator()
  {
//     fNormalizationIntegrals = NULL;
    delete fNormalizationIntegrals;
    
    fCosThetaPDFFlat.clear();
    fPhiPDFFlat.clear();
    fPDF1Flat.clear();
    fOmegaFlat.clear();
    fOmega2Flat.clear();
    fCosThetaPDF.clear();
    fPhiPDF.clear();
    fPDF1.clear();
    fOmega.clear();
    fOmega2.clear();
//     std::cout << "GEventsGenerator destructor calling. . ." << std::endl;
  }
  std::vector<double> GEventsGenerator::GetSliceWithinBounds(const std::vector<double>& toCutVec, const std::vector<double>& controlVec,
                                             const double& lowBound, const double& upBound) const
  {
       std::vector<double> res;
       for(auto&& i = 0; i < static_cast<int>(toCutVec.size()); i++)
       {
//            std::cout << controlVec[i] << "  " << lowBound << "  " << upBound << '\n';
           if( (controlVec[i] <= upBound) && (controlVec[i] >= lowBound) )
           {
//                std::cout << toCutVec[i] << '\n';
               res.push_back(toCutVec[i]);
           }
       }
       return res;
  }
  
  void GEventsGenerator::GSetThetaLimits(const double& lowCosTheta, const double& upCosTheta)
  {
    fLowCosThetaLimit = lowCosTheta;
    fUpCosThetaLimit = upCosTheta;
    std::cout << this << "  " << fLowVarsLimits[0] << '\n';
  }
  void GEventsGenerator::GSetPhiLimits(const double& lowPhi, const double& upPhi)
  {
    fLowPhiLimit = lowPhi;
    fUpPhiLimit = upPhi;
    std::cout << this << "  " << fLowVarsLimits[0] << '\n';
  }
  void GEventsGenerator::GSetSLimits(const double& lowsij, const double& upsij, const double& lowsjk, const double& upsjk, 
      const double& lowsik, const double& upsik, const double& lows, const double& ups)
  {
    fLowSijLimit = lowsij;
    fUpSijLimit = upsij;
    fLowSjkLimit = lowsjk;
    fUpSjkLimit = upsjk;
    fLowSikLimit = lowsik;
    fUpSikLimit = upsik;
    fLowSLimit = lows;
    fUpSLimit = ups;
  }
  
  void GEventsGenerator::GGenerateEvents(const double& low, const double& up)
  {
    srand(time(0));
    int nEvents = 0, nAllEvents = 0;
    GMinimizeWithMonteCarlo();
    std::cout.precision(16);
    double maxPDFValue = -fMinFunctionValue;
    std::cout << low << "  " << up << '\n';
    while(nEvents < fNEvents)
    {
      auto&& inVars = fMf->GGetPhaseSpace();
      auto&& ui = double(rand())/RAND_MAX;
      auto&& pdf = -fMf->GProcessingComputationOfPDF(&inVars[0], &fParsValues[0]);
      if( (ui*maxPDFValue < pdf) && (std::sqrt(inVars[4]) <= up) && (std::sqrt(inVars[4]) >= 0) )
      {
        nEvents++;
        if(nEvents%1000 == 0 )
                  std::cout << nEvents << std::endl;
        
//         std::cout << "Truth: " << inVars[0] << "  " << inVars[1] << "  " << inVars[2] << "  " << inVars[3] << "  " << inVars[4] << "  " << inVars[5] << '\n';
        auto s = inVars[4];
        auto sKpi1 = inVars[2];
        auto spi1pi2 = inVars[3];
        auto sKpi2 = fMf->GSij(s, sKpi1, spi1pi2);
    
        auto cosTh = ( (fCharge == 0) && (sKpi2>sKpi1) )?(-inVars[0]):inVars[0];
    
        fMf->GPhaseSpaceTo4Vectors();
        (*fKaon_4V).emplace_back(fMf->GGet4VecK());
        (*fPion1_4V).emplace_back(fMf->GGet4VecPi1());
        (*fPion2_4V).emplace_back(fMf->GGet4VecPi2());
        (*fPhoton_4V).emplace_back(fMf->GGet4VecGamma());
        fCosThetaPDF.emplace_back(cosTh);
        fPhiPDF.emplace_back(inVars[1]);
        fSij.emplace_back(inVars[2]);
        fSjk.emplace_back(inVars[3]);
        fM.emplace_back(inVars[4]);
        fPDF1Flat.emplace_back(pdf);
      }
    }
  }
  void GEventsGenerator::GCoeffsForIS()
  {
    fK1coeffOfij = kMK0star_892 * kGammaK0star_892 
    / ( atan2( fUpSijLimit - kMK0star_892 * kMK0star_892, kMK0star_892 * kGammaK0star_892)-
        atan2(fLowSijLimit - kMK0star_892 * kMK0star_892, kMK0star_892 * kGammaK0star_892) );
    fK1coeffOfik = kMK0star_892 * kGammaK0star_892 
    / ( atan2( fUpSikLimit - kMK0star_892 * kMK0star_892, kMK0star_892 * kGammaK0star_892)-
        atan2(fLowSikLimit - kMK0star_892 * kMK0star_892, kMK0star_892 * kGammaK0star_892) );
    fK1coeffOfjk = kMRho0_775 * kGammaRho0_775 
    / ( atan2( fUpSjkLimit - kMRho0_775 * kMRho0_775, kMRho0_775 * kGammaRho0_775)-
        atan2(fLowSjkLimit - kMRho0_775 * kMRho0_775, kMRho0_775 * kGammaRho0_775) );
    fAtancoeffOfij = atan2(fLowSijLimit - kMK0star_892 * kMK0star_892, kMK0star_892 * kGammaK0star_892);
    fAtancoeffOfik = atan2(fLowSikLimit - kMK0star_892 * kMK0star_892, kMK0star_892 * kGammaK0star_892); 
    fAtancoeffOfjk = atan2(fLowSjkLimit - kMRho0_775 * kMRho0_775, kMRho0_775 * kGammaRho0_775); 
    std::cout << fK1coeffOfij << "  " << fK1coeffOfik << "  " << fK1coeffOfjk << "   " << 
    fAtancoeffOfij << "  " << fAtancoeffOfik << "   " << fAtancoeffOfjk << std::endl;
    std::cout << fLowSijLimit << "  " << fUpSijLimit << std::endl;
  }
}

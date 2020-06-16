#include "EvtGenModels/EvtGCouplingConstants.hh"
#include "EvtGenModels/EvtGKinematics.hh"
namespace Gamapola{
  GCouplingConstants::GCouplingConstants(): GInterfaceForMathFunctions(),
  fArg(0)
  {
//     std::cout << "GCouplingConstants constructor calling. . ." << std::endl;
  }
  GCouplingConstants::~GCouplingConstants()
  {
//     std::cout << "GCouplingConstants destructor calling. . ." << std::endl;
  }
  void GCouplingConstants::GCouplingConstantsCalculation()
  {
    double momKst = GKinematics::GMomentaA( kMK0star_892 * kMK0star_892, kMPion1 * kMPion1, kMKaon);
    double momRho = GKinematics::GMomentaA( kMRho0_775 * kMRho0_775, kMPion1 * kMPion1, kMPion2);
    
    fgKstarKpi = sqrt( 3. * 2. * TMath::Pi() * kMK0star_892 * kMK0star_892 * kGammaK0star_892 / ( momKst * momKst * momKst ) );
    fgRhoPiPi = -sqrt( 3. * 2. * TMath::Pi() * kMRho0_775 * kMRho0_775 * kGammaRho0_775 / ( momRho * momRho * momRho ) );
    fgK1KappaPi[0] = 43.0666;
//     sqrt(0.28 * kGammaK1_1270 / GCouplingKappa(sLow, sUp));
    fgK1KappaPi[1] = 0.;
//     std::cout << fgKstarKpi << "   " << fgRhoPiPi << "   " << fgK1KappaPi[0] << std::endl;
  }
}

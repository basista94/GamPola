#include "EvtGenModels/EvtGBreitWigner.hh"
#include "EvtGenModels/EvtGKinematics.hh"

namespace Gamapola{
  GBreitWigner::GBreitWigner():GInterfaceForMathFunctions(),
  fArg(0)
  {
    std::cout << "GBreitWigner constructor calling. . ." << std::endl;
  }
  GBreitWigner::~GBreitWigner()
  {
//     std::cout << "GBreitWigner destructor calling. . ." << std::endl;
  }
  std::complex<double> GBreitWigner::GBWKstr(double sij, int nSwap)
  {
    double widthKstr = GWidthBarrierKstr(sij);
//     double kMPion1 = kMPion;
    double q = GKinematics::GMomentaA(sij,kMKaon*kMKaon, kMPion1);
    double bwf = GBlattWeisskopfFactor(q,0.);
    
    double gamma0 = kMK0star_892 * std::sqrt(kMK0star_892 * kMK0star_892 + kGammaK0star_892 * kGammaK0star_892);
    double c = kMK0star_892 * kGammaK0star_892 * gamma0 / std::sqrt(kMK0star_892 * kMK0star_892 + gamma0);
    
    std::complex<double> denom = std::complex<double>(sij - kMK0star_892 * kMK0star_892, kGammaK0star_892 * kMK0star_892);
    double numer = std::sqrt(c*bwf);
    return 1./denom;
//     std::complex<double> denom(sij - kMK0star_892 * kMK0star_892, kGammaK0star_892 * kMK0star_892);
//     return 1./denom;
  }
  std::complex<double> GBreitWigner::GBWRho(double sij)
  {
    double rel = GKinematics::GMomentaA(sij, kMPion1 * kMPion1, kMPion2) / GKinematics::GMomentaA(kMRho0_775 * kMRho0_775, kMPion1 * kMPion1, kMPion2);
    std::complex<double> denom(sij - kMRho0_775 * kMRho0_775, kGammaRho0_775 * kMRho0_775 * rel * rel * rel);
    return  1./denom;
  }
  std::complex<double> GBreitWigner::GBWK1(double s, int nRes)
  {
//     std::cout << fResMasses[nRes] << "  " << fResWidths[nRes] << "   " << nRes << std::endl;
    std::complex<double> denom(s - fResMasses[nRes] * fResMasses[nRes], fResWidths[nRes] * fResMasses[nRes]);
    return 1./denom;
  }
  std::complex<double> GBreitWigner::GBWKappa(double sij)
  {
    std::complex<double> denom(sij - kMKappa_800 * kMKappa_800, kGammaKappa_800 * kMKappa_800);
    return 1./denom;
  }
  std::complex<double> GBreitWigner::GBWK1bw(const double& s, const double& mass, const double& width)
  {    
    std::complex<double> denom = std::complex<double>(s - mass * mass, width * mass);
    return 1./denom;
  }
  double GBreitWigner::GWidth(const double& s, const double& sbc, int nRes)
  {
//     double kMPion1 = kMPion;
    double q0 = GKinematics::GMomentaA(fResMasses[nRes]*fResMasses[nRes],kMK0star_892*kMK0star_892, kMPion1);
    double q = GKinematics::GMomentaA(s,kMK0star_892*kMK0star_892, kMPion1);
    double bwf = GBlattWeisskopfFactor(q,q0);
    double width = fResWidths[nRes] * (q / q0) * (q / q0) * (q / q0) * (fResMasses[nRes] / sqrt(s)) * bwf;
//     if(q != q)
//       std::cout << s << "  " << (kMK0star_892 + kMPion1)*(kMK0star_892 + kMPion1) << '\n';
//       std::cout << q << "  " << q0 << "  " << (fResMasses[nRes] / sqrt(s)) << bwf << '\n';
    return width;
  }
  double GBreitWigner::GBlattWeisskopfFactor(const double& q, const double& q0)
  {
    double R = 1.5;
    double f2 = (1+R*R*q0*q0)/(1+R*R*q*q);
    return f2;
  }
  double GBreitWigner::GWidthBarrierKstr(const double& sbc)
  {
//     double kMPion1 = kMPion;
    double q0 = GKinematics::GMomentaA(kMK0star_892*kMK0star_892,kMKaon*kMKaon, kMPion1);
    double q = GKinematics::GMomentaA(sbc,kMKaon*kMKaon, kMPion1);
    double bwf = GBlattWeisskopfFactor(q,q0);
    double width = kGammaK0star_892 * (q / q0) * (q / q0) * (q / q0) * (kMK0star_892 / sqrt(sbc)) * bwf;
//     std::cout << q << "  " << width << " "  << kGammaK0star_892 << '\n';
    return width;
  }
}

#ifndef GPHYSICSPDG_HH
#define GPHYSICSPDG_HH

#include <iostream>
#include <cmath>
#include <unordered_map>
typedef const double cd;
namespace Gamapola
{
  template<class T>
  using Dict1D = std::unordered_map<std::string, T>;
  template<class T>
  using Dict2D = std::unordered_map<std::string, Dict1D<T>>;
  class Object
    {
    public:
        virtual ~Object() = 0;
    };
    inline Object::~Object() {}

    template<class T>
    class AnyObj : public Object
    {
    public:
        typedef T Type;
        explicit AnyObj(const Type& data) : data(data) {}
        AnyObj() {}
        Type data;
    };

  cd pi = 4*std::atan(1.);
  cd kPhaseSpaceFactor = 1./(32.* pi * pi * sqrt(pi));
  cd ka = 10000;
  cd kMBmeson = 5.27925999999999984169107847264968;
//   cd kMKaon = 0.49367700;
//   cd kMPion1 = 0.13957018;
//   cd kMPion2 = 0.13957018;//0.134976;
  
//   static const int fMaxNRes = 7;
//   static const int fMaxNDecays = 3;
//   static const int fNSpacialDims = 3;
//   static const int fMaxNCharges = 3;
//   static const int fMaxNWaves = 3;
//   cd pi = 4*atan(1.);
//   cd ka = 10000;
//   cd kMBmeson = 5.27925999999999984169107847264968;
//   cd kMKaon = 0.498;
//   cd kMPion = 0.140;
  cd kMK1_1270 = 1.272;
  cd kGammaK1_1270 = 0.090;
  cd kBr1K1_1270 = 0.16;
  cd kBr2K1_1270 = 0.42;
  
  cd kMK1_1400 = 1.403;
  cd kGammaK1_1400 = 0.174;
  cd kBr1K1_1400 = 0.94;
  cd kBr2K1_1400 = 0.03;
  
  cd kMKst_1410 = 1.414;
  cd kGammaKst_1410 = 0.232;
  cd kFF1Kst_1410 = 1;
  cd kFF2Kst_1410 = 1;
  cd kBr1Kst_1410 = 0.93;
  cd kBr2Kst_1410 = 0.07;
  
  cd kMKst_1680 = 1.717;
  cd kGammaKst_1680 = 0.32;
  cd kFF1Kst_1680 = 1.;
  cd kFF2Kst_1680 = 1.;
  cd kBr1Kst_1680 = 0.299;
  cd kBr2Kst_1680 = 0.314;
  
  cd kMK2_1770 = 1.773;
  cd kGammaK2_1770 = 0.186;
  cd kFF1K2_1770 = 1.;
  cd kFF2K2_1770 = 1.;
  
  cd kMK2_1600 = 1.605;
  cd kGammaK2_1600 = 0.115;
//   cd kMK2_1770 = 1.773;
//   cd kGammaK2_1770 = 0.186;
  
  cd kMK2_1430 = 1.4256;
  cd kGammaK2_1430 = 0.0985;
  cd kFF1K2_1430 = 1.;
  cd kFF2K2_1430 = 1.;
  cd kBr1K2_1430 = 0.247;
  cd kBr2K2_1430 = 0.087;
  
//   Kaons resonances
  cd kMK0star_892 = 0.89581;
  cd kGammaK0star_892 = 0.0474;
  cd kMK0star_1430 = 1.43;
  cd kGammaK0star_1430 = 0.27;
  
//   Kappa resonance
  cd kMKappa_800 = 0.8;
  cd kGammaKappa_800 = 0.5;
  
//   Rho resonances
  cd kMRho0_775 = 0.775;
  cd kGammaRho0_775 = 0.149;
  cd kMRhoPlus = 0.775;
  cd kGammaRhoPlus = 0.150;
  cd kMRhoMinus = 0.775;
  cd kGammaRhoMinus = 0.150;
  
//   
  cd kQPC = 4;
  cd kRpi = 1/0.4;
  cd kRK = 1/0.4;
  cd kRKst = 1/0.4;
  cd kRrho = 1/0.4;
  cd kRomega = 1/0.4;
  cd kRK1 = 1/0.4;
  cd kRb1 = 1/0.4;
}
#endif

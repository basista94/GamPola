#ifndef GEVENTSGENERATOR_H
#define GEVENTSGENERATOR_H

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <limits>
#include <cmath>
#include <TF1.h>
#include <TMath.h>
#include "ginac/ginac.h"
#include "GInterfaceForMinimization.h"

using namespace GiNaC;
namespace Gamapola
{
  class GEventsGenerator: virtual public GInterfaceForMinimization{
  public:
    GEventsGenerator();
    ~GEventsGenerator();
//     void GSetEventsNumber(const int& nEvents);
    void GGenerateEvents(const double& low = std::numeric_limits<float>::min(), const double& up = std::numeric_limits<float>::max());
    void GSetThetaLimits(const double& lowTheta, const double& upTheta);
    void GSetPhiLimits(const double& lowPhi, const double& upPhi);
    void GSetSLimits(const double& lowsij, const double& upsij, const double& lowsjk, const double& upsjk, 
      const double& lowsik, const double& upsik, const double& lows, const double& ups);
//     void GSetDecayMode(const char* decMode[], const int* charges, const int& nRes);
    void GSetCharge(const int& charge);
    void GCoeffsForIS();
    std::vector<double> GetSliceWithinBounds(const std::vector<double>& toCutVec, const std::vector<double>& controlVec,
                                             const double& lowBound, const double& upBound) const;
    const std::vector<double>& GGetCosThetaPDF() const
    {
      return fCosThetaPDF;
    }
    const std::vector<double>& GGetPhiPDF() const
    {
      return fPhiPDF;
    }
    const std::vector<double>& GGetPhiPDF2() const
    {
      return fPhiPDF2;
    }
    const std::vector<double>& GGetPDF1() const
    {
      return fPDF1;
    }
    const std::vector<double>& GGetOmega() const
    {
      return fOmega;
    }
    const std::vector<double>& GGetOmega2() const
    {
      return fOmega2;
    }
    
    const double& GGetIntegralOfPDF() const
    {
      return fIntegralOfPDF;
    }
    const int& GGetEventsNumber() const
    {
      return fNEvents;
    }
    const std::vector<double>& GGetCosThetaPDFFlat() const
    {
      return fCosThetaPDFFlat;
    }
    const std::vector<double>& GGetPhiPDFFlat() const
    {
      return fPhiPDFFlat;
    }
    const std::vector<double>& GGetPDF1Flat() const
    {
      return fPDF1Flat;
    }
    const std::vector<double>& GGetOmegaFlat() const
    {
      return fOmegaFlat;
    }
    const std::vector<double>& GGetOmega2Flat() const
    {
      return fOmega2Flat;
    }
    const std::vector<double>& GGetSijFlat() const
    {
      return fSijFlat;
    }
    const std::vector<double>& GGetSjkFlat() const
    {
      return fSjkFlat;
    }
    const std::vector<double>& GGetSikFlat() const
    {
      return fSikFlat;
    }
    const std::vector<double>& GGetSij() const
    {
      return fSij;
    }
    const std::vector<double>& GGetSjk() const
    {
      return fSjk;
    }
    const std::vector<double>& GGetSik() const
    {
      return fSik;
    }
    const std::vector<double>& GGetMFlat() const
    {
      return fMFlat;
    }
    const std::vector<double>& GGetM() const
    {
      return fM;
    }
    const std::complex<double>* GGetNormalizationIntegrals()
    {
      return fNormalizationIntegrals;
    }
    const std::vector<GVector4D>& GGetKaonV() const
    {
        return *fKaon_4V;
    }
    const std::vector<GVector4D>& GGetPion1V() const
    {
        return *fPion1_4V;
    }
    const std::vector<GVector4D>& GGetPion2V() const
    {
        return *fPion2_4V;
    }
    const std::vector<GVector4D>& GGetPhotonV() const
    {
        return *fPhoton_4V;
    }
  private:
    int fArg;
//     limits on generation of statistical parameters
    double fLowCosThetaLimit;
    double fUpCosThetaLimit;
    double fLowPhiLimit;
    double fUpPhiLimit;
    double fLowSijLimit;
    double fUpSijLimit;
    double fLowSjkLimit;
    double fUpSjkLimit;
    double fLowSikLimit;
    double fUpSikLimit;
    double fLowSLimit;
    double fUpSLimit;
//     generated kinematic variables
    std::vector<double> fCosThetaPDF;
    std::vector<double> fPhiPDF;
    std::vector<double> fPhiPDF2;
    std::vector<double> fPDF1;
    std::vector<double> fOmega;
    std::vector<double> fOmega2;
//     generated flat kinematic variables
    std::vector<double> fCosThetaPDFFlat;
    std::vector<double> fPhiPDFFlat;
    std::vector<double> fPDF1Flat;
    std::vector<double> fOmegaFlat;
    std::vector<double> fOmega2Flat;
    std::vector<double> fSijFlat;
    std::vector<double> fSjkFlat;
    std::vector<double> fSikFlat;
    std::vector<double> fSij;
    std::vector<double> fSjk;
    std::vector<double> fSik;
    std::vector<double> fMFlat;
    std::vector<double> fM;
    //     Coefficients for importance sampling
    double fK1coeffOfij;
    double fAtancoeffOfij;
    double fK1coeffOfik;
    double fAtancoeffOfik;
    double fK1coeffOfjk;
    double fAtancoeffOfjk;
    
    std::shared_ptr<std::vector<GVector4D>> fKaon_4V;
    std::shared_ptr<std::vector<GVector4D>> fPion1_4V;
    std::shared_ptr<std::vector<GVector4D>> fPion2_4V;
    std::shared_ptr<std::vector<GVector4D>> fPhoton_4V;
  };
}
#endif

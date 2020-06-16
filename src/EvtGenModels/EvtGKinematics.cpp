#include "EvtGenModels/EvtGKinematics.hh"
#include "EvtGenModels/EvtGSimpleRandomEngine.hh"
#include "EvtGenModels/EvtGPhysicsPDG.hh"

namespace Gamapola{
  GKinematics::GKinematics():GInterfaceForMathFunctions(),
  fArg(0)
  {
    std::cout << "GKinematics constructor calling. . ." << std::endl;
//     if(!fRandomEngine)
    fRandomEngine = new GSimpleRandomEngine();
  }
  GKinematics::~GKinematics()
  {
      if(fRandomEngine)
          delete fRandomEngine;
//     std::cout << "GKinematics destructor calling. . ." << std::endl;
  }
  void GKinematics::GKinematicsConstanstsCalculation(const double& mRes, const int& nRes)
  {
      
    fMomentaConst[fIndex[nRes]][0] = GMomentaA(mRes * mRes, kMK0star_892*kMK0star_892,kMPion1);
    fMomentaConst[fIndex[nRes]][1] = GMomentaA(mRes * mRes, kMRho0_775*kMRho0_775,kMKaon);
    fMomentaConst[fIndex[nRes]][2] = GMomentaA(mRes * mRes, kMK0star_1430*kMK0star_1430,kMPion1);
  }
  double GKinematics::GSjkMin(const double& sij, const double& mi, const double& mj, const double& mk, const double& M)
  {
    auto&& Ejstar = GEnergyA(sij, mj*mj, mi);
//     double Ejstar = ( sij - mi * mi + mj * mj ) / ( 2 * sqrt ( sij ) );
    double Ekstar = ( M * M - sij - mk * mk ) / ( 2 * sqrt ( sij ) );
    double sjk = ( Ejstar + Ekstar ) * ( Ejstar + Ekstar ) - 
                 ( sqrt( Ejstar * Ejstar - mj * mj ) + sqrt( Ekstar * Ekstar - mk * mk ) ) * 
                 ( sqrt( Ejstar * Ejstar - mj * mj ) + sqrt( Ekstar * Ekstar - mk * mk ) );
    return  sjk ;
  }
  double GKinematics::GSjkMax(const double& sij, const double& mi, const double& mj, const double& mk, const double& M)
  {
    auto&& Ejstar = GEnergyA(sij, mj*mj, mi);
//     double Ejstar = ( sij - mi * mi + mj * mj ) / ( 2 * sqrt ( sij ) );
    double Ekstar = ( M * M - sij - mk * mk ) / ( 2 * sqrt ( sij ) );
    double sjk = ( Ejstar + Ekstar ) * ( Ejstar + Ekstar ) - 
                 ( sqrt( Ejstar * Ejstar - mj * mj ) - sqrt( Ekstar * Ekstar - mk * mk ) ) * 
                 ( sqrt( Ejstar * Ejstar - mj * mj ) - sqrt( Ekstar * Ekstar - mk * mk ) );
    return  sjk ;
  } 
  double GKinematics::GMomentaA(double s, double sij, double mb)
  { 
    return sqrt ( ( s - ( sqrt( sij ) + mb) * (sqrt( sij ) + mb) ) * 
    (s - (sqrt( sij ) - mb) * (sqrt( sij ) - mb) ) / ( 4. * s ) );
  }
  double GKinematics::GEnergyA(double s, double sij, double mb)
  {
//     std::cout << sqrt(s) << "  " << mb << "  " << sqrt(sij) << std::endl;
    return ( s - mb * mb + sij ) / ( 2 * sqrt( s ) );
  } 
  double GKinematics::GPi1Pi2Vec4(double s, double sij, double sjk)
  {
//     std::cout << s << " " << sij << "  " << sjk << "  " <<  kMKaon << std::endl;
    return ( s - sij - sjk + kMKaon * kMKaon ) / 2. ;
  }
  void GKinematics::GPi1Pi2Vec3()
  {
    fPiPjVec[0] = fMomiV[1] * fMomjV[2] - fMomiV[2] * fMomjV[1];
    fPiPjVec[1] = fMomiV[2] * fMomjV[0] - fMomiV[0] * fMomjV[2];
    fPiPjVec[2] = fMomiV[0] * fMomjV[1] - fMomiV[1] * fMomjV[0];
//     std::cout << fPiPjVec[0] << "  " << fMomiV[1] * fMomjV[2] << "  " << fMomiV[2] * fMomjV[1] << std::endl;
  }
  double GKinematics::GSij(double s, double sij, double sjk)
  {
    return ( s - sij - sjk + kMKaon * kMKaon + kMPion1 * kMPion1 + kMPion2 * kMPion2) ;
  }
  std::complex<double> GKinematics::GScPrd(const double* pi, const std::complex<double>* pj)
  {
    return pi[0] * pj[0] + pi[1] * pj[1] + pi[2] * pj[2]; 
  }
  double GKinematics::GPiPj()
  {
    return fEnergy1 * fEnergy2 + (kMPion1 * kMPion1 + kMPion2 * kMPion2 - f_Sij) / 2.;
  }
  void GKinematics::GCalculateKinematics(double* x, const int& charge)
  {
    fCosTheta = x[0];
    auto s = x[4];
    auto sKpi1 = x[2];
    auto spi1pi2 = x[3];
    auto sKpi2 = GSij(s, sKpi1, spi1pi2);
    
    double cosTh = ( (charge == 0) && (sKpi2>sKpi1) )?(-x[0]):x[0];
    double sinTh = sin(acos(cosTh));

    fEnergy1 = GEnergyA(s, kMPion1*kMPion1, sqrt(sKpi2));
    fEnergy2 = GEnergyA(s, kMPion2*kMPion2, sqrt(sKpi1));
    fEnergy3 = GEnergyA(s, kMKaon*kMKaon, sqrt(spi1pi2));
    fEnergy4 = (kMBmeson*kMBmeson-s)/(2*sqrt(s));

    fPiPjVec4 = GPi1Pi2Vec4(s, sKpi1, sKpi2);
    f_Sij = GSij(s, sKpi1, sKpi2);
    fMomij = GPiPj();
    
    fMomi = sqrt(fEnergy1 * fEnergy1 - kMPion1 * kMPion1);
    fMomj = sqrt(fEnergy2 * fEnergy2 - kMPion2 * kMPion2);
    fMomk = sqrt(fEnergy3 * fEnergy3 - kMKaon * kMKaon);
    fMoml = fEnergy4; 

    double delta = acos(fMomij / (fMomi * fMomj));
    auto&& alpha = atan2(fMomj*sin(delta),fMomi+fMomj*cos(delta));
    
    fPhii = (2 * x[1] - delta) / 2;
    fPhij = (2 * x[1] + delta) / 2;
    fPhik = pi+alpha+x[1]-delta/2.;

    fMomiV[0] = fMomi * cosTh * cos(fPhii);
    fMomiV[1] = fMomi * sin(fPhii);
    fMomiV[2] = -fMomi * sinTh * cos(fPhii);
    
    fMomjV[0] = fMomj * cosTh * cos(fPhij);
    fMomjV[1] = fMomj * sin(fPhij);
    fMomjV[2] = -fMomj * sinTh * cos(fPhij);
    
    fMomkV[0] = fMomk * cosTh * cos(fPhik);
    fMomkV[1] = fMomk * sin(fPhik);
    fMomkV[2] = -fMomk * sinTh * cos(fPhik);
    
    fqK.set(fEnergy3, fMomk * cosTh * cos(fPhik), fMomk * sin(fPhik), -fMomk * sinTh * cos(fPhik));
    fqPi1.set(fEnergy1, fMomi * cosTh * cos(fPhii), fMomi * sin(fPhii), -fMomi * sinTh * cos(fPhii));
    fqPi2.set(fEnergy2, fMomj * cosTh * cos(fPhij), fMomj * sin(fPhij), -fMomj * sinTh * cos(fPhij));
    
    GPi1Pi2Vec3();
    fE0PiPj = GScPrd(fPiPjVec, fEpsilon0);
    fE0Pi = GScPrd(fMomiV, fEpsilon0);
    fE0Pj = GScPrd(fMomjV, fEpsilon0);
    for(size_t iAxis = 0; iAxis < 3; iAxis++)
    {
      fK2KinC4combo[iAxis] = ( fE0Pi * fPiPjVec[iAxis] + fE0PiPj * fMomiV[iAxis] );
      fK2KinC5combo[iAxis] = ( fE0Pj * fPiPjVec[iAxis] + fE0PiPj * fMomjV[iAxis] );
//       std::cout << "Check kin: " << fK2KinC4combo[iAxis] << "   " << fK2KinC5combo[iAxis] << '\n';
//       if(fK2KinC5combo[iAxis] != fK2KinC5combo[iAxis])
//           std::cout << fE0Pj << "  " << fPiPjVec[iAxis] << "  " << fE0PiPj << "  " << fMomjV[iAxis] << '\n';
    }
  }
  
  std::vector<double> GKinematics::GGetPhaseSpace()
  {
      auto&& M = 1.4;
      auto&& s12Min = (kMKaon + kMPion1) * (kMKaon + kMPion1);
      auto&& s23Min = (kMPion1 + kMPion2) * (kMPion1 + kMPion2);
      auto&& s13Min = (kMKaon + kMPion2) * (kMKaon + kMPion2);
      auto&& sMin = 1./*(m1 + m2 + m3) * (m1 + m2 + m3)*/;
      auto&& sMax = (M+0.6)*(M+0.6);
      while(1)
      {
        auto&& s = sMin + ( sMax - sMin ) * double(rand())/RAND_MAX;
        auto&& M = sqrt(s);
        auto&& s12Max = ( M - kMPion2 ) * ( M - kMPion2 );
        auto&& s23Max = ( M - kMKaon ) * ( M - kMKaon );
        auto&& sjk = s23Min + (s23Max - s23Min)*double(rand())/RAND_MAX;
        auto&& sij = s12Min + (s12Max - s12Min)*double(rand())/RAND_MAX;
        auto&& lowsjk = GSjkMin(sij, kMKaon, kMPion1, kMPion2, M);
        auto&& upsjk = GSjkMax(sij, kMKaon, kMPion1, kMPion2, M);
        if((sjk < lowsjk) || (sjk > upsjk))
            continue;
//         auto&& sik = kMKaon*kMKaon+kMPion1*kMPion1+kMPion2*kMPion2+M*M-sjk-sij;
        auto&& costheta = -1 + 2 * double(rand())/RAND_MAX;
        auto&& phi = 0. + (2 * pi - 0.) * double(rand())/RAND_MAX;
        
//         std::cout << M << '\n';
//         std::cout << costheta << "  " << phi << "  " << sij << "  " << sjk << "   " << M*M << '\n';
        return std::vector<double>{costheta, phi, sij, sjk,  M*M};
      }
  }
  void GKinematics::GPhaseSpaceTo4Vectors()
  {
      auto&& cosTh = fCosTheta;
      auto&& sinTh = sin(acos(cosTh));
    
      GVector4D q_K(fEnergy3, 
                    fMomk * cosTh * cos(fPhik), 
                    fMomk * sin(fPhik), 
                    -fMomk * sinTh * cos(fPhik));
      
      GVector4D q_piplus(fEnergy1, 
                         fMomi * cosTh * cos(fPhii), 
                         fMomi * sin(fPhii), 
                         -fMomi * sinTh * cos(fPhii));
      
      GVector4D q_pizero(fEnergy2, 
                         fMomj * cosTh * cos(fPhij), 
                         fMomj * sin(fPhij), 
                         -fMomj * sinTh * cos(fPhij));
      
      GVector4D q_gamma(fEnergy4,0,0,-fMoml); 
      
      
      auto twoPi = 2 * pi;
      auto mB2 = kMBmeson*kMBmeson;
      auto M2 = (fEnergy1+fEnergy2+fEnergy3) * (fEnergy1+fEnergy2+fEnergy3);
      auto coeff = (mB2-M2)/(mB2+M2); 
//     
//     std::cout << q_K+q_piplus+q_pizero << '\n'; 
    
//       std::cout << q_K << "  " << q_piplus << "  " << q_pizero << "  " << q_gamma << '\n';
      auto alphaAngle = twoPi * fRandomEngine->random();
      auto beta = std::acos(-1 + 2*fRandomEngine->random());
      GVector3D boost(0,0,coeff);
      
      q_gamma = boostTo(q_gamma,boost,false);
      q_K = boostTo(q_K,boost,false);
      q_piplus = boostTo(q_piplus,boost,false);
      q_pizero = boostTo(q_pizero,boost,false);
      
      
      
      q_gamma.applyRotateEuler(alphaAngle,beta,-alphaAngle);
      q_K.applyRotateEuler(alphaAngle,beta,-alphaAngle);
      q_piplus.applyRotateEuler(alphaAngle,beta,-alphaAngle);
      q_pizero.applyRotateEuler(alphaAngle,beta,-alphaAngle);
      
      
      fq_gamma = q_gamma;
      fq_K = q_K;
      fq_pi1 = q_piplus;
      fq_pi2 = q_pizero;
  }
  
  std::vector<double> GKinematics::G4VectorsToPhaseSpace(const GVector4D& qK, 
                                                         const GVector4D& qpi1, 
                                                         const GVector4D& qpi2, 
                                                         const GVector4D& qGamma, 
                                                         const int& charge)
  {
      auto&& k1new = qK+qpi1+qpi2;
      auto&& mK = qK.mass();
      auto&& mPi1 = qpi1.mass();
      auto&& mPi2 = qpi2.mass();
      auto&& sKpi1 = (qK+qpi1)*(qK+qpi1);
      auto&& spi1pi2 = (qpi1+qpi2)*(qpi1+qpi2);
      auto&& M = k1new.mass();
//       std::cout << sKpi1 << "  " << spi1pi2 << "  " << M*M << '\n';
      auto&& sKpi2 = mPi1*mPi1+mPi2*mPi2+mK*mK+M*M-spi1pi2-sKpi1;
      GVector3D boost2Kres(-k1new.get(1)/k1new.get(0),-k1new.get(2)/k1new.get(0),-k1new.get(3)/k1new.get(0));
      auto&& qK_Kres = boostTo(qK, boost2Kres, false);
      auto&& qGamma_Kres = boostTo(qGamma,boost2Kres,false);
      auto&& qpi1_Kres = boostTo(qpi1,boost2Kres,false);
      auto&& qpi2_Kres = boostTo(qpi2,boost2Kres,false);
      GVector4D ez(1,-qGamma_Kres.get(1),-qGamma_Kres.get(2),-qGamma_Kres.get(3));
      ez = ez/qGamma_Kres.d3mag();
      auto&& n = qpi1_Kres.cross(qpi2_Kres)/qpi1_Kres.cross(qpi2_Kres).d3mag();
      auto&& costheta = n.dot(ez);
      auto&& ey = ez.cross(n)/ez.cross(n).d3mag();
      auto&& ex = ey.cross(ez);
      auto&& theta = acos(costheta);
      auto&& exprime = ex*cos(theta)-ez*sin(theta);
      auto&& phi1 = atan2(ey.dot(qpi1_Kres),exprime.dot(qpi1_Kres));
      auto&& phi2 = atan2(ey.dot(qpi2_Kres),exprime.dot(qpi2_Kres));      
      auto&& phi = (phi1+phi2)/2.;
      auto&& delta = acos(qpi1_Kres.dot(qpi2_Kres)/(qpi1_Kres.d3mag()*qpi2_Kres.d3mag()));

      if(fabs(delta-(phi2 - phi1)) <= 1e-3)
          phi = (phi < 0)?(phi+2*pi):phi;
      else
          phi+=pi;
      
      return std::vector<double>{costheta, phi, sKpi1, spi1pi2, M * M};
  }
}

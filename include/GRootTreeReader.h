#ifndef GROOTTREEREADER_H
#define GROOTTREEREADER_H

#include "TTree.h"
#include "TFile.h"
#include "EvtGenModels/EvtGVector4D.hh"
#include <memory>

namespace Gamapola
{
    class GRootTreeReader
    {
    public:
        explicit GRootTreeReader(const std::string& dataFile, const int& charge)
        {
            double pi_const = 4*std::atan(1);
            hfile2 = std::make_shared<TFile>(dataFile.c_str());
            if(hfile2)
                std::cout << "HEMA" << '\n';
            auto&& raw_tree = (TTree*)hfile2->Get("DalitzEventList");
            T2.reset(raw_tree);
            if(charge == 1)
            {
                T2->SetBranchAddress("_1_K#_E",&(kaon4Vec[0]));
                T2->SetBranchAddress("_1_K#_Px",&(kaon4Vec[1]));
                T2->SetBranchAddress("_1_K#_Py",&(kaon4Vec[2]));
                T2->SetBranchAddress("_1_K#_Pz",&(kaon4Vec[3]));
  
                T2->SetBranchAddress("_2_pi#_E",&(pi14Vec[0]));
                T2->SetBranchAddress("_2_pi#_Px",&(pi14Vec[1]));
                T2->SetBranchAddress("_2_pi#_Py",&(pi14Vec[2]));
                T2->SetBranchAddress("_2_pi#_Pz",&(pi14Vec[3]));
                
                T2->SetBranchAddress("_3_pi~_E",&(pi24Vec[0]));
                T2->SetBranchAddress("_3_pi~_Px",&(pi24Vec[1]));
                T2->SetBranchAddress("_3_pi~_Py",&(pi24Vec[2]));
                T2->SetBranchAddress("_3_pi~_Pz",&(pi24Vec[3]));
  
                T2->SetBranchAddress("_4_gamma0_E",&(gamma4Vec[0]));
                T2->SetBranchAddress("_4_gamma0_Px",&(gamma4Vec[1]));
                T2->SetBranchAddress("_4_gamma0_Py",&(gamma4Vec[2]));
                T2->SetBranchAddress("_4_gamma0_Pz",&(gamma4Vec[3]));
            }   
            else
            {
                T2->SetBranchAddress("_1_K#_E",&(kaon4Vec[0]));
                T2->SetBranchAddress("_1_K#_Px",&(kaon4Vec[1]));
                T2->SetBranchAddress("_1_K#_Py",&(kaon4Vec[2]));
                T2->SetBranchAddress("_1_K#_Pz",&(kaon4Vec[3]));
                
                T2->SetBranchAddress("_2_pi0_E",&(pi14Vec[0]));
                T2->SetBranchAddress("_2_pi0_Px",&(pi14Vec[1]));
                T2->SetBranchAddress("_2_pi0_Py",&(pi14Vec[2]));
                T2->SetBranchAddress("_2_pi0_Pz",&(pi14Vec[3]));
                
                T2->SetBranchAddress("_3_pi~_E",&(pi24Vec[0]));
                T2->SetBranchAddress("_3_pi~_Px",&(pi24Vec[1]));
                T2->SetBranchAddress("_3_pi~_Py",&(pi24Vec[2]));
                T2->SetBranchAddress("_3_pi~_Pz",&(pi24Vec[3]));
  
                T2->SetBranchAddress("_4_gamma0_E",&(gamma4Vec[0]));
                T2->SetBranchAddress("_4_gamma0_Px",&(gamma4Vec[1]));
                T2->SetBranchAddress("_4_gamma0_Py",&(gamma4Vec[2]));
                T2->SetBranchAddress("_4_gamma0_Pz",&(gamma4Vec[3]));
            }
            
            T2->SetBranchAddress("cosTheta",&(fPhaseSpace[0]));
            T2->SetBranchAddress("phi",&(fPhaseSpace[1]));
            T2->SetBranchAddress("sKPi1",&(fPhaseSpace[2]));
            T2->SetBranchAddress("sPi1Pi2",&(fPhaseSpace[3]));
            T2->SetBranchAddress("Mkpipi",&(fPhaseSpace[4]));
            T2->SetBranchAddress("decrate",&(decrate));
//             std::cout << T2->GetEntries() << '\n';
        }
        
        void GGetEntry(const int& entry)
        {
            T2->GetEntry(entry);
            GVector4D q_K(kaon4Vec[0], kaon4Vec[1], kaon4Vec[2], kaon4Vec[3]);
            GVector4D q_pi1(pi14Vec[0], pi14Vec[1], pi14Vec[2], pi14Vec[3]);
            GVector4D q_pi2(pi24Vec[0], pi24Vec[1], pi24Vec[2], pi24Vec[3]);
            GVector4D q_gamma(gamma4Vec[0], gamma4Vec[1], gamma4Vec[2], gamma4Vec[3]);
            fKaon = q_K;
            fPion1 = q_pi1;
            fPion2 = q_pi2;
            fGamma = q_gamma;
        }
        
        const GVector4D& GGetKaonFourVec() const { return fKaon; }
        const GVector4D& GGetPion1FourVec() const { return fPion1; }
        const GVector4D& GGetPion2FourVec() const { return fPion2; }
        const GVector4D& GGetGammaFourVec() const { return fGamma; }
        const std::array<double, 5>& GGetPhaseSpace() const {return fPhaseSpace;}
        const double& GGetDecRate() const { return decrate; }
        int GGetNEntries() const { return T2->GetEntries(); }
    private:
        const static int vecLength = 4;
        double kaon4Vec[vecLength] = {0};
        double pi14Vec[vecLength] = {0};
        double pi24Vec[vecLength] = {0};
        double gamma4Vec[vecLength] = {0};
        double decrate;
        std::array<double, 5> fPhaseSpace;
        GVector4D fKaon;
        GVector4D fPion1;
        GVector4D fPion2;
        GVector4D fGamma;
        std::shared_ptr<TFile> hfile2;
        std::shared_ptr<TTree> T2;
    };
}

#endif

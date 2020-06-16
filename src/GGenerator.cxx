#include "GGenerator.h"
#include "EvtGenModels/EvtGPhysicsPDG.hh"
#include "GamPolaParser/GGeneratorParser.h"
#include "GamPolaPlotter/GGampolaPlotter.h"
#include "GamPolaParser/GCutsParser.h"
#include "TTree.h"
#include "TFile.h"

namespace Gamapola
{
    GGenerator::GGenerator(const int& charge):
    mf(new GMathFunctions),
    eg(new GEventsGenerator),
    modelPars(0),
//     fDictionary(0),
    fCharge(charge)
    {
    }
    
    GGenerator::~GGenerator()
    {
    }
    void GGenerator::GSetCouplings(const std::vector<double>& couplings)
    {
        eg->GSetFunctionParameters(couplings.size(), const_cast<double*>(&couplings[0]));
        eg->GSetFunction(&Gamapola::GMathFunctions::GProcessingComputationOfPDF);
        eg->GSetDecayMode(fCharge);
    }
    void GGenerator::GSetCouplings(const std::string& fileName)
    {
        std::shared_ptr<GGeneratorParser> gp(new GGeneratorParser(fileName));
        auto&& model = gp->GParse();
     eg->GSetFunctionParameters(model.size(), &model[0]);
     eg->GSetFunction(&Gamapola::GMathFunctions::GProcessingComputationOfPDF);
     eg->GSetDecayMode(fCharge);
    }
    
    void GGenerator::GGenerate(const int& nEvents, const std::string& cut)
    {
        std::shared_ptr<GCutsParser> fp(new GCutsParser(cut));
        auto dictionary = fp->GParse();
        auto&& low = dictionary["sKpipi"].first;
        auto&& up = dictionary["sKpipi"].second;
        if(up == 0)
        {
//             low = std::numeric_limits<float>::min();
            up = std::numeric_limits<float>::max();
        }
        eg->GSetEventsNumber(nEvents);
        eg->GGenerateEvents(low, up);
    }
    void GGenerator::GWriteToContainer()
    {
        auto&& nEvents = static_cast<int>(eg->GGetCosThetaPDF().size());
        for(auto iEv = 0; iEv < nEvents; ++iEv)
        {
            fEvents[0].emplace_back(eg->GGetCosThetaPDF()[iEv]);
            fEvents[1].emplace_back(eg->GGetPhiPDF()[iEv]);
            fEvents[2].emplace_back(eg->GGetSij()[iEv]);
            fEvents[3].emplace_back(eg->GGetSjk()[iEv]);
            fEvents[4].emplace_back(eg->GGetM()[iEv]);
        }
    }
    void GGenerator::GWriteToFile(const std::string& outFile)
    {
        double Bmeson[4];
        double Kaon[4];
        double pi1[4];
        double pi2[4];
        double gamma[4];
   
        int pdgB = 0;
        int pdgK = 0;
        int pdgpi1 = 0;
        int pdgpi2 = 0;
        int pdggamma = 0;
        double decrate = 0.;
        double mkpipi, costheta, phi, sKPi1, sPi1Pi2;
        
        auto&& nEvents = static_cast<int>(eg->GGetKaonV().size());
        TFile f(outFile.c_str(),"recreate");
        TTree t2("DalitzEventList","DalitzEventList");
        if(fCharge == 0)
        {
            pdgB = 511;
            pdgK = 321;
            pdgpi1 = 111;
            pdgpi2 = -211;
            pdggamma = 22;
            
            t2.Branch("_0_B0_pdg",&pdgB,"_0_B0_pdg/I");
            t2.Branch("_0_B0_E",&Bmeson[0],"_0_B0_E/D");
            t2.Branch("_0_B0_Px",&Bmeson[1],"_0_B0_Px/D");
            t2.Branch("_0_B0_Py",&Bmeson[2],"_0_B0_Py/D");
            t2.Branch("_0_B0_Pz",&Bmeson[3],"_0_B0_Pz/D");
   
            t2.Branch("_1_K#_pdg",&pdgK,"_1_K#_pdg/I");
            t2.Branch("_1_K#_E",&Kaon[0],"_1_K#_E/D");
            t2.Branch("_1_K#_Px",&Kaon[1],"_1_K#_Px/D");
            t2.Branch("_1_K#_Py",&Kaon[2],"_1_K#_Py/D");
            t2.Branch("_1_K#_Pz",&Kaon[3],"_1_K#_Pz/D");
   
            t2.Branch("_2_pi0_pdg",&pdgpi1,"_2_pi0_pdg/I");
            t2.Branch("_2_pi0_E",&pi1[0],"_2_pi0_E/D");
            t2.Branch("_2_pi0_Px",&pi1[1],"_2_pi0_Px/D");
            t2.Branch("_2_pi0_Py",&pi1[2],"_2_pi0_Py/D");
            t2.Branch("_2_pi0_Pz",&pi1[3],"_2_pi0_Pz/D");
   
            t2.Branch("_3_pi~_pdg",&pdgpi2,"_3_pi~_pdg/I");
            t2.Branch("_3_pi~_E",&pi2[0],"_3_pi~_E/D");
            t2.Branch("_3_pi~_Px",&pi2[1],"_3_pi~_Px/D");
            t2.Branch("_3_pi~_Py",&pi2[2],"_3_pi~_Py/D");
            t2.Branch("_3_pi~_Pz",&pi2[3],"_3_pi~_Pz/D");
   
            t2.Branch("_4_gamma0_pdg",&pdggamma,"_4_gamma0_pdg/I");
            t2.Branch("_4_gamma0_E",&gamma[0],"_4_gamma0_E/D");
            t2.Branch("_4_gamma0_Px",&gamma[1],"_4_gamma0_Px/D");
            t2.Branch("_4_gamma0_Py",&gamma[2],"_4_gamma0_Py/D");
            t2.Branch("_4_gamma0_Pz",&gamma[3],"_4_gamma0_Pz/D");
            
            t2.Branch("cosTheta", &costheta, "cosTheta/D");
            t2.Branch("phi", &phi, "phi/D");
            t2.Branch("sKPi1", &sKPi1, "sKPi1/D");
            t2.Branch("sPi1Pi2", &sPi1Pi2, "sPi1Pi2/D");
            t2.Branch("Mkpipi", &mkpipi, "Mkpipi/D");
            t2.Branch("decrate", &decrate, "decrate/D");
        }
        else
        {
            pdgB = 521;
            pdgK = 321;
            pdgpi1 = 211;
            pdgpi2 = -211;
            pdggamma = 22;
            
            t2.Branch("_0_B#_pdg",&pdgB,"_0_B#_pdg/I");
            t2.Branch("_0_B#_E",&Bmeson[0],"_0_B#_E/D");
            t2.Branch("_0_B#_Px",&Bmeson[1],"_0_B#_Px/D");
            t2.Branch("_0_B#_Py",&Bmeson[2],"_0_B#_Py/D");
            t2.Branch("_0_B#_Pz",&Bmeson[3],"_0_B#_Pz/D");
   
            t2.Branch("_1_K#_pdg",&pdgK,"_1_K#_pdg/I");
            t2.Branch("_1_K#_E",&Kaon[0],"_1_K#_E/D");
            t2.Branch("_1_K#_Px",&Kaon[1],"_1_K#_Px/D");
            t2.Branch("_1_K#_Py",&Kaon[2],"_1_K#_Py/D");
            t2.Branch("_1_K#_Pz",&Kaon[3],"_1_K#_Pz/D");
   
            t2.Branch("_2_pi#_pdg",&pdgpi1,"_2_pi#_pdg/I");
            t2.Branch("_2_pi#_E",&pi1[0],"_2_pi#_E/D");
            t2.Branch("_2_pi#_Px",&pi1[1],"_2_pi#_Px/D");
            t2.Branch("_2_pi#_Py",&pi1[2],"_2_pi#_Py/D");
            t2.Branch("_2_pi#_Pz",&pi1[3],"_2_pi#_Pz/D");
   
            t2.Branch("_3_pi~_pdg",&pdgpi2,"_3_pi~_pdg/I");
            t2.Branch("_3_pi~_E",&pi2[0],"_3_pi~_E/D");
            t2.Branch("_3_pi~_Px",&pi2[1],"_3_pi~_Px/D");
            t2.Branch("_3_pi~_Py",&pi2[2],"_3_pi~_Py/D");
            t2.Branch("_3_pi~_Pz",&pi2[3],"_3_pi~_Pz/D");
   
            t2.Branch("_4_gamma0_pdg",&pdggamma,"_4_gamma0_pdg/I");
            t2.Branch("_4_gamma0_E",&gamma[0],"_4_gamma0_E/D");
            t2.Branch("_4_gamma0_Px",&gamma[1],"_4_gamma0_Px/D");
            t2.Branch("_4_gamma0_Py",&gamma[2],"_4_gamma0_Py/D");
            t2.Branch("_4_gamma0_Pz",&gamma[3],"_4_gamma0_Pz/D");
            
            t2.Branch("cosTheta", &costheta, "cosTheta/D");
            t2.Branch("phi", &phi, "phi/D");
            t2.Branch("sKPi1", &sKPi1, "sKPi1/D");
            t2.Branch("sPi1Pi2", &sPi1Pi2, "sPi1Pi2/D");
            t2.Branch("Mkpipi", &mkpipi, "Mkpipi/D");
            t2.Branch("decrate", &decrate, "decrate/D");
        }
   
        
        for(auto iEv = 0; iEv < nEvents; ++iEv)
        {
            for(auto iComp = 0; iComp < 4; ++iComp)
            {
                Kaon[iComp] = eg->GGetKaonV()[iEv].get(iComp);
                pi1[iComp] = eg->GGetPion1V()[iEv].get(iComp);
                pi2[iComp] = eg->GGetPion2V()[iEv].get(iComp);
                gamma[iComp] = eg->GGetPhotonV()[iEv].get(iComp);
                Bmeson[iComp] = Kaon[iComp] + pi1[iComp] + pi2[iComp] + gamma[iComp];
            }
//             auto inv = std::sqrt((Kaon[0] + pi1[0] + pi2[0]) * (Kaon[0] + pi1[0] + pi2[0]) - 
//             (Kaon[1] + pi1[1] + pi2[1]) * (Kaon[1] + pi1[1] + pi2[1]) -
//             (Kaon[2] + pi1[2] + pi2[2]) * (Kaon[2] + pi1[2] + pi2[2]) - 
//             (Kaon[3] + pi1[3] + pi2[3]) * (Kaon[3] + pi1[3] + pi2[3]));
//             
//             if(inv > 1.8 || inv < 1)
//                 std::cout << inv << "  " << iEv << '\n';
            costheta = eg->GGetCosThetaPDF()[iEv];
            phi = eg->GGetPhiPDF()[iEv];
            sKPi1 = eg->GGetSij()[iEv];
            sPi1Pi2 = eg->GGetSjk()[iEv];
            mkpipi = std::sqrt(eg->GGetM()[iEv]);
            decrate = eg->GGetPDF1Flat()[iEv];
            t2.Fill();
        }
//         P(sick|been in Italy) = P(been in Italy|sick)*P(A)/P(B)
        
        t2.Write();
            
    }

}

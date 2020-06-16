#include "GFitter.h"
#include "GamPolaParser/GFitterParser.h"
#include "GamPolaParser/GCutsParser.h"
#include "GGOFKolmogorovTest.h"
#include "TTree.h"
#include "TFile.h"

using namespace GiNaC;

namespace Gamapola{
    GFitter::GFitter(const int& charge, const std::string& tosave):
    ea(new GEventsAnalyzer(charge)),
    mf(new GMathFunctions()),
    finVars(0),
    fkindOfVars(0),
    fnamesOfVars(0),
    flowVars(0),
    fupVars(0),
    fstepsVars(0),
    pars(0),
    fDictionary(0),
    fCharge(charge)
    {
        srand(time(0));
        fDictionary[0] = "gammaQPC"; fDictionary[1] = "f2"; fDictionary[2] = "thetaK1";
        fDictionary[3] = "phiDK*"; fDictionary[4] = "phiSRho"; fDictionary[5] = "phiDRho";
        fDictionary[6] = "ggKappaRe"; fDictionary[7] = "lambda"; fDictionary[8] = "ffRe_1400";
        fDictionary[9] = "ffIm_1400"; fDictionary[10] = "ggKappaIm"; fDictionary[11] = "ff2_Re_1410";
        fDictionary[12] = "ff2_Im_1410"; fDictionary[13] = "ff3_Re_1680"; fDictionary[14] = "ff3_Im_1680";
        fDictionary[15] = "ff4_Re_1430"; fDictionary[16] = "ff4_Im_1430"; fDictionary[17] = "gg3Kstr_1410";
        fDictionary[18] = "gg3ImKstr_1410"; fDictionary[19] = "gg3Rho_1410"; fDictionary[20] = "gg3ImRho_1410";
        fDictionary[21] = "gg4Kstr_1680"; fDictionary[22] = "gg4ImKstr_1680"; fDictionary[23] = "gg4Rho_1680";
        fDictionary[24] = "gg4ImRho_1680"; fDictionary[25] = "gg5Kstr_1430"; fDictionary[26] = "gg5ImKstr_1430";
        fDictionary[27] = "gg5Rho_1430"; fDictionary[28] = "gg5ImRho_1430"; fDictionary[29] = "kstr1270S_Re";
        fDictionary[30] = "kstr1270S_Im"; fDictionary[31] = "kstr1270D_Re"; fDictionary[32] = "kstr1270D_Im";
        fDictionary[33] = "kstr1400S_Re"; fDictionary[34] = "kstr1400S_Im"; fDictionary[35] = "kstr1400D_Re";
        fDictionary[36] = "kstr1400D_Im"; fDictionary[37] = "rho1270S_Re"; fDictionary[38] = "rho1270S_Im";
        fDictionary[39] = "rho1270D_Re"; fDictionary[40] = "rho1270D_Im"; fDictionary[41] = "rho1400S_Re";
        fDictionary[42] = "rho1400S_Im"; fDictionary[43] = "rho1400D_Re"; fDictionary[44] = "rho1400D_Im";
        GSetDefaultParameters();
        fDefVars = finVars;
        fGOF = std::make_shared<GGOFKolmogorovTest>();
        csvfile.open (tosave);
        csvfile << "min val," << "status," << "pValueCosTheta," << "pValuePhi," << "pValueSKPi," << "pValueSPiPi," << "pValueMKPiPi,";
        for(auto i = 0; i < fDictionary.size(); ++i)
            csvfile << fDictionary[i]+"," << fDictionary[i]+"_err,";
        csvfile << '\n';
    }
    
    GFitter::~GFitter()
    {
        csvfile.close();
        delete ea;
        delete mf;
    }
    void GFitter::GSetDefaultParameters()
    {
        std::shared_ptr<GFitterParser> fp(new GFitterParser("../additionalFiles/model_default.fit"));
        auto dictionary = fp->GParse();
        auto valPars = static_cast<int>(fDictionary.size());
        for(auto i = 0; i < valPars; ++i)
        {
            fnamesOfVars.push_back(fDictionary[i]);
            if(dictionary.find(fDictionary[i]) != dictionary.end())
            {
                auto record = dictionary[fDictionary[i]];
                finVars.push_back(std::get<0>(record));
                fkindOfVars.push_back((std::get<1>(record)));
                flowVars.push_back(std::get<2>(record));
                fupVars.push_back(std::get<3>(record));
                fstepsVars.push_back(std::get<4>(record));
//                 std::cout << varType << '\n';
            }
            else
            {
                
                throw "Check model_default.fit file!";
            }
        }
    }
    void GFitter::GSetModelParameters(const std::string& fileName)
    {
        
        std::shared_ptr<GFitterParser> fp(new GFitterParser(fileName));
        auto dictionary = fp->GParse();
        auto valPars = static_cast<int>(fDictionary.size());
        for(auto i = 0; i < valPars; ++i)
        {
            fnamesOfVars.push_back(fDictionary[i]);
            if(dictionary.find(fDictionary[i]) != dictionary.end())
            {
                auto record = dictionary[fDictionary[i]];
                finVars[i] = std::get<0>(record);
                fkindOfVars[i] = std::get<1>(record);
                flowVars[i] = std::get<2>(record);
                fupVars[i] = std::get<3>(record);
                fstepsVars[i] = std::get<4>(record);
//                 std::cout << varType << '\n';
            }
        }
        ea->GSetFunctionVariables(valPars, &finVars[0], &flowVars[0], &fupVars[0], &fstepsVars[0], fkindOfVars, fnamesOfVars, true);
    }
    
    void GFitter::GSetData(const std::string& dataFile, const std::string& cut)
    {
        pars.clear();
        fCut = cut;
        std::shared_ptr<GCutsParser> fp(new GCutsParser(cut));
        auto dictionary = fp->GParse();
        std::cout << dictionary["sKpipi"].first << "  " << dictionary["sKpipi"].second << '\n';
        double pi_const = 4*std::atan(1);
        auto hfile2 = std::make_shared<TFile>(dataFile.c_str());
        if(hfile2)
            std::cout << "HEMA" << '\n';
        TTree *T2 = (TTree*)hfile2->Get("DalitzEventList");
        double kaon4Vec[4] = {0};
        double pi14Vec[4] = {0};
        double pi24Vec[4] = {0};
        double gamma4Vec[4] = {0};
        if(fCharge == 1)
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
        auto nentries = T2->GetEntries();
        std::cout << nentries << '\n';
        
        int nPars = 5;
//         pars.resize(nPars*nentries+2);
        pars.resize(2);
//         pars[0] = nentries;
        pars[1] = nPars;
        
        auto nEvts = 0;
        for(auto ii = 0; ii < nentries; ii++)
        {
            T2->GetEntry(ii);
            GVector4D q_K(kaon4Vec[0], kaon4Vec[1], kaon4Vec[2], kaon4Vec[3]);
            GVector4D q_piplus(pi14Vec[0], pi14Vec[1], pi14Vec[2], pi14Vec[3]);
            GVector4D q_pizero(pi24Vec[0], pi24Vec[1], pi24Vec[2], pi24Vec[3]);
            GVector4D q_gamma(gamma4Vec[0], gamma4Vec[1], gamma4Vec[2], gamma4Vec[3]);
            auto&& phaseSpace = GKinematics::G4VectorsToPhaseSpace(q_K, q_piplus, q_pizero, q_gamma, fCharge);
            
//             std::cout << phaseSpace[4] << '\n';
            if( (std::sqrt(phaseSpace[4]) <= dictionary["sKpipi"].second) && (dictionary["sKpipi"].first <= std::sqrt(phaseSpace[4])) )
            {
                fGOF->GSetTruthEvent(phaseSpace);
                for(auto dim = 0; dim < nPars; ++dim)
                {
                    pars.push_back(phaseSpace[dim]);
//                  pars[2+ii*nPars+dim] = phaseSpace[dim];
                }   
                ++nEvts;
            }
        }
        
        std::cout << nEvts << "  " << nentries << "  " << pars.size() <<  '\n';
        pars[0] = nEvts;
    }
    void GFitter::GGenerateNormalizationIntegrals(const char* file, const std::string& cut)
    {
        std::shared_ptr<GCutsParser> fp(new GCutsParser(cut));
        auto dictionary = fp->GParse();
        auto&& low = dictionary["sKpipi"].first;
        auto&& up = dictionary["sKpipi"].second;
        if(up == 0)
        {
            up = std::numeric_limits<float>::max();
        }
        auto&& symPars = 26;
        const char* resNames[] = {"K1_1400", "K*_1410","K*_1680", "K2_1430","K1_1270"};
        int charges[] = {fCharge,fCharge,fCharge,fCharge,fCharge};
        
        ex p[symPars] =
        {
            symbol("gQPMC"), symbol("f2"),   symbol("thK1"), symbol("phDK"), 
            symbol("phSR"),  symbol("phDR"), symbol("gg"),       symbol("lambda"),
            symbol("ff"),    symbol("ff2"),  symbol("ff3"),      symbol("ff4"),
            symbol("gg3Kstr"), symbol("gg3Rho"), 
            symbol("gg4Kstr"), symbol("gg4Rho"),
            symbol("gg5Kstr"), symbol("gg5Rho"),
    
            symbol("kstr1270S"),    symbol("kstr1270D"),  symbol("kstr1400S"),      symbol("kstr1400D"),
            symbol("rho1270S"),    symbol("rho1270D"),  symbol("rho1400S"),      symbol("rho1400D")
        };
        
         auto valPars = 45;
        double modelPars3[valPars] = {
      4.,       4,    TMath::Pi()/3.,         0.,              0.,              0,    
   // gammaQPC(0), f2(1),     thetaK1(2),    phiDK*(3),      phiSRho(4),     phiDRho(5), 
         0.0,            0.7,          0.47,               0.,              0.0,
   //   ggKappaRe(6),   lambda(7),   ffRe_1400(8) ,    ffIm_1400(9),    ggKappaIm(10)
         0.78,             0.,                1.24,                   0.,               0.44,           0.,
//      ff2_Re_1410(11)  ff2_Im_1410(12)  ff3_Re_1680(13)   ff3_Im_1680(14)    ff4_Re_1430(15)  ff4_Im_1430(16)
       1.,                      0.,            1.,              0.,      
//     gg3Kstr_1410(17)  gg3ImKstr_1410(18)   gg3Rho_1410(19)   gg3ImRho_1410(20)  
       1.,                      0.,            1.,              0.,
//     gg4Kstr_1680(21)  gg4ImKstr_1680(22)    gg4Rho_1680(23)   gg4ImRho_1680(24)
       1.,                      0.,            1.,              0.,
//     gg5Kstr_1430(25)  gg5ImKstr_1430(26)    gg5Rho_1430(27)   gg5ImRho_1430(28)
       
         1.0,          0.0,              1.0,           0.0,
         1.0,          0.0,              1.0,           0.0,
         1.0,          0.0,              1.0,           0.0,
         1.0,          0.0,              1.0,           0.0
     };
     
     int charge = fCharge;
        auto ea = std::make_shared<GEventsAnalyzer>(charge);
        ea->GWriteNormalizationIntegrals(file, "rearchivate", symPars, resNames, p, &modelPars3[0], charges);        
        ea->GSetFunctionParameters(valPars, &modelPars3[0]);
        ea->GSetFunction(&Gamapola::GMathFunctions::GProcessingComputationOfPDF);
        ea->GSetDecayMode(charge);
        ea->GSetEventsNumber(5000);
     
        ea->GGenerateEvents(low, up);
    }
    
    void GFitter::GSetNormalizationIntegrals(const std::string& fileName)
    {
        auto&& symPars = 26;
        const int nRes = 5;
        const char* resNames[] = {"K1_1400", "K*_1410","K*_1680", "K2_1430","K1_1270"};
        int charges[] = {fCharge,fCharge,fCharge,fCharge,fCharge};
        ex p[symPars] =
        {
            symbol("gQPMC"), symbol("f2"),   symbol("thK1"), symbol("phDK"), 
            symbol("phSR"),  symbol("phDR"), symbol("gg"),       symbol("lambda"),
            symbol("ff"),    symbol("ff2"),  symbol("ff3"),      symbol("ff4"),
            symbol("gg3Kstr"), symbol("gg3Rho"), 
            symbol("gg4Kstr"), symbol("gg4Rho"),
            symbol("gg5Kstr"), symbol("gg5Rho"),
    
            symbol("kstr1270S"),    symbol("kstr1270D"),  symbol("kstr1400S"),      symbol("kstr1400D"),
            symbol("rho1270S"),    symbol("rho1270D"),  symbol("rho1400S"),      symbol("rho1400D")
        };
        ea->GReadNormalizationIntegrals(fileName.c_str(), symPars, resNames, p, &finVars[0], charges);
    }
    
    void GFitter::GFit(const bool& isGOF)
    {
        auto valPars = 45;
        int charges[] = {fCharge,fCharge,fCharge,fCharge,fCharge};
        ea->GSetFunctionParameters((int)pars.size(), &pars[0]);
//         ea->GSetFunctionVariables(valPars, &finVars[0], &flowVars[0], &fupVars[0], &fstepsVars[0], fkindOfVars, fnamesOfVars, true);
        ea->GSetDecayMode(charges[0]);
        ea->GAnalyzeEvents();
        if(isGOF)
        {
            auto&& fp = ea->GGetFitParameters();
            auto&& fpe = ea->GGetFitParametersErrors();
            auto&& minVal = ea->GGetMinFuncValue();
            auto&& status = ea->GGetStatus();
            fGOF->GSetGeneratedData(fp, 5000, fCharge, fCut);
            fGOF->GMakeTest();
            auto&& res = fGOF->GGetpValues();
            csvfile << minVal << "," << status << ",";
            for(auto el: res)
                csvfile << el << ",";
            for(auto&& i = 0; i < valPars; ++i)
            {
                csvfile << fp[i] << "," << fpe[i] << ",";
                if(fkindOfVars[i] == "fixed")
                    continue;
                std::cout << std::left << std::setw(15) << fDictionary[i] << std::left << std::setw(15) << ":    truth: " << fDefVars[i] << std::left << std::setw(15) << "     fit: " << fp[i] << std::left << std::setw(15) << "    error: " << fpe[i] << '\n';
            }
            csvfile << '\n';
        }
    }
}

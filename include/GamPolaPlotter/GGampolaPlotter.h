#ifndef GGAMPOLAPLOTTER_H
#define GGAMPOLAPLOTTER_H

#include <unordered_map>
#include "TH1D.h"
#include "TH2D.h"
#include "EvtGenModels/EvtGVector4D.hh"
#include "GRootTreeReader.h"
#include "GamPolaParser/GPlotterParser.h"
#include "EvtGenModels/EvtGVector3D.hh"
#include <functional>

namespace Gamapola
{
    class GGampolaPlotter
    {
    public:
        explicit GGampolaPlotter(const std::string& input_filename)
        {
            auto&& pp = std::make_shared<GPlotterParser>(input_filename);
            pp->GParse();
            fDict1D = pp->GGetDictionary1D();
            for(const auto& el : fDict1D)
            {
                GDistributeFuncPtrs(el.first);
                fHisto1D.insert(std::make_pair(std::string(el.first), std::move(std::make_shared<TH1D>(el.first.c_str(), el.first.c_str(), std::get<1>(el.second), std::get<2>(el.second),std::get<3>(el.second)))));
            }
            fDict1DCut = pp->GGetDictionaryCut1D();
            for(const auto& el : fDict1DCut)
            {
                GDistributeFuncPtrs(std::get<2>(el.second));
                fHisto1D.insert(std::make_pair(std::string(el.first), std::move(std::make_shared<TH1D>(el.first.c_str(), el.first.c_str(), std::get<3>(el.second), std::get<4>(el.second),std::get<5>(el.second)))));
                
            }
            
            fDict2D = pp->GGetDictionary2D();
            for(const auto& el : fDict2D)
            {
                GDistributeFuncPtrs(el.first);
                fHisto2D.insert(std::make_pair(std::string(el.first), std::move(std::make_shared<TH2D>(el.first.c_str(), el.first.c_str(), std::get<1>(el.second), std::get<2>(el.second),std::get<3>(el.second), std::get<4>(el.second), std::get<5>(el.second),std::get<6>(el.second)))));
            }
        }
        void GDistributeFuncPtrs(const std::string& distr)
        {
            std::cout << distr << '\n';
            if(distr == "sKpipi")
                fFuncPtrMap.insert(std::make_pair(distr, std::bind(&GGampolaPlotter::GFillSKPiPi, this)));
            if(distr == "sKpi1")
                fFuncPtrMap.insert(std::make_pair(distr, std::bind(&GGampolaPlotter::GFillSKPi1, this)));
            if(distr == "cosTheta")
                fFuncPtrMap.insert(std::make_pair(distr, std::bind(&GGampolaPlotter::GFillCosTheta, this)));
            if(distr == "sKpi1Kpi2")
                fFuncPtrMap.insert(std::make_pair(distr, std::bind(&GGampolaPlotter::GFillKPi1vsKPi2, this)));
        }
        void plot(const std::string& data, const int& charge, const std::string& out)
        {
            fTreeReader = std::make_shared<GRootTreeReader>(data, charge);
            auto&& nEntries = fTreeReader->GGetNEntries();
            auto&& ffile = std::make_shared<TFile>(out.c_str(), "recreate") ;
            for(auto&& entry = 0; entry < nEntries; ++entry)
            {
                fTreeReader->GGetEntry(entry);
                auto&& kaon = fTreeReader->GGetKaonFourVec();
                auto&& pion1 = fTreeReader->GGetPion1FourVec();
                auto&& pion2 = fTreeReader->GGetPion2FourVec();
                auto&& kpipi = kaon+pion1+pion2;
                fControlValue = std::sqrt(kpipi*kpipi);
                for(auto& el : fFuncPtrMap)
                    el.second();
            }
            for(auto& el : fHisto1D)
                el.second->Write();
            for(auto& el : fHisto2D)
                el.second->Write();
            
            ffile->Close();
        }
        void GFillSKPiPi()
        {
            auto&& kaon = fTreeReader->GGetKaonFourVec();
            auto&& pion1 = fTreeReader->GGetPion1FourVec();
            auto&& pion2 = fTreeReader->GGetPion2FourVec();
            auto&& kpipi = kaon+pion1+pion2;
            Update1D("sKpipi", std::sqrt(kpipi*kpipi));
        }
        void GFillSKPi1()
        {
            auto&& kaon = fTreeReader->GGetKaonFourVec();
            auto&& pion1 = fTreeReader->GGetPion1FourVec();
            auto&& kpi1 = kaon+pion1;
            Update1D("sKpi1", std::sqrt(kpi1*kpi1));
        }
        void GFillCosTheta()
        {
            auto&& qK = fTreeReader->GGetKaonFourVec();
            auto&& qpi1 = fTreeReader->GGetPion1FourVec();
            auto&& qpi2 = fTreeReader->GGetPion2FourVec();
            auto&& k1new = qK+qpi1+qpi2;
            GVector3D boost2Kres(-k1new.get(1)/k1new.get(0),-k1new.get(2)/k1new.get(0),-k1new.get(3)/k1new.get(0));
            auto&& qGamma_Kres = boostTo(fTreeReader->GGetGammaFourVec(),boost2Kres,false);
            auto&& qpi1_Kres = boostTo(qpi1,boost2Kres,false);
            auto&& qpi2_Kres = boostTo(qpi2,boost2Kres,false);
            
            GVector4D ez(1,-qGamma_Kres.get(1),-qGamma_Kres.get(2),-qGamma_Kres.get(3));
            ez = ez/qGamma_Kres.d3mag();
            auto&& n = qpi1_Kres.cross(qpi2_Kres)/qpi1_Kres.cross(qpi2_Kres).d3mag();
            Update1D("cosTheta", n.dot(ez));
        }
        void GFillKPi1vsKPi2()
        {
            auto&& kaon = fTreeReader->GGetKaonFourVec();
            auto&& pion1 = fTreeReader->GGetPion1FourVec();
            auto&& pion2 = fTreeReader->GGetPion2FourVec();
            auto&& kpi1 = (kaon+pion1);
            auto&& kpi2 = (kaon+pion2);
            auto&& skpi1 = kpi1*kpi1;
            auto&& skpi2 = kpi2*kpi2;
            Update2D("sKpi1Kpi2", skpi1, skpi2);
        }
        
        void Update1D(const std::string& entry, const double& value)
        {
            std::cout << entry << '\n';
            for(const auto& el : fDict1D)
                if(std::get<0>(el.second) == entry)
                    fHisto1D[el.first]->Fill(value);
            for(const auto& el : fDict1DCut)
                if( (std::get<2>(el.second) == entry) && (fControlValue >= std::get<0>(el.second)) && (fControlValue <= std::get<1>(el.second)))
                    fHisto1D[el.first]->Fill(value);
        }
        void Update2D(const std::string& entry, const double& valueX, const double& valueY)
        {
            for(const auto& el : fDict2D)
                if(std::get<0>(el.second) == entry)
                    fHisto2D[el.first]->Fill(valueX, valueY);
        }
        ~GGampolaPlotter()
        {}
    private: 
        std::shared_ptr<GRootTreeReader> fTreeReader;
        std::unordered_map<std::string, std::shared_ptr<TH1D>> fHisto1D;
        std::unordered_map<std::string, std::shared_ptr<TH2D>> fHisto2D;
        std::unordered_map<std::string, std::function<void()>> fFuncPtrMap;
        Dict1DPlain fDict1D;
        Dict1DCut fDict1DCut;
        Dict2DPlain fDict2D;
        double fControlValue;
    };
}

#endif

#include "GGOFKolmogorovTest.h"
#include "GGenerator.h"
#include "TMath.h"
namespace Gamapola{
    void GGOFKolmogorovTest::GMakeTest()
    {
        fpValues.clear();
        auto&& nTruth = static_cast<int>(fTruthData[0].size());
        auto&& nGenerated = static_cast<int>(fGeneratedData[0].size());
        std::array<double, 5> stats;
        for(auto&& iTest = 0; iTest < 5; ++iTest)
        {
            
            auto&& truthDistr = fTruthData[iTest];
            auto&& genData = fGeneratedData[iTest];
            std::sort (truthDistr.begin(), truthDistr.end());
            std::sort (genData.begin(), genData.end());
            auto&& stat = TMath::KolmogorovTest(nTruth, &truthDistr[0], nGenerated, &genData[0], "D");
            fpValues.emplace_back(stat);
        }
    }
    void GGOFKolmogorovTest::GSetTruthEvent(const std::vector<double>& phaseSpace)
    {
        for(auto i = 0; i < static_cast<int>(phaseSpace.size()); ++i)
        {
            fTruthData[i].emplace_back(phaseSpace[i]);
        }
//         std::cout << "fTruthData: " << fTruthData[0].size() << '\n';
    }
    void GGOFKolmogorovTest::GSetGeneratedData(const std::vector<double>& modelPars, const int& nEvts, const int& charge, const std::string& cut)
    {
        auto&& nEvents = static_cast<int>(fTruthData[0].size());
        auto g1 = std::make_shared<GGenerator>(charge);
        g1->GSetCouplings(modelPars);
        g1->GGenerate(nEvents, cut);
        g1->GWriteToContainer();
        g1->GWriteToFile("neutral_gen_test.root");
        fGeneratedData = g1->GGetData();
    }
}

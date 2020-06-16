#ifndef GGOFKOLMOGOROVTEST
#define GGOFKOLMOGOROVTEST

#include "GGoodnessOfFit.h"
#include <vector>
#include <array>
#include <memory>
namespace Gamapola{
    class GGOFKolmogorovTest: public GGoodnessOfFit
    {
    public: 
        virtual void GMakeTest();
        virtual void GSetTruthEvent(const std::vector<double>& phaseSpace);
        virtual void GSetGeneratedData(const std::vector<double>& modelPars, const int& nEvts, const int& charge, const std::string& cut);
        const std::vector<double>& GGetpValues() const {return fpValues;}
    private:
        std::array<std::vector<double>,5> fTruthData;
        std::array<std::vector<double>,5> fGeneratedData;
        std::vector<double> fpValues;
    };
}

#endif

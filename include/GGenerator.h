#ifndef GGENERATOR_H
#define GGENERATOR_H

#include <memory>
#include <unordered_map>
#include <array>
#include "GEventsGenerator.h"
#include "EvtGenModels/EvtGMathFunctions.hh"

namespace Gamapola
{
    class GGenerator
    {
    public:
        GGenerator(const int& charge);
        ~GGenerator();
        void GSetCouplings(const std::string& fileName);
        void GSetCouplings(const std::vector<double>& couplings);
        void GGenerate(const int& nEvts, const std::string& cut = "");
        void GWriteToFile(const std::string& outFile);
        void GWriteToContainer();
        const std::array<std::vector<double>, 5>& GGetData() const { return fEvents; }
        std::shared_ptr<GEventsGenerator> GGetEvtGen()
        {
            return eg;
        }
        std::shared_ptr<GMathFunctions> GGetMath()
        {
            return mf;
        }
    private:
        std::shared_ptr<GMathFunctions> mf;
        std::shared_ptr<GEventsGenerator> eg;
        std::vector<double> modelPars;
        int fCharge;
        std::array<std::vector<double>, 5> fEvents;
    };
}

#endif

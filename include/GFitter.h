#ifndef GFITTER_H
#define GFITTER_H

#include <memory>
#include "GEventsAnalyzer.h"
#include "EvtGenModels/EvtGMathFunctions.hh"
#include "GGenerator.h"
#include "GGoodnessOfFit.h"

namespace Gamapola
{
    class GFitter
    {
    public:
        GFitter(const int& charge, const std::string& tosave);
        ~GFitter();
        void GSetModelParameters(const std::string& fileName);
        void GSetData(const std::string& inputFileName, const std::string& cut);
        void GSetNormalizationIntegrals(const std::string& fileName);
        void GGenerateNormalizationIntegrals(const char* file, const std::string& cut);
        void GFit(const bool& isGOF = false);
        GEventsAnalyzer* GGetEvtAnalyzer() { return ea; }
    private:
        void GSetDefaultParameters();
        GEventsAnalyzer* ea;
        GMathFunctions* mf;
        std::vector<double> finVars;
        std::vector<double> fDefVars;
        std::vector<std::string> fkindOfVars;
        std::vector<std::string> fnamesOfVars;
        std::vector<double> flowVars;
        std::vector<double> fupVars;
        std::vector<double> fstepsVars;
        std::vector<double> pars;
        std::unordered_map<int, std::string> fDictionary; 
        int fCharge;
        std::shared_ptr<GGoodnessOfFit> fGOF;
        std::string fCut;
        std::ofstream csvfile;
    };
}

#endif

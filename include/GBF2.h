#ifndef GBF2_HH
#define GBF2_HH

#include <vector>
#include <string>
#include <memory>
#include <unordered_map>
#include "EvtGenModels/EvtGMathFunctions.hh"
#include "GEventsAnalyzer.h"
#include "GFitter.h"

namespace Gamapola{
    class GBF2
    {
    public: 
        GBF2(const int& charge, const std::string& modelPars, const std::string& normInts);
        virtual ~GBF2();
        void GSetModelParameters(const std::string& fileName);
        void GHistogramize(const int& nBins, const double& low, const double& up, const std::string& infile, const std::string& outfile);
        void GHistogramize(const std::unordered_map<std::string, std::tuple<int, double, double>>&, const std::string& infile, const std::string& outfile, const int& nEvts);
        double GChi2(double* x, double* p);
        void GSetData(const std::string& data, const std::string& key);
        void GFit(const std::string& filename);
    private:
        void GSetDefaultParameters();
        std::vector<std::string> GParseCSV(std::ifstream& data);
        void GFillMC(const std::string& str);
        void GFillData();
        
        double GRatio(double* x, double*p);
        int GGetBin(const double& s, const double& sMin, const double& sMax, const int& nBins);
        std::shared_ptr<GMathFunctions> fMf;
        int fcharge;
        std::vector<double> model;
        std::unordered_map<int, std::string> fDictionary; 
        std::shared_ptr<GFitter> fFitter;
        std::unordered_map<int,std::vector<std::vector<std::vector<double>>>> fPhaseSpace;
        std::unordered_map<int,std::vector<std::vector<double>>> fPDF; 
        GEventsAnalyzer* ea;
        double fTempPDF;
        std::unordered_map<int,std::vector<double>> fHeights;
        std::unordered_map<int,std::vector<double>> fHeightErrs;
        std::unordered_map<int,double> fNdata;
        
        std::vector<std::string> fvalidKeys;
        std::unordered_map<std::string, int> fwordMap;
        std::unordered_map<std::string, std::string> fdataMap;
        
        std::vector<double> finVars;
        std::vector<double> fDefVars;
        std::vector<std::string> fkindOfVars;
        std::vector<std::string> fnamesOfVars;
        std::vector<double> flowVars;
        std::vector<double> fupVars;
        std::vector<double> fstepsVars;
        int fNDataTotal;
    };
}
#endif

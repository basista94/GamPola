#ifndef GPLOTTERPARSER_H
#define GPLOTTERPARSER_H

#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_map>

namespace Gamapola
{
    using Dict1DPlain = std::unordered_map<std::string, std::tuple<std::string, int, double, double>>;
    using Dict1DCut = std::unordered_map<std::string, std::tuple<double, double, std::string, int, double, double>>;
    using Dict2DPlain = std::unordered_map<std::string, std::tuple<std::string, int, double, double, int, double, double>>;
    class GPlotterParser
    {
    public:
        explicit GPlotterParser(const std::string& filename)
        {
            fFileName = filename;
        }
        virtual void GParse()
        {
            std::ifstream file;
            file.open(fFileName);
            std::string name;
            while (std::getline(file, name))
            {
                name.erase(remove_if(name.begin(), name.end(), isspace), name.end());
                if(name == "" || name[0] == std::string("#"))
                    continue;
                GParseHeader(name);
            }
        }
        const Dict1DPlain& GGetDictionary1D() const { return fDictionary; }
        const Dict1DCut& GGetDictionaryCut1D() const { return fDictionaryCut; }
        const Dict2DPlain& GGetDictionary2D() const { return fDictionary2D; }
    private:
        void GParseHeader(const std::string& str)
        {
            auto&& header = GExtractSubstringBefore(str, ":");
            auto&& dimensionlity = GExtractSubstringBefore(header, ",");
            auto&& rangeFlag = GExtractSubstringAfter(header, ",");
            if(dimensionlity == "1D" && rangeFlag == "All")
                GParseString(GExtractSubstringAfter(str, ":"));
            if(dimensionlity == "2D" && rangeFlag == "All")
                GParseString2D(GExtractSubstringAfter(str, ":"));
            if(dimensionlity == "1D" && rangeFlag == "cut")
                GParseCutString(GExtractSubstringAfter(str, ":"));
        }
        void GParseString(const std::string& str)
        {
            auto&& label = GExtractSubstringBefore(str, ":");
            auto&& rest0 = GExtractSubstringAfter(str, ":");
//             std::cout << str << '\n';
            auto&& parName = GExtractSubstringBefore(rest0, ",");
            auto&& rest1 = GExtractSubstringAfter(rest0, ",");
            
            auto&& parBins = stoi(GExtractSubstringBefore(rest1, ","));
            auto&& rest2 = GExtractSubstringAfter(rest1, ",");

            auto&& parLow = stod(GExtractSubstringBefore(rest2, ","));
            auto&& rest3 = GExtractSubstringAfter(rest2, ",");
            
            auto&& parUp = stod(GExtractSubstringBefore(rest3, ","));
//             std::cout << label << "   " << parName << "  " << parBins << "  " << parLow << "  " << parUp << '\n';
            fDictionary[label] = std::make_tuple(parName, parBins, parLow, parUp);
        }
        void GParseCutString(const std::string& str)
        {
            auto&& cutRest = GExtractSubstringBefore(str, ":");
            auto&& low = stod(GExtractSubstringBefore(cutRest, "<sKpipi<"));
            auto&& up = stod(GExtractSubstringAfter(cutRest, "<sKpipi<"));
            
            auto&& substr = GExtractSubstringAfter(str, ":");
            
            auto&& label = GExtractSubstringBefore(substr, ":");
            auto&& rest0 = GExtractSubstringAfter(substr, ":");
//             std::cout << str << '\n';
            auto&& parName = GExtractSubstringBefore(rest0, ",");
            auto&& rest1 = GExtractSubstringAfter(rest0, ",");
            
            auto&& parBins = stoi(GExtractSubstringBefore(rest1, ","));
            auto&& rest2 = GExtractSubstringAfter(rest1, ",");

            auto&& parLow = stod(GExtractSubstringBefore(rest2, ","));
            auto&& rest3 = GExtractSubstringAfter(rest2, ",");
            
            auto&& parUp = stod(GExtractSubstringBefore(rest3, ","));
            
            std::cout << low << "  " << up << "  " << cutRest << "   " << parName << "  " << parBins << "  " << parLow << "  " << parUp << '\n';
            fDictionaryCut[label] = std::make_tuple(low, up, parName, parBins, parLow, parUp);
        }
        void GParseString2D(const std::string& str)
        {
            auto&& label = GExtractSubstringBefore(str, ":");
            auto&& rest0 = GExtractSubstringAfter(str, ":");
//             std::cout << str << '\n';
            auto&& parName = GExtractSubstringBefore(rest0, ",");
            auto&& rest1X = GExtractSubstringAfter(rest0, ",");
            
            auto&& parBinsX = stoi(GExtractSubstringBefore(rest1X, ","));
            auto&& rest2X = GExtractSubstringAfter(rest1X, ",");

            auto&& parLowX = stod(GExtractSubstringBefore(rest2X, ","));
            auto&& rest3X = GExtractSubstringAfter(rest2X, ",");
            
            auto&& parUpX = stod(GExtractSubstringBefore(rest3X, ","));
            auto&& rest1Y = GExtractSubstringAfter(rest3X, ",");
            
            auto&& parBinsY = stoi(GExtractSubstringBefore(rest1Y, ","));
            auto&& rest2Y = GExtractSubstringAfter(rest1Y, ",");

            auto&& parLowY = stod(GExtractSubstringBefore(rest2Y, ","));
            auto&& rest3Y = GExtractSubstringAfter(rest2Y, ",");
            
            auto&& parUpY = stod(GExtractSubstringBefore(rest3Y, ","));
//             fDictionary[parName] = std::make_tuple(parBins, parLow, parUp);
//             std::cout << label << "   " << parName << "  " << parBinsX << "  " << parLowX << "  " << parUpX << "  " << parBinsY << "  " << parLowY << "  " << parUpY << "  " << '\n';
            fDictionary2D[label] = std::make_tuple(parName, parBinsX, parLowX, parUpX, parBinsY, parLowY, parUpY);
        }
        std::string GExtractSubstringBefore(const std::string& str, const std::string& spec)
        {
            auto&& parName = str.substr(0, str.find(spec));
            parName.erase(remove_if(parName.begin(), parName.end(), isspace), parName.end());
            return parName;
        }
        std::string GExtractSubstringAfter(const std::string& str, const std::string& spec)
        {
            auto&& parStrVal = str.substr(str.find(spec) + spec.length());
            return parStrVal;
        }
        std::string fFileName; 
        std::vector<std::string> fParNames;
        std::vector<double> fParValues;
        Dict1DPlain fDictionary;
        Dict1DCut fDictionaryCut;
        Dict2DPlain fDictionary2D;
    };
}
#endif

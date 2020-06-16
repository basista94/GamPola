#ifndef GCUTSPARSER_H
#define GCUTSPARSER_H

namespace Gamapola
{
    class GCutsParser
    {
    public:
        explicit GCutsParser(const std::string& cut)
        {
            fCut = cut;
        }
        virtual std::unordered_map<std::string, std::pair<double, double>> GParse()
        {
            fCut.erase(remove_if(fCut.begin(), fCut.end(), isspace), fCut.end());
            if(fCut == "")
                return fDictionary;
            auto&& parLow = stod(GExtractSubstringBefore(fCut, "<"));

            auto&& rest1 = GExtractSubstringAfter(fCut, "<");
            
            auto&& parName = GExtractSubstringBefore(rest1, "<");
            auto&& parUp = stod(GExtractSubstringAfter(rest1, "<"));
            
            fDictionary[parName] = std::make_pair(parLow, parUp);
            return fDictionary;
        }
        
    private:
        std::string GExtractSubstringBefore(const std::string& str, const std::string& spec)
        {
            auto&& parName = str.substr(0, str.find(spec));
            parName.erase(remove_if(parName.begin(), parName.end(), isspace), parName.end());
            return parName;
        }
        std::string GExtractSubstringAfter(const std::string& str, const std::string& spec)
        {
            auto&& parStrVal = str.substr(str.find(spec) + 1);
            return parStrVal;
        }
        std::string fCut; 
        std::unordered_map<std::string, std::pair<double, double>> fDictionary;
    };
}

#endif

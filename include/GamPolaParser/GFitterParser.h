#ifndef GFITTERPARSER_H
#define GFITTERPARSER_H

namespace Gamapola
{
    class GFitterParser
    {
    public:
        explicit GFitterParser(const std::string& filename)
        {
            fFileName = filename;
        }
        virtual std::unordered_map<std::string, std::tuple<double, std::string, double, double, double>> GParse()
        {
            std::ifstream file;
            file.open(fFileName);
            std::string name;
            while (std::getline(file, name))
            {
                name.erase(remove_if(name.begin(), name.end(), isspace), name.end());
                if(name == "" || name[0] == std::string("#"))
                    continue;
                GParseString(name);
            }
            
            return fDictionary;
        }
        
    private:
        void GParseString(const std::string& str)
        {
//             std::cout << str << '\n';
            auto&& parName = GExtractSubstringBefore(str, ":");
//             if(parName == "cut")
//                 std::cout << "Goodbye" << '\n';
//                 return;
            auto&& rest1 = GExtractSubstringAfter(str, ":");
            
            auto&& parVal = stod(GExtractSubstringBefore(rest1, ","));
            auto&& rest2 = GExtractSubstringAfter(rest1, ",");
            
            auto&& parType = GExtractSubstringBefore(rest2, ",");
            auto&& rest3 = GExtractSubstringAfter(rest2, ",");
            
            auto&& parLow = stod(GExtractSubstringBefore(rest3, ","));
            auto&& rest4 = GExtractSubstringAfter(rest3, ",");
            
            auto&& parUp = stod(GExtractSubstringBefore(rest4, ","));
            auto&& rest5 = GExtractSubstringAfter(rest4, ",");
            
            auto&& parStep = stod(GExtractSubstringBefore(rest5, ";"));
//             std::cout << parName << "  " << parVal << "  " << parType << "  " << parLow << "  " << parUp << "   " << parStep << '\n';
            fDictionary[parName] = std::make_tuple(parVal, parType, parLow, parUp, parStep);
        }
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
        std::string fFileName; 
        std::vector<std::string> fParNames;
        std::vector<double> fParValues;
        std::unordered_map<std::string, std::tuple<double, std::string, double, double, double>> fDictionary;
    };
}

#endif

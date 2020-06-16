#ifndef GGENERATORPARSER_h
#define GGENERATORPARSER_h

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_map>

namespace Gamapola
{    
    class GGeneratorParser
    {
    public:
        explicit GGeneratorParser(const std::string& filename)
        {
            fFileName = filename;
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
        }
        virtual std::vector<double> GParse()
        {
            std::ifstream file;
            file.open(fFileName);
            std::string name;
            while (std::getline(file, name))
            {
                if(name == "")
                    continue;
                GParseString(name);
            }
            
            auto valPars = static_cast<int>(fDictionary.size());
            std::vector<double> model;
            for(auto i = 0; i < valPars; ++i)
            {
                if(fDictionaryy.find(fDictionary[i]) != fDictionaryy.end())
                    model.emplace_back(fDictionaryy[fDictionary[i]]);
                else
                {
                    std::cout << fDictionary[i] << '\n';
                    throw "Check parameters!";
                }
            }
     
            return model;
        }
        
    private:
        void GParseString(const std::string& str)
        {
            auto&& parName = GExtractParName(str);
            auto&& parVal = GExtractParValue(str);
            fDictionaryy[parName] = parVal;
        }
        std::string GExtractParName(const std::string& str)
        {
            auto&& parName = str.substr(0, str.find(":"));
            parName.erase(remove_if(parName.begin(), parName.end(), isspace), parName.end());
            return parName;
        }
        double GExtractParValue(const std::string& str)
        {
            auto&& parStrVal = str.substr(str.find(":") + 1);
            return stod(parStrVal);
        }
        std::string fFileName; 
        std::vector<std::string> fParNames;
        std::vector<double> fParValues;
        std::unordered_map<std::string, double> fDictionaryy;
        std::unordered_map<int, std::string> fDictionary; 
    };
}

#endif

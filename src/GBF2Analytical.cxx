#include "GBF2Analytical.h"
#include "GRootTreeReader.h"
#include "GamPolaParser/GGeneratorParser.h"
#include "GamPolaParser/GFitterParser.h"
#include "GGenerator.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMinuitMinimizer.h"
#include "Math/WrappedMultiTF1.h"
#include <string>
#include <assert.h>

using namespace GiNaC;

namespace Gamapola{

GBF2Analytical::GBF2Analytical(const int& charge, const std::string& modelPars, const std::string& normInts):
fMf(new GMathFunctions()),
fFitter(new GFitter(charge, normInts+".csv")),
fcharge(charge),
fvalidKeys({"cosTheta", "phi", "Mkpipi", "Mkpi", "Mpipi", "Mkpi2"}),
fwordMap({{"cosTheta", 0}, {"phi", 1}, {"Mpipi", 3}, {"Mkpi", 2}, {"Mkpipi",4}, {"Mkpi2",5}}),
fNDataTotal(0)
{   
//     fcharge = charge;
//     fMf = std::make_shared<GMathFunctions>();
    fMf->GSetResonances(fcharge);
    
    fFitter->GSetModelParameters(modelPars);    
    fFitter->GSetNormalizationIntegrals(normInts);
    fFitter->GGetEvtAnalyzer()->GSetDecayMode(fcharge);
    ea = fFitter->GGetEvtAnalyzer();

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
        GSetDefaultParameters();
        fDefVars = finVars;
}

GBF2Analytical::~GBF2Analytical()
{
}
void GBF2Analytical::GSetModelParameters(const std::string& fileName)
{
        
        std::shared_ptr<GFitterParser> fp(new GFitterParser(fileName));
        auto dictionary = fp->GParse();
        auto valPars = static_cast<int>(fDictionary.size());
        for(auto i = 0; i < valPars; ++i)
        {
            fnamesOfVars.push_back(fDictionary[i]);
            if(dictionary.find(fDictionary[i]) != dictionary.end())
            {
                auto record = dictionary[fDictionary[i]];
                finVars[i] = std::get<0>(record);
                fkindOfVars[i] = std::get<1>(record);
                flowVars[i] = std::get<2>(record);
                fupVars[i] = std::get<3>(record);
                fstepsVars[i] = std::get<4>(record);
//                 std::cout << varType << '\n';
            }
        }
        ea->GSetFunctionVariables(valPars, &finVars[0], &flowVars[0], &fupVars[0], &fstepsVars[0], fkindOfVars, fnamesOfVars, true);
}
    
void GBF2Analytical::GSetDefaultParameters()
{
        std::shared_ptr<GFitterParser> fp(new GFitterParser("../additionalFiles/model_default.fit"));
        auto dictionary = fp->GParse();
        auto valPars = static_cast<int>(fDictionary.size());
        for(auto i = 0; i < valPars; ++i)
        {
            fnamesOfVars.push_back(fDictionary[i]);
            if(dictionary.find(fDictionary[i]) != dictionary.end())
            {
                auto record = dictionary[fDictionary[i]];
                finVars.push_back(std::get<0>(record));
                fkindOfVars.push_back((std::get<1>(record)));
                flowVars.push_back(std::get<2>(record));
                fupVars.push_back(std::get<3>(record));
                fstepsVars.push_back(std::get<4>(record));
//                 std::cout << varType << '\n';
            }
            else
            {
                
                throw "Check model_default.fit file!";
            }
        }
}

int GBF2Analytical::GGetBin(const double& s, const double& sMin, const double& sMax, const int& nBins)
{
    return static_cast<int>((s - sMin) / (sMax - sMin) * nBins);
}
int GBF2Analytical::GGetLHCbBin(const double& mKpipi)
{
    if(mKpipi < 1.3 && mKpipi >= 1.1)
        return 1;
    else if(mKpipi < 1.4 && mKpipi >= 1.3)
        return 2;
    else if(mKpipi < 1.6 && mKpipi >= 1.4)
        return 3;
    else if(mKpipi < 1.9 && mKpipi >= 1.6)
        return 4;
    else
        return -1;
}
void GBF2Analytical::GHistogramize(const std::unordered_map<std::string, std::tuple<int, double, double>>& dict, const std::string& infile, const std::string& outfile, const int& nEvts)
{
    auto g1 = std::make_shared<GGenerator>(fcharge);
    g1->GSetCouplings(infile);
    std::shared_ptr<GGeneratorParser> gp(new GGeneratorParser(infile));
    auto&& model = gp->GParse();
    g1->GGenerate(nEvts);
    auto&& eg = g1->GGetEvtGen();
    std::ofstream csvfile;
    csvfile.open (outfile);
    csvfile << "pdf," << "cosTheta," << "phi," << "sKPi," << "sPiPi," << "sKpipi,";
    for(auto el: dict)
    {
        if(std::find(fvalidKeys.begin(), fvalidKeys.end(), el.first) != fvalidKeys.end())
        {
            std::cout << el.first << '\n';
            csvfile << std::string("bin_")+el.first+std::string(",");
        }
    }
    csvfile << "bin_LHCb" << '\n';
    for(auto&& entry = 0; entry < nEvts; ++entry)
    {
        auto&& s = eg->GGetM()[entry];
        auto&& sjk = eg->GGetSjk()[entry];
        auto&& sik = eg->GGetSij()[entry];
        auto&& phi = eg->GGetPhiPDF()[entry];
        auto&& cosTheta = eg->GGetCosThetaPDF()[entry];
        auto sKpi2 = fMf->GSij(s, sik, sjk);
        
        auto&& phaseSpace = std::vector<double>{cosTheta, phi, sik, sjk, s};
        auto&& pdf = ea->GPDF(&model[0], &phaseSpace[0]);
        
//         phaseSpace[5] = std::sqrt(phaseSpace[5]);
        phaseSpace[4] = std::sqrt(phaseSpace[4]);
        phaseSpace[3] = std::sqrt(phaseSpace[3]);
        phaseSpace[2] = std::sqrt(phaseSpace[2]);
        phaseSpace.emplace_back(std::sqrt(sKpi2));
        
        auto&& isBroken = false;
        std::vector<int> binsVec;
        for(auto el: dict)
        {
            if(std::find(fvalidKeys.begin(), fvalidKeys.end(), el.first) != fvalidKeys.end())
            {
                auto iBin = GGetBin(phaseSpace[fwordMap[el.first]], std::get<1>(el.second), std::get<2>(el.second), std::get<0>(el.second));
                if(iBin < 0 || iBin >= std::get<0>(el.second))
                {
                    isBroken = true;
                    break;
                }
                binsVec.emplace_back(iBin);
            }
        }
        if(!isBroken)
        {
            auto&& lhcbBinning = GGetLHCbBin(phaseSpace[4]);
//             auto iBin = GGetBin(phaseSpace[0], -1., 1., 20);
            csvfile << pdf << "," << cosTheta << "," << phi << "," << sik << "," << sjk << "," << s << ",";
            for(auto el: binsVec)
                csvfile << el << ",";
             csvfile << lhcbBinning << "," << '\n';
        }
    }
    csvfile.close();
}
void GBF2Analytical::GSetData(const std::string& data, const std::string& key)
{
    fdataMap[key] = data;
//     auto mkpipiH = GFillData(data);
//     GFillMC(mc);
}
void GBF2Analytical::GFillData()
{
    fNDataTotal = 0;
    for(auto el : fdataMap)
    {
        fNdata[fwordMap[el.first]] = 0;
        std::ifstream data_if;
        data_if.open(el.second);
        auto&& row1 = GParseCSV(data_if);
    
        auto&& nCols = static_cast<int>(row1.size());
        while(!data_if.eof())
        {
            auto&& row1 = GParseCSV(data_if);
            if(static_cast<int>(row1.size()) == nCols)
            {
                fHeights[fwordMap[el.first]].emplace_back(std::stod(row1[0]));
                fHeightErrs[fwordMap[el.first]].emplace_back(std::stod(row1[1]));
                fNdata[fwordMap[el.first]] += std::stod(row1[0]);
            }
        }
        data_if.close();
//         std::cout << fHeights[fwordMap[el.first]].size() << "   " << fHeightErrs[fwordMap[el.first]].size() << '\n';
        fNDataTotal += fNdata[fwordMap[el.first]];
    }
    for(auto i = 0; i < 4; ++i)
    {
        std::ifstream data_if;
        data_if.open("../lhcb_bin"+std::to_string(i+1)+".csv");
        auto&& row1 = GParseCSV(data_if);
        auto&& nCols = static_cast<int>(row1.size());
        while(!data_if.eof())
        {
            auto&& row1 = GParseCSV(data_if);
            if(static_cast<int>(row1.size()) == nCols)
            {
                fHeights[100+i].emplace_back(std::stod(row1[0]));
                fHeightErrs[100+i].emplace_back(std::stod(row1[1]));
                fNdata[100+i] += std::stod(row1[0]);
            }
        }
        data_if.close();
    }
}

void GBF2Analytical::GFillMC(const std::string& mc)
{
    fLHCbBins.clear();
    std::ifstream data_if;
    data_if.open(mc);
    auto&& header = GParseCSV(data_if);
    std::vector<std::string> keys;
    auto&& cntr = 0;
    for(auto& el: header)
    {
//         std::cout << el << '\n';
            if(cntr > 5)
            {
                el.erase(el.begin(), el.begin()+4);
                keys.push_back(el);
//                 std::cout << el << '\n';
            }
            ++cntr;
    }
    auto&& nCols = static_cast<int>(header.size());
// //     assert(nCols == 0);
//     std::cout << nCols << '\n';
    for(auto key: keys)
    {
        auto&& nBins = static_cast<int>(fHeightErrs[fwordMap[key]].size());
        if(nBins <= 0)
            continue;
//         std::cout << key << '\n';
        fPhaseSpace[fwordMap[key]].resize(nBins);
        fPDF[fwordMap[key]].resize(nBins);
//         fPDF[fwordMap[key]][2].size();
    }
    fPhaseSpace[0].resize(20);
    fPDF[0].resize(20);
        
    while(!data_if.eof())
    {
        auto&& row1 = GParseCSV(data_if);
        if(row1.size() == nCols)
        {
            auto&& pdf = std::stod(row1[0]);
            auto&& phaseSpace = std::vector<double>{std::stod(row1[1]), std::stod(row1[2]), std::stod(row1[3]), std::stod(row1[4]), std::stod(row1[5])};
            auto&& binsCntr = 6;
            for(auto key: keys)
            {
                if(binsCntr == 10)
                    continue;
//                 std::cout << key << '\n';
//                 std::cout << "************************************************" << '\n';
                auto&& iBin = std::stoi(row1[binsCntr]);
                if(fPhaseSpace[fwordMap[key]].size() > 0)
                {
                    fPhaseSpace[fwordMap[key]][iBin].emplace_back(phaseSpace);
                    fPDF[fwordMap[key]][iBin].emplace_back(pdf);
                }
                ++binsCntr;
            }
            fLHCbBins.emplace_back(std::stoi(row1[binsCntr-1]));
        }
    }
    data_if.close();
}
std::vector<std::string> GBF2Analytical::GParseCSV(std::ifstream& data)
{
    std::vector<std::string>   result;
    std::string                line;
    std::getline(data,line);

    std::stringstream          lineStream(line);
    std::string                cell;

    while(std::getline(lineStream,cell, ','))
    {
        result.push_back(cell);
    }
    // This checks for a trailing comma with no data after it.
//     if (!lineStream && cell.empty())
//     {
//         // If there was a trailing comma then add an empty element.
//         result.push_back("");
//     }
    return result;
}

void GBF2Analytical::GFit(const std::string& filename)
{
    GFillData();
    GFillMC("histogram.csv");
    std::shared_ptr<GGeneratorParser> gp(new GGeneratorParser("../additionalFiles/model.ge"));
    auto&& fRandomize = true;
//     g1->GSetCouplings("../additionalFiles/model.ge");
    auto&& modelStar = gp->GParse();
    auto&& fNVariables = static_cast<int>(modelStar.size());
//     double*p = nullptr;
//     auto&& chi2 = GChi2(&modelStar[0], p);
//     std::cout << chi2 << '\n';
    GSetModelParameters(filename);
    auto&& funcPointer = new TF1("GChi2Total", this, &Gamapola::GBF2Analytical::GChi2Total, 0, 0, 0, "GBF2Analytical", "GLogLikelihoodPDF1");
    TMinuitMinimizer globalMin(ROOT::Minuit::kMigrad);
    ROOT::Math::WrappedMultiTF1 g1(*funcPointer, fNVariables);
    globalMin.SetFunction(g1);
    
    for(int i = 0; i < fNVariables; ++i)
    {
      if(std::strcmp(fkindOfVars[i].c_str(), "fixed") != 0 && fRandomize)
          finVars[i] = flowVars[i] + (fupVars[i] - flowVars[i]) * (double)rand()/RAND_MAX;
      if(std::strcmp(fkindOfVars[i].c_str(), "limited") == 0)
        globalMin.SetLimitedVariable(i, fDictionary[i].c_str(), modelStar[i],fstepsVars[i], flowVars[i], fupVars[i]);
      else if(std::strcmp(fkindOfVars[i].c_str(), "fixed") == 0)
        globalMin.SetFixedVariable(i, fDictionary[i].c_str(), modelStar[i]);
      else
        globalMin.SetVariable(i, fDictionary[i].c_str(), modelStar[i],fstepsVars[i]);
      
//       std::cout << fDictionary[i] << "   !" << fkindOfVars[i] << "!  " << finVars[i] << std::endl;
    }
    
//     globalMin.SetVariable(21, fDictionary[21].c_str(), 0.71,fstepsVars[21]);
//     globalMin.SetVariable(22, fDictionary[22].c_str(), 1.65,fstepsVars[22]);
//     
    globalMin.Minimize();
    globalMin.PrintResults();
}
double GBF2Analytical::GRatio(double* x, double*p)
{
//     std::cout << ea->GPDF(x, p) << "  " << fTempPDF << '\n';
   return ea->GPDF(x, p) / fTempPDF;
}

double GBF2Analytical::GChi2Total(double* x, double* p)
{
    return GChi2CosTheta(x,p)+GChi2(x,p);
}
double GBF2Analytical::GChi2CosTheta(double* x, double* p)
{
    auto&& chi2Total = 0.;
    auto&& cosTheta = fPhaseSpace[0];
    auto&& nBins = static_cast<int>(cosTheta.size());
    std::array<int,4> HiSum = {0};
    std::array<std::vector<double>,4> HiVec;
//     std::cout << nBins << '\n';
    
    for(auto&& iBin = 0; iBin < nBins; ++iBin)
    {
        auto&& phaseSpaceVec = cosTheta[iBin];
        auto&& tempPDFVec = fPDF[0][iBin];
        auto&& nRecords = static_cast<int>(tempPDFVec.size());
        std::array<double,4> nominator = {0.};
        std::array<double,4> denominator = {0.};
        std::array<double,4> Hi = {0.};
        std::array<int,4> nRecDiscr = {0};
        for(auto&& iRecord = 0; iRecord < nRecords; ++iRecord)
        {
            auto&& phaseSpace = phaseSpaceVec[iRecord];
            auto&& mkpipi = std::sqrt(phaseSpace[4]);
//             std::cout << mkpipi << '\n';
            if(1.1 <= mkpipi && mkpipi < 1.3)
            {
                nominator[0] += ea->GPDF(x, &phaseSpace[0]);
                denominator[0] += tempPDFVec[iRecord];
                ++nRecDiscr[0];
            }
            else if(1.3 <= mkpipi && mkpipi < 1.4)
            {
                nominator[1] += ea->GPDF(x, &phaseSpace[0]);
                denominator[1] += tempPDFVec[iRecord];
                ++nRecDiscr[1];
            }
            else if(1.4 <= mkpipi && mkpipi < 1.6)
            {
                nominator[2] += ea->GPDF(x, &phaseSpace[0]);
                denominator[2] += tempPDFVec[iRecord];
                ++nRecDiscr[2];
            }
            else if(1.6 <= mkpipi && mkpipi < 1.9)
            {
                nominator[3] += ea->GPDF(x, &phaseSpace[0]);
                denominator[3] += tempPDFVec[iRecord];
                ++nRecDiscr[3];
            }
//                 nominator += ea->GPDF(x, &phaseSpace[0]);
//                 denominator += tempPDFVec[iRecord];
        }
        for(auto i = 0; i < 4; ++i)
        {
            if(nRecDiscr[i] == 0)
                Hi[i] = 0.;
            else
                Hi[i] = nRecDiscr[i] * nominator[i] / denominator[i];
            
            HiSum[i] += Hi[i];
            HiVec[i].emplace_back(Hi[i]);
        }
    }
    for(auto iHist = 0; iHist < 4; ++iHist)
    {
        auto&& graph1 = std::make_shared<TGraph>();
        auto&& graph2 = std::make_shared<TGraph>();
        auto&& key = 100+iHist;
        auto&& cAlpha = fNdata[key] / HiSum[iHist];
        auto&& chi2 = 0.;
        for(auto&& iBin = 0; iBin < nBins; ++iBin)
        {
//             if( (el.first == 4) && (iBin < 10 || iBin > 50))
//                 continue;
            auto&& comp = (fHeights[key][iBin] - cAlpha * HiVec[iHist][iBin]) * (fHeights[key][iBin] - cAlpha * HiVec[iHist][iBin]) / (fHeightErrs[key][iBin] * fHeightErrs[key][iBin] + cAlpha * cAlpha * HiVec[iHist][iBin]);
            chi2 += comp;
            
            graph1->SetPoint(iBin, iBin, fHeights[key][iBin]);
            graph2->SetPoint(iBin, iBin, cAlpha * HiVec[iHist][iBin]);
        }
        std::cout << key << "  " << chi2 << '\n';
        chi2Total += chi2;
        auto&& c1 = std::make_shared<TCanvas>();
        graph1->SetMarkerStyle(20);
        graph1->SetMarkerColor(kRed);
        graph2->SetMarkerStyle(20);
        graph2->SetMarkerColor(kBlack);
        graph2->Draw("AP");
        graph1->Draw("P");
        c1->SaveAs((std::string("comp")+std::to_string(key)+std::string(".pdf")).c_str());
    }
//     std::cout << HiVec[0].size() << "  " << HiVec[1].size() << "  " << HiVec[2].size() << "  " << HiVec[3].size() << '\n';
    return chi2Total;
}
double GBF2Analytical::GChi2(double* x, double*p)
{
    auto&& chi2Total = 0.;
    for(auto el: fPhaseSpace)
    {
        if(el.first == 0)
        {
//             std::cout << el.second.size() << '\n';
            continue;
        }
        auto&& nBins = static_cast<int>(el.second.size());
        auto&& HiSum = 0.;
        std::vector<double> HiVec;
        for(auto&& iBin = 0; iBin < nBins; ++iBin)
        {
            auto&& phaseSpaceVec = el.second[iBin];
            auto&& tempPDFVec = fPDF[el.first][iBin];
            auto&& nRecords = static_cast<int>(tempPDFVec.size());
            auto&& nominator = 0.;
            auto&& denominator = 0.;
            auto&& Hi = 0.;
            for(auto&& iRecord = 0; iRecord < nRecords; ++iRecord)
            {
                auto&& phaseSpace = phaseSpaceVec[iRecord];
                fTempPDF = tempPDFVec[iRecord];
                nominator += ea->GPDF(x, &phaseSpace[0]);
                denominator += fTempPDF;
//                 Hi += GRatio(x, &phaseSpace[0]);
//                 std::cout << Hi << '\n';
            }
            if(nRecords == 0)
                Hi = 0;
            else
                Hi = nRecords * nominator / denominator;
            HiSum += Hi;
            HiVec.emplace_back(Hi);
//             std::cout << nRecords << " " << nominator << "  " << denominator << "  " << Hi << '\n';
        }
//         std::cout << HiSum << "  " << fNdata[el.first] << '\n';
        auto&& chi2 = 0.;
        auto&& cAlpha = fNdata[el.first] / HiSum;
//         auto&& cBeta = 1. / fNdata[el.first];
        auto&& graph1 = std::make_shared<TGraph>();
        auto&& graph2 = std::make_shared<TGraph>();
        auto&& graph3 = std::make_shared<TGraph>();
        for(auto&& iBin = 0; iBin < nBins; ++iBin)
        {
//             if( (el.first == 4) && (iBin < 10 || iBin > 50))
//                 continue;
            auto&& comp = (fHeights[el.first][iBin] - cAlpha * HiVec[iBin]) * (fHeights[el.first][iBin] - cAlpha * HiVec[iBin]) / (fHeightErrs[el.first][iBin] * fHeightErrs[el.first][iBin] + cAlpha * cAlpha * HiVec[iBin]);
            chi2 += comp;
//             std::cout << fHeights[el.first][iBin] << "  " << HiVec[iBin] << "  " << iBin << '\n';
            graph1->SetPoint(iBin, iBin, fHeights[el.first][iBin]);
            graph2->SetPoint(iBin, iBin, cAlpha * HiVec[iBin]);
            graph3->SetPoint(iBin, iBin, comp);
        }
        auto&& c1 = std::make_shared<TCanvas>();
        graph1->SetMarkerStyle(20);
        graph1->SetMarkerColor(kRed);
        graph2->SetMarkerStyle(20);
        graph2->SetMarkerColor(kBlack);
        graph2->Draw("AP");
        graph1->Draw("P");
        c1->SaveAs((std::string("comp")+std::to_string(el.first)+std::string(".pdf")).c_str());
        auto&& c2 = std::make_shared<TCanvas>();
        graph3->SetMarkerStyle(20);
        graph3->SetMarkerColor(kBlack);
        graph3->Draw("AP");
        c2->SaveAs((std::string("diff")+std::to_string(el.first)+std::string(".pdf")).c_str());
        chi2Total += chi2;
        std::cout << el.first << "  " << chi2 << '\n';
    }
    std::cout << chi2Total << '\n';
    return chi2Total;
}
}

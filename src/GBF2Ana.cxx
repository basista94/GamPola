#include "GBF2Ana.h"
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

GBF2Ana::GBF2Ana(const int& charge, const std::string& modelPars, const std::string& normInts):
fMf(new GMathFunctions()),
fcharge(charge),
fFitter(new GFitter(charge, normInts+".csv")),
fPDF(0),
fvalidKeys({"cosTheta", "phi", "Mkpipi", "Mkpi", "Mpipi", "Mkpi2"}),
fwordMap({{"cosTheta", 0}, {"phi", 1}, {"Mpipi", 3}, {"Mkpi", 2}, {"Mkpipi",4}, {"Mkpi2",5}}),
fNDataTotal(0),
analyticalCoeffs(0),
binsRecs(0),
fChi2TotVec(0),
fParVec(0)
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

GBF2Ana::~GBF2Ana()
{
}
void GBF2Ana::GSetModelParameters(const std::string& fileName)
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
    
void GBF2Ana::GSetDefaultParameters()
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

int GBF2Ana::GGetBin(const double& s, const double& sMin, const double& sMax, const int& nBins)
{
    return static_cast<int>((s - sMin) / (sMax - sMin) * nBins);
}
int GBF2Ana::GGetLHCbBin(const double& mKpipi)
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
void GBF2Ana::GHistogramize(const std::unordered_map<std::string, std::tuple<int, double, double>>& dict, const std::string& infile, const std::string& outfile, const int& nEvts)
{
    auto g1 = std::make_shared<GGenerator>(fcharge);
    g1->GSetCouplings(infile);
    std::shared_ptr<GGeneratorParser> gp(new GGeneratorParser(infile));
    auto&& model = gp->GParse();
    auto&& modelCoeffs = ea->EvaluateModelCoeffs(&model[0]);
    
    g1->GGenerate(nEvts);
    g1->GWriteToFile("gampola_input.root");
    auto&& eg = g1->GGetEvtGen();
    for(auto&& entry = 0; entry < nEvts; ++entry)
    {
        auto&& s = eg->GGetM()[entry];
        auto&& sjk = eg->GGetSjk()[entry];
        auto&& sik = eg->GGetSij()[entry];
        auto&& phi = eg->GGetPhiPDF()[entry];
        auto&& cosTheta = eg->GGetCosThetaPDF()[entry];
        auto sKpi2 = fMf->GSij(s, sik, sjk);
        auto&& phaseSpace = std::vector<double>{cosTheta, phi, sik, sjk, s};
        auto kino = ea->EvaluateKinematicsCoeffs(&phaseSpace[0],&model[0]);
        auto pdfnum = -fMf->GProcessingComputationOfPDF(&phaseSpace[0],&model[0])/*/ea->GDenominator(&model[0])*/;
//        std::complex<double> pdfC;
//         for(auto i = 0; i < kino.size(); ++i)
//         {
//             pdfC += kino[i] * modelCoeffs[i];
//         }
//         std::cout << " CHECK: " << pdfC << "   " << pdf << '\n';
        
        for(auto& el : kino)
            el /= pdfnum;
        
        auto&& pdf = 0.;
        
        for(auto i = 0; i < static_cast<int>(kino.size()); ++i)
            pdf += std::real(kino[i] * modelCoeffs[i]);
        
        phaseSpace[4] = std::sqrt(phaseSpace[4]);
        phaseSpace[3] = std::sqrt(phaseSpace[3]);
        phaseSpace[2] = std::sqrt(phaseSpace[2]);
        phaseSpace.emplace_back(std::sqrt(sKpi2));
//         if(phaseSpace[4] >= 1.6)
//             continue;
        auto&& isBroken = false;
        std::vector<int> binsVec;
        auto&& cosThetaCntr = 0;
        auto&& cntr = 0;
        for(auto el: dict)
        {
            if(el.first == "cosTheta")
                cosThetaCntr = cntr;
            if(std::find(fvalidKeys.begin(), fvalidKeys.end(), el.first) != fvalidKeys.end())
            {
                auto iBin = GGetBin(phaseSpace[fwordMap[el.first]], std::get<1>(el.second), std::get<2>(el.second), std::get<0>(el.second));
                if( (phaseSpace[fwordMap[el.first]] < std::get<1>(el.second)) ||
                    (phaseSpace[fwordMap[el.first]] > std::get<2>(el.second)) ||
                    iBin < 0 || iBin >= std::get<0>(el.second))
                {
                    isBroken = true;
                    break;
                }
                binsVec.emplace_back(iBin);
            }
            ++cntr;
        }
        if(!isBroken)
        {
            auto&& cntr = 0;
            for(auto el: dict)
            {
//                 if(fwordMap[el.first] == 0)
//                     continue;
                if(fwordMap[el.first] == 0)
                {
                    std::vector<double> low = {1.1,1.3,1.4,1.6};
                    std::vector<double> up = {1.3,1.4,1.6,1.9};
                    for(auto i = 0; i < 4; ++i)
                    {
                        if(low[i] <= phaseSpace[4] && phaseSpace[4] < up[i])
                        {
                            auto&& key = 100+i;
                            if(analyticalCoeffs[key][binsVec[cosThetaCntr]].size() == 0)
                            {
                                analyticalCoeffs[key][binsVec[cosThetaCntr]] = kino;
                                binsRecs[key][binsVec[cosThetaCntr]] = 1;
                                fPDF[key][binsVec[cosThetaCntr]] += pdf;
                            }
                            else
                            {
                                ++binsRecs[key][binsVec[cosThetaCntr]];
//                             std::cout << analyticalCoeffs[key][binsVec[cosThetaCntr]].size() << '\n';
                                for(auto i = 0; i < static_cast<int>(kino.size()); ++i)
                                {
                                    analyticalCoeffs[key][binsVec[cosThetaCntr]][i] += kino[i];
                                }
                                fPDF[key][binsVec[cosThetaCntr]] += pdf;
                            }
                        }
                    }
                }
                if(analyticalCoeffs[fwordMap[el.first]][binsVec[cntr]].size() == 0)
                {
                    analyticalCoeffs[fwordMap[el.first]][binsVec[cntr]] = kino;
                    binsRecs[fwordMap[el.first]][binsVec[cntr]] = 1;
                }
                else
                {
//                     std::cout << analyticalCoeffs[fwordMap[el.first]][binsVec[cntr]].size() << '\n';
                    ++binsRecs[fwordMap[el.first]][binsVec[cntr]];
                    for(auto i = 0; i < static_cast<int>(kino.size()); ++i)
                    {
                        analyticalCoeffs[fwordMap[el.first]][binsVec[cntr]][i] += kino[i];
                    }
                }
                fPDF[fwordMap[el.first]][binsVec[cntr]] += pdf;
                ++cntr;
            }
        }
    }
    auto&& lhcb = binsRecs[4];
    auto&& graphRatio = std::make_shared<TGraph>();
        auto&& graphHiaplhaStr = std::make_shared<TGraph>();
        auto&& graphdenom = std::make_shared<TGraph>();
        for(auto&& iBin = 0; iBin < lhcb.size(); ++iBin)
        {
            auto&& denominator = fPDF[4][iBin];
            auto&& Hialpha_str = binsRecs[4][iBin];
            auto&& Ai = 0.;
            if(Hialpha_str == 0)
            {
                graphRatio->SetPoint(iBin, iBin, Hialpha_str);
            }
            else
            {
                graphRatio->SetPoint(iBin, iBin, Hialpha_str / denominator);
            }
            graphHiaplhaStr->SetPoint(iBin, iBin, Hialpha_str);
            graphdenom->SetPoint(iBin, iBin, denominator);
        }
        auto&& cRatio = std::make_shared<TCanvas>();
        graphRatio->SetMarkerStyle(20);
        graphRatio->SetMarkerColor(kRed);
        graphRatio->Draw("AP");
        cRatio->SaveAs("ratio.pdf");
        
        auto&& cHalpha = std::make_shared<TCanvas>();
        graphHiaplhaStr->SetMarkerStyle(20);
        graphHiaplhaStr->SetMarkerColor(kRed);
        graphHiaplhaStr->Draw("AP");
        cHalpha->SaveAs("Hi_alphaStr.pdf");
        
        auto&& cDenom = std::make_shared<TCanvas>();
        graphdenom->SetMarkerStyle(20);
        graphdenom->SetMarkerColor(kRed);
        graphdenom->Draw("AP");
        cDenom->SaveAs("denominator.pdf");
//     
//     auto&& graph1 = std::make_shared<TGraph>();
//     auto&& evtsN = 0.;
//     for(auto i = 0; i < lhcb.size(); ++i)
//     {
//         graph1->SetPoint(i, i, binsRecs[4][i]);
//         evtsN += binsRecs[4][i];
//     }
//     std::cout << "Events number: " << evtsN << '\n';
//     auto&& c1 = std::make_shared<TCanvas>();
//     graph1->SetMarkerStyle(20);
//     graph1->SetMarkerColor(kRed);
//     graph1->Draw("AP");
//     c1->SaveAs("bf2_input.pdf");
//     for(auto el: analyticalCoeffs)
//     {
//         for(auto el2 : el.second)
//         {
//             auto&& kino = el2.second;
//             std::complex<double> pdfC;
//             for(auto ii = 0; ii < kino.size(); ++ii)
//                 pdfC += kino[ii] * modelCoeffs[ii];
//             std::cout << std::real(pdfC) << "  " << fPDF[el.first][el2.first] << '\n';
//         }
//     }
    
}
void GBF2Ana::GSetData(const std::string& data, const std::string& key)
{
    fdataMap[key] = data;
//     auto mkpipiH = GFillData(data);
//     GFillMC(mc);
}
void GBF2Ana::GFillData()
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

void GBF2Ana::GFillMC(const std::string& mc)
{
}
std::vector<std::string> GBF2Ana::GParseCSV(std::ifstream& data)
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

void GBF2Ana::GFit(const std::string& filename)
{
    GFillData();
    GFillMC("histogram.csv");
    std::shared_ptr<GGeneratorParser> gp(new GGeneratorParser("../additionalFiles/model.ge"));
    auto&& fRandomize = true;
//     g1->GSetCouplings("../additionalFiles/model.ge");
    auto&& modelStar = gp->GParse();
    auto&& fNVariables = static_cast<int>(modelStar.size());
    double*p = nullptr;
//     auto&& chi2 = GChi2(&modelStar[0], p);
//     std::cout << chi2 << '\n';
    GSetModelParameters(filename);
    auto&& funcPointer = new TF1("GChi2Total", this, &Gamapola::GBF2Ana::GChi2Total, 0, 0, 0, "GBF2Ana", "GLogLikelihoodPDF1");
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
        globalMin.SetFixedVariable(i, fDictionary[i].c_str(), finVars[i]);
      else
      {
        globalMin.SetVariable(i, fDictionary[i].c_str(), modelStar[i],fstepsVars[i]);
        fVarID = i;
      }
      
//       std::cout << fDictionary[i] << "   !" << fkindOfVars[i] << "!  " << finVars[i] << std::endl;
    }
    
//     globalMin.SetVariable(21, fDictionary[21].c_str(), 0.71,fstepsVars[21]);
//     globalMin.SetVariable(22, fDictionary[22].c_str(), 1.65,fstepsVars[22]);
//     

    globalMin.SetMaxIterations(500000);
    globalMin.Minimize();
    globalMin.PrintResults();
    
//     auto&& nPoints = 10000;
//     for(auto i = 0; i < nPoints; ++i)
//     {
//         auto&& parval = -3 + 6. * i / nPoints;
//         modelStar[fVarID] = parval;
//         fChi2TotVec.push_back(GChi2Total(&modelStar[0], p));
//         fParVec.push_back(modelStar[fVarID]);
//     }
//     
//     auto&& graph1 = std::make_shared<TGraph>();
//     for(auto i = 0; i < fChi2TotVec.size(); ++i)
//     {
//         graph1->SetPoint(i, fParVec[i], fChi2TotVec[i]);
//     }
//     auto&& c1 = std::make_shared<TCanvas>();
//     graph1->SetMarkerStyle(20);
//     graph1->SetMarkerSize(0.3);
//     graph1->SetMarkerColor(kRed);
//     graph1->Draw("AP");
//     c1->SaveAs("dependency.pdf");
    std::vector<double> pars(fNVariables,0);
    for(auto i = 0; i < fNVariables; ++i)
    {
        pars[i]=globalMin.X()[i];
//         std::cout << fDictionary[i] << "  " << globalMin.X()[i] << '\n';
    }
    auto g11 = std::make_shared<GGenerator>(fcharge);
    g11->GSetCouplings(pars);
    g11->GGenerate(10000);
    g11->GWriteToFile("charged_10k.root");

//     modelStar[8] = 0.;
//     modelStar[9] = 0.;
//     GChi2(&modelStar[0], p);
//     
//     auto g11 = std::make_shared<GGenerator>(fcharge);
//     g11->GSetCouplings(modelStar);
//     g11->GGenerate(10000);
//     g11->GWriteToFile("gampola.root");
}
double GBF2Ana::GRatio(double* x, double*p)
{
//     std::cout << ea->GPDF(x, p) << "  " << fTempPDF << '\n';
   return ea->GPDF(x, p) / fTempPDF;
}

double GBF2Ana::GChi2Total(double* x, double* p)
{
    auto&& chi2_curr = GChi2CosTheta(x,p) + GChi2(x,p);
    return chi2_curr;
}
double GBF2Ana::GChi2CosTheta(double* x, double* p)
{
    auto&& chi2Total = 0.;
    auto&& cosTheta = analyticalCoeffs[0];
    auto&& nBins = static_cast<int>(cosTheta.size());
    auto&& modelCoeffs = ea->EvaluateModelCoeffs(x);
    std::array<double,4> HiSum = {0};
    std::array<std::vector<double>,4> HiVec;
//     std::cout << nBins << '\n';
//     
    
    for(auto&& iBin = 0; iBin < nBins; ++iBin)
    {
        std::array<double,4> nominator = {0.};
        std::array<double,4> denominator = {0.};
        std::array<double,4> Hi = {0.};
        std::array<int,4> nRecDiscr = {0};
        
        for(auto i = 0; i < 4; ++i)
        {
            auto&& kinCoeffs = analyticalCoeffs[100+i][iBin];
            auto&& nIntegrals = static_cast<int>(kinCoeffs.size());
            for(auto iInt = 0; iInt < nIntegrals; ++iInt)
            {
                nominator[i] += std::real(kinCoeffs[iInt] * modelCoeffs[iInt]);
            }
            denominator[i] = fPDF[100+i][iBin];
            auto&& nRecords = binsRecs[100+i][iBin];
            if(nRecords == 0)
                Hi[i] = 0;
            else
            {
                Hi[i] = nRecords * nominator[i] / denominator[i]/* / ea->GDenominator(x)*/;
//                 std::cout << ea->GDenominator(x) << '\n';
            }
//             std::cout << "LHCb bin: " << i << " Bin # " << iBin << " Hi: " << Hi[i] << " nominator: " <<  nominator[i] << " denominator: " << denominator[i] << " Hi_alpha*: " << nRecords << '\n';
            HiSum[i] += Hi[i];
            HiVec[i].emplace_back(Hi[i]);
//             std::cout << Hi[i] << " ggggggggggggggggggggggggggggggggggggg " << '\n';
        }
    }
    for(auto iHist = 0; iHist < 4; ++iHist)
    {
//         std::cout << iHist << '\n';
        auto&& graph1 = std::make_shared<TGraph>();
        auto&& graph2 = std::make_shared<TGraph>();
        auto&& graph3 = std::make_shared<TGraph>();
        auto&& key = 100+iHist;
        
//         std::cout << "Ai = " << HiSum[iHist] << " Ndatai = " << fNdata[key] << '\n';
        
        auto&& cAlpha = fNdata[key] / HiSum[iHist];
        auto&& chi2 = 0.;
        for(auto&& iBin = 0; iBin < nBins; ++iBin)
        {
//             if( (el.first == 4) && (iBin < 10 || iBin > 50))
//                 continue;
            auto&& comp = (fHeights[key][iBin] - cAlpha * HiVec[iHist][iBin]) * (fHeights[key][iBin] - cAlpha * HiVec[iHist][iBin]) / (fHeightErrs[key][iBin] * fHeightErrs[key][iBin]/* + cAlpha * cAlpha * HiVec[iHist][iBin]*/);
            chi2 += comp;
            
            graph1->SetPoint(iBin, iBin, fHeights[key][iBin]);
            graph2->SetPoint(iBin, iBin, cAlpha * HiVec[iHist][iBin]);
            graph3->SetPoint(iBin, iBin, comp);
        }
//         std::cout << key << "  " << chi2 << '\n';
        chi2Total += chi2;
        auto&& c1 = std::make_shared<TCanvas>();
        graph1->SetMarkerStyle(20);
        graph1->SetMarkerColor(kRed);
        graph2->SetMarkerStyle(20);
        graph1->SetMinimum(0);
        graph2->SetMinimum(0);
        graph2->SetMarkerColor(kBlack);
        graph2->Draw("AP");
        graph1->Draw("P");
//         std::cout << key << '\n';
        c1->SaveAs((std::string("comp")+std::to_string(key)+std::string(".pdf")).c_str());
        
        auto&& c2 = std::make_shared<TCanvas>();
        graph3->SetMarkerStyle(20);
        graph3->SetMarkerColor(kBlack);
        graph3->Draw("AP");
        c2->SaveAs((std::string("chi2")+std::to_string(100+iHist)+std::string(".pdf")).c_str());
    }
//     std::cout << HiVec[0].size() << "  " << HiVec[1].size() << "  " << HiVec[2].size() << "  " << HiVec[3].size() << '\n';
    return chi2Total;
}
double GBF2Ana::GChi2(double* x, double*p)
{
    auto&& chi2Total = 0.;
    auto&& modelCoeffs = ea->EvaluateModelCoeffs(x); 
    for(auto el: analyticalCoeffs)
    {
        if(el.first == 0 || el.first == 100 || el.first == 101 || el.first == 102 || el.first == 103)
        {
            continue;
        }
        auto&& nBins = static_cast<int>(el.second.size());
        auto&& AiSum = 0.;
        std::vector<double> HiVec;
//         auto&& graphRatio = std::make_shared<TGraph>();
//         auto&& graphHiaplhaStr = std::make_shared<TGraph>();
//         auto&& graphdenom = std::make_shared<TGraph>();
        for(auto&& iBin = 0; iBin < nBins; ++iBin)
        {
// //             analyticalCoeffs
            auto&& kinCoeffs = el.second[iBin];
            auto&& nIntegrals = static_cast<int>(kinCoeffs.size());
            auto&& numerator = 0.;
            for(auto a = 0; a < nIntegrals; ++a)
            {
                numerator += std::real(kinCoeffs[a] * modelCoeffs[a]);
            }
            auto&& denominator = fPDF[el.first][iBin];
            auto&& Hialpha_str = binsRecs[el.first][iBin];
            auto&& Ai = 0.;
            if(Hialpha_str == 0)
            {
                Ai = 0;
//                 graphRatio->SetPoint(iBin, iBin, Hialpha_str);
            }
            else
            {
                Ai = Hialpha_str * numerator / denominator/* / ea->GDenominator(x)*/;
//                 graphRatio->SetPoint(iBin, iBin, Hialpha_str / denominator);
            }
            AiSum += Ai;
            HiVec.emplace_back(Ai);
//             graphHiaplhaStr->SetPoint(iBin, iBin, Hialpha_str);
//             graphdenom->SetPoint(iBin, iBin, denominator);
        }
//         auto&& cRatio = std::make_shared<TCanvas>();
//         graphRatio->SetMarkerStyle(20);
//         graphRatio->SetMarkerColor(kRed);
//         graphRatio->Draw("AP");
//         cRatio->SaveAs("ratio.pdf");
//         
//         auto&& cHalpha = std::make_shared<TCanvas>();
//         graphHiaplhaStr->SetMarkerStyle(20);
//         graphHiaplhaStr->SetMarkerColor(kRed);
//         graphHiaplhaStr->Draw("AP");
//         cHalpha->SaveAs("Hi_alphaStr.pdf");
//         
//         auto&& cDenom = std::make_shared<TCanvas>();
//         graphdenom->SetMarkerStyle(20);
//         graphdenom->SetMarkerColor(kRed);
//         graphdenom->Draw("AP");
//         cDenom->SaveAs("denominator.pdf");
// //         std::cout << HiSum << "  " << fNdata[el.first] << '\n';
        auto&& chi2 = 0.;
        auto&& cAlpha = fNdata[el.first] / AiSum;
//         auto&& cBeta = 1. / fNdata[el.first];
        auto&& graph1 = std::make_shared<TGraph>();
        auto&& graph2 = std::make_shared<TGraph>();
        auto&& graph3 = std::make_shared<TGraph>();
        for(auto&& iBin = 0; iBin < nBins; ++iBin)
        {
//             if( (el.first == 4) && (iBin > 50))
//                 continue;
            auto&& comp = (fHeights[el.first][iBin] - cAlpha * HiVec[iBin]) * (fHeights[el.first][iBin] - cAlpha * HiVec[iBin]) / (fHeightErrs[el.first][iBin] * fHeightErrs[el.first][iBin]/* + cAlpha * cAlpha * HiVec[iBin]*/);
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
        graph1->SetMinimum(0);
        graph2->SetMinimum(0);
        graph2->Draw("AP");
        graph1->Draw("P");
//         std::cout << std::string("comp")+std::to_string(el.first)+std::string(".pdf") << '\n';
        
        c1->SaveAs((std::string("comp")+std::to_string(el.first)+std::string(".pdf")).c_str());
        auto&& c2 = std::make_shared<TCanvas>();
        graph3->SetMarkerStyle(20);
        graph3->SetMarkerColor(kBlack);
        graph3->Draw("AP");
        c2->SaveAs((std::string("diff")+std::to_string(el.first)+std::string(".pdf")).c_str());
        chi2Total += chi2;
//         std::cout << el.first << "  " << chi2 << '\n';
    }
//     std::cout << chi2Total << '\n';
    return chi2Total;
}
}

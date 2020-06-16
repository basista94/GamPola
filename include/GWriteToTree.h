#ifndef GWRITETOTREE_H
#define GWRITETOTREE_H

#include <iostream>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TDirectory.h>

namespace Gamapola
{
  class GWriteToTree{
  public:
    GWriteToTree();
    virtual ~GWriteToTree();
    void GCreateFile(const char* filename);
    void GAddHistogram(const char* histName, const std::vector<double>& content, 
                       const double& lowValue, const double& upValue, const int& nBins);
    void GAddHistogram2D(const char* histName, const std::vector<double>& content1, 
                       const double& lowValue1, const double& upValue1, const int& nBins1,
                       const std::vector<double>& content2, 
                       const double& lowValue2, const double& upValue2, const int& nBins2);
    
    void GWriteHistograms(const int& closeFile = 1);
    void GMakeDirectory(const char* dirName);
    void GChangeDirectory(const char* dirName);
  private:
    int fArg;
    TFile * ffile ;
    std::vector<const double*> fAddressesOfPars;
    std::vector<const double*> fAddressesOfPars2D[2];
    std::vector<TH1D*> fHistograms;
    std::vector<TH2D*> fHistograms2D;
    std::vector<unsigned int> fSizeOfParsVectors;
    std::vector<unsigned int> fSizeOfParsVectors2D;
    TDirectory* fCurrentDirectory;
  };
}
#endif

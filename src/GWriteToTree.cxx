#include "GWriteToTree.h"

namespace Gamapola{
  GWriteToTree::GWriteToTree():
  fArg(0),
  ffile(0),
  fAddressesOfPars(0),
  fHistograms(0),
  fSizeOfParsVectors(0),
  fCurrentDirectory(0)
  {
    std::cout << "GWriteToTree constructor calling. . ." << std::endl;
  }
  GWriteToTree::~GWriteToTree()
  {
    ffile->Close();
    ffile = NULL;
    fCurrentDirectory = NULL;
    for(std::vector<const double*>::iterator it = fAddressesOfPars.begin(); it != fAddressesOfPars.end(); ++it)
      delete *it;
    for(std::vector<const double*>::iterator it = fAddressesOfPars2D[0].begin(); it != fAddressesOfPars2D[0].end(); ++it)
      delete *it; 
    for(std::vector<const double*>::iterator it = fAddressesOfPars2D[1].begin(); it != fAddressesOfPars2D[1].end(); ++it)
      delete *it;
    for(std::vector<TH1D*>::iterator it = fHistograms.begin(); it != fHistograms.end(); ++it)
      delete *it;
    for(std::vector<TH2D*>::iterator it = fHistograms2D.begin(); it != fHistograms2D.end(); ++it)
      delete *it;
    delete ffile;
    delete fCurrentDirectory;
//     std::cout << "GWriteToTree destructor calling. . ." << std::endl;
  }
  void GWriteToTree::GCreateFile(const char* filename)
  {
    ffile = new TFile( filename, "recreate") ;
    fCurrentDirectory = gDirectory;
  }
  void GWriteToTree::GAddHistogram(const char* histName, const std::vector<double>& content, 
                       const double& lowValue, const double& upValue, const int& nBins)
  {
    TH1D* h = new TH1D(histName, histName, nBins, lowValue, upValue);
    auto&& temp = content;
    fAddressesOfPars.push_back(&temp[0]);
//     std::cout << fAddressesOfPars.size()-1 << '\n';
//     if(histName == "cosTheta")
//     for(int i = 0; i < content.size(); i++)
//     {
//       std::cout << fAddressesOfPars[fAddressesOfPars.size()-1][i] << std::endl;
//     }
    fHistograms.push_back(h);
    fSizeOfParsVectors.push_back(content.size());
  }
  void GWriteToTree::GAddHistogram2D(const char* histName, const std::vector<double>& contentx, 
                       const double& lowValuex, const double& upValuex, const int& nBinsx,
                       const std::vector<double>& contenty, 
                       const double& lowValuey, const double& upValuey, const int& nBinsy)
  {
    TH2D* h2d = new TH2D(histName, histName, nBinsx, lowValuex, upValuex, nBinsy, lowValuey, upValuey);
    fAddressesOfPars2D[0].push_back(&contentx[0]);
    fAddressesOfPars2D[1].push_back(&contenty[0]);
    fHistograms2D.push_back(h2d);
    fSizeOfParsVectors2D.push_back(contentx.size());
  }
  void GWriteToTree::GWriteHistograms(const int& closeFile)
  {
//     for(int i = 0; i < fSizeOfParsVectors[1]; i++)
//     {
//       std::cout << fAddressesOfPars[1][i] << std::endl;
//     }
//     std::cout << fAddressesOfPars.size() << '\n';
    for(unsigned int j = 0; j < fHistograms.size(); j++)
    {
        
      for(unsigned int i = 0; i < fSizeOfParsVectors[j]; i++)
      {
        fHistograms[j]->Fill(fAddressesOfPars[j][i]);
//         if(j == 0)
//             std::cout << fAddressesOfPars[j][i] << std::endl;
      }
      fHistograms[j]->Write();
    }
    
    for(unsigned int j = 0; j < fHistograms2D.size(); j++)
    {
      for(unsigned int iP = 0; iP < fSizeOfParsVectors2D[j]; iP++)
        fHistograms2D[j]->Fill(fAddressesOfPars2D[0][j][iP], fAddressesOfPars2D[1][j][iP]);
      fHistograms2D[j]->Write();
    }
    fHistograms.clear();
    fHistograms2D.clear();
    fAddressesOfPars.clear();
    fAddressesOfPars2D[0].clear();
    fAddressesOfPars2D[1].clear();
    fSizeOfParsVectors.clear();
    fSizeOfParsVectors2D.clear();
    if(closeFile == 1)
        ffile->Close();
  }
  void GWriteToTree::GMakeDirectory(const char* dirName)
  {
    gDirectory->mkdir(dirName);
  }
  void GWriteToTree::GChangeDirectory(const char* dirName)
  {
    gDirectory->cd(dirName);
  }
}

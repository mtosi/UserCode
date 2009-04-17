// Class: VBFHZZllbbValidationPlots.h
// Description:  Some Basic validation plots for jets.
// Author: K. Kousouris
// Date:  27 - August - 2008
//
#ifndef VBFHZZllbbValidationPlots_h
#define VBFHZZllbbValidationPlots_h
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TFile.h>
#include "TNamed.h"
#include <vector>
#include <map>
#include "FWCore/Framework/interface/EDAnalyzer.h"

class VBFHZZllbbValidationPlots : public edm::EDAnalyzer 
   {
     public:
       VBFHZZllbbValidationPlots(edm::ParameterSet const& cfg);
     private:
       void beginJob(edm::EventSetup const& iSetup);
       void analyze(edm::Event const& e, edm::EventSetup const& iSetup);
       void endJob();

       void fillHist1D (const TString& histName, const Double_t& x);
       void fillHist2D (const TString& histName, const Double_t& x, const Double_t& y);
       void fillProfile(const TString& histName, const Double_t& x, const Double_t& y);

       std::string outputFileName_; 
       double      minPtCut_;
       double      dRmatchCut_;
       bool        MC_;

       std::map<TString, TH1*> histNames1DMap_;  
       std::map<TString, TH2*> histNames2DMap_;
       std::map<TString, TProfile*> profileNamesMap_; 

       TFile* outputFile_;
   

  };

#endif

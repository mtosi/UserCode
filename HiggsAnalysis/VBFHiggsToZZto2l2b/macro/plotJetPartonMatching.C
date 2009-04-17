#include <string.h>
#include <sstream>
#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
#include "TLine.h"
#include "Riostream.h"

#include "AnalysisExamples/AnalysisClasses/interface/mergeObj.h"
#include "AnalysisExamples/AnalysisClasses/interface/findObj.h"
#include "AnalysisExamples/AnalysisClasses/interface/THStackLegend.h"
#include "AnalysisExamples/AnalysisClasses/interface/ListFashionAttributedHisto.h"
#include "AnalysisExamples/AnalysisObjects/interface/xSecLO.h"
#include "AnalysisExamples/AnalysisObjects/interface/xSecNLO.h"
#include "AnalysisExamples/AnalysisObjects/interface/BR.h"

using namespace anaobj;

void plotJetPartonMatching(TString outputFileName = "plotJetPartonMatchingOutput.root") {

  // general root setting
  gROOT->Reset(); 
  //  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.1);
  gStyle->SetLabelSize(0.04,"X");
  gStyle->SetLabelSize(0.04,"Y");
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);

  // ************************************************************
  std::vector<TString> sgnFileSuffix;
  //  sgnFileSuffix.push_back("VBFH160ZZ");
  sgnFileSuffix.push_back("VBFH200ZZ");
  //  sgnFileSuffix.push_back("VBFH400ZZ");
  //  sgnFileSuffix.push_back("VBFH800ZZ");
  // ************************************************************
  // List of Files
  TList * sgnFileList = new TList();
  std::vector<TString>::const_iterator sgnFileSuffix_itr = sgnFileSuffix.begin();
  for ( ; sgnFileSuffix_itr != sgnFileSuffix.end(); ++sgnFileSuffix_itr )
    sgnFileList->Add( TFile::Open("jetPartonMatch"+*sgnFileSuffix_itr+"output.root"  ) );
    // ************************************************************
  std::vector<double> sgnCrossSections;
  sgnCrossSections.push_back(1000*xSec_VBFH160_*BR_H160ZZ_);
  sgnCrossSections.push_back(1000*xSec_VBFH200_*BR_H200ZZ_);
  sgnCrossSections.push_back(1000*xSec_VBFH400_*BR_H400ZZ_);
  sgnCrossSections.push_back(1000*xSec_VBFH800_*BR_H800ZZ_);
  // ************************************************************
  std::vector<TString> varNameVector;
  varNameVector.push_back("Bothok_05");
  varNameVector.push_back("Bothok_05_20pc");
  varNameVector.push_back("Bothok_02");
  varNameVector.push_back("Bothok_02_20pc");
  // ************************************************************

  int indexSample = 0;
  std::vector<TH1*> sgnEventsNumberTH1vector;
  if (sgnFileSuffix.size() != 0) 
    findObj( sgnEventsNumberTH1vector, 
	     sgnFileList, 
	     "Allevents", 
	     &sgnFileSuffix
	     );
  
  int nbins_   = sgnEventsNumberTH1vector[0]->GetNbinsX();
  double xmin_ = sgnEventsNumberTH1vector[0]->GetXaxis()->GetXmin();
  double xmax_ = sgnEventsNumberTH1vector[0]->GetXaxis()->GetXmax();

  TFile * outputFile = TFile::Open( outputFileName, "RECREATE" );
  std::vector<TString>::const_iterator varNameVector_itr = varNameVector.begin();

  for ( ; varNameVector_itr != varNameVector.end(); ++varNameVector_itr ) {
    std::vector<TH1*> sgnTH1vector;
    
    TString varName = *varNameVector_itr;
    std::cout << "varName: " << varName << std::endl;
    findObj( sgnTH1vector, 
	     sgnFileList, 
	     *varNameVector_itr, 
	     &sgnFileSuffix
	     );

    outputFile->cd();
    TDirectory * directory = outputFile->mkdir( *varNameVector_itr, *varNameVector_itr );
    directory->cd();

    std::vector<TH1*>::const_iterator sgnTH1vector_itr = sgnTH1vector.begin();
    for ( ; sgnTH1vector_itr != sgnTH1vector.end();++sgnTH1vector_itr) {

      TString histoName  = (*sgnTH1vector_itr)->GetName();
      TString histoTitle = (*sgnTH1vector_itr)->GetTitle();

      TH1D* fracHisto = new TH1D(histoName+"frac",histoTitle,nbins_,xmin_,xmax_);
      fracHisto->Sumw2();
      fracHisto->Divide(sgnTH1vector[0],sgnEventsNumberTH1vector[0]);
  
      TCanvas * canvas = new TCanvas(*varNameVector_itr,*varNameVector_itr,1000,800);
      canvas->cd();
      fracHisto->Draw();

      fracHisto->Write();
      canvas->Write();
      canvas->Print(*varNameVector_itr+"JetPartonMatchingPlot.jpg");
      
    }
  }
  
}


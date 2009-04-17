#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TMath.h"
#include "TLegend.h"
#include "TString.h"
#include "TPaveLabel.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>


void mHtoZZplot( TString inputFileName1,
		 TString inputFileName2,
		 TString inputFileName3,
		 TString inputFileName4) {
  // general root setting
  gROOT->Reset(); 
  //  gROOT->SetBatch(kTRUE);
  //  gStyle->SetOptStat(0);
  gStyle->SetTitleSize(0.1);
  gStyle->SetLabelSize(0.04,"X");
  gStyle->SetLabelSize(0.04,"Y");
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);

  const double xSec_ggH160_  = 21.0; // in pb	  
  const double xSecPosUncert_ggH160_ =  18.; // %   
  const double xSecNegUncert_ggH160_ = -14.; // %	  
  const double xSec_VBFH160_ = 3.32; // in pb	  
  const double xSecPosUncert_VBFH160_ =  0.5; // % 
  const double xSecNegUncert_VBFH160_ = -1.0; // % 
  const double BR_H160ZZ_ = 0.4334;                 

  const double xSec_ggH200_  = 14.8; // in pb	  
  const double xSecPosUncert_ggH200_ =  17.; // %   
  const double xSecNegUncert_ggH200_ = -14.; // %	  
  const double xSec_VBFH200_ = 2.53; // in pb	  
  const double xSecPosUncert_VBFH200_ =  0.4; // % 
  const double xSecNegUncert_VBFH200_ = -1.1; // % 
  const double BR_H200ZZ_ = 0.2613;                 

  const double xSec_ggH800_  = 0.397; // in pb	  
  const double xSecPosUncert_ggH800_ =  16.; // %   
  const double xSecNegUncert_ggH800_ = -14.; // %	  
  const double xSec_VBFH800_ = 0.196; // in pb	  
  const double xSecPosUncert_VBFH800_ =  0.7; // % 
  const double xSecNegUncert_VBFH800_ = -1.5; // % 
  const double BR_H800ZZ_ = 0.2928;                 

  TFile * inputFile1 = new TFile(inputFileName1);
  TFile * inputFile2 = new TFile(inputFileName2);
  TFile * inputFile3 = new TFile(inputFileName3);
  TFile * inputFile4 = new TFile(inputFileName4);
 
  //    TString HistoName("MHtoZZ");
    TString HistoName("recHlljj_AFTERZdileptonMass");
  TH1D * histo1 = dynamic_cast<TH1D*>(inputFile1->Get(HistoName));
  TString title = histo1->GetTitle();
  histo1->SetTitle(title+ " for m_{H} = 160GeV/c^{2}");
  histo1->GetXaxis()->SetTitle("m_{H} GeV/c^{2}");
  TH1D * histo2 = dynamic_cast<TH1D*>(inputFile2->Get(HistoName));
  histo2->SetTitle(title+ " for m_{H} = 200GeV");
  histo2->GetXaxis()->SetTitle("m_{H} GeV/c^{2}");
  TH1D * histo3 = dynamic_cast<TH1D*>(inputFile3->Get(HistoName));
  histo3->SetTitle(title+ " for m_{H} = 400GeV");
  histo3->GetXaxis()->SetTitle("m_{H} GeV/c^{2}");
  TH1D * histo4 = dynamic_cast<TH1D*>(inputFile4->Get(HistoName));
  histo4->SetTitle(title+ " for m_{H} = 800GeV");
  histo4->GetXaxis()->SetTitle("m_{H} GeV/c^{2}");

  TCanvas * canvas = new TCanvas ( "diffxSec", "differential xSec", 1200, 800 );
  gStyle->SetOptStat(0);
  canvas->UseCurrentStyle();
  canvas->Divide(2,2);
  canvas->cd(1);
  histo1->SetMarkerColor(215);
  histo1->Draw();
  canvas->cd(2);
  histo2->SetMarkerColor(215);
  histo2->Draw();
  canvas->cd(3);
  histo3->SetMarkerColor(215);
  histo3->Draw();
  canvas->cd(4);
  histo4->SetMarkerColor(215);
  histo4->Draw();

}

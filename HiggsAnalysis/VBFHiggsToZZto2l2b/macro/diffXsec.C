#include "TH1.h"
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


void diffXsec( TString inputFile1Name,
	       TString inputFile2Name,
	       TString inputFile3Name,
	       TString inputFile4Name,
	       TString eventsNumberHistoName,
	       TString HistoName) {
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

  const double xSec_ggH400_  = ; // in pb	  
  const double xSecPosUncert_ggH400_ =  ; // %   
  const double xSecNegUncert_ggH400_ = ; // %	  
  const double xSec_VBFH400_ = ; // in pb	  
  const double xSecPosUncert_VBFH400_ =  // % 
  const double xSecNegUncert_VBFH400_ =  // % 
  const double BR_H400ZZ_ =                  

  const double xSec_ggH800_  = 0.397; // in pb	  
  const double xSecPosUncert_ggH800_ =  16.; // %   
  const double xSecNegUncert_ggH800_ = -14.; // %	  
  const double xSec_VBFH800_ = 0.196; // in pb	  
  const double xSecPosUncert_VBFH800_ =  0.7; // % 
  const double xSecNegUncert_VBFH800_ = -1.5; // % 
  const double BR_H800ZZ_ = 0.2928;                 

  TFile * inputFile1 = new TFile(inputFile1Name);
  TFile * inputFile2 = new TFile(inputFile2Name);
  TFile * inputFile3 = new TFile(inputFile3Name);
  TFile * inputFile4 = new TFile(inputFile4Name);
 
  TH1D * eventsNumberHisto1 = dynamic_cast<TH1D*>(inputFile1->Get(eventsNumberHistoName));
  TH1D * histo1 = dynamic_cast<TH1D*>(inputFile1->Get(HistoName));
  TH1D * eventsNumberHisto2 = dynamic_cast<TH1D*>(inputFile2->Get(eventsNumberHistoName));
  TH1D * histo2 = dynamic_cast<TH1D*>(inputFile2->Get(HistoName));
  TH1D * eventsNumberHisto3 = dynamic_cast<TH1D*>(inputFile3->Get(eventsNumberHistoName));
  TH1D * histo3 = dynamic_cast<TH1D*>(inputFile3->Get(HistoName));
  TH1D * eventsNumberHisto4 = dynamic_cast<TH1D*>(inputFile4->Get(eventsNumberHistoName));
  TH1D * histo4 = dynamic_cast<TH1D*>(inputFile4->Get(HistoName));

  TString HistoTitle = histo1->GetTitle();
  histo1->SetTitle(HistoTitle+ " for m_{H} = 160GeV");
  if(histo1->GetName()->Contains("bbbar"))         histo1->GetXaxis()->SetTitle("m_{bb} (GeV/c^{2})");
  else if(histo1->GetName()->Contains("muon"))     histo1->GetXaxis()->SetTitle("m_{#mu#mu} (GeV/c^{2})");
  else if(histo1->GetName()->Contains("electron")) histo1->GetXaxis()->SetTitle("m_{ee} (GeV/c^{2})");

  histo2->SetTitle(HistoTitle+ " for m_{H} = 200GeV");
  histo2->GetXaxis()->SetTitle("m_{bb} (GeV/c^{2})");
  histo3->SetTitle(HistoTitle+ " for m_{H} = 400GeV");
  histo3->GetXaxis()->SetTitle("m_{bb} (GeV/c^{2})");
  histo4->SetTitle(HistoTitle+ " for m_{H} = 800GeV");
  histo4->GetXaxis()->SetTitle("m_{bb} (GeV/c^{2})");


  double events1 = (double)eventsNumberHisto1->GetMaximum(); 
  double events2 = (double)eventsNumberHisto2->GetMaximum();
  double events3 = (double)eventsNumberHisto3->GetMaximum();
  double events4 = (double)eventsNumberHisto4->GetMaximum();
  double xSecEff1 = (xSec_ggH160_+xSec_VBFH160_)*BR_H160ZZ_;
  double xSecEff2 = (xSec_ggH200_+xSec_VBFH200_)*BR_H200ZZ_;
  double xSecEff3 = (xSec_ggH400_+xSec_VBFH400_)*BR_H400ZZ_; 
  double xSecEff4 = (xSec_ggH800_+xSec_VBFH800_)*BR_H800ZZ_; 

  double factor1 = xSecEff1/events1;
  double factor2 = xSecEff2/events2;
  double factor3 = xSecEff3/events3;
  double factor4 = xSecEff4/events4;
  std::cout << "events1: "  << events1  << " events2: "  << events2   << " events3: "  << events3  << " events4: "  << events4  << std::endl; 
  std::cout << "xSecEff1: " << xSecEff1 << " xSecEff2: " << xSecEff2  << " xSecEff3: " << xSecEff3 << " xSecEff4: " << xSecEff4 << std::endl;
  std::cout << "factor1: "  << factor1  << " factor2: "  << factor2   << " factor3: "  << factor3  << " factor4: "  << factor4  << std::endl;

  TLegend * legend1 = new TLegend(0.5,0.6,0.99,0.99);
  legend1->SetFillColor(0);
  legend1->SetBorderSize(0);  
  legend1->SetTextFont(72);
  legend1->SetTextSize(0.027);
  legend1->SetFillColor(0);
  TLegend * legend2 = new TLegend(0.5,0.6,0.99,0.99);
  legend2->SetFillColor(0);
  legend2->SetBorderSize(0);  
  legend2->SetTextFont(72);
  legend2->SetTextSize(0.027);
  legend2->SetFillColor(0);
  TLegend * legend3 = new TLegend(0.5,0.6,0.99,0.99);
  legend3->SetFillColor(0);
  legend3->SetBorderSize(0);  
  legend3->SetTextFont(72);
  legend3->SetTextSize(0.027);
  legend3->SetFillColor(0);
  TLegend * legend4 = new TLegend(0.5,0.6,0.99,0.99);
  legend4->SetFillColor(0);
  legend4->SetBorderSize(0);  
  legend4->SetTextFont(72);
  legend4->SetTextSize(0.027);
  legend4->SetFillColor(0);
  char nev[50];

  TCanvas * canvas = new TCanvas ( "diffxSec", "differential xSec", 1200, 800 );
  gStyle->SetOptStat(0);
  canvas->UseCurrentStyle();
  canvas->Divide(2,2);
  canvas->cd(1);
  histo1->Scale(factor1);
  histo1->Draw();
  sprintf(nev,"effective luminosity:%.2f",factor1);
  legend1->AddEntry(histo1,nev,"");
  legend1->Draw();
  canvas->cd(2);
  histo2->Scale(factor2);
  histo2->Draw();
  legend2->Clear();
  sprintf(nev,"effective luminosity:%.2f",factor2);
  legend2->AddEntry(histo2,nev,"");
  legend2->Draw();
  canvas->cd(3);
  histo3->Scale(factor3);
  histo3->Draw();
  legend3->Clear();
  sprintf(nev,"effective luminosity:%.2f",factor3);
  legend3->AddEntry(histo3,nev,"");
  legend3->Draw();
  canvas->cd(4);
  histo4->Scale(factor4);
  histo4->Draw();
  legend4->Clear();
  sprintf(nev,"effective luminosity:%.2f",factor4);
  legend4->AddEntry(histo4,nev,"");
  legend4->Draw();

}

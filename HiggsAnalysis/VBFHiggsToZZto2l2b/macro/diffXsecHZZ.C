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


void diffXsecHZZ( TString inputFile1Name,
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

  const double luminosityFactor = 100000.; // 100fb^-1 of luminosity assumed in histograms

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

  const double xSec_ggH400_  = 7.88; // in pb	  
  const double xSecPosUncert_ggH400_ =  17.0; // %   
  const double xSecNegUncert_ggH400_ = -14.0; // %	  
  const double xSec_VBFH400_ = 0.869; // in pb	  
  const double xSecPosUncert_VBFH400_ = -0.2; // % 
  const double xSecNegUncert_VBFH400_ = -1.1; // % 
  const double BR_H400ZZ_ = 0.2742;                

  const double xSec_ggH800_  = 0.397; // in pb	  
  const double xSecPosUncert_ggH800_ =  16.; // %   
  const double xSecNegUncert_ggH800_ = -14.; // %	  
  const double xSec_VBFH800_ = 0.196; // in pb	  
  const double xSecPosUncert_VBFH800_ =  0.7; // % 
  const double xSecNegUncert_VBFH800_ = -1.5; // % 
  const double BR_H800ZZ_ = 0.2928;                 

  const double idiotFactor = 2*0.01513*0.03363;

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
  int nbin = histo1->GetNbinsX();
  double xmin = histo1->GetXaxis()->GetXmin();
  double xmax = histo1->GetXaxis()->GetXmax();
  TH1D * histoTot = new TH1D(HistoName+"TOT","inclusive "+HistoTitle,nbin,xmin,xmax);

  histo1->SetTitle(HistoTitle+ " for m_{H} = 160GeV");
  histo2->SetTitle(HistoTitle+ " for m_{H} = 200GeV");
  histo3->SetTitle(HistoTitle+ " for m_{H} = 400GeV");
  histo4->SetTitle(HistoTitle+ " for m_{H} = 800GeV");

  double histo1Entries = histo1->GetEntries();
  double histo2Entries = histo2->GetEntries();
  double histo3Entries = histo3->GetEntries();
  double histo4Entries = histo4->GetEntries();
  double events1 = (double)eventsNumberHisto1->GetMaximum(); 
  double events2 = (double)eventsNumberHisto2->GetMaximum();
  double events3 = (double)eventsNumberHisto3->GetMaximum();
  double events4 = (double)eventsNumberHisto4->GetMaximum();
  double xSecEff1 = (xSec_ggH160_+xSec_VBFH160_)*BR_H160ZZ_;
  double xSecEff2 = (xSec_ggH200_+xSec_VBFH200_)*BR_H200ZZ_;
  double xSecEff3 = (xSec_ggH400_+xSec_VBFH400_)*BR_H400ZZ_; 
  double xSecEff4 = (xSec_ggH800_+xSec_VBFH800_)*BR_H800ZZ_; 
  double invLuminosityEff1 = xSecEff1/events1;
  double invLuminosityEff2 = xSecEff2/events2;
  double invLuminosityEff3 = xSecEff3/events3;
  double invLuminosityEff4 = xSecEff4/events4;

  TLegend * legend1 = new TLegend(0.5,0.6,0.89,0.8);
  legend1->SetFillColor(0);
  legend1->SetBorderSize(0);  
  legend1->SetTextFont(72);
  legend1->SetTextSize(0.035);
  legend1->SetFillColor(0);
  TLegend * legend2 = new TLegend(0.5,0.6,0.89,0.8);
  legend2->SetFillColor(0);
  legend2->SetBorderSize(0);  
  legend2->SetTextFont(72);
  legend2->SetTextSize(0.035);
  legend2->SetFillColor(0);
  TLegend * legend3 = new TLegend(0.5,0.6,0.89,0.8);
  legend3->SetFillColor(0);
  legend3->SetBorderSize(0);  
  legend3->SetTextFont(72);
  legend3->SetTextSize(0.035);
  legend3->SetFillColor(0);
  TLegend * legend4 = new TLegend(0.5,0.6,0.89,0.8);
  legend4->SetFillColor(0);
  legend4->SetBorderSize(0);  
  legend4->SetTextFont(72);
  legend4->SetTextSize(0.035);
  legend4->SetFillColor(0);
  TLegend * legendTot = new TLegend(0.7,0.65,0.89,0.85);
  legendTot->SetFillColor(0);
  legendTot->SetBorderSize(0);  
  legendTot->SetTextFont(72);
  legendTot->SetTextSize(0.035);
  legendTot->SetFillColor(0);
  char nev[50];

  TCanvas * canvas = new TCanvas ( "diffxSec", "differential xSec", 1200, 800 );
  gStyle->SetOptStat(0);
  //  canvas->UseCurrentStyle();
  canvas->Divide(2,2);
  canvas->cd(1);
  canvas->SetLogy();
  //histo1->Sumw2();
  histo1->SetNormFactor(double(histo1Entries*invLuminosityEff1*luminosityFactor*idiotFactor));
  histo1->Draw("histo");
  sprintf(nev,"number of events in 100 fb^{-1}:");
  legend1->AddEntry(histo1,nev,"");
  sprintf(nev,"%.2f",histo1Entries*invLuminosityEff1*luminosityFactor*idiotFactor);
  legend1->AddEntry(histo1,nev,"");
  legend1->Draw();
  canvas->cd(2);
  histo2->SetNormFactor(double(histo2Entries*invLuminosityEff2*luminosityFactor*idiotFactor));
  histo2->Draw();
  legend2->Clear();
  sprintf(nev,"number of events in 100 fb^{-1}:");
  legend2->AddEntry(histo2,nev,"");
  sprintf(nev,"%.2f",histo2Entries*invLuminosityEff2*luminosityFactor*idiotFactor);
  legend2->AddEntry(histo2,nev,"");
  legend2->Draw();
  canvas->cd(3);
  histo3->SetNormFactor(double(histo3Entries*invLuminosityEff3*luminosityFactor*idiotFactor));
  histo3->Draw();
  legend3->Clear();
  sprintf(nev,"number of events in 100 fb^{-1}:");
  legend3->AddEntry(histo3,nev,"");
  sprintf(nev,"%.2f",histo3Entries*invLuminosityEff3*luminosityFactor*idiotFactor);
  legend3->AddEntry(histo3,nev,"");
  legend3->Draw();
  canvas->cd(4);
  histo4->SetNormFactor(double(histo4Entries*idiotFactor*invLuminosityEff4*luminosityFactor));
  histo4->Draw();
  legend4->Clear();
  sprintf(nev,"number of events in 100 fb^{-1}:");
  legend4->AddEntry(histo4,nev,"");
  std::cout << "integral: " << histo4Entries*invLuminosityEff4*luminosityFactor*idiotFactor << std::endl;
  sprintf(nev,"%.2f",histo4Entries*invLuminosityEff4*luminosityFactor*idiotFactor);
  legend4->AddEntry(histo4,nev,"");
  legend4->Draw();

  histoTot->Sumw2();
  histoTot->Add(histo1);
  histoTot->Add(histo2);
  histoTot->Add(histo3);
  histoTot->Add(histo4);

  TCanvas * canvasTot = new TCanvas ( "inclusivediffxSec", "differential xSec", 1200, 400 );
  gStyle->SetOptStat(0);
  //  canvasTot->UseCurrentStyle();
  canvasTot->cd();
  gPad->SetLogy();
  histo1->SetTitle("ZZ (Z->ll,Z->b\\bar{b}) mass (100 fb^{-1})");
  histo1->GetXaxis()->SetTitle("m_{ZZ} (GeV)");
  histo1->GetYaxis()->SetTitle("d#sigma/dm_{ZZ} events/10GeV");
  histo1->SetMinimum(0.001);
  histo1->SetLineColor(70);
  histo1->SetFillColor(70);
  histo1->Draw();
  histo2->SetLineColor(76);
  histo2->SetFillColor(76);
  histo2->Draw("same");
  histo3->SetLineColor(82);
  histo3->SetFillColor(82);
  histo3->Draw("same");
  histo4->SetLineColor(86);
  histo4->SetFillColor(86);
  histo4->Draw("same");
  histoTot->Draw("same");

  sprintf(nev,"m_{H}=160 GeV");
  legendTot->AddEntry(histo1,nev,"l");
  sprintf(nev,"m_{H}=200 GeV");
  legendTot->AddEntry(histo2,nev,"l");
  sprintf(nev,"m_{H}=400 GeV");
  legendTot->AddEntry(histo3,nev,"l");
  sprintf(nev,"m_{H}=800 GeV");
  legendTot->AddEntry(histo4,nev,"l");
  sprintf(nev,"inclusive");
  legendTot->AddEntry(histoTot,nev,"lp");
  legendTot->Draw();

  //  TFile * outputFile = new TFile("diffXsecHZZ.root","RECREATE");
  TFile * outputFile = new TFile("output.root","RECREATE");
  histo1->SetName(HistoName+"160");
  histo2->SetName(HistoName+"200");
  histo3->SetName(HistoName+"400");
  histo4->SetName(HistoName+"800");
  histo1->Write();
  histo2->Write();
  histo3->Write();
  histo4->Write();
  histoTot->Write();

  //  canvas->Print("HZZ.jpg");
  //  canvasTot->Print("HZZtot.jpg");

  outputFile->Close();
}

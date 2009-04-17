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


void diffXsecZbb( TString inputFile1Name,
		  TString inputFile2Name,
		  TString inputFile3Name,
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

  const double xSec_Zbb0jets_  = 1.66; // in pb	  
  const double xSec_Zbb1jets_  = 0.29; // in pb	  
  const double xSec_Zbb2jets_  = 0.05; // in pb	  
  const double xSec_Zbb3jets_  = 0.01; // in pb	  

  const double leptonFactor = 2.0;

  const double luminosityFactor = 100000.; // 100fb^-1 of luminosity assumed in histograms

  TFile * inputFile1 = new TFile(inputFile1Name);
  TFile * inputFile2 = new TFile(inputFile2Name);
  TFile * inputFile3 = new TFile(inputFile3Name);
 
  TH1D * eventsNumberHisto1 = dynamic_cast<TH1D*>(inputFile1->Get(eventsNumberHistoName));
  TH1D * histo1 = dynamic_cast<TH1D*>(inputFile1->Get(HistoName));
  TH1D * histoirrBkg1 = dynamic_cast<TH1D*>(inputFile1->Get(HistoName+"irrBkg"));
  TH1D * eventsNumberHisto2 = dynamic_cast<TH1D*>(inputFile2->Get(eventsNumberHistoName));
  TH1D * histo2 = dynamic_cast<TH1D*>(inputFile2->Get(HistoName));
  TH1D * histoirrBkg2 = dynamic_cast<TH1D*>(inputFile2->Get(HistoName+"irrBkg"));
  TH1D * eventsNumberHisto3 = dynamic_cast<TH1D*>(inputFile3->Get(eventsNumberHistoName));
  TH1D * histo3 = dynamic_cast<TH1D*>(inputFile3->Get(HistoName));
  TH1D * histoirrBkg3 = dynamic_cast<TH1D*>(inputFile3->Get(HistoName+"irrBkg"));

  TString HistoTitle = histo1->GetTitle();
  int nbin = histo1->GetNbinsX();
  double xmin = histo1->GetXaxis()->GetXmin();
  double xmax = histo1->GetXaxis()->GetXmax();
  TH1D * histoTot       = new TH1D(HistoName+"TOT",      "inclusive "+HistoTitle,nbin,xmin,xmax);
  TH1D * histoirrBkgTot = new TH1D(HistoName+"irrBkgTOT","inclusive irrudicible background "+HistoTitle,nbin,xmin,xmax);

  histo1->SetTitle(HistoTitle+ " for Zbb 0 jets");
  histo2->SetTitle(HistoTitle+ " for Zbb 1 jets");
  histo3->SetTitle(HistoTitle+ " for Zbb 2 jets");
  histoirrBkg1->SetTitle("irriducible background "+HistoTitle+ " for Zbb 0 jets");
  histoirrBkg2->SetTitle("irriducible background "+HistoTitle+ " for Zbb 1 jets");
  histoirrBkg3->SetTitle("irriducible background "+HistoTitle+ " for Zbb 2 jets");

  double histo1Entries = histo1->GetEntries();
  double histo2Entries = histo2->GetEntries();
  double histo3Entries = histo3->GetEntries();
  double histoirrBkg1Entries = histoirrBkg1->GetEntries();
  double histoirrBkg2Entries = histoirrBkg2->GetEntries();
  double histoirrBkg3Entries = histoirrBkg3->GetEntries();
  double events1 = (double)eventsNumberHisto1->GetMaximum(); 
  double events2 = (double)eventsNumberHisto2->GetMaximum();
  double events3 = (double)eventsNumberHisto3->GetMaximum();
  double xSecEff1  = leptonFactor*xSec_Zbb0jets_;
  double xSecEff2  = leptonFactor*xSec_Zbb1jets_;
  double xSecEff3  = leptonFactor*xSec_Zbb2jets_;
  double invLuminosityEff1 = xSecEff1/events1;
  double invLuminosityEff2 = xSecEff2/events2;
  double invLuminosityEff3 = xSecEff3/events3;


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
  TLegend * legendTot = new TLegend(0.7,0.65,0.89,0.85);
  legendTot->SetFillColor(0);
  legendTot->SetBorderSize(0);  
  legendTot->SetTextFont(72);
  legendTot->SetTextSize(0.035);
  legendTot->SetFillColor(0);
  char nev[50];

  TCanvas * canvas = new TCanvas ( "diffxSec", "differential xSec", 1200, 400 );
  gStyle->SetOptStat(0);
  //  canvas->UseCurrentStyle();
  canvas->Divide(3,1);
  canvas->cd(1);
  histo1->SetNormFactor(histo1Entries*invLuminosityEff1*luminosityFactor);
  histo1->Draw();
  sprintf(nev,"number of events in 100 fb^{-1}:");
  sprintf(nev,"%.2f",histo1Entries*invLuminosityEff1*luminosityFactor);
  legend1->AddEntry(histo1,nev,"");
  legend1->Draw();
  canvas->cd(2);
  histo2->SetNormFactor(histo2Entries*invLuminosityEff2*luminosityFactor);
  histo2->Draw();
  legend2->Clear();
  sprintf(nev,"number of events in 100 fb^{-1}:");
  sprintf(nev,"%.2f",histo2Entries*invLuminosityEff2*luminosityFactor);
  legend2->AddEntry(histo2,nev,"");
  legend2->Draw();
  canvas->cd(3);
  histo3->SetNormFactor(histo3Entries*invLuminosityEff3*luminosityFactor);
  histo3->Draw();
  legend3->Clear();
  sprintf(nev,"number of events in 100 fb^{-1}:");
  sprintf(nev,"%.2f",histo3Entries*invLuminosityEff3*luminosityFactor);
  legend3->AddEntry(histo3,nev,"");
  legend3->Draw();

  histoTot->Sumw2();
  histoTot->Add(histo1);
  histoTot->Add(histo2);
  histoTot->Add(histo3);

  TCanvas * canvasTot = new TCanvas ( "inclusivediffxSec", "differential xSec", 1200, 400 );
  gStyle->SetOptStat(0);
  canvasTot->UseCurrentStyle();
  canvasTot->cd();
  histo1->SetTitle("4-body mass (100 fb^{-1})");
  histo1->GetXaxis()->SetTitle("m_{llbb} (GeV)");
  histo1->GetYaxis()->SetTitle("d#sigma/dm_{llbb} events/10GeV");
  histo1->SetMinimum(0.001);
  histo1->SetLineColor(51);
  histo1->Draw();
  histo2->SetLineColor(56);
  histo2->Draw("same");
  histo3->SetLineColor(60);
  histo3->Draw("same");
  histoTot->Draw("same");

  legendTot->Clear();
  sprintf(nev,"Zb\\bar{b}+0jets");
  legendTot->AddEntry(histo1,nev,"l");
  sprintf(nev,"Zb\\bar{b}+1jets");
  legendTot->AddEntry(histo2,nev,"l");
  sprintf(nev,"Zb\\bar{b}+2jets");
  legendTot->AddEntry(histo3,nev,"l");
  sprintf(nev,"inclusive");
  legendTot->AddEntry(histoTot,nev,"l");
  legendTot->Draw();



  TCanvas * canvasirrBkg = new TCanvas ( "diffxSecirrbkg", "differential xSec for irriducible background", 1200, 400 );
  gStyle->SetOptStat(0);
  canvasirrBkg->UseCurrentStyle();
  canvasirrBkg->Divide(3,1);
  canvasirrBkg->cd(1);
  histoirrBkg1->SetNormFactor(histoirrBkg1Entries*invLuminosityEff1*luminosityFactor);
  histoirrBkg1->Draw();
  sprintf(nev,"number of events in 100 fb^{-1}:");
  sprintf(nev,"%.2f",histoirrBkg1Entries*invLuminosityEff1*luminosityFactor);
  legend1->AddEntry(histoirrBkg1,nev,"");
  legend1->Draw();
  canvasirrBkg->cd(2);
  histoirrBkg2->SetNormFactor(histoirrBkg2Entries*invLuminosityEff2*luminosityFactor);
  histoirrBkg2->Draw();
  legend2->Clear();
  sprintf(nev,"number of events in 100 fb^{-1}:");
  sprintf(nev,"%.2f",histoirrBkg2Entries*invLuminosityEff2*luminosityFactor);
  legend2->AddEntry(histoirrBkg2,nev,"");
  legend2->Draw();
  canvasirrBkg->cd(3);
  histoirrBkg3->SetNormFactor(histoirrBkg3Entries*invLuminosityEff3*luminosityFactor);
  histoirrBkg3->Draw();
  legend3->Clear();
  sprintf(nev,"number of events in 100 fb^{-1}:");
  sprintf(nev,"%.2f",histoirrBkg3Entries*invLuminosityEff3*luminosityFactor);
  legend3->AddEntry(histoirrBkg3,nev,"");
  legend3->Draw();

  histoirrBkgTot->Sumw2();
  histoirrBkgTot->Add(histoirrBkg1);
  histoirrBkgTot->Add(histoirrBkg2);
  histoirrBkgTot->Add(histoirrBkg3);

  TCanvas * canvasirrBkgTot = new TCanvas ( "inclusivediffxSecirrBkg", "differential xSec", 1200, 400 );
  gStyle->SetOptStat(0);
  canvasirrBkgTot->UseCurrentStyle();
  canvasirrBkgTot->cd();
  gPad->SetLogy();
  histoirrBkg1->SetTitle("4-body mass (100 fb^{-1})");
  histoirrBkg1->GetXaxis()->SetTitle("m_{llbb} (GeV)");
  histoirrBkg1->GetYaxis()->SetTitle("d#sigma/dm_{llbb} events/10GeV");
  histoirrBkg1->SetMinimum(0.001);
  histoirrBkg1->SetLineColor(51);
  histoirrBkg1->Draw();
  histoirrBkg2->SetLineColor(56);
  histoirrBkg2->Draw("same");
  histoirrBkg3->SetLineColor(60);
  histoirrBkg3->Draw("same");
  histoirrBkgTot->Draw("same");

  legendTot->Clear();
  sprintf(nev,"Zb\\bar{b}+0jets");
  legendTot->AddEntry(histoirrBkg1,nev,"l");
  sprintf(nev,"Zb\\bar{b}+1jets");
  legendTot->AddEntry(histoirrBkg2,nev,"l");
  sprintf(nev,"Zb\\bar{b}+2jets");
  legendTot->AddEntry(histoirrBkg3,nev,"l");
  sprintf(nev,"inclusive");
  legendTot->AddEntry(histoirrBkgTot,nev,"l");
  legendTot->Draw();

  canvas->Print("Zbb.jpg");
  canvasTot->Print("Zbbtot.jpg");
  canvasirrBkg->Print("ZbbIrr.jpg");
  canvasirrBkgTot->Print("ZbbIrrtot.jpg");

  TFile * outputFile = new TFile("diffXsecZbb.root","RECREATE");

  histo1->Write();
  histo2->Write();
  histo3->Write();
  histoTot->Write();
  histoirrBkg1->Write();
  histoirrBkg2->Write();
  histoirrBkg3->Write();
  histoirrBkgTot->Write();

  outputFile->Close();

}


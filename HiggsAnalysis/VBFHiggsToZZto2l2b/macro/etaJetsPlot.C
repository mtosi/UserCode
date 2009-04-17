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


void etaJetsPlot( TString inputFileName1,
		  TString inputFileName2,
		  TString inputFileName3,
		  TString inputFileName4
		  ) {

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
  gStyle->SetPalette(216);

  TFile * inputFile1 = new TFile(inputFileName1);
  TFile * inputFile2 = new TFile(inputFileName2);
  TFile * inputFile3 = new TFile(inputFileName3);
  TFile * inputFile4 = new TFile(inputFileName4);
 
  TString HistoJetName(       "OfflinejetEta"       );
  TString HistoCentralJetName("OfflinecentraljetEta");
  TString HistoForwardJetName("OfflineforwardjetEta");
  TH1D * histo1Jet = dynamic_cast<TH1D*>(inputFile1->Get(HistoJetName));
  histo1Jet->SetTitle("jet |#eta| distribution");
  histo1Jet->GetXaxis()->SetTitle("jet |#eta|");
  TH1D * histo1CentralJets = dynamic_cast<TH1D*>(inputFile1->Get(HistoCentralJetName));
  histo1CentralJets->SetTitle("central jet |#eta| distribution");
  histo1CentralJets->GetXaxis()->SetTitle("central jet |#eta|");
  TH1D * histo1ForwardJets = dynamic_cast<TH1D*>(inputFile1->Get(HistoForwardJetName));
  histo1ForwardJets->SetTitle("forward jet jet |#eta| distribution");
  histo1ForwardJets->GetXaxis()->SetTitle("forward jet jet |#eta|");
  TH1D * histo2Jet = dynamic_cast<TH1D*>(inputFile2->Get(HistoJetName));
  histo2Jet->SetTitle("jet |#eta| distribution");
  histo2Jet->GetXaxis()->SetTitle("jet |#eta|");
  TH1D * histo2CentralJets = dynamic_cast<TH1D*>(inputFile2->Get(HistoCentralJetName));
  histo2CentralJets->SetTitle("central jet |#eta| distribution");
  histo2CentralJets->GetXaxis()->SetTitle("central jet |#eta|");
  TH1D * histo2ForwardJets = dynamic_cast<TH1D*>(inputFile2->Get(HistoForwardJetName));
  histo2ForwardJets->SetTitle("forward jet jet |#eta| distribution");
  histo2ForwardJets->GetXaxis()->SetTitle("forward jet jet |#eta|");
  TH1D * histo3Jet = dynamic_cast<TH1D*>(inputFile3->Get(HistoJetName));
  histo3Jet->SetTitle("jet |#eta| distribution");
  histo3Jet->GetXaxis()->SetTitle("jet |#eta|");
  TH1D * histo3CentralJets = dynamic_cast<TH1D*>(inputFile3->Get(HistoCentralJetName));
  histo3CentralJets->SetTitle("central jet |#eta| distribution");
  histo3CentralJets->GetXaxis()->SetTitle("central jet |#eta|");
  TH1D * histo3ForwardJets = dynamic_cast<TH1D*>(inputFile3->Get(HistoForwardJetName));
  histo3ForwardJets->SetTitle("forward jet jet |#eta| distribution");
  histo3ForwardJets->GetXaxis()->SetTitle("forward jet jet |#eta|");
  TH1D * histo4Jet = dynamic_cast<TH1D*>(inputFile4->Get(HistoJetName));
  histo4Jet->SetTitle("jet |#eta| distribution");
  histo4Jet->GetXaxis()->SetTitle("jet |#eta|");
  TH1D * histo4CentralJets = dynamic_cast<TH1D*>(inputFile4->Get(HistoCentralJetName));
  histo4CentralJets->SetTitle("central jet |#eta| distribution");
  histo4CentralJets->GetXaxis()->SetTitle("central jet |#eta|");
  TH1D * histo4ForwardJets = dynamic_cast<TH1D*>(inputFile4->Get(HistoForwardJetName));
  histo4ForwardJets->SetTitle("forward jet jet |#eta| distribution");
  histo4ForwardJets->GetXaxis()->SetTitle("forward jet jet |#eta|");


  TLegend * legend = new TLegend(0.85,0.8,0.99,0.99);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);  
  legend->SetTextFont(72);
  legend->SetTextSize(0.027);
  legend->SetFillColor(0);
  legend->AddEntry(histo1Jet,        "all jet",    "lpe");
  legend->AddEntry(histo1CentralJets,"central jet","lf");
  legend->AddEntry(histo1ForwardJets,"forward jet","lf");


  TCanvas * canvas = new TCanvas ( "jetEta", "jet |#eta|)", 1200, 800 );
  gStyle->SetOptStat(0);
  canvas->UseCurrentStyle();
  canvas->Divide(2,2);
  canvas->cd(1);
  histo1CentralJets->SetLineColor(62);
  histo1CentralJets->SetFillColor(62);
  histo1CentralJets->SetFillStyle(3002);
  //  histo1CentralJets->Scale(1./(double)histo1CentralJets->GetEntries());
  histo1CentralJets->Draw();
  histo1ForwardJets->SetLineColor(95);
  histo1ForwardJets->SetFillColor(95);
  histo1ForwardJets->SetFillStyle(3002);
  //  histo1ForwardJets->Scale(1./(double)histo1ForwardJets->GetEntries());
  histo1ForwardJets->Draw("same");
  histo1Jet->SetMarkerStyle(kFullDotMedium);
  histo1Jet->SetMarkerColor(51);
  //  histo1Jet->Scale(1./(double)histo1Jet->GetEntries());
  histo1Jet->Draw("pesame");
  legend->Draw();

  canvas->cd(2);
  histo2CentralJets->SetLineColor(62);
  histo2CentralJets->SetFillColor(62);
  histo2CentralJets->SetFillStyle(3002);
  //  histo2CentralJets->Scale(1./(double)histo2CentralJets->GetEntries());
  histo2CentralJets->Draw();
  histo2ForwardJets->SetLineColor(95);
  histo2ForwardJets->SetFillColor(95);
  histo2ForwardJets->SetFillStyle(3002);
  //  histo2ForwardJets->Scale(1./(double)histo2ForwardJets->GetEntries());
  histo2ForwardJets->Draw("same");
  histo2Jet->SetMarkerStyle(kFullDotMedium);
  histo2Jet->SetMarkerColor(51);
  //  histo2Jet->Scale(1./(double)histo2Jet->GetEntries());
  histo2Jet->Draw("pesame");
  legend->Draw();

  canvas->cd(3);
  histo3CentralJets->SetLineColor(62);
  histo3CentralJets->SetFillColor(62);
  histo3CentralJets->SetFillStyle(3002);
  //  histo3CentralJets->Scale(1./(double)histo3CentralJets->GetEntries());
  histo3CentralJets->Draw();
  histo3ForwardJets->SetLineColor(95);
  histo3ForwardJets->SetFillColor(95);
  histo3ForwardJets->SetFillStyle(3002);
  //  histo3ForwardJets->Scale(1./(double)histo3ForwardJets->GetEntries());
  histo3ForwardJets->Draw("same");
  histo3Jet->SetMarkerStyle(kFullDotMedium);
  histo3Jet->SetMarkerColor(51);
  //  histo3Jet->Scale(1./(double)histo3Jet->GetEntries());
  histo3Jet->Draw("pesame");
  legend->Draw();

  canvas->cd(4);
  histo4CentralJets->SetLineColor(62);
  histo4CentralJets->SetFillColor(62);
  histo4CentralJets->SetFillStyle(3002);
  //  histo4CentralJets->Scale(1./(double)histo4CentralJets->GetEntries());
  histo4CentralJets->Draw();
  histo4ForwardJets->SetLineColor(95);
  histo4ForwardJets->SetFillColor(95);
  histo4ForwardJets->SetFillStyle(3002);
  //  histo4ForwardJets->Scale(1./(double)histo4ForwardJets->GetEntries());
  histo4ForwardJets->Draw("same");
  histo4Jet->SetMarkerStyle(kFullDotMedium);
  histo4Jet->SetMarkerColor(51);
  //  histo4Jet->Scale(1./(double)histo4Jet->GetEntries());
  histo4Jet->Draw("pesame");
  legend->Draw();


}

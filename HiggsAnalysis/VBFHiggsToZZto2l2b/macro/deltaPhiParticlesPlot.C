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


void delatPhiParticlesPlot( TString inputFileName1,
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

  TFile * inputFile1 = new TFile(inputFileName1);
  TFile * inputFile2 = new TFile(inputFileName2);
  TFile * inputFile3 = new TFile(inputFileName3);
  TFile * inputFile4 = new TFile(inputFileName4);
 
  TString HistobbbarName(    "DeltaPhiZbbar"     );
  TString HistoElectronsName("DeltaPhiZelectrons");
  TString HistoMuonsName(    "DeltaPhiZmuons"    );
  TString HistoZZName(       "DeltaPhiHZZ"       );
  TH1D * histo1bbbar = dynamic_cast<TH1D*>(inputFile1->Get(HistobbbarName));
  histo1bbbar->SetTitle("#Delta#Phi between Zdecay particles");
  histo1bbbar->GetXaxis()->SetTitle("#Delta#Phi (rad)");
  TH1D * histo1Electrons = dynamic_cast<TH1D*>(inputFile1->Get(HistoElectronsName));
  histo1Electrons->SetTitle("#Delta#Phi between Zdecay particles");
  histo1Electrons->GetXaxis()->SetTitle("#Delta#Phi (rad)");
  TH1D * histo1Muons = dynamic_cast<TH1D*>(inputFile1->Get(HistoMuonsName));
  histo1Muons->SetTitle("#Delta#Phi between Zdecay particles");
  histo1Muons->GetXaxis()->SetTitle("#Delta#Phi (rad)");
  TH1D * histo1ZZ =  dynamic_cast<TH1D*>(inputFile1->Get(HistoZZName));
  histo1ZZ->SetTitle("#Delta#Phi between ZZ");
  histo1ZZ->GetXaxis()->SetTitle("#Delta#Phi (rad)");
  TH1D * histo2bbbar = dynamic_cast<TH1D*>(inputFile2->Get(HistobbbarName));
  histo2bbbar->SetTitle("#Delta#Phi between Zdecay particles");
  histo2bbbar->GetXaxis()->SetTitle("#Delta#Phi (rad)");
  TH1D * histo2Electrons = dynamic_cast<TH1D*>(inputFile2->Get(HistoElectronsName));
  histo2Electrons->SetTitle("#Delta#Phi between Zdecay particles");
  histo2Electrons->GetXaxis()->SetTitle("#Delta#Phi (rad)");
  TH1D * histo2Muons = dynamic_cast<TH1D*>(inputFile2->Get(HistoMuonsName));
  histo2Muons->SetTitle("#Delta#Phi between Zdecay particles");
  histo2Muons->GetXaxis()->SetTitle("#Delta#Phi (rad)");
  TH1D * histo2ZZ =  dynamic_cast<TH1D*>(inputFile2->Get(HistoZZName));
  histo2ZZ->SetTitle("#Delta#Phi between ZZ");
  histo2ZZ->GetXaxis()->SetTitle("#Delta#Phi (rad)");
  TH1D * histo3bbbar = dynamic_cast<TH1D*>(inputFile3->Get(HistobbbarName));
  histo3bbbar->SetTitle("#Delta#Phi between Zdecay particles");
  histo3bbbar->GetXaxis()->SetTitle("#Delta#Phi (rad)");
  TH1D * histo3Electrons = dynamic_cast<TH1D*>(inputFile3->Get(HistoElectronsName));
  histo3Electrons->SetTitle("#Delta#Phi between Zdecay particles");
  histo3Electrons->GetXaxis()->SetTitle("#Delta#Phi (rad)");
  TH1D * histo3Muons = dynamic_cast<TH1D*>(inputFile3->Get(HistoMuonsName));
  histo3Muons->SetTitle("#Delta#Phi between Zdecay particles");
  histo3Muons->GetXaxis()->SetTitle("#Delta#Phi (rad)");
  TH1D * histo3ZZ =  dynamic_cast<TH1D*>(inputFile3->Get(HistoZZName));
  histo3ZZ->SetTitle("#Delta#Phi between ZZ");
  histo3ZZ->GetXaxis()->SetTitle("#Delta#Phi (rad)");
  TH1D * histo4bbbar = dynamic_cast<TH1D*>(inputFile4->Get(HistobbbarName));
  histo4bbbar->SetTitle("#Delta#Phi between Zdecay particles");
  histo4bbbar->GetXaxis()->SetTitle("#Delta#Phi (rad)");
  TH1D * histo4Electrons = dynamic_cast<TH1D*>(inputFile4->Get(HistoElectronsName));
  histo4Electrons->SetTitle("#Delta#Phi between Zdecay particles");
  histo4Electrons->GetXaxis()->SetTitle("#Delta#Phi (rad)");
  TH1D * histo4Muons = dynamic_cast<TH1D*>(inputFile4->Get(HistoMuonsName));
  histo4Muons->SetTitle("#Delta#Phi between Zdecay particles");
  histo4Muons->GetXaxis()->SetTitle("#Delta#Phi (rad)");
  TH1D * histo4ZZ =  dynamic_cast<TH1D*>(inputFile4->Get(HistoZZName));
  histo4ZZ->SetTitle("#Delta#Phi between ZZ");
  histo4ZZ->GetXaxis()->SetTitle("#Delta#Phi (rad)");

  TLegend * legend = new TLegend(0.85,0.8,0.99,0.99);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);  
  legend->SetTextFont(72);
  legend->SetTextSize(0.027);
  legend->SetFillColor(0);
  legend->AddEntry(histo1bbbar,    "b\\bar{b}",     "l");
  legend->AddEntry(histo1Electrons,"e^{=}e^{-}",    "l");
  legend->AddEntry(histo1Muons,    "#mu^{+}#mu^{-}","l");

  TCanvas * canvasZ = new TCanvas ( "deltaPhiZ", "#Delta#Phi (rad)", 1200, 800 );
  gStyle->SetOptStat(0);
  canvasZ->UseCurrentStyle();
  canvasZ->Divide(2,2);
  canvasZ->cd(1);
  histo1Muons->SetLineColor(51);
  histo1Muons->Scale(1./histo1Muons->GetEntries());
  histo1Muons->Draw();
  histo1bbbar->SetLineColor(79);
  histo1bbbar->Scale(1./(double)histo1bbbar->GetEntries());
  histo1bbbar->Draw("same");
  histo1Electrons->SetLineColor(62);
  histo1Electrons->Scale(1./(double)histo1Electrons->GetEntries());
  histo1Electrons->Draw("same");
  legend->Draw();

  canvasZ->cd(2);
  histo2Muons->SetLineColor(51);
  histo2Muons->Scale(1./histo2Muons->GetEntries());
  histo2Muons->Draw();
  histo2bbbar->SetLineColor(79);
  histo2bbbar->Scale(1./(double)histo2bbbar->GetEntries());
  histo2bbbar->Draw("same");
  histo2Electrons->SetLineColor(62);
  histo2Electrons->Scale(1./(double)histo2Electrons->GetEntries());
  histo2Electrons->Draw("same");
  legend->Draw();

  canvasZ->cd(3);
  histo3Muons->SetLineColor(51);
  histo3Muons->Scale(1./histo3Muons->GetEntries());
  histo3Muons->Draw();
  histo3bbbar->SetLineColor(79);
  histo3bbbar->Scale(1./(double)histo3bbbar->GetEntries());
  histo3bbbar->Draw("same");
  histo3Electrons->SetLineColor(62);
  histo3Electrons->Scale(1./(double)histo3Electrons->GetEntries());
  histo3Electrons->Draw("same");
  legend->Draw();

  canvasZ->cd(4);
  histo4Muons->SetLineColor(51);
  histo4Muons->Scale(1./histo4Muons->GetEntries());
  histo4Muons->Draw();
  histo4bbbar->SetLineColor(79);
  histo4bbbar->Scale(1./(double)histo4bbbar->GetEntries());
  histo4bbbar->Draw("same");
  histo4Electrons->SetLineColor(62);
  histo4Electrons->Scale(1./(double)histo4Electrons->GetEntries());
  histo4Electrons->Draw("same");
  legend->Draw();

  TCanvas * canvasH = new TCanvas ( "deltaPhiH", "#Delta#Phi (rad)", 1200, 800 );
  gStyle->SetOptStat(0);
  canvasH->UseCurrentStyle();
  canvasH->Divide(2,2);
  canvasH->cd(1);
  histo1ZZ->SetLineColor(215);
  histo1ZZ->Scale(1./histo1ZZ->GetEntries());
  histo1ZZ->Draw();

  canvasH->cd(2);
  histo2ZZ->SetLineColor(215);
  histo2ZZ->Scale(1./histo2ZZ->GetEntries());
  histo2ZZ->Draw();

  canvasH->cd(3);
  histo3ZZ->SetLineColor(215);
  histo3ZZ->Scale(1./histo3ZZ->GetEntries());
  histo3ZZ->Draw();

  canvasH->cd(4);
  histo4ZZ->SetLineColor(215);
  histo4ZZ->Scale(1./histo4ZZ->GetEntries());
  histo4ZZ->Draw();

}

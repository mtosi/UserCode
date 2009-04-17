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


void diffXsecplot( TString inputSgnFileName,
		   TString inputBkgFileName,
		   TString HistoSgnName,
		   TString HistoBkgName) {
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

  TFile * inputSgnFile = new TFile(inputSgnFileName);
  TFile * inputBkgFile = new TFile(inputBkgFileName);
 
  TH1D * histoSgn1   = dynamic_cast<TH1D*>(inputSgnFile->Get(HistoSgnName+"160"));
  TH1D * histoSgn2   = dynamic_cast<TH1D*>(inputSgnFile->Get(HistoSgnName+"200"));
  TH1D * histoSgn3   = dynamic_cast<TH1D*>(inputSgnFile->Get(HistoSgnName+"400"));
  TH1D * histoSgn4   = dynamic_cast<TH1D*>(inputSgnFile->Get(HistoSgnName+"800"));
  TH1D * histoSgnTot = dynamic_cast<TH1D*>(inputSgnFile->Get(HistoSgnName+"TOT"));
  TH1D * histoBkg    = dynamic_cast<TH1D*>(inputBkgFile->Get(HistoBkgName+"TOT"));
  TH1D * histoirrBkg = dynamic_cast<TH1D*>(inputBkgFile->Get(HistoBkgName+"irrBkgTOT"));

  TLegend * legend = new TLegend(0.5,0.6,0.99,0.99);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);  
  legend->SetTextFont(72);
  legend->SetTextSize(0.027);
  legend->SetFillColor(0);

  histoSgn1->SetLineColor(70);
  histoSgn1->SetFillColor(70);
  histoSgn2->SetLineColor(76);
  histoSgn2->SetFillColor(76);
  histoSgn3->SetLineColor(82);
  histoSgn3->SetFillColor(82);
  histoSgn4->SetLineColor(86);
  histoSgn4->SetFillColor(86);
  histoBkg->SetLineColor(52);
  histoBkg->SetFillColor(52);
  histoBkg->SetFillStyle(3003);
  histoBkg->SetTitle("4-body mass (100 fb^{-1})");
  histoBkg->GetXaxis()->SetTitle("m_{llbb} (GeV)");
  histoBkg->GetYaxis()->SetTitle("d#sigma/dm_{llbb} events/10GeV");
  histoirrBkg->SetLineColor(52);
  histoirrBkg->SetFillColor(52);
  histoirrBkg->SetFillStyle(3003);
  histoirrBkg->SetTitle("4-body mass (100 fb^{-1})");
  histoirrBkg->GetXaxis()->SetTitle("m_{llbb} (GeV)");
  histoirrBkg->GetYaxis()->SetTitle("d#sigma/dm_{llbb} events/10GeV");
  legend->AddEntry(histoSgn1,"HZZ w/ m_{H}=160 GeV","f");
  legend->AddEntry(histoSgn2,"HZZ w/ m_{H}=200 GeV","f");
  legend->AddEntry(histoSgn3,"HZZ w/ m_{H}=400 GeV","f");
  legend->AddEntry(histoSgn4,"HZZ w/ m_{H}=800 GeV","f");
  legend->AddEntry(histoBkg,"inclusive Zbb+njets","f");

  TCanvas * canvas = new TCanvas ( "diffxSec", "differential xSec", 1200, 600 );
  gStyle->SetOptStat(0);
  canvas->UseCurrentStyle();
  canvas->Divide(2,1);
  canvas->cd(1);
  gPad->SetLogy();
  histoBkg->SetMinimum(0.001);
  histoBkg->Draw("histo");
  histoSgn1->Draw("same");
  histoSgn2->Draw("same");
  histoSgn3->Draw("same");
  histoSgn4->Draw("same");
  legend->Draw();

  legend->Clear();
  legend->AddEntry(histoSgn1,"HZZ w/ m_{H}=160 GeV","f");
  legend->AddEntry(histoSgn2,"HZZ w/ m_{H}=200 GeV","f");
  legend->AddEntry(histoSgn3,"HZZ w/ m_{H}=400 GeV","f");
  legend->AddEntry(histoSgn4,"HZZ w/ m_{H}=800 GeV","f");
  legend->AddEntry(histoirrBkg,"irrudicible inclusive Zbb","f");
  canvas->cd(2);
  gPad->SetLogy();
  histoirrBkg->SetMinimum(0.001);
  histoirrBkg->Draw("histo");
  histoSgn1->Draw("same");
  histoSgn2->Draw("same");
  histoSgn3->Draw("same");
  histoSgn4->Draw("same");
  legend->Draw();


  canvas->Print("HZZandZbb.jpg");
  canvas->Print("HZZandZbb.eps");

}

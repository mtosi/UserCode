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


void signalPlots( ) {
  // general root setting
  gROOT->Reset(); 
  //  gROOT->SetBatch(kTRUE);
  //  gStyle->SetOptStat("nemrou");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(111);
  gStyle->SetTitleSize(0.1);
  gStyle->SetLabelSize(0.04,"X");
  gStyle->SetLabelSize(0.04,"Y");
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);

  TFile * inputSignalFile1 = new TFile("VBFH160ZZoutput.root");
  TFile * inputSignalFile2 = new TFile("VBFH200ZZoutput.root");
  //TFile * inputSignalFile3 = new TFile("VBFH400ZZoutput.root");
  TFile * inputSignalFile4 = new TFile("VBFH800ZZoutput.root");
 
  TFile * inputBackgroundFile1 = new TFile("Zbb0Jetsoutput.root");
  TFile * inputBackgroundFile2 = new TFile("Zbb1Jetsoutput.root");
  TFile * inputBackgroundFile3 = new TFile("Zbb2Jetsoutput.root");

  TFile * inputBackgroundFile4 = new TFile("ZZ0Jetsoutput.root");
  TFile * inputBackgroundFile5 = new TFile("ZZ1Jetsoutput.root");
  TFile * inputBackgroundFile6 = new TFile("ZZ2Jetsoutput.root");
  TFile * inputBackgroundFile7 = new TFile("ZZINCLJetsoutput.root");

 
  TString HistoNameMC("DeltaRZZ");
  TString xAxisTitle("#DeltaR(Z_{1},Z_{2})");

  TH1D * sgnHistoMC1 = dynamic_cast<TH1D*>(inputSignalFile1->Get(HistoNameMC));
  TString title = sgnHistoMC1->GetTitle();
  sgnHistoMC1->SetTitle(title+ " for m_{H} = 160GeV/c^{2}");
  sgnHistoMC1->GetXaxis()->SetTitle(xAxisTitle);
  TH1D * sgnHistoMC2 = dynamic_cast<TH1D*>(inputSignalFile2->Get(HistoNameMC));
  sgnHistoMC2->SetTitle(title+ " for m_{H} = 200GeV");
  sgnHistoMC2->GetXaxis()->SetTitle(xAxisTitle);
//  TH1D * sgnHistoMC3 = dynamic_cast<TH1D*>(inputSignalFile3->Get(HistoNameMC));
//  sgnHistoMC3->SetTitle(title+ " for m_{H} = 400GeV");
//  sgnHistoMC3->GetXaxis()->SetTitle(xAxisTitle);
  TH1D * sgnHistoMC4 = dynamic_cast<TH1D*>(inputSignalFile4->Get(HistoNameMC));
  sgnHistoMC4->SetTitle(title+ " for m_{H} = 800GeV");
  sgnHistoMC4->GetXaxis()->SetTitle(xAxisTitle);

  TString HistoNameREC("recHZZdeltaR");
  TH1D * sgnHistoREC1 = dynamic_cast<TH1D*>(inputSignalFile1->Get(HistoNameREC));
  title = sgnHistoREC1->GetTitle();
  sgnHistoREC1->SetTitle(title+ " for m_{H} = 160GeV/c^{2}");
  sgnHistoREC1->GetXaxis()->SetTitle(xAxisTitle);
  TH1D * sgnHistoREC2 = dynamic_cast<TH1D*>(inputSignalFile2->Get(HistoNameREC));
  sgnHistoREC2->SetTitle(title+ " for m_{H} = 200GeV");
  sgnHistoREC2->GetXaxis()->SetTitle(xAxisTitle);
//  TH1D * sgnHistoREC3 = dynamic_cast<TH1D*>(inputSignalFile3->Get(HistoNameREC));
//  sgnHistoREC3->SetTitle(title+ " for m_{H} = 400GeV");
//  sgnHistoREC3->GetXaxis()->SetTitle(xAxisTitle);
  TH1D * sgnHistoREC4 = dynamic_cast<TH1D*>(inputSignalFile4->Get(HistoNameREC));
  sgnHistoREC4->SetTitle(title+ " for m_{H} = 800GeV");
  sgnHistoREC4->GetXaxis()->SetTitle(xAxisTitle);

  TH1D * bkgZbb0HistoREC = dynamic_cast<TH1D*>(inputBackgroundFile1->Get(HistoNameREC));
  title = bkgZbb0HistoREC->GetTitle();
  bkgZbb0HistoREC->SetTitle(title+ " for Zbb + 0 jets");
  bkgZbb0HistoREC->GetXaxis()->SetTitle(xAxisTitle);
  TH1D * bkgZbb1HistoREC = dynamic_cast<TH1D*>(inputBackgroundFile2->Get(HistoNameREC));
  title = bkgZbb1HistoREC->GetTitle();
  bkgZbb1HistoREC->SetTitle(title+ " for Zbb + 1 jets");
  bkgZbb1HistoREC->GetXaxis()->SetTitle(xAxisTitle);
  TH1D * bkgZbb2HistoREC = dynamic_cast<TH1D*>(inputBackgroundFile3->Get(HistoNameREC));
  title = bkgZbb2HistoREC->GetTitle();
  bkgZbb2HistoREC->SetTitle(title+ " for Zbb + 2 jets");
  bkgZbb2HistoREC->GetXaxis()->SetTitle(xAxisTitle);

  TH1D * bkgZZ0HistoREC = dynamic_cast<TH1D*>(inputBackgroundFile4->Get(HistoNameREC));
  title = bkgZZ0HistoREC->GetTitle();
  bkgZZ0HistoREC->SetTitle(title+ " for ZZ + 0 jets");
  bkgZZ0HistoREC->GetXaxis()->SetTitle(xAxisTitle);
  TH1D * bkgZZ1HistoREC = dynamic_cast<TH1D*>(inputBackgroundFile5->Get(HistoNameREC));
  title = bkgZZ1HistoREC->GetTitle();
  bkgZZ1HistoREC->SetTitle(title+ " for ZZ + 1 jets");
  bkgZZ1HistoREC->GetXaxis()->SetTitle(xAxisTitle);
  TH1D * bkgZZ2HistoREC = dynamic_cast<TH1D*>(inputBackgroundFile6->Get(HistoNameREC));
  title = bkgZZ2HistoREC->GetTitle();
  bkgZZ2HistoREC->SetTitle(title+ " for ZZ + 2 jets");
  bkgZZ2HistoREC->GetXaxis()->SetTitle(xAxisTitle);
 
  // output file
  TFile *outputfile;
  outputfile = dynamic_cast<TFile*>(gROOT->FindObject("signalPlots.root")); 
  if (outputfile) outputfile->Close();
  outputfile = new TFile("signalPlots.root","RECREATE","signal histograms");

  char nev[50];
  TLegend * legend1 = new TLegend(0.85,0.8,0.99,0.99);
  legend1->SetFillColor(0);
  legend1->SetBorderSize(0);  
  legend1->SetTextFont(72);
  legend1->SetTextSize(0.027);
  legend1->SetFillColor(0);
  TLegend * legend2 = new TLegend(0.85,0.8,0.99,0.99);
  legend2->SetFillColor(0);
  legend2->SetBorderSize(0);  
  legend2->SetTextFont(72);
  legend2->SetTextSize(0.027);
  legend2->SetFillColor(0);
  TLegend * legend3 = new TLegend(0.85,0.8,0.99,0.99);
  legend3->SetFillColor(0);
  legend3->SetBorderSize(0);  
  legend3->SetTextFont(72);
  legend3->SetTextSize(0.027);
  legend3->SetFillColor(0);
  TLegend * legend4 = new TLegend(0.85,0.8,0.99,0.99);
  legend4->SetFillColor(0);
  legend4->SetBorderSize(0);  
  legend4->SetTextFont(72);
  legend4->SetTextSize(0.027);
  legend4->SetFillColor(0);

  TCanvas * canvas1 = new TCanvas ( "c1", "c1", 1200, 800 );
  canvas1->UseCurrentStyle();
  canvas1->Divide(2,2);
  canvas1->cd(1);
  sgnHistoMC1->SetLineColor(70);
  sgnHistoMC1->Draw();
  sgnHistoREC1->SetLineColor(70);
  sgnHistoREC1->SetFillColor(70);
  sgnHistoREC1->SetFillStyle(3002);
  sgnHistoREC1->Draw("same");
  sprintf(nev,"MC: %d",int(sgnHistoMC1->GetEntries()));
  legend1->AddEntry(sgnHistoMC1,nev, "bf");
  sprintf(nev,"rec: %d",int(sgnHistoREC1->GetEntries()));
  legend1->AddEntry(sgnHistoREC1,nev,"bf");
  legend1->Draw();
  canvas1->cd(2);
  sgnHistoMC2->SetLineColor(76);
  sgnHistoMC2->Draw();
  sgnHistoREC2->SetLineColor(76);
  sgnHistoREC2->SetFillColor(76);
  sgnHistoREC2->SetFillStyle(3002);
  sgnHistoREC2->Draw("same");
  sprintf(nev,"MC: %d",int(sgnHistoMC2->GetEntries()));
  legend2->AddEntry(sgnHistoMC2,nev, "lf");
  sprintf(nev,"rec: %d",int(sgnHistoREC2->GetEntries()));
  legend2->AddEntry(sgnHistoREC2,nev,"lf");
  legend2->Draw();
  canvas1->cd(3);
//  sgnHistoMC3->SetLineColor(82);
//  sgnHistoMC3->Draw();
//  sgnHistoREC3->SetLineColor(82);
//  sgnHistoREC3->SetFillColor(82);
//  sgnHistoREC3->SetFillStyle(3002);
//  sgnHistoREC3->Draw("same");
//  sprintf(nev,"MC: %d",int(sgnHistoMC3->GetEntries()));
//  legend3->AddEntry(sgnHistoMC3,nev, "lf");
//  sprintf(nev,"rec: %d",int(sgnHistoREC3->GetEntries()));
//  legend3->AddEntry(sgnHistoREC3,nev,"lf");
//  legend3->Draw();
  canvas1->cd(4);
  sgnHistoMC4->SetLineColor(86);
  sgnHistoMC4->Draw();
  sgnHistoREC4->SetLineColor(86);
  sgnHistoREC4->SetFillColor(86);
  sgnHistoREC4->SetFillStyle(3002);
  sgnHistoREC4->Draw("same");
  sprintf(nev,"MC: %d",int(sgnHistoMC4->GetEntries()));
  legend4->AddEntry(sgnHistoMC4,nev, "lf");
  sprintf(nev,"rec: %d",int(sgnHistoREC4->GetEntries()));
  legend4->AddEntry(sgnHistoREC4,nev,"lf");
  legend4->Draw();


  TLegend * legend = new TLegend(0.85,0.8,0.99,0.99);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);  
  legend->SetTextFont(72);
  legend->SetTextSize(0.027);
  legend->SetFillColor(0);
  TCanvas * canvas2 = new TCanvas ( "c2", "c2", 1200, 800 );
  canvas2->UseCurrentStyle();
  canvas2->cd();
  sgnHistoREC4->SetLineColor(86);
  sgnHistoREC4->SetFillColor(86);
  sgnHistoREC4->Draw();
  bkgZZ0HistoREC->SetLineColor(51);
  bkgZZ0HistoREC->SetFillColor(51);
  bkgZZ0HistoREC->SetFillStyle(3002);
  bkgZZ0HistoREC->Draw("same");
  bkgZZ1HistoREC->SetLineColor(55);
  bkgZZ1HistoREC->SetFillColor(55);
  bkgZZ1HistoREC->SetFillStyle(3004);
  bkgZZ1HistoREC->Draw("same");
  bkgZZ2HistoREC->SetLineColor(59);
  bkgZZ2HistoREC->SetFillColor(59);
  bkgZZ2HistoREC->SetFillStyle(3005);
  bkgZZ2HistoREC->Draw("same");
  sprintf(nev,"rec: %d",int(sgnHistoREC4->GetEntries()));
  legend->AddEntry(sgnHistoREC4,nev,"bf");
  sprintf(nev,"rec: %d",int(bkgZZ0HistoREC->GetEntries()));
  legend->AddEntry(bkgZZ0HistoREC,nev,"bf");
  sprintf(nev,"rec: %d",int(bkgZZ1HistoREC->GetEntries()));
  legend->AddEntry(bkgZZ1HistoREC,nev,"bf");
  sprintf(nev,"rec: %d",int(bkgZZ2HistoREC->GetEntries()));
  legend->AddEntry(bkgZZ2HistoREC,nev,"bf");
  legend->Draw();

}

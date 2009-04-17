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


void mZ1vsmZ2plot( TString inputFileName1,
		   TString inputFileName2,
		   TString inputFileName3,
		   TString inputFileName4){
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
 
  TString HistoName("MZ1vsMZ2");
  TH2D * histo1 = dynamic_cast<TH2D*>(inputFile1->Get(HistoName));
  TString title = histo1->GetTitle();
  histo1->SetTitle(title+ " for m_{H} = 160GeV/c^{2}");
  histo1->GetXaxis()->SetTitle("m_{Z1} GeV/c^{2}");
  histo1->GetYaxis()->SetTitle("m_{Z2} GeV/c^{2}");
  TH2D * histo2 = dynamic_cast<TH2D*>(inputFile2->Get(HistoName));
  histo2->SetTitle(title+ " for m_{H} = 200GeV");
  histo2->GetXaxis()->SetTitle("m_{Z1} GeV/c^{2}");
  histo2->GetYaxis()->SetTitle("m_{Z2} GeV/c^{2}");
  TH2D * histo3 = dynamic_cast<TH2D*>(inputFile3->Get(HistoName));
  histo3->SetTitle(title+ " for m_{H} = 400GeV");
  histo3->GetXaxis()->SetTitle("m_{Z1} GeV/c^{2}");
  histo3->GetYaxis()->SetTitle("m_{Z2} GeV/c^{2}");
  TH2D * histo4 = dynamic_cast<TH2D*>(inputFile4->Get(HistoName));
  histo4->SetTitle(title+ " for m_{H} = 800GeV");
  histo4->GetXaxis()->SetTitle("m_{Z1} GeV/c^{2}");
  histo4->GetYaxis()->SetTitle("m_{Z2} GeV/c^{2}");

  TCanvas * canvas = new TCanvas ( "diffxSec", "differential xSec", 1200, 800 );
  gStyle->SetOptStat(0);
  canvas->UseCurrentStyle();
  canvas->Divide(2,2);
  canvas->cd(1);
  histo1->SetMarkerColor(70);
  histo1->Draw();
  canvas->cd(2);
  histo2->SetMarkerColor(76);
  histo2->Draw();
  canvas->cd(3);
  histo3->SetMarkerColor(82);
  histo3->Draw();
  canvas->cd(4);
  histo4->SetMarkerColor(86);
  histo4->Draw();

  canvas->Print("mZ1vsmZ2.jpg");
  canvas->Print("mZ1vsmZ2.eps");

}

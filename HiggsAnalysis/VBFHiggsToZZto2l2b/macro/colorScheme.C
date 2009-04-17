#include <string.h>
#include <sstream>
#include "TStyle.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "Riostream.h"

#include "AnalysisExamples/AnalysisClasses/interface/ListFashionAttributedHisto.h"
//#include "AnalysisExamples/AnalysisClasses/interface/FashionAttributedHisto.h"

void colorScheme() {

  using namespace anaobj;

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

 TLegend *lLegend0 = new TLegend(0.05,0.05,0.95,0.95,"VBF H->ZZ");
 lLegend0->SetFillColor(0);
 lLegend0->SetTextFont(72);
 lLegend0->SetTextSize(0.09);
 TLegend *lLegend1 = new TLegend(0.05,0.05,0.95,0.95,"ZZ+nJETS");
 lLegend1->SetFillColor(0);
 lLegend1->SetTextFont(72);
 lLegend1->SetTextSize(0.09);
 TLegend *lLegend2 = new TLegend(0.05,0.05,0.95,0.95,"Zbb+nJETS");
 lLegend2->SetFillColor(0);
 lLegend2->SetTextFont(72);
 lLegend2->SetTextSize(0.09);
 TLegend *lLegend3 = new TLegend(0.05,0.05,0.95,0.95,"WZ+nJETS");
 lLegend3->SetFillColor(0);
 lLegend3->SetTextFont(72);
 lLegend3->SetTextSize(0.09);
 TLegend *lLegend4 = new TLegend(0.05,0.05,0.95,0.95,"tt+nJETS");
 lLegend4->SetFillColor(0);
 lLegend4->SetTextFont(72);
 lLegend4->SetTextSize(0.09);

 std::ostringstream name; 
 name << "VBFH" << 160 << "ZZ";
 TH1D *pH1 = new TH1D(name.str().c_str(),name.str().c_str(),0,0,0);
 name.str("");
 name << "VBFH" << 200 << "ZZ";
 TH1D *pH2 = new TH1D(name.str().c_str(),name.str().c_str(),0,0,0);
 name.str("");
 name << "VBFH" << 400 << "ZZ";
 TH1D *pH3 = new TH1D(name.str().c_str(),name.str().c_str(),0,0,0);
 name.str("");
 name << "VBFH" << 800 << "ZZ";
 TH1D *pH4 = new TH1D(name.str().c_str(),name.str().c_str(),0,0,0);
 name.str("");

 name << "ZZnJETS";
 TH1D *pZZ = new TH1D(name.str().c_str(),name.str().c_str(),0,0,0);
 name.str("");
 name << "ZZ" << 0 << "JETS";
 TH1D *pZZ1 = new TH1D(name.str().c_str(),name.str().c_str(),0,0,0);
 name.str("");
 name << "ZZ" << 1 << "JETS";
 TH1D *pZZ2 = new TH1D(name.str().c_str(),name.str().c_str(),0,0,0);
 name.str("");
 name << "ZZ" << 2 << "JETS";
 TH1D *pZZ3 = new TH1D(name.str().c_str(),name.str().c_str(),0,0,0);
 name.str("");

 name << "ZbbnJETS";
 TH1D *pZbb = new TH1D(name.str().c_str(),name.str().c_str(),0,0,0);
 name.str("");
 name << "Zbb" << 0 << "JETS";
 TH1D *pZbb1 = new TH1D(name.str().c_str(),name.str().c_str(),0,0,0);
 name.str("");
 name << "Zbb" << 1 << "JETS";
 TH1D *pZbb2 = new TH1D(name.str().c_str(),name.str().c_str(),0,0,0);
 name.str("");
 name << "Zbb" << 2 << "JETS";
 TH1D *pZbb3 = new TH1D(name.str().c_str(),name.str().c_str(),0,0,0);
 name.str("");

 name << "WZnJETS";
 TH1D *pWZ = new TH1D(name.str().c_str(),name.str().c_str(),0,0,0);
 name.str("");
 name << "WZ" << 0 << "JETS";
 TH1D *pWZ1 = new TH1D(name.str().c_str(),name.str().c_str(),0,0,0);
 name.str("");
 name << "WZ" << 1 << "JETS";
 TH1D *pWZ2 = new TH1D(name.str().c_str(),name.str().c_str(),0,0,0);
 name.str("");
 name << "WZ" << 2 << "JETS";
 TH1D *pWZ3 = new TH1D(name.str().c_str(),name.str().c_str(),0,0,0);
 name.str("");

 name << "ttnJETS";
 TH1D *ptt = new TH1D(name.str().c_str(),name.str().c_str(),0,0,0);
 name.str("");
 name << "tt" << 0 << "JETS";
 TH1D *ptt1 = new TH1D(name.str().c_str(),name.str().c_str(),0,0,0);
 name.str("");
 name << "tt" << 1 << "JETS";
 TH1D *ptt2 = new TH1D(name.str().c_str(),name.str().c_str(),0,0,0);
 name.str("");
 name << "tt" << 2 << "JETS";
 TH1D *ptt3 = new TH1D(name.str().c_str(),name.str().c_str(),0,0,0);
 name.str("");

 pH1->SetLineColor(86); 
 pH1->SetFillColor(86);
 // pH1->SetFillStyle(3003);
 pH1->SetMarkerColor(86);
 pH1->SetMarkerStyle(20);
 pH1->SetMarkerSize(0.4);
 pH2->SetLineColor(84); 
 pH2->SetFillColor(84);
 pH2->SetFillStyle(3001);
 pH2->SetMarkerColor(84);
 pH2->SetMarkerStyle(20);
 pH2->SetMarkerSize(0.6);
 pH3->SetLineColor(211); 
 pH3->SetFillColor(211);
 pH3->SetFillStyle(3009);
 pH3->SetMarkerColor(211);
 pH3->SetMarkerStyle(20);
 pH3->SetMarkerSize(0.8);
 pH4->SetLineColor(210); 
 pH4->SetFillColor(210);
 pH4->SetFillStyle(3015);
 pH4->SetMarkerColor(210);
 pH4->SetMarkerStyle(20);
 pH4->SetMarkerSize(1.0);

 ptt->SetLineColor(92); 
 ptt->SetFillColor(92);
 // ptt->SetFillStyle(3001);
 ptt->SetMarkerColor(92);
 ptt->SetMarkerStyle(21);
 // ptt->SetMarkerSize(0.4);
 ptt1->SetLineColor(91); 
 ptt1->SetFillColor(91);
 ptt1->SetFillStyle(3001);
 ptt1->SetMarkerColor(91);
 ptt1->SetMarkerStyle(21);
 ptt1->SetMarkerSize(0.4);
 ptt2->SetLineColor(92); 
 ptt2->SetFillColor(92);
 ptt2->SetFillStyle(3009);
 ptt2->SetMarkerColor(92);
 ptt2->SetMarkerStyle(21);
 ptt2->SetMarkerSize(0.6);
 ptt3->SetLineColor(93); 
 ptt3->SetFillColor(93);
 ptt3->SetFillStyle(3015);
 ptt3->SetMarkerColor(93);
 ptt3->SetMarkerStyle(21);
 ptt3->SetMarkerSize(0.8);

 pZbb->SetLineColor(215); 
 pZbb->SetFillColor(215);
 // pZbb->SetFillStyle(3001);
 pZbb->SetMarkerColor(215);
 pZbb->SetMarkerStyle(22);
 // pZbb->SetMarkerSize(0.4);
 pZbb1->SetLineColor(216); 
 pZbb1->SetFillColor(216);
 pZbb1->SetFillStyle(3001);
 pZbb1->SetMarkerColor(216);
 pZbb1->SetMarkerStyle(22);
 pZbb1->SetMarkerSize(0.4);
 pZbb2->SetLineColor(215); 
 pZbb2->SetFillColor(215);
 pZbb2->SetFillStyle(3009);
 pZbb2->SetMarkerColor(215);
 pZbb2->SetMarkerStyle(22);
 pZbb2->SetMarkerSize(0.6);
 pZbb3->SetLineColor(214); 
 pZbb3->SetFillColor(214);
 pZbb3->SetFillStyle(3015);
 pZbb3->SetMarkerColor(214);
 pZbb3->SetMarkerStyle(22);
 pZbb3->SetMarkerSize(0.8);

 pZZ->SetLineColor(65); 
 pZZ->SetFillColor(65);
 // pZZ->SetFillStyle(3001);
 pZZ->SetMarkerColor(65);
 pZZ->SetMarkerStyle(23);
 // pZZ->SetMarkerSize(0.4);
 pZZ1->SetLineColor(67); 
 pZZ1->SetFillColor(67);
 pZZ1->SetFillStyle(3001);
 pZZ1->SetMarkerColor(67);
 pZZ1->SetMarkerStyle(23);
 pZZ1->SetMarkerSize(0.4);
 pZZ2->SetLineColor(65); 
 pZZ2->SetFillColor(65);
 pZZ2->SetFillStyle(3009);
 pZZ2->SetMarkerColor(65);
 pZZ2->SetMarkerStyle(23);
 pZZ2->SetMarkerSize(0.6);
 pZZ3->SetLineColor(63); 
 pZZ3->SetFillColor(63);
 pZZ3->SetFillStyle(3015);
 pZZ3->SetMarkerColor(63);
 pZZ3->SetMarkerStyle(23);
 pZZ3->SetMarkerSize(0.8);

 pWZ->SetLineColor(98); 
 pWZ->SetFillColor(98);
 // pWZ->SetFillStyle(3001);
 pWZ->SetMarkerColor(98);
 pWZ->SetMarkerStyle(29);
 // pWZ->SetMarkerSize(0.4);
 pWZ1->SetLineColor(96); 
 pWZ1->SetFillColor(96);
 pWZ1->SetFillStyle(3001);
 pWZ1->SetMarkerColor(96);
 pWZ1->SetMarkerStyle(29);
 pWZ1->SetMarkerSize(0.4);
 pWZ2->SetLineColor(98); 
 pWZ2->SetFillColor(98);
 pWZ2->SetFillStyle(3009);
 pWZ2->SetMarkerColor(98);
 pWZ2->SetMarkerStyle(29);
 pWZ2->SetMarkerSize(0.6);
 pWZ3->SetLineColor(100); 
 pWZ3->SetFillColor(100);
 pWZ3->SetFillStyle(3015);
 pWZ3->SetMarkerColor(100);
 pWZ3->SetMarkerStyle(29);
 pWZ3->SetMarkerSize(0.8);

 lLegend0->AddEntry(pH1->Clone(),pH1->GetTitle(),"lfp");
 lLegend0->AddEntry(pH2->Clone(),pH2->GetTitle(),"lfp");
 lLegend0->AddEntry(pH3->Clone(),pH3->GetTitle(),"lfp");
 lLegend0->AddEntry(pH4->Clone(),pH4->GetTitle(),"lfp");

 lLegend1->AddEntry(pZZ->Clone(), pZZ->GetTitle(), "lfp");
 lLegend1->AddEntry(pZZ1->Clone(),pZZ1->GetTitle(),"lfp");
 lLegend1->AddEntry(pZZ2->Clone(),pZZ2->GetTitle(),"lfp");
 lLegend1->AddEntry(pZZ3->Clone(),pZZ3->GetTitle(),"lfp");

 lLegend2->AddEntry(pZbb->Clone(), pZbb->GetTitle(), "lfp");
 lLegend2->AddEntry(pZbb1->Clone(),pZbb1->GetTitle(),"lfp");
 lLegend2->AddEntry(pZbb2->Clone(),pZbb2->GetTitle(),"lfp");
 lLegend2->AddEntry(pZbb3->Clone(),pZbb3->GetTitle(),"lfp");

 lLegend3->AddEntry(pWZ->Clone(), pWZ->GetTitle(), "lfp");
 lLegend3->AddEntry(pWZ1->Clone(),pWZ1->GetTitle(),"lfp");
 lLegend3->AddEntry(pWZ2->Clone(),pWZ2->GetTitle(),"lfp");
 lLegend3->AddEntry(pWZ3->Clone(),pWZ3->GetTitle(),"lfp");

 lLegend4->AddEntry(ptt, ptt->GetTitle(), "lfp");
 lLegend4->AddEntry(ptt1,ptt1->GetTitle(),"lfp");
 lLegend4->AddEntry(ptt2,ptt2->GetTitle(),"lfp");
 lLegend4->AddEntry(ptt3,ptt3->GetTitle(),"lfp");


 TCanvas *lC0 = new TCanvas("temp0","temp0",1000,600);
 lC0->Divide(5,1);
 lC0->cd(1);
 lLegend0->Draw();
 lC0->cd(2);
 lLegend1->Draw();
 lC0->cd(3);
 lLegend2->Draw();
 lC0->cd(4);
 lLegend3->Draw();
 lC0->cd(5);
 lLegend4->Draw();

 lC0->Print("COLORscheme.pdf");
// lC0->Print("COLORscheme.jpg");
// lC0->Print("COLORscheme.eps");
// lC0->Print("COLORscheme.ps");


}

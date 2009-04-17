#include <string>
#include <vector>
#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TObjString.h>
#include <TPRegexp.h>
#include <TCanvas.h>
#include <TH1F.h>

using namespace std;

TObject* get(TFile & file, const string & name) {
  TObject* object = file.Get( name.c_str() );
  if (object == 0) {
    string message = "ERROR: cannot find object " + name;
    cerr << message << endl;
    throw( message );
  } 
  return object;
}

vector<string> get_list(TFile& file) {
  TPRegexp pcreJetTag("JetTag_(.*)BJetTags_GLOBAL");

  vector<string> list;
  for (int i = 0; i < file.GetListOfKeys()->GetSize(); ++i) {
    TString name( file.GetListOfKeys()->At(i)->GetName() );
    if (pcreJetTag.MatchB( name ))
      list.push_back( string( ((TObjString *) pcreJetTag.MatchS( name )->At(1))->GetString() ) );
  }

  return list;
}

void list(TFile & file) {
  vector<string> list = get_list(file);

  cout << "BJetTags" << endl;
  for (unsigned int i = 0; i < list.size(); i++)
      cout << '\t' << list[i] << endl;
  cout << endl;
}

enum eta_range {
  GLOBAL  = 0,
  BARREL  = 1,
  ENDCAP  = 2,
  FORWARD = 2
};

enum energy_range {
  PT_GLOBAL = 0,
  PT_50_80  = 1,
  PT_80_120 = 2
};

const char* global = "_GLOBAL";
const char* eta_range_tag[]            = { "", "_ETA_0-1.4",  "_ETA_1.4-2.4" };
const char* eta_range_description[]    = { "", " (barrel)",   " (forward)"   };
const char* energy_range_tag[]         = { "", "_PT_50-80",   "_PT_80-120"   };
const char* energy_range_description[] = { "", " (pT 50-80)", " (pT 80-120)" };

enum flavour_color {
  colorBottom = kRed,
  colorCharm  = kBlue,
  colorLight  = kGreen,
  colorGluon  = kBlack
};

void plot(TFile & file, const string & name, eta_range eta = GLOBAL, energy_range energy = PT_GLOBAL, bool keep = false) {
  string tag;
  if (eta != GLOBAL || energy != PT_GLOBAL) {
    tag += eta_range_tag[eta];
    tag += energy_range_tag[energy];
  } else {
    tag = global;
  }

  string title_tag;
  title_tag += eta_range_description[eta];
  title_tag += energy_range_description[energy];

  string folder = "JetTag_" + name + "BJetTags" + tag;
  string name_discriminant  = name + "_discriminant";
  string name_efficiency    = name + "_efficiency";
  string name_mistag        = name + "_mistag";
  string title_discriminant = name + ": discriminant by flavour"          + title_tag;
  string title_efficiency   = name + ": efficiency vs. discriminator cut" + title_tag;
  string title_mistag       = name + ": mistag vs. b tag efficiency"      + title_tag;
  string name_b_discr = folder + "/" + "discr_"                   + name + "BJetTags" + tag + "B";
  string name_c_discr = folder + "/" + "discr_"                   + name + "BJetTags" + tag + "C";
  string name_x_discr = folder + "/" + "discr_"                   + name + "BJetTags" + tag + "DUS";
  string name_g_discr = folder + "/" + "discr_"                   + name + "BJetTags" + tag + "G";
  string name_b_eff   = folder + "/" + "effVsDiscrCut_discr_"     + name + "BJetTags" + tag + "B";
  string name_c_eff   = folder + "/" + "effVsDiscrCut_discr_"     + name + "BJetTags" + tag + "C";
  string name_x_eff   = folder + "/" + "effVsDiscrCut_discr_"     + name + "BJetTags" + tag + "DUS";
  string name_g_eff   = folder + "/" + "effVsDiscrCut_discr_"     + name + "BJetTags" + tag + "G";
  string name_c_vs_b  = folder + "/" + "FlavEffVsBEff_C_discr_"   + name + "BJetTags" + tag;
  string name_x_vs_b  = folder + "/" + "FlavEffVsBEff_DUS_discr_" + name + "BJetTags" + tag;
  string name_g_vs_b  = folder + "/" + "FlavEffVsBEff_G_discr_"   + name + "BJetTags" + tag;

  // added by mia
  string name_jetEta      = name + "_jetEta";
  string name_jetMomentum = name + "_jetMomentum";
  string name_jetPt       = name + "_jetPt";
  string title_jetEta      = name + ": jet eta"      + title_tag;
  string title_jetMomentum = name + ": jet momentum" + title_tag;
  string title_jetPt       = name + ": jet pt"       + title_tag;

  string name_b_jetEta = folder + "/" + "jetEta_" + name + "BJetTags" + tag + "B";
  string name_c_jetEta = folder + "/" + "jetEta_" + name + "BJetTags" + tag + "C";
  string name_x_jetEta = folder + "/" + "jetEta_" + name + "BJetTags" + tag + "DUS";
  string name_g_jetEta = folder + "/" + "jetEta_" + name + "BJetTags" + tag + "G";

  string name_b_jetMomentum = folder + "/" + "jetMomentum_" + name + "BJetTags" + tag + "B";
  string name_c_jetMomentum = folder + "/" + "jetMomentum_" + name + "BJetTags" + tag + "C";
  string name_x_jetMomentum = folder + "/" + "jetMomentum_" + name + "BJetTags" + tag + "DUS";
  string name_g_jetMomentum = folder + "/" + "jetMomentum_" + name + "BJetTags" + tag + "G";

  string name_b_jetPt = folder + "/" + "jetPt_" + name + "BJetTags" + tag + "B";
  string name_c_jetPt = folder + "/" + "jetPt_" + name + "BJetTags" + tag + "C";
  string name_x_jetPt = folder + "/" + "jetPt_" + name + "BJetTags" + tag + "DUS";
  string name_g_jetPt = folder + "/" + "jetPt_" + name + "BJetTags" + tag + "G";
  
  // discriminant distribution, by flavour
  TCanvas* discriminant = new TCanvas(name_discriminant.c_str(), title_discriminant.c_str());
  float max = 0;

  TH1F* plot_b_discr = (TH1F*) get( file, name_b_discr );
  if (plot_b_discr->Integral() > 0.0)
    plot_b_discr->Scale( 1.0 / plot_b_discr->Integral() );
  plot_b_discr->SetLineColor(colorBottom);
  plot_b_discr->SetLineWidth(2);
  if (plot_b_discr->GetMaximum() > max) max = plot_b_discr->GetMaximum();

  TH1F* plot_c_discr = (TH1F*) get( file, name_c_discr );
  if (plot_c_discr->Integral() > 0.0)
    plot_c_discr->Scale( 1.0 / plot_c_discr->Integral() );
  plot_c_discr->SetLineColor(colorCharm);
  plot_c_discr->SetLineWidth(2);
  if (plot_c_discr->GetMaximum() > max) max = plot_c_discr->GetMaximum();

  TH1F* plot_x_discr = (TH1F*) get( file, name_x_discr );
  if (plot_x_discr->Integral() > 0.0)
    plot_x_discr->Scale( 1.0 / plot_x_discr->Integral() );
  plot_x_discr->SetLineColor(colorLight);
  plot_x_discr->SetLineWidth(2);
  if (plot_x_discr->GetMaximum() > max) max = plot_x_discr->GetMaximum();

  TH1F* plot_g_discr = (TH1F*) get( file, name_g_discr );
  if (plot_g_discr->Integral() > 0.0)
    plot_g_discr->Scale( 1.0 / plot_g_discr->Integral() );
  plot_g_discr->SetLineColor(colorGluon);
  plot_g_discr->SetLineWidth(2);
  if (plot_g_discr->GetMaximum() > max) max = plot_g_discr->GetMaximum();

  plot_b_discr->SetMinimum(0.);
  plot_b_discr->SetMaximum(max * 1.1);
  plot_b_discr->SetTitle(title_discriminant.c_str());
  plot_b_discr->Draw("");
  plot_c_discr->Draw("same");
  plot_x_discr->Draw("same");
  plot_g_discr->Draw("same");
  discriminant->SetFillColor( kWhite );
  discriminant->SetGridx( true );
  discriminant->SetGridy( true );
  discriminant->SetFrameBorderMode( 0 );
  discriminant->SetFrameLineWidth( 2 );
  discriminant->SetCanvasSize(1280, 960);
  discriminant->SaveAs((name_discriminant + ".png").c_str());
  discriminant->SaveAs((name_discriminant + ".eps").c_str());
  discriminant->SetCanvasSize(400, 300);
  discriminant->SaveAs((name_discriminant + "_small.png").c_str());
  discriminant->SaveAs((name_discriminant + "_small.eps").c_str());

  // efficiency vs. discriminator cut, by flavour
  TCanvas* efficiency = new TCanvas(name_efficiency.c_str(), title_efficiency.c_str());

  TH1F* plot_b_eff = (TH1F*) get( file, name_b_eff );
  plot_b_eff->SetMarkerColor(colorBottom);
  plot_b_eff->SetMarkerSize(0.2);
  plot_b_eff->SetMarkerStyle(kFullDotMedium);
  if (plot_b_eff->GetMaximum() > max) max = plot_b_eff->GetMaximum();

  TH1F* plot_c_eff = (TH1F*) get( file, name_c_eff );
  plot_c_eff->SetMarkerColor(colorCharm);
  plot_c_eff->SetMarkerSize(0.2);
  plot_c_eff->SetMarkerStyle(kFullDotMedium);
  if (plot_c_eff->GetMaximum() > max) max = plot_c_eff->GetMaximum();

  TH1F* plot_x_eff = (TH1F*) get( file, name_x_eff );
  plot_x_eff->SetMarkerColor(colorLight);
  plot_x_eff->SetMarkerSize(0.2);
  plot_x_eff->SetMarkerStyle(kFullDotMedium);
  if (plot_x_eff->GetMaximum() > max) max = plot_x_eff->GetMaximum();

  TH1F* plot_g_eff = (TH1F*) get( file, name_g_eff );
  plot_g_eff->SetMarkerColor(colorGluon);
  plot_g_eff->SetMarkerSize(0.2);
  plot_g_eff->SetMarkerStyle(kFullDotMedium);
  if (plot_g_eff->GetMaximum() > max) max = plot_g_eff->GetMaximum();

  plot_b_eff->SetMinimum(0.);
  if (plot_b_eff->GetMaximum() < 0.2) {
    // soft lepton algorithms
    plot_b_eff->SetMaximum(0.2);
  } else {
    // other algorithms
    plot_b_eff->SetMaximum(1.);
  }
  plot_b_eff->SetTitle(title_efficiency.c_str());
  plot_b_eff->Draw("");
  plot_c_eff->Draw("same");
  plot_x_eff->Draw("same");
  plot_g_eff->Draw("same");
  efficiency->SetFillColor( kWhite );
  efficiency->SetGridx( true );
  efficiency->SetGridy( true );
  efficiency->SetFrameBorderMode( 0 );
  efficiency->SetFrameLineWidth( 2 );
  efficiency->SetCanvasSize(1280, 960);
  efficiency->SaveAs((name_efficiency + ".png").c_str());
  efficiency->SaveAs((name_efficiency + ".eps").c_str());
  efficiency->SetCanvasSize(400, 300);
  efficiency->SaveAs((name_efficiency + "_small.png").c_str());
  efficiency->SaveAs((name_efficiency + "_small.eps").c_str());

  // mistag vs. efficency plot
  TCanvas* mistag = new TCanvas(name_mistag.c_str(), title_mistag.c_str());
  mistag->SetLogy(true);

  TH1F* plot_c_vs_b = (TH1F*) get( file, name_c_vs_b );
  plot_c_vs_b->SetMarkerColor(colorCharm);
  plot_c_vs_b->SetMarkerSize(0.2);
  plot_c_vs_b->SetMarkerStyle(kFullDotMedium);
  
  TH1F* plot_x_vs_b = (TH1F*) get( file, name_x_vs_b );
  plot_x_vs_b->SetMarkerColor(colorLight);
  plot_x_vs_b->SetMarkerSize(0.2);
  plot_x_vs_b->SetMarkerStyle(kFullDotMedium);
  
  TH1F* plot_g_vs_b = (TH1F*) get( file, name_g_vs_b );
  plot_g_vs_b->SetMarkerColor(colorGluon);
  plot_g_vs_b->SetMarkerSize(0.2);
  plot_g_vs_b->SetMarkerStyle(kFullDotMedium);
  
  plot_c_vs_b->SetMinimum(0.0001);
  plot_c_vs_b->SetMaximum(1.);
  plot_c_vs_b->SetTitle(title_mistag.c_str());
  plot_c_vs_b->Draw("");
  plot_x_vs_b->Draw("same");
  plot_g_vs_b->Draw("same");
  
  mistag->SetFillColor( kWhite );
  mistag->SetGridx( true );
  mistag->SetGridy( true );
  mistag->SetFrameBorderMode( 0 );
  mistag->SetFrameLineWidth( 2 );
  mistag->SetCanvasSize(1280, 960);
  mistag->SaveAs((name_mistag + ".png").c_str());
  mistag->SaveAs((name_mistag + ".eps").c_str());
  mistag->SetCanvasSize(400, 300);
  mistag->SaveAs((name_mistag + "_small.png").c_str());
  mistag->SaveAs((name_mistag + "_small.eps").c_str());

  // jet eta distribution, by flavour
  TCanvas* jetEta = new TCanvas(name_jetEta.c_str(), title_jetEta.c_str());
  float max = 0;

  TH1F* plot_b_jetEta = (TH1F*) get( file, name_b_jetEta );
  if (plot_b_jetEta->Integral() > 0.0)
    plot_b_jetEta->Scale( 1.0 / plot_b_jetEta->Integral() );
  plot_b_jetEta->SetLineColor(colorBottom);
  plot_b_jetEta->SetLineWidth(2);
  if (plot_b_jetEta->GetMaximum() > max) max = plot_b_jetEta->GetMaximum();

  TH1F* plot_c_jetEta = (TH1F*) get( file, name_c_jetEta );
  if (plot_c_jetEta->Integral() > 0.0)
    plot_c_jetEta->Scale( 1.0 / plot_c_jetEta->Integral() );
  plot_c_jetEta->SetLineColor(colorCharm);
  plot_c_jetEta->SetLineWidth(2);
  if (plot_c_jetEta->GetMaximum() > max) max = plot_c_jetEta->GetMaximum();

  TH1F* plot_x_jetEta = (TH1F*) get( file, name_x_jetEta );
  if (plot_x_jetEta->Integral() > 0.0)
    plot_x_jetEta->Scale( 1.0 / plot_x_jetEta->Integral() );
  plot_x_jetEta->SetLineColor(colorLight);
  plot_x_jetEta->SetLineWidth(2);
  if (plot_x_jetEta->GetMaximum() > max) max = plot_x_jetEta->GetMaximum();

  TH1F* plot_g_jetEta = (TH1F*) get( file, name_g_jetEta );
  if (plot_g_jetEta->Integral() > 0.0)
    plot_g_jetEta->Scale( 1.0 / plot_g_jetEta->Integral() );
  plot_g_jetEta->SetLineColor(colorGluon);
  plot_g_jetEta->SetLineWidth(2);
  if (plot_g_jetEta->GetMaximum() > max) max = plot_g_jetEta->GetMaximum();

  plot_b_jetEta->SetMinimum(0.);
  plot_b_jetEta->SetMaximum(max * 1.1);
  plot_b_jetEta->SetTitle(title_jetEta.c_str());
  plot_b_jetEta->Draw("");
  plot_c_jetEta->Draw("same");
  plot_x_jetEta->Draw("same");
  plot_g_jetEta->Draw("same");
  jetEta->SetFillColor( kWhite );
  jetEta->SetGridx( true );
  jetEta->SetGridy( true );
  jetEta->SetFrameBorderMode( 0 );
  jetEta->SetFrameLineWidth( 2 );
  jetEta->SetCanvasSize(1280, 960);
  jetEta->SaveAs((name_jetEta + ".png").c_str());
  jetEta->SaveAs((name_jetEta + ".eps").c_str());
  jetEta->SetCanvasSize(400, 300);
  jetEta->SaveAs((name_jetEta + "_small.png").c_str());
  jetEta->SaveAs((name_jetEta + "_small.eps").c_str());

  // jet momentum distribution, by flavour
  TCanvas* jetMomentum = new TCanvas(name_jetMomentum.c_str(), title_jetMomentum.c_str());
  float max = 0;

  TH1F* plot_b_jetMomentum = (TH1F*) get( file, name_b_jetMomentum );
  if (plot_b_jetMomentum->Integral() > 0.0)
    plot_b_jetMomentum->Scale( 1.0 / plot_b_jetMomentum->Integral() );
  plot_b_jetMomentum->SetLineColor(colorBottom);
  plot_b_jetMomentum->SetLineWidth(2);
  if (plot_b_jetMomentum->GetMaximum() > max) max = plot_b_jetMomentum->GetMaximum();

  TH1F* plot_c_jetMomentum = (TH1F*) get( file, name_c_jetMomentum );
  if (plot_c_jetMomentum->Integral() > 0.0)
    plot_c_jetMomentum->Scale( 1.0 / plot_c_jetMomentum->Integral() );
  plot_c_jetMomentum->SetLineColor(colorCharm);
  plot_c_jetMomentum->SetLineWidth(2);
  if (plot_c_jetMomentum->GetMaximum() > max) max = plot_c_jetMomentum->GetMaximum();

  TH1F* plot_x_jetMomentum = (TH1F*) get( file, name_x_jetMomentum );
  if (plot_x_jetMomentum->Integral() > 0.0)
    plot_x_jetMomentum->Scale( 1.0 / plot_x_jetMomentum->Integral() );
  plot_x_jetMomentum->SetLineColor(colorLight);
  plot_x_jetMomentum->SetLineWidth(2);
  if (plot_x_jetMomentum->GetMaximum() > max) max = plot_x_jetMomentum->GetMaximum();

  TH1F* plot_g_jetMomentum = (TH1F*) get( file, name_g_jetMomentum );
  if (plot_g_jetMomentum->Integral() > 0.0)
    plot_g_jetMomentum->Scale( 1.0 / plot_g_jetMomentum->Integral() );
  plot_g_jetMomentum->SetLineColor(colorGluon);
  plot_g_jetMomentum->SetLineWidth(2);
  if (plot_g_jetMomentum->GetMaximum() > max) max = plot_g_jetMomentum->GetMaximum();

  plot_b_jetMomentum->SetMinimum(0.);
  plot_b_jetMomentum->SetMaximum(max * 1.1);
  plot_b_jetMomentum->SetTitle(title_jetMomentum.c_str());
  plot_b_jetMomentum->Draw("");
  plot_c_jetMomentum->Draw("same");
  plot_x_jetMomentum->Draw("same");
  plot_g_jetMomentum->Draw("same");
  jetMomentum->SetFillColor( kWhite );
  jetMomentum->SetGridx( true );
  jetMomentum->SetGridy( true );
  jetMomentum->SetFrameBorderMode( 0 );
  jetMomentum->SetFrameLineWidth( 2 );
  jetMomentum->SetCanvasSize(1280, 960);
  jetMomentum->SaveAs((name_jetMomentum + ".png").c_str());
  jetMomentum->SaveAs((name_jetMomentum + ".eps").c_str());
  jetMomentum->SetCanvasSize(400, 300);
  jetMomentum->SaveAs((name_jetMomentum + "_small.png").c_str());
  jetMomentum->SaveAs((name_jetMomentum + "_small.eps").c_str());

  // jet pt distribution, by flavour
  TCanvas* jetPt = new TCanvas(name_jetPt.c_str(), title_jetPt.c_str());
  float max = 0;

  TH1F* plot_b_jetPt = (TH1F*) get( file, name_b_jetPt );
  if (plot_b_jetPt->Integral() > 0.0)
    plot_b_jetPt->Scale( 1.0 / plot_b_jetPt->Integral() );
  plot_b_jetPt->SetLineColor(colorBottom);
  plot_b_jetPt->SetLineWidth(2);
  if (plot_b_jetPt->GetMaximum() > max) max = plot_b_jetPt->GetMaximum();

  TH1F* plot_c_jetPt = (TH1F*) get( file, name_c_jetPt );
  if (plot_c_jetPt->Integral() > 0.0)
    plot_c_jetPt->Scale( 1.0 / plot_c_jetPt->Integral() );
  plot_c_jetPt->SetLineColor(colorCharm);
  plot_c_jetPt->SetLineWidth(2);
  if (plot_c_jetPt->GetMaximum() > max) max = plot_c_jetPt->GetMaximum();

  TH1F* plot_x_jetPt = (TH1F*) get( file, name_x_jetPt );
  if (plot_x_jetPt->Integral() > 0.0)
    plot_x_jetPt->Scale( 1.0 / plot_x_jetPt->Integral() );
  plot_x_jetPt->SetLineColor(colorLight);
  plot_x_jetPt->SetLineWidth(2);
  if (plot_x_jetPt->GetMaximum() > max) max = plot_x_jetPt->GetMaximum();

  TH1F* plot_g_jetPt = (TH1F*) get( file, name_g_jetPt );
  if (plot_g_jetPt->Integral() > 0.0)
    plot_g_jetPt->Scale( 1.0 / plot_g_jetPt->Integral() );
  plot_g_jetPt->SetLineColor(colorGluon);
  plot_g_jetPt->SetLineWidth(2);
  if (plot_g_jetPt->GetMaximum() > max) max = plot_g_jetPt->GetMaximum();

  plot_b_jetPt->SetMinimum(0.);
  plot_b_jetPt->SetMaximum(max * 1.1);
  plot_b_jetPt->SetTitle(title_jetPt.c_str());
  plot_b_jetPt->Draw("");
  plot_c_jetPt->Draw("same");
  plot_x_jetPt->Draw("same");
  plot_g_jetPt->Draw("same");
  jetPt->SetFillColor( kWhite );
  jetPt->SetGridx( true );
  jetPt->SetGridy( true );
  jetPt->SetFrameBorderMode( 0 );
  jetPt->SetFrameLineWidth( 2 );
  jetPt->SetCanvasSize(1280, 960);
  jetPt->SaveAs((name_jetPt + ".png").c_str());
  jetPt->SaveAs((name_jetPt + ".eps").c_str());
  jetPt->SetCanvasSize(400, 300);
  jetPt->SaveAs((name_jetPt + "_small.png").c_str());
  jetPt->SaveAs((name_jetPt + "_small.eps").c_str());

  if (! keep) {
    delete discriminant;
    delete efficiency;
    delete mistag;

    delete jetEta;
    delete jetMomentum;
    delete jetPt;
  }
}

void make_all_plots(TFile & file) {
  vector<string> list = get_list(file);
  
  for (unsigned int i = 0; i < list.size(); i++) {
    // softElectron btag is defined only in the barrel region
    if (list[i] == "softElectron")
      plot( file, list[i], BARREL, PT_GLOBAL );
    else
      plot( file, list[i], GLOBAL, PT_GLOBAL );
  }
}

// convenience functions with implicit file 
//#ifdef __CINT__

TFile * _file0;

void list(void) {
  list( *_file0 );
}

void plot(const string & name, eta_range eta = GLOBAL, energy_range energy = PT_GLOBAL, bool keep = false) {
  plot( *_file0, name, eta, energy, keep );
}

void make_all_plots() {
  make_all_plots( *_file0 );
}

//#endif // __CINT__

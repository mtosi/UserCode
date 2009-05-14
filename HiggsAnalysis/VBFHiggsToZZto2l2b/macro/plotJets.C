/*
   created by mia tosi on the Sep 1st, 2008
    - macro skeleton for the jets analysis:
      - input file w/ one of its branches [offlineJets]
      - loop over the events
      - loop over the jets collection per event
      - get jets information per jets
      - fill some histograms on jets information
      - write histograms into the output file
      - save output file
*/

#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TFile.h"
#include "TStyle.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

// This is needed to avoid errors from CINT
#if !defined(__CINT__) && !defined(__MAKECINT__)

// Classes to be stored
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/BaseJet.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/BaseMEt.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/OfflineMEt.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/OfflineJet.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/MCParticle.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/GlobalMuon.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/SimpleElectron.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/SimpleTau.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/Summary.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/PythiaParticleIndex.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/ParticlesMass.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/ParticlesCharge.h"

#endif

int plotJets ( int sample = 0 ) {

  using namespace std;
  using namespace edm;
  using namespace vbfhzz2l2b;

  // general root setting
  gROOT->Reset(); 
  //  gROOT->SetBatch(kTRUE);
  gStyle->SetTitleSize(0.1);
  gStyle->SetLabelSize(0.04,"X");
  gStyle->SetLabelSize(0.04,"Y");
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  // setting overflow, underflow and other options
  gStyle->SetOptStat("nemrou");
  gStyle->SetOptFit(111);

  gStyle->SetPalette(1);
  const Int_t NCont = 255;
  /*
  const Int_t NRGBs = 5;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  gStyle->CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  */
  gStyle->SetNumberContours(NCont);

  const double PI_= TMath::Pi();
  
  TChain events("Events");
  // set the buffers for the branches
  std::vector<GlobalMuon> gloMuon_vec;
  TBranch* gloMuon_B;
  std::vector<OfflineJet> offJets_vec;
  TBranch* offJets_B;


  TString sampleName = "";
  TString castorSubDirectory = "";
  switch (sample) {
  case(1):
    sampleName = "VBFH160ZZ";
    castorSubDirectory = "HIGGS_ZZ/";
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS160_ZZ_1.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS160_ZZ_2.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS160_ZZ_3.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS160_ZZ_4.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS160_ZZ_5.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS160_ZZ_6.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS160_ZZ_7.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS160_ZZ_8.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS160_ZZ_9.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS160_ZZ_10.root");
    //    events.SetBranchAddress("vbfhzz2l2bGlobalMuons_offlineProd_globalMuons_PROD.obj",&gloMuon_vec,&gloMuon_B);
    events.SetBranchAddress("vbfhzz2l2bOfflineJets_offlineProd_offlineJets_PROD.obj",&offJets_vec,&offJets_B);
    break;
  case(2):
    sampleName = "VBFH200ZZ";
    castorSubDirectory = "HIGGS_ZZ/";
    events.Add("eventRootFile/VBFHIGGS200_ZZ_1.root");
    events.Add("eventRootFile/VBFHIGGS200_ZZ_2.root");
    events.Add("eventRootFile/VBFHIGGS200_ZZ_3.root");
    events.Add("eventRootFile/VBFHIGGS200_ZZ_4.root");
    events.Add("eventRootFile/VBFHIGGS200_ZZ_5.root");
    events.Add("eventRootFile/VBFHIGGS200_ZZ_6.root");
    events.Add("eventRootFile/VBFHIGGS200_ZZ_7.root");
    events.Add("eventRootFile/VBFHIGGS200_ZZ_8.root");
    events.Add("eventRootFile/VBFHIGGS200_ZZ_9.root");
    events.Add("eventRootFile/VBFHIGGS200_ZZ_10.root");
    //    events.SetBranchAddress("vbfhzz2l2bGlobalMuons_offlineProd_globalMuons_PROD.obj",&gloMuon_vec,&gloMuon_B);
    events.SetBranchAddress("vbfhzz2l2bOfflineJets_offlineProd_offlineJets_PROD.obj",&offJets_vec,&offJets_B);
    break;
  case(3):
    sampleName = "VBFH400ZZ";
    castorSubDirectory = "HIGGS_ZZ/";
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS400_ZZ_1.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS400_ZZ_2.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS400_ZZ_3.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS400_ZZ_4.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS400_ZZ_5.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS400_ZZ_6.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS400_ZZ_7.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS400_ZZ_8.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS400_ZZ_9.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS400_ZZ_10.root");
    //    events.SetBranchAddress("vbfhzz2l2bGlobalMuons_offlineProd_globalMuons_PROD.obj",&gloMuon_vec,&gloMuon_B);
    events.SetBranchAddress("vbfhzz2l2bOfflineJets_offlineProd_offlineJets_PROD.obj",&offJets_vec,&offJets_B);
    break;
  case(4):
    sampleName = "VBFH800ZZ";
    castorSubDirectory = "HIGGS_ZZ/";
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS800_ZZ_1.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS800_ZZ_2.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS800_ZZ_3.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS800_ZZ_4.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS800_ZZ_5.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS800_ZZ_6.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS800_ZZ_7.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS800_ZZ_8.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS800_ZZ_9.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS800_ZZ_10.root");
    //    events.SetBranchAddress("vbfhzz2l2bGlobalMuons_offlineProd_globalMuons_PROD.obj",&gloMuon_vec,&gloMuon_B);
    events.SetBranchAddress("vbfhzz2l2bOfflineJets_offlineProd_offlineJets_PROD.obj",&offJets_vec,&offJets_B);
    break;
  case(10):
    sampleName = "ZZ_0JETS";
    castorSubDirectory = sampleName+"/";
    // input files
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/ZZ_0JETS/ZZ_0JETS_1.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/ZZ_0JETS/ZZ_0JETS_2.root");
    //    events.SetBranchAddress("vbfhzz2l2bGlobalMuons_offlineProd_globalMuons_PROD1.obj",&gloMuon_vec,&gloMuon_B);
    events.SetBranchAddress("vbfhzz2l2bOfflineJets_offlineProd_offlineJets_PROD1.obj",&offJets_vec,&offJets_B);
    break;
  case(11):
    sampleName = "ZZ_1JETS";
    castorSubDirectory = sampleName+"/";
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/ZZ_1JETS/ZZ_1JETS_1.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/ZZ_1JETS/ZZ_1JETS_2.root");
    //    events.SetBranchAddress("vbfhzz2l2bGlobalMuons_offlineProd_globalMuons_PROD1.obj",&gloMuon_vec,&gloMuon_B);
    events.SetBranchAddress("vbfhzz2l2bOfflineJets_offlineProd_offlineJets_PROD1.obj",&offJets_vec,&offJets_B);
    break;
  case(12):
    sampleName = "ZZ_2JETS";
    castorSubDirectory = sampleName+"/";
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/ZZ_2JETS/ZZ_2JETS_1.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/ZZ_2JETS/ZZ_2JETS_2.root");
    //    events.SetBranchAddress("vbfhzz2l2bGlobalMuons_offlineProd_globalMuons_PROD1.obj",&gloMuon_vec,&gloMuon_B);
    events.SetBranchAddress("vbfhzz2l2bOfflineJets_offlineProd_offlineJets_PROD1.obj",&offJets_vec,&offJets_B);
    break;
  case(13):
    sampleName = "ZZ_3JETS";
    castorSubDirectory = sampleName+"/";
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/ZZ_3JETS/ZZ_3JETS_1.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/ZZ_3JETS/ZZ_3JETS_2.root");
    //    events.SetBranchAddress("vbfhzz2l2bGlobalMuons_offlineProd_globalMuons_PROD1.obj",&gloMuon_vec,&gloMuon_B);
    events.SetBranchAddress("vbfhzz2l2bOfflineJets_offlineProd_offlineJets_PROD1.obj",&offJets_vec,&offJets_B);
    break;
  default:
    sampleName = "VBFH200ZZ";
    castorSubDirectory = "HIGGS_ZZ/";
    events.Add("eventRootFile/VBFHIGGS200_ZZ_1.root");
    events.Add("eventRootFile/VBFHIGGS200_ZZ_2.root");
    events.Add("eventRootFile/VBFHIGGS200_ZZ_3.root");
    events.Add("eventRootFile/VBFHIGGS200_ZZ_4.root");
    events.Add("eventRootFile/VBFHIGGS200_ZZ_5.root");
    events.Add("eventRootFile/VBFHIGGS200_ZZ_6.root");
    events.Add("eventRootFile/VBFHIGGS200_ZZ_7.root");
    events.Add("eventRootFile/VBFHIGGS200_ZZ_8.root");
    events.Add("eventRootFile/VBFHIGGS200_ZZ_9.root");
    events.Add("eventRootFile/VBFHIGGS200_ZZ_10.root");
    //    events.SetBranchAddress("vbfhzz2l2bGlobalMuons_offlineProd_globalMuons_PROD.obj",&gloMuon_vec,&gloMuon_B);
    events.SetBranchAddress("vbfhzz2l2bOfflineJets_offlineProd_offlineJets_PROD.obj",&offJets_vec,&offJets_B);
  }


  // output file
  TFile *outputfile;
  outputfile = dynamic_cast<TFile*>(gROOT->FindObject(sampleName+"plotJets.root")); 
  if (outputfile) outputfile->Close();
  outputfile = new TFile(sampleName+"plotJets.root","RECREATE","jets histograms");

  // jets information histogram
  TH1F * histo_jetsNumber      = new TH1F("jetsNumber",     "number of jets",                        50,     0.,   50. ); 
  TH1F * histo_jetE            = new TH1F("jetE",           "jet energy",                           520,     0. ,5200. ); 
  TH1F * histo_jetEMfrac       = new TH1F("jetEMfrac",      "jet EM energy fraction",                22,     0.,    1.1);
  TH1F * histo_jetUncorrEt     = new TH1F("jetUncorrEt",    "jet uncorrected transverse energy",    300,     0.,  300. );
  TH1F * histo_jetEt           = new TH1F("jetEt",          "jet corrected transverse energy",      300,     0.,  300. );
  TH1F * histo_jetEx           = new TH1F("jetEx",          "jet energy projection along x-axis",   500,  -500.,  500. ); 
  TH1F * histo_jetEy           = new TH1F("jetEy",          "jet energy projection along y-axis",   500,  -500. , 500. ); 
  TH1F * histo_jetEz           = new TH1F("jetEz",          "jet energy projection along z-axis",  2000, -2000. ,2000. ); 
  TH1F * histo_jetEta          = new TH1F("jetEta",         "jet pseudo-rapidity",                  200,    -5.,    5. );
  TH1F * histo_jetPhi          = new TH1F("jetPhi",         "jet azimutal angle",                   120,   -PI_,    PI_);
  TH1F * histo_jetMass         = new TH1F("jetMass",        "jet mass",                             200,     0.,  100. );
  TH1F * histo_jetTagTkMassS1  = new TH1F("jetTagTkMassS1", "tag jet: invariant mass w/ S1 tracks",  50,     0.,   10. );
  TH1F * histo_jetTagTkMassS2  = new TH1F("jetTagTkMassS2", "tag jet: invariant mass w/ S2 tracks",  50,     0.,   10. );
  TH1F * histo_jetTagTkMassS3  = new TH1F("jetTagTkMassS3", "tag jet: invariant mass w/ S3 tracks",  50,     0.,   10. );
  TH1F * histo_jetTagTkNumS1   = new TH1F("jetTagTknums1",  "tag jet: number of S1 tracks",          30,     0.,   15. );
  TH1F * histo_jetTagTkNumS2   = new TH1F("jetTagTknums2",  "tag jet: number of S2 tracks",          30,     0.,   15. );
  TH1F * histo_jetTagTkNumS3   = new TH1F("jetTagTknums3",  "tag jet: number of S3 tracks",          30,     0.,   15. );
  TH1F * histo_jetTagTkSumPtS1 = new TH1F("jetTagTkSumPtS1","tag jet: #Sigma p_{T} of S1 tracks",    700,     0.,  350. );
  TH1F * histo_jetTagTkSumPtS2 = new TH1F("jetTagTkSumPtS2","tag jet: #Sigma p_{T} of S2 tracks",    700,     0.,  350. );
  TH1F * histo_jetTagTkSumPtS3 = new TH1F("jetTagTkSumPtS3","tag jet: #Sigma p_{T} of S3 tracks",    700,     0.,  350. );

  TH2F * histo_jetEtaVSEt       = new TH2F("jetEtaVSEt",       "jet #eta VS transverse energys",         120, 0.,60.,200,-5.,5.);
  TH2F * histo_jetEtaVSTagMass1 = new TH2F("jetEtaVSTagMass1", "jet #eta VS invariant mass w/ S1 tracks", 50, 0.,10.,200,-5.,5.);
  TH2F * histo_jetEtaVSTagMass2 = new TH2F("jetEtaVSTagMass2", "jet #eta VS invariant mass w/ S2 tracks", 50, 0.,10.,200,-5.,5.);
  TH2F * histo_jetEtaVSTagMass3 = new TH2F("jetEtaVSTagMass3", "jet #eta VS invariant mass w/ S3 tracks", 50, 0.,10.,200,-5.,5.);

  TH1F * histo_jetsEt25Number  = new TH1F("jetsEt25Number", "fraction of jet w/ corrected E_{T} #ge 25 GeV",500,0.,1.1);
  TH1F * histo_jetsEt30Number  = new TH1F("jetsEt30Number", "fraction of jet w/ corrected E_{T} #ge 30 GeV",500,0.,1.1);
  TH1F * histo_jetsEt35Number  = new TH1F("jetsEt35Number", "fraction of jet w/ corrected E_{T} #ge 35 GeV",500,0.,1.1);
  TH1F * histo_jetsEt40Number  = new TH1F("jetsEt40Number", "fraction of jet w/ corrected E_{T} #ge 40 GeV",500,0.,1.1);
  TH1F * histo_jetsEt45Number  = new TH1F("jetsEt45Number", "fraction of jet w/ corrected E_{T} #ge 45 GeV",500,0.,1.1);
  TH1F * histo_jetsE100Number  = new TH1F("jetsE100Number", "fraction of jet w/ E #ge 100 GeV",             500,0.,1.1);
  TH1F * histo_jetsEta15Number = new TH1F("jetsEta15Number","fraction of jet w/ |#eta| #le 1.5",            500,0.,1.1);

  TH2F * histo_jetsEtVSjetsNumber        = new TH2F("jetsEtVSjetsNumber",       "jet E_{T} VS number of jets",                                     100, 0.,100.,50,0.,50.);
  TH2F * histo_jetMaxEtVSjetsNumber      = new TH2F("jetMaxEtVSjetsNumber",     "maximum jet E_{T} VS number of jets",                             100, 0.,100.,50,0.,50.);
  TH2F * histo_jetMaxEtVSjetsEt25Number  = new TH2F("jetMaxEtVSjetsEt25Number", "maximum jet E_{T} VS number of jets w/ corrected E_{T} >= 25 GeV",100, 0.,100.,50,0.,50.);
  TH2F * histo_jetMaxEtVSjetsEt30Number  = new TH2F("jetMaxEtVSjetsEt30Number", "maximum jet E_{T} VS number of jets w/ corrected E_{T} >= 30 GeV",100, 0.,100.,50,0.,50.);
  TH2F * histo_jetMaxEtVSjetsEt35Number  = new TH2F("jetMaxEtVSjetsEt35Number", "maximum jet E_{T} VS number of jets w/ corrected E_{T} >= 35 GeV",100, 0.,100.,50,0.,50.);
  TH2F * histo_jetMaxEtVSjetsEt40Number  = new TH2F("jetMaxEtVSjetsEt40Number", "maximum jet E_{T} VS number of jets w/ corrected E_{T} >= 40 GeV",100, 0.,100.,50,0.,50.);
  TH2F * histo_jetMaxEtVSjetsEt45Number  = new TH2F("jetMaxEtVSjetsEt45Number", "maximum jet E_{T} VS number of jets w/ corrected E_{T} >= 45 GeV",100, 0.,100.,50,0.,50.);
  TH2F * histo_jetMaxEtVSjetsE100Number  = new TH2F("jetMaxEtVSjetsE100Number", "maximum jet E_{T} VS number of jets w/ corrected E >= 100 GeV",   100, 0.,100.,50,0.,50.);
  TH2F * histo_jetMaxEtVSjetsEta15Number = new TH2F("jetMaxEtVSjetsEta15Number","maximum jet E_{T} VS number of jets w/ corrected |#eta|< 1.5",    100, 0.,100.,50,0.,50.);
  TH2F * histo_jetMinEtVSjetsNumber      = new TH2F("jetMinEtVSjetsNumber",     "minimum jet E_{T} VS number of jets",                             50, 0.,50.,50,0.,50.);
  TH2F * histo_jetMinEtVSjetsEt25Number  = new TH2F("jetMinEtVSjetsEt25Number", "minimum jet E_{T} VS number of jets w/ corrected E_{T} >= 25 GeV",50, 0.,50.,50,0.,50.);
  TH2F * histo_jetMinEtVSjetsEt30Number  = new TH2F("jetMinEtVSjetsEt30Number", "minimum jet E_{T} VS number of jets w/ corrected E_{T} >= 30 GeV",50, 0.,50.,50,0.,50.);
  TH2F * histo_jetMinEtVSjetsEt35Number  = new TH2F("jetMinEtVSjetsEt35Number", "minimum jet E_{T} VS number of jets w/ corrected E_{T} >= 35 GeV",50, 0.,50.,50,0.,50.);
  TH2F * histo_jetMinEtVSjetsEt40Number  = new TH2F("jetMinEtVSjetsEt40Number", "minimum jet E_{T} VS number of jets w/ corrected E_{T} >= 40 GeV",50, 0.,50.,50,0.,50.);
  TH2F * histo_jetMinEtVSjetsEt45Number  = new TH2F("jetMinEtVSjetsEt45Number", "minimum jet E_{T} VS number of jets w/ corrected E_{T} >= 45 GeV",50, 0.,50.,50,0.,50.);
  TH2F * histo_jetMinEtVSjetsE100Number  = new TH2F("jetMinEtVSjetsE100Number", "minimum jet E_{T} VS number of jets w/ corrected E >= 100 GeV",   50, 0.,50.,50,0.,50.);
  TH2F * histo_jetMinEtVSjetsEta15Number = new TH2F("jetMinEtVSjetsEta15Number","minimum jet E_{T} VS number of jets w/ corrected |#eta|< 1.5",    50, 0.,50.,50,0.,50.);

  TH1F * histo_jetsVtxZ  = new TH1F("jetsVtxZ", "jets longitudinal vertex position",200,-20.,20.);
  TH2F * histo_jetsVtxXY = new TH2F("jetsVtxXY","jets transverse vertex position",  300,-0.15,0.15,300,-0.15,0.15);

  TH2F * histo_jetFirstMaxEtaVSjetSecondMaxEta = new TH2F("jetFirstMaxEtaVSjetSecondMaxEta","",200,-5.,5.,100,0.,5.);
  // global muons
  math::XYZTLorentzVector null_XYZTLorentzVector = math::XYZTLorentzVector(0.,0.,0.,0.);
  GlobalMuon null_globalMuon = GlobalMuon( 0.,0.,0.,0,
					   null_XYZTLorentzVector,
					   math::XYZPoint(0.,0.,0.),
					   0.,0.,0.,0.,0.,0.,0.,0,0,0.,0,0.,0.);
  OfflineJet null_offlinejet = OfflineJet( 0.,0.,0.,0.,0.,
					   null_XYZTLorentzVector,
					   math::XYZPoint(0.,0.,0.),
					   0.,0.,0.,0, 0., 0.,0, 0., 0.,0, 0., 0. );
  unsigned int numberOfEvents = events.GetEntries();

  int jets_number = 0;
  unsigned int et25Counter  = 0;
  unsigned int et30Counter  = 0;
  unsigned int et35Counter  = 0;
  unsigned int et40Counter  = 0;
  unsigned int et45Counter  = 0;
  unsigned int e100Counter  = 0;
  unsigned int eta15Counter = 0;
      
  // loop over the events
  for( unsigned int evt_iter = 0;
       evt_iter < 1000;
       //       evt_iter < numberOfEvents;
       ++evt_iter) {

    // need to call SetAddress since TBranch's change for each file read
    offJets_B->SetAddress(&offJets_vec);
    offJets_B->GetEntry(evt_iter);
    events.GetEntry(evt_iter,0);

    // now can access data
    if ( evt_iter/100 == float(evt_iter)/100. ) std::cout << "event = " << evt_iter << std::endl;
    jets_number = offJets_vec.size();
    
    if ( jets_number != 0 ) {
      //      if ( evt_iter/100 == float(evt_iter)/100. ) std::cout << "number of jets = " << jets_number << std::endl;

      et25Counter  = 0;
      et30Counter  = 0;
      et35Counter  = 0;
      et40Counter  = 0;
      et45Counter  = 0;
      e100Counter  = 0;
      eta15Counter = 0;
      double jetMaxEt   = -99.;
      double jetMinEt   =  999.;
      double jetMinEt25 =  999.;
      double jetMinEt30 =  999.;
      double jetMinEt35 =  999.;
      double jetMinEt40 =  999.;
      double jetMinEt45 =  999.;
      double eta_max[2] = {0.,0.};

      // loop on all the jets
      for ( unsigned int offJets_iter = 0; offJets_iter < offJets_vec.size(); ++offJets_iter ) { 
	// take the jet
	OfflineJet * offlineJets = &( offJets_vec[offJets_iter] );
	
	math::XYZTLorentzVector jet_p4       = offlineJets->p4();
	math::XYZPoint          jet_vertex   = offlineJets->vertex();
	double                  jet_e        = offlineJets->e();
	double                  jet_EMfrac   = offlineJets->emEnergyFraction(); 
        double                  jet_uncorrEt = offlineJets->uncorrEt();
	double                  jet_et       = offlineJets->et(); 
	double                  jet_ex       = offlineJets->ex();
	double                  jet_ey       = offlineJets->ey();
	double                  jet_ez       = offlineJets->ez();
	double                  jet_eta      = offlineJets->eta();
	double                  jet_phi      = offlineJets->phi();
	double                  jet_mass     = offlineJets->jetMass();
	float                   jet_discriminatorHE = offlineJets->discriminatorHighEff();
        float                   jet_discriminatorHP = offlineJets->discriminatorHighPur();
	double                  jet_tagTkMassS1  = offlineJets->tagTkMassS1(); 
        double                  jet_tagTkMassS2  = offlineJets->tagTkMassS2();
        double                  jet_tagTkMassS3  = offlineJets->tagTkMassS3();
        int                     jet_tagTkNumS1   = offlineJets->tkNumS1();
        int                     jet_tagTkNumS2   = offlineJets->tkNumS2();
        int                     jet_tagTkNumS3   = offlineJets->tkNumS3();
        double                  jet_tagTkSumPtS1 = offlineJets->tkSumPtS1();
        double                  jet_tagTkSumPtS2 = offlineJets->tkSumPtS2();
        double                  jet_tagTkSumPtS3 = offlineJets->tkSumPtS3();

	if (jet_et  >=  25. ) { et25Counter++; if (jet_et <= jetMinEt) jetMinEt25 = jet_et;}
	if (jet_et  >=  30. ) { et30Counter++; if (jet_et <= jetMinEt) jetMinEt30 = jet_et;}
	if (jet_et  >=  35. ) { et35Counter++; if (jet_et <= jetMinEt) jetMinEt35 = jet_et;}
	if (jet_et  >=  40. ) { et40Counter++; if (jet_et <= jetMinEt) jetMinEt40 = jet_et;}
	if (jet_et  >=  45. ) { et45Counter++; if (jet_et <= jetMinEt) jetMinEt45 = jet_et;}
	if (jet_e   >= 100. ) e100Counter++;
	if (jet_eta >=   1.5) eta15Counter++;

	if (jet_et >= jetMaxEt) jetMaxEt = jet_et;
	if (jet_et <= jetMinEt) jetMinEt = jet_et;


	if (fabs(jet_eta) >= fabs(eta_max[0])) eta_max[0]=jet_eta;
	else if(fabs(jet_eta) >= fabs(eta_max[1])) eta_max[1]=jet_eta;

	// fill jets information per jets
        histo_jetE->Fill(jet_e);
        histo_jetEMfrac->Fill(jet_EMfrac);
	histo_jetUncorrEt->Fill(jet_uncorrEt);
	histo_jetEt->Fill(jet_et);
	histo_jetEx->Fill(jet_ex);
        histo_jetEy->Fill(jet_ey);
	histo_jetEz->Fill(jet_ez);
	histo_jetEta->Fill(jet_eta);
	histo_jetPhi->Fill(jet_phi);
	if ( jet_mass         != 0. ) histo_jetMass->Fill(jet_mass);
	if ( jet_tagTkMassS1  != 0. ) histo_jetTagTkMassS1->Fill(jet_tagTkMassS1);
	if ( jet_tagTkMassS2  != 0. ) histo_jetTagTkMassS2->Fill(jet_tagTkMassS2);
	if ( jet_tagTkMassS3  != 0. ) histo_jetTagTkMassS3->Fill(jet_tagTkMassS3);
	if ( jet_tagTkNumS1   != 0  ) histo_jetTagTkNumS1->Fill(jet_tagTkNumS1);
	if ( jet_tagTkNumS2   != 0  ) histo_jetTagTkNumS2->Fill(jet_tagTkNumS2);
	if ( jet_tagTkNumS3   != 0  ) histo_jetTagTkNumS3->Fill(jet_tagTkNumS3);
	if ( jet_tagTkSumPtS1 != 0. ) histo_jetTagTkSumPtS1->Fill(jet_tagTkSumPtS1);
	if ( jet_tagTkSumPtS2 != 0. ) histo_jetTagTkSumPtS2->Fill(jet_tagTkSumPtS2);
	if ( jet_tagTkSumPtS3 != 0. ) histo_jetTagTkSumPtS3->Fill(jet_tagTkSumPtS3);

	histo_jetEtaVSEt->Fill(jet_et,jet_eta);
	if ( jet_tagTkMassS1  != 0. ) histo_jetEtaVSTagMass1->Fill(jet_tagTkMassS1,jet_eta);	
	if ( jet_tagTkMassS2  != 0. ) histo_jetEtaVSTagMass2->Fill(jet_tagTkMassS2,jet_eta);	
	if ( jet_tagTkMassS3  != 0. ) histo_jetEtaVSTagMass3->Fill(jet_tagTkMassS3,jet_eta);	

	histo_jetsEtVSjetsNumber->Fill(jets_number,jet_et);
	histo_jetsVtxZ->Fill(jet_vertex.z());
	histo_jetsVtxXY->Fill(jet_vertex.x(),jet_vertex.y());

      } // end loop on jets
      histo_jetsEt25Number->Fill(double(et25Counter)/double(jets_number));
      histo_jetsEt30Number->Fill(double(et30Counter)/double(jets_number));
      histo_jetsEt35Number->Fill(double(et35Counter)/double(jets_number));
      histo_jetsEt40Number->Fill(double(et40Counter)/double(jets_number));
      histo_jetsEt45Number->Fill(double(et45Counter)/double(jets_number));
      histo_jetsE100Number->Fill(double(e100Counter)/double(jets_number));  
      histo_jetsEta15Number->Fill(double(eta15Counter)/double(jets_number)); 
    
      histo_jetMaxEtVSjetsNumber->Fill(     jets_number, jetMaxEt);
      histo_jetMaxEtVSjetsEt25Number->Fill( et25Counter, jetMaxEt);
      histo_jetMaxEtVSjetsEt30Number->Fill( et30Counter, jetMaxEt);
      histo_jetMaxEtVSjetsEt35Number->Fill( et35Counter, jetMaxEt);
      histo_jetMaxEtVSjetsEt40Number->Fill( et40Counter, jetMaxEt);
      histo_jetMaxEtVSjetsEt45Number->Fill( et45Counter, jetMaxEt);
      histo_jetMaxEtVSjetsE100Number->Fill( e100Counter, jetMaxEt);
      histo_jetMaxEtVSjetsEta15Number->Fill(eta15Counter,jetMaxEt);
      histo_jetMinEtVSjetsNumber->Fill(     jets_number, jetMinEt);
      histo_jetMinEtVSjetsEt25Number->Fill( et25Counter, jetMinEt25);
      histo_jetMinEtVSjetsEt30Number->Fill( et30Counter, jetMinEt30);
      histo_jetMinEtVSjetsEt35Number->Fill( et35Counter, jetMinEt35);
      histo_jetMinEtVSjetsEt40Number->Fill( et40Counter, jetMinEt40);
      histo_jetMinEtVSjetsEt45Number->Fill( et45Counter, jetMinEt45);
      histo_jetMinEtVSjetsE100Number->Fill( e100Counter, jetMinEt);
      histo_jetMinEtVSjetsEta15Number->Fill(eta15Counter,jetMinEt);

      if (eta_max[0]*eta_max[1] < 0) histo_jetFirstMaxEtaVSjetSecondMaxEta->Fill(-fabs(eta_max[0]),fabs(eta_max[1]));
      else histo_jetFirstMaxEtaVSjetSecondMaxEta->Fill(fabs(eta_max[0]),fabs(eta_max[1]));
      
    } // end if (jets_number != 0)
    
    // fill jets information per event
    histo_jetsNumber->Fill(jets_number);
  } // end loop on events
  /*
  TCanvas * offJets_canvas1 = new TCanvas("c1","c1",1000,800);
  offJets_canvas1->cd();
  histo_jetsNumber->GetXaxis()->SetTitle("jets number");
  histo_jetsNumber->Draw();
  
  TCanvas * offJets_canvas2 = new TCanvas("c2", "c2", 1000, 800);
  offJets_canvas2->Divide(2,1);
  offJets_canvas2->cd(1);
  histo_jetEt->GetXaxis()->SetTitle("e_{T} (GeV)");
  histo_jetEt->Draw();
  offJets_canvas2->cd(2);
  histo_jetUncorrEt->GetXaxis()->SetTitle("uncorrected e_{T} (GeV)");
  histo_jetUncorrEt->Draw();


  TCanvas * offJets_canvas3 = new TCanvas("c3", "c3", 1000, 800);
  offJets_canvas3->Divide(2,2);
  offJets_canvas3->cd(1);
  histo_jetEx->GetXaxis()->SetTitle("e_{x} (GeV)");
  histo_jetEx->Draw();
  offJets_canvas3->cd(2);
  histo_jetEy->GetXaxis()->SetTitle("e_{y} (GeV)");
  histo_jetEy->Draw();
  offJets_canvas3->cd(3);
  histo_jetEz->GetXaxis()->SetTitle("e_{z} (GeV)");
  histo_jetEz->Draw();

  TCanvas * offJets_canvas4 = new TCanvas("c4", "c4", 1000, 800);
  offJets_canvas4->Divide(2,2);
  offJets_canvas4->cd(1);
  histo_jetMass->GetXaxis()->SetTitle("jet mass (GeV/c^{2})");
  histo_jetMass->Draw();
  offJets_canvas4->cd(2);
  histo_jetPhi->GetXaxis()->SetTitle("#phi (rad)");
  histo_jetPhi->Draw();
  offJets_canvas4->cd(4);
  histo_jetEta->GetXaxis()->SetTitle("#eta");
  histo_jetEta->Draw();


  TCanvas * offJets_canvas5 = new TCanvas("c5", "c5", 1000, 800);
  offJets_canvas5->Divide(2,2);
  offJets_canvas5->cd(1);
  histo_jetTagTkMassS1->GetXaxis()->SetTitle("tagTkMassS1");
  histo_jetTagTkMassS1->Draw();
  offJets_canvas5->cd(2);
  histo_jetTagTkMassS2->GetXaxis()->SetTitle("tagTkMassS2");
  histo_jetTagTkMassS2->Draw();
  offJets_canvas5->cd(3);
  histo_jetTagTkMassS3->GetXaxis()->SetTitle("tagTkMassS3");
  histo_jetTagTkMassS3->Draw();
*/
  TCanvas * offJets_canvas5 = new TCanvas("c5", "c5", 1000, 800);
  offJets_canvas5->Divide(4,2);
  offJets_canvas5->cd(1);
  histo_jetsNumber->GetXaxis()->SetTitle("fraction of jets");
  histo_jetsNumber->GetXaxis()->SetTitleSize(0.04);
  histo_jetsNumber->Draw();
  offJets_canvas5->cd(2);
  histo_jetsEt25Number->GetXaxis()->SetTitle("fraction of jets");
  histo_jetsEt25Number->GetXaxis()->SetTitleSize(0.04);
  histo_jetsEt25Number->Draw();
  offJets_canvas5->cd(3);
  histo_jetsEt30Number->GetXaxis()->SetTitle("fraction of jets");
  histo_jetsEt30Number->GetXaxis()->SetTitleSize(0.04);
  histo_jetsEt30Number->Draw();
  offJets_canvas5->cd(4);
  histo_jetsEt35Number->GetXaxis()->SetTitle("fraction of jets");
  histo_jetsEt35Number->GetXaxis()->SetTitleSize(0.04);
  histo_jetsEt35Number->Draw();
  offJets_canvas5->cd(5);
  histo_jetsEt40Number->GetXaxis()->SetTitle("fraction of jets");
  histo_jetsEt40Number->GetXaxis()->SetTitleSize(0.04);
  histo_jetsEt40Number->Draw();
  offJets_canvas5->cd(6);
  histo_jetsEt45Number->GetXaxis()->SetTitle("fraction of jets");
  histo_jetsEt45Number->GetXaxis()->SetTitleSize(0.04);
  histo_jetsEt45Number->Draw();
  offJets_canvas5->cd(7);
  histo_jetsE100Number->GetXaxis()->SetTitle("fraction of jets");
  histo_jetsE100Number->GetXaxis()->SetTitleSize(0.04);
  histo_jetsE100Number->Draw();  
  offJets_canvas5->cd(8);
  histo_jetsEta15Number->GetXaxis()->SetTitle("fraction of jets");
  histo_jetsEta15Number->GetXaxis()->SetTitleSize(0.04);
  histo_jetsEta15Number->Draw();

  TCanvas * offJets_canvas6 = new TCanvas("c6", "c6", 1000, 800);
  offJets_canvas6->Divide(2,2);
  offJets_canvas6->cd(1);
  histo_jetEtaVSEt->GetYaxis()->SetTitle("#eta ");
  histo_jetEtaVSEt->GetXaxis()->SetTitle("E_{T}(GeV)");
  histo_jetEtaVSEt->SetStats(kFALSE);
  histo_jetEtaVSEt->Draw("colz");
  offJets_canvas6->cd(2);  
  histo_jetEtaVSTagMass1->GetYaxis()->SetTitle("#eta ");
  histo_jetEtaVSTagMass1->GetXaxis()->SetTitle("invariant mass w/ S1 tracks (GeV/c^{2})");
  histo_jetEtaVSTagMass1->SetStats(kFALSE);
  histo_jetEtaVSTagMass1->Draw("colz");
  offJets_canvas6->cd(3);
  histo_jetEtaVSTagMass2->GetYaxis()->SetTitle("#eta ");
  histo_jetEtaVSTagMass2->GetXaxis()->SetTitle("invariant mass w/ S2 tracks (GeV/c^{2})");
  histo_jetEtaVSTagMass2->SetStats(kFALSE);
  histo_jetEtaVSTagMass2->Draw("colz");
  offJets_canvas6->cd(4);
  histo_jetEtaVSTagMass3->GetYaxis()->SetTitle("#eta ");
  histo_jetEtaVSTagMass3->GetXaxis()->SetTitle("invariant mass w/ S3 tracks (GeV/c^{2})");
  histo_jetEtaVSTagMass3->SetStats(kFALSE);
  histo_jetEtaVSTagMass3->Draw("colz");

  TCanvas * offJets_canvas7 = new TCanvas("c7", "c7", 1000, 800);
  offJets_canvas7->Divide(2,2);
  offJets_canvas7->cd(1);
  histo_jetsEtVSjetsNumber->GetXaxis()->SetTitle("number of jets");
  histo_jetsEtVSjetsNumber->GetXaxis()->SetTitleSize(0.04);
  histo_jetsEtVSjetsNumber->GetYaxis()->SetTitle("E_{T} (GeV)");
  histo_jetsEtVSjetsNumber->GetYaxis()->SetTitleSize(0.04);
  histo_jetsEtVSjetsNumber->Draw();
  offJets_canvas7->cd(3);
  histo_jetMaxEtVSjetsNumber->GetXaxis()->SetTitle("number of jets");
  histo_jetMaxEtVSjetsNumber->GetXaxis()->SetTitleSize(0.04);
  histo_jetMaxEtVSjetsNumber->GetYaxis()->SetTitle("E_{T}^{max} (GeV)");
  histo_jetMaxEtVSjetsNumber->GetYaxis()->SetTitleSize(0.04);
  histo_jetMaxEtVSjetsNumber->Draw();
  offJets_canvas7->cd(4);
  histo_jetMinEtVSjetsNumber->GetXaxis()->SetTitle("number of jets");
  histo_jetMinEtVSjetsNumber->GetXaxis()->SetTitleSize(0.04);
  histo_jetMinEtVSjetsNumber->GetYaxis()->SetTitle("E_{T}^{min} (GeV)");
  histo_jetMinEtVSjetsNumber->GetYaxis()->SetTitleSize(0.04);
  histo_jetMinEtVSjetsNumber->Draw();

  TCanvas * offJets_canvas8 = new TCanvas("c8", "c8", 1000, 800);
  offJets_canvas8->Divide(3,3);
  offJets_canvas8->cd(1);
  histo_jetMaxEtVSjetsEt25Number->GetXaxis()->SetTitle("number of jets");
  histo_jetMaxEtVSjetsEt25Number->GetXaxis()->SetTitleSize(0.04);
  histo_jetMaxEtVSjetsEt25Number->GetYaxis()->SetTitle("E_{T} (GeV)");
  histo_jetMaxEtVSjetsEt25Number->GetYaxis()->SetTitleSize(0.04);
  histo_jetMaxEtVSjetsEt25Number->Draw();
  histo_jetMinEtVSjetsEt25Number->SetMarkerColor(kRed);
  histo_jetMinEtVSjetsEt25Number->Draw("same");
  offJets_canvas8->cd(2);
  histo_jetMaxEtVSjetsEt30Number->GetXaxis()->SetTitle("number of jets");
  histo_jetMaxEtVSjetsEt30Number->GetXaxis()->SetTitleSize(0.04);
  histo_jetMaxEtVSjetsEt30Number->GetYaxis()->SetTitle("E_{T} (GeV)");
  histo_jetMaxEtVSjetsEt30Number->GetYaxis()->SetTitleSize(0.04);
  histo_jetMaxEtVSjetsEt30Number->Draw();
  histo_jetMinEtVSjetsEt30Number->SetMarkerColor(kRed);
  histo_jetMinEtVSjetsEt30Number->Draw("same");
  offJets_canvas8->cd(3);
  histo_jetMaxEtVSjetsEt35Number->GetXaxis()->SetTitle("number of jets");
  histo_jetMaxEtVSjetsEt35Number->GetXaxis()->SetTitleSize(0.04);
  histo_jetMaxEtVSjetsEt35Number->GetYaxis()->SetTitle("E_{T} (GeV)");
  histo_jetMaxEtVSjetsEt35Number->GetYaxis()->SetTitleSize(0.04);
  histo_jetMaxEtVSjetsEt35Number->Draw();
  histo_jetMinEtVSjetsEt35Number->SetMarkerColor(kRed);
  histo_jetMinEtVSjetsEt35Number->Draw("same");
  offJets_canvas8->cd(4);
  histo_jetMaxEtVSjetsEt40Number->GetXaxis()->SetTitle("number of jets");
  histo_jetMaxEtVSjetsEt40Number->GetXaxis()->SetTitleSize(0.04);
  histo_jetMaxEtVSjetsEt40Number->GetYaxis()->SetTitle("E_{T} (GeV)");
  histo_jetMaxEtVSjetsEt40Number->GetYaxis()->SetTitleSize(0.04);
  histo_jetMaxEtVSjetsEt40Number->Draw();
  histo_jetMinEtVSjetsEt40Number->SetMarkerColor(kRed);
  histo_jetMinEtVSjetsEt40Number->Draw("same");
  offJets_canvas8->cd(5);
  histo_jetMaxEtVSjetsEt45Number->GetXaxis()->SetTitle("number of jets");
  histo_jetMaxEtVSjetsEt45Number->GetXaxis()->SetTitleSize(0.04);
  histo_jetMaxEtVSjetsEt45Number->GetYaxis()->SetTitle("E_{T} (GeV)");
  histo_jetMaxEtVSjetsEt45Number->GetYaxis()->SetTitleSize(0.04);
  histo_jetMaxEtVSjetsEt45Number->Draw();
  histo_jetMinEtVSjetsEt45Number->SetMarkerColor(kRed);
  histo_jetMinEtVSjetsEt45Number->Draw("same");
  offJets_canvas8->cd(6);
  histo_jetMaxEtVSjetsE100Number->GetXaxis()->SetTitle("number of jets");
  histo_jetMaxEtVSjetsE100Number->GetXaxis()->SetTitleSize(0.04);
  histo_jetMaxEtVSjetsE100Number->GetYaxis()->SetTitle("E_{T} (GeV)");
  histo_jetMaxEtVSjetsE100Number->GetYaxis()->SetTitleSize(0.04);
  histo_jetMaxEtVSjetsE100Number->Draw();
  histo_jetMinEtVSjetsE100Number->SetMarkerColor(kRed);
  histo_jetMinEtVSjetsE100Number->Draw("same");
  offJets_canvas8->cd(7);
  histo_jetMaxEtVSjetsEta15Number->GetXaxis()->SetTitle("number of jets");
  histo_jetMaxEtVSjetsEta15Number->GetXaxis()->SetTitleSize(0.04);
  histo_jetMaxEtVSjetsEta15Number->GetYaxis()->SetTitle("E_{T} (GeV)");
  histo_jetMaxEtVSjetsEta15Number->GetYaxis()->SetTitleSize(0.04);
  histo_jetMaxEtVSjetsEta15Number->Draw();
  histo_jetMinEtVSjetsEt35Number->SetMarkerColor(kRed);
  histo_jetMinEtVSjetsEta15Number->Draw("same");

  TCanvas * offJets_canvas9 = new TCanvas("c9", "c9", 1000, 800);
  offJets_canvas9->Divide(2,2);
  offJets_canvas9->cd(1);
  histo_jetsVtxXY->GetXaxis()->SetTitle("vertex x (cm)");
  histo_jetsVtxXY->GetXaxis()->SetTitleSize(0.04);
  histo_jetsVtxXY->GetYaxis()->SetTitle("vertex y (cm)");
  histo_jetsVtxXY->GetYaxis()->SetTitleSize(0.04);
  histo_jetsVtxXY->Draw();
  offJets_canvas9->cd(2);
  histo_jetsVtxZ->GetXaxis()->SetTitle("vertex z (cm)");
  histo_jetsVtxZ->GetXaxis()->SetTitleSize(0.04);
  histo_jetsVtxZ->Draw();

  TCanvas * offJets_canvas10 = new TCanvas("c10", "c10", 1000, 800);
  offJets_canvas10->cd();
  histo_jetFirstMaxEtaVSjetSecondMaxEta->GetXaxis()->SetTitle("sign(#eta_{1^{st}},#eta_{2^{nd}})*|#eta_{1^{st}}|");
  histo_jetFirstMaxEtaVSjetSecondMaxEta->GetXaxis()->SetTitleSize(0.04);
  histo_jetFirstMaxEtaVSjetSecondMaxEta->GetYaxis()->SetTitle("|#eta_{2^{nd}}|");
  histo_jetFirstMaxEtaVSjetSecondMaxEta->GetYaxis()->SetTitleSize(0.04);
  histo_jetFirstMaxEtaVSjetSecondMaxEta->Draw();


  /*
  histo_jetsNumber->Write();
  histo_jetE->Write();
  histo_jetEMfrac->Write();
  histo_jetUncorrEt->Write();
  histo_jetEt->Write();
  histo_jetEx->Write();
  histo_jetEy->Write();
  histo_jetEz->Write();
  histo_jetEta->Write();
  histo_jetPhi->Write();
  histo_jetMass->Write();
  histo_jetTagTkMassS1->Write();
  histo_jetTagTkMassS2->Write();
  histo_jetTagTkMassS3->Write();
  histo_jetTagTkNumS1->Write();
  histo_jetTagTkNumS2->Write();
  histo_jetTagTkNumS3->Write();
  histo_jetTagTkSumPtS1->Write();
  histo_jetTagTkSumPtS2->Write();
  histo_jetTagTkSumPtS3->Write();
  histo_jetEtaVSEt->Write();
  histo_jetEtaVSTagMass1->Write();
  histo_jetEtaVSTagMass2->Write();
  histo_jetEtaVSTagMass3->Write();
  */
  
//  offJets_canvas1->Write();
//  offJets_canvas2->Write();
//  offJets_canvas3->Write();
//  offJets_canvas4->Write();
  offJets_canvas5->Write();
  offJets_canvas6->Write();
  offJets_canvas7->Write();
  offJets_canvas8->Write();
  offJets_canvas9->Write();
  offJets_canvas10->Write();

//  offJets_canvas1->Print("jetsCanvas1.eps");
//  offJets_canvas2->Print("jetsCanvas2.eps");
//  offJets_canvas3->Print("jetsCanvas3.eps");
//  offJets_canvas4->Print("jetsCanvas4.eps");
//  offJets_canvas5->Print("jetsCanvas5.eps");
//  offJets_canvas6->Print("jetsCanvas6.eps");
  gStyle->SetOptStat(0);
  offJets_canvas10->Print("fwdJetsEtaPlot.eps");

  outputfile->Write();

  return 0;
};

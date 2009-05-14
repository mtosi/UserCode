/*
  created by mia tosi on the Aug 26th, 2008
  - macro skeleton for the HZZ analysis:
  - input file w/ one of its branches [global muons]
  - loop over the events
  - loop over the muon collection per event
  - get muon information per muon
  - fill some histograms on muon information
  - write histograms into the output file
  - save output file
*/

#include "TF1.h"
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
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/BaseJet.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/BaseMEt.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/OfflineMEt.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/OfflineJet.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/MCParticle.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/GlobalMuon.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/SimpleElectron.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/SimpleTau.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/Summary.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/PythiaParticleIndex.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/ParticlesMass.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/ParticlesCharge.h"

#endif

int plotGlobalMuons ( int sample = 0 ) {

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
    events.SetBranchAddress("vbfhzz2l2bGlobalMuons_offlineProd_globalMuons_PROD.obj",&gloMuon_vec,&gloMuon_B);
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
    events.SetBranchAddress("vbfhzz2l2bGlobalMuons_offlineProd_globalMuons_PROD.obj",&gloMuon_vec,&gloMuon_B);
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
    events.SetBranchAddress("vbfhzz2l2bGlobalMuons_offlineProd_globalMuons_PROD.obj",&gloMuon_vec,&gloMuon_B);
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
    events.SetBranchAddress("vbfhzz2l2bGlobalMuons_offlineProd_globalMuons_PROD.obj",&gloMuon_vec,&gloMuon_B);
    break;
  case(10):
    sampleName = "ZZ_0JETS";
    castorSubDirectory = sampleName+"/";
    // input files
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/ZZ_0JETS/ZZ_0JETS_1.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/ZZ_0JETS/ZZ_0JETS_2.root");
    events.SetBranchAddress("vbfhzz2l2bGlobalMuons_offlineProd_globalMuons_PROD1.obj",&gloMuon_vec,&gloMuon_B);
    break;
  case(11):
    sampleName = "ZZ_1JETS";
    castorSubDirectory = sampleName+"/";
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/ZZ_1JETS/ZZ_1JETS_1.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/ZZ_1JETS/ZZ_1JETS_2.root");
    events.SetBranchAddress("vbfhzz2l2bGlobalMuons_offlineProd_globalMuons_PROD1.obj",&gloMuon_vec,&gloMuon_B);
    break;
  case(12):
    sampleName = "ZZ_2JETS";
    castorSubDirectory = sampleName+"/";
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/ZZ_2JETS/ZZ_2JETS_1.root");
    //    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/ZZ_2JETS/ZZ_2JETS_2.root");
    events.SetBranchAddress("vbfhzz2l2bGlobalMuons_offlineProd_globalMuons_PROD1.obj",&gloMuon_vec,&gloMuon_B);
    break;
  case(13):
    sampleName = "ZZ_3JETS";
    castorSubDirectory = sampleName+"/";
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/ZZ_3JETS/ZZ_3JETS_1.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/ZZ_3JETS/ZZ_3JETS_2.root");
    events.SetBranchAddress("vbfhzz2l2bGlobalMuons_offlineProd_globalMuons_PROD1.obj",&gloMuon_vec,&gloMuon_B);
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
    events.SetBranchAddress("vbfhzz2l2bGlobalMuons_offlineProd_globalMuons_PROD.obj",&gloMuon_vec,&gloMuon_B);
  }

  // output file
  TFile *outputfile;
  outputfile = dynamic_cast<TFile*>(gROOT->FindObject(sampleName+"plotGlobalMuons.root")); 
  if (outputfile) outputfile->Close();
  outputfile = new TFile(sampleName+"plotGlobalMuons.root","RECREATE","global muons histograms");

  // muon information histogram 
  TH1F * histo_muonsNumber   = new TH1F("muonsNumber",  "number of global muons",                                20,   0. ,   10.  );
  TH1F * histo_muonPt        = new TH1F("muonPt",       "global muon transverse momentum",                      150,   0. ,  200.  );
  TH1F * histo_muonPz        = new TH1F("muonPz",       "global muon momentum along z-axis",                    600,-300. ,  300.  );
  TH1F * histo_muonEta       = new TH1F("muonEta",      "global muon pseudo-rapidity",                          100,  -5. ,    5.  );
  TH1F * histo_muonPhi       = new TH1F("muonPhi",      "global muon azimutal angle",                   int(PI_*20),  -PI_,    PI_ );
  TH1F * histo_muonCharge    = new TH1F("muonCharge",   "global muon charge",                                    11,  -1.1,    1.1 ); 
  TH1F * histo_muonEM        = new TH1F("muonEM",       "global muon EM energy deposit",                         20,  -0.5,    0.5 );
  TH1F * histo_muonHAD       = new TH1F("muonHAD",      "global muon hadronic energy deposit",                   20,  -0.5,    0.5 );
  TH1F * histo_muonHO        = new TH1F("muonHO",       "global muon Ho energy deposit",                         20,  -0.5,    0.5 );
  TH1F * histo_muonIsoSumPt  = new TH1F("muonIsoSumPt", "global muon isolation: #Sigma p_{T}",                   20,  -0.5,    0.5 );
  TH1F * histo_muonIsoEMet   = new TH1F("muonIsoEMet",  "global muon isolation: EM transverse energy",           20,  -0.5,    0.5 );
  TH1F * histo_muonIsoHADet  = new TH1F("muonIsoHADet", "global muon isolation: hadronic transverse energy",     20,  -0.5,    0.5 );
  TH1F * histo_muonIsoHOet   = new TH1F("muonIsoHOet",  "global muon isolation: ho transverse energy",           20,  -0.5,    0.5 );
  TH1F * histo_muonIsoNjets  = new TH1F("muonIsoNjets", "global muon isolation: number of jets",                 20,  -0.5,    1.5 );
  TH1F * histo_muonIsoNtk    = new TH1F("muonIsonTk",   "global muon isolation: number of tracks",               20,  -0.5,    1.5 );
  TH1F * histo_muonNormChi2  = new TH1F("muonNormChi2", "global muon isolation: #chi^{2}/n.o.f.",               100,   0. ,   10.  );
  TH1F * histo_muonHitsNum   = new TH1F("muonHitsNum",  "global muon number of hits",                            40,   5.,    30.  );
  TH1F * histo_muonDxy       = new TH1F("muonDxy",      "global muon transverse impact parameter",              100,  -0.1,    0.1 );
  TH1F * histo_muonDxyError  = new TH1F("muonDxyError", "global muon transverse impact parameter error",        100,   0. ,    0.01);
  TH1F * histo_muonDxySignif = new TH1F("muonDxySignif","global muon transverse impact parameter significance", 200, -100. , 100.  );  


  TH1F * histo_muonsPt[4];
  TH1F * histo_muonsPt2[4];
  TH1F * histo_muonsDxy[4];
  TH1F * histo_muonsDxyError[4];
  TH1F * histo_muonsDxySignif[4];
  TH2F * histo_posMuPzVSnegMuPz[4];
  TH2F * histo_muonsPhiVSmuonsEta[4];
  TH1F * histo_dimuonDeltaPhi[4];
  TH1F * histo_dimuonDeltaEta[4];
  TH1F * histo_dimuonDeltaR[4];
  TH1F * histo_dimuonMass[4];
  TH1F * histo_dimuonPt[4];
  TH1F * histo_dimuonPhi[4];
  TH1F * histo_dimuonEta[4];
  TH2F * histo_muonsPtVSdimuonMass[4];
  TH2F * histo_muonsEtaVSdimuonMass[4];
  TH2F * histo_muonsDxySignifVSdimuonMass[4];
  TH2F * histo_dimuonDeltaPhiVSdimuonMass[4];
  TH2F * histo_dimuonDeltaEtaVSdimuonMass[4];
  TH2F * histo_dimuonDeltaRVSdimuonMass[4];
  TH2F * histo_dimuonPtVSdimuonMass[4];

  TH1F * histo_ZmuonsDxy[4];
  TH1F * histo_ZmuonsDxyError[4];
  TH1F * histo_ZmuonsDxySignif[4];
  TH2F * histo_ZposMuPzVSZnegMuPz[4];
  TH2F * histo_ZmuonsPhiVSZmuonsEta[4];
  TH1F * histo_ZdimuonDeltaPhi[4];
  TH1F * histo_ZdimuonDeltaEta[4];
  TH1F * histo_ZdimuonDeltaR[4];
  TH1F * histo_ZmuonsPt[4];
  TH1F * histo_ZmuonsPt2[4];
  TH1F * histo_ZMass[4];
  TH1F * histo_ZMass_pt10Cut[4];
  TH1F * histo_ZMass_pt15Cut[4];
  TH1F * histo_ZMass_pt20Cut[4];
  TH1F * histo_ZPt[4];
  TH1F * histo_ZPhi[4];
  TH1F * histo_ZEta[4];
  TH2F * histo_ZmuonsPtVSZMass[4];
  TH2F * histo_ZmuonsEtaVSZMass[4];
  TH2F * histo_ZmuonsDxySignifVSZMass[4];
  TH2F * histo_ZdimuonDeltaPhiVSZMass[4];
  TH2F * histo_ZdimuonDeltaEtaVSZMass[4];
  TH2F * histo_ZdimuonDeltaRVSZMass[4];
  TH2F * histo_ZPtVSZMass[4];
  char name[50];
  for ( unsigned int index = 0; index < 4; index++ ) {
    sprintf(name,"muonsPt_%d",       index); histo_muonsPt[index]        = new TH1F(name,"transverse momentum of the muons w/ opposite charge",                      200,   0. , 200.  );
    sprintf(name,"muonsPt2_%d",      index); histo_muonsPt2[index]       = new TH1F(name,"transverse momentum of the muons w/ opposite charge squared",              200,   0. ,4000.  );
    sprintf(name,"muonsDxy_%d",      index); histo_muonsDxy[index]       = new TH1F(name,"transverse impact parameter of the muons w/ opposite charge",              100,  -0.1,   0.1 );  
    sprintf(name,"muonsDxyError_%d", index); histo_muonsDxyError[index]  = new TH1F(name,"transverse impact parameter error of the muons w/ opposite charge",        100,   0.,    0.01);  
    sprintf(name,"muonsDxySignif_%d",index); histo_muonsDxySignif[index] = new TH1F(name,"transverse impact parameter significance of the muons w/ opposite charge", 200,-100.,  100.  );  
    sprintf(name,"muonPosPzVSmuonNegPz_%d", index); histo_posMuPzVSnegMuPz[index]     = new TH2F(name,"global muon VS global anti-muon momentum along z-axis",        600,-300.,300.,600,-300.,300.);
    sprintf(name,"muonsPhiVSmuonsEta_%d",   index); histo_muonsPhiVSmuonsEta[index]   = new TH2F(name,"muons azimutal angle VS muons pseudo-rapidity",        int(PI_*20), -PI_, PI_,100,  -5.,  5.);
    sprintf(name,"dimuonDeltaPhi_%d",index); histo_dimuonDeltaPhi[index] = new TH1F(name,"#Delta #phi between 2 muons w/ opposite charge",                   int(PI_*10),   0., PI_);
    sprintf(name,"dimuonDeltaEta_%d",index); histo_dimuonDeltaEta[index] = new TH1F(name,"#Delta #eta between 2 muons w/ opposite charge",                            50,   0.,  5.);
    sprintf(name,"dimuonDeltaR_%d",  index); histo_dimuonDeltaR[index]   = new TH1F(name,"#Delta R between 2 muons w/ opposite charge",                               50,   0.,  5.);
    sprintf(name,"dimuonMass_%d",    index); histo_dimuonMass[index]     = new TH1F(name,"invariant mass of 2 muons system w/ opposite charge",                      120,   0.,120.);
    sprintf(name,"dimuonPt_%d",      index); histo_dimuonPt[index]       = new TH1F(name,"transverse momentum of 2 muons system w/ opposite charge",                 200,   0.,200.);
    sprintf(name,"dimuonEta_%d",     index); histo_dimuonEta[index]      = new TH1F(name,"pseudo-rapidity of 2 muons system w/ opposite charge",                     100,  -5.,  5.);
    sprintf(name,"dimuonPhi_%d",     index); histo_dimuonPhi[index]      = new TH1F(name,"azimutal angle of 2 muons system w/ opposite charge",              int(PI_*20), -PI_, PI_);
    sprintf(name,"muonsPtVSdimuonMass_%d",       index); histo_muonsPtVSdimuonMass[index]        = new TH2F(name,"muons transverse momentum VS invariant mass of 2 muons system",
							 						    120,   0.,120.,        200,   0.,200.);
    sprintf(name,"muonsEtaVSdimuonMass_%d",      index); histo_muonsEtaVSdimuonMass[index]       = new TH2F(name,"muons pseudo-rapidity VS invariant mass of 2 muons system",
							 						    120,   0.,120.,        100,  -5.,  5.);
    sprintf(name,"muonsDxySignifVSdimuonMass_%d",index); histo_muonsDxySignifVSdimuonMass[index] = new TH2F(name,"muons transverse impact parameter significance VS invariant mass of 2 muons system",
							 						    120,   0.,120.,        200,-100.,100.);
    sprintf(name,"dimuonDeltaPhiVSdimuonMass_%d",index); histo_dimuonDeltaPhiVSdimuonMass[index] = new TH2F(name,"#Delta#phi VS invariant mass of 2 muons system w/ opposite charge",
							 						    120,   0.,120.,int(PI_*10),   0., PI_);
    sprintf(name,"dimuonDeltaEtaVSdimuonMass_%d",index); histo_dimuonDeltaEtaVSdimuonMass[index] = new TH2F(name,"#Delta#eta VS invariant mass of 2 muons system",
							 						    120,   0.,120.,         50,   0.,  5.);
    sprintf(name,"dimuonDeltaRVSdimuonMass_%d",  index); histo_dimuonDeltaRVSdimuonMass[index]   = new TH2F(name,"#DeltaR VS invariant mass of 2 muons system",
							 						    120,   0.,120.,         50,   0.,  5.);
    sprintf(name,"dimuonPtVSdimuonMass_%d",      index); histo_dimuonPtVSdimuonMass[index]       = new TH2F(name,"muons transverse momentum VS invariant mass of 2 muons system",
													    120,   0.,120.,        200,   0.,200.);

    sprintf(name,"ZmuonsPt_%d",       index); histo_ZmuonsPt[index]        = new TH1F(name,"transverse momentum of the muons from the reconstructed Z",                      200,   0. , 200.  );
    sprintf(name,"ZmuonsPt2_%d",      index); histo_ZmuonsPt2[index]       = new TH1F(name,"transverse momentum of the muons from the reconstructed Z squared",              200,   0. ,4000.  );
    sprintf(name,"ZmuonsDxy_%d",      index); histo_ZmuonsDxy[index]       = new TH1F(name,"transverse impact parameter of the muons from the reconstructed Z",              100,  -0.1,   0.1 );
    sprintf(name,"ZmuonsDxyErrror_%d",index); histo_ZmuonsDxyError[index]  = new TH1F(name,"transverse impact parameter error of the muons from the reconstructed Z",        100,   0. ,   0.01);
    sprintf(name,"ZmuonsDxySignif_%d",index); histo_ZmuonsDxySignif[index] = new TH1F(name,"transverse impact parameter significance of the muons from the reconstructed Z", 200,-100. , 100.  );  
    sprintf(name,"ZdimuonDeltaPhi_%d",index); histo_ZdimuonDeltaPhi[index] = new TH1F(name,"#Delta #phi between the 2 muons from the reconstructed Z",               int(PI_*10),   0.,     PI_);
    sprintf(name,"ZdimuonDeltaEta_%d",index); histo_ZdimuonDeltaEta[index] = new TH1F(name,"#Delta #phi between the 2 muons from the reconstructed Z",                        50,   0.,    5.  );
    sprintf(name,"ZdimuonDeltaR_%d",  index); histo_ZdimuonDeltaR[index]   = new TH1F(name,"#Delta R between the 2 muons from the reconstructed Z",                           50,   0.,    5.  );
    sprintf(name,"ZMass_%d",          index); histo_ZMass[index]           = new TH1F(name,"reconstructed Z mass",                                                           120,   0.,  120.  );
    sprintf(name,"ZMass_pt10Cut_%d",  index); histo_ZMass_pt10Cut[index]   = new TH1F(name,"reconstructed Z mass w/ muons p_{T} > 10.GeV/c",                                 120,   0.,  120.  );
    sprintf(name,"ZMass_pt15Cut_%d",  index); histo_ZMass_pt15Cut[index]   = new TH1F(name,"reconstructed Z mass w/ muons p_{T} > 15.GeV/c",                                 120,   0.,  120.  );
    sprintf(name,"ZMass_pt20Cut_%d",  index); histo_ZMass_pt20Cut[index]   = new TH1F(name,"reconstructed Z mass w/ muons p_{T} > 20.GeV/c",                                 120,   0.,  120.  );
    sprintf(name,"ZPt_%d",            index); histo_ZPt[index]             = new TH1F(name,"reconstructed Z transverse momentum",                                            200,   0.,  200.  );
    sprintf(name,"ZEta_%d",           index); histo_ZEta[index]            = new TH1F(name,"reconstructed Z pseudo-rapidity",                                                100,  -5.,    5.  );
    sprintf(name,"ZPhi_%d",           index); histo_ZPhi[index]            = new TH1F(name,"reconstructed Z azimutal angle",                                         int(PI_*20), -PI_,     PI_);
    sprintf(name,"ZmuonPosPzVSZmuonNegPz_%d",index); histo_ZposMuPzVSZnegMuPz[index]   = new TH2F(name,"reconstructed Z muon VS anti-muon momentum along z-axis",        600,-300.,300.,600,-300.,300.);
    sprintf(name,"ZmuonsPhiVSZmuonsEta_%d",  index); histo_ZmuonsPhiVSZmuonsEta[index] = new TH2F(name,"azimutal angle VS pseudo-rapidity of muons from Z",      int(PI_*20), -PI_, PI_,100,  -5.,  5.);
    sprintf(name,"ZmuonsPtVSZMass_%d",       index); histo_ZmuonsPtVSZMass[index]        = new TH2F(name,"transverse momentum of muons from Z VS reconstructed Z mass",    
												    120,   0.,120.,        200,   0.,    200.);
    sprintf(name,"ZmuonsEtaVSZMass_%d",      index); histo_ZmuonsEtaVSZMass[index]       = new TH2F(name,"pseudo-rapidity of muons from Z VS reconstructed Z mass",        
												    120,   0.,120.,        100,  -5.,      5.);
    sprintf(name,"ZmuonsDxySignifVSZMass_%d",index); histo_ZmuonsDxySignifVSZMass[index] = new TH2F(name,"transverse impact parameter significance of muons from Z VS reconstructed Z mass",
												    120,   0.,120.,        200,-100.,    100.);
    sprintf(name,"ZdimuonDeltaPhiVSZMass_%d",index); histo_ZdimuonDeltaPhiVSZMass[index] = new TH2F(name,"#Delta#phi between muons from Z VS reconstructed Z mass",        
												    120,   0.,120.,int(PI_*10),   0.,     PI_);
    sprintf(name,"ZdimuonDeltaEtaVSZMass_%d",index); histo_ZdimuonDeltaEtaVSZMass[index] = new TH2F(name,"#Delta#eta between muons from Z VS reconstructed Z mass",        
												    120,   0.,120.,         50,   0.,    5.  );
    sprintf(name,"ZdimuonDeltaRVSZMass_%d",  index); histo_ZdimuonDeltaRVSZMass[index]   = new TH2F(name,"#DeltaR between muons from Z VS reconstructed Z mass",    
												    120,   0.,120.,	    50,   0.,    5.  );
    sprintf(name,"ZPtVSZMass_%d",            index); histo_ZPtVSZMass[index]             = new TH2F(name,"transverse momentum VS reconstructed Z mass",        
												    120,   0.,120.,        200,   0.,    200.);
  }

  TH1F * histo_muonsVtxZ  = new TH1F("muonsVtxZ", "muons longitudinal vertex position",200,-20.,20.);
  TH2F * histo_muonsVtxXY = new TH2F("muonsVtxXY","muons transverse vertex position",  300,-0.15,0.15,300,-0.15,0.15);
  TH1F * histo_ZmuonsVtxZ  = new TH1F("ZmuonsVtxZ", "longitudinal vertex position of muons from reconstructed Z",200,-20.,20.);
  TH2F * histo_ZmuonsVtxXY = new TH2F("ZmuonsVtxXY","transverse vertex position of muons from reconstructed Z",  300,-0.15,0.15,300,-0.15,0.15);
  TH1F * histo_ZmuonsDeltaVtxX = new TH1F("ZmuonsDeltaVtxX","difference in vertex x-position between the muons from reconstructed Z",50,0.,0.05);
  TH1F * histo_ZmuonsDeltaVtxY = new TH1F("ZmuonsDeltaVtxY","difference in vertex y-position between the muons from reconstructed Z",50,0.,0.05);
  TH1F * histo_ZmuonsDeltaVtxZ = new TH1F("ZmuonsDeltaVtxZ","difference in vertex z-position between the muons from reconstructed Z",50,0.,0.05);

  // global muons
  math::XYZTLorentzVector null_XYZTLorentzVector = math::XYZTLorentzVector(0.,0.,0.,0.);
  GlobalMuon null_globalMuon = GlobalMuon( 0.,0.,0.,0,
					   null_XYZTLorentzVector,
					   math::XYZPoint(0.,0.,0.),
					   0.,0.,0.,0.,0.,0.,0.,0,0,0.,0,0.,0.);

  unsigned int numberOfEvents = events.GetEntries();

  int muons_number = 0;
  // loop over the events
  for( unsigned int evt_iter = 0;
       evt_iter < numberOfEvents;
       evt_iter++) {

    // need to call SetAddress since TBranch's change for each file read
    gloMuon_B->SetAddress(&gloMuon_vec);
    gloMuon_B->GetEntry(evt_iter);
    events.GetEntry(evt_iter,0);

    // now can access data
    if ( evt_iter/100 == float(evt_iter)/100. ) std::cout << "event = " << evt_iter << std::endl;
    muons_number = gloMuon_vec.size();
    
    //    if ( muons_number != 0 && muons_number != 4 ) {
    if ( muons_number != 0 ) {
      if ( evt_iter/100 == float(evt_iter)/100. ) std::cout << "number of global muons = " << muons_number << std::endl;
      
      GlobalMuon Z_muon1 = null_globalMuon;
      GlobalMuon Z_muon2 = null_globalMuon;
      // loop on all the global muons
      for ( unsigned int gloMuon_iter = 0; gloMuon_iter < gloMuon_vec.size(); gloMuon_iter++ ) { 
	// take the muon
	GlobalMuon * globalMuon = &( gloMuon_vec[gloMuon_iter] );
	
	double                  muon_pt     = globalMuon->pt();
	double                  muon_eta    = globalMuon->eta(); 
	double                  muon_phi    = globalMuon->phi(); 
	int                     muon_charge = globalMuon->charge();
	math::XYZTLorentzVector muon_p4     = globalMuon->p4();
	math::XYZPoint          muon_vertex = globalMuon->vertex();
	double                  muon_em     = globalMuon->caloEm();
	double                  muon_had    = globalMuon->caloHad();
	double                  muon_ho     = globalMuon->caloHo();
	float                   muon_isosumpt   = globalMuon->iso03SumPt();
	float                   muon_isoemet    = globalMuon->iso03emEt();
	float                   muon_isohadet   = globalMuon->iso03hadEt();
	float                   muon_isohoet    = globalMuon->iso03hoEt();
	int                     muon_isonjets   = globalMuon->iso03nJets();
	int                     muon_isontracks = globalMuon->iso03nTracks();
	double                  muon_normchi2   = globalMuon->normChi2();
	int                     muon_hitsnum    = globalMuon->hitsNum();
	double                  muon_dxy        = globalMuon->impactParamXY();
	double                  muon_dxyerror   = globalMuon->errorImpactParamXY();

	double muon_dxysignif = muon_dxy/muon_dxyerror;
	double muon_pz = muon_p4.pz();

	double leptonicZmassUncertainty = 0.3*ZMass_;
	for ( unsigned int gloMuon_iter2 = (gloMuon_vec.size()-1); gloMuon_iter2 > gloMuon_iter; gloMuon_iter2--) {
	  GlobalMuon * globalMuon2 = &( gloMuon_vec[gloMuon_iter2] );
	  int muon2_charge = globalMuon2->charge();
	    
	  // Z->mumu reconstruction through p4 of 2 muons w/ opposite charge
	  if ( muon_charge*muon2_charge < 0 ){
	    double                  muon2_pt  = globalMuon2->pt(); 
	    double                  muon2_eta = globalMuon2->eta(); 
	    double                  muon2_phi = globalMuon2->phi(); 
	    math::XYZTLorentzVector muon2_p4  = globalMuon2->p4();
	    double                  muon2_pz  = muon2_p4.pz();
	    double                  muon2_dxy        = globalMuon->impactParamXY();
	    double                  muon2_dxyerror   = globalMuon->errorImpactParamXY();

	    double muon2_dxysignif = muon2_dxy/muon2_dxyerror;

	    double dimuonDeltaEta = std::fabs(muon_eta-muon2_eta);
	    double dimuonDeltaPhi = PI_ - std::fabs( std::fabs( muon_phi-muon2_phi ) - PI_ );
	    double dimuonDeltaR = TMath::Sqrt(pow(dimuonDeltaEta,2)+pow(dimuonDeltaPhi,2));
	    double dimuonPt   = (muon_p4+muon2_p4).pt();
	    double dimuonEta  = (muon_p4+muon2_p4).eta();
	    double dimuonPhi  = (muon_p4+muon2_p4).phi();
	    double dimuonMass = (muon_p4+muon2_p4).M();
	    if (std::fabs(dimuonMass-ZMass_) < leptonicZmassUncertainty) {
	      leptonicZmassUncertainty = std::fabs(dimuonMass-ZMass_);
	      Z_muon1 = *globalMuon;
	      Z_muon2 = *globalMuon2;
	    }

	    if (muon_charge < 0) histo_posMuPzVSnegMuPz[0]->Fill(muon_pz,muon2_pz);
	    else histo_posMuPzVSnegMuPz[0]->Fill(muon2_pz,muon_pz);

	    histo_muonsPtVSdimuonMass[0]->Fill(dimuonMass,muon_pt);
	    histo_muonsPtVSdimuonMass[0]->Fill(dimuonMass,muon2_pt);
	    histo_muonsPhiVSmuonsEta[0]->Fill(muon_phi,muon_eta);
	    histo_muonsPhiVSmuonsEta[0]->Fill(muon2_phi,muon2_eta);
	    histo_muonsEtaVSdimuonMass[0]->Fill(dimuonMass,muon_eta);
	    histo_muonsEtaVSdimuonMass[0]->Fill(dimuonMass,muon2_eta);
	    histo_muonsDxySignifVSdimuonMass[0]->Fill(dimuonMass,muon_dxysignif);
	    histo_muonsDxySignifVSdimuonMass[0]->Fill(dimuonMass,muon2_dxysignif);
	    histo_dimuonDeltaPhiVSdimuonMass[0]->Fill(dimuonMass,dimuonDeltaPhi);
	    histo_dimuonDeltaEtaVSdimuonMass[0]->Fill(dimuonMass,dimuonDeltaEta);
	    histo_dimuonDeltaRVSdimuonMass[0]->Fill(dimuonMass,dimuonDeltaR);
	    histo_dimuonPtVSdimuonMass[0]->Fill(dimuonMass,dimuonPt);

	    histo_muonsPt[0]->Fill( muon_pt  );
	    histo_muonsPt[0]->Fill( muon2_pt );
	    histo_muonsPt2[0]->Fill( pow(muon_pt,2)  );
	    histo_muonsPt2[0]->Fill( pow(muon2_pt,2) );
	    histo_muonsDxy[0]->Fill(muon_dxy);
	    histo_muonsDxy[0]->Fill(muon2_dxy);
	    histo_muonsDxyError[0]->Fill(muon_dxyerror);
	    histo_muonsDxyError[0]->Fill(muon2_dxyerror);
	    histo_muonsDxySignif[0]->Fill(muon_dxysignif);
	    histo_muonsDxySignif[0]->Fill(muon2_dxysignif);
	    histo_dimuonMass[0]->Fill(     dimuonMass     );
	    histo_dimuonDeltaPhi[0]->Fill( dimuonDeltaPhi );
	    histo_dimuonDeltaEta[0]->Fill( dimuonDeltaEta );
	    histo_dimuonDeltaR[0]->Fill(   dimuonDeltaR   );
	    histo_dimuonPt[0]->Fill(  dimuonPt  );
	    histo_dimuonEta[0]->Fill( dimuonEta );
	    histo_dimuonPhi[0]->Fill( dimuonPhi );
	    int index = -9;
	    switch(muons_number){
	    case(2):
	      index = 1;
	      break;
	    case(3):	
	      index = 2;
	      break;
	    case(4):	
	      index = 3;
	      break;
	    default:
	      index = 3;
	    }
	    if (muon_charge < 0) histo_posMuPzVSnegMuPz[index]->Fill(muon_pz,muon2_pz);
	    else histo_posMuPzVSnegMuPz[index]->Fill(muon2_pz,muon_pz);

	    histo_muonsPtVSdimuonMass[index]->Fill(dimuonMass,muon_pt);
	    histo_muonsPtVSdimuonMass[index]->Fill(dimuonMass,muon2_pt);
	    histo_muonsEtaVSdimuonMass[index]->Fill(dimuonMass,muon_eta);
	    histo_muonsEtaVSdimuonMass[index]->Fill(dimuonMass,muon2_eta);
	    histo_muonsPhiVSmuonsEta[index]->Fill(muon_phi,muon_eta);
	    histo_muonsPhiVSmuonsEta[index]->Fill(muon2_phi,muon2_eta);
	    histo_muonsDxySignifVSdimuonMass[index]->Fill(dimuonMass,muon_dxysignif);
	    histo_muonsDxySignifVSdimuonMass[index]->Fill(dimuonMass,muon2_dxysignif);
	    histo_dimuonDeltaPhiVSdimuonMass[index]->Fill(dimuonMass,dimuonDeltaPhi);
	    histo_dimuonDeltaEtaVSdimuonMass[index]->Fill(dimuonMass,dimuonDeltaEta);
	    histo_dimuonDeltaRVSdimuonMass[index]->Fill(dimuonMass,dimuonDeltaR);
	    histo_dimuonPtVSdimuonMass[index]->Fill(dimuonMass,dimuonPt);
	    histo_muonsPt[index]->Fill( muon_pt  );
	    histo_muonsPt[index]->Fill( muon2_pt );
	    histo_muonsPt2[index]->Fill( pow(muon_pt,2)  );
	    histo_muonsPt2[index]->Fill( pow(muon2_pt,2) );
	    histo_muonsDxy[index]->Fill(       muon_dxy        );
	    histo_muonsDxy[index]->Fill(       muon2_dxy       );
	    histo_muonsDxyError[index]->Fill(  muon_dxyerror   );
	    histo_muonsDxyError[index]->Fill(  muon2_dxyerror  );
	    histo_muonsDxySignif[index]->Fill( muon_dxysignif  );
	    histo_muonsDxySignif[index]->Fill( muon2_dxysignif );
	    histo_dimuonMass[index]->Fill(     dimuonMass     );
	    histo_dimuonDeltaPhi[index]->Fill( dimuonDeltaPhi );
	    histo_dimuonDeltaEta[index]->Fill( dimuonDeltaEta );
	    histo_dimuonDeltaR[index]->Fill(   dimuonDeltaR   );
	    histo_dimuonPt[index]->Fill(  dimuonPt  );
	    histo_dimuonEta[index]->Fill( dimuonEta );
	    histo_dimuonPhi[index]->Fill( dimuonPhi );
	    histo_muonsVtxZ->Fill(muon_vertex.z());
	    histo_muonsVtxXY->Fill(muon_vertex.x(),muon_vertex.y());
	  }
	}
	
	// fill global muon information per muon
	histo_muonCharge->Fill(muon_charge);
	histo_muonPt->Fill(muon_pt);
	histo_muonPz->Fill(muon_pz);
	histo_muonEta->Fill(muon_eta);
	histo_muonPhi->Fill(muon_phi);
	histo_muonDxy->Fill(muon_dxy);
	histo_muonDxyError->Fill(muon_dxyerror);
	histo_muonDxySignif->Fill(muon_dxysignif);
	histo_muonEM->Fill(muon_em);
	histo_muonHAD->Fill(muon_had);
	histo_muonHO->Fill(muon_ho);
	histo_muonIsoSumPt->Fill(muon_isosumpt);
	histo_muonIsoEMet->Fill(muon_isoemet);
	histo_muonIsoHADet->Fill(muon_isohadet);
	histo_muonIsoHOet->Fill(muon_isohoet);
	histo_muonIsoNjets->Fill(muon_isonjets);
	histo_muonIsoNtk->Fill(muon_isontracks);
	histo_muonNormChi2->Fill(muon_normchi2);
	histo_muonHitsNum->Fill(muon_hitsnum);
	
      } // end loop on muons

      if ( muons_number >= 2 &&
	   Z_muon1.p4() != Z_muon2.p4() ) {
	double Zmuon1_eta       = Z_muon1.eta();
	double Zmuon1_phi       = Z_muon1.phi();
	double Zmuon1_pt        = Z_muon1.pt();
	double Zmuon1_pz        = (Z_muon1.p4()).pz();
	double Zmuon1_dxy       = Z_muon1.impactParamXY();
	double Zmuon1_dxyerror  = Z_muon1.errorImpactParamXY();
	double Zmuon1_dxysignif = Zmuon1_dxy/Zmuon1_dxyerror;
	double Zmuon2_eta       = Z_muon2.eta();
	double Zmuon2_phi       = Z_muon2.phi();
	double Zmuon2_pt        = Z_muon2.pt();
	double Zmuon2_pz        = (Z_muon2.p4()).pz();
	double Zmuon2_dxy       = Z_muon2.impactParamXY();
	double Zmuon2_dxyerror  = Z_muon2.errorImpactParamXY();
	double Zmuon2_dxysignif = Zmuon2_dxy/Zmuon2_dxyerror;

	double Zdimuon_deltaEta = std::fabs(Z_muon1.eta()-Z_muon2.eta());
	double Zdimuon_deltaPhi = PI_ - std::fabs( std::fabs( Z_muon1.phi()-Z_muon2.phi() ) - PI_ );
	double Zdimuon_deltaR   = TMath::Sqrt(pow(Zdimuon_deltaEta,2)+pow(Zdimuon_deltaPhi,2));
      
	math::XYZTLorentzVector Z_p4 = Z_muon1.p4()+Z_muon2.p4();
	double Z_mass = Z_p4.M();
	double Z_pt   = Z_p4.pt();
	double Z_phi  = Z_p4.phi();
	double Z_eta  = Z_p4.eta();

	if (Z_muon1.charge() < 0) histo_ZposMuPzVSZnegMuPz[0]->Fill(Zmuon1_pz,Zmuon2_pz);
	else histo_ZposMuPzVSZnegMuPz[0]->Fill(Zmuon2_pz,Zmuon1_pz);
	histo_ZmuonsPtVSZMass[0]->Fill(Z_mass,Zmuon1_pt);
	histo_ZmuonsPtVSZMass[0]->Fill(Z_mass,Zmuon2_pt);
	histo_ZmuonsEtaVSZMass[0]->Fill(Z_mass,Zmuon1_eta);
	histo_ZmuonsEtaVSZMass[0]->Fill(Z_mass,Zmuon2_eta);
	histo_ZmuonsPhiVSZmuonsEta[0]->Fill(Zmuon1_phi,Zmuon1_eta);
	histo_ZmuonsPhiVSZmuonsEta[0]->Fill(Zmuon2_phi,Zmuon2_eta);
	histo_ZmuonsDxySignifVSZMass[0]->Fill(Z_mass,Zmuon1_dxysignif);
	histo_ZmuonsDxySignifVSZMass[0]->Fill(Z_mass,Zmuon2_dxysignif);
	histo_ZdimuonDeltaPhiVSZMass[0]->Fill(Z_mass,Zdimuon_deltaPhi);
	histo_ZdimuonDeltaEtaVSZMass[0]->Fill(Z_mass,Zdimuon_deltaEta);
	histo_ZdimuonDeltaRVSZMass[0]->Fill(Z_mass,Zdimuon_deltaR);
	histo_ZPtVSZMass[0]->Fill(Z_mass,Z_pt);

	histo_ZmuonsPt[0]->Fill(        Zmuon1_pt       );
	histo_ZmuonsPt[0]->Fill(        Zmuon2_pt       );
	histo_ZmuonsPt2[0]->Fill(   pow(Zmuon1_pt,2)    );
	histo_ZmuonsPt2[0]->Fill(   pow(Zmuon2_pt,2)    );
	histo_ZmuonsDxy[0]->Fill(       Zmuon1_dxy      );
	histo_ZmuonsDxy[0]->Fill(       Zmuon2_dxy      );
	histo_ZmuonsDxyError[0]->Fill(  Zmuon1_dxyerror );
	histo_ZmuonsDxyError[0]->Fill(  Zmuon2_dxyerror );
	histo_ZmuonsDxySignif[0]->Fill( Zmuon1_dxysignif);
	histo_ZmuonsDxySignif[0]->Fill( Zmuon2_dxysignif);

	histo_ZdimuonDeltaPhi[0]->Fill( Zdimuon_deltaPhi );
	histo_ZdimuonDeltaEta[0]->Fill( Zdimuon_deltaEta );
	histo_ZdimuonDeltaR[0]->Fill(   Zdimuon_deltaR   );
	histo_ZMass[0]->Fill( Z_mass );
	if (Zmuon1_pt>=10. && Zmuon2_pt>=10.) histo_ZMass_pt10Cut[0]->Fill( Z_mass );
	if (Zmuon1_pt>=15. && Zmuon2_pt>=15.) histo_ZMass_pt15Cut[0]->Fill( Z_mass );
	if (Zmuon1_pt>=20. && Zmuon2_pt>=20.) histo_ZMass_pt20Cut[0]->Fill( Z_mass );
	histo_ZPt[0]->Fill(   Z_pt   );
	histo_ZPhi[0]->Fill(  Z_phi  );
	histo_ZEta[0]->Fill(  Z_eta  );
	int index = -9;
	switch(muons_number){
	case(2):
	  index = 1;
	  histo_ZmuonsVtxZ->Fill((Z_muon1.vertex()).z()); 
	  histo_ZmuonsVtxZ->Fill((Z_muon2.vertex()).z()); 
	  histo_ZmuonsVtxXY->Fill((Z_muon1.vertex()).x(),(Z_muon1.vertex()).y());
	  histo_ZmuonsVtxXY->Fill((Z_muon2.vertex()).x(),(Z_muon2.vertex()).y());
	  histo_ZmuonsDeltaVtxX->Fill(std::fabs((Z_muon1.vertex()).x()-(Z_muon2.vertex()).x()));
	  histo_ZmuonsDeltaVtxY->Fill(std::fabs((Z_muon1.vertex()).y()-(Z_muon2.vertex()).y()));
	  histo_ZmuonsDeltaVtxZ->Fill(std::fabs((Z_muon1.vertex()).z()-(Z_muon2.vertex()).z()));
	  break;
	case(3):	
	  index = 2;
	  break;
	case(4):	
	  index = 3;
	  break;
	default:
	  index = 3;
	}
	if (Z_muon1.charge() < 0) histo_ZposMuPzVSZnegMuPz[index]->Fill(Zmuon1_pz,Zmuon2_pz);
	else histo_ZposMuPzVSZnegMuPz[index]->Fill(Zmuon2_pz,Zmuon1_pz);
	histo_ZmuonsPtVSZMass[index]->Fill(Z_mass,Zmuon1_pt);
	histo_ZmuonsPtVSZMass[index]->Fill(Z_mass,Zmuon2_pt);
	histo_ZmuonsEtaVSZMass[index]->Fill(Z_mass,Zmuon1_eta);
	histo_ZmuonsEtaVSZMass[index]->Fill(Z_mass,Zmuon2_eta);
	histo_ZmuonsPhiVSZmuonsEta[index]->Fill(Zmuon1_phi,Zmuon1_eta);
	histo_ZmuonsPhiVSZmuonsEta[index]->Fill(Zmuon2_phi,Zmuon2_eta);
	histo_ZmuonsDxySignifVSZMass[index]->Fill(Z_mass,Zmuon1_dxysignif);
	histo_ZmuonsDxySignifVSZMass[index]->Fill(Z_mass,Zmuon2_dxysignif);
	histo_ZdimuonDeltaPhiVSZMass[index]->Fill(Z_mass,Zdimuon_deltaPhi);
	histo_ZdimuonDeltaEtaVSZMass[index]->Fill(Z_mass,Zdimuon_deltaEta);
	histo_ZdimuonDeltaRVSZMass[index]->Fill(Z_mass,Zdimuon_deltaR);
	histo_ZPtVSZMass[index]->Fill(Z_mass,Z_pt);

	histo_ZmuonsPt[index]->Fill(        Zmuon1_pt       );
	histo_ZmuonsPt[index]->Fill(        Zmuon2_pt       );
	histo_ZmuonsPt2[index]->Fill(   pow(Zmuon1_pt,2)    );
	histo_ZmuonsPt2[index]->Fill(   pow(Zmuon2_pt,2)    );
	histo_ZmuonsDxy[index]->Fill(       Zmuon1_dxy       );
	histo_ZmuonsDxy[index]->Fill(       Zmuon2_dxy       );
	histo_ZmuonsDxyError[index]->Fill(  Zmuon1_dxyerror  );
	histo_ZmuonsDxyError[index]->Fill(  Zmuon2_dxyerror  );
	histo_ZmuonsDxySignif[index]->Fill( Zmuon1_dxysignif );
	histo_ZmuonsDxySignif[index]->Fill( Zmuon2_dxysignif );
	histo_ZdimuonDeltaPhi[index]->Fill( Zdimuon_deltaPhi );
	histo_ZdimuonDeltaEta[index]->Fill( Zdimuon_deltaEta );
	histo_ZdimuonDeltaR[index]->Fill(   Zdimuon_deltaR   );
	histo_ZMass[index]->Fill( Z_mass );
	if (Zmuon1_pt>=10. && Zmuon2_pt>=10.) histo_ZMass_pt10Cut[index]->Fill( Z_mass );
	if (Zmuon1_pt>=15. && Zmuon2_pt>=15.) histo_ZMass_pt15Cut[index]->Fill( Z_mass );
	if (Zmuon1_pt>=20. && Zmuon2_pt>=20.) histo_ZMass_pt20Cut[index]->Fill( Z_mass );
	histo_ZPt[index]->Fill(   Z_pt   );
	histo_ZPhi[index]->Fill(  Z_phi  );
	histo_ZEta[index]->Fill(  Z_eta  );
      }
    } // end if (muons_number != 0)
    
  
    // fill global muon information per event
    histo_muonsNumber->Fill(muons_number);
    
  } // end loop on events


  // preliminary fit to reconstructed Z mass
  // ---------------------------------------
  TF1 * Zgaus;
  double prob        = -99.;
  double mean_value  = -99.; 
  double mean_error  = -99.;  
  double sigma_value = -99.; 
  double sigma_error = -99.;
  // output latex file
  ofstream latexfile;
  latexfile.open(sampleName+"_leptonicZfit.tex", ofstream::out | ofstream::app);
  char fitResult[50];
  histo_ZMass[1]->Fit("gaus","","",70.,120.);
  Zgaus = histo_ZMass[1]->GetFunction("gaus");
  prob        = Zgaus->GetProb();
  mean_value  = Zgaus->GetParameter(1); 
  mean_error  = Zgaus->GetParError(1);  
  sigma_value = Zgaus->GetParameter(2); 
  sigma_error = Zgaus->GetParError(2);
  std::cout << "reconstructed Z mass: " << mean_value << " +- " << sigma_value << " [prob: " << prob << std::endl;
  latexfile << "\\hline" << endl;
  sprintf(fitResult,"reconstructed Z mass");                  latexfile << fitResult << "&\\\\" << endl;
  sprintf(fitResult,"%.2f \\pm %.2f",mean_value,mean_error);  latexfile << "mean: & $" << fitResult << "$\\\\" << endl;
  sprintf(fitResult,"%.2f \\pm %.2f",sigma_value,sigma_error);latexfile << "sigma: & $" << fitResult << "$\\\\" << endl;
  latexfile << "\\hline" << endl;

  histo_ZMass_pt10Cut[1]->Fit("gaus","","",70.,120.);
  Zgaus = histo_ZMass_pt10Cut[1]->GetFunction("gaus");
  prob        = Zgaus->GetProb();
  mean_value  = Zgaus->GetParameter(1); 
  mean_error  = Zgaus->GetParError(1);  
  sigma_value = Zgaus->GetParameter(2); 
  sigma_error = Zgaus->GetParError(2);
  std::cout << "reconstructed Z mass w/ muon pt>10GeV: " << mean_value << " +- " << sigma_value << " [prob: " << prob << std::endl;
  latexfile << "\\hline" << endl;
  sprintf(fitResult,"reconstructed Z mass w/ muon $p_{T}>10GeV/c$");latexfile << fitResult << "&\\\\" << endl;
  sprintf(fitResult,"%.2f \\pm %.2f",mean_value,mean_error);        latexfile << "mean: & $" << fitResult << "$\\\\" << endl;
  sprintf(fitResult,"%.2f \\pm %.2f",sigma_value,sigma_error);      latexfile << "sigma: & $" << fitResult << "$\\\\" << endl;
  latexfile << "\\hline" << endl;

  histo_ZMass_pt15Cut[1]->Fit("gaus","","",70.,120.);
  Zgaus = histo_ZMass_pt15Cut[1]->GetFunction("gaus");
  prob        = Zgaus->GetProb();
  mean_value  = Zgaus->GetParameter(1); 
  mean_error  = Zgaus->GetParError(1);  
  sigma_value = Zgaus->GetParameter(2); 
  sigma_error = Zgaus->GetParError(2);
  std::cout << "reconstructed Z mass w/ muon pt>15GeV: " << mean_value << " +- " << sigma_value << " [prob: " << prob << std::endl;
  latexfile << "\\hline" << endl;
  sprintf(fitResult,"reconstructed Z mass w/ muon $p_{T}>15GeV/c$");latexfile << fitResult << "&\\\\" << endl;
  sprintf(fitResult,"%.2f \\pm %.2f",mean_value,mean_error);
  latexfile << "mean: & $" << fitResult << "$\\\\" << endl;
  sprintf(fitResult,"%.2f \\pm %.2f",sigma_value,sigma_error);
  latexfile << "sigma: & $" << fitResult << "$\\\\" << endl;
  latexfile << "\\hline" << endl;

  histo_ZMass_pt20Cut[1]->Fit("gaus","","",70.,120.);
  Zgaus = histo_ZMass_pt20Cut[1]->GetFunction("gaus");
  prob        = Zgaus->GetProb();
  mean_value  = Zgaus->GetParameter(1); 
  mean_error  = Zgaus->GetParError(1);  
  sigma_value = Zgaus->GetParameter(2); 
  sigma_error = Zgaus->GetParError(2);
  std::cout << "reconstructed Z mass w/ muon pt>20GeV: " << mean_value << " +- " << sigma_value << " [prob: " << prob << std::endl;
  latexfile << "\\hline" << endl;
  sprintf(fitResult,"reconstructed Z mass w/ muon $p_{T}>20GeV/c$");latexfile << fitResult << "&\\\\" << endl;
  sprintf(fitResult,"%.2f \\pm %.2f",mean_value,mean_error);        latexfile << "mean: & $" << fitResult << "$\\\\" << endl;
  sprintf(fitResult,"%.2f \\pm %.2f",sigma_value,sigma_error);      latexfile << "sigma: & $" << fitResult << "$\\\\" << endl;
  latexfile << "\\hline" << endl;

  latexfile.close();

  
  //  TCanvas * gloMuon_canvas1 = new TCanvas("c1","c1",1000,800);
  //  gloMuon_canvas1->cd();
  //  histo_muonsNumber->GetXaxis()->SetTitle("number of muons");
  //  histo_muonsNumber->GetXaxis()->SetTitleSize(0.04);
  //  histo_muonsNumber->Draw();

  //  TCanvas * gloMuon_canvas2 = new TCanvas("c2","c2",1000,800);
  //  gloMuon_canvas2->Divide(2,2);
  //  gloMuon_canvas2->cd(1);
  //  histo_muonPt->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  //  histo_muonPt->GetXaxis()->SetTitleSize(0.04);
  //  histo_muonPt->Draw();
  //  gloMuon_canvas2->cd(2);
  //  histo_muonEta->GetXaxis()->SetTitle("#eta");
  //  histo_muonEta->GetXaxis()->SetTitleSize(0.04);
  //  histo_muonEta->Draw();
  //  gloMuon_canvas2->cd(3);
  //  histo_muonPhi->GetXaxis()->SetTitle("#phi (rad)");
  //  histo_muonPhi->GetXaxis()->SetTitleSize(0.04);
  //  histo_muonPhi->Draw();
  //  gloMuon_canvas2->cd(4);
  //  histo_muonPz->GetXaxis()->SetTitle("p_{z} (GeV/c)");
  //  histo_muonPz->GetXaxis()->SetTitleSize(0.04);
  //  histo_muonPz->Draw();

  TCanvas * gloMuon_canvas3 = new TCanvas("c3","c3",1000,800);
  gloMuon_canvas3->Divide(3,2);
  gloMuon_canvas3->cd(1);
  histo_muonsDxy[1]->GetXaxis()->SetTitle("d_{xy} (cm)");
  histo_muonsDxy[1]->GetXaxis()->SetTitleSize(0.04);
  histo_muonsDxy[1]->Draw();
  gloMuon_canvas3->cd(2);
  histo_muonsDxyError[1]->GetXaxis()->SetTitle("#sigma_{d_{xy}} (cm)");
  histo_muonsDxyError[1]->GetXaxis()->SetTitleSize(0.04);
  histo_muonsDxyError[1]->Draw();
  gloMuon_canvas3->cd(3);
  histo_muonsDxySignif[1]->GetXaxis()->SetTitle("#frac{d_{xy}}{#sigma_{d_{xy}}}");
  histo_muonsDxySignif[1]->GetXaxis()->SetTitleSize(0.04);
  histo_muonsDxySignif[1]->Draw();
  gloMuon_canvas3->cd(4);
  histo_ZmuonsDxy[1]->GetXaxis()->SetTitle("d_{xy} (cm)");
  histo_ZmuonsDxy[1]->GetXaxis()->SetTitleSize(0.04);
  histo_ZmuonsDxy[1]->Draw();
  gloMuon_canvas3->cd(5);
  histo_ZmuonsDxyError[1]->GetXaxis()->SetTitle("#sigma_{d_{xy}} (cm)");
  histo_ZmuonsDxyError[1]->GetXaxis()->SetTitleSize(0.04);
  histo_ZmuonsDxyError[1]->Draw();
  gloMuon_canvas3->cd(6);
  histo_ZmuonsDxySignif[1]->GetXaxis()->SetTitle("#frac{d_{xy}}{#sigma_{d_{xy}}}");
  histo_ZmuonsDxySignif[1]->GetXaxis()->SetTitleSize(0.04);
  histo_ZmuonsDxySignif[1]->Draw();

  TCanvas * gloMuon_canvas4 = new TCanvas("c4","c4",1000,800);
  gloMuon_canvas4->Divide(2,2);
  gloMuon_canvas4->cd(1);
  histo_muonsDxySignif[2]->SetLineColor(kRed);
  histo_muonsDxySignif[2]->SetFillColor(kRed);
  histo_muonsDxySignif[2]->SetFillStyle(3001);
  histo_ZmuonsDxySignif[1]->GetXaxis()->SetTitle("#frac{d_{xy}}{#sigma_{d_{xy}}}");
  histo_ZmuonsDxySignif[1]->GetXaxis()->SetTitleSize(0.04);
  histo_ZmuonsDxySignif[1]->Draw();
  histo_muonsDxySignif[2]->Draw("same");
  gloMuon_canvas4->cd(2);
  histo_muonsPt[2]->SetFillStyle(3001);
  histo_muonsPt[2]->SetFillColor(kRed);
  histo_muonsPt[2]->SetLineColor(kRed);
  histo_ZmuonsPt[1]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  histo_ZmuonsPt[1]->GetXaxis()->SetTitleSize(0.04);
  histo_ZmuonsPt[1]->Draw();
  histo_muonsPt[2]->Draw("same");
  gloMuon_canvas4->cd(3);
  histo_dimuonMass[2]->SetFillStyle(3001);
  histo_dimuonMass[2]->SetFillColor(kRed);
  histo_dimuonMass[2]->SetLineColor(kRed);
  histo_ZMass[1]->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  histo_ZMass[1]->GetXaxis()->SetTitleSize(0.04);
  histo_ZMass[1]->Draw();
  histo_dimuonMass[2]->Draw("same");
  gloMuon_canvas4->cd(4);
  histo_dimuonDeltaR[2]->SetFillStyle(3001);
  histo_dimuonDeltaR[2]->SetFillColor(kRed);
  histo_dimuonDeltaR[2]->SetLineColor(kRed);
  histo_ZdimuonDeltaR[1]->GetXaxis()->SetTitle("#DeltaR");
  histo_ZdimuonDeltaR[1]->GetXaxis()->SetTitleSize(0.04);
  histo_ZdimuonDeltaR[1]->Draw();
  histo_dimuonDeltaR[2]->Draw("same");

  TCanvas * gloMuon_canvas5 = new TCanvas("c5","c5",1000,800);
  gloMuon_canvas5->Divide(2,2);
  gloMuon_canvas5->cd(1);
  histo_ZmuonsPtVSZMass[1]->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  histo_ZmuonsPtVSZMass[1]->GetXaxis()->SetTitleSize(0.04);
  histo_ZmuonsPtVSZMass[1]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  histo_ZmuonsPtVSZMass[1]->GetYaxis()->SetTitleSize(0.04);
  histo_ZmuonsPtVSZMass[1]->Draw();
  histo_muonsPtVSdimuonMass[2]->SetMarkerColor(kRed);
  histo_muonsPtVSdimuonMass[2]->Draw("same");
  gloMuon_canvas5->cd(2);
  histo_ZmuonsEtaVSZMass[1]->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  histo_ZmuonsEtaVSZMass[1]->GetXaxis()->SetTitleSize(0.04);
  histo_ZmuonsEtaVSZMass[1]->GetYaxis()->SetTitle("#eta_{#mu}");
  histo_ZmuonsEtaVSZMass[1]->GetYaxis()->SetTitleSize(0.04);
  histo_ZmuonsEtaVSZMass[1]->Draw();
  histo_muonsEtaVSdimuonMass[2]->SetMarkerColor(kRed);
  histo_muonsEtaVSdimuonMass[2]->Draw("same");
  gloMuon_canvas5->cd(3);
  histo_ZmuonsPhiVSZmuonsEta[1]->GetXaxis()->SetTitle("#phi (rad)");
  histo_ZmuonsPhiVSZmuonsEta[1]->GetXaxis()->SetTitleSize(0.04);
  histo_ZmuonsPhiVSZmuonsEta[1]->GetYaxis()->SetTitle("#eta_{#mu}");
  histo_ZmuonsPhiVSZmuonsEta[1]->GetYaxis()->SetTitleSize(0.04);
  histo_ZmuonsPhiVSZmuonsEta[1]->Draw();
  histo_muonsPhiVSmuonsEta[2]->SetMarkerColor(kRed);
  histo_muonsPhiVSmuonsEta[2]->Draw("same");
  gloMuon_canvas5->cd(4);

  TCanvas * gloMuon_canvas6 = new TCanvas("c6","c6",1000,800);
  gloMuon_canvas6->Divide(3,2);
  gloMuon_canvas6->cd(1);
  histo_ZmuonsDxySignifVSZMass[1]->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  histo_ZmuonsDxySignifVSZMass[1]->GetXaxis()->SetTitleSize(0.04);
  histo_ZmuonsDxySignifVSZMass[1]->GetYaxis()->SetTitle("#frac{d_{xy}}{#sigma_{d_{xy}}}(#mu^{+}#mu^{-})");
  histo_ZmuonsDxySignifVSZMass[1]->GetYaxis()->SetTitleSize(0.04);
  histo_ZmuonsDxySignifVSZMass[1]->Draw();
  histo_muonsDxySignifVSdimuonMass[2]->SetMarkerColor(kRed);
  histo_muonsDxySignifVSdimuonMass[2]->Draw("same");
  gloMuon_canvas6->cd(2);
  histo_ZdimuonDeltaPhiVSZMass[1]->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  histo_ZdimuonDeltaPhiVSZMass[1]->GetXaxis()->SetTitleSize(0.04);
  histo_ZdimuonDeltaPhiVSZMass[1]->GetYaxis()->SetTitle("#Delta#phi(#mu^{+}#mu^{-})");
  histo_ZdimuonDeltaPhiVSZMass[1]->GetYaxis()->SetTitleSize(0.04);
  histo_ZdimuonDeltaPhiVSZMass[1]->Draw();
  histo_dimuonDeltaPhiVSdimuonMass[2]->SetMarkerColor(kRed);
  histo_dimuonDeltaPhiVSdimuonMass[2]->Draw("same");
  gloMuon_canvas6->cd(3);
  histo_ZdimuonDeltaRVSZMass[1]->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  histo_ZdimuonDeltaRVSZMass[1]->GetXaxis()->SetTitleSize(0.04);
  histo_ZdimuonDeltaRVSZMass[1]->GetYaxis()->SetTitle("#DeltaR(#mu^{+}#mu^{-})");
  histo_ZdimuonDeltaRVSZMass[1]->GetYaxis()->SetTitleSize(0.04);
  histo_ZdimuonDeltaRVSZMass[1]->Draw();
  histo_dimuonDeltaRVSdimuonMass[2]->SetMarkerColor(kRed);
  histo_dimuonDeltaRVSdimuonMass[2]->Draw("same");
  gloMuon_canvas6->cd(4);
  histo_ZPtVSZMass[1]->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  histo_ZPtVSZMass[1]->GetXaxis()->SetTitleSize(0.04);
  histo_ZPtVSZMass[1]->GetYaxis()->SetTitle("reconstructed Z p_{T} (GeV/c)");
  histo_ZPtVSZMass[1]->GetYaxis()->SetTitleSize(0.04);
  histo_ZPtVSZMass[1]->Draw();
  histo_dimuonPtVSdimuonMass[2]->SetMarkerColor(kRed);
  histo_dimuonPtVSdimuonMass[2]->Draw("same");
  gloMuon_canvas6->cd(5);
  histo_ZdimuonDeltaEtaVSZMass[1]->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  histo_ZdimuonDeltaEtaVSZMass[1]->GetXaxis()->SetTitleSize(0.04);
  histo_ZdimuonDeltaEtaVSZMass[1]->GetYaxis()->SetTitle("#Delta#eta(#mu^{+}#mu^{-})");
  histo_ZdimuonDeltaEtaVSZMass[1]->GetYaxis()->SetTitleSize(0.04);
  histo_ZdimuonDeltaEtaVSZMass[1]->Draw();
  histo_dimuonDeltaEtaVSdimuonMass[2]->SetMarkerColor(kRed);
  histo_dimuonDeltaEtaVSdimuonMass[2]->Draw("same");


  TCanvas * gloMuon_canvas7 = new TCanvas("c7","c7",1000,800);
  gloMuon_canvas7->Divide(2,2);
  gloMuon_canvas7->cd(1);
  histo_ZMass[1]->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  histo_ZMass[1]->GetXaxis()->SetTitleSize(0.04);
  histo_ZMass[1]->Draw();
  gloMuon_canvas7->cd(2);
  histo_ZMass_pt10Cut[1]->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  histo_ZMass_pt10Cut[1]->GetXaxis()->SetTitleSize(0.04);
  histo_ZMass_pt10Cut[1]->Draw();
  gloMuon_canvas7->cd(3);
  histo_ZMass_pt15Cut[1]->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  histo_ZMass_pt15Cut[1]->GetXaxis()->SetTitleSize(0.04);
  histo_ZMass_pt15Cut[1]->Draw();
  gloMuon_canvas7->cd(4);
  histo_ZMass_pt20Cut[1]->GetXaxis()->SetTitle("m_{#mu^{+}#mu^{-}} (GeV/c^{2})");
  histo_ZMass_pt20Cut[1]->GetXaxis()->SetTitleSize(0.04);
  histo_ZMass_pt20Cut[1]->Draw();


  TCanvas * gloMuon_canvas8 = new TCanvas("c8","c8",1000,800);
  gloMuon_canvas8->Divide(3,3);
  gloMuon_canvas8->cd(1);
  histo_muonsVtxXY->GetXaxis()->SetTitle("vertex x (cm)");
  histo_muonsVtxXY->GetXaxis()->SetTitleSize(0.04);
  histo_muonsVtxXY->GetYaxis()->SetTitle("vertex y (cm)");
  histo_muonsVtxXY->GetYaxis()->SetTitleSize(0.04);
  histo_muonsVtxXY->Draw();
  gloMuon_canvas8->cd(3);
  histo_muonsVtxZ->GetXaxis()->SetTitle("vertex z (cm)");
  histo_muonsVtxZ->GetXaxis()->SetTitleSize(0.04);
  histo_muonsVtxZ->Draw();
  gloMuon_canvas8->cd(4);
  histo_ZmuonsVtxXY->GetXaxis()->SetTitle("vertex x (cm)");
  histo_ZmuonsVtxXY->GetXaxis()->SetTitleSize(0.04);
  histo_ZmuonsVtxXY->GetYaxis()->SetTitle("vertex y (cm)");
  histo_ZmuonsVtxXY->GetYaxis()->SetTitleSize(0.04);
  histo_ZmuonsVtxXY->Draw();
  gloMuon_canvas8->cd(6);
  histo_ZmuonsVtxZ->GetXaxis()->SetTitle("vertex z (cm)");
  histo_ZmuonsVtxZ->GetXaxis()->SetTitleSize(0.04);
  histo_ZmuonsVtxZ->Draw();
  gloMuon_canvas8->cd(7);
  histo_ZmuonsDeltaVtxX->GetXaxis()->SetTitle("#Deltax (cm)");
  histo_ZmuonsDeltaVtxX->GetXaxis()->SetTitleSize(0.04);
  histo_ZmuonsDeltaVtxX->Draw();
  gloMuon_canvas8->cd(8);
  histo_ZmuonsDeltaVtxY->GetXaxis()->SetTitle("#Deltay (cm)");
  histo_ZmuonsDeltaVtxY->GetXaxis()->SetTitleSize(0.04);
  histo_ZmuonsDeltaVtxY->Draw();
  gloMuon_canvas8->cd(9);
  histo_ZmuonsDeltaVtxZ->GetXaxis()->SetTitle("#Deltaz (cm)");
  histo_ZmuonsDeltaVtxZ->GetXaxis()->SetTitleSize(0.04);
  histo_ZmuonsDeltaVtxZ->Draw();

  TCanvas * gloMuon_canvas9 = new TCanvas("c9","c9",1000,800);
  gloMuon_canvas9->cd();
  histo_muonsVtxXY->GetXaxis()->SetTitle("vertex x (cm)");
  histo_muonsVtxXY->GetXaxis()->SetTitleSize(0.04);
  histo_muonsVtxXY->GetYaxis()->SetTitle("vertex y (cm)");
  histo_muonsVtxXY->GetYaxis()->SetTitleSize(0.04);
  histo_muonsVtxXY->Draw();
  TCanvas * gloMuon_canvas10 = new TCanvas("c10","c10",1000,800);
  gloMuon_canvas10->cd();
  histo_muonsDxy[1]->GetXaxis()->SetTitle("d_{xy} (cm)");
  histo_muonsDxy[1]->GetXaxis()->SetTitleSize(0.04);
  histo_muonsDxy[1]->Draw();


  histo_muonsNumber->Write();
  histo_muonPt->Write();
  histo_muonEta->Write();
  histo_muonPhi->Write();
  histo_muonCharge->Write();
  histo_muonEM->Write();
  histo_muonHAD->Write();
  histo_muonHO->Write();
  histo_muonIsoSumPt->Write();
  histo_muonIsoEMet->Write();
  histo_muonIsoHADet->Write();
  histo_muonIsoHOet->Write();
  histo_muonIsoNjets->Write();
  histo_muonIsoNtk->Write();
  histo_muonNormChi2->Write();
  histo_muonHitsNum->Write();
  histo_muonDxy->Write();
  histo_muonDxyError->Write();
  histo_muonDxySignif->Write();
  histo_muonsVtxZ->Write();
  histo_muonsVtxXY->Write();
  histo_ZmuonsVtxZ->Write(); 
  histo_ZmuonsVtxZ->Write(); 
  histo_ZmuonsVtxXY->Write();
  histo_ZmuonsVtxXY->Write();
  histo_ZmuonsDeltaVtxX->Write();
  histo_ZmuonsDeltaVtxY->Write();
  histo_ZmuonsDeltaVtxZ->Write();

  for ( unsigned int index = 0; index < 4; index++ ) {
    histo_posMuPzVSnegMuPz[index]->Write();
    histo_muonsPtVSdimuonMass[index]->Write();
    histo_muonsEtaVSdimuonMass[index]->Write();
    histo_muonsDxySignifVSdimuonMass[index]->Write();
    histo_dimuonDeltaPhiVSdimuonMass[index]->Write();
    histo_dimuonDeltaRVSdimuonMass[index]->Write();
    histo_dimuonPtVSdimuonMass[index]->Write();
    histo_muonsPhiVSmuonsEta[index]->Write();
    histo_dimuonDeltaPhi[index]->Write();
    histo_dimuonDeltaR[index]->Write();
    histo_muonsDxy[index]->Write();
    histo_muonsDxyError[index]->Write();
    histo_muonsDxySignif[index]->Write();
    histo_muonsPt[index]->Write();
    histo_muonsPt2[index]->Write();
    histo_dimuonMass[index]->Write();
    histo_dimuonPt[index]->Write();
    histo_dimuonEta[index]->Write();
    histo_dimuonPhi[index]->Write();
    histo_ZmuonsDxy[index]->Write();
    histo_ZmuonsDxyError[index]->Write();
    histo_ZmuonsDxySignif[index]->Write();
    histo_ZmuonsPt[index]->Write();
    histo_ZmuonsPt2[index]->Write();
    histo_ZdimuonDeltaPhi[index]->Write();
    histo_ZdimuonDeltaR[index]->Write();
    histo_ZMass[index]->Write();
    histo_ZMass_pt10Cut[index]->Write();
    histo_ZMass_pt15Cut[index]->Write();
    histo_ZMass_pt20Cut[index]->Write();
    histo_ZPt[index] ->Write();
    histo_ZEta[index]->Write();
    histo_ZPhi[index]->Write();
    histo_ZposMuPzVSZnegMuPz[index]->Write();
    histo_ZmuonsPtVSZMass[index]->Write();
    histo_ZmuonsEtaVSZMass[index]->Write();
    histo_ZmuonsDxySignifVSZMass[index]->Write();
    histo_ZdimuonDeltaPhiVSZMass[index]->Write();
    histo_ZdimuonDeltaRVSZMass[index]->Write();
    histo_ZPtVSZMass[index]->Write();
    histo_ZmuonsPhiVSZmuonsEta[index]->Write();
  }

  //  gloMuon_canvas1->Write();
  //  gloMuon_canvas2->Write();
  gloMuon_canvas3->Write();
  gloMuon_canvas4->Write();
  gloMuon_canvas5->Write();
  gloMuon_canvas6->Write();
  gloMuon_canvas7->Write();
  gloMuon_canvas8->Write();

  gStyle->SetOptStat("ne");
  gloMuon_canvas10->Print("plotGlobalMuonsImpactParameter.jpg");
  gloMuon_canvas9->Print("plotGlobalMuonsPrimaryVertex.jpg");
  //  gloMuon_canvas2->Print("plotGlobalMuons2.eps");
  //  gloMuon_canvas3->Print("plotGlobalMuons3.eps");

  outputfile->Write();

  return 0;
};

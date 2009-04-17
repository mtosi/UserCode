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
#include "AnalysisExamples/AnalysisObjects/interface/BaseJet.h"
#include "AnalysisExamples/AnalysisObjects/interface/BaseMEt.h"
#include "AnalysisExamples/AnalysisObjects/interface/OfflineMEt.h"
#include "AnalysisExamples/AnalysisObjects/interface/OfflineJet.h"
#include "AnalysisExamples/AnalysisObjects/interface/MCParticle.h"
#include "AnalysisExamples/AnalysisObjects/interface/GlobalMuon.h"
#include "AnalysisExamples/AnalysisObjects/interface/SimpleElectron.h"
#include "AnalysisExamples/AnalysisObjects/interface/SimpleTau.h"
#include "AnalysisExamples/AnalysisObjects/interface/Summary.h"

#include "AnalysisExamples/AnalysisObjects/interface/PythiaParticleIndex.h"
#include "AnalysisExamples/AnalysisObjects/interface/ParticlesMass.h"
#include "AnalysisExamples/AnalysisObjects/interface/ParticlesCharge.h"

#endif

int plotMuonIsolation ( int sample = 0 ) {

  using namespace std;
  using namespace edm;
  using namespace anaobj;

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
    events.SetBranchAddress("anaobjGlobalMuons_offlineProd_globalMuons_PROD.obj",&gloMuon_vec,&gloMuon_B);
    events.SetBranchAddress("anaobjOfflineJets_offlineProd_offlineJets_PROD.obj",&offJets_vec,&offJets_B);
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
    events.SetBranchAddress("anaobjGlobalMuons_offlineProd_globalMuons_PROD.obj",&gloMuon_vec,&gloMuon_B);
    events.SetBranchAddress("anaobjOfflineJets_offlineProd_offlineJets_PROD.obj",&offJets_vec,&offJets_B);
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
    events.SetBranchAddress("anaobjGlobalMuons_offlineProd_globalMuons_PROD.obj",&gloMuon_vec,&gloMuon_B);
    events.SetBranchAddress("anaobjOfflineJets_offlineProd_offlineJets_PROD.obj",&offJets_vec,&offJets_B);
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
    events.SetBranchAddress("anaobjGlobalMuons_offlineProd_globalMuons_PROD.obj",&gloMuon_vec,&gloMuon_B);
    events.SetBranchAddress("anaobjOfflineJets_offlineProd_offlineJets_PROD.obj",&offJets_vec,&offJets_B);
    break;
  case(10):
    sampleName = "ZZ_0JETS";
    castorSubDirectory = sampleName+"/";
    // input files
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/ZZ_0JETS/ZZ_0JETS_1.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/ZZ_0JETS/ZZ_0JETS_2.root");
    events.SetBranchAddress("anaobjGlobalMuons_offlineProd_globalMuons_PROD1.obj",&gloMuon_vec,&gloMuon_B);
    events.SetBranchAddress("anaobjOfflineJets_offlineProd_offlineJets_PROD1.obj",&offJets_vec,&offJets_B);
    break;
  case(11):
    sampleName = "ZZ_1JETS";
    castorSubDirectory = sampleName+"/";
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/ZZ_1JETS/ZZ_1JETS_1.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/ZZ_1JETS/ZZ_1JETS_2.root");
    events.SetBranchAddress("anaobjGlobalMuons_offlineProd_globalMuons_PROD1.obj",&gloMuon_vec,&gloMuon_B);
    events.SetBranchAddress("anaobjOfflineJets_offlineProd_offlineJets_PROD1.obj",&offJets_vec,&offJets_B);
    break;
  case(12):
    sampleName = "ZZ_2JETS";
    castorSubDirectory = sampleName+"/";
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/ZZ_2JETS/ZZ_2JETS_1.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/ZZ_2JETS/ZZ_2JETS_2.root");
    events.SetBranchAddress("anaobjGlobalMuons_offlineProd_globalMuons_PROD1.obj",&gloMuon_vec,&gloMuon_B);
    events.SetBranchAddress("anaobjOfflineJets_offlineProd_offlineJets_PROD1.obj",&offJets_vec,&offJets_B);
    break;
  case(13):
    sampleName = "ZZ_3JETS";
    castorSubDirectory = sampleName+"/";
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/ZZ_3JETS/ZZ_3JETS_1.root");
    events.Add("castor:/castor/cern.ch/user/t/tosi/FastSim/NOPU/ZZ_3JETS/ZZ_3JETS_2.root");
    events.SetBranchAddress("anaobjGlobalMuons_offlineProd_globalMuons_PROD1.obj",&gloMuon_vec,&gloMuon_B);
    events.SetBranchAddress("anaobjOfflineJets_offlineProd_offlineJets_PROD1.obj",&offJets_vec,&offJets_B);
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
    events.SetBranchAddress("anaobjGlobalMuons_offlineProd_globalMuons_PROD.obj",&gloMuon_vec,&gloMuon_B);
    events.SetBranchAddress("anaobjOfflineJets_offlineProd_offlineJets_PROD.obj",&offJets_vec,&offJets_B);
  }


  // output file
  TFile *outputfile;
  outputfile = dynamic_cast<TFile*>(gROOT->FindObject(sampleName+"plotMuonIsolation.root")); 
  if (outputfile) outputfile->Close();
  outputfile = new TFile(sampleName+"plotMuonIsolation.root","RECREATE","jets histograms");

  TH1F * histo_muonJetDeltaR    = new TH1F("muonJetDeltaR",   "distance in #eta-#phi space between muon and jet",        100,0.,10.);
  TH1F * histo_muonJetDeltaRmin = new TH1F("muonJetDeltaRmin","minimum distance in #eta-#phi space between muon and jet",100,0., 3.);

  TH2F * histo_muonJetDeltaRminVSmuonsNumber = new TH2F("muonJetDeltaRminVSmuonsNumber","minimum distance between muon and jet VS number of global muons",           100,0.,3.,         20,   0., 10.);
  TH2F * histo_muonJetDeltaRminVSmuonPt      = new TH2F("muonJetDeltaRminVSmuonPt",     "minimum distance between muon and jet VS global muon transverse momentum",  100,0.,3.,        150,   0.,200.);
  TH2F * histo_muonJetDeltaRminVSmuonPz      = new TH2F("muonJetDeltaRminVSmuonPz",     "minimum distance between muon and jet VS global muon momentum along z-axis",100,0.,3.,        600,-300.,300.);
  TH2F * histo_muonJetDeltaRminVSmuonEta     = new TH2F("muonJetDeltaRminVSmuonEta",    "minimum distance between muon and jet VS global muon pseudo-rapidity",      100,0.,3.,        100,  -5.,  5.);
  TH2F * histo_muonJetDeltaRminVSmuonPhi     = new TH2F("muonJetDeltaRminVSmuonPhi",    "minimum distance between muon and jet VS global muon azimutal angle",       100,0.,3.,int(PI_*20), -PI_, PI_);
  TH2F * histo_muonJetDeltaRminVSdimuonMass  = new TH2F("muonJetDeltaRminVSdimuonMass", "minimum distance between muon and jet VS invariant mass of 2 muons system", 100,0.,3.,        120,   0.,120.);
	 					     
  TH2F * histo_muonJetDeltaRminVSjetsNumber  = new TH2F("muonJetDeltaRminVSjetsNumber", "minimum distance between muon and jet VS number of jets",                                              100,0.,3.,  50,   0.,  50.); 
  TH2F * histo_muonJetDeltaRminVSjetUncorrEt = new TH2F("muonJetDeltaRminVSjetUncorrEt","minimum distance between muon and jet VS jet uncorrected transverse energy",                           100,0.,3., 300,   0., 300.);
  TH2F * histo_muonJetDeltaRminVSjetEt       = new TH2F("muonJetDeltaRminVSjetEt",      "minimum distance between muon and jet VS jet corrected transverse energy",                             100,0.,3., 300,   0., 300.);
  TH2F * histo_muonJetDeltaRminVSjetEta      = new TH2F("muonJetDeltaRminVSjetEta",     "minimum distance between muon and jet VS jet pseudo-rapidity",                                         100,0.,3., 200,  -5.,   5.);
  TH2F * histo_muonJetDeltaRminVSjetPhi      = new TH2F("muonJetDeltaRminVSjetPhi",     "minimum distance between muon and jet VS jet azimutal angle",                                          100,0.,3., 120, -PI_,  PI_);
  TH2F * histo_muonJetDeltaRminVSjetMass     = new TH2F("muonJetDeltaRminVSjetMass",    "minimum distance between muon and jet VS jet mass",                                                    100,0.,3., 200,   0., 100.);
  TH2F * histo_muonJetDeltaRminVSdijetMass   = new TH2F("muonJetDeltaRminVSdijetMass",  "minimum distance between muon and jet VS invariant mass of 2 jets w/ E_{T}#ge 40.GeV && |#eta|#le 2.5",100,0.,3.,1000,   0.,1000.);


  math::XYZTLorentzVector null_XYZTLorentzVector = math::XYZTLorentzVector(0.,0.,0.,0.);
  GlobalMuon null_globalMuon = GlobalMuon( 0.,0.,0.,0,
					   null_XYZTLorentzVector,
					   math::XYZPoint(0.,0.,0.),
					   0.,0.,0.,0.,0.,0.,0.,0,0,0.,0,0.,0.);
  OfflineJet null_offlineJet = OfflineJet( 0.,0.,0.,0.,0.,
					   null_XYZTLorentzVector,
					   math::XYZPoint(0.,0.,0.),
					   0.,0.,0.,0, 0., 0.,0, 0., 0.,0, 0., 0. );
  unsigned int numberOfEvents = events.GetEntries();

  int jets_number  = -99;
  int muons_number = -99;
  // loop over the events
  for( unsigned int evt_iter = 0;
       evt_iter < numberOfEvents;
       ++evt_iter) {

    // need to call SetAddress since TBranch's change for each file read
    offJets_B->SetAddress(&offJets_vec);
    offJets_B->GetEntry(evt_iter);
    events.GetEntry(evt_iter,0);
    gloMuon_B->SetAddress(&gloMuon_vec);
    gloMuon_B->GetEntry(evt_iter);
    events.GetEntry(evt_iter,0);

    // now can access data
    if ( evt_iter/100 == float(evt_iter)/100. ) std::cout << "event = " << evt_iter << std::endl;
    jets_number = offJets_vec.size();
    muons_number = gloMuon_vec.size();

    if ( jets_number != 0 ) {
      double muon_jetDeltaRmin = 3.;
      if ( muons_number != 0 ) {
	if ( evt_iter/100 == float(evt_iter)/100. ) std::cout << "number of global muons = " << muons_number << std::endl;
	if ( evt_iter/100 == float(evt_iter)/100. ) std::cout << "number of jets = " << jets_number << std::endl;
      
	GlobalMuon Z_muon1 = null_globalMuon;
	GlobalMuon Z_muon2 = null_globalMuon;
	// loop on all the global muons
	for ( unsigned int gloMuon_iter = 0; gloMuon_iter < gloMuon_vec.size(); gloMuon_iter++ ) { 
	  muon_jetDeltaRmin = 3.;
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

	  double dimuonMass = -99.;
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
	      dimuonMass = (muon_p4+muon2_p4).M();
	      if (std::fabs(dimuonMass-ZMass_) < leptonicZmassUncertainty) {
		leptonicZmassUncertainty = std::fabs(dimuonMass-ZMass_);
		Z_muon1 = *globalMuon;
		Z_muon2 = *globalMuon2;
	      }

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
	      }
	    }
	  }

	  OfflineJet closestJet = null_offlineJet;
	  // loop on all the jets
	  for ( unsigned int offJets_iter = 0; offJets_iter < offJets_vec.size(); ++offJets_iter ) { 
	    // take the jet
	    OfflineJet * offlineJet = &( offJets_vec[offJets_iter] );
	    
	    math::XYZTLorentzVector jet_p4       = offlineJet->p4();
	    double                  jet_uncorrEt = offlineJet->uncorrEt();
	    double                  jet_et       = offlineJet->et(); 
	    double                  jet_eta      = offlineJet->eta();
	    double                  jet_phi      = offlineJet->phi();
	    double                  jet_mass     = offlineJet->jetMass();

	    double muon_jetDeltaPhi =  PI_ - std::fabs( std::fabs( muon_phi-jet_phi ) - PI_ );
	    double muon_jetDeltaEta = std::fabs(muon_eta-jet_eta);
	    double muon_jetDeltaR   = TMath::Sqrt(pow(muon_jetDeltaPhi,2)+pow(muon_jetDeltaEta,2));
	    if ( muon_jetDeltaR <= muon_jetDeltaRmin ) {
	      muon_jetDeltaRmin = muon_jetDeltaR;
	      closestJet = *offlineJet;
	    }
	    histo_muonJetDeltaR->Fill(muon_jetDeltaR);

	  } // end loop on jets

	  if ( closestJet.et() >= 40. && std::fabs(closestJet.eta()) <= 2.5  ) {
	    for ( unsigned int offJets_iter = 0; offJets_iter < offJets_vec.size(); offJets_iter++ ) { 
	      if ( offJets_vec[offJets_iter].et() >= 40. &&
		   std::fabs(offJets_vec[offJets_iter].eta()) <= 2.5 ) {
		math::XYZTLorentzVector jet2_p4 = offJets_vec[offJets_iter].p4();
		if (closestJet.p4() != jet2_p4)	{
		  double dijetMass = (closestJet.p4()+jet2_p4).M();
		  histo_muonJetDeltaRminVSdijetMass->Fill(muon_jetDeltaRmin,dijetMass);
		}
	      } 
	    }
	  }

	  histo_muonJetDeltaRmin->Fill(muon_jetDeltaRmin);
	  
	  histo_muonJetDeltaRminVSmuonPt->Fill(muon_jetDeltaRmin,muon_pt);
	  histo_muonJetDeltaRminVSmuonPz->Fill(muon_jetDeltaRmin,muon_pz);
	  histo_muonJetDeltaRminVSmuonEta->Fill(muon_jetDeltaRmin,muon_eta);
	  histo_muonJetDeltaRminVSmuonPhi->Fill(muon_jetDeltaRmin,muon_phi);
	  if (dimuonMass != -99.) histo_muonJetDeltaRminVSdimuonMass->Fill(muon_jetDeltaRmin,dimuonMass);

	  histo_muonJetDeltaRminVSjetsNumber->Fill(muon_jetDeltaRmin,jets_number);
	  histo_muonJetDeltaRminVSjetUncorrEt->Fill(muon_jetDeltaRmin,closestJet.uncorrEt());
	  histo_muonJetDeltaRminVSjetEt->Fill(muon_jetDeltaRmin,  closestJet.et());
	  histo_muonJetDeltaRminVSjetEta->Fill(muon_jetDeltaRmin, closestJet.eta());
	  histo_muonJetDeltaRminVSjetPhi->Fill(muon_jetDeltaRmin, closestJet.phi());
	  histo_muonJetDeltaRminVSjetMass->Fill(muon_jetDeltaRmin,closestJet.jetMass());

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

	}

      } // end if (muons_number != 0)
      histo_muonJetDeltaRminVSmuonsNumber->Fill(muon_jetDeltaRmin,muons_number); 
    } // end if (jets_number != 0)
    // fill jets information per event
    
  } // end loop on events

  TCanvas * muIso_canvas1 = new TCanvas("c1", "c1", 1000, 800);
  muIso_canvas1->Divide(2,1);
  muIso_canvas1->cd(1);
  histo_muonJetDeltaR->GetXaxis()->SetTitle("#DeltaR(#mu,jet)");
  histo_muonJetDeltaR->GetXaxis()->SetTitleSize(0.04);
  histo_muonJetDeltaR->Draw();
  muIso_canvas1->cd(2);
  histo_muonJetDeltaRmin->GetXaxis()->SetTitle("#DeltaR(#mu,jet)_{min}");
  histo_muonJetDeltaRmin->GetXaxis()->SetTitleSize(0.04);
  histo_muonJetDeltaRmin->Draw();

  TCanvas * muIso_canvas2 = new TCanvas("c2", "c2", 1000, 800);
  muIso_canvas2->Divide(3,2);
  muIso_canvas2->cd(1);
  histo_muonJetDeltaRminVSmuonsNumber->GetXaxis()->SetTitle("#DeltaR(#mu,jet)_{min}");
  histo_muonJetDeltaRminVSmuonsNumber->GetXaxis()->SetTitleSize(0.04);
  histo_muonJetDeltaRminVSmuonsNumber->GetYaxis()->SetTitle("muons number");
  histo_muonJetDeltaRminVSmuonsNumber->GetYaxis()->SetTitleSize(0.04);
  histo_muonJetDeltaRminVSmuonsNumber->Draw();
  muIso_canvas2->cd(2);
  histo_muonJetDeltaRminVSmuonPt->GetXaxis()->SetTitle("#DeltaR(#mu,jet)_{min}");
  histo_muonJetDeltaRminVSmuonPt->GetXaxis()->SetTitleSize(0.04);
  histo_muonJetDeltaRminVSmuonPt->GetYaxis()->SetTitle("muons p_{T} (GeV/c)");
  histo_muonJetDeltaRminVSmuonPt->GetYaxis()->SetTitleSize(0.04);
  histo_muonJetDeltaRminVSmuonPt->Draw();
  muIso_canvas2->cd(3);
  histo_muonJetDeltaRminVSmuonPz->GetXaxis()->SetTitle("#DeltaR(#mu,jet)_{min}");
  histo_muonJetDeltaRminVSmuonPz->GetXaxis()->SetTitleSize(0.04);
  histo_muonJetDeltaRminVSmuonPz->GetYaxis()->SetTitle("muon p_{z} (GeV/c)");
  histo_muonJetDeltaRminVSmuonPz->GetYaxis()->SetTitleSize(0.04);
  histo_muonJetDeltaRminVSmuonPz->Draw();
  muIso_canvas2->cd(4);
  histo_muonJetDeltaRminVSmuonEta->GetXaxis()->SetTitle("#DeltaR(#mu,jet)_{min}");
  histo_muonJetDeltaRminVSmuonEta->GetXaxis()->SetTitleSize(0.04);
  histo_muonJetDeltaRminVSmuonEta->GetYaxis()->SetTitle("muon #eta");
  histo_muonJetDeltaRminVSmuonEta->GetYaxis()->SetTitleSize(0.04);
  histo_muonJetDeltaRminVSmuonEta->Draw();
  muIso_canvas2->cd(5);
  histo_muonJetDeltaRminVSmuonPhi->GetXaxis()->SetTitle("#DeltaR(#mu,jet)_{min}");
  histo_muonJetDeltaRminVSmuonPhi->GetXaxis()->SetTitleSize(0.04);
  histo_muonJetDeltaRminVSmuonPhi->GetYaxis()->SetTitle("muon #phi (rad)");
  histo_muonJetDeltaRminVSmuonPhi->GetYaxis()->SetTitleSize(0.04);
  histo_muonJetDeltaRminVSmuonPhi->Draw(); 					     
  muIso_canvas2->cd(6);
  histo_muonJetDeltaRminVSdimuonMass->GetXaxis()->SetTitle("#DeltaR(#mu,jet)_{min}");
  histo_muonJetDeltaRminVSdimuonMass->GetXaxis()->SetTitleSize(0.04);
  histo_muonJetDeltaRminVSdimuonMass->GetYaxis()->SetTitle("m_{#mu^{+}#mu{-}} (GeV/c^{2})");
  histo_muonJetDeltaRminVSdimuonMass->GetYaxis()->SetTitleSize(0.04);
  histo_muonJetDeltaRminVSdimuonMass->Draw(); 					     

  TCanvas * muIso_canvas3 = new TCanvas("c3", "c3", 1000, 800);
  muIso_canvas3->Divide(4,2);
  muIso_canvas3->cd(1);
  histo_muonJetDeltaRminVSjetsNumber->GetXaxis()->SetTitle("#DeltaR(#mu,jet)_{min}"); 
  histo_muonJetDeltaRminVSjetsNumber->GetXaxis()->SetTitleSize(0.04); 
  histo_muonJetDeltaRminVSjetsNumber->GetYaxis()->SetTitle("jets number"); 
  histo_muonJetDeltaRminVSjetsNumber->GetYaxis()->SetTitleSize(0.04); 
  histo_muonJetDeltaRminVSjetsNumber->Draw(); 
  muIso_canvas3->cd(2);
  histo_muonJetDeltaRminVSjetUncorrEt->GetXaxis()->SetTitle("#DeltaR(#mu,jet)_{min}"); 
  histo_muonJetDeltaRminVSjetUncorrEt->GetXaxis()->SetTitleSize(0.04); 
  histo_muonJetDeltaRminVSjetUncorrEt->GetYaxis()->SetTitle("jet E_{T}^{uncorr} (GeV)"); 
  histo_muonJetDeltaRminVSjetUncorrEt->GetYaxis()->SetTitleSize(0.04); 
  histo_muonJetDeltaRminVSjetUncorrEt->Draw();
  muIso_canvas3->cd(3);
  histo_muonJetDeltaRminVSjetEt->GetXaxis()->SetTitle("#DeltaR(#mu,jet)_{min}"); 
  histo_muonJetDeltaRminVSjetEt->GetXaxis()->SetTitleSize(0.04); 
  histo_muonJetDeltaRminVSjetEt->GetYaxis()->SetTitle("jet E_{T} (GeV)"); 
  histo_muonJetDeltaRminVSjetEt->GetYaxis()->SetTitleSize(0.04); 
  histo_muonJetDeltaRminVSjetEt->Draw();
  muIso_canvas3->cd(4);
  histo_muonJetDeltaRminVSjetEta->GetXaxis()->SetTitle("#DeltaR(#mu,jet)_{min}"); 
  histo_muonJetDeltaRminVSjetEta->GetXaxis()->SetTitleSize(0.04); 
  histo_muonJetDeltaRminVSjetEta->GetYaxis()->SetTitle("jet #eta"); 
  histo_muonJetDeltaRminVSjetEta->GetYaxis()->SetTitleSize(0.04); 
  histo_muonJetDeltaRminVSjetEta->Draw();
  muIso_canvas3->cd(5);
  histo_muonJetDeltaRminVSjetPhi->GetXaxis()->SetTitle("#DeltaR(#mu,jet)_{min}"); 
  histo_muonJetDeltaRminVSjetPhi->GetXaxis()->SetTitleSize(0.04); 
  histo_muonJetDeltaRminVSjetPhi->GetYaxis()->SetTitle("jet #phi (rad)"); 
  histo_muonJetDeltaRminVSjetPhi->GetYaxis()->SetTitleSize(0.04); 
  histo_muonJetDeltaRminVSjetPhi->Draw();
  muIso_canvas3->cd(6);
  histo_muonJetDeltaRminVSjetMass->GetXaxis()->SetTitle("#DeltaR(#mu,jet)_{min}"); 
  histo_muonJetDeltaRminVSjetMass->GetXaxis()->SetTitleSize(0.04); 
  histo_muonJetDeltaRminVSjetMass->GetYaxis()->SetTitle("jet mass (GeV/c^{2}"); 
  histo_muonJetDeltaRminVSjetMass->GetYaxis()->SetTitleSize(0.04); 
  histo_muonJetDeltaRminVSjetMass->Draw();
  muIso_canvas3->cd(7);
  histo_muonJetDeltaRminVSdijetMass->GetXaxis()->SetTitle("#DeltaR(#mu,jet)_{min}"); 
  histo_muonJetDeltaRminVSdijetMass->GetXaxis()->SetTitleSize(0.04); 
  histo_muonJetDeltaRminVSdijetMass->GetYaxis()->SetTitle("m_{jj} (GeV/c^{2}"); 
  histo_muonJetDeltaRminVSdijetMass->GetYaxis()->SetTitleSize(0.04); 
  histo_muonJetDeltaRminVSdijetMass->Draw();


  muIso_canvas1->Write();
  muIso_canvas2->Write();
  muIso_canvas3->Write();

  histo_muonJetDeltaR->Write();
  histo_muonJetDeltaRmin->Write();

  histo_muonJetDeltaRminVSmuonsNumber->Write();
  histo_muonJetDeltaRminVSmuonPt->Write();
  histo_muonJetDeltaRminVSmuonPz->Write();
  histo_muonJetDeltaRminVSmuonEta->Write();
  histo_muonJetDeltaRminVSmuonPhi->Write();
  histo_muonJetDeltaRminVSdimuonMass->Write();

  histo_muonJetDeltaRminVSjetsNumber->Write(); 
  histo_muonJetDeltaRminVSjetUncorrEt->Write();
  histo_muonJetDeltaRminVSjetEt->Write();
  histo_muonJetDeltaRminVSjetEta->Write();
  histo_muonJetDeltaRminVSjetPhi->Write();
  histo_muonJetDeltaRminVSjetMass->Write();
  histo_muonJetDeltaRminVSdijetMass->Write();

  outputfile->Write();

  return 0;
};

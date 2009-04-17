/*
  created by mia tosi on the Sep 2st, 2008
  - macro skeleton for the jets analysis:
  - input file w/ one of its branches [offlineJets]
  - loop over the events
  - loop over the jets collection per event
  - get jets information per jets
  - fill some histograms on jets 
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

#endif

int plotMCParticles () {

  using namespace std;
  using namespace edm;
  using namespace anaobj;

  // setting overflow, underflow and other options
  gStyle->SetOptStat("nemrou");


  const double PI_= TMath::Pi();

  // output file
  TFile *outputfile;
  outputfile = dynamic_cast<TFile*>(gROOT->FindObject("plotMCParticles.root")); 
  if (outputfile) outputfile->Close();
  outputfile = new TFile("plotMCParticles.root","RECREATE","MC particles histograms");

  TChain events("Events");

  // input files
  //  events.Add("rfio:/castor/cern.ch/user/t/tosi/FastSim/NOPU/HIGGS_ZZ/VBFHIGGS200_ZZ_1.root");
  events.Add("VBFHIGGS200_ZZ_1.root");
//  events.Add("VBFHIGGS200_ZZ_2.root");
//  events.Add("VBFHIGGS200_ZZ_3.root");
//  events.Add("VBFHIGGS200_ZZ_4.root");
//  events.Add("VBFHIGGS200_ZZ_5.root");
//  events.Add("VBFHIGGS200_ZZ_6.root");
//  events.Add("VBFHIGGS200_ZZ_7.root");
//  events.Add("VBFHIGGS200_ZZ_8.root");
//  events.Add("VBFHIGGS200_ZZ_9.root");
//  events.Add("VBFHIGGS200_ZZ_10.root");
  
  // set the buffers for the branches
  //MCParticles
  std::vector<MCParticle> mcparticle_vec;
  TBranch* mcparticle_B;
  events.SetBranchAddress("anaobjMCParticles_offlineProd_MCParticles_PROD.obj",&mcparticle_vec,&mcparticle_B);

  // MC particles  histogram
  TH1F * histo_mcPid  = new TH1F ("mcPid",  "MC particle id",                                42, -16.,   26.);
  TH1F * histo_mcMPid = new TH1F ("mcMPid", "mother's MC particle id",                       42, -16.,   26.);
  TH1F * histo_mcMass = new TH1F ("mcMass", "MC particle mass",                             250,   0.,  250.);
  TH1F * histo_mcPt   = new TH1F ("mcPt",   "MC particle transerve momentum",               400,   0.,  700.);
  TH1F * histo_mcPx   = new TH1F ("mcPx",   "MC particle transerve momentum along x-axis",  500, -500., 500.);
  TH1F * histo_mcPy   = new TH1F ("mcPy",   "MC particle transerve momentum along y-axis",  500, -500., 500.);
  TH1F * histo_mcPz   = new TH1F ("mcPz",   "MC particle transerve momentum along z-axis", 2000,-2000.,2000.);
  TH1F * histo_mcP    = new TH1F ("mcP",    "MC particle momentum",                         500 ,   0.,3000.);
  TH1F * histo_mcEta  = new TH1F ("mcEta",  "MC particle pseudo-rapidity",                  200,  -10.,  10.);
  TH1F * histo_mcPhi  = new TH1F ("mcPhi",  "MC particle azimutal angle",                   120,  -PI_,  PI_);

  TH1F * histo_mceta2 = new TH1F ("eta_quarks_daZ", "eta_quarks_daZ", 200, -7., 7.);
  TH1F * histo_mcphi2 = new TH1F ("phi_quarks_daZ", "phi_quarks_daZ", 100 , -PI_, PI_);
  TH1F * histo_mcpt2 = new TH1F ("pt_quarks_daZ", "pt quarks daZ", 200 , 0., 450.);
  TH1F * histo_mcInvv = new TH1F ("massa_invariante_2quarks_daZ", "massa invariante 2quarks daZ", 200 , 0., 400.);

  TH1F * histo_Qpt = new TH1F ("pt_quark_antiquark","pt_quark_antiquark", 200, 0., 400.); 
  TH1F * histo_Qphi = new TH1F ("phi_quark_antiquark","phi_quark_antiquark", 100, -PI_, PI_); 
  TH1F * histo_Qeta = new TH1F ("eta_quark_antiquark","eta_quark_antiquark", 200, -9., 12.); 
  TH1F * histo_Qmass = new TH1F ("mass_quark_antiquark","mass_quark_antiquark", 200,-10., 250.); 
  TH1F * histo_Qp = new TH1F ("p_quark_antiquark","p_quark_antiquark", 300, 0., 3000.); 
  TH1F * histo_Qpx = new TH1F ("px_quark_antiquark","px_quark_antiquark", 300,-500., 500.); 
  TH1F * histo_Qpy = new TH1F ("py_quark_antiquark","py_quark_antiquark", 300,-700., 700.); 
  TH1F * histo_Qpz = new TH1F ("pz_quark_antiquark","pz_quark_antiquark", 300, 3000., 3000.); 

  int mcparticle_number = 0;
  int jets_number = 0;
  int muon_number = 0;

  unsigned int nCuts = 4;
  unsigned int lightQcounter[nCuts];
  unsigned int cQcounter[nCuts];   
  unsigned int bQcounter[nCuts];   
  unsigned int hadronicZcounter[nCuts];
  unsigned int Hcounter[nCuts];
  for (unsigned int index = 0; index < nCuts; ++index ) {
    lightQcounter[index]    = 0;
    cQcounter[index]        = 0;   
    bQcounter[index]        = 0;
    hadronicZcounter[index] = 0;
    Hcounter[index]         = 0;
  }
  unsigned int leptonicZcounter = 0;

  // loop over the events
  unsigned int numberOfEvents = events.GetEntries();
  for( unsigned int evt_iter = 0;
       evt_iter < numberOfEvents;
       ++evt_iter) {

    // need to call SetAddress since TBranch's change for each file read
    mcparticle_B->SetAddress(&mcparticle_vec);
    mcparticle_B->GetEntry(evt_iter);

    events.GetEntry(evt_iter,0);

    // now can access data
    if ( evt_iter/100 == float(evt_iter)/100. ) std::cout << "event = " << evt_iter << std::endl;
    mcparticle_number = mcparticle_vec.size();

    //use vector to store quarks information
    vector<MCParticle> info;

    if ( evt_iter/100 == float(evt_iter)/100. ) std::cout << "number of MCParticle = " << mcparticle_number << std::endl;
    if ( mcparticle_number != 0 ) {
      if ( evt_iter/100 == float(evt_iter)/100. ) std::cout << "number of MCParticle = " << mcparticle_number << std::endl;
      
      math::XYZTLorentzVector p4Zleptonic = math::XYZTLorentzVector(0.,0.,0.,0.);
      math::XYZTLorentzVector p4Zhadronic[nCuts];
      math::XYZTLorentzVector p4Higgs[nCuts];
      for (unsigned int index = 0; index < nCuts; ++index ) { 
	p4Zhadronic[index] = math::XYZTLorentzVector(0.,0.,0.,0.);
	p4Higgs[index] = math::XYZTLorentzVector(0.,0.,0.,0.);
      }
      // loop on all the MCParticle
      for ( unsigned int mcparticle_iter1 = 0; mcparticle_iter1 < mcparticle_vec.size(); ++mcparticle_iter1 ) { 
	// take the MCParticle
	MCParticle * MCParticle1 = &( mcparticle_vec[mcparticle_iter1] );
	
	double mc_eta = MCParticle1->eta();
	int mc_mPid=MCParticle1->mPid();
	double mc_mass=MCParticle1->mass();
	double mc_p=MCParticle1->p();
	double mc_phi=MCParticle1->phi();
	int mc_pid=MCParticle1->pid();
	double mc_pt=MCParticle1->pt();
	double mc_px=MCParticle1->px();
	double mc_py=MCParticle1->py();
	double mc_pz=MCParticle1->pz();

	if ( mc_mPid == pythiaZ_ &&  std::fabs(mc_pid) < pythiat_ ) {
	  histo_MCeta2->Fill(mc_eta); 
	  histo_MCphi2->Fill(mc_phi);
	  histo_MCpt2->Fill(mc_pt);

	  info.push_back(*MCParticle1);

	}
	histo_MCeta->Fill(mc_eta);
	histo_MCmPid->Fill(mc_mPid);
	histo_MCmass->Fill(mc_mass);
	histo_MCp->Fill(mc_p);
	histo_MCphi->Fill(mc_phi);
	histo_MCpid->Fill(mc_pid);
	histo_MCpt->Fill(mc_pt);
	histo_MCpx->Fill(mc_px);
	histo_MCpy->Fill(mc_py);
	histo_MCpz->Fill(mc_pz);
	  
      } // end loop on particles
      // if ( evt_iter/100 == float(evt_iter)/100. ) cout<< "numero quarks da Z " << info.size() << endl;  
    
      if ( info.size() == 2 &&
	   ( info[0].pid() == -info[1].pid()) ) {
	
	double e0 = sqrt(info[0].mass()*info[0].mass() +info[0].p()*info[0].p() ); 
	double e1 = sqrt(info[1].mass()*info[1].mass() +info[1].p()*info[1].p() ); 
	double minv = sqrt( pow(e0+e1,2) - pow(info[0].px()+info[1].px(),2) - pow(info[0].py()+info[1].py(),2) - pow(info[0].pz()+info[1].pz(),2) );

	for( unsigned int info_iter =0; info_iter < info.size(); ++info_iter ) {
	  histo_Qpt->Fill( info[info_iter].pt() );
	  histo_Qphi->Fill( info[info_iter].phi() );
	  histo_Qeta->Fill(info[info_iter].eta() );
	  histo_Qmass->Fill( info[info_iter].mass() );
	  histo_Qp->Fill( info[info_iter].p() );
	  histo_Qpx->Fill( info[info_iter].px() );
	  histo_Qpy->Fill( info[info_iter].py() );
	  histo_Qpz->Fill( info[info_iter].pz() );
	}

      } // end if (info.size() == 2 && ( info[0].pid() != info[1].pid() ) )

    } // end if ( jets_number >= 2 && muon_number >= 2 && mcparticle_number != 0 ) 

  }  // end loop on events
  
  std::cout << "number of Z->mu+mu-: " << leptonicZcounter << std::endl;
  for (unsigned int index = 0; index < nCuts; ++index ) {
    std::cout << "nCuts: " << index << std::endl;
    std::cout << "number of 'good' jet matched to light quark: " << lightQcounter[index]    << std::endl;
    std::cout << "number of 'good' jet matched to c quark: "     << cQcounter[index]        << std::endl;
    std::cout << "number of 'good' jet matched to b quark: "     << bQcounter[index]        << std::endl;
    std::cout << "number of Z->jj: "                             << hadronicZcounter[index] << std::endl;
    std::cout << "number of 'H->ZZ->mumujj: "                    << Hcounter[index]         << std::endl;
  }

  //histo_muon_pt->GetXaxis()->SetTitle("muon p_{T} (GeV/c)");

 TCanvas * mcparticle_can = new TCanvas("c1","c1",1000,800);
 mcparticle_can->Divide(1,2);
 mcparticle_can->cd(1);
 histo_MCphi->GetXaxis()->SetTitle("#phi (radianti)");
 histo_MCphi->Draw();
 mcparticle_can->cd(2);
 histo_MCeta->GetXaxis()->SetTitle("#eta");
 histo_MCeta->Draw();
 
  outputfile->Write();
  
  return 0;
};

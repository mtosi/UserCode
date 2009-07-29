// -*- C++ -*-
//
// Package:    VBFHZZllbbAnalyzer
// Class:      VBFHZZllbbAnalyzer
// 
/**\class VBFHZZllbbAnalyzer VBFHZZllbbAnalyzer.cc HiggsAnalysis/VBFHiggsToZZto2l2b/src/VBFHZZllbbAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Mia TOSI
//         Created:  Mon Feb  2 17:31:44 CET 2009
// $Id: VBFHZZllbbAnalyzer.h,v 1.2 2009/04/28 16:45:01 tosi Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//  to access TFileService within a Module
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/VBFHZZllbbElectronTrackIsolationAlgos.h"

// class decleration
//

class VBFHZZllbbAnalyzer : public edm::EDAnalyzer {
   public:
      explicit VBFHZZllbbAnalyzer(const edm::ParameterSet&);
      ~VBFHZZllbbAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

  // ----------member data ---------------------------
  TTree *mytree_;
  // isolation algorithm
  VBFHZZllbbElectronTrackIsolationAlgos * eleTrkIsoAlgo_;

  // input parameter and tag label
  int           whichSim_;
  edm::InputTag vertexLabel_;
  edm::InputTag trackLabel_;
  edm::InputTag muonLabel_;
  edm::InputTag electronLabel_;
  bool          eleTrkIsoAlgoFlag_;
  edm::InputTag metLabel_;
  edm::InputTag caloJetLabel_;
  std::string   corIC5CaloJetsWithBTagLabel_;
  std::string   corIC5PFJetsWithBTagLabel_;
  bool          corIC5PFJetsWithBTagFlag_;
  edm::InputTag genParticleLabel_;
  edm::InputTag genJetLabel_;
  edm::InputTag genMetLabel_;

  double        leptonPtCut_;	          
  double        jetEtCut_;
  double        jetPartonDeltaR2Cut_;
  double        jetLeptonDeltaRCut_;
  double        jetEMfracCut_;

  
  unsigned nleptons_;
  unsigned njets_;
  unsigned nZjets_;
  unsigned nforwardjets_;

  int eventcounter_;
  int badeventcounter_;
  int goodLeptonicEventCounter_;
  int twoJetsAboveEtCutEventsCounter_;
  int twoJetsAboveEtCut_DeltaR2016EventsCounter_;
  int twoJetsAboveEtCut_DeltaR2CutEventsCounter_;
  int fourJetsAboveEtCutEventsCounter_;

  int nbin_;

  edm::Service<TFileService> fs;
  
  TH1D * h_jet_N, * h_jet_et, * h_jet_pt, * h_jet_eta, * h_jet_phi, * h_jet_bDiscr;
  TH1D * h_jet_coret, * h_jet_corpt;
  TH1D * h_jetmuon_deltaR, * h_jetmet_deltaR, * h_jetelectron_deltaR;
  TH1D * h_bjet_N, * h_bjetmuon_deltaR, * h_bjetmet_deltaR, * h_bjetelectron_deltaR;
  TH1D * h_muon_N, * h_muon_pt, * h_muon_eta, * h_muon_phi;
  TH1D * h_dimuon_mass, * h_dimuon_pt, * h_dimuon_eta, * h_dimuon_phi;
  TH1D * h_dimuon_deltaEta, * h_dimuon_deltaR, * h_dimuon_deltaPhi;
  TH1D * h_dimuon_pz1pz2, * h_dimuon_eta1eta2;
  TH2D * h_dimuon_massVSpt, * h_dimuon_massVSeta, * h_dimuon_massVSphi;
  TH2D * h_dimuon_massVSdeltaEta, * h_dimuon_massVSdeltaR, * h_dimuon_massVSdeltaPhi;
  TH2D * h_dimuon_massVSpz1pz2, * h_dimuon_massVSeta1eta2;
  TH1D * h_electron_N, * h_electron_pt, * h_electron_eta, * h_electron_phi;
  TH1D * h_met_N, * h_met_pt, * h_met_eta, * h_met_phi;

};

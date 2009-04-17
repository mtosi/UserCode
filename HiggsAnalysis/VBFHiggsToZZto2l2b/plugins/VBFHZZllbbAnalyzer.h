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
// $Id$
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

#include "TH1F.h"

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
  edm::InputTag genParticleLabel_;
  edm::InputTag caloJetLabel_;
  edm::InputTag muonLabel_;
  edm::InputTag electronLabel_;
  edm::InputTag bTagLabel_;
  edm::InputTag metLabel_;
  double ptMax_;
  
  edm::Service<TFileService> fs;
  
  TH1D * h_jet_pt,  * h_jet_eta,  * h_jet_phi, * h_jet_bDiscr;
  TH1D * h_jetmuon_deltaR, * h_jetmet_deltaR, * h_jetelectron_deltaR;
  TH1D * h_bjetmuon_deltaR, * h_bjetmet_deltaR, * h_bjetelectron_deltaR;
  TH1D * h_muon_pt, * h_muon_eta, * h_muon_phi;
  TH1D * h_electron_pt, * h_electron_eta, * h_electron_phi;
  TH1D * h_met_pt, * h_met_eta, * h_met_phi;
};

// -*- C++ -*-
//
// Package:    VBFHZZllbbRAWCORBTAGtesting
// Class:      VBFHZZllbbRAWCORBTAGtesting
// 
/**\class VBFHZZllbbRAWCORBTAGtesting VBFHZZllbbRAWCORBTAGtesting.cc HiggsAnalysis/VBFHiggsToZZto2l2b/src/VBFHZZllbbRAWCORBTAGtesting.cc

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

class VBFHZZllbbRAWCORBTAGtesting : public edm::EDAnalyzer {
   public:
      explicit VBFHZZllbbRAWCORBTAGtesting(const edm::ParameterSet&);
      ~VBFHZZllbbRAWCORBTAGtesting();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

  // ----------member data ---------------------------
  edm::InputTag caloJetLabel_;
  edm::InputTag corrCaloJetLabel_;
  edm::InputTag muonLabel_;
  edm::InputTag electronLabel_;
  edm::InputTag metLabel_;
  //  edm::InputTag bTagLabel_;
  std::vector<edm::ParameterSet> bTagConfigLabel_;
  std::string jetCorrectionServiceLabel_;

  std::vector<edm::InputTag> bJetTagInputTags_;

  edm::Service<TFileService> fs;
  
  TH1D * h_jet_pt,  * h_jet_eta,  * h_jet_phi, * h_jet_emFrac;
  TH1D * h_jetmuon_deltaR, * h_jetmet_deltaR, * h_jetelectron_deltaR;
  TH1D * h_corjet_pt,  * h_corjet_eta,  * h_corjet_phi, * h_corjet_emFrac;
  TH1D * h_corjetmuon_deltaR, * h_corjetmet_deltaR, * h_corjetelectron_deltaR;

  std::vector<TH1D*> hvec_bjet_pt, hvec_bjet_eta, hvec_bjet_phi, hvec_bjet_emFrac, hvec_bjet_bDiscr;
  std::vector<TH1D*> hvec_bjetmuon_deltaR, hvec_bjetelectron_deltaR, hvec_bjetmet_deltaR;

  TH1D * h_coronflyjet_pt;
  TH1D * h_bcoronflyjet_pt;


};

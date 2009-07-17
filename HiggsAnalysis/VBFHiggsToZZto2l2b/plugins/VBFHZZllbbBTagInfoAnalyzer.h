// -*- C++ -*-
//
// Package:    VBFHZZllbbBTagInfoAnalyzer
// Class:      VBFHZZllbbBTagInfoAnalyzer
// 
/**\class VBFHZZllbbBTagInfoAnalyzer VBFHZZllbbBTagInfoAnalyzer.cc HiggsAnalysis/VBFHiggsToZZto2l2b/src/VBFHZZllbbBTagInfoAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Mia TOSI
//         Created:  Mon Feb  2 17:31:44 CET 2009
// $Id: VBFHZZllbbBTagInfoAnalyzer.h,v 1.2 2009/07/06 13:15:48 tosi Exp $
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

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"


#include "TH1F.h"

// class decleration
//

class VBFHZZllbbBTagInfoAnalyzer : public edm::EDAnalyzer {
   public:
      explicit VBFHZZllbbBTagInfoAnalyzer(const edm::ParameterSet&);
      ~VBFHZZllbbBTagInfoAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze (const edm::Event&, const edm::EventSetup&);
      virtual void endJob  () ;

  double tracksInvariantMass(const reco::TrackRefVector & selectedTracks ) const;

  // ----------member data ---------------------------
  edm::InputTag impactParameterTagInfosLabel_;
  edm::InputTag secondaryVertexTagInfosLabel_;

  edm::Service<TFileService> fs;

  TH1D * muonpvX_;
  TH1D * muonpvY_;
  TH1D * muonpvZ_;

  TH1D * pvX_;
  TH1D * pvY_;
  TH1D * pvZ_;
  TH1D * BSpvX_;
  TH1D * BSpvY_;
  TH1D * BSpvZ_;
  TH1D * PXpvX_;
  TH1D * PXpvY_;
  TH1D * PXpvZ_;

  TH1D * tagpvX_;
  TH1D * tagpvY_;
  TH1D * tagpvZ_;
  TH1D * tagsvX_;
  TH1D * tagsvY_;
  TH1D * tagsvZ_;
  TH1D * distanceTo1track_;
  TH1D * distanceToJetAxis_;
  TH1D * ip2d_;
  TH1D * ip3d_;
  TH1D * ip2dSig_;
  TH1D * ip3dSig_;

  TH1D * selTrkNum_;
  TH1D * selTrkNumS1_;
  TH1D * selTrkNumS2_;
  TH1D * selTrkNumS3_;
  TH1D * selTrkSumPt_;
  TH1D * selTrkSumPtS1_;
  TH1D * selTrkSumPtS2_;
  TH1D * selTrkSumPtS3_;
  TH1D * selTrkMass_;
  TH1D * selTrkMassS1_;
  TH1D * selTrkMassS2_;
  TH1D * selTrkMassS3_;

};

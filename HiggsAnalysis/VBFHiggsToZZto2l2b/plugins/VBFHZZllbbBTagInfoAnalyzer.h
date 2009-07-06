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
// $Id: VBFHZZllbbBTagInfoAnalyzer.h,v 1.1 2009/05/14 10:52:01 tosi Exp $
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
  
};

// -*- C++ -*-
//
// Package:    VBFHZZllbbBhadronReconstruction
// Class:      VBFHZZllbbBhadronReconstruction
// 
/**\class VBFHZZllbbBhadronReconstruction VBFHZZllbbBhadronReconstruction.cc HiggsAnalysis/VBFHiggsToZZto2l2b/src/VBFHZZllbbBhadronReconstruction.cc

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

// class decleration
//

class VBFHZZllbbBhadronReconstruction : public edm::EDAnalyzer {
   public:
      explicit VBFHZZllbbBhadronReconstruction(const edm::ParameterSet&);
      ~VBFHZZllbbBhadronReconstruction();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

  // ----------member data ---------------------------
  std::string corJetWithBTagLabel_;

  double MB_;   // B hadron mass (GeV/c^2)
  double MB2_;  // B hadron mass (GeV/c^2)
  


};

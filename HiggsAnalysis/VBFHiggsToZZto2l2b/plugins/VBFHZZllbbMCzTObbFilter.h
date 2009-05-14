#ifndef VBFHZZllbbMCzTObbFilter_h
#define VBFHZZllbbMCzTObbFilter_h
// -*- C++ -*-
//
// Package:    VBFHZZllbbMCzTObbFilter
// Class:      VBFHZZllbbMCzTObbFilter
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


// class decleration
//
class VBFHZZllbbMCzTObbFilter : public edm::EDFilter {
   public:
      explicit VBFHZZllbbMCzTObbFilter(const edm::ParameterSet&);
      ~VBFHZZllbbMCzTObbFilter();


      virtual bool filter(edm::Event&, const edm::EventSetup&);
   private:
      // ----------member data ---------------------------

  edm::InputTag genParticleLabel_;
  int signal_;

};

#endif


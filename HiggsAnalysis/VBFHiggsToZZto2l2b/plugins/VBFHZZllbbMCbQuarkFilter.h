#ifndef VBFHZZllbbMCbQuarkFilter_h
#define VBFHZZllbbMCbQuarkFilter_h
// -*- C++ -*-
//
// Package:    VBFHZZllbbMCbQuarkFilter
// Class:      VBFHZZllbbMCbQuarkFilter
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
class VBFHZZllbbMCbQuarkFilter : public edm::EDFilter {
   public:
      explicit VBFHZZllbbMCbQuarkFilter(const edm::ParameterSet&);
      ~VBFHZZllbbMCbQuarkFilter();


      virtual bool filter(edm::Event&, const edm::EventSetup&);
   private:
      // ----------member data ---------------------------

  edm::InputTag genParticleLabel_;
  int signal_;

};

#endif


#ifndef VBFHZZllbbMCprocessFilter_h
#define VBFHZZllbbMCprocessFilter_h
// -*- C++ -*-
//
// Package:    VBFHZZllbbMCprocessFilter
// Class:      VBFHZZllbbMCprocessFilter
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
class VBFHZZllbbMCprocessFilter : public edm::EDFilter {
   public:
      explicit VBFHZZllbbMCprocessFilter(const edm::ParameterSet&);
      ~VBFHZZllbbMCprocessFilter();


      virtual bool filter(edm::Event&, const edm::EventSetup&);
   private:
      // ----------member data ---------------------------
      
  int  whichSim_;
  int  signal_;
  bool vbfFlag_;
};

#endif


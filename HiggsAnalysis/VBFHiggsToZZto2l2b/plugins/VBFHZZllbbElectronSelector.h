#ifndef  VBFHZZllbbElectronSelector_h
#define  VBFHZZllbbElectronSelector_h

/**\class VBFHZZllbbElectronSelector
 *
 *
 * Original Author:  Dominique Fortin
 * Modified by Nicola De Filippis
 * Refine electron collection to begin with
 *
 */

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

class VBFHZZllbbElectronSelector : public edm::EDProducer {
 public:
  explicit VBFHZZllbbElectronSelector(const edm::ParameterSet& );
  ~VBFHZZllbbElectronSelector();

 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  edm::InputTag sourceLabel_;
  edm::InputTag sourceIDLabel_;

  float sourcePtMinCut_;
  float sourceEtaMaxCut_;
  
  int counterelectron,counterelectronbefore;
 
};
  
#endif


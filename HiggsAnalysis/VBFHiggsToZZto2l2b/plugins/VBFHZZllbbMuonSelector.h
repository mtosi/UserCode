#ifndef  VBFHZZllbbMuonSelector_h
#define  VBFHZZllbbMuonnSelector_h

/**\class VBFHZZllbbMuonSelector
 *
 * Refine muon collection to begin with
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

class VBFHZZllbbMuonSelector : public edm::EDProducer {
 public:
  explicit VBFHZZllbbMuonSelector(const edm::ParameterSet& );
  ~VBFHZZllbbMuonSelector();

 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  edm::InputTag sourceLabel_;

  double sourceMinPtBarrelCut_;
  double sourceMinPtEndcapCut_;
  double sourceMinPEndcapCut_;
  double sourceMaxEtaCut_;

};

#endif

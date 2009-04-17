#ifndef VBFHZZllbbElectronIsolationProducer_h
#define VBFHZZllbbElectronIsolationProducer_h

/**\class VBFHZZllbbElectronIsolationProducer
 *
 * Compute isolation for cones around electron candidates
 */
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"


class VBFHZZllbbElectronIsolationProducer : public edm::EDProducer {

 public:
  explicit VBFHZZllbbElectronIsolationProducer(const edm::ParameterSet&);
  ~VBFHZZllbbElectronIsolationProducer();

 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  edm::InputTag electronLabel_;
  edm::InputTag electronVetoLabel_;
  edm::InputTag muonsLabel_;
  edm::InputTag tracksLabel_;
  double isolationConeCut_;
  double isolationVetoCut_;
  double isolationCut_;
  //  double ptMin;
  // doudle maxDz;

};

#endif

#ifndef VBFHZZllbbMuonIsolationProducer_h
#define VBFHZZllbbMuonIsolationProducer_h

/**\class VBFHZZllbbMuonIsolationProducer
 *
 * Compute isolation for cones around electron candidates
 */
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

class VBFHZZllbbMuonIsolationProducer : public edm::EDProducer {

 public:
  explicit VBFHZZllbbMuonIsolationProducer(const edm::ParameterSet&);
  ~VBFHZZllbbMuonIsolationProducer();

 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  edm::InputTag muonLabel_;
  edm::InputTag electronLabel_;
  edm::InputTag trackLabel_;
  edm::InputTag theECALIsoDepositLabel_;    // EM calorimeter Isolation deposit label
  edm::InputTag theHCALIsoDepositLabel_;    // Hadron calorimeter Isolation deposit label
  edm::InputTag theHOCALIsoDepositLabel_;   // Outer calorimeter Isolation deposit label
  edm::InputTag theTrackerIsoDepositLabel_; // Tracker Isolation deposit label 
  double isoCut_;
  double vetoConeSize_;
  double mainConeSize_;
  //  double ptMin;
  // double maxDz;

};

#endif

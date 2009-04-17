#ifndef VBFHZZllbbPRESELECTION_h
#define VBFHZZllbbPRESELECTION_h

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

using namespace edm;
using namespace std;

class VBFHZZllbbPreSelection : public edm::EDFilter {
  
 public:
  // Constructor
  explicit VBFHZZllbbPreSelection(const edm::ParameterSet&);

  // Destructor
  ~VBFHZZllbbPreSelection();

  /// Get event properties to send to builder to fill seed collection
  virtual bool filter(edm::Event&, const edm::EventSetup& );


 private:
  int eventsCounter_, selectedEventsCounter_;

  double tightLeptonMinPt_;
  double softLeptonMinPt_;
  int tightLeptonMinNumber_;
  int softLeptonMinNumber_;
  double tightJetMinPt_;
  double softJetMinPt_;
  int tightJetMinNumber_;
  int softJetMinNumber_;

  // Reco samples
  edm::InputTag electronLabel_;
  edm::InputTag muonLabel_;
  edm::InputTag metLabel_;
  edm::InputTag jetLabel_;
  std::string   corJetsWithBTagLabel_;
  edm::InputTag mcParticleLabel_;
  edm::InputTag genJetLabel_;
  edm::InputTag genMetLabel_;

};

#endif

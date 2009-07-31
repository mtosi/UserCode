#ifndef HiggsAnalysis_VBFHZZllbbSkim
#define HiggsAnalysis_VBFHZZllbbSkim

/* \class VBFHZZllbbSkim
 *
 *
 * Filter to select 4 lepton events based on the
 * 1 or 2 electron or 1 or 2 muon HLT trigger, 
 * and four leptons (no flavour requirement).
 * No charge requirements are applied on event.
 *
 * \author Dominique Fortin - UC Riverside
 *
 */

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

class VBFHZZllbbSkim : public edm::EDFilter {
  
 public:
  // Constructor
  explicit VBFHZZllbbSkim(const edm::ParameterSet&);

  // Destructor
  ~VBFHZZllbbSkim();

  /// Get event properties to send to builder to fill seed collection
  virtual bool filter(edm::Event&, const edm::EventSetup& );


 private:
  int nEvents_, nSelectedEvents_;

  bool debug_;
  double tightLeptonMinPt_;
  double softLeptonMinPt_;
  int tightLeptonMinNumber_;
  int softLeptonMinNumber_;

  // Reco samples
  edm::InputTag muonLabel_;
  edm::InputTag electronLabel_;
};

#endif

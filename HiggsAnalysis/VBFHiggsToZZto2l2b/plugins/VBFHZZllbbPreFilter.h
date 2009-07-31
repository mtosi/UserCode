#ifndef HiggsAnalysis_VBFHZZllbbPreFilter
#define HiggsAnalysis_VBFHZZllbbPreFilter

/* \class HiggsTo4LeptonsSkim
 *
 *
 * Filter to select 4 lepton events (4e, 4mu, 2e2mu) within
 * fiducial volume (|eta| < 2.4)
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

class VBFHZZllbbPreFilter : public edm::EDFilter {
  
 public:
  // Constructor
  explicit VBFHZZllbbPreFilter(const edm::ParameterSet&);

  // Destructor
  ~VBFHZZllbbPreFilter();

  /// Get event properties to send to builder to fill seed collection
  virtual bool filter(edm::Event&, const edm::EventSetup& );


 private:
  int nEvents_, nSelectedEvents_;
  
  bool debug_;
  int  leptonFlavour_;

  edm::InputTag MCParticleLabel_;

};

#endif

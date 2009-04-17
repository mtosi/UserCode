#ifndef HiggsAnalysis_VBFHZZllbbSkimEff
#define HiggsAnalysis_VBFHZZllbbSkimEff

/* \class HiggsTo4LeptonsSkimEff
 *
 * EDAnalyzer to study the HLT and skim efficiency for signal
 * A preselection on the generaged event is built in
 *
 * \author Dominique Fortin - UC Riverside
 *
 */

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

using namespace edm;
using namespace std;

class VBFHZZllbbSkimEff : public edm::EDAnalyzer {
  
 public:
  // Constructor
  explicit VBFHZZllbbSkimEff(const edm::ParameterSet&);

  // Destructor
  ~VBFHZZllbbSkimEff();

  /// Get event properties to send to builder to fill seed collection
  virtual void analyze(const edm::Event&, const edm::EventSetup& );


 private:
  // input tag
  bool debug_;
  double tightLeptonMinPt_;
  double softLeptonMinPt_;
  int tightLeptonMinNumber_;
  int softLeptonMinNumber_;

  // event counters
  int nEvents_;
  int nSelTwoE_, nSelTwoM_, nSelTwoL_, nSelTau_;
  int nTwoE_, nTwoM_, nTwoL_, nTau_;
  int nOutE_, nOutM_, nOutTau_;

  // gen particle
  edm::InputTag MCParticleLabel_;
  // reco samples
  edm::InputTag muonLabel_;
  edm::InputTag electronLabel_;
};

#endif

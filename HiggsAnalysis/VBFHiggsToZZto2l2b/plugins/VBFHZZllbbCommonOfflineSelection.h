#ifndef VBFHZZllbbCommonOfflineSelection_h
#define VBFHZZllbbCommonOfflineSelection_h

/* \class VBFHZZllbbCommonOfflineSelection
 *
 *
 * Analysis selection:
 * - tight isolation
 * - vertexing
 */

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"


// Class declaration
class VBFHZZllbbCommonOfflineSelection : public edm::EDProducer {
  
 public:
  // Constructor
  explicit VBFHZZllbbCommonOfflineSelection(const edm::ParameterSet&);

  // Destructor
  ~VBFHZZllbbCommonOfflineSelection();

 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  bool match(double mass, double pt, int charge, const reco::CandidateCollection *c1Coll);

  std::string decaychannel;
  bool useBestCandidate;

  edm::InputTag leptonLabel_;
  edm::InputTag bestCandidatesLeptonsLabel_;
  edm::InputTag leptonMapLabel_;
  edm::InputTag leptonIsoVarLabel_;
  std::vector<double> leptonIsoVarCut_;

  edm::InputTag leptonVertexLabel_;
  edm::InputTag leptonVertexMapLabel_;
  std::vector<double> vertexVarCut_;

};

#endif

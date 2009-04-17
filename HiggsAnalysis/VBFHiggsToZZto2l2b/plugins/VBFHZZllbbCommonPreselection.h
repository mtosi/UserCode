#ifndef VBFHZZllbbCommonPreselection_h
#define VBFHZZllbbCommonPreselection_h

/* \class VBFHZZllbbCommonPreSelection
 *
 *
 * Analysis preselection:
 * m_ll > 12 GeV
 * m_H  > 100 GeV
 *
 */

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

// Class declaration
class VBFHZZllbbCommonPreselection : public edm::EDProducer {
  
 public:
  // Constructor
  explicit VBFHZZllbbCommonPreselection(const edm::ParameterSet&);

  // Destructor
  ~VBFHZZllbbCommonPreselection();

 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  std::string   decaychannel;
  edm::InputTag hLabel_;
  edm::InputTag leptonLabel_;
  edm::InputTag zToLLLabel_;
  edm::InputTag leptonLooseIsolLabel_;

  double llMassCut_;
  double fourBodyMassCut_;
  int  nLeptonCut_;
  int  numberOfLeptonCombsCut_;
  int  numberOf2l2bCombsCut_;
  int  nLooseLeptonCut_;

};

#endif

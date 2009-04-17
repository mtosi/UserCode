// -*- C++ -*-
//
// Package:    VBFHZZllbbCorJetWithBTagProducer
// Class:      VBFHZZllbbCorJetWithBTagProducer
// 
/**\class VBFHZZllbbCorJetWithBTagProducer VBFHZZllbbCorJetWithBTagProducer.cc HiggsAnalysis/VBFHZZllbbCorJetWithBTagProducer/src/VBFHZZllbbCorJetWithBTagProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Mia TOSI
//         Created:  Wed Mar 18 16:33:51 CET 2009
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class decleration
//
class VBFHZZllbbCorJetWithBTagProducer : public edm::EDProducer {
public:
  
  explicit VBFHZZllbbCorJetWithBTagProducer (const edm::ParameterSet& fParameters);
  virtual ~VBFHZZllbbCorJetWithBTagProducer ();
  virtual void produce(edm::Event&, const edm::EventSetup&);
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void endJob() ;
  
  // ----------member data ---------------------------
  // Jet algorithm: it can be any Calo or PF algorithm
  std::vector<edm::ParameterSet> bTagConfigLabel_;
  std::string jetCorrectionService_;
  
  std::vector<edm::InputTag> bJetTagInputTags_;
  
};


// -*- C++ -*-
//
// Package:    VBFHZZllbbCorJetWithBTagAnalyzer
// Class:      VBFHZZllbbCorJetWithBTagAnalyzer
// 
/**\class VBFHZZllbbCorJetWithBTagAnalyzer VBFHZZllbbCorJetWithBTagAnalyzer.cc HiggsAnalysis/VBFHiggsToZZto2l2b/src/VBFHZZllbbCorJetWithBTagAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Mia TOSI
//         Created:  Mon Feb  2 17:31:44 CET 2009
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbCorJetWithBTagAnalyzer.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/CorJetWithBTagDiscr.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace vbfhzz2l2b;

// class decleration
//

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
VBFHZZllbbCorJetWithBTagAnalyzer::VBFHZZllbbCorJetWithBTagAnalyzer(const edm::ParameterSet& iConfig) :
  corJetWithBTagLabel_( iConfig.getUntrackedParameter<std::string>("corJetWithBTagLabel") )
{
   //now do what ever initialization is needed

}

VBFHZZllbbCorJetWithBTagAnalyzer::~VBFHZZllbbCorJetWithBTagAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
VBFHZZllbbCorJetWithBTagAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<vbfhzz2l2b::CorJetWithBTagDiscrCollection> corJetWithBTagHandle;
  //  edm::Handle<vbfhzz2l2b::CorJetWithBTagsDiscrCollection> corJetWithBTagHandle;
  iEvent.getByLabel(corJetWithBTagLabel_,"corJetWithBTagDiscr",corJetWithBTagHandle);
  std::cout << "corJetWithBTagHandle->size(): " << corJetWithBTagHandle->size() << std::endl;

  int jetIndex = 0;    
  std::vector<reco::JetBaseRef> jets = vbfhzz2l2b::CorJetBTagDiscrAssociation::allJets(*corJetWithBTagHandle);
  for ( std::vector<reco::JetBaseRef>::const_iterator jet = jets.begin(); 
	jet != jets.end();
	++jet, jetIndex++ ) {
    std::cout << "jetIndex: " << jetIndex << std::endl;
    
    double corrEt = (*corJetWithBTagHandle)[*jet].corEt_;
    double uncorrEt = (*jet)->et();
    std::cout << "uncorrEt: " << uncorrEt << std::endl;
    std::cout << "corrEt: " << corrEt << std::endl;
    double emFrac = (dynamic_cast<const reco::CaloJet*>(&**jet))->emEnergyFraction();
    std::cout << "emFrac: " << emFrac << std::endl;

  }

  /*
  vbfhzz2l2b::CorJetWithBTagDiscrCollection::const_iterator corJetWithBTag = corJetWithBTagHandle->begin();
  int index = 0;
  for ( ; corJetWithBTag != corJetWithBTagHandle->end(); ++corJetWithBTag, index++ ) {
    std::cout << "index: " << index;
    double uncorEt = corJetWithBTag->corrEt_;
    std::cout << " --> uncorEt: " << uncorEt;
    double corEt = corJetWithBTag->et();
    std::cout << " <-- corEt: " << corEt;
//    //    std::cout << " discriminatorsVec.size: " << (corJetWithBTag->discriminators())->size() << std::endl;
//    const reco::Jet * refenrenceJet = corJetWithBTag->refJet();
//    
//    double emFrac = (dynamic_cast<const reco::CaloJet*>(refenrenceJet))->emEnergyFraction();
//    std::cout << " === emFrac: " << emFrac << std::endl;
  }
  */
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
VBFHZZllbbCorJetWithBTagAnalyzer::beginJob(const edm::EventSetup&)
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
VBFHZZllbbCorJetWithBTagAnalyzer::endJob() {
}

DEFINE_FWK_MODULE(VBFHZZllbbCorJetWithBTagAnalyzer);


/* \class VBFHZZllbbSkim 
 *
 * Consult header file for description
 *
 * author:  Dominique Fortin - UC Riverside
 *
 */


// system include files
#include "HiggsAnalysis/Skimming/interface/VBFHZZllbbSkim.h"

// User include files
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Muons:
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"

// Electrons
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

// C++
#include <iostream>
#include <vector>

using namespace std;
using namespace edm;
using namespace reco;

// Constructor
VBFHZZllbbSkim::VBFHZZllbbSkim(const edm::ParameterSet& iConfig) :

  // local debug flag
  debug_( iConfig.getParameter<bool>("VBFHZZllbbDebug") ),
  // reconstructed objects
  muonLabel_(     iConfig.getParameter<edm::InputTag>("muonLabel"    ) ),
  electronLabel_( iConfig.getParameter<edm::InputTag>("electronLabel") ),
  
  // minimum Pt for leptons for skimming
  tightLeptonMinPt_( iConfig.getParameter<double>("tightMinimumPt" ) ),
  softLeptonMinPt_(  iConfig.getParameter<double>("softMinimumPt"  ) ),
  tightLeptonMinNumber_( iConfig.getParameter<int>(   "tightLeptonMinimumNumber") ),
  softLeptonMinNumber_(  iConfig.getParameter<int>(   "softLeptonMinimumNumber" ) ) {

  std::cout << "[VBFHZZllbbSkim::VBFHZZllbbSkim]" << std::endl;
  // initializing event counter
  nEvents_         = 0;
  nSelectedEvents_ = 0;

}


// Destructor
VBFHZZllbbSkim::~VBFHZZllbbSkim() {
  std::cout << "*****************************************" << std::endl;
  std::cout << " number_events_read " << nEvents_ 
	    << " number_events_kept " << nSelectedEvents_ 
	    << " => efficiency:     " << ((double)nSelectedEvents_)/((double)nEvents_)
	    << std::endl;
  edm::LogVerbatim("VBFHZZllbbSkim") 
    << " number_events_read " << nEvents_
    << " number_events_kept " << nSelectedEvents_ 
    << " => efficiency:     " << ((double)nSelectedEvents_)/((double)nEvents_)
    << std::endl;
  std::cout << "*****************************************" << std::endl;

  std::cout << "[VBFHZZllbbSkim::~VBFHZZllbbSkim]" << std::endl;
}



// Filter event
bool VBFHZZllbbSkim::filter(edm::Event& iEvent, const edm::EventSetup& setup ) {

  nEvents_++;

  bool keepEvent     = false;
  int  nTightLeptons = 0;
  int  nSoftLeptons  = 0;
  

  // look at muons: get the muon track collection from the event
  edm::Handle<reco::MuonCollection> muonHandle;
  iEvent.getByLabel(muonLabel_,muonHandle);

  // loop over muon collections and count how many muons there are, 
  // and how many are above threshold
  reco::MuonCollection::const_iterator muon_itr = muonHandle->begin(); 
  for ( ; muon_itr != muonHandle->end(); ++muon_itr ) {      
    if ( muon_itr->pt() >= tightLeptonMinPt_ ) nTightLeptons++; 
    if ( muon_itr->pt() >= softLeptonMinPt_  ) nSoftLeptons++; 
  } 
  
  // look at electrons: get the electron track collection from the event
  edm::Handle<reco::GsfElectronCollection> electronHandle;
  iEvent.getByLabel(electronLabel_,electronHandle);

  reco::GsfElectronCollection::const_iterator electron_itr = electronHandle->begin();
  // Loop over electron collections and count how many muons there are, 
  // and how many are above threshold
  for ( ; electron_itr != electronHandle->end(); ++electron_itr ) {
    if ( electron_itr->pt() >= tightLeptonMinPt_ ) nTightLeptons++;
    if ( electron_itr->pt() >= softLeptonMinPt_  ) nSoftLeptons++; 
  }
  
  // make decision:
  if ( nTightLeptons >= tightLeptonMinNumber_ && 
       nSoftLeptons  >= softLeptonMinNumber_     ) keepEvent = true;

  if ( keepEvent ) nSelectedEvents_++;

  return keepEvent;
}



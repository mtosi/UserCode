
/* \class VBFHZZllbbPreSelection 
 *
 */


// system include files
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbPreSelection.h"

// user include files
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// muons:
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"

// electrons
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

// jets:
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

// met:
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"

// b-tagging
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/CorJetWithBTagDiscr.h"

// C++
#include <iostream>
#include <vector>

using namespace std;
using namespace edm;
using namespace reco;
using namespace vbfhzz2l2b;


// Constructor
VBFHZZllbbPreSelection::VBFHZZllbbPreSelection(const edm::ParameterSet& iConfig) :

  // reconstructed objects
  electronLabel_        ( iConfig.getParameter<edm::InputTag> ( "electronLabel"        ) ),
  muonLabel_            ( iConfig.getParameter<edm::InputTag> ( "muonLabel"            ) ),
  metLabel_             ( iConfig.getParameter<edm::InputTag> ( "metLabel"             ) ),
  jetLabel_             ( iConfig.getParameter<edm::InputTag> ( "jetLabel"             ) ),
  corJetsWithBTagLabel_ ( iConfig.getParameter<std::string>   ( "corJetsWithBTagLabel" ) ),
  mcParticleLabel_      ( iConfig.getParameter<edm::InputTag> ( "mcParticleLabel"      ) ),
  genJetLabel_          ( iConfig.getParameter<edm::InputTag> ( "genJetLabel"          ) ),
  genMetLabel_          ( iConfig.getParameter<edm::InputTag> ( "genMetLabel"          ) ),
  // minimum Pt for leptons
  tightLeptonMinPt_ ( iConfig.getParameter<double> ( "tightLeptonMinPt" ) ),
  softLeptonMinPt_  ( iConfig.getParameter<double> ( "softLeptonMinPt"  ) ),
  // minimum number for leptons
  tightLeptonMinNumber_ ( iConfig.getParameter<int>( "tightLeptonMinNumber" ) ),
  softLeptonMinNumber_  ( iConfig.getParameter<int>( "softLeptonMinNumber"  ) ),
  // minimum Pt for jets
  tightJetMinPt_ ( iConfig.getParameter<double> ( "tightJetMinPt" ) ),
  softJetMinPt_  ( iConfig.getParameter<double> ( "softJetMinPt"  ) ),
  // minimum number for jets
  tightJetMinNumber_ ( iConfig.getParameter<int>( "tightJetMinNumber" ) ),
  softJetMinNumber_  ( iConfig.getParameter<int>( "softJetMinNumber"  ) )
 {

  std::cout << "[VBFHZZllbbPreSelection::VBFHZZllbbPreSelection]" << std::endl;
  // initializing event counter
  eventsCounter_         = 0;
  selectedEventsCounter_ = 0;

}


// Destructor
VBFHZZllbbPreSelection::~VBFHZZllbbPreSelection() {
  std::cout << "*****************************************" << std::endl;
  std::cout << " number_events_read " << eventsCounter_ 
	    << " number_events_kept " << selectedEventsCounter_ 
	    << " => efficiency:     " << ((double)selectedEventsCounter_)/((double)eventsCounter_)
	    << std::endl;
  edm::LogVerbatim("VBFHZZllbbPreSelection") 
    << " number_events_read " << eventsCounter_
    << " number_events_kept " << selectedEventsCounter_ 
    << " => efficiency:     " << ((double)selectedEventsCounter_)/((double)eventsCounter_)
    << std::endl;
  std::cout << "*****************************************" << std::endl;

  std::cout << "[VBFHZZllbbPreSelection::~VBFHZZllbbPreSelection]" << std::endl;
}



// Filter event
bool VBFHZZllbbPreSelection::filter(edm::Event& iEvent, const edm::EventSetup& setup ) {

  eventsCounter_++;

  bool keepEvent = true;

  int nTightLeptons = 0;
  int nSoftLeptons  = 0;
  int nTightJets    = 0;
  int nSoftJets     = 0;
  
  // look at muons: get the muon track collection from the event
  edm::Handle<reco::MuonCollection> muonHandle;
  iEvent.getByLabel(muonLabel_,muonHandle);

  // loop over muon collections and count how many muons there are, 
  // and how many are above threshold
  for ( reco::MuonCollection::const_iterator muon_itr = muonHandle->begin(); 
	muon_itr != muonHandle->end(); 
	++muon_itr ) {      
    if ( muon_itr->pt() >= tightLeptonMinPt_ ) nTightLeptons++; 
    if ( muon_itr->pt() >= softLeptonMinPt_  ) nSoftLeptons++; 
  } 
  
  // look at electrons: get the electron track collection from the event
  edm::Handle<reco::GsfElectronCollection> electronHandle;
  iEvent.getByLabel(electronLabel_,electronHandle);

  // Loop over electron collections and count how many electrons there are, 
  // and how many are above threshold
  for ( reco::GsfElectronCollection::const_iterator electron_itr = electronHandle->begin();
	electron_itr != electronHandle->end(); 
	++electron_itr ) {
    if ( electron_itr->pt() >= tightLeptonMinPt_ ) nTightLeptons++;
    if ( electron_itr->pt() >= softLeptonMinPt_  ) nSoftLeptons++; 
  }

  // look at jets:
  edm::Handle<reco::CaloJetCollection> jetHandle ;
  iEvent.getByLabel (jetLabel_, jetHandle) ;

  // look at corrected jets:
  edm::Handle<vbfhzz2l2b::CorJetWithBTagDiscrCollection> corJetWithBTagHandle;
  iEvent.getByLabel(corJetsWithBTagLabel_,"corJetWithBTagDiscr",corJetWithBTagHandle);
  // Loop over jet collections and count how many jets there are, 
  // and how many are above threshold
  int jetIndex = 0;    
  std::vector<reco::JetBaseRef> jets = vbfhzz2l2b::CorJetBTagDiscrAssociation::allJets(*corJetWithBTagHandle);
  for ( std::vector<reco::JetBaseRef>::const_iterator jet_itr = jets.begin(); 
	jet_itr != jets.end();
	++jet_itr, jetIndex++ ) {
    std::cout << "jetIndex: " << jetIndex << std::endl;
    
    std::vector<double> discrVec = (*corJetWithBTagHandle)[*jet_itr].discrVec_;
    double corrEt = (*corJetWithBTagHandle)[*jet_itr].corEt_;
    double corrPt = (*corJetWithBTagHandle)[*jet_itr].corPt_;
    double uncorrEt = (*jet_itr)->et();
    double uncorrPt = (*jet_itr)->pt();
    std::cout << "uncorrEt: " << uncorrEt
	      << " --> corrEt: " << corrEt << std::endl;
    std::cout << "uncorrPt: " << uncorrPt
	      << " --> corrPt: " << corrPt << std::endl;
    for ( unsigned int index = 0; index != discrVec.size(); index++) 
      std::cout << "discr[" << index << "]: " << discrVec[index] << std::endl;

    double emFrac = (dynamic_cast<const reco::CaloJet*>(&**jet_itr))->emEnergyFraction();
    std::cout << "emFrac: " << emFrac << std::endl;
    if ( corrPt >= tightJetMinPt_ ) nTightJets++;
    else if ( corrEt >= tightJetMinPt_ ) std::cout << "SILLY BANANA!!! tightJetMinPt CUT VALUE IS REFERED TO ET" << std::cout;
    if ( corrPt >= softJetMinPt_  ) nSoftJets++; 
    else if ( corrEt >= softJetMinPt_ ) std::cout << "SILLY BANANA!!! softJetMinPt CUT VALUE IS REFERED TO ET" << std::cout;
  }

  // make decision:
  if ( nTightLeptons < tightLeptonMinNumber_ || 
       nSoftLeptons  < softLeptonMinNumber_     ) keepEvent = false;
  if ( nTightJets < tightJetMinNumber_ ||
       nSoftJets  < softJetMinNumber_     ) keepEvent = false;


  if ( keepEvent ) selectedEventsCounter_++;

  return keepEvent;
}



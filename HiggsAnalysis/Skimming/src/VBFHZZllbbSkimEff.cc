
/* \class HiggsTo4LeptonsSkimEff 
 *
 * Consult header file for description
 *
 * author:  Dominique Fortin - UC Riverside
 *
 */


// system include files
#include "HiggsAnalysis/Skimming/interface/VBFHZZllbbSkimEff.h"

// User include files
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Muons:
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"

// Electrons
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

// Candidate handling
#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "AnalysisExamples/AnalysisObjects/interface/PythiaParticleIndex.h"

// C++
#include <iostream>
#include <vector>

using namespace std;
using namespace edm;
using namespace reco;
using namespace anaobj;


// Constructor
VBFHZZllbbSkimEff::VBFHZZllbbSkimEff(const edm::ParameterSet& iConfig) :

  // local debug flag
  debug_( iConfig.getParameter<bool>("VBFHZZllbbDebug") ),
  // gen particle
  MCParticleLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "genParticleLabel" ) ),
  // reconstructed objects
  muonLabel_(     iConfig.getParameter<edm::InputTag>("muonLabel"    ) ),
  electronLabel_( iConfig.getParameter<edm::InputTag>("electronLabel") ),

  // minimum Pt for leptons for skimming
  tightLeptonMinPt_( iConfig.getParameter<double>("tightMinimumPt" ) ),
  softLeptonMinPt_(  iConfig.getParameter<double>("softMinimumPt"  ) ),
  tightLeptonMinNumber_( iConfig.getParameter<int>(   "tightLeptonMinimumNumber") ),
  softLeptonMinNumber_(  iConfig.getParameter<int>(   "softLeptonMinimumNumber" ) ) {

  // initializing event counter
  nEvents_  = 0;
  nSelTwoE_ = nSelTwoM_ = nSelTwoL_ = nSelTau_ = 0;
  nTwoE_    = nTwoM_    = nTwoL_    = nTau_    = 0;
  nOutE_ = nOutM_ = nOutTau_ = 0;

}


// Destructor
VBFHZZllbbSkimEff::~VBFHZZllbbSkimEff() {
  std::cout << "************************************************************" << std::endl;
  std::cout << "Number of events read " << nEvents_ << std::endl;
  std::cout << "*** Efficiency for the various subsamples *** " <<  std::endl;

  std::cout << "Two leptons: " 
	    << " preselected " << nTwoL_ << " [out: " << nOutM_+nOutE_ << "]"
	    << " kept        " << nSelTwoL_
	    << " efficiency  " << ((double)nSelTwoL_)/((double)nTwoL_) << std::endl;
  std::cout << "Two muons:   "
	    << " preselected " << nTwoM_ << " [out: " << nOutM_ << "]"        
	    << " kept        " << nSelTwoM_
	    << " efficiency  " << ((double)nSelTwoM_)/((double)nTwoM_) << std::endl;
  std::cout << "Two elecs:   "
	    << " preselected " << nTwoE_  << " [out: " << nOutE_ << "]"
	    << " kept        " << nSelTwoE_
	    << " efficiency  " << ((double)nSelTwoE_)/((double)nTwoE_) << std::endl;
  std::cout << "with taus:   "
	    << " preselected " << nTau_  << " [out: " << nOutTau_ << "]"
	    << " kept        " << nSelTau_
	    << " efficiency  " << ((double)nSelTau_)/((double)nTau_) << std::endl;
  std::cout << "************************************************************" << std::endl;


}



// Filter event
void VBFHZZllbbSkimEff::analyze(const edm::Event& iEvent, const edm::EventSetup& setup ) {

  nEvents_++;

  bool keepEvent = false;

  // First, pre-selection:
  int nMuon = 0;
  int nElec = 0;
  int nTau  = 0;
  int nOutMuon = 0;
  int nOutElec = 0;
  int nOutTau  = 0;

  bool isTwoE = false;
  bool isTwoM = false;
  bool isTwoL = false;
  bool isTau  = false;

  // get gen particle candidates 
  edm::Handle < GenParticleCollection > MCparticles;
  iEvent.getByLabel( MCParticleLabel_, MCparticles );

  GenParticleCollection::const_iterator MCparticle = MCparticles->begin();
  for ( ; MCparticle != MCparticles->end(); ++MCparticle ) {    
    if ( MCparticle->mother() != 0 ) {
      // mother is a Z
      if ( MCparticle->mother()->pdgId() == pythiaZ_ ) {
	// muons:
	if ( fabs(MCparticle->pdgId()) == pythiamu_ ) {
	  // in fiducial volume:
	  if ( MCparticle->pt() < 3. ) continue;
	  if ( fabs(MCparticle->eta()) <= 2.4 ) nMuon++;
	  else nOutMuon++;
	}
	// electrons:
	else if ( fabs(MCparticle->pdgId()) == pythiae_ ) {
	  // In fiducial volume:
	  if ( MCparticle->pt() < 3. ) continue;
	  if ( fabs(MCparticle->eta()) <= 2.5 ) nElec++;
	  else nOutElec++;
	}
	// taus:
	else if ( fabs(MCparticle->pdgId()) == pythiatau_ ) {
	  // In fiducial volume:
	  //  if ( MCparticle->pt() < 3. ) continue; // why do they not apply this request?!?
	  if ( fabs(MCparticle->eta()) <= 2.5 ) nTau++;
	}
      }
    }
  }

  if (nOutMuon > 0) nOutM_++;
  if (nOutElec > 0) nOutE_++;
  if (nOutTau  > 0) nOutTau_++;

  if (nMuon >= 2) {
    isTwoM = true;
    nTwoM_++;
  }
  if (nElec >= 2) {
    isTwoE = true;
    nTwoE_++;
  }
  if (nTau >= 2) {
    isTau = true;
    nTau_++;
  }  
  if (isTwoE || isTwoM) {
    isTwoL = true;
    nTwoL_++;
  }
  if (isTwoL) keepEvent = true;
  else 
    return;

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
  else keepEvent = false;

  if ( keepEvent ) {
    if (isTwoE) nSelTwoE_++;
    if (isTwoM) nSelTwoM_++;
    if (isTwoL) nSelTwoL_++;
    if (isTau)  nSelTau_++;
  }

}

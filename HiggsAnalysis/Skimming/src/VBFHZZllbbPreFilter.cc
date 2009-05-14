
/* \class HiggsTo4LeptonsPreFilter 
 *
 * Consult header file for description
 *
 * author:  Dominique Fortin - UC Riverside
 *
 */


// system include files
#include "HiggsAnalysis/Skimming/interface/VBFHZZllbbPreFilter.h"

// User include files
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Candidate handling
//#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/PythiaParticleIndex.h"

// C++
#include <iostream>
#include <vector>

using namespace std;
using namespace edm;
using namespace reco;
using namespace vbfhzz2l2b;


// Constructor
VBFHZZllbbPreFilter::VBFHZZllbbPreFilter(const edm::ParameterSet& iConfig) :
  MCParticleLabel_(iConfig.getUntrackedParameter<edm::InputTag>( "genParticleLabel" ) ),
  // Local Debug flag
  leptonFlavour_( iConfig.getParameter<int> ("VBFHZZllbbPreFilterLeptonFlavour") ) {
  
  // initializing event counter
  nSelectedEvents_ = 0;
  nEvents_         = 0;
}


// Destructor
VBFHZZllbbPreFilter::~VBFHZZllbbPreFilter() {
  std::cout << "*****************************************" << std::endl;
  std::cout << "number of events processed: " << nEvents_ << std::endl;
  std::cout << "number of events kept:      " << nSelectedEvents_ << std::endl;
  std::cout << " => efficiency:             " << ((double)nSelectedEvents_)/((double)nEvents_) << std::endl;
  std::cout << "*****************************************" << std::endl;
}


// Filter event
bool VBFHZZllbbPreFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  nEvents_++;

  bool keepEvent = false;

  bool twoL = false;
  bool twoE = false;
  bool twoM = false;

  int nElec = 0;
  int nMuon = 0;

  // get gen particle candidates 
  edm::Handle < GenParticleCollection > MCparticles;
  //  edm::Handle < CandidateCollection > MCparticles;
  iEvent.getByLabel( MCParticleLabel_, MCparticles );


  GenParticleCollection::const_iterator MCparticle = MCparticles->begin();
  //  CandidateCollection::const_iterator MCparticle = MCparticles->begin();
  for ( ; MCparticle != MCparticles->end(); ++MCparticle ) {    
    if ( MCparticle->mother() != 0 ) {
    // mother is a Z
      if ( MCparticle->mother()->pdgId() == pythiaZ_ ) {
	// muons:
	if ( fabs(MCparticle->pdgId()) == pythiamu_ ) 
	  // in fiducial volume:
	  if ( MCparticle->pt() >= 3. && fabs(MCparticle->eta()) <= 2.4 ) nMuon++;
	// electrons:
	if ( fabs(MCparticle->pdgId()) == pythiae_ ) 
	  // In fiducial volume:
	  if ( MCparticle->pt() >= 3. && fabs(MCparticle->eta()) <= 2.5 ) nElec++;
      }
    }
  }

  if (nElec >= 2)   twoE = true;
  if (nMuon >= 2)   twoM = true;
  if (twoE || twoM) twoL = true;
  
  if ( leptonFlavour_ == 0 && twoL ) keepEvent = true;    
  if ( leptonFlavour_ == 1 && twoM ) keepEvent = true;    
  if ( leptonFlavour_ == 2 && twoE ) keepEvent = true;    
  
  if (keepEvent ) nSelectedEvents_++;
  
  return keepEvent;
  
}



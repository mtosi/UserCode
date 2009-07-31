#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbMCzTObbFilter.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/PythiaParticleIndex.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/ParticlesMass.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/ParticlesCharge.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <iostream>

using namespace vbfhzz2l2b;
using namespace edm;
using namespace std;

//! constructor
VBFHZZllbbMCzTObbFilter::VBFHZZllbbMCzTObbFilter(const edm::ParameterSet& iConfig) :
  genParticleLabel_ ( iConfig.getParameter<edm::InputTag> ( "genParticleLabel" ) ),
  signal_           ( iConfig.getParameter<int>           ( "signal"           ) )  // 1:Signal,  0:Background
{
}


// ------------------------------------------------------------------------------------


//! destructor
VBFHZZllbbMCzTObbFilter::~VBFHZZllbbMCzTObbFilter()
{}


// ------------------------------------------------------------------------------------


//! filtering method
bool 
VBFHZZllbbMCzTObbFilter::filter (edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  bool zTObb = false;

  edm::Handle < reco::GenParticleCollection > genParticlesHandle;
  iEvent.getByLabel( genParticleLabel_, genParticlesHandle ) ;
  
  for ( reco::GenParticleCollection::const_iterator particle_itr = genParticlesHandle->begin();
	particle_itr != genParticlesHandle->end(); ++particle_itr ) {

    int tmpParticleId = particle_itr->pdgId();
    // Z boson
    if ( tmpParticleId != pythiaZ_ ) continue;
    int tmpMotherId = 0;
    if ( particle_itr->mother() != 0 ) tmpMotherId = particle_itr->mother()->pdgId();
    
    // Z to bbbar
    if ( tmpMotherId == tmpParticleId ) continue;
    int tmpZdaughterId = particle_itr->daughter(0)->pdgId();
    if ( fabs(tmpZdaughterId) == pythiab_ ) zTObb = true; 
    
  } // end loop over particle_itr
   
  unsigned int bQuarkExtraCounter = 0; 
  for ( reco::GenParticleCollection::const_iterator particle_itr = genParticlesHandle->begin();
	particle_itr != genParticlesHandle->end(); ++particle_itr ) {    
    int tmpParticleId = particle_itr->pdgId();
    if ( fabs(tmpParticleId) != pythiab_ ) continue;
    
    int tmpMotherId = 0;
    if ( particle_itr->mother() != 0 ) tmpMotherId = particle_itr->mother()->pdgId();
    if ( tmpMotherId == tmpParticleId && tmpMotherId == pythiaZ_ ) continue;
    
    bQuarkExtraCounter++;
  }

  if ( bQuarkExtraCounter ) 
    std::cout << "number of b quark not from Z decay: " << bQuarkExtraCounter << std::endl;
  
  return zTObb ;

}
	


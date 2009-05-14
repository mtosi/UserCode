#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbMCbQuarkFilter.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/PythiaParticleIndex.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/ParticlesMass.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/ParticlesCharge.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <iostream>

using namespace vbfhzz2l2b;
using namespace edm;
using namespace std;

//! constructor
VBFHZZllbbMCbQuarkFilter::VBFHZZllbbMCbQuarkFilter(const edm::ParameterSet& iConfig) :
  genParticleLabel_ ( iConfig.getParameter<edm::InputTag> ( "genParticleLabel" ) ),
  signal_           ( iConfig.getParameter<int>           ( "signal"           ) )  // 0:Signal,  1:Background
{
}


// ------------------------------------------------------------------------------------


//! destructor
VBFHZZllbbMCbQuarkFilter::~VBFHZZllbbMCbQuarkFilter()
{}


// ------------------------------------------------------------------------------------


//! filtering method
bool 
VBFHZZllbbMCbQuarkFilter::filter (edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  bool bQuark = false;

  edm::Handle < reco::GenParticleCollection > genParticlesHandle;
  iEvent.getByLabel( genParticleLabel_, genParticlesHandle ) ;
  
  for ( reco::GenParticleCollection::const_iterator particle_itr = genParticlesHandle->begin();
	particle_itr != genParticlesHandle->end(); ++particle_itr ) {
    int tmpParticleId = particle_itr->pdgId();
    if ( fabs(tmpParticleId) != pythiab_ ) continue;
    
    int tmpMotherId = 0;
    if ( particle_itr->mother() != 0 ) tmpMotherId = particle_itr->mother()->pdgId();
    if ( tmpMotherId == tmpParticleId ) continue;
    
    bQuark = true; 
  }

  return bQuark ;

}
	


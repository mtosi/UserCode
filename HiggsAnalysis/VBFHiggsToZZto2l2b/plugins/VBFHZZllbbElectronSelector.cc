/**\class VBFHZZllbbElectronSelector.cc
 *
 *
 */

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbElectronSelector.h"

#include "FWCore/Framework/interface/ESHandle.h"

// Electrons:
#include <DataFormats/EgammaCandidates/interface/GsfElectron.h>
#include <DataFormats/EgammaCandidates/interface/GsfElectronFwd.h>

// Candidate handling
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include <Math/VectorUtil.h>
#include <memory>
#include <vector>

using namespace edm;
using namespace std;
using namespace reco;


// constructor
VBFHZZllbbElectronSelector::VBFHZZllbbElectronSelector(const edm::ParameterSet& pset) {

  sourceLabel_   = pset.getParameter<edm::InputTag>("sourceCollection");
  sourceIDLabel_ = pset.getParameter<edm::InputTag>("sourceIDLabel");
  sourcePtMinCut_  = pset.getParameter<double>("sourcePtMin");
  sourceEtaMaxCut_ = pset.getParameter<double>("sourceEtaMax");

  string iName = "selectedElectron";
  produces<reco::GsfElectronCollection>(); 

  counterelectron=0;	
  counterelectronbefore=0;

}


VBFHZZllbbElectronSelector::~VBFHZZllbbElectronSelector() {

  cout << "Size of electrons selected before EleID" << counterelectronbefore 
       << " and after EleID" << counterelectron << endl;
}

void VBFHZZllbbElectronSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  auto_ptr<reco::GsfElectronCollection> Gelectron( new reco::GsfElectronCollection );

  // Get all pixel match GSF electron candidates
  edm::Handle<edm::View<GsfElectron> > electrons;
  iEvent.getByLabel(sourceLabel_, electrons);

  // Get the valuemap of electron from electron ID
  edm::Handle< edm::ValueMap<float> >  eIDValueMap;
  iEvent.getByLabel( sourceIDLabel_ , eIDValueMap);
  const edm::ValueMap<float> & eIdmapCuts = * eIDValueMap ;
  
  // Loop over GsfElectrons
  for (unsigned int i = 0; i < electrons->size(); ++i) {	  
      Ref<edm::View<reco::GsfElectron> > electronRef(electrons,i);

      if ( electronRef->pt() < sourcePtMinCut_ ) continue;
      if ( fabs( electronRef->eta() ) > sourceEtaMaxCut_ ) continue;

      ++counterelectronbefore;
      if ( eIdmapCuts[electronRef] > 0 ) {
           Gelectron->push_back( *electronRef );
           ++counterelectron;
      }
   } 

  const string iName = "";
  iEvent.put( Gelectron, iName );

}


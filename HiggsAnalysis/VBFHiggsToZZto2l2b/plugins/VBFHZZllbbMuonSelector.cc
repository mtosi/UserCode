/**\class VBFHZZllbbMuonSelector.cc
 *
 * Original Author:  Dominique Fortin
 * Modified by       Nicola De Filippis 
 */

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbMuonSelector.h"

#include "FWCore/Framework/interface/ESHandle.h"

// Muons:
#include <DataFormats/MuonReco/interface/Muon.h>
#include <DataFormats/MuonReco/interface/MuonFwd.h>
#include <DataFormats/TrackReco/interface/Track.h>

// Candidate handling
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include "DataFormats/Common/interface/AssociationVector.h"

// system include files
#include <Math/VectorUtil.h>
#include <memory>
#include <vector>

using namespace edm;
using namespace std;
using namespace reco;

// constructor
VBFHZZllbbMuonSelector::VBFHZZllbbMuonSelector(const edm::ParameterSet& pset) {

  sourceLabel_ = pset.getParameter<edm::InputTag>("sourceLabel");
  sourceMinPtBarrelCut_ = pset.getParameter<double>("sourceMinPtBarrelCut");
  sourceMinPtEndcapCut_ = pset.getParameter<double>("sourceMinPtEndcapCut");
  sourceMinPEndcapCut_  = pset.getParameter<double>("sourceMinPEndcapCut");
  sourceMaxEtaCut_      = pset.getParameter<double>("sourceMaxEtaCut");

  string alias;
  string iName = "selectedMuon";
  produces<reco::MuonCollection>(); 

}


VBFHZZllbbMuonSelector::~VBFHZZllbbMuonSelector() {
 
}


//
// member functions
//
void VBFHZZllbbMuonSelector::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // muons
  auto_ptr<reco::MuonCollection> Gmuon( new reco::MuonCollection );

  edm::Handle<edm::View<Muon> > muons;
  iEvent.getByLabel(sourceLabel_, muons);

  edm::View<reco::Muon>::const_iterator muon_itr = muons->begin();    

  // Loop over muons
  for ( ; muon_itr != muons->end(); ++muon_itr ) {
    bool isEndcap = false;
    bool isBarrel = false;

    if ( fabs( muon_itr->eta() ) > sourceMaxEtaCut_ ) continue;
    if ( fabs( muon_itr->eta() ) < 1.1 ) {
      isBarrel = true;
    } else {
      isEndcap = true;
    }
    // Other criteria here: 

    if ( isEndcap && 
	 muon_itr->pt() > sourceMinPtEndcapCut_ && 
	 muon_itr->p()  > sourceMinPEndcapCut_     ) Gmuon->push_back( *muon_itr );
    if ( isBarrel && 
	 muon_itr->pt() > sourceMinPtBarrelCut_ ) Gmuon->push_back( *muon_itr );   
  }

  
  const string iName = "";
  iEvent.put( Gmuon, iName );

}


/**\class VBFHZZllbbElectronIsolationProducer
 *
 * Compute isolation for cones around electron candidates
 */

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbElectronIsolationProducer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Candidate handling
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"

// Electrons
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

// Muons
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"


// Tracker tracks
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/CandidateTkIsolation.h"

//
#include <memory>
#include <vector>

using namespace edm;
using namespace std;
using namespace reco;

// constructor
VBFHZZllbbElectronIsolationProducer::VBFHZZllbbElectronIsolationProducer(const edm::ParameterSet& pset) {
  
  electronLabel_     = pset.getParameter<edm::InputTag>("ElectronsLabel"    );
  electronVetoLabel_ = pset.getParameter<edm::InputTag>("ElectronsVetoLabel");
  muonsLabel_        = pset.getParameter<edm::InputTag>("MuonsLabel"        );
  tracksLabel_       = pset.getParameter<edm::InputTag>("TracksLabel"       );
  isolationConeCut_ = pset.getParameter<double>("isolationConeCut"    );
  isolationVetoCut_ = pset.getParameter<double>("isolationConeVetoCut");
  isolationCut_     = pset.getParameter<double>("isolationCut"        );
  //ptMin          = pset.getParameter<double>("minTrackerTrackPt");
  //maxDz          = pset.getParameter<double>("isolationDeltaZ");

  string alias;
  string iName = "ElectronIsolation";

  produces<vector<float> >( alias = iName + "SumpT"       ).setBranchAlias( alias );
  produces<vector<float> >( alias = iName + "SumpToverpT" ).setBranchAlias( alias );
  produces<vector<float> >( alias = iName + "SumpT2"      ).setBranchAlias( alias );
  produces<reco::GsfElectronCollection>();

}


// destructor
VBFHZZllbbElectronIsolationProducer::~VBFHZZllbbElectronIsolationProducer() {

}


void VBFHZZllbbElectronIsolationProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  auto_ptr<vector<float> > electronIsoSumpT(    new vector<float> );
  auto_ptr<vector<float> > electronIsoSumpT_pT( new vector<float> );
  auto_ptr<vector<float> > electronIsoSumpT2(   new vector<float> );

  auto_ptr<reco::GsfElectronCollection> isolatedElectron( new reco::GsfElectronCollection );

 
  // Get the tracker tracks
  Handle<edm::View<reco::Track> > tracks;
  iEvent.getByLabel(tracksLabel_.label(), tracks);
  const edm::View<reco::Track>* trackCollection = tracks.product () ;

  //Get electron candidates  
  Handle<edm::View<reco::GsfElectron> > electronCandidates;
  iEvent.getByLabel(electronLabel_, electronCandidates);
  cout << "Size of input collection of electrons" << electronCandidates->size() << endl;

  //Get electron veto candidates
  Handle<edm::View<reco::GsfElectron> > electronVetoCandidates;
  iEvent.getByLabel(electronVetoLabel_, electronVetoCandidates);
  const edm::View<reco::GsfElectron>* Resolved_Collection = electronVetoCandidates.product () ;
 
   // get muons
  Handle<edm::View<reco::Muon> > allMuons;
  iEvent.getByLabel(muonsLabel_, allMuons);
  const edm::View<reco::Muon>* muonCollection =  allMuons.product () ;

  edm::View<reco::GsfElectron>::const_iterator electronCandidates_itr = electronCandidates->begin();
  for ( ; electronCandidates_itr != electronCandidates->end(); ++electronCandidates_itr) {
    CandidateTkIsolation myTkIsolation(&(*electronCandidates_itr),trackCollection,Resolved_Collection,muonCollection,1) ;
    
    myTkIsolation.setIntRadius (isolationVetoCut_) ;
    myTkIsolation.setExtRadius (isolationConeCut_) ;
    
    // isolation cut
    if ( myTkIsolation.getPtTracksCorr()[0]/electronCandidates_itr->pt() < isolationCut_ ){
      electronIsoSumpT->push_back(myTkIsolation.getPtTracksCorr()[0]);
      electronIsoSumpT_pT->push_back(myTkIsolation.getPtTracksCorr()[0]/electronCandidates_itr->pt());
      electronIsoSumpT2->push_back(myTkIsolation.getPtTracksCorr()[0]*myTkIsolation.getPtTracksCorr()[0]);

      isolatedElectron->push_back(*electronCandidates_itr);
    }
  }

  sort(electronIsoSumpT->begin(),electronIsoSumpT->end());
  sort(electronIsoSumpT_pT->begin(),electronIsoSumpT_pT->end());        
  sort(electronIsoSumpT2->begin(),electronIsoSumpT2->end()); 
  
  const string & isoName = "ElectronIsolation";
  iEvent.put( electronIsoSumpT,    isoName + "SumpT"       );
  iEvent.put( electronIsoSumpT_pT, isoName + "SumpToverpT" );
  iEvent.put( electronIsoSumpT2,   isoName + "SumpT2"      );
  const string iName = "";
  iEvent.put( isolatedElectron, iName );
}



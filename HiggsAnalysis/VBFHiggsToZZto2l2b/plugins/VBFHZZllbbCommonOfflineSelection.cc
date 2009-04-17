/* \class VBFHZZllbbCommonOfflineSelection
 *
 *
 * Tight isolation: electron and muons
 */


// system include files
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbCommonOfflineSelection.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/VBFHZZllbbElectronAssociationMap.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/VBFHZZllbbMuonAssociationMap.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <memory>

// namespaces
using namespace edm;
using namespace std;
using namespace reco;

struct SortCandByDecreasingPt {
  bool operator()( const Candidate &c1, const Candidate &c2) const {
    return c1.pt() > c2.pt();
  }
};

bool cmp( float a, float b ) {
  return a > b;
}
 
// Constructor
VBFHZZllbbCommonOfflineSelection::VBFHZZllbbCommonOfflineSelection(const edm::ParameterSet& pset) {

  // Decay Channel
  decaychannel = pset.getParameter<std::string>("decaychannel");
  useBestCandidate = pset.getParameter<bool> ("useBestCandidate" );
  bestCandidatesLeptonsLabel_ = pset.getParameter<edm::InputTag>("bestCandidatesLeptons");

  // tight isolation
  leptonLabel_    = pset.getParameter<edm::InputTag>("leptonLabel"   );
  leptonMapLabel_ = pset.getParameter<edm::InputTag>("leptonMapLabel");
  // tight isolation cuts
  leptonIsoVarLabel_ = pset.getParameter<edm::InputTag>  ("leptonIsoVarLabel");
  leptonIsoVarCut_   = pset.getParameter<vector<double> >("leptonIsoVarCut"  );

  // vertexing
  leptonVertexLabel_    = pset.getParameter<edm::InputTag>("leptonVertexLabel"   );
  leptonVertexMapLabel_ = pset.getParameter<edm::InputTag>("leptonVertexMapLabel");

  // vertexing cuts
  vertexVarCut_ = pset.getParameter<vector<double> >("vertexVarCut");

  std::cout << "Starting Offline selection for channel " << decaychannel << std::endl;

  string alias;
  if (decaychannel=="2e2b"){
    produces<bool> (alias = decaychannel + "OffselTightCombIsolEle").setBranchAlias( alias );
  }
  if (decaychannel=="2mu2b"){
    produces<bool> (alias = decaychannel + "OffselTightCombIsolMu").setBranchAlias( alias );
  }

  produces<bool> (alias = decaychannel + "OffselVertComb").setBranchAlias( alias );
  produces<bool> (alias = decaychannel + "Offsel"  ).setBranchAlias( alias );

}


// Destructor
VBFHZZllbbCommonOfflineSelection::~VBFHZZllbbCommonOfflineSelection() {

}


// Filter event (event preselection)
void VBFHZZllbbCommonOfflineSelection::produce(edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  bool matched=false;
  
  // flags
  auto_ptr<bool> boolTightCombIsoElectron ( new bool );
  *boolTightCombIsoElectron  = false;
  auto_ptr<bool> boolTightCombIsoMuon ( new bool );
  *boolTightCombIsoMuon   = false;
  
  auto_ptr<bool> boolVertexComb ( new bool );
  *boolVertexComb = false;
  
  // Offline selection flag
  auto_ptr<bool> boolOffSel ( new bool );
  *boolOffSel = false;  
  
  //Tight isolation for electrons
  int nTightCombIsolElectronCounter = 0;
  if (decaychannel=="2e2b"){
    
    Handle<edm::View<GsfElectron> > electronHandle;
    iEvent.getByLabel(leptonLabel_, electronHandle);
    
    Handle<VBFHZZllbbElectronAssociationMap> isoElectronMapHandle;
    iEvent.getByLabel(leptonMapLabel_, isoElectronMapHandle);
    
    
    vector<float> isoElectronVector;
    int index=0;
    edm::View<reco::GsfElectron>::const_iterator electron_itr = electronHandle->begin();
    for ( ; electron_itr != electronHandle->end(); ++electron_itr) {
      edm::Ref<edm::View<reco::GsfElectron> > electronTrackRef(electronHandle,index);
      std::cout << "Isolation value from electron map " << (*isoElectronMapHandle)[electronTrackRef] << std::endl;
      
      if (useBestCandidate) {
	Handle<CandidateCollection> bestLeptonCandidatesHandle;
	iEvent.getByLabel(bestCandidatesLeptonsLabel_, bestLeptonCandidatesHandle);
	const reco::CandidateCollection* bestLeptons = bestLeptonCandidatesHandle.product();
	matched=match(electron_itr->mass(),electron_itr->pt(),electron_itr->charge(),bestLeptons);
      }
      else {
	matched=true;
      }
      
      if (matched) isoElectronVector.push_back((*isoElectronMapHandle)[electronTrackRef]);
      index ++;
    }
    // sorting in decreasing order
    sort(isoElectronVector.begin(),isoElectronVector.end(),cmp);
    // cut on the two least isolated electrons; sum of iso;
    // if the least are tight isolated --> all are tight isolated
    
    if (isoElectronVector.size() >= 2){
      if ( (isoElectronVector.at(0)+isoElectronVector.at(1)) < leptonIsoVarCut_.at(0) ) nTightCombIsolElectronCounter++;
    }
    else {
      nTightCombIsolElectronCounter=-999;
      std::cout << "Warning: there are no 2 values for isolation variables" << std::endl;
    }     
  }
  
  // Tight isolation for muons
  int nTightCombIsolMuonCounter=0;
  if (decaychannel=="2mu2b"){    
    
    Handle<edm::View<Muon> > muonHandle;
    iEvent.getByLabel(leptonLabel_, muonHandle);
    
    Handle<VBFHZZllbbMuonAssociationMap> isoMuonMapHandle;
    iEvent.getByLabel(leptonMapLabel_, isoMuonMapHandle);
    
    vector<float> isoMuonVector;
    int index=0;
    edm::View<reco::Muon>::const_iterator muon_itr = muonHandle->begin(); 
    for ( ; muon_itr != muonHandle->end(); ++muon_itr) {
      edm::Ref<edm::View<reco::Muon> > muonTrackRef(muonHandle,index);
      std::cout << "Isolation value from muon map " << (*isoMuonMapHandle)[muonTrackRef] << std::endl;
      
      if (useBestCandidate) {
	Handle<CandidateCollection> bestLeptonCandidatesHandle;
	iEvent.getByLabel(bestCandidatesLeptonsLabel_, bestLeptonCandidatesHandle);
	const reco::CandidateCollection* bestLeptons = bestLeptonCandidatesHandle.product () ;
	matched=match(muon_itr->mass(),muon_itr->pt(),muon_itr->charge(),bestLeptons);
      }
      else {
	matched=true;
      }
      
      if (matched) isoMuonVector.push_back((*isoMuonMapHandle)[muonTrackRef]);
      index++;
    }
    // sorting in decreasing order
    sort(isoMuonVector.begin(),isoMuonVector.end(),cmp);
    // cut on the two least isolated muons; sum of iso; if the least are tight isolated --> all are tight isolated
    if (isoMuonVector.size() >= 2){
      if ( (isoMuonVector.at(0)+isoMuonVector.at(1)) < leptonIsoVarCut_.at(0) ) nTightCombIsolMuonCounter++;
    }
    else {
      nTightCombIsolMuonCounter=-999;
      std::cout << "Warning: there are no 2 values for isolation variables" << std::endl;
    }     
  }

  // Vertexing      
  vector<float> vertexLeptonVector;
  int nVertexCombCounter=0;
  
  if (decaychannel=="2mu2b"){
    Handle<edm::View<Muon> > muonVertexHandle;
    iEvent.getByLabel(leptonVertexLabel_, muonVertexHandle);
    
    Handle<VBFHZZllbbMuonAssociationMap> vertexMuonMapHandle;
    iEvent.getByLabel(leptonVertexMapLabel_, vertexMuonMapHandle);
    
    int index=0;
    edm::View<reco::Muon> ::const_iterator muonVertex_itr = muonVertexHandle->begin(); 
    for ( ; muonVertex_itr != muonVertexHandle->end(); ++muonVertex_itr) {
      edm::Ref<edm::View<reco::Muon> > muonTrackRef(muonVertexHandle,index);
      std::cout << "Vertexing value from muon map " << (*vertexMuonMapHandle)[muonTrackRef] << std::endl;

      if (useBestCandidate) {
	Handle<CandidateCollection> bestLeptonCandidatesHandle;
	iEvent.getByLabel(bestCandidatesLeptonsLabel_, bestLeptonCandidatesHandle);
	const reco::CandidateCollection* bestLeptons = bestLeptonCandidatesHandle.product();
	matched=match(muonVertex_itr->mass(),muonVertex_itr->pt(),muonVertex_itr->charge(),bestLeptons);
      }
      else {
	matched=true;
      }

      if (matched) vertexLeptonVector.push_back(fabs((*vertexMuonMapHandle)[muonTrackRef]));
      index++;
    }    
  }
  
  if (decaychannel=="2e2b"){
    Handle<edm::View<GsfElectron> > electronVertexHandle;
    iEvent.getByLabel(leptonVertexLabel_, electronVertexHandle);
    
    Handle<VBFHZZllbbElectronAssociationMap> vertexElectronMapHandle;
    iEvent.getByLabel(leptonVertexMapLabel_, vertexElectronMapHandle);

    int index=0;
    edm::View<reco::GsfElectron>::const_iterator electronVertex_itr = electronVertexHandle->begin(); 
    for ( ; electronVertex_itr != electronVertexHandle->end(); ++electronVertex_itr) {
      edm::Ref<edm::View<reco::GsfElectron> > electronTrackRef(electronVertexHandle,index);
      std::cout << "Vertexing value from electron map " << (*vertexElectronMapHandle)[electronTrackRef] << std::endl;
    
      if (useBestCandidate) {
	Handle<CandidateCollection> bestLeptonCandidatesHandle;
	iEvent.getByLabel(bestCandidatesLeptonsLabel_, bestLeptonCandidatesHandle);
	const reco::CandidateCollection* bestLeptons = bestLeptonCandidatesHandle.product();
	matched=match(electronVertex_itr->mass(),electronVertex_itr->pt(),electronVertex_itr->charge(),bestLeptons);
      }
      else {
	matched=true;
      }

      if (matched) vertexLeptonVector.push_back(fabs((*vertexElectronMapHandle)[electronTrackRef]));
      index++;
    }
  }
  
  // sorting in decreasing order
  sort(vertexLeptonVector.begin(),vertexLeptonVector.end(),cmp);
  // cut on the two least IP leptons;  
  if (vertexLeptonVector.size() >= 2){
    if ( vertexLeptonVector.at(0) < vertexVarCut_.at(0) && 
	 vertexLeptonVector.at(1) < vertexVarCut_.at(1) ) nVertexCombCounter++;
  }
  else {
    nVertexCombCounter=-999;
    std::cout << "Warning: there are no 2 values for vertex significance" << std::endl;
  }
  
  
  //  2e2b, 2mu2b channels
  if (decaychannel=="2e2b"){
    if (nTightCombIsolElectronCounter >= 1) {
      *boolTightCombIsoElectron = true;
      if ( nVertexCombCounter>=1 ) 
	*boolOffSel=true;
    }
    iEvent.put(boolTightCombIsoElectron, decaychannel + "OffSelTightCombIsolEle");

  }
  if (decaychannel=="2mu2b"){
    if (nTightCombIsolMuonCounter >= 1) {
      *boolTightCombIsoMuon  = true;
      if ( nVertexCombCounter>=1 ) 
	*boolOffSel=true;
    }
    iEvent.put(boolTightCombIsoMuon, decaychannel + "OffSelTightCombIsolMu" );
  }
  
  if (nVertexCombCounter>= 1) *boolVertexComb=true; 
  iEvent.put(boolVertexComb, decaychannel + "OffSelVertComb");
  iEvent.put(boolOffSel, decaychannel + "OffSel");
  
}

void VBFHZZllbbCommonOfflineSelection::beginJob(const edm::EventSetup& iSetup) {
}

void VBFHZZllbbCommonOfflineSelection::endJob() {
  std::cout << "Created offline selection variables" << std::endl;
}


bool VBFHZZllbbCommonOfflineSelection::match(double mass, 
					     double pt, 
					     int charge, 
					     const reco::CandidateCollection *c1Coll){

  bool found=false;

  for( CandidateCollection::const_iterator pp = c1Coll->begin();pp != c1Coll->end(); ++pp ) {

    if ((abs(pp->p4().mass()-mass)  <0.001 ) &&
        (abs(pp->p4().pt()  -pt)    <0.001 ) &&
        (abs(pp->charge()   -charge)<0.001 )  ){
      found=true;
      std::cout << "Found lepton in the best leptons collection" << std::endl;
    }
  }
  return found;
}


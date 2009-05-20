// -*- C++ -*-
//
// Package:    SimpleNtple
// Class:      SimpleNtple
// 
/**\class SimpleNtple SimpleNtple.cc Analysis/SimpleNtple/src/SimpleNtple.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/RefToBase.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/SimpleNtple.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

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

// jet correction
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/CorJetWithBTagDiscr.h"

// utilities
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/VBFHZZllbbUtils.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/ProcessIndex.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace vbfhzz2l2b;


enum { FASTSIM = 0,
       FULLSIM = 1
}; // to be added to VBFHZZllbbUtils

SimpleNtple::SimpleNtple(const edm::ParameterSet& iConfig) :
  whichSim_ ( iConfig.getParameter<int> ( "whichSim"  ) ), // 0:FastSim, 1:FullSim
  trackLabel_    ( iConfig.getParameter<edm::InputTag> ( "trackLabel"    ) ),
  muonLabel_     ( iConfig.getParameter<edm::InputTag> ( "muonLabel"     ) ),
  electronLabel_ ( iConfig.getParameter<edm::InputTag> ( "electronLabel" ) ),
  metLabel_      ( iConfig.getParameter<edm::InputTag> ( "metLabel"      ) ),
  tagJetLabel_   ( iConfig.getParameter<edm::InputTag> ( "tagJetLabel"   ) ),
  corIC5CaloJetsWithBTagLabel_ ( iConfig.getParameter<std::string> ( "corIC5CaloJetsWithBTagLabel" ) ),
  corIC5PFJetsWithBTagFlag_    ( iConfig.getParameter<bool>        ( "corIC5PFJetsWithBTagFlag"    ) ),
  genParticleLabel_ ( iConfig.getParameter<edm::InputTag> ( "genParticleLabel" ) ),
  genJetLabel_      ( iConfig.getParameter<edm::InputTag> ( "genJetLabel"      ) ),
  genMetLabel_      ( iConfig.getParameter<edm::InputTag> ( "genMetLabel"      ) ) {

  if ( corIC5PFJetsWithBTagFlag_ ) 
    corIC5PFJetsWithBTagLabel_ = iConfig.getParameter<std::string> ( "corIC5PFJetsWithBTagLabel" );
 
  //now do what ever initialization is needed
  edm::Service<TFileService> fs ;
  mytree_  = fs->make <TTree>("VBFSimpleTree","VBFSimpleTree"); 
  
  std::cout << "[SimpleNtple::SimpleNtple] DONE" << std::endl;

}


// --------------------------------------------------------------------


SimpleNtple::~SimpleNtple()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

  delete tagjetInvMass_;    // depends on the tag jet definition
  delete tagjetDeltaEta_;
  delete tagjetZeppenfeld_;
  delete zjetInvMass_;    // depends on the z jet definition => btagger?
  delete zjetDeltaEta_;
  delete zjetZeppenfeld_;

  delete eleP4_ ;
  delete eleVtxP3_;
  delete eleEt_;
  delete elePt_;
  delete eleIsoSumPt_;
  delete eleIsoNtrack_;
  delete eleD0_;
  delete eleDxy_;
  delete eleDxyError_;
  delete eleID_;

  delete muP4_ ;
  delete muVtxP3_;
  delete muEt_;
  delete muPt_;
  delete muIsoSumPt_;
  delete muIsoNtrack_;
  delete muD0_;
  delete muDxy_;
  delete muDxyError_;
  delete muID_;

  delete tagjetP4_;
  delete tagjetVtxP3_;
  delete tagjetEmFrac_;
  delete tagjetChFrac_;
  delete tagjetCorEt_;
  delete tagjetCorPt_;

  delete btagjetP4_;
  delete btagjetVtxP3_; 
  delete btagjetEmFrac_;
  delete btagjetChFrac_;
  delete btagjetCorEt_;
  delete btagjetCorPt_;
  delete btagjetCompoSVbTagDiscr_;
  delete btagjetHighEFFbTagDiscr_;
  delete btagjetHighPURbTagDiscr_;

  delete metP4_ ;
  delete metSig_ ;

  delete trackP4_ ;
  delete genparticleP4_ ;
  delete genjetP4_;
  delete genmetP4_;

  std::cout << "[SimpleNtple::~SimpleNtple]" << std::endl;
  
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
SimpleNtple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  InitObjs();
  std::cout << "SimpleNtple::analyze] InitObjs DONE" << std::endl;
  FillEvent                  (iEvent, iSetup);
  FillcorIC5CaloJetsWithBTag (iEvent, iSetup);
  FillMuon                   (iEvent, iSetup);
  FillElectron               (iEvent, iSetup);
  FillZhad                   (iEvent, iSetup);
  FillZlep                   (iEvent, iSetup);
  FillMet                    (iEvent, iSetup);
  FillTagJet                 (iEvent, iSetup); // not implemented yet
  if ( corIC5PFJetsWithBTagFlag_ )
    FillcorIC5PFJetsWithBTag (iEvent, iSetup);   
  if ( whichSim_ == FULLSIM )
    FillTrack (iEvent, iSetup);
  FillGenParticle  (iEvent, iSetup); // got an error message in execution
  FillGenJet       (iEvent, iSetup);
  FillGenMet       (iEvent, iSetup);
  
  mytree_->Fill();

  std::cout << "[SimpleNtple::analyze] DONE" << std::endl;
  
}


// --------------------------------------------------------------------

void SimpleNtple::InitObjs() {

  // event obj
  evtID_    = 0;
  evtRun_   = 0;
  evtEvent_ = 0;
  tagjetInvMass_    -> clear ();    // depends on the tag jet definition
  tagjetDeltaEta_   -> clear ();
  tagjetZeppenfeld_ -> clear ();
  zjetInvMass_      -> clear ();    // depends on the z jet definition => btagger?
  zjetDeltaEta_     -> clear ();
  zjetZeppenfeld_   -> clear ();
  //electrons;
  eleN_ = 0;
  eleP4_        -> Clear ();
  eleEt_        -> clear ();
  elePt_        -> clear ();
  eleIsoSumPt_  -> clear ();
  eleIsoNtrack_ -> clear ();
  eleD0_        -> clear ();
  eleDxy_       -> clear ();
  eleDxyError_  -> clear ();
  eleID_        -> clear ();
  eleVtxP3_     -> Clear ();
  //muons
  muN_ = 0;
  muP4_        -> Clear ();
  muEt_        -> clear ();
  muPt_        -> clear ();
  muIsoSumPt_  -> clear ();
  muIsoNtrack_ -> clear ();
  muD0_        -> clear ();
  muDxy_       -> clear ();
  muDxyError_  -> clear ();
  muID_        -> clear ();
  muVtxP3_     -> Clear ();
  // tag jets
  tagjetN_      = 0;
  tagjetNtrack_ = 0;
  tagjetP4_     -> Clear ();
  tagjetEmFrac_ -> clear ();
  tagjetChFrac_ -> clear ();
  tagjetCorEt_  -> clear ();
  tagjetCorPt_  -> clear ();
  tagjetEmFrac_ -> clear ();
  tagjetVtxP3_  -> Clear ();
  // other jets with b tag
  btagjetN_      = 0;
  btagjetNtrack_ = 0;
  btagjetP4_               -> Clear ();
  btagjetEmFrac_           -> clear ();
  btagjetChFrac_           -> clear ();
  btagjetCorEt_            -> clear ();
  btagjetCorPt_            -> clear ();
  btagjetCompoSVbTagDiscr_ -> clear ();
  btagjetHighEFFbTagDiscr_ -> clear ();
  btagjetHighPURbTagDiscr_ -> clear ();
  btagjetVtxP3_            -> Clear ();
 
  // met
  metP4_  -> Clear ();
  metSig_ -> clear ();

  trackP4_       -> Clear () ;
  genparticleP4_ -> Clear () ;
  genjetP4_	 -> Clear () ;
  genmetP4_	 -> Clear () ;
  
}

// --------------------------------------------------------------------
void SimpleNtple::FillEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  //  std::cout << "SimpleNtple::FillEvent" << std::endl;

  if ( whichSim_ == FULLSIM ) {
    edm::Handle<edm::HepMCProduct> evtMC;
    try {
      iEvent.getByLabel("source", evtMC); }
    catch(...) {
      std::cerr << "[SimpleNtple::FillKindEvent] defined as FullSim, but HepMCProduct::source not found" << std::endl; }
    const HepMC::GenEvent * mcEv = evtMC->GetEvent();
    evtID_ = mcEv->signal_process_id();
  }
  else if ( whichSim_ == FASTSIM ) {
    edm::Handle<int> genProcessID;
    try {
      iEvent.getByLabel( "genEventProcID", genProcessID ); }
    catch(...) {
      std::cerr << "[SimpleNtple::FillKindEvent] defined as FastSim, but genEventProcID not found" << std::endl; }

    evtID_ = *genProcessID;
  }
  else {
    std::cout << "--> WARNING: simulation not specificied!!" << std::endl;
  }


}

// --------------------------------------------------------------------
void SimpleNtple::FillZhad(const edm::Event& iEvent, const edm::EventSetup& iSetup){
}

// --------------------------------------------------------------------
void SimpleNtple::FillZlep(const edm::Event& iEvent, const edm::EventSetup& iSetup){
}

// --------------------------------------------------------------------
void SimpleNtple::FillcorIC5CaloJetsWithBTag(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  std::cout << "[SimpleNtple::FillcorIC5CaloJetsWithBTag]" << std::endl;

  edm::Handle<vbfhzz2l2b::CorJetWithBTagDiscrCollection> corJetWithBTagHandle;
  iEvent.getByLabel(corIC5CaloJetsWithBTagLabel_,"corJetWithBTagDiscr",corJetWithBTagHandle);
  //  std::cout << "corJetWithBTagHandle->size(): " << corJetWithBTagHandle->size() << std::endl;

  TClonesArray &jetP4 = *btagjetP4_;
  int jetIndex = 0;
  std::vector<reco::JetBaseRef> jets = vbfhzz2l2b::CorJetBTagDiscrAssociation::allJets(*corJetWithBTagHandle);  
  for ( std::vector<reco::JetBaseRef>::const_iterator jet = jets.begin(); 
	jet != jets.end(); ++jet, jetIndex++ ) {
    //    std::cout << "jetIndex: " << jetIndex << std::endl;
    
    std::vector<double> discrVec = (*corJetWithBTagHandle)[*jet].discrVec_;
//    std::cout << "corJetWithBTag highEffDiscr: "       << discrVec[vbfhzz2l2b::HIGHEFF]    << std::endl;
//    std::cout << "corJetWithBTag highPurDiscr: "       << discrVec[vbfhzz2l2b::HIGHPUR]    << std::endl;
//    std::cout << "corJetWithBTag combSecVtxDiscr: "    << discrVec[vbfhzz2l2b::COMBSECVTX] << std::endl;
    double corrEt = (*corJetWithBTagHandle)[*jet].corEt_;
    double uncorrEt = (*jet)->et();
    double uncorrPt = (*jet)->pt();
    double corPt    = (corrEt/uncorrEt)*uncorrPt;
    double emFrac = (dynamic_cast<const reco::CaloJet*>(&**jet))->emEnergyFraction();

    vbfhzz2l2b::setMomentum (myvector_, (*jet)->p4());
    new (jetP4[jetIndex]) TLorentzVector (myvector_);
    //    btagjetNtrack_ = ;
    btagjetEmFrac_  -> push_back (emFrac);
    //    btagjetChFrac_ -> push_back ();
    btagjetCorEt_            -> push_back (corrEt);
    btagjetCorPt_	     -> push_back (corPt);
    btagjetHighPURbTagDiscr_ -> push_back (discrVec[vbfhzz2l2b::HIGHEFF]   );
    btagjetHighEFFbTagDiscr_ -> push_back (discrVec[vbfhzz2l2b::HIGHPUR]   );
    btagjetCompoSVbTagDiscr_ -> push_back (discrVec[vbfhzz2l2b::COMBSECVTX]);
    //  btagjetVtxP3_;  

  }
  btagjetN_ = jets.size();
}

// --------------------------------------------------------------------
void SimpleNtple::FillcorIC5PFJetsWithBTag(const edm::Event& iEvent, const edm::EventSetup& iSetup){
}

// --------------------------------------------------------------------
void SimpleNtple::FillMuon(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  //  std::cout << "SimpleNtple::FillMu" << std::endl;
  edm::Handle<reco::MuonCollection> muonHandle;
  iEvent.getByLabel (muonLabel_,muonHandle);


  TClonesArray &muonP4 = *muP4_;
  TClonesArray &muonVtxP3 = *muVtxP3_;

  int muonIndex = 0;
  for ( reco::MuonCollection::const_iterator muon_itr = muonHandle->begin();
	muon_itr != muonHandle->end(); ++muon_itr, muonIndex++ ) {
    vbfhzz2l2b::setMomentum (myvector_, muon_itr->p4());
    //    setVertex(myvertex_,muon_itr->Vertex());
    new (muonP4[muonIndex]) TLorentzVector (myvector_);
    new (muonVtxP3[muonIndex]) TVector3 (myvertex_);

    //    muEt_        -> push_back ();
    muPt_        -> push_back (muon_itr->pt());
    //    muIsoSumPt_  -> push_back ();
    //    muIsoNtrack_ -> push_back ();
    //    muD0_        -> push_back ();
    //    muDxy_       -> push_back ();
    //    muDxyError_  -> push_back ();
    //    muID_        -> push_back ();

  }
  muN_ = muonHandle->size(); 

}

// --------------------------------------------------------------------
void SimpleNtple::FillElectron(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::cout << "[SimpleNtple::FillElectron]" << std::endl;

  edm::Handle<reco::PixelMatchGsfElectronCollection> electronHandle ;
  iEvent.getByLabel (electronLabel_,electronHandle) ;

  TClonesArray &electronP4 = *eleP4_;
  TClonesArray &electronVtxP3 = *eleVtxP3_;
  int electronIndex = 0;
  for ( reco::PixelMatchGsfElectronCollection::const_iterator electron_itr = electronHandle->begin();
	electron_itr != electronHandle->end(); ++electron_itr, electronIndex++ ) {
    vbfhzz2l2b::setMomentum (myvector_, electron_itr->p4());
    new (electronP4[electronIndex]) TLorentzVector (myvector_);
  }

  eleN_ = electronHandle->size(); 
  
}

// --------------------------------------------------------------------
void SimpleNtple::FillMet(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  //  std::cout << "SimpleNtple::FillMet" << std::endl;
  edm::Handle<reco::CaloMETCollection> metCollectionHandle;
  iEvent.getByLabel (metLabel_ , metCollectionHandle);
  const CaloMETCollection *calometcol = metCollectionHandle.product();
  const CaloMET *calomet = &(calometcol->front());

  TClonesArray &MET = *metP4_;
  vbfhzz2l2b::setMomentum (myvector_, calomet->p4());
  new (MET[0]) TLorentzVector (myvector_);
  
}


// --------------------------------------------------------------------
void SimpleNtple::FillTagJet(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  //  std::cout << "SimpleNtple::FillTagJet" << std::endl;

  // edm::Handle<reco::RecoChargedCandidateCollection> tagJetHandle;
  // edm::Handle<reco::CandidateCollection> jetTagsHandle ;
  edm::Handle<reco::CaloJetCollection> tagJetHandle;
  iEvent.getByLabel (tagJetLabel_, tagJetHandle) ;

  typedef reco::CaloJetCollection::const_iterator tagJetItr;  

  // looking for the highest invariant mass jets pair
  std::pair<tagJetItr,tagJetItr> maxInvMassPair = 
    vbfhzz2l2b::findPair_maxInvMass_ptMinCut<tagJetItr>(tagJetHandle->begin(), tagJetHandle->end(),
							20., 15.);

  double invMass    = -99.;
  double deltaEta   = -99.;
  double zeppenfeld = -999.;
  if (maxInvMassPair.first != maxInvMassPair.second) {
    invMass    =  (     (maxInvMassPair.first)->p4() + ((maxInvMassPair.second)->p4()) ).M();
    deltaEta   =  fabs( (maxInvMassPair.first)->eta() - (maxInvMassPair.second)->eta() );
    zeppenfeld =  (     (maxInvMassPair.first)->pz() * (maxInvMassPair.second)->pz() );
  }
  tagjetInvMass_->push_back(invMass);
  tagjetDeltaEta_->push_back(deltaEta);
  tagjetZeppenfeld_->push_back(zeppenfeld);

  // looking for the highest delta eta jets pair
  std::pair<tagJetItr,tagJetItr> maxDeltaEtaPair = 
    vbfhzz2l2b::findPair_maxDeltaEta_ptMinCut<tagJetItr>(tagJetHandle->begin(), tagJetHandle->end(),
							 20., 15.);
  invMass    = -99.;
  deltaEta   = -99.;
  zeppenfeld = -999.;
  if(maxDeltaEtaPair.first != maxDeltaEtaPair.second) {
    invMass    =  (     (maxDeltaEtaPair.first)->p4() + ((maxDeltaEtaPair.second)->p4()) ).M();
    deltaEta   =  fabs( (maxDeltaEtaPair.first)->eta() - (maxDeltaEtaPair.second)->eta() );
    zeppenfeld =  (     (maxDeltaEtaPair.first)->pz() * (maxDeltaEtaPair.second)->pz() );
  }
  tagjetInvMass_->push_back(invMass);
  tagjetDeltaEta_->push_back(deltaEta);
  tagjetZeppenfeld_->push_back(zeppenfeld);


  TClonesArray &jetTag = *tagjetP4_;
  vbfhzz2l2b::setMomentum (myvector_, (maxInvMassPair.first)->p4());
  new (jetTag[0]) TLorentzVector (myvector_);
  vbfhzz2l2b::setMomentum (myvector_, (maxInvMassPair.second)->p4());
  new (jetTag[1]) TLorentzVector (myvector_);
  vbfhzz2l2b::setMomentum (myvector_, (maxDeltaEtaPair.first)->p4());
  new (jetTag[2]) TLorentzVector (myvector_);
  vbfhzz2l2b::setMomentum (myvector_, (maxDeltaEtaPair.second)->p4());
  new (jetTag[3]) TLorentzVector (myvector_);

}

// --------------------------------------------------------------------
void  SimpleNtple::FillTrack(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::cout << "[SimpleNtple::FillTrack]" << std::endl;
  edm::Handle<reco::TrackCollection> trackHandle ;
  iEvent.getByLabel (trackLabel_, trackHandle) ;

  TClonesArray &track = *trackP4_;
  int trackIndex = 0;
  for (reco::TrackCollection::const_iterator track_itr = trackHandle->begin (); 
       track_itr != trackHandle->end (); ++track_itr, trackIndex++ ) { 

    math::XYZVector mom = track_itr->innerMomentum () ; 
    myvector_.SetPx (mom.x ()) ;
    myvector_.SetPy (mom.y ()) ;
    myvector_.SetPz (mom.z ()) ;
    
    new (track[trackIndex]) TLorentzVector (myvector_);
  }
}


// --------------------------------------------------------------------
void SimpleNtple::FillGenParticle(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::cout << "[SimpleNtple::FillGenParticle]" << std::endl;
  edm::Handle<reco::GenParticleCollection> genParticleHandle; 
  iEvent.getByLabel (genParticleLabel_,genParticleHandle);
 
  TClonesArray &genParticle = *genparticleP4_;
  int particleIndex = 0;
  if ( evtID_ == HWWFusion_ || evtID_ == HZZFusion_ ){ //---- only if VBF
    for (reco::GenParticleCollection::const_iterator particle_itr = genParticleHandle->begin(); 
	 particle_itr != genParticleHandle->end(); ++particle_itr, particleIndex++ ) {
      
      int pdg = particle_itr->pdgId();
      int status = particle_itr->status();
      int mother1 = 0;
      if ( particle_itr->numberOfMothers() > 0 ) mother1 = particle_itr->mother(0)->pdgId();
      int mother2 = 0;
      if ( particle_itr->numberOfMothers() > 1 ) mother2 = particle_itr->mother(1)->pdgId();
      int daughter1 = 0;
      if ( particle_itr->numberOfDaughters() > 0 ) daughter1 = particle_itr->daughter(0)->pdgId();
      int daughter2 = 0;
      if ( particle_itr->numberOfDaughters() > 1 ) daughter2 = particle_itr->daughter(1)->pdgId();
      double px = particle_itr->px();
      double py = particle_itr->py();
      double pz = particle_itr->pz();
      double e  = particle_itr->energy();
      double vx = 0.;
      double vy = 0.;
      double vz = 0.;
      double time = 0.;
  
      new (genParticle[particleIndex]) TParticle (pdg, status, 
						  mother1, mother2,
						  daughter1, daughter2,
						  px, py, pz, e,
						  vx, vy, vz, time);
    }
  }
}




// --------------------------------------------------------------------


void SimpleNtple::FillGenJet(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::cout << "[SimpleNtple::FillGenJet]" << std::endl;
  edm::Handle< reco::GenJetCollection > genJetHandle ;
  iEvent.getByLabel( genJetLabel_, genJetHandle ) ;
  
  TClonesArray &genJet = *genjetP4_;
  int genjetIndex = 0;
  for (reco::GenJetCollection::const_iterator genjet_itr = genJetHandle->begin (); 
       genjet_itr != genJetHandle->end (); ++genjet_itr, genjetIndex++ ) { 
  
    myvector_.SetPx ( genjet_itr->px() );
    myvector_.SetPy ( genjet_itr->py() );
    myvector_.SetPz ( genjet_itr->pz() );
    myvector_.SetE  ( genjet_itr->emEnergy() + genjet_itr->hadEnergy() );
    new (genJet[genjetIndex]) TLorentzVector (myvector_);

  }

}



// --------------------------------------------------------------------
void SimpleNtple::FillGenMet(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  std::cout << "[SimpleNtple::FillGenMet]" << std::endl;
  edm::Handle< reco::GenMETCollection > genMetHandle ;
  iEvent.getByLabel( genMetLabel_, genMetHandle ) ;
  
  TClonesArray &genMets = *genmetP4_;
  int metIndex = 0;
  for (reco::GenMETCollection::const_iterator met_itr = genMetHandle->begin (); 
       met_itr != genMetHandle->end (); ++met_itr, metIndex++ ) { 
    myvector_.SetPx ( met_itr->px () );
    myvector_.SetPy ( met_itr->py () );
    myvector_.SetPz ( met_itr->pz () );
    myvector_.SetE ( met_itr->emEnergy () + met_itr->hadEnergy () );
    new (genMets[metIndex]) TLorentzVector (myvector_);
  }
  std::cout << "[SimpleNtple::FillGenMet] DONE" << std::endl;
}


// ------------ method called once each job just before starting event loop  ------------
void SimpleNtple::beginJob(const edm::EventSetup& iSetup)
{

  std::cout << "[SimpleNtple::beginJob]" << std::endl;
  tagjetInvMass_    = new std::vector<double>;    // depends on the tag jet definition
  tagjetDeltaEta_   = new std::vector<double>;
  tagjetZeppenfeld_ = new std::vector<double>;
  zjetInvMass_    = new std::vector<double>;    // depends on the z jet definition => btagger?
  zjetDeltaEta_   = new std::vector<double>;
  zjetZeppenfeld_ = new std::vector<double>;
  mytree_->Branch("evtID",&evtID_,"evtID_/I");
 
  // vector of the TLorentz Vectors of electron
  eleP4_    = new TClonesArray ("TLorentzVector");
  eleVtxP3_ = new TClonesArray ("TVector3");
  eleEt_        = new std::vector<double>;
  elePt_        = new std::vector<double>;
  eleD0_        = new std::vector<double>;
  eleDxy_       = new std::vector<double>;
  eleDxyError_  = new std::vector<double>;
  eleIsoSumPt_  = new std::vector<double>;
  eleIsoNtrack_ = new std::vector<double>;
  eleID_        = new std::vector<int>;
  mytree_->Branch("electronP4",     "TClonesArray",       &eleP4_,    256000,0);
  mytree_->Branch("electronVtxP3",  "TClonesArray",       &eleVtxP3_, 256000,0);
  mytree_->Branch("eleN",                                 &eleN_,"eleN_/I");
  mytree_->Branch("eleIsoSumPt",    "std::vector<double>",&eleIsoSumPt_);
  mytree_->Branch("eleIsoEleNtrack","std::vector<double>",&eleIsoNtrack_);
  mytree_->Branch("eleID",          "std::vector<int>",   &eleID_);     

  // vector of the TLorentz Vectors of muon
  muP4_    = new TClonesArray ("TLorentzVector");
  muVtxP3_ = new TClonesArray ("TVector3");
  muEt_        = new std::vector<double>;
  muPt_        = new std::vector<double>;
  muD0_        = new std::vector<double>;
  muDxy_       = new std::vector<double>;
  muDxyError_  = new std::vector<double>;
  muIsoSumPt_  = new std::vector<double>;
  muIsoNtrack_ = new std::vector<double>;
  muID_        = new std::vector<int>;
  mytree_->Branch("muonP4",     "TClonesArray",       &muP4_,    256000,0);
  mytree_->Branch("muonVtxP3",  "TClonesArray",       &muVtxP3_, 256000,0);
  mytree_->Branch("muN",                              &muN_,      "muN_/I");
  mytree_->Branch("muIsoSumPt", "std::vector<double>",&muIsoSumPt_);
  mytree_->Branch("muIsoNtrack","std::vector<double>",&muIsoNtrack_);
  mytree_->Branch("muID",       "std::vector<int>",   &muID_);

  // vector with the 2 tag TLorentzVectors
  tagjetP4_    = new TClonesArray ("TLorentzVector");
  tagjetVtxP3_ = new TClonesArray ("TVector3");
  tagjetEmFrac_ = new std::vector<double>;
  tagjetChFrac_ = new std::vector<double>;
  tagjetCorEt_  = new std::vector<double>;
  tagjetCorPt_  = new std::vector<double>;
  mytree_->Branch ("tagjetP4",    "TClonesArray",       &tagjetP4_,    256000,0);
  mytree_->Branch ("tagjetVtxP3", "TClonesArray",       &tagjetVtxP3_, 256000,0);
  mytree_->Branch ("tagjetEmFrac","std::vector<double>",&tagjetEmFrac_);
  mytree_->Branch ("tagjetChFrac","std::vector<double>",&tagjetChFrac_);
  mytree_->Branch ("tagjetCorEt", "std::vector<double>",&tagjetCorEt_);
  mytree_->Branch ("tagjetCorPt", "std::vector<double>",&tagjetCorPt_);

  // vector of the TLorentz Vectors of other jets with b tag
  btagjetP4_    = new TClonesArray ("TLorentzVector");
  btagjetVtxP3_ = new TClonesArray ("TVector3");
  btagjetEmFrac_           = new std::vector<double>;
  btagjetChFrac_           = new std::vector<double>;
  btagjetCorEt_            = new std::vector<double>;
  btagjetCorPt_            = new std::vector<double>;
  btagjetCompoSVbTagDiscr_ = new std::vector<double>;
  btagjetHighEFFbTagDiscr_ = new std::vector<double>;
  btagjetHighPURbTagDiscr_ = new std::vector<double>;
  mytree_->Branch ("btagjetP4",              "TClonesArray",       &btagjetP4_,    256000,0);
  mytree_->Branch ("btagjetVtxP3",           "TClonesArray",       &btagjetVtxP3_, 256000,0);
  mytree_->Branch ("btagjetEmFrac",          "std::vector<double>",&btagjetEmFrac_);
  mytree_->Branch ("btagjetChFrac",          "std::vector<double>",&btagjetChFrac_);
  mytree_->Branch ("btagjetCorEt",           "std::vector<double>",&btagjetCorEt_);
  mytree_->Branch ("btagjetCorPt",           "std::vector<double>",&btagjetCorPt_);
  mytree_->Branch ("btagjetcompoSVbTagDiscr","std::vector<double>",&btagjetCompoSVbTagDiscr_);
  mytree_->Branch ("btagjethighEFFbTagDiscr","std::vector<double>",&btagjetHighEFFbTagDiscr_);
  mytree_->Branch ("btagjethighPURbTagDiscr","std::vector<double>",&btagjetHighPURbTagDiscr_);

  // vector of the TLorentz Vectors of met
  metP4_ = new TClonesArray ("TLorentzVector");
  metSig_ = new std::vector<double>;
  mytree_->Branch ("metP4", "TClonesArray", &metP4_, 256000,0);

  // vector of the TLorentz Vectors of tracks
  trackP4_ = new TClonesArray ("TLorentzVector");
  mytree_->Branch ("trackP4", "TClonesArray", &trackP4_, 256000,0);

  // vector of the TLorentz Vectors of other genParticle
  genparticleP4_ = new TClonesArray ("TParticle");
  mytree_->Branch ("genParticleP4", "TClonesArray", &genparticleP4_, 256000,0);

  // vector of the TLorentz Vectors of other genJets
  genjetP4_ = new TClonesArray ("TLorentzVector");
  mytree_->Branch ("genJetP4", "TClonesArray", &genjetP4_, 256000,0);

  // vector of the TLorentz Vectors of other genMet
  genmetP4_ = new TClonesArray ("TLorentzVector");
  mytree_->Branch ("genMetP4", "TClonesArray", &genmetP4_, 256000,0);

}


// ------------ method called once each job just after ending the event loop  ------------


void 
SimpleNtple::endJob() {
  std::cout << "[SimpleNtple::endJob]" << std::endl;
}


// --------------------------------------------------------------------


void SimpleNtple::setVertex (TVector3 &myvector, const TVector3 & mom) {
  myvector.SetX (mom.X());
  myvector.SetY (mom.Y());
  myvector.SetZ (mom.Z());
}


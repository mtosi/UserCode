// -*- C++ -*-
//
// Package:    newSimpleNtple
// Class:      newSimpleNtple
// 
/*

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

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/newSimpleNtple.h"

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
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/PythiaParticleIndex.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace vbfhzz2l2b;

newSimpleNtple::newSimpleNtple(const edm::ParameterSet& iConfig) :
  whichSim_          ( iConfig.getParameter<int> ( "whichSim"  ) ), // 0:FastSim, 1:FullSim
  vertexLabel_       ( iConfig.getParameter<edm::InputTag> ( "vertexLabel"      ) ),
  trackLabel_        ( iConfig.getParameter<edm::InputTag> ( "trackLabel"       ) ),
  muonLabel_         ( iConfig.getParameter<edm::InputTag> ( "muonLabel"        ) ),
  electronLabel_     ( iConfig.getParameter<edm::InputTag> ( "electronLabel"    ) ),
  eleTrkIsoAlgoFlag_ ( iConfig.getParameter<bool>          ( "eleTrkIsoAlgoFlag") ),
  metLabel_          ( iConfig.getParameter<edm::InputTag> ( "metLabel"         ) ),
  tagJetLabel_       ( iConfig.getParameter<edm::InputTag> ( "tagJetLabel"      ) ),
  corIC5CaloJetsWithBTagLabel_ ( iConfig.getParameter<std::string> ( "corIC5CaloJetsWithBTagLabel" ) ),
  corIC5PFJetsWithBTagFlag_    ( iConfig.getParameter<bool>        ( "corIC5PFJetsWithBTagFlag"    ) ),
  genParticleLabel_  ( iConfig.getParameter<edm::InputTag> ( "genParticleLabel" ) ),
  genJetLabel_       ( iConfig.getParameter<edm::InputTag> ( "genJetLabel"      ) ),
  genMetLabel_       ( iConfig.getParameter<edm::InputTag> ( "genMetLabel"      ) ) {

  if ( corIC5PFJetsWithBTagFlag_ ) 
    corIC5PFJetsWithBTagLabel_ = iConfig.getParameter<std::string> ( "corIC5PFJetsWithBTagLabel" );
  

  if ( eleTrkIsoAlgoFlag_ )
    eleTrkIsoAlgo_ = new VBFHZZllbbElectronTrackIsolationAlgos(
		    iConfig.getParameter         <double> ("coneRadius") ,
		    iConfig.getParameter         <double> ("vetoRadius") ,
		    iConfig.getParameter         <double> ("otherVetoRadius") ,
		    iConfig.getParameter         <double> ("ptMin") ,
		    iConfig.getParameter         <double> ("lipMax") ,
		    iConfig.getUntrackedParameter<bool>   ("useTkQuality",true)
		    );

  //now do what ever initialization is needed
  edm::Service<TFileService> fs ;
  mytree_  = fs->make <TTree>("VBFSimpleTree","VBFSimpleTree"); 
  
  //  std::cout << "[newSimpleNtple::newSimpleNtple] DONE" << std::endl;
  
}


// --------------------------------------------------------------------


newSimpleNtple::~newSimpleNtple()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  delete eleTrkIsoAlgo_;

  std::cout << "[newSimpleNtple::~newSimpleNtple]" << std::endl;
  
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
newSimpleNtple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::cout << "[newSimpleNtple::analyze]" << std::endl;

  InitObjs();
  std::cout << "newSimpleNtple::analyze] InitObjs DONE" << std::endl;
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
  if ( whichSim_ == vbfhzz2l2b::FULLSIM )
    FillTrack (iEvent, iSetup);
  FillGenParticle            (iEvent, iSetup); // got an error message in execution
  FillGenJet                 (iEvent, iSetup);
  FillGenMet                 (iEvent, iSetup);
  
  mytree_->Fill();

}


// --------------------------------------------------------------------

void newSimpleNtple::InitObjs() {
  evtID_ = 0;

  // event obj
  CloneEvt_->Clear();
  evt_ = new ((*CloneEvt_)[0]) EVT();
  
  // jet objs
  CloneJet_->Clear();
  for ( unsigned int index = 0; index < 100; index++ ) 
    jet_ = new ((*CloneJet_)[index]) JET();

  // muon objs
  CloneMuon_->Clear();     
  for ( unsigned int index = 0; index < 30; index++ ) 
    muon_ = new ((*CloneMuon_)[index]) MUON();

  // electron objs
  CloneElectron_->Clear();     
  for ( unsigned int index = 0; index < 30; index++ ) 
    electron_ = new ((*CloneElectron_)[index]) ELECTRON();

  // reconstructed hadronic Z objs
  CloneZhad_->Clear();
  for ( unsigned int index = 0; index < 10; index++ ) 
    Zhad_ =new ((*CloneZhad_)[index]) ZHAD();
  
  // tag jets
  invmasstagjetP4_               -> Clear () ;
  invmasstagjetEmEnergyFraction_ -> clear () ;
  invmasstagjetChFrac_           -> clear () ;
  invmasstagjetCorEt_            -> clear () ;
  invmasstagjetCorPt_            -> clear () ;
  invmasstagjetPrimVtxP3_        -> Clear () ;
  deltaetatagjetP4_               -> Clear () ;
  deltaetatagjetEmEnergyFraction_ -> clear () ;
  deltaetatagjetChFrac_           -> clear () ; 
  deltaetatagjetCorEt_            -> clear () ;
  deltaetatagjetCorPt_            -> clear () ;
  deltaetatagjetPrimVtxP3_        -> Clear () ;
  zeptagjetP4_               -> Clear () ;
  zeptagjetEmEnergyFraction_ -> clear () ;
  zeptagjetChFrac_           -> clear () ;
  zeptagjetCorEt_            -> clear () ;
  zeptagjetCorPt_            -> clear () ;
  zeptagjetPrimVtxP3_        -> Clear () ;

  // track
  trackP4_       -> Clear () ;

  // gen particle
  genparticleP4_ 	  -> Clear ();
  genparticlePrimVtxP3_	  -> Clear ();
  genparticlePdgID_	  -> clear ();
  genparticleStatus_	  -> clear ();
  genparticleIndex_	  -> clear ();
  genparticleMomN_	  -> clear ();
  genparticleMomPdgID_	  -> clear ();
  //  genparticleMomPdgIndex_ -> clear ();
  genparticleKidN_	  -> clear ();
  genparticleKidPdgID_	  -> clear ();
  //  genparticleKidPdgIndex_ -> clear ();

  genjetP4_	   -> Clear () ;
  genjetPrimVtxP3_ -> Clear () ;
  genjetKidN_      -> clear () ;
  genjetKidPdgID_  -> clear () ;

  genmetP4_	   -> Clear () ;
  genmetPrimVtxP3_ -> Clear () ;

}

void newSimpleNtple::FillEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  //  std::cout << "[newSimpleNtple::FillEvent]" << std::endl;

  if ( whichSim_ == vbfhzz2l2b::FULLSIM ) {
    edm::Handle<edm::HepMCProduct> evtMC;
    try {
      iEvent.getByLabel("source", evtMC); }
    catch(...) {
      std::cerr << "[newSimpleNtple::FillKindEvent] defined as FullSim, but HepMCProduct::source not found" << std::endl; }
    const HepMC::GenEvent * mcEv = evtMC->GetEvent();
    evtID_ = mcEv->signal_process_id();
  }
  else if ( whichSim_ == vbfhzz2l2b::FASTSIM ) {
    edm::Handle<int> genProcessID;
    try {
      iEvent.getByLabel( "genEventProcID", genProcessID ); }
    catch(...) {
      std::cerr << "[newSimpleNtple::FillKindEvent] defined as FastSim, but genEventProcID not found" << std::endl; }

    evtID_ = *genProcessID;
  }
  else {
    std::cout << "--> WARNING: simulation not specificied!!" << std::endl;
  }

  int run   = iEvent.id().run();
  int event = iEvent.id().event();
  std::cout << "event id: " << iEvent.id() << " -->run: " << run << " and event: " << event <<std::endl;
  evt_->Run     = iEvent.id().run();             // run number
  evt_->Event   = iEvent.id().event();           // event number
  std::cout << "evt_->Run: " << evt_->Run << " and evt_->Event: " << evt_->Event << std::endl;
//  evt_->Ilum    = eventAuxiliaryHandle->luminosityBlock(); // instantaneous luminosity (e30)
  evt_->eventID = evtID_;
//  evt_->nPV;            // number of primary vertex
//  evt_->trigpath;       // Z_BB Trigger Path: 1*main + 10*test1 + 100*test2 (if exists). Ex 101 means main + 2nd test trigger where fired.
//  
//  evt_->P4bquark1;      // first  b-quark
//  evt_->P4bquark2;      // second b-quark
//  evt_->indjetb1;       // first b-quark associated jet index
//  evt_->indjetb2;       // second b-quark associated jet index
//  
//  evt_->pthat;          // pthat
//  
//  evt_->njet;           // number of jets in the event
//  evt_->ngoodjet;       // number of "good" jets in the event
//  evt_->nbtag;          // number of tag in the event
//  evt_->Zvertex;        // Z of the reconstructed primary vertex
//  evt_->P2met;          // Missing Et vector
//
//  evt_->nmuon;          // number of muons in the event	 
//  evt_->nelectron;      // number of electrons in the event
//  evt_->nZhad;          // number of hadronic Z            
  


}


// --------------------------------------------------------------------

void newSimpleNtple::FillElectron(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::cout << "[newSimpleNtple::FillElectron]" << std::endl;

  edm::Handle<reco::PixelMatchGsfElectronCollection> electronHandle ;
  iEvent.getByLabel (electronLabel_,electronHandle) ;

  int nElectron = 0;
  if(electronHandle->size() < 30 ) nElectron = electronHandle->size();
  else nElectron = 30;

  int counter = 0;
  for ( reco::PixelMatchGsfElectronCollection::const_iterator electron_itr = electronHandle->begin();
	electron_itr != electronHandle->end(); ++electron_itr ) {
    vbfhzz2l2b::setMomentum (myvector_, electron_itr->p4());
  }
  
}


// --------------------------------------------------------------------


void newSimpleNtple::FillMuon(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::cout << "[newSimpleNtple::FillMu]" << std::endl;

  edm::Handle<reco::MuonCollection> muonHandle;
  iEvent.getByLabel (muonLabel_,muonHandle);
  
  int nMuon = 0;
  if(muonHandle->size() < 30 ) nMuon = muonHandle->size(); 
  else nMuon = 30;

  int counter = 0;
  for ( reco::MuonCollection::const_iterator muon_itr = muonHandle->begin();
	muon_itr != muonHandle->end(); ++muon_itr ) {
    setMomentum (myvector_, muon_itr->p4());
  }
}


// --------------------------------------------------------------------


void newSimpleNtple::FillMet(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::cout << "[newSimpleNtple::FillMet]" << std::endl;
  edm::Handle<reco::CaloMETCollection> metCollectionHandle;
  iEvent.getByLabel (metLabel_ , metCollectionHandle);
  const CaloMETCollection *calometcol = metCollectionHandle.product();
  const CaloMET *calomet = &(calometcol->front());

//  TClonesArray &MET = *m_MET;
//  setMomentum (myvector_, calomet->p4());
//  new (MET[0]) TLorentzVector (myvector_);
  
}


// --------------------------------------------------------------------


void newSimpleNtple::FillTagJet(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  std::cout << "[newSimpleNtple::FillTagJet]" << std::endl;

  //  std::cout << "[SimpleNtple::FillTagJet]" << std::endl;

  edm::Handle<reco::CaloJetCollection> tagJetHandle;
  iEvent.getByLabel (tagJetLabel_, tagJetHandle) ;

  typedef reco::CaloJetCollection::const_iterator tagJetItr;  

  // looking for the highest invariant mass jets pair
  std::pair<tagJetItr,tagJetItr> maxInvMassPair = 
    vbfhzz2l2b::findPair_maxInvMass_ptMinCut<tagJetItr>(tagJetHandle->begin(), tagJetHandle->end(),
							20., 15.);

  double invMass   = -99.;
  double deltaEta  = -99.;
  double zeppenfeld = -999.;
  int    njets     = 0;
  if (maxInvMassPair.first != maxInvMassPair.second) {
    invMass    =  (     (maxInvMassPair.first)->p4() + ((maxInvMassPair.second)->p4()) ).M();
    deltaEta   =  fabs( (maxInvMassPair.first)->eta() - (maxInvMassPair.second)->eta() );
    zeppenfeld =  (     (maxInvMassPair.first)->pz() * (maxInvMassPair.second)->pz() );
    njets = 2;

    TClonesArray &invmassjetTag = *invmasstagjetP4_;
    vbfhzz2l2b::setMomentum (myvector_, (maxInvMassPair.first)->p4());
    new (invmassjetTag[0]) TLorentzVector (myvector_);
    vbfhzz2l2b::setMomentum (myvector_, (maxInvMassPair.second)->p4());
    new (invmassjetTag[1]) TLorentzVector (myvector_);
    invmasstagjetEmEnergyFraction_->push_back((maxInvMassPair.first)->emEnergyFraction());
    invmasstagjetEmEnergyFraction_->push_back((maxInvMassPair.second)->emEnergyFraction());
  }
  //  invmasstagjetInvMass_->push_back(invMass);
  //  invmasstagjetDeltaEta_->push_back(deltaEta);
  //  invmasstagjetZeppenfeld_->push_back(zeppenfeld);
  //  invmasstagjetN_ = njets;

  // looking for the highest delta eta jets pair
  std::pair<tagJetItr,tagJetItr> maxDeltaEtaPair = 
    vbfhzz2l2b::findPair_maxDeltaEta_ptMinCut<tagJetItr>(tagJetHandle->begin(), tagJetHandle->end(),
							 20., 15.);
  invMass   = -99.;
  deltaEta  = -99.;
  zeppenfeld = -999.;
  njets     = 0;
  if(maxDeltaEtaPair.first != maxDeltaEtaPair.second) {
    invMass    =  (     (maxDeltaEtaPair.first)->p4() + ((maxDeltaEtaPair.second)->p4()) ).M();
    deltaEta   =  fabs( (maxDeltaEtaPair.first)->eta() - (maxDeltaEtaPair.second)->eta() );
    zeppenfeld =  (     (maxDeltaEtaPair.first)->pz() * (maxDeltaEtaPair.second)->pz() );
    njets = 2;

    TClonesArray &deltaetajetTag = *deltaetatagjetP4_;
    vbfhzz2l2b::setMomentum (myvector_, (maxDeltaEtaPair.first)->p4());
    new (deltaetajetTag[0]) TLorentzVector (myvector_);
    vbfhzz2l2b::setMomentum (myvector_, (maxDeltaEtaPair.second)->p4());
    new (deltaetajetTag[1]) TLorentzVector (myvector_);
    deltaetatagjetEmEnergyFraction_->push_back((maxDeltaEtaPair.first)->emEnergyFraction());
    deltaetatagjetEmEnergyFraction_->push_back((maxDeltaEtaPair.second)->emEnergyFraction());
  }
  //  deltaetatagjetInvMass_->push_back(invMass);
  //  deltaetatagjetDeltaEta_->push_back(deltaEta);
  //  deltaetatagjetZeppenfeld_->push_back(zeppenfeld);
  //  deltaetatagjetN_ = njets;


  // looking for the highest zeppenfeld variable value jets pair
  std::pair<tagJetItr,tagJetItr> maxZepPair = 
    vbfhzz2l2b::findPair_maxZeppenfeld_ptMinCut<tagJetItr>(tagJetHandle->begin(), tagJetHandle->end(),
							   20., 15.);

  invMass   = -99.;
  deltaEta  = -99.;
  zeppenfeld = -999.;
  njets     = 0;
  if (maxZepPair.first != maxZepPair.second) {
    invMass    =  (     (maxInvMassPair.first)->p4() + ((maxInvMassPair.second)->p4()) ).M();
    deltaEta   =  fabs( (maxInvMassPair.first)->eta() - (maxInvMassPair.second)->eta() );
    zeppenfeld =  (     (maxInvMassPair.first)->pz() * (maxInvMassPair.second)->pz() );
    njets = 2;

    TClonesArray &zepjetTag = *zeptagjetP4_;
    vbfhzz2l2b::setMomentum (myvector_, (maxZepPair.first)->p4());
    new (zepjetTag[0]) TLorentzVector (myvector_);
    vbfhzz2l2b::setMomentum (myvector_, (maxZepPair.second)->p4());
    new (zepjetTag[1]) TLorentzVector (myvector_);
    zeptagjetEmEnergyFraction_->push_back((maxZepPair.first)->emEnergyFraction());
    zeptagjetEmEnergyFraction_->push_back((maxZepPair.second)->emEnergyFraction());
  }

  //  zeptagjetInvMass_->push_back(invMass);
  //  zeptagjetDeltaEta_->push_back(deltaEta);
  //  zeptagjetZeppenfeld_->push_back(zeppenfeld);
  //  zeptagjetN_ = njets;
}
// --------------------------------------------------------------------

void newSimpleNtple::FillcorIC5CaloJetsWithBTag(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  std::cout << "[newSimpleNtple::FillcorIC5CaloJetsWithBTag]" << std::endl;

  edm::Handle<vbfhzz2l2b::CorJetWithBTagDiscrCollection> corJetWithBTagHandle;
  iEvent.getByLabel(corIC5CaloJetsWithBTagLabel_,"corJetWithBTagDiscr",corJetWithBTagHandle);
  //  std::cout << "corJetWithBTagHandle->size(): " << corJetWithBTagHandle->size() << std::endl;

  int jetIndex = 0;    
  std::vector<reco::JetBaseRef> jets = vbfhzz2l2b::CorJetBTagDiscrAssociation::allJets(*corJetWithBTagHandle);
  for ( std::vector<reco::JetBaseRef>::const_iterator jet = jets.begin(); 
	jet != jets.end();
	++jet, jetIndex++ ) {
    //    std::cout << "jetIndex: " << jetIndex << std::endl;
    
    std::vector<double> discrVec = (*corJetWithBTagHandle)[*jet].discrVec_;
    double corrEt = (*corJetWithBTagHandle)[*jet].corEt_;
    double uncorrEt = (*jet)->et();
    //    std::cout << "uncorrEt: " << uncorrEt << std::endl;
    //    std::cout << "corrEt: " << corrEt << std::endl;
    //    for ( unsigned int index = 0; index != discrVec.size(); index++) 
    //      std::cout << "discr[" << index << "]: " << discrVec[index] << std::endl;

    double emFrac = (dynamic_cast<const reco::CaloJet*>(&**jet))->emEnergyFraction();
    //    std::cout << "emFrac: " << emFrac << std::endl;

//     std::cout << "corJetWithBTag highEffDiscr: "       << corJet->discriminators().discriminatorHighEff()       << std::endl;
//     std::cout << "corJetWithBTag highPurDiscr: "       << corJet->discriminators().discriminatorHighPur()       << std::endl;
//     std::cout << "corJetWithBTag combSecVtxDiscr: "    << corJet->discriminators().discriminatorCombSecVtx()    << std::endl;
//     //  std::cout << "corJetWithBTag combSecVtxMVADiscr: " << corJet->discriminators().discriminatorCombSecVtxMVA() << std::endl;
//     //  std::cout << "corJetWithBTag softMuonDiscr: "      << corJet->discriminators().discriminatorSoftMuon()      << std::endl;
//     //  std::cout << "corJetWithBTag softElectronDiscr: "  << corJet->discriminators().discriminatorSoftElectron()  << std::endl;
//     std::cout << "corJetWithBTag jetProbDiscr: "       << corJet->discriminators().discriminatorJetProb()       << std::endl;
//     

    setMomentum (myvector_, (*jet)->p4());      
  }

}

void newSimpleNtple::FillcorIC5PFJetsWithBTag(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
}
// --------------------------------------------------------------------

void newSimpleNtple::FillZhad(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
}

void newSimpleNtple::FillZlep(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
}

void  newSimpleNtple::FillTrack(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::cout << "[newSimpleNtple::FillTracks]" << std::endl;

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


void newSimpleNtple::FillGenParticle(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::cout << "[newSimpleNtple::FillGenParticles]" << std::endl;
  edm::Handle<reco::GenParticleCollection> genParticleHandle; 
  iEvent.getByLabel (genParticleLabel_,genParticleHandle);
 
  TClonesArray &genParticleP4        = *genparticleP4_;
  TClonesArray &genparticlePrimVtxP3 = *genparticlePrimVtxP3_;
  if ( evtID_ == HWWFusion_ || evtID_ == HZZFusion_ ){ //---- only if VBF

    int particleIndex = 0;
    for (reco::GenParticleCollection::const_iterator particle_itr = genParticleHandle->begin(); 
	 particle_itr != genParticleHandle->end(); ++particle_itr, particleIndex++ ) {

      if ( fabs(particle_itr->pdgId()) <= pythiaH_ ) {
	int pdgID = particle_itr->pdgId();
	int status = particle_itr->status();
	int mother1pdgID = 0;
	if ( particle_itr->numberOfMothers() > 0 ) mother1pdgID = particle_itr->mother(0)->pdgId();
	int mother2pdgID = 0;
	if ( particle_itr->numberOfMothers() > 1 ) mother2pdgID = particle_itr->mother(1)->pdgId();
	int daughter1pdgID = 0;
	if ( particle_itr->numberOfDaughters() > 0 ) daughter1pdgID = particle_itr->daughter(0)->pdgId();
	int daughter2pdgID = 0;
	if ( particle_itr->numberOfDaughters() > 1 ) daughter2pdgID = particle_itr->daughter(1)->pdgId();
	double px = particle_itr->px();
	double py = particle_itr->py();
	double pz = particle_itr->pz();
	double e  = particle_itr->energy();
	double vx = particle_itr->vx();
	double vy = particle_itr->vy();
	double vz = particle_itr->vz();
	double time = 0.;
	
	new (genParticleP4[particleIndex]) TParticle (pdgID, status, 
						      mother1pdgID, mother2pdgID,
						      daughter1pdgID, daughter2pdgID,
						      px, py, pz, e,
						      vx, vy, vz, time);
	vbfhzz2l2b::setVertex (myvertex_, particle_itr->vertex());
	new (genparticlePrimVtxP3[particleIndex]) TVector3 (myvertex_);

	genparticlePdgID_       -> push_back (pdgID);
	genparticleStatus_      -> push_back (status);
	genparticleIndex_       -> push_back (particleIndex);
	genparticleMomN_        -> push_back (particle_itr->numberOfMothers());
	genparticleMomPdgID_    -> push_back (std::pair<int,int>(mother1pdgID,mother2pdgID));
	//genparticleMomPdgIndex_ -> push_back ();
	genparticleKidN_	      -> push_back (particle_itr->numberOfDaughters());
	genparticleKidPdgID_    -> push_back (std::pair<int,int>(mother1pdgID,mother2pdgID));
	//genparticleKidPdgIndex_ -> push_back (); 

      } // if pdgId <= pythiaH_
    } // loop over genParticles
  } // if VBF
}




// --------------------------------------------------------------------


void newSimpleNtple::FillGenJet(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::cout << "[newSimpleNtple::FillGenJet]" << std::endl;
  edm::Handle< reco::GenJetCollection > genJetHandle ;
  iEvent.getByLabel( genJetLabel_, genJetHandle ) ;
  
  TClonesArray &genjetP4        = *genjetP4_;
  TClonesArray &genjetPrimVtxP3 = *genjetPrimVtxP3_;

  int genjetIndex = 0;
  for (reco::GenJetCollection::const_iterator genjet_itr = genJetHandle->begin (); 
       genjet_itr != genJetHandle->end (); ++genjet_itr, genjetIndex++ ) { 

    std::cout << "genjet_itr->numberOfDaughters(): " << genjet_itr->numberOfDaughters() << std::endl;
    std::vector<int> daughterPdgIDVec;
    for (std::vector<edm::Ptr< reco::Candidate> >::const_iterator daughter_itr = (genjet_itr->daughterPtrVector()).begin();
	 daughter_itr != (genjet_itr->daughterPtrVector()).end(); ++daughter_itr ) {

      daughterPdgIDVec.push_back ((*daughter_itr)->pdgId());
    }

    vbfhzz2l2b::setMomentum (myvector_, genjet_itr->p4());
    vbfhzz2l2b::setVertex   (myvertex_, genjet_itr->vertex());
    new (genjetP4[genjetIndex])        TLorentzVector (myvector_);
    new (genjetPrimVtxP3[genjetIndex]) TVector3       (myvertex_);

    genjetKidN_     -> push_back (genjet_itr->numberOfDaughters());
    genjetKidPdgID_ -> push_back (daughterPdgIDVec);

  }

  //  std::cout << "[SimpleNtple::FillGenJet] DONE" << std::endl;
}



// --------------------------------------------------------------------


void newSimpleNtple::FillGenMet(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::cout << "[newSimpleNtple::FillGenMet]" << std::endl;

  edm::Handle< reco::GenMETCollection > genMetCollectionHandle ;
  iEvent.getByLabel( genMetLabel_, genMetCollectionHandle ) ;
  
  const reco::GenMETCollection *genmetcol = genMetCollectionHandle.product();
  const reco::GenMET *genmet = &(genmetcol->front());

  TClonesArray &genmetP4        = *genmetP4_;
  TClonesArray &genmetPrimVtxP3 = *genmetPrimVtxP3_;
  vbfhzz2l2b::setMomentum (myvector_, genmet->p4());
  vbfhzz2l2b::setVertex   (myvertex_, genmet->vertex());
  new (genmetP4[0])        TLorentzVector (myvector_);
  new (genmetPrimVtxP3[0]) TVector3       (myvertex_);

}

// ------------ method called once each job just before starting event loop  ------------


void newSimpleNtple::beginJob(const edm::EventSetup& iSetup)
{

  std::cout << "[newSimpleNtple::beginJob]" << std::endl;

  CloneEvt_ = new TClonesArray("EVT");
  mytree_->Branch("EVT","TClonesArray",&CloneEvt_,256000,0);
  CloneJet_ = new TClonesArray("JET");
  mytree_->Branch("JET","TClonesArray",&CloneJet_,256000,0);
  CloneMuon_ = new TClonesArray("MUON");
  mytree_->Branch("MUON","TClonesArray",&CloneMuon_,256000,0);
  CloneElectron_ = new TClonesArray("ELECTRON");
  mytree_->Branch("ELECTRON","TClonesArray",&CloneElectron_,256000,0);
  CloneZhad_ = new TClonesArray("ZHAD");
  mytree_->Branch("ZHAD","TClonesArray",&CloneZhad_,256000,0);
  

  std::cout << "[newSimpleNtple::beginJob] DONE" << std::endl;

  // vector of the TLorentz Vectors of tag jets with inv mass criteria
  invmasstagjetP4_        = new TClonesArray ("TLorentzVector");
  invmasstagjetPrimVtxP3_ = new TClonesArray ("TVector3");
  invmasstagjetEmEnergyFraction_ = new std::vector<double>;
  invmasstagjetChFrac_           = new std::vector<double>;
  invmasstagjetCorEt_            = new std::vector<double>;
  invmasstagjetCorPt_            = new std::vector<double>;
  invmasstagjetCompoSVbTagDiscr_ = new std::vector<double>;
  invmasstagjetHighEFFbTagDiscr_ = new std::vector<double>;
  invmasstagjetHighPURbTagDiscr_ = new std::vector<double>;
  mytree_->Branch ("invmasstagjetP4",              "TClonesArray",       &invmasstagjetP4_,        256000,0);
  mytree_->Branch ("invmasstagjetPrimVtxP3",       "TClonesArray",       &invmasstagjetPrimVtxP3_, 256000,0);
  mytree_->Branch ("invmasstagjetEmFrac",          "std::vector<double>",&invmasstagjetEmEnergyFraction_);
  mytree_->Branch ("invmasstagjetChFrac",          "std::vector<double>",&invmasstagjetChFrac_);
  mytree_->Branch ("invmasstagjetCorEt",           "std::vector<double>",&invmasstagjetCorEt_);
  mytree_->Branch ("invmasstagjetCorPt",           "std::vector<double>",&invmasstagjetCorPt_);
  mytree_->Branch ("invmasstagjetcompoSVbTagDiscr","std::vector<double>",&invmasstagjetCompoSVbTagDiscr_);
  mytree_->Branch ("invmasstagjethighEFFbTagDiscr","std::vector<double>",&invmasstagjetHighEFFbTagDiscr_);
  mytree_->Branch ("invmasstagjethighPURbTagDiscr","std::vector<double>",&invmasstagjetHighPURbTagDiscr_);

  // vector of the TLorentz Vectors of tag jets with inv mass criteria
  deltaetatagjetP4_        = new TClonesArray ("TLorentzVector");
  deltaetatagjetPrimVtxP3_ = new TClonesArray ("TVector3");
  deltaetatagjetEmEnergyFraction_ = new std::vector<double>;
  deltaetatagjetChFrac_           = new std::vector<double>;
  deltaetatagjetCorEt_            = new std::vector<double>;
  deltaetatagjetCorPt_            = new std::vector<double>;
  deltaetatagjetCompoSVbTagDiscr_ = new std::vector<double>;
  deltaetatagjetHighEFFbTagDiscr_ = new std::vector<double>;
  deltaetatagjetHighPURbTagDiscr_ = new std::vector<double>;
  mytree_->Branch ("deltaetatagjetP4",              "TClonesArray",       &deltaetatagjetP4_,        256000,0);
  mytree_->Branch ("deltaetatagjetPrimVtxP3",       "TClonesArray",       &deltaetatagjetPrimVtxP3_, 256000,0);
  mytree_->Branch ("deltaetatagjetEmFrac",          "std::vector<double>",&deltaetatagjetEmEnergyFraction_);
  mytree_->Branch ("deltaetatagjetChFrac",          "std::vector<double>",&deltaetatagjetChFrac_);
  mytree_->Branch ("deltaetatagjetCorEt",           "std::vector<double>",&deltaetatagjetCorEt_);
  mytree_->Branch ("deltaetatagjetCorPt",           "std::vector<double>",&deltaetatagjetCorPt_);
  mytree_->Branch ("deltaetatagjetcompoSVbTagDiscr","std::vector<double>",&deltaetatagjetCompoSVbTagDiscr_);
  mytree_->Branch ("deltaetatagjethighEFFbTagDiscr","std::vector<double>",&deltaetatagjetHighEFFbTagDiscr_);
  mytree_->Branch ("deltaetatagjethighPURbTagDiscr","std::vector<double>",&deltaetatagjetHighPURbTagDiscr_);


  // vector of the TLorentz Vectors of tag jets with inv mass criteria
  zeptagjetP4_        = new TClonesArray ("TLorentzVector");
  zeptagjetPrimVtxP3_ = new TClonesArray ("TVector3");
  zeptagjetEmEnergyFraction_ = new std::vector<double>;
  zeptagjetChFrac_           = new std::vector<double>;
  zeptagjetCorEt_            = new std::vector<double>;
  zeptagjetCorPt_            = new std::vector<double>;
  zeptagjetCompoSVbTagDiscr_ = new std::vector<double>;
  zeptagjetHighEFFbTagDiscr_ = new std::vector<double>;
  zeptagjetHighPURbTagDiscr_ = new std::vector<double>;
  mytree_->Branch ("zeptagjetP4",              "TClonesArray",       &zeptagjetP4_,        256000,0);
  mytree_->Branch ("zeptagjetPrimVtxP3",       "TClonesArray",       &zeptagjetPrimVtxP3_, 256000,0);
  mytree_->Branch ("zeptagjetEmFrac",          "std::vector<double>",&zeptagjetEmEnergyFraction_);
  mytree_->Branch ("zeptagjetChFrac",          "std::vector<double>",&zeptagjetChFrac_);
  mytree_->Branch ("zeptagjetCorEt",           "std::vector<double>",&zeptagjetCorEt_);
  mytree_->Branch ("zeptagjetCorPt",           "std::vector<double>",&zeptagjetCorPt_);
  mytree_->Branch ("zeptagjetcompoSVbTagDiscr","std::vector<double>",&zeptagjetCompoSVbTagDiscr_);
  mytree_->Branch ("zeptagjethighEFFbTagDiscr","std::vector<double>",&zeptagjetHighEFFbTagDiscr_);
  mytree_->Branch ("zeptagjethighPURbTagDiscr","std::vector<double>",&zeptagjetHighPURbTagDiscr_);

  // vector of the TLorentz Vectors of tracks
  trackP4_ = new TClonesArray ("TLorentzVector");
  mytree_->Branch ("trackP4", "TClonesArray", &trackP4_, 256000,0);


  genparticleP4_          = new TClonesArray ("TParticle");
  //  genparticleP4_          = new TClonesArray ("TLorentzVector");
  genparticlePrimVtxP3_   = new TClonesArray ("TVector3");
  genparticlePdgID_       = new std::vector<int>;
  genparticleStatus_      = new std::vector<int>;
  genparticleIndex_	  = new std::vector<int>;
  genparticleMomN_	  = new std::vector<int>;
  genparticleMomPdgID_	  = new std::vector< std::pair<int,int> >;
  //  genparticleMomPdgIndex_ = new std::vector< std::pair<int,int> >;
  genparticleKidN_	  = new std::vector<int>;
  genparticleKidPdgID_	  = new std::vector< std::pair<int,int> >;
  //  genparticleKidPdgIndex_ = new std::vector< std::pair<int,int> >;
  mytree_->Branch("genparticleP4",         "TClonesArray",                     &genparticleP4_,       256000,0);
  mytree_->Branch("genparticlePrimVtxP3",  "TClonesArray",                     &genparticlePrimVtxP3_,256000,0);
  mytree_->Branch("genparticlePdgID",      "std::vector<int>",                 &genparticlePdgID_);
  mytree_->Branch("genparticleStatus",     "std::vector<int>",                 &genparticleStatus_);
  mytree_->Branch("genparticleIndex",      "std::vector<int>",                 &genparticleIndex_);
  mytree_->Branch("genparticleMomN",       "std::vector<int>",                 &genparticleMomN_);
  mytree_->Branch("genparticleMomPdgID",   "std::vector< std::pair<int,int> >",&genparticleMomPdgID_);
  //mytree_->Branch("genparticleMomPdgIndex","std::vector< std::pair<int,int> >",&genparticleMomPdgIndex_);
  mytree_->Branch("genparticleKidN",       "std::vector<int>",                 &genparticleKidN_);
  mytree_->Branch("genparticleKidPdgID",   "std::vector< std::pair<int,int> >",&genparticleKidPdgID_);
  //  mytree_->Branch("genparticleKidPdgIndex","std::vector< std::pair<int,int> >",&genparticleKidPdgIndex_);


  // vector of the TLorentz Vectors of other genJets
  genjetP4_        = new TClonesArray ("TLorentzVector");
  genjetPrimVtxP3_ = new TClonesArray ("TVector3");
  genjetKidN_     = new std::vector<int>;
  genjetKidPdgID_ = new std::vector< std::vector<int> >;
  mytree_->Branch ("genjetP4",        "TClonesArray",                   &genjetP4_,        256000,0);
  mytree_->Branch ("genjetPrimVtxP3", "TClonesArray",                   &genjetPrimVtxP3_, 256000,0);
  mytree_->Branch ("genjetKidN",      "std::vector<int>",               &genjetKidN_    );   
  mytree_->Branch ("genjetKidPdgID",  "std::vector< std::vector<int> >",&genjetKidPdgID_);

  // vector of the TLorentz Vectors of other genMet
  genmetP4_        = new TClonesArray ("TLorentzVector");
  genmetPrimVtxP3_ = new TClonesArray ("TVector3");
  mytree_->Branch ("genmetP4",       "TClonesArray", &genmetP4_,        256000,0);
  mytree_->Branch ("genmetPrimVtxP3","TClonesArray", &genmetPrimVtxP3_, 256000,0);


}


// ------------ method called once each job just after ending the event loop  ------------


void 
newSimpleNtple::endJob() {
  //  mytree_->Write();
  std::cout << "[newSimpleNtple::endJob]" << std::endl;

  delete invmasstagjetP4_;
  delete invmasstagjetPrimVtxP3_;
  delete invmasstagjetEmEnergyFraction_;
  delete invmasstagjetChFrac_;
  delete invmasstagjetCorEt_;
  delete invmasstagjetCorPt_;
  delete invmasstagjetCompoSVbTagDiscr_;
  delete invmasstagjetHighEFFbTagDiscr_;
  delete invmasstagjetHighPURbTagDiscr_;

  delete deltaetatagjetP4_;
  delete deltaetatagjetPrimVtxP3_;
  delete deltaetatagjetEmEnergyFraction_;
  delete deltaetatagjetChFrac_;
  delete deltaetatagjetCorEt_;
  delete deltaetatagjetCorPt_;
  delete deltaetatagjetCompoSVbTagDiscr_;
  delete deltaetatagjetHighEFFbTagDiscr_;
  delete deltaetatagjetHighPURbTagDiscr_;

  delete zeptagjetP4_;
  delete zeptagjetPrimVtxP3_;
  delete zeptagjetEmEnergyFraction_;
  delete zeptagjetChFrac_;
  delete zeptagjetCorEt_;
  delete zeptagjetCorPt_;
  delete zeptagjetCompoSVbTagDiscr_;
  delete zeptagjetHighEFFbTagDiscr_;
  delete zeptagjetHighPURbTagDiscr_;

  delete trackP4_ ;

  delete genparticleP4_;         
  delete genparticlePrimVtxP3_;      
  delete genparticlePdgID_;      
  delete genparticleStatus_;      
  delete genparticleIndex_;	 
  delete genparticleMomN_;	 
  delete genparticleMomPdgID_;	 
  //  delete genparticleMomPdgIndex_;
  delete genparticleKidN_;	 
  delete genparticleKidPdgID_;	 
  //  delete genparticleKidPdgIndex_;

  delete genjetP4_;
  delete genjetPrimVtxP3_;
  delete genjetKidN_;
  delete genjetKidPdgID_;

  delete genmetP4_;
  delete genmetPrimVtxP3_;
  
  //  std::cout << "[newSimpleNtple::endJob]" << std::endl;


}


// --------------------------------------------------------------------

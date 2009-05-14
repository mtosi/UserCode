// -*- C++ -*-
//
// Package:    SimpleNtple
// Class:      SimpleNtple
// 

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

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/VBFHZZllbbUtils.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/ProcessIndex.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace vbfhzz2l2b;
using namespace vbfhzz2l2b::SimpleNtpleObj;


enum { FASTSIM = 0,
       FULLSIM = 1
};

newSimpleNtple::newSimpleNtple(const edm::ParameterSet& iConfig) :
  whichSim_ ( iConfig.getParameter<int> ( "whichSim"  ) ), // 0:FastSim, 1:FullSim
  tracksLabel_   ( iConfig.getParameter<edm::InputTag> ( "tracksLabel"   ) ),
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
 
  std::cout << "[newSimpleNtple::newSimpleNtple]" << std::endl;
  //now do what ever initialization is needed
  /*
  edm::Service<TFileService> fs ;
  mytree_  = fs->make <TTree>("VBFSimpleTree","VBFSimpleTree"); 
  */
  mytree_  = new TTree("VBFSimpleTree","VBFSimpleTree"); 
  std::cout << "[newSimpleNtple::newSimpleNtple] DONE" << std::endl;

}


// --------------------------------------------------------------------


newSimpleNtple::~newSimpleNtple()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

//  delete CloneEvt_;
//  delete CloneJet_;
//  delete CloneMuon_;
//  delete CloneElectron_;
//  delete CloneZhad_;

  delete m_tagJets;
  delete m_MET;
  delete m_tracks;
  delete m_genParticles;
  delete m_genJets;
  delete m_genMet;

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
  FillTagJet                 (iEvent, iSetup);
  if ( corIC5PFJetsWithBTagFlag_ )
    FillcorIC5PFJetsWithBTag (iEvent, iSetup);   
  if ( whichSim_ == FULLSIM )
    FillTracks (iEvent, iSetup);
  FillGenParticles (iEvent, iSetup);
  FillGenJet       (iEvent, iSetup);
  FillGenMet       (iEvent, iSetup);
  
  mytree_->Fill();
  
}


// --------------------------------------------------------------------

void newSimpleNtple::InitObjs() {
  evtID_ = 0;

  /*
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
  */
  m_tagJets      -> Clear () ;
  m_MET          -> Clear () ;
  m_tracks       -> Clear () ;
  m_genParticles -> Clear () ;
  m_genJets      -> Clear () ;
  m_genMet       -> Clear () ;
  
}

void newSimpleNtple::FillEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  edm::Handle<edm::EventAuxiliary> eventAuxiliaryHandle;
  iEvent.getByLabel("EventAuxiliary",eventAuxiliaryHandle);
  std::cout << "eventAuxiliaryHandle->run(): " << eventAuxiliaryHandle->run() << std::endl;

  std::cout << "[newSimpleNtple::FillEvent]" << std::endl;
  if ( whichSim_ == FULLSIM ) {
    edm::Handle<edm::HepMCProduct> evtMC;
    try {
      iEvent.getByLabel("source", evtMC); }
    catch(...) {
      std::cerr << "[newSimpleNtple::FillKindEvent] defined as FullSim, but HepMCProduct::source not found" << std::endl; }
    const HepMC::GenEvent * mcEv = evtMC->GetEvent();
    evtID_ = mcEv->signal_process_id();
  }
  else if ( whichSim_ == FASTSIM ) {
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


//  evt_->Run     = eventAuxiliaryHandle->run();             // run number
//  evt_->Event   = eventAuxiliaryHandle->event();           // event number
//  evt_->Ilum    = eventAuxiliaryHandle->luminosityBlock(); // instantaneous luminosity (e30)
//  evt_->eventID = evtID;
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

  TClonesArray &MET = *m_MET;
  setMomentum (myvector_, calomet->p4());
  new (MET[0]) TLorentzVector (myvector_);
  
}


// --------------------------------------------------------------------


void newSimpleNtple::FillTagJet(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  std::cout << "[newSimpleNtple::FillTagJet]" << std::endl;

  edm::Handle<reco::RecoChargedCandidateCollection> tagJetHandle;
  iEvent.getByLabel (tagJetLabel_, tagJetHandle) ;
  
  math::XYZTLorentzVector sumLV = (*tagJetHandle)[0].p4() + (*tagJetHandle)[1].p4() ;
  invMassTagJet_ = sumLV.M();

  TClonesArray &jetTag = *m_tagJets;
  int counter = 0;
  for (RecoChargedCandidateCollection::const_iterator jet_itr = tagJetHandle->begin (); 
       jet_itr != tagJetHandle->end (); ++jet_itr ) { 
    
    vbfhzz2l2b::setMomentum (myvector_, jet_itr->p4());
    new (jetTag[counter]) TLorentzVector (myvector_);
    counter++;
  }
}
// --------------------------------------------------------------------

void newSimpleNtple::FillcorIC5CaloJetsWithBTag(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  std::cout << "[newSimpleNtple::FillcorIC5CaloJetsWithBTag]" << std::endl;

  edm::Handle<vbfhzz2l2b::CorJetWithBTagDiscrCollection> corJetWithBTagHandle;
  iEvent.getByLabel(corIC5CaloJetsWithBTagLabel_,"corJetWithBTagDiscr",corJetWithBTagHandle);
  std::cout << "corJetWithBTagHandle->size(): " << corJetWithBTagHandle->size() << std::endl;

  int jetIndex = 0;    
  std::vector<reco::JetBaseRef> jets = vbfhzz2l2b::CorJetBTagDiscrAssociation::allJets(*corJetWithBTagHandle);
  for ( std::vector<reco::JetBaseRef>::const_iterator jet = jets.begin(); 
	jet != jets.end();
	++jet, jetIndex++ ) {
    std::cout << "jetIndex: " << jetIndex << std::endl;
    
    std::vector<double> discrVec = (*corJetWithBTagHandle)[*jet].discrVec_;
    double corrEt = (*corJetWithBTagHandle)[*jet].corEt_;
    double uncorrEt = (*jet)->et();
    std::cout << "uncorrEt: " << uncorrEt << std::endl;
    std::cout << "corrEt: " << corrEt << std::endl;
    for ( unsigned int index = 0; index != discrVec.size(); index++) 
      std::cout << "discr[" << index << "]: " << discrVec[index] << std::endl;

    double emFrac = (dynamic_cast<const reco::CaloJet*>(&**jet))->emEnergyFraction();
    std::cout << "emFrac: " << emFrac << std::endl;

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

void  newSimpleNtple::FillTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::cout << "[newSimpleNtple::FillTracks]" << std::endl;
  edm::Handle<reco::TrackCollection> tracksHandle ;
  iEvent.getByLabel (tracksLabel_, tracksHandle) ;

  TClonesArray &tracks = *m_tracks;
  int counter = 0;
  for (reco::TrackCollection::const_iterator track_itr = tracksHandle->begin (); 
       track_itr != tracksHandle->end (); ++track_itr ) { 

    math::XYZVector mom = track_itr->innerMomentum () ; 
    myvector_.SetPx (mom.x ()) ;
    myvector_.SetPy (mom.y ()) ;
    myvector_.SetPz (mom.z ()) ;
    
    new (tracks[counter]) TLorentzVector (myvector_);
    counter++;
  }
}


// --------------------------------------------------------------------


void newSimpleNtple::FillGenParticles(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::cout << "[newSimpleNtple::FillGenParticles]" << std::endl;
  edm::Handle<reco::GenParticleCollection> genParticlesHandle; 
  iEvent.getByLabel (genParticleLabel_,genParticlesHandle);
 
  TClonesArray &genParticles = *m_genParticles;
  int counter = 0;
 
  if ( evtID_ == HWWFusion_ || evtID_ == HZZFusion_ ){ //---- only if VBF
    for (reco::GenParticleCollection::const_iterator particle_itr = genParticlesHandle->begin(); 
	 particle_itr != genParticlesHandle->end(); ++particle_itr ) {
   
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
  
   new (genParticles[counter]) TParticle (pdg, status, 
					  mother1, mother2,
					  daughter1, daughter2,
					  px, py, pz, e,
					  vx, vy, vz, time);
  
   counter++;
  }
 }
}




// --------------------------------------------------------------------


void newSimpleNtple::FillGenJet(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::cout << "[newSimpleNtple::FillGenJet]" << std::endl;
  edm::Handle< reco::GenJetCollection > genJetsHandle ;
  iEvent.getByLabel( genJetLabel_, genJetsHandle ) ;

 TClonesArray &genJets = *m_genJets;
 int counter = 0;
 for (reco::GenJetCollection::const_iterator genjet_itr = genJetsHandle->begin (); 
      genjet_itr != genJetsHandle->end (); ++genjet_itr ) { 
  
  myvector_.SetPx ( genjet_itr->px() );
  myvector_.SetPy ( genjet_itr->py() );
  myvector_.SetPz ( genjet_itr->pz() );
  myvector_.SetE  ( genjet_itr->emEnergy() + genjet_itr->hadEnergy() );
  new (genJets[counter]) TLorentzVector (myvector_);
  counter++;
 }
}



// --------------------------------------------------------------------


void newSimpleNtple::FillGenMet(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::cout << "newSimpleNtple::FillGenMet" << std::endl;
 edm::Handle< reco::GenMETCollection > genMetHandle ;
 iEvent.getByLabel( genMetLabel_, genMetHandle ) ;

 TClonesArray &genMets = *m_genMet;
 int counter = 0;
 for (reco::GenMETCollection::const_iterator gMIt = genMetHandle->begin (); 
      gMIt != genMetHandle->end (); 
      ++gMIt ) 
 { 
  myvector_.SetPx ((*gMIt).px ()) ;
  myvector_.SetPy ((*gMIt).py ()) ;
  myvector_.SetPz ((*gMIt).pz ()) ;
  myvector_.SetE ((*gMIt).emEnergy () + (*gMIt).hadEnergy ()) ;
  new (genMets[counter]) TLorentzVector (myvector_);
  counter++;
 }
}

// ------------ method called once each job just before starting event loop  ------------


void newSimpleNtple::beginJob(const edm::EventSetup& iSetup)
{

  std::cout << "[newSimpleNtple::beginJob]" << std::endl;

  //  CloneEvt_ = new TClonesArray("EVT");
  std::cout << "[newSimpleNtple::beginJob] DONE w/ CloneEvt_" << std::endl;
  /*
  mytree_->Branch("EVT","TClonesArray",&CloneEvt_,256000,0);
  std::cout << "[newSimpleNtple::beginJob] DONE w/ CloneEvt_" << std::endl;
  CloneJet_ = new TClonesArray("JET");
  mytree_->Branch("JET","TClonesArray",&CloneJet_,256000,0);
  CloneMuon_ = new TClonesArray("MUON");
  mytree_->Branch("MUON","TClonesArray",&CloneMuon_,256000,0);
  CloneElectron_ = new TClonesArray("ELECTRON");
  mytree_->Branch("ELECTRON","TClonesArray",&CloneElectron_,256000,0);
  CloneZhad_ = new TClonesArray("ZHAD");
  mytree_->Branch("ZHAD","TClonesArray",&CloneZhad_,256000,0);
  */

  // vector with the 2 tag TLorentzVectors
  m_tagJets = new TClonesArray ("TLorentzVector");
  mytree_->Branch ("tagJets", "TClonesArray", &m_tagJets, 256000,0);

  // vector of the TLorentz Vectors of other jets
  m_MET = new TClonesArray ("TLorentzVector");
  mytree_->Branch ("MET", "TClonesArray", &m_MET, 256000,0);

  // vector of the TLorentz Vectors of other jets
  m_tracks = new TClonesArray ("TLorentzVector");
  mytree_->Branch ("tracks", "TClonesArray", &m_tracks, 256000,0);

  // vector of the TLorentz Vectors of other genParticles
  m_genParticles = new TClonesArray ("TParticle");
  mytree_->Branch ("genParticles", "TClonesArray", &m_genParticles, 256000,0);

  // vector of the TLorentz Vectors of other genJets
  m_genJets = new TClonesArray ("TLorentzVector");
  mytree_->Branch ("genJets", "TClonesArray", &m_genJets, 256000,0);

  // vector of the TLorentz Vectors of other genMet
  m_genMet = new TClonesArray ("TLorentzVector");
  mytree_->Branch ("genMet", "TClonesArray", &m_genMet, 256000,0);


  std::cout << "[newSimpleNtple::beginJob] DONE" << std::endl;



}


// ------------ method called once each job just after ending the event loop  ------------


void 
newSimpleNtple::endJob() {
  std::cout << "[newSimpleNtple::endJob]" << std::endl;
}


// --------------------------------------------------------------------

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


SimpleNtple::SimpleNtple(const edm::ParameterSet& iConfig) :
  whichSim_          ( iConfig.getParameter<int> ( "whichSim"  ) ), // 0:FastSim, 1:FullSim
  vertexLabel_       ( iConfig.getParameter<edm::InputTag> ( "vertexLabel"      ) ),
  trackLabel_        ( iConfig.getParameter<edm::InputTag> ( "trackLabel"       ) ),
  muonLabel_         ( iConfig.getParameter<edm::InputTag> ( "muonLabel"        ) ),
  electronLabel_     ( iConfig.getParameter<edm::InputTag> ( "electronLabel"    ) ),
  metLabel_          ( iConfig.getParameter<edm::InputTag> ( "metLabel"         ) ),
  tagJetLabel_       ( iConfig.getParameter<edm::InputTag> ( "tagJetLabel"      ) ),
  corIC5CaloJetsWithBTagLabel_ ( iConfig.getParameter<std::string> ( "corIC5CaloJetsWithBTagLabel" ) ),
  corIC5PFJetsWithBTagFlag_    ( iConfig.getParameter<bool>        ( "corIC5PFJetsWithBTagFlag"    ) ),
  genParticleLabel_  ( iConfig.getParameter<edm::InputTag> ( "genParticleLabel" ) ),
  genJetLabel_       ( iConfig.getParameter<edm::InputTag> ( "genJetLabel"      ) ),
  genMetLabel_       ( iConfig.getParameter<edm::InputTag> ( "genMetLabel"      ) ) {

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
  std::cout << "[SimpleNtple::analyze] InitObjs DONE" << std::endl;
  
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

    
  std::cout << "[SimpleNtple::analyze] DONE" << std::endl;
}


// --------------------------------------------------------------------

void SimpleNtple::InitObjs() {

  // event obj
  evtID_             = 0;
  evtRun_            = 0;
  evtEvent_          = 0;
  jetN_              = 0;
  btagjetN_          = 0;
  eleN_              = 0;
  muN_               = 0;
  glbmuN_            = 0;
  glbmuPromptTightN_ = 0;
  invmasstagjetN_    = 0;
  deltaetatagjetN_   = 0;
  zeptagjetN_        = 0;
  invmasstagjetInvMass_    -> clear ();    // depends on the tag jet definition
  invmasstagjetDeltaEta_   -> clear ();
  invmasstagjetZeppenfeld_ -> clear ();
  deltaetatagjetInvMass_    -> clear ();    // depends on the tag jet definition
  deltaetatagjetDeltaEta_   -> clear ();
  deltaetatagjetZeppenfeld_ -> clear ();
  zeptagjetInvMass_    -> clear ();    // depends on the tag jet definition
  zeptagjetDeltaEta_   -> clear ();
  zeptagjetZeppenfeld_ -> clear ();
  zjetInvMass_      -> clear ();    // depends on the z jet definition => btagger?
  zjetDeltaEta_     -> clear ();
  zjetZeppenfeld_   -> clear ();
}

// --------------------------------------------------------------------
void SimpleNtple::FillEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  std::cout << "SimpleNtple::FillEvent" << std::endl;

  if ( whichSim_ == vbfhzz2l2b::FULLSIM ) {
    edm::Handle<edm::HepMCProduct> evtMC;
    try {
      iEvent.getByLabel("source", evtMC); }
    catch(...) {
      std::cerr << "[SimpleNtple::FillKindEvent] defined as FullSim, but HepMCProduct::source not found" << std::endl; }
    const HepMC::GenEvent * mcEv = evtMC->GetEvent();
    evtID_ = mcEv->signal_process_id();
  }
  else if ( whichSim_ == vbfhzz2l2b::FASTSIM ) {
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
}

// --------------------------------------------------------------------
void SimpleNtple::FillcorIC5PFJetsWithBTag(const edm::Event& iEvent, const edm::EventSetup& iSetup){
}

// --------------------------------------------------------------------
void SimpleNtple::FillMuon(const edm::Event& iEvent, const edm::EventSetup& iSetup){
}

// --------------------------------------------------------------------
void SimpleNtple::FillElectron(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
}

// --------------------------------------------------------------------
void SimpleNtple::FillMet(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  
}


// --------------------------------------------------------------------
void SimpleNtple::FillTagJet(const edm::Event& iEvent, const edm::EventSetup& iSetup){
}

// --------------------------------------------------------------------
void  SimpleNtple::FillTrack(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
}

// --------------------------------------------------------------------
void SimpleNtple::FillGenParticle(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
}




// --------------------------------------------------------------------


void SimpleNtple::FillGenJet(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

}



// --------------------------------------------------------------------
void SimpleNtple::FillGenMet(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
}


// ------------ method called once each job just before starting event loop  ------------
void SimpleNtple::beginJob(const edm::EventSetup& iSetup)
{

  std::cout << "[SimpleNtple::beginJob]" << std::endl;

  invmasstagjetInvMass_    = new std::vector<double>;    // depends on the tag jet definition
  invmasstagjetDeltaEta_   = new std::vector<double>;
  invmasstagjetZeppenfeld_ = new std::vector<double>;
  deltaetatagjetInvMass_    = new std::vector<double>;    // depends on the tag jet definition
  deltaetatagjetDeltaEta_   = new std::vector<double>;
  deltaetatagjetZeppenfeld_ = new std::vector<double>;
  zeptagjetInvMass_    = new std::vector<double>;    // depends on the tag jet definition
  zeptagjetDeltaEta_   = new std::vector<double>;
  zeptagjetZeppenfeld_ = new std::vector<double>;
  zjetInvMass_    = new std::vector<double>;    // depends on the z jet definition => btagger?
  zjetDeltaEta_   = new std::vector<double>;
  zjetZeppenfeld_ = new std::vector<double>;

  mytree_->Branch("evtID",&evtID_,"evtID_/I");
  mytree_->Branch("whichSim",&whichSim_,"whichSim_/I");
  mytree_->Branch("jetN",             &jetN_,             "jetN_/I"             );
  mytree_->Branch("btagjetN",         &btagjetN_,         "btagjetN_/I"         );
  mytree_->Branch("invmasstagjetN",   &invmasstagjetN_,   "invmasstagjetN_/I"   );
  mytree_->Branch("deltaetatagjetN",  &deltaetatagjetN_,  "deltaetatagjetN_/I"  );
  mytree_->Branch("zeptagjetN",       &zeptagjetN_,       "zeptagjetN_/I"       );
  mytree_->Branch("muN",              &muN_,              "muN_/I"              );
  mytree_->Branch("glbmuN",           &glbmuN_,           "glbmumuN_/I"         );
  mytree_->Branch("glbmuPromptTightN",&glbmuPromptTightN_,"glbmuPromptTightN_/I");
  mytree_->Branch("eleN",             &eleN_,             "&eleN_/I"            );
  mytree_->Branch("invmasstagjetInvMass",   "std::vector<double>",&invmasstagjetInvMass_   );
  mytree_->Branch("invmasstagjetDeltaEta",  "std::vector<double>",&invmasstagjetDeltaEta_  );
  mytree_->Branch("invmasstagjetZeppenfeld","std::vector<double>",&invmasstagjetZeppenfeld_);
  mytree_->Branch("deltaetatagjetInvMass",   "std::vector<double>",&deltaetatagjetInvMass_   );
  mytree_->Branch("deltaetatagjetDeltaEta",  "std::vector<double>",&deltaetatagjetDeltaEta_  );
  mytree_->Branch("deltaetatagjetZeppenfeld","std::vector<double>",&deltaetatagjetZeppenfeld_);
  mytree_->Branch("zeptagjetInvMass",   "std::vector<double>",&zeptagjetInvMass_   );
  mytree_->Branch("zeptagjetDeltaEta",  "std::vector<double>",&zeptagjetDeltaEta_  );
  mytree_->Branch("zeptagjetZeppenfeld","std::vector<double>",&zeptagjetZeppenfeld_);
  mytree_->Branch("zjetInvMass",     "std::vector<double>",&zjetInvMass_     );
  mytree_->Branch("zjetDeltaEta",    "std::vector<double>",&zjetDeltaEta_    );
  mytree_->Branch("zjetZeppenfeld",  "std::vector<double>",&zjetZeppenfeld_  );

  std::cout << "DONE w/ event branch" << std::endl;
}


// ------------ method called once each job just after ending the event loop  ------------


void 
SimpleNtple::endJob() {

  delete invmasstagjetInvMass_;    // depends on the tag jet definition
  delete invmasstagjetDeltaEta_;
  delete invmasstagjetZeppenfeld_;
  delete deltaetatagjetInvMass_;    // depends on the tag jet definition
  delete deltaetatagjetDeltaEta_;
  delete deltaetatagjetZeppenfeld_;
  delete zeptagjetInvMass_;    // depends on the tag jet definition
  delete zeptagjetDeltaEta_;
  delete zeptagjetZeppenfeld_;
  delete zjetInvMass_;    // depends on the z jet definition => btagger?
  delete zjetDeltaEta_;
  delete zjetZeppenfeld_;

  std::cout << "[SimpleNtple::endJob]" << std::endl;

}


// --------------------------------------------------------------------




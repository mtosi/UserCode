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
//
// Original Author:  Alessio Ghezzi
//         Created:  Tue Jun  5 19:34:31 CEST 2007
// $Id: SimpleNtple.cc,v 1.3 2009/02/23 15:42:58 amassiro Exp $
//
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

using namespace std;
using namespace edm;
using namespace reco;
using namespace vbfhzz2l2b;


enum { FASTSIM = 0,
       FULLSIM = 1
};

SimpleNtple::SimpleNtple(const edm::ParameterSet& iConfig) :
  whichSim_   ( iConfig.getParameter<int>           ( "whichSim"  ) ), // 0:FastSim, 1:FullSim
  TracksTag_  ( iConfig.getParameter<edm::InputTag> ( "TracksTag" ) ),
  EleTag_     ( iConfig.getParameter<edm::InputTag> ( "EleTag"    ) ),
  IsolEleTag_ ( iConfig.getParameter<edm::InputTag> ( "IsolEleTag") ),
  MuTag_      ( iConfig.getParameter<edm::InputTag> ( "MuTag"     ) ),
  IsolMuTag_  ( iConfig.getParameter<edm::InputTag> ( "IsolMuTag" ) ),
  MetTag_     ( iConfig.getParameter<edm::InputTag> ( "MetTag"    ) ),
  TagJetTag_  ( iConfig.getParameter<edm::InputTag> ( "TagJetTag" ) ),
  JetTag_     ( iConfig.getParameter<edm::InputTag> ( "JetTag"    ) ),
  bool_IterativeCone5CaloJetsTag_         (iConfig.getParameter<bool>( "bool_IterativeCone5CaloJetsTag"         ) ),
  bool_IterativeCone5PFJetsTag_           (iConfig.getParameter<bool>( "bool_IterativeCone5PFJetsTag"           ) ),
  bool_corIterativeCone5CaloJetsWithBTag_ (iConfig.getParameter<bool>( "bool_corIterativeCone5CaloJetsWithBTag" ) ),
  bool_corIterativeCone5PFJetsWithBTag_   (iConfig.getParameter<bool>( "bool_corIterativeCone5PFJetsWithBTag"   ) ),
  bool_SisCone5CaloJetsTag_               (iConfig.getParameter<bool>( "bool_SisCone5CaloJetsTag"               ) ),
  bool_SisCone5PFJetsTag_                 (iConfig.getParameter<bool>( "bool_SisCone5PFJetsTag"                 ) ),
  bool_corSisCone5CaloJetsWithBTag_       (iConfig.getParameter<bool>( "bool_corSisCone5CaloJetsWithBTag"       ) ),
  bool_corSisCone5PFJetsWithBTag_         (iConfig.getParameter<bool>( "bool_corSisCone5PFJetsWithBTag"         ) ),

  MCtruthTag_ ( iConfig.getParameter<edm::InputTag> ( "MCtruthTag" ) ),
  genJetTag_  ( iConfig.getParameter<edm::InputTag> ( "genJetTag"  ) ),
  genMetTag_  ( iConfig.getParameter<edm::InputTag> ( "genMetTag"  ) ) {

  //  corJetWithBTag_( iConfig.getParameter<std::string>( "corJetWithBTag" ) )
  if ( bool_IterativeCone5CaloJetsTag_         ) IterativeCone5CaloJetsTag_         = iConfig.getParameter<edm::InputTag>( "IterativeCone5CaloJetsTag"         );
  if ( bool_IterativeCone5PFJetsTag_           ) IterativeCone5PFJetsTag_           = iConfig.getParameter<edm::InputTag>( "IterativeCone5PFJetsTag"           );
  if ( bool_corIterativeCone5CaloJetsWithBTag_ ) corIterativeCone5CaloJetsWithBTag_ = iConfig.getParameter<std::string>  ( "corIterativeCone5CaloJetsWithBTag" );
  if ( bool_corIterativeCone5PFJetsWithBTag_   ) corIterativeCone5PFJetsWithBTag_   = iConfig.getParameter<std::string>  ( "corIterativeCone5PFJetsWithBTag"   );
  if ( bool_SisCone5CaloJetsTag_               ) SisCone5CaloJetsTag_               = iConfig.getParameter<edm::InputTag>( "SisCone5CaloJetsTag"               );
  if ( bool_SisCone5PFJetsTag_                 ) SisCone5PFJetsTag_                 = iConfig.getParameter<edm::InputTag>( "SisCone5PFJetsTag"                 );
  if ( bool_corSisCone5CaloJetsWithBTag_       ) corSisCone5CaloJetsWithBTag_       = iConfig.getParameter<std::string>  ( "corSisCone5CaloJetsWithBTag"       );
  if ( bool_corSisCone5PFJetsWithBTag_         ) corSisCone5PFJetsWithBTag_         = iConfig.getParameter<std::string>  ( "corSisCone5PFJetsWithBTag"         );
 
  std::cout << "SimpleNtple::SimpleNtple" << std::endl;
  //now do what ever initialization is needed
  edm::Service<TFileService> fs ;
  mytree_  = fs->make <TTree>("VBFSimpleTree","VBFSimpleTree"); 
  
}


// --------------------------------------------------------------------


SimpleNtple::~SimpleNtple()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

  delete emFrac_;

  delete emFracWithBTag_;
  delete corEtWithBTag_;
  delete compoSVbTagDiscrWithBTag_;
  delete highEFFbTagDiscrWithBTag_;
  delete highPURbTagDiscrWithBTag_;
  delete discriminatorVecWithBTag_;

  delete m_tagJets ;
  delete m_otherJets ;
  delete m_otherJets_IterativeCone5CaloJets ;
  delete m_otherJets_IterativeCone5PFJets;
  delete m_otherJets_corIterativeCone5CaloJetsWithBTag ;
  delete m_otherJets_corIterativeCone5PFJetsWithBtag;
  delete m_otherJets_SisCone5CaloJets ;
  delete m_otherJets_SisCone5PFJets;
  delete m_otherJets_corSisCone5CaloJetsWithBTag ;
  delete m_otherJets_corSisCone5PFJetsWithBtag;            

  delete m_electrons ;
  delete m_muons ;
  delete m_MET ;
  delete m_tracks ;
  delete m_genParticles ;
  delete m_genJets;
  delete m_genMet;

  std::cout << "SimpleNtple::~SimpleNtple" << std::endl;
  
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
SimpleNtple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  Init();
  FillKindEvent (iEvent, iSetup);

  FillEle (iEvent, iSetup);
  FillMu (iEvent, iSetup);
  FillMet (iEvent, iSetup);
  //   FillTagJet (iEvent, iSetup); //---- AM --- not now!
  FillJet (iEvent, iSetup);
  if ( bool_IterativeCone5CaloJetsTag_         ) FillIC5CaloJets            (iEvent, iSetup);	       
  if ( bool_IterativeCone5PFJetsTag_           ) FillIC5PFJets              (iEvent, iSetup);	       	 
  if ( bool_corIterativeCone5CaloJetsWithBTag_ ) FillcorIC5CaloJetsWithBTag (iEvent, iSetup); 
  if ( bool_corIterativeCone5PFJetsWithBTag_   ) FillcorIC5PFJetsWithBTag   (iEvent, iSetup);   
  if ( bool_SisCone5CaloJetsTag_               ) FillSC5CaloJets            (iEvent, iSetup);	       
  if ( bool_SisCone5PFJetsTag_                 ) FillSC5PFJets              (iEvent, iSetup);     
  if ( bool_corSisCone5CaloJetsWithBTag_       ) FillcorSC5CaloJetsWithBTag (iEvent, iSetup); 
  if ( bool_corSisCone5PFJetsWithBTag_         ) FillcorSC5PFJetsWithBTag   (iEvent, iSetup);
  if ( whichSim_ == FULLSIM ) FillTracks (iEvent, iSetup);
  FillGenParticles (iEvent, iSetup); //---- AM --- to call after FillKindEvent
  FillGenJet (iEvent, iSetup);
  FillGenMet (iEvent, iSetup);
  
  mytree_->Fill();
  
  emFrac_->clear();

  emFracWithBTag_->clear();
  corEtWithBTag_->clear();
  compoSVbTagDiscrWithBTag_->clear();
  highEFFbTagDiscrWithBTag_->clear();
  highPURbTagDiscrWithBTag_->clear();
  discriminatorVecWithBTag_->clear();

  m_tagJets      -> Clear () ;
  m_otherJets    -> Clear () ;  
  m_otherJets_IterativeCone5CaloJets            -> Clear();
  m_otherJets_IterativeCone5PFJets              -> Clear();
  m_otherJets_corIterativeCone5CaloJetsWithBTag -> Clear();
  m_otherJets_corIterativeCone5PFJetsWithBtag   -> Clear();
  m_otherJets_SisCone5CaloJets                  -> Clear();
  m_otherJets_SisCone5PFJets                    -> Clear();
  m_otherJets_corSisCone5CaloJetsWithBTag       -> Clear();
  m_otherJets_corSisCone5PFJetsWithBtag         -> Clear();            
  m_electrons    -> Clear () ;
  m_muons        -> Clear () ;
  m_MET          -> Clear () ;
  m_tracks       -> Clear () ;
  m_genParticles -> Clear () ;
  m_genJets      -> Clear () ;
  m_genMet       -> Clear () ;
  
}


// --------------------------------------------------------------------

void SimpleNtple::FillKindEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  std::cout << "SimpleNtple::FillKindEvent" << std::endl;

  if ( whichSim_ == FULLSIM ) {
    edm::Handle<edm::HepMCProduct> evtMC;
    try {
      iEvent.getByLabel("source", evtMC); }
    catch(...) {
      std::cerr << "[SimpleNtple::FillKindEvent] defined as FullSim, but HepMCProduct::source not found" << std::endl; }
    const HepMC::GenEvent * mcEv = evtMC->GetEvent();
    IdEvent = mcEv->signal_process_id();
  }
  else if ( whichSim_ == FASTSIM ) {
    edm::Handle<int> genProcessID;
    try {
      iEvent.getByLabel( "genEventProcID", genProcessID ); }
    catch(...) {
      std::cerr << "[SimpleNtple::FillKindEvent] defined as FastSim, but genEventProcID not found" << std::endl; }

    IdEvent = *genProcessID;
  }
  else {
    std::cout << "--> WARNING: simulation not specificied!!" << std::endl;
  }


}


// --------------------------------------------------------------------

void SimpleNtple::FillEle(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  std::cout << "SimpleNtple::FillEle" << std::endl;
  edm::Handle<edm::View<reco::PixelMatchGsfElectron> > EleHandle ;
  //  edm::Handle<reco::PixelMatchGsfElectronCollection> EleHandle ;
  iEvent.getByLabel (EleTag_,EleHandle) ;

  if(EleHandle->size() < 30 ) nEle = EleHandle->size();
  else nEle = 30;

  TClonesArray &electrons = *m_electrons;
  int counter = 0;
  for (int i=0; i< nEle; i++)
    {
      setMomentum (myvector, (*EleHandle)[i].p4());
      new (electrons[counter]) TLorentzVector (myvector);
      counter++;
    }
  
}


// --------------------------------------------------------------------


void SimpleNtple::FillMu(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  std::cout << "SimpleNtple::FillMu" << std::endl;
  edm::Handle<reco::MuonCollection> MuHandle;
  iEvent.getByLabel (MuTag_,MuHandle);
  
  if(MuHandle->size() < 30 ) nMu = MuHandle->size(); 
  else nMu = 30;

  TClonesArray &muons = *m_muons;
  int counter = 0;
  for (int i=0; i< nMu; i++)
    {
      setMomentum (myvector, (*MuHandle)[i].p4());
      new (muons[counter]) TLorentzVector (myvector);
      counter++;
    }


}


// --------------------------------------------------------------------


void SimpleNtple::FillMet(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  std::cout << "SimpleNtple::FillMet" << std::endl;
  edm::Handle<reco::CaloMETCollection> metCollectionHandle;
  iEvent.getByLabel (MetTag_ , metCollectionHandle);
  const CaloMETCollection *calometcol = metCollectionHandle.product();
  const CaloMET *calomet = &(calometcol->front());

  TClonesArray &MET = *m_MET;
  setMomentum (myvector, calomet->p4());
  new (MET[0]) TLorentzVector (myvector);
  
}


// --------------------------------------------------------------------


void SimpleNtple::FillTagJet(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  std::cout << "SimpleNtple::FillTagJet" << std::endl;

  edm::Handle<reco::RecoChargedCandidateCollection> jetTagsHandle;
  // edm::Handle<reco::CandidateCollection> jetTagsHandle ;
  iEvent.getByLabel (TagJetTag_, jetTagsHandle) ;
  
  if (jetTagsHandle->size () != 2) return ;
  math::XYZTLorentzVector sumLV = (*jetTagsHandle)[0].p4() + (*jetTagsHandle)[1].p4() ;
  MinvTags = sumLV.M();

  TClonesArray &jetTag = *m_tagJets;
  int counter = 0;
  for (RecoChargedCandidateCollection::const_iterator jet = jetTagsHandle->begin (); 
       jet != jetTagsHandle->end (); 
       ++jet )
    { 
      setMomentum (myvector, jet->p4());
      new (jetTag[counter]) TLorentzVector (myvector);
      counter++;
    }
}


// --------------------------------------------------------------------


void SimpleNtple::FillJet(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  std::cout << "SimpleNtple::FillJet" << std::endl;

  edm::Handle<reco::CaloJetCollection> JetHandle ;
  iEvent.getByLabel (JetTag_, JetHandle) ;

  TClonesArray &jetOther = *m_otherJets;

  FillCaloJet(JetHandle,jetOther);

//  int counter = 0;
//  for (CaloJetCollection::const_iterator jet = JetHandle->begin (); 
//       jet != JetHandle->end (); 
//       ++jet, counter++ ) { 
//      setMomentum (myvector, jet->p4());        
//      new (jetOther[counter]) TLorentzVector (myvector);
//
//      emFrac_->push_back(jet->emEnergyFraction());
//  }

}

void SimpleNtple::FillCaloJet(const edm::Handle<reco::CaloJetCollection>& caloJetHandle, TClonesArray &caloJetClonesArray) {

  std::cout << "SimpleNtple::FillCaloJet" << std::endl;

  int counter = 0;
  for (CaloJetCollection::const_iterator jet = caloJetHandle->begin (); 
       jet != caloJetHandle->end (); 
       ++jet, counter++ ) { 
    setMomentum (myvector, jet->p4());        
    new (caloJetClonesArray[counter]) TLorentzVector (myvector);
    
    emFrac_->push_back(jet->emEnergyFraction());
  }
  
  
}

void SimpleNtple::FillPFJet(const edm::Handle<reco::PFJetCollection>& pfJetHandle, TClonesArray &pfJetClonesArray) {

  std::cout << "SimpleNtple::FillPFJet" << std::endl;

  int counter = 0;
  for (PFJetCollection::const_iterator jet = pfJetHandle->begin (); 
       jet != pfJetHandle->end (); 
       ++jet, counter++ ) { 
    setMomentum (myvector, jet->p4());        
    new (pfJetClonesArray[counter]) TLorentzVector (myvector);
    
  }
  
  
}

void SimpleNtple::FillIC5CaloJets(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  std::cout << "SimpleNtple::FillIC5CaloJets" << std::endl;

  edm::Handle<reco::CaloJetCollection> JetHandle ;
  iEvent.getByLabel (IterativeCone5CaloJetsTag_, JetHandle) ;

  TClonesArray &jetOther = *m_otherJets_IterativeCone5CaloJets;

  FillCaloJet(JetHandle,jetOther);

}            

void SimpleNtple::FillIC5PFJets(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  std::cout << "SimpleNtple::FillIC5PFJets" << std::endl;

  edm::Handle<reco::PFJetCollection> JetHandle ;
  iEvent.getByLabel (IterativeCone5PFJetsTag_, JetHandle) ;

  TClonesArray &jetOther = *m_otherJets_IterativeCone5PFJets;

  FillPFJet(JetHandle,jetOther);

}
// --------------------------------------------------------------------

void SimpleNtple::FillcorIC5CaloJetsWithBTag(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  std::cout << "SimpleNtple::FillcorIC5CaloJetsWithBTag" << std::endl;

  TClonesArray &jetOtherWithBTag = *m_otherJets_corIterativeCone5CaloJetsWithBTag;

  int counter = 0;

  edm::Handle<vbfhzz2l2b::CorJetWithBTagDiscrCollection> corJetWithBTagHandle;
  iEvent.getByLabel(corIterativeCone5CaloJetsWithBTag_,"corJetWithBTagDiscr",corJetWithBTagHandle);
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

    setMomentum (myvector, (*jet)->p4());      
    new (jetOtherWithBTag[counter]) TLorentzVector (myvector);
    std::cout << "jetOtherWithBTag set" << std::endl;

    compoSVbTagDiscrWithBTag_->push_back(discrVec[0]);
    std::cout << "compoSVbTagDiscrWithBTag set" << std::endl;
    highPURbTagDiscrWithBTag_->push_back(discrVec[1]);
    std::cout << "highPURbTagDiscrWithBTag set" << std::endl;
    highEFFbTagDiscrWithBTag_->push_back(discrVec[2]);
    std::cout << "highEFFbTagDiscrWithBTag set" << std::endl;
    discriminatorVecWithBTag_->push_back(discrVec);
    std::cout << "discriminatorVecWithBTag set" << std::endl;
    corEtWithBTag_->push_back(corrEt);
    std::cout << "corEtWithBTag set" << std::endl;
    emFracWithBTag_->push_back(emFrac);
    std::cout << "emFrac set" << std::endl;
    
    counter++;    

  }

}

void SimpleNtple::FillcorIC5PFJetsWithBTag(const edm::Event& iEvent, const edm::EventSetup& iSetup){
}
void SimpleNtple::FillSC5CaloJets(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  std::cout << "SimpleNtple::FillSC5CaloJets" << std::endl;

  edm::Handle<reco::CaloJetCollection> JetHandle ;
  iEvent.getByLabel (SisCone5CaloJetsTag_, JetHandle) ;

  TClonesArray &jetOther = *m_otherJets_SisCone5CaloJets;

  FillCaloJet(JetHandle,jetOther);
}
void SimpleNtple::FillSC5PFJets(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  std::cout << "SimpleNtple::FillIC5PFJets" << std::endl;

  edm::Handle<reco::PFJetCollection> JetHandle ;
  iEvent.getByLabel (SisCone5PFJetsTag_, JetHandle) ;

  TClonesArray &jetOther = *m_otherJets_SisCone5PFJets;

  FillPFJet(JetHandle,jetOther);

}
void SimpleNtple::FillcorSC5CaloJetsWithBTag(const edm::Event& iEvent, const edm::EventSetup& iSetup){
}
void SimpleNtple::FillcorSC5PFJetsWithBTag(const edm::Event& iEvent, const edm::EventSetup& iSetup){
}
// --------------------------------------------------------------------


void  SimpleNtple::FillTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::cout << "SimpleNtple::FillTracks" << std::endl;
  //  edm::Handle<trackCollection> TracksHandle ;
  edm::Handle<reco::TrackCollection> TracksHandle ;
  iEvent.getByLabel (TracksTag_, TracksHandle) ;
  std::cout << "TracksHandle DONE" << std::endl;

  TClonesArray &tracks = *m_tracks;
  int counter = 0;
  int index = 0;
  //  for (trackCollection::const_iterator tkIt = TracksHandle->begin (); 
  for (reco::TrackCollection::const_iterator tkIt = TracksHandle->begin (); 
       tkIt != TracksHandle->end (); 
       ++tkIt, index++ ) 
    { 
      std::cout << "index: " << index << std::endl;
      //      math::XYZVector mom = (*tkIt).innerMomentum () ; 
      math::XYZVector mom = tkIt->innerMomentum () ; 
      myvector.SetPx (mom.x ()) ;
      myvector.SetPy (mom.y ()) ;
      myvector.SetPz (mom.z ()) ;
      new (tracks[counter]) TLorentzVector (myvector);
      counter++;
     }
}


// --------------------------------------------------------------------


void SimpleNtple::FillGenParticles(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::cout << "SimpleNtple::FillGenParticles" << std::endl;
 edm::Handle<reco::GenParticleCollection> genParticlesHandle; 
 iEvent.getByLabel (MCtruthTag_,genParticlesHandle);
 
 TClonesArray &genParticles = *m_genParticles;
 int counter = 0;
 
 if (IdEvent==123 || IdEvent==124){//---- only if VBF
  for (reco::GenParticleCollection::const_iterator genIt = genParticlesHandle->begin (); 
       genIt != genParticlesHandle->end (); 
       ++genIt ) 
  {
   
   Int_t pdg = genIt->pdgId();
   Int_t status = genIt->status();
   Int_t mother1 = 0;
   if (genIt->numberOfMothers()>0) mother1 = genIt->mother(0)->pdgId();
   Int_t mother2 = 0;
   if (genIt->numberOfMothers()>1) mother2 = genIt->mother(1)->pdgId();
   Int_t daughter1 = 0;
   if (genIt->numberOfDaughters()>0) daughter1 = genIt->daughter(0)->pdgId();
   Int_t daughter2 = 0;
   if (genIt->numberOfDaughters()>1) daughter2 = genIt->daughter(1)->pdgId();
   Double_t px = genIt->px();
   Double_t py = genIt->py();
   Double_t pz = genIt->pz();
   Double_t etot = genIt->energy();
   Double_t vx = 0;
   Double_t vy = 0;
   Double_t vz = 0;
   Double_t time = 0;
  
   new (genParticles[counter]) TParticle (pdg, status, mother1, mother2, daughter1, daughter2, px, py, pz, etot, vx, vy, vz, time);
  
    
   counter++;
  }
 }
}




// --------------------------------------------------------------------


void SimpleNtple::FillGenJet(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::cout << "SimpleNtple::FillGenJet" << std::endl;
 edm::Handle< reco::GenJetCollection > genJetsHandle ;
 iEvent.getByLabel( genJetTag_, genJetsHandle ) ;

 TClonesArray &genJets = *m_genJets;
 int counter = 0;
 for (reco::GenJetCollection::const_iterator gJIt = genJetsHandle->begin (); 
      gJIt != genJetsHandle->end (); 
      ++gJIt ) 
 { 
  
  myvector.SetPx ((*gJIt).px ()) ;
  myvector.SetPy ((*gJIt).py ()) ;
  myvector.SetPz ((*gJIt).pz ()) ;
  myvector.SetE ((*gJIt).emEnergy () + (*gJIt).hadEnergy ()) ;
  new (genJets[counter]) TLorentzVector (myvector);
  counter++;
 }
}



// --------------------------------------------------------------------


void SimpleNtple::FillGenMet(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::cout << "SimpleNtple::FillGenMet" << std::endl;
 edm::Handle< reco::GenMETCollection > genMetHandle ;
 iEvent.getByLabel( genMetTag_, genMetHandle ) ;

 TClonesArray &genMets = *m_genMet;
 int counter = 0;
 for (reco::GenMETCollection::const_iterator gMIt = genMetHandle->begin (); 
      gMIt != genMetHandle->end (); 
      ++gMIt ) 
 { 
  myvector.SetPx ((*gMIt).px ()) ;
  myvector.SetPy ((*gMIt).py ()) ;
  myvector.SetPz ((*gMIt).pz ()) ;
  myvector.SetE ((*gMIt).emEnergy () + (*gMIt).hadEnergy ()) ;
  new (genMets[counter]) TLorentzVector (myvector);
  counter++;
 }
}

// --------------------------------------------------------------------


void SimpleNtple::Init(){

  std::cout << "SimpleNtple::Init" << std::endl;
  nEle = 0; 
  nMu = 0;
  for (int i=0;i<30;i++){
    IsolEleSumPt[i]=0;IsolEleNTracks[i]=0;EleId[i]=0;
    IsolMuSumPt[i]=0;IsolMuNTracks[i]=0;
  }

  MinvTags = -1;
}


// ------------ method called once each job just before starting event loop  ------------


void SimpleNtple::beginJob(const edm::EventSetup& iSetup)
{

  std::cout << "SimpleNtple::beginJob" << std::endl;
  mytree_->Branch("IdEvent",&IdEvent,"IdEvent/I");
 
  mytree_->Branch("nEle",&nEle,"nEle/I");
  mytree_->Branch("IsolEleSumPt",IsolEleSumPt,"IsolEleSumPt[30]/F");
  mytree_->Branch("IsolEleNTracks",IsolEleNTracks,"IsolEleNTracks[30]/F");
  mytree_->Branch("EleId",EleId,"EleId[30]/I");

  mytree_->Branch("nMu",&nMu,"nMu/I");
  mytree_->Branch("IsolMuSumPt",IsolMuSumPt,"IsolMuSumPt[30]/F");
  mytree_->Branch("IsolMuNTracks",IsolMuNTracks,"IsolMuNTracks[30]/F");

  mytree_->Branch("MinvTags",&MinvTags,"MinvTags/F");

  // vector with the 2 tag TLorentzVectors
  m_tagJets = new TClonesArray ("TLorentzVector");
  mytree_->Branch ("tagJets", "TClonesArray", &m_tagJets, 256000,0);

  // vector of the TLorentz Vectors of other jets
  m_otherJets = new TClonesArray ("TLorentzVector");
  emFrac_     = new std::vector<double>;
  mytree_->Branch ("otherJets", "TClonesArray", &m_otherJets, 256000,0);
  mytree_->Branch ("jetEMfrac", "std::vector<double>",&emFrac_);

  m_otherJets_IterativeCone5CaloJets = new TClonesArray ("TLorentzVector");
  mytree_->Branch ("otherJets_IterativeCone5CaloJets", "TClonesArray", &m_otherJets_IterativeCone5CaloJets, 256000,0);
  m_otherJets_IterativeCone5PFJets = new TClonesArray ("TLorentzVector");
  mytree_->Branch ("otherJets_IterativeCone5PFJets", "TClonesArray", &m_otherJets_IterativeCone5PFJets, 256000,0);

  // vector of the TLorentz Vectors of other jets with b tag
  m_otherJets_corIterativeCone5CaloJetsWithBTag = new TClonesArray ("TLorentzVector");
  emFracWithBTag_           = new std::vector<double>;
  corEtWithBTag_            = new std::vector<double>;
  compoSVbTagDiscrWithBTag_ = new std::vector<double>;
  highEFFbTagDiscrWithBTag_ = new std::vector<double>;
  highPURbTagDiscrWithBTag_ = new std::vector<double>;
  discriminatorVecWithBTag_ = new std::vector<std::vector<double> >;
  mytree_->Branch ("otherJets_corIterativeCone5CaloJetsWithBTag", "TClonesArray", &m_otherJets_corIterativeCone5CaloJetsWithBTag, 256000,0);
  mytree_->Branch ("jetEMfracWithBTag",       "std::vector<double>",&emFracWithBTag_);
  mytree_->Branch ("jetCorEtWithBTag",        "std::vector<double>",&corEtWithBTag_);
  mytree_->Branch ("compoSVbTagDiscrWithBTag","std::vector<double>",&compoSVbTagDiscrWithBTag_);
  mytree_->Branch ("highEFFbTagDiscrWithBTag","std::vector<double>",&highEFFbTagDiscrWithBTag_);
  mytree_->Branch ("highPURbTagDiscrWithBTag","std::vector<double>",&highPURbTagDiscrWithBTag_);
  mytree_->Branch ("discriminatorVecWithBTag","std::vector<std::vector<double> >",&discriminatorVecWithBTag_);

  m_otherJets_corIterativeCone5PFJetsWithBtag = new TClonesArray ("TLorentzVector");
  mytree_->Branch ("otherJets_corIterativeCone5PFJetsWithBtag", "TClonesArray", &m_otherJets_corIterativeCone5PFJetsWithBtag, 256000,0);
  m_otherJets_SisCone5CaloJets  = new TClonesArray ("TLorentzVector");
  mytree_->Branch ("otherJets_SisCone5CaloJets", "TClonesArray", &m_otherJets_SisCone5CaloJets, 256000,0);
  m_otherJets_SisCone5PFJets = new TClonesArray ("TLorentzVector");
  mytree_->Branch ("otherJets_SisCone5PFJets", "TClonesArray", &m_otherJets_SisCone5PFJets, 256000,0);
  m_otherJets_corSisCone5CaloJetsWithBTag  = new TClonesArray ("TLorentzVector");
  mytree_->Branch ("otherJets_corSisCone5CaloJetsWithBTag", "TClonesArray", &m_otherJets_corSisCone5CaloJetsWithBTag, 256000,0);
  m_otherJets_corSisCone5PFJetsWithBtag = new TClonesArray ("TLorentzVector");            
  mytree_->Branch ("otherJets_corSisCone5PFJetsWithBtag", "TClonesArray", &m_otherJets_corSisCone5PFJetsWithBtag, 256000,0);

  // vector of the TLorentz Vectors of other jets
  m_electrons = new TClonesArray ("TLorentzVector");
  mytree_->Branch ("electrons", "TClonesArray", &m_electrons, 256000,0);

  // vector of the TLorentz Vectors of other jets
  m_muons = new TClonesArray ("TLorentzVector");
  mytree_->Branch ("muons", "TClonesArray", &m_muons, 256000,0);

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

}


// ------------ method called once each job just after ending the event loop  ------------


void 
SimpleNtple::endJob() {
  std::cout << "SimpleNtple::endJob" << std::endl;
}


// --------------------------------------------------------------------


void SimpleNtple::setMomentum (TLorentzVector &myvector, const LorentzVector & mom)
{
  myvector.SetPx (mom.Px());
  myvector.SetPy (mom.Py());
  myvector.SetPz (mom.Pz());
  myvector.SetE (mom.E());
}


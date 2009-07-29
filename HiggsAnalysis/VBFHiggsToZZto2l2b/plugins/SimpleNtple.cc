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
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/PythiaParticleIndex.h"

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
  
  //  std::cout << "[SimpleNtple::SimpleNtple] DONE" << std::endl;
  
}


// --------------------------------------------------------------------


SimpleNtple::~SimpleNtple()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

  delete eleTrkIsoAlgo_;

  //  std::cout << "[SimpleNtple::~SimpleNtple]" << std::endl;
  
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
SimpleNtple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  InitObjs();
  //  std::cout << "[SimpleNtple::analyze] InitObjs DONE" << std::endl;

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

  //  std::cout << "[SimpleNtple::analyze] DONE" << std::endl;
}


// --------------------------------------------------------------------

void SimpleNtple::InitObjs() {

  //  std::cout << "[SimpleNtple::InitObjs]" << std::endl;

  // event obj
  evtID_    = 0;
  evtRun_   = 0;
  evtEvent_ = 0;
  invmasstagjetN_    = 0;
  deltaetatagjetN_   = 0;
  zeptagjetN_        = 0;
  jetN_              = 0;
  btagjetN_          = 0;
  eleN_              = 0;
  muN_               = 0;
  glbmuN_            = 0;
  glbmuPromptTightN_ = 0;

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

  //electrons;
  eleP4_              -> Clear() ;
  eleLongLived_       -> clear() ;
  elePdgID_           -> clear() ;
  eleCaloEnergy_      -> clear() ;
  eleCaloEnergyError_ -> clear() ;
  eleCaloEnergySig_   -> clear() ;
  eleEnergy_          -> clear() ;
  eleHadOverEm_       -> clear() ;
  eleEt_              -> clear() ;
  elePt_              -> clear() ;
  eleIsoVal_          -> clear() ;
  eleIsoSumPt_        -> clear() ;
  eleIsoNtrack_       -> clear() ;
  eleD0_              -> clear() ;
  eleD0Error_         -> clear() ;
  eleDxy_             -> clear() ;
  eleDxyError_        -> clear() ;
  eleDz_              -> clear() ;
  eleDzError_         -> clear() ;
  eleEtaError_        -> clear() ;
  elePhiError_        -> clear() ;
  eleNormChi2_        -> clear() ;
  eleQoverP_          -> clear() ;
  eleQoverPError_     -> clear() ;
  eleCharge_          -> clear() ;
  elePrimVtxP3_       -> Clear() ;

  //muons
  glbmuPromptTightFlag_ -> clear() ;
  glbmuP4_              -> Clear() ;
  glbmuPrimVtxP3_       -> Clear() ;
  glbmuCharge_          -> clear() ;
  glbmuPdgID_           -> clear() ;
  glbmuEmEnergy_        -> clear() ;   
  glbmuEmS9Energy_      -> clear() ; 
  glbmuHadEnergy_       -> clear() ;  
  glbmuHadS9Energy_     -> clear() ;
  glbmuHoEnergy_        -> clear() ;   
  glbmuHoS9Energy_      -> clear() ; 
  glbmuIso03emEt_       -> clear() ;
  glbmuIso03hadEt_      -> clear() ;
  glbmuIso03hoEt_       -> clear() ;
  glbmuIso03nJets_      -> clear() ;
  glbmuIso03nTracks_    -> clear() ;
  glbmuIso03sumPt_      -> clear() ;
  glbmuIso05emEt_       -> clear() ;
  glbmuIso05hadEt_      -> clear() ;
  glbmuIso05hoEt_       -> clear() ;
  glbmuIso05nJets_      -> clear() ;
  glbmuIso05nTracks_    -> clear() ;
  glbmuIso05sumPt_      -> clear() ;
  glbmuChi2_            -> clear() ;
  glbmuNdof_            -> clear() ;
  glbmud0_              -> clear() ;
  glbmud0Err_           -> clear() ;
  glbmudz_              -> clear() ;
  glbmudzErr_           -> clear() ;

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

  // other jets with b tag
  btagjetNtrack_           -> clear () ;
  btagjetP4_               -> Clear () ;
  btagjetEmEnergyFraction_ -> clear () ;
  btagjetChFrac_           -> clear () ;
  btagjetCorEt_            -> clear () ;
  btagjetCorPt_            -> clear () ;
  btagjetCompoSVbTagDiscr_ -> clear () ;
  btagjetHighEFFbTagDiscr_ -> clear () ;
  btagjetHighPURbTagDiscr_ -> clear () ; 
  btagjetPrimVtxP3_        -> Clear () ; 
  btagjetSecVtxP3_         -> Clear () ;
  btagjetIP2d_             -> clear () ; 
  btagjetIP2dErr_	   -> clear () ;
  btagjetIP2dSig_	   -> clear () ;
  btagjetIP3d_	 	   -> clear () ;
  btagjetIP3dErr_	   -> clear () ;
  btagjetIP3dSig_          -> clear () ;
  btagjetEtaetaMoment_     -> clear () ;
  btagjetPhiphiMoment_     -> clear () ;
  btagjetEtaphiMoment_     -> clear () ;
  btagjetMaxEInEmTowers_   -> clear () ;	
  btagjetMaxEInHadTowers_  -> clear () ;

  // met
  metP4_        -> Clear ();
  metPrimVtxP3_ -> Clear();
  metSig_       -> clear ();
  metSumEt_     -> clear ();

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

// --------------------------------------------------------------------
void SimpleNtple::FillEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  //  std::cout << "[SimpleNtple::FillEvent]" << std::endl;

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

  //  std::cout << "[SimpleNtple::FillcorIC5CaloJetsWithBTag]" << std::endl;

  edm::Handle<vbfhzz2l2b::CorJetWithBTagDiscrCollection> corJetWithBTagHandle;
  iEvent.getByLabel(corIC5CaloJetsWithBTagLabel_,"corJetWithBTagDiscr",corJetWithBTagHandle);
  //  std::cout << "corJetWithBTagHandle->size(): " << corJetWithBTagHandle->size() << std::endl;

  edm::Handle<reco::SecondaryVertexTagInfoCollection> secondaryVtxTagInfosHandle;
  iEvent.getByLabel( "secondaryVertexTagInfos", secondaryVtxTagInfosHandle );  


  
  const TrackIPTagInfoRef & trkIPtagInfoRef = secondaryVtxTagInfosHandle -> at(0) . trackIPTagInfoRef ();
  //  reco::TrackIPTagInfoCollection::const_iterator trackIPTagInfos = 

  std::cout << "****************************************************************" << std::endl;
  std::cout << "************************ primary vertex ************************" << std::endl;
  edm::Ref< reco::VertexCollection > primVtxColl = (trkIPtagInfoRef.product()) -> at(0) . primaryVertex ();
  std::cout << "primVtxColl.size(): " << (primVtxColl.product())->size() << std::endl;

  double chi2 = primVtxColl->chi2(); // chi-squares
  bool   isFake = primVtxColl->isFake(); //	Tells whether a Vertex is fake, i.e.
  bool   isValid = primVtxColl->isValid(); // Tells whether the vertex is valid.
  double ndof = primVtxColl->ndof(); // Number of degrees of freedom
  // Meant to be Double32_t for soft-assignment fitters:
  // tracks may contribute to the vertex with fractional weights. 
  double normalizedChi2 = primVtxColl->normalizedChi2(); //	chi-squared divided by n.d.o.f.
  const reco::Vertex::Point & position = primVtxColl->position (); // position
  size_t tracksSize = primVtxColl->tracksSize(); //	number of tracks
  double pvx = primVtxColl->x(); // x coordinate
  double pvy = primVtxColl->y(); // y coordinate
  double pvz = primVtxColl->z(); // z coordinate
  double pvxError = primVtxColl->xError(); //	error on x
  double pvyError = primVtxColl->yError(); //	error on y
  double pvzError = primVtxColl->zError(); //	error on z
  
  std::cout << "x: " << pvx << "+o-" << pvxError << std::endl;
  std::cout << "y: " << pvy << "+o-" << pvyError << std::endl;
  std::cout << "z: " << pvz << "+o-" << pvzError << std::endl;
  std::cout << "isFake: " << isFake << " <--> isValid: " << isValid << std::endl;
  std::cout << "chi2: " << chi2 << " ndof: " << ndof << " => normalizedChi2: " << normalizedChi2 << "(" << chi2/ndof << ")" << std::endl;
  std::cout << "tracksSize: " << tracksSize << std::endl;
  std::cout << "****************************************************************" << std::endl;
    
      
  
  TClonesArray &jetP4        = *btagjetP4_;
  TClonesArray &jetPrimVtxP3 = *btagjetPrimVtxP3_;
  TClonesArray &jetSecVtxP3  = *btagjetSecVtxP3_;
  std::vector<double> ip2dValVec;
  std::vector<double> ip2dErrVec;
  std::vector<double> ip2dSigVec;
  std::vector<double> ip3dValVec;
  std::vector<double> ip3dErrVec;
  std::vector<double> ip3dSigVec;
  int jetIndex = 0;
  std::vector<reco::JetBaseRef> jets = vbfhzz2l2b::CorJetBTagDiscrAssociation::allJets(*corJetWithBTagHandle);  
  reco::SecondaryVertexTagInfoCollection::const_iterator secondaryVtxTagInfos_itr = secondaryVtxTagInfosHandle->begin(); 
  std::cout << "///////////// loop over jet collection ///////////" << std::endl;
  for ( std::vector<reco::JetBaseRef>::const_iterator jet = jets.begin();
	jet != jets.end(); ++jet, jetIndex++, ++secondaryVtxTagInfos_itr ) {
    std::cout << "   jetIndex: " << jetIndex << std::endl;

    //    RefToBase<Jet> sVjetRef = secondaryVtxTagInfos_itr -> jet();
    //    const JetTracksAssociationRef & jetTrackAssRef = secondaryVtxTagInfos_itr -> jtaRef();
    unsigned int   nSelTrks   = secondaryVtxTagInfos_itr -> nSelectedTracks ();
    TrackRefVector selTrkColl = secondaryVtxTagInfos_itr ->trackIPTagInfoRef() -> selectedTracks ();
  
    //    TrackRefVector trkColl = secondaryVtxTagInfos_itr -> tracks (); // returns a list of tracks associated to the jet [belongs to BaseTagInfo]
    //    std::cout << "trkColl.size():    " << trkColl.size() << std::endl;
    std::cout << "   selTrkColl.size(): " << selTrkColl.size() << std::endl;

    std::cout << "   /-/-/-/-/-/-/-/ loop over selTrkColl /-/-/-/-/-/-/-/" << std::endl;
    for ( unsigned int index = 0; index != selTrkColl.size(); index++ ) { // of course the maximum number of tracks [for which one gets IP data] is the amount of selected tracks
      //      TrackRef selTrkRef = secondaryVtxTagInfos_itr -> track (index);
      //      const reco::SecondaryVertexTagInfo::TrackData selTrkDataByRef  = secondaryVtxTagInfos_itr -> trackData (selTrkRef);
      //      const reco::TrackIPTagInfo::TrackIPData selTrkIPDataByRef = secondaryVtxTagInfos_itr -> trackIPData (selTrkRef);

      const reco::SecondaryVertexTagInfo::TrackData selTrkData  = secondaryVtxTagInfos_itr -> trackData (index); // what can I do w/ this?!?!?
      const reco::TrackIPTagInfo::TrackIPData selTrkIPData = secondaryVtxTagInfos_itr -> trackIPData (index);

      std::cout << "      --- nSelTrks index: " << index << std::endl;
      std::cout << "      selTrkData.usedForVertexFit():     " << selTrkData.usedForVertexFit()     << std::endl;
      std::cout << "      selTrkData.associatedToVertex():   " << selTrkData.associatedToVertex()   << std::endl; // associatedToVertex() referers to secondary vertex
      if (selTrkData.associatedToVertex()) {
	std::cout << "        nVtxCand: " << secondaryVtxTagInfos_itr -> nVertexCandidates () << std::endl;
	std::cout << "        nVtx:     " << secondaryVtxTagInfos_itr -> nVertices ()         << std::endl;
      }
      std::cout << "      selTrkIPData.distanceToFirstTrack: " << selTrkIPData.distanceToFirstTrack << std::endl;
      std::cout << "      selTrkIPData.distanceToJetAxis:    " << selTrkIPData.distanceToJetAxis    << std::endl;
      std::cout << "      selTrkIPData.ip2d.significance():  " << selTrkIPData.ip2d.significance()  << std::endl;   
      std::cout << "      selTrkIPData.ip3d.significance():  " << selTrkIPData.ip3d.significance()  << std::endl;
      ip2dValVec.push_back ( selTrkIPData.ip2d.value()        );
      ip2dErrVec.push_back ( selTrkIPData.ip2d.error()        );
      ip2dSigVec.push_back ( selTrkIPData.ip2d.significance() );
      ip3dValVec.push_back ( selTrkIPData.ip3d.value()        );
      ip3dErrVec.push_back ( selTrkIPData.ip3d.error()        );
      ip3dSigVec.push_back ( selTrkIPData.ip3d.significance() );
      
    }
    std::cout << "   ----------------------------------------------------------------" << std::endl;

    unsigned int trackIPTagInfosIndex = 0;
    const TrackIPTagInfoRef & trkIPtagInfoRef = secondaryVtxTagInfos_itr -> trackIPTagInfoRef ();
    std::cout << "   (trkIPtagInfoRef.product())->size(): " << (trkIPtagInfoRef.product())->size() << std::endl;
    std::cout << "   /---/---/---/---/ loop over trkIPtagInfos /---/---/---/---/" << std::endl;
    for ( reco::TrackIPTagInfoCollection::const_iterator trackIPTagInfos_itr = (trkIPtagInfoRef.product())->begin(); 
	 trackIPTagInfos_itr !=  (trkIPtagInfoRef.product())->end(); ++trackIPTagInfos_itr, trackIPTagInfosIndex++ ) {
      std::cout << "     ---> trackIPTagInfosIndex: " << trackIPTagInfosIndex << std::endl;

      // look @ TrackIPData
      const std::vector<reco::TrackIPTagInfo::TrackIPData> & ipColl = trackIPTagInfos_itr -> impactParameterData (); // vectors of TrackIPData orderd as the selectedTracks()

      int selTrkNum    = trackIPTagInfos_itr->selectedTracks().size();
      int selTrkNum_S1 = 0;
      int selTrkNum_S2 = 0;
      int selTrkNum_S3 = 0;
      
      double selTrkSumPt    = 0.;
      double selTrkSumPt_S1 = 0.;
      double selTrkSumPt_S2 = 0.;
      double selTrkSumPt_S3 = 0.;

      // additional vectors to store subgroups of selectedTracks
      reco::TrackRefVector selTrkColl_S1; // track collection w/ significance > 1.
      reco::TrackRefVector selTrkColl_S2; // track collection w/ significance > 2.
      reco::TrackRefVector selTrkColl_S3; // track collection w/ significance > 3.

      std::cout << "      ipColl.size():     " << ipColl.size()     << std::endl;
      std::cout << "      selTrkColl.size(): " << selTrkColl.size() << std::endl;
    
      unsigned int trkIndex = 0;
      std::cout << "      /-*-/-*-/-*-/ loop over selTrkColl /-*-/-*-/-*-/" << std::endl;
      reco::TrackRefVector::const_iterator selTrkColl_itr = selTrkColl.begin();
      for ( reco::TrackRefVector::const_iterator selTrkColl_itr = selTrkColl.begin();
	    selTrkColl_itr != selTrkColl.end();
	    ++selTrkColl_itr, trkIndex++ ) {

	const reco::TrackIPTagInfo::TrackIPData selTrkIP = secondaryVtxTagInfos_itr -> trackIPData (trkIndex);

	++selTrkNum;
	selTrkSumPt += (*selTrkColl_itr)->pt();

	std::cout << "         ----- trkIndex: " << trkIndex << std::endl;
	std::cout << "         selTrkColl_itr -> charge ():    " << (*selTrkColl_itr) -> charge()   << std::endl;
	std::cout << "         selTrkIP.ip2d.significance:     " << selTrkIP.ip2d.significance()    << std::endl;
	std::cout << "         selTrkIP.ip3d.significance:     " << selTrkIP.ip3d.significance()    << std::endl;
	//	std::cout << "         ipColl_itr->ip2d.significance:  " << ipColl_itr->ip2d.significance() << std::endl;
	//	std::cout << "         ipColl_itr->ip3d.significance:  " << ipColl_itr->ip3d.significance() << std::endl;

	// take tracks with minimum IP significance
	if ( fabs(selTrkIP.ip3d.significance()) > 1. ) {
	  //	if ( fabs(ipColl_itr->ip3d.significance()) > 1. ) {
	  selTrkColl_S1.push_back(*selTrkColl_itr);
	  selTrkSumPt_S1 += (*selTrkColl_itr)->pt();
	  ++selTrkNum_S1;
	}
	if ( fabs(selTrkIP.ip3d.significance()) > 2. ) {
	  //	if ( fabs(ipColl_itr->ip3d.significance()) > 2. ) {
	  selTrkColl_S2.push_back(*selTrkColl_itr);
	  selTrkSumPt_S2 += (*selTrkColl_itr)->pt();
	  ++selTrkNum_S2;
	}
	if ( fabs(selTrkIP.ip3d.significance()) > 3. ) {
	  //	if ( fabs(ipColl_itr->ip3d.significance()) > 3. ) {
	  selTrkColl_S3.push_back(*selTrkColl_itr);
	  selTrkSumPt_S3 += (*selTrkColl_itr)->pt();
	  ++selTrkNum_S3;
	}
      }
      std::cout << "      /-*-/-*-/-*-/-*-/-*-/-*-/-*-/-*-/-*-/-*-/-*-/-*-/" << std::endl;

      // evaluate tag tracks invariant mass
      double bTagTkInvMass    = vbfhzz2l2b::tracksInvariantMass( selTrkColl    );
      double bTagTkInvMass_S1 = vbfhzz2l2b::tracksInvariantMass( selTrkColl_S1 );
      double bTagTkInvMass_S2 = vbfhzz2l2b::tracksInvariantMass( selTrkColl_S2 );
      double bTagTkInvMass_S3 = vbfhzz2l2b::tracksInvariantMass( selTrkColl_S3 );
      
      std::cout << "      bTagTkInvMass   : " << bTagTkInvMass   
		<< "      bTagTkInvMass_S1: " << bTagTkInvMass_S1
		<< "      bTagTkInvMass_S2: " << bTagTkInvMass_S2
		<< "      bTagTkInvMass_S3: " << bTagTkInvMass_S3 << std::endl;
      
      std::cout << "      selTrkSumPt   : " << selTrkSumPt    
		<< "      selTrkSumPt_S1: " << selTrkSumPt_S1 
		<< "      selTrkSumPt_S2: " << selTrkSumPt_S2 
		<< "      selTrkSumPt_S3: " << selTrkSumPt_S3 << std::endl;
      
    }
    std::cout << "   /---/---/---/---/---/---/---/---/---/---/---/---/" << std::endl;
    
    
  
    unsigned int nVtxCand = secondaryVtxTagInfos_itr -> nVertexCandidates ();
    unsigned int nVtx     = secondaryVtxTagInfos_itr -> nVertices ();
    std::cout << "   nVtxCand: " << nVtxCand << std::endl;
    std::cout << "   nVtx:     " << nVtx << std::endl;
    for ( unsigned int secVtxCandIndex = 0; secVtxCandIndex != nVtxCand; secVtxCandIndex++ ) {
      std::cout << "      ****************************************************************" << std::endl;
      std::cout << "      *********************** secondary vertex ***********************" << std::endl;
      
      unsigned int nVtxTrks = secondaryVtxTagInfos_itr -> nVertexTracks (secVtxCandIndex);
      std::cout << "      nVtxTrks: " << nVtxTrks << std::endl;

      const Vertex   secVtx = secondaryVtxTagInfos_itr -> secondaryVertex (secVtxCandIndex);
      double x = secVtx.x(); // x coordinate
      double y = secVtx.y(); // y coordinate
      double z = secVtx.z(); // z coordinate
      double xError = secVtx.xError(); //	error on x
      double yError = secVtx.yError(); //	error on y
      double zError = secVtx.zError(); //	error on z
      std::cout << "      x: " << x << "+o-" << xError << std::endl;
      std::cout << "      y: " << y << "+o-" << yError << std::endl;
      std::cout << "      z: " << z << "+o-" << zError << std::endl;
      reco::TrackRefVector vtxTrkColl = secondaryVtxTagInfos_itr -> vertexTracks (secVtxCandIndex);
      std::cout << "      vtxTrks.size(): " << vtxTrkColl.size() << std::endl;
      unsigned int trkIndex = secondaryVtxTagInfos_itr -> findTrack (vtxTrkColl[0]);
      std::cout << "      trkIndex: " << trkIndex << " <--> 0" << std::endl;
      trkIndex = secondaryVtxTagInfos_itr -> findTrack (vtxTrkColl[1]);
      std::cout << "      trkIndex: " << trkIndex << " <--> 1" << std::endl;
      std::cout << "      ****************************************************************" << std::endl;
    }
 
    std::cout << "   ****************************************************************" << std::endl;
    std::cout << "   ********************** b tag dicriminator **********************" << std::endl;
  
    std::vector<double> discrVec = (*corJetWithBTagHandle)[*jet].discrVec_;
//    std::cout << "corJetWithBTag highEffDiscr: "       << discrVec[vbfhzz2l2b::HIGHEFF]    << std::endl;
//    std::cout << "corJetWithBTag highPurDiscr: "       << discrVec[vbfhzz2l2b::HIGHPUR]    << std::endl;
//    std::cout << "corJetWithBTag combSecVtxDiscr: "    << discrVec[vbfhzz2l2b::COMBSECVTX] << std::endl;
    double corrEt = (*corJetWithBTagHandle)[*jet].corEt_;
    double uncorrEt = (*jet)->et();
    double uncorrPt = (*jet)->pt();
    double corPt    = (corrEt/uncorrEt)*uncorrPt;
    double etaeta = (*jet)->etaetaMoment ();
    double phiphi = (*jet)->phiphiMoment ();
    double etaphi = (*jet)->etaphiMoment ();
    double emFrac       = (dynamic_cast<const reco::CaloJet*>(&**jet))->emEnergyFraction();
    double maxEinEmTow  = (dynamic_cast<const reco::CaloJet*>(&**jet))->maxEInEmTowers();
    double maxEinHadTow = (dynamic_cast<const reco::CaloJet*>(&**jet))->maxEInHadTowers();

    vbfhzz2l2b::setMomentum (myvector_, (*jet)->p4());
    vbfhzz2l2b::setVertex   (myvertex_, pvx, pvy, pvz);
    new (jetP4[jetIndex])        TLorentzVector (myvector_);
    new (jetPrimVtxP3[jetIndex]) TVector3       (myvertex_);
    if ( secondaryVtxTagInfos_itr -> nVertexCandidates() )
      vbfhzz2l2b::setVertex (myvertex_, secondaryVtxTagInfos_itr -> secondaryVertex (0) );
    else 
      vbfhzz2l2b::setVertex (myvertex_, 0.,0.,0.);
    new (jetSecVtxP3[jetIndex])  TVector3 (myvertex_);
    btagjetNtrack_            -> push_back (nSelTrks);
    btagjetEmEnergyFraction_  -> push_back (emFrac);
    btagjetChFrac_            -> push_back (corPt/corrEt);
    btagjetCorEt_             -> push_back (corrEt);
    btagjetCorPt_	      -> push_back (corPt);
    btagjetHighPURbTagDiscr_  -> push_back (discrVec[vbfhzz2l2b::HIGHEFF]   );
    btagjetHighEFFbTagDiscr_  -> push_back (discrVec[vbfhzz2l2b::HIGHPUR]   );
    btagjetCompoSVbTagDiscr_  -> push_back (discrVec[vbfhzz2l2b::COMBSECVTX]);
    btagjetIP2d_              -> push_back (ip2dValVec);
    btagjetIP2dErr_	      -> push_back (ip2dErrVec);
    btagjetIP2dSig_	      -> push_back (ip2dSigVec);
    btagjetIP3d_	      -> push_back (ip3dValVec);
    btagjetIP3dErr_	      -> push_back (ip3dErrVec);
    btagjetIP3dSig_           -> push_back (ip3dSigVec);
    btagjetEtaetaMoment_      -> push_back (etaeta);
    btagjetPhiphiMoment_      -> push_back (phiphi);
    btagjetEtaphiMoment_      -> push_back (etaphi);
    btagjetMaxEInEmTowers_    -> push_back (maxEinEmTow );	
    btagjetMaxEInHadTowers_   -> push_back (maxEinHadTow);
  }
  std::cout << "///////////////////////////////////////////////////////////////////////" << std::endl;
  btagjetN_ = jets.size();
}

// --------------------------------------------------------------------
void SimpleNtple::FillcorIC5PFJetsWithBTag(const edm::Event& iEvent, const edm::EventSetup& iSetup){
}

// --------------------------------------------------------------------
void SimpleNtple::FillMuon(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  //  std::cout << "[SimpleNtple::FillMuon]" << std::endl;

  edm::Handle<reco::MuonCollection> muonHandle;
  iEvent.getByLabel (muonLabel_,muonHandle);


  TClonesArray &muonP4    = *glbmuP4_;
  TClonesArray &muonPrimVtxP3 = *glbmuPrimVtxP3_;

  int muonIndex = 0;
  int globalmuonCounter = 0;
  int globalmuonprompttightCounter = 0;
  for ( reco::MuonCollection::const_iterator muon_itr = muonHandle->begin();
	muon_itr != muonHandle->end(); ++muon_itr, muonIndex++ ) {

    bool glbmuPromptTightFlag = false;
    if ( !muon_itr->isGlobalMuon() ) continue;
    globalmuonCounter++;
    if ( muon_itr->isGood(reco::Muon::GlobalMuonPromptTight) ) {
      globalmuonprompttightCounter++;
      glbmuPromptTightFlag = true;
    }

    glbmuPromptTightFlag_ -> push_back(glbmuPromptTightFlag);
    vbfhzz2l2b::setMomentum (myvector_, muon_itr->p4()    );
    vbfhzz2l2b::setVertex   (myvertex_, muon_itr->vertex());
    new (muonP4[muonIndex])        TLorentzVector (myvector_);
    new (muonPrimVtxP3[muonIndex]) TVector3       (myvertex_);

    glbmuCharge_       -> push_back (muon_itr -> charge()          );
    glbmuPdgID_        -> push_back (muon_itr -> pdgId()           );
    glbmuEmEnergy_     -> push_back (muon_itr -> calEnergy().em    );   
    glbmuEmS9Energy_   -> push_back (muon_itr -> calEnergy().emS9  ); 
    glbmuHadEnergy_    -> push_back (muon_itr -> calEnergy().had   );  
    glbmuHadS9Energy_  -> push_back (muon_itr -> calEnergy().hadS9 );
    glbmuHoEnergy_     -> push_back (muon_itr -> calEnergy().ho    );   
    glbmuHoS9Energy_   -> push_back (muon_itr -> calEnergy().hoS9  ); 
    glbmuIso03emEt_    -> push_back (muon_itr -> isolationR03().emEt    );
    glbmuIso03hadEt_   -> push_back (muon_itr -> isolationR03().hadEt   );
    glbmuIso03hoEt_    -> push_back (muon_itr -> isolationR03().hoEt    );
    glbmuIso03nJets_   -> push_back (muon_itr -> isolationR03().nJets   );
    glbmuIso03nTracks_ -> push_back (muon_itr -> isolationR03().nTracks );
    glbmuIso03sumPt_   -> push_back (muon_itr -> isolationR03().sumPt   );
    glbmuIso05emEt_    -> push_back (muon_itr -> isolationR05().emEt    );
    glbmuIso05hadEt_   -> push_back (muon_itr -> isolationR05().hadEt   );
    glbmuIso05hoEt_    -> push_back (muon_itr -> isolationR05().hoEt    );
    glbmuIso05nJets_   -> push_back (muon_itr -> isolationR05().nJets   );
    glbmuIso05nTracks_ -> push_back (muon_itr -> isolationR05().nTracks );
    glbmuIso05sumPt_   -> push_back (muon_itr -> isolationR05().sumPt   );
    glbmuChi2_         -> push_back (muon_itr -> globalTrack()->chi2()    );
    glbmuNdof_         -> push_back (muon_itr -> globalTrack()->ndof()    );
    glbmud0_           -> push_back (muon_itr -> globalTrack()->d0()      );
    glbmud0Err_        -> push_back (muon_itr -> globalTrack()->d0Error() );
    glbmudz_           -> push_back (muon_itr -> globalTrack()->dz()      );
    glbmudzErr_        -> push_back (muon_itr -> globalTrack()->dzError() );

  }
  muN_               = muonHandle->size(); 
  glbmuN_            = globalmuonCounter;
  glbmuPromptTightN_ = globalmuonprompttightCounter;
}

// --------------------------------------------------------------------
void SimpleNtple::FillElectron(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //  std::cout << "[SimpleNtple::FillElectron]" << std::endl;

  edm::Handle<reco::PixelMatchGsfElectronCollection> electronHandle ;
  iEvent.getByLabel (electronLabel_,electronHandle) ;

  edm::Handle<reco::TrackCollection> trackHandle ;
  iEvent.getByLabel (trackLabel_, trackHandle) ;

  TClonesArray &electronP4 = *eleP4_;
  TClonesArray &electronPrimVtxP3 = *elePrimVtxP3_;
  int electronIndex = 0;
  for ( reco::PixelMatchGsfElectronCollection::const_iterator electron_itr = electronHandle->begin();
	electron_itr != electronHandle->end(); ++electron_itr ) {
    vbfhzz2l2b::setMomentum (myvector_, electron_itr->p4());
    vbfhzz2l2b::setVertex   (myvertex_, electron_itr->vertex());
    new (electronP4[electronIndex])        TLorentzVector (myvector_);
    new (electronPrimVtxP3[electronIndex]) TVector3       (myvertex_);

    eleLongLived_       -> push_back (electron_itr->longLived());
    elePdgID_           -> push_back (electron_itr->pdgId());
    eleCaloEnergy_      -> push_back (electron_itr->caloEnergy());
    eleCaloEnergyError_ -> push_back (electron_itr->caloEnergyError());
    eleCaloEnergySig_   -> push_back (fabs(electron_itr->caloEnergy())/
				      electron_itr->caloEnergyError());
    eleEnergy_          -> push_back (electron_itr->energy());
    eleHadOverEm_       -> push_back (electron_itr->hadronicOverEm());
    eleEt_              -> push_back (electron_itr->et());
    elePt_              -> push_back (electron_itr->pt());
    eleIsoVal_          -> push_back (eleTrkIsoAlgo_->isoValue (electronHandle,trackHandle,*electron_itr));
    eleIsoSumPt_        -> push_back (eleTrkIsoAlgo_->sumPt    (electronHandle,trackHandle,*electron_itr));
    eleIsoNtrack_       -> push_back (eleTrkIsoAlgo_->numOfTrks(electronHandle,trackHandle,*electron_itr));
    eleD0_              -> push_back (electron_itr->bestTrack()->d0());
    eleD0Error_         -> push_back (electron_itr->bestTrack()->d0Error());
    eleD0Sig_           -> push_back (fabs(electron_itr->bestTrack()->d0())/
				      electron_itr->bestTrack()->d0Error());
    eleDxy_             -> push_back (electron_itr->bestTrack()->dxy());
    eleDxyError_        -> push_back (electron_itr->bestTrack()->dxyError());
    eleDxySig_          -> push_back (fabs(electron_itr->bestTrack()->dxy())/
				      electron_itr->bestTrack()->dxyError());
    eleDz_              -> push_back (electron_itr->bestTrack()->dz());
    eleDzError_         -> push_back (electron_itr->bestTrack()->dzError());
    eleDzSig_           -> push_back (fabs(electron_itr->bestTrack()->dz())/
				      electron_itr->bestTrack()->dzError());
    eleEtaError_        -> push_back (electron_itr->bestTrack()->etaError());
    elePhiError_        -> push_back (electron_itr->bestTrack()->phiError());
    eleNormChi2_        -> push_back (electron_itr->bestTrack()->normalizedChi2());
    eleQoverP_          -> push_back (electron_itr->bestTrack()->qoverp());
    eleQoverPError_     -> push_back (electron_itr->bestTrack()->qoverpError());
    eleCharge_          -> push_back (electron_itr->charge());

  }

  eleN_ = electronHandle->size(); 
  
}

// --------------------------------------------------------------------
void SimpleNtple::FillMet(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  //  std::cout << "[SimpleNtple::FillMet]" << std::endl;

  edm::Handle<reco::CaloMETCollection> metCollectionHandle;
  iEvent.getByLabel (metLabel_ , metCollectionHandle);
  const reco::CaloMETCollection *calometcol = metCollectionHandle.product();
  const reco::CaloMET *calomet = &(calometcol->front());

  TClonesArray &metP4        = *metP4_;
  TClonesArray &metPrimVtxP3 = *metPrimVtxP3_;
  vbfhzz2l2b::setMomentum (myvector_, calomet->p4());
  vbfhzz2l2b::setVertex   (myvertex_, calomet->vertex());
  new (metP4[0])        TLorentzVector (myvector_);
  new (metPrimVtxP3[0]) TVector3       (myvertex_);

  metSig_   -> push_back (calomet->metSignificance());
  metSumEt_ -> push_back (calomet->sumEt());
  
}


// --------------------------------------------------------------------
void SimpleNtple::FillTagJet(const edm::Event& iEvent, const edm::EventSetup& iSetup){

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
  invmasstagjetInvMass_->push_back(invMass);
  invmasstagjetDeltaEta_->push_back(deltaEta);
  invmasstagjetZeppenfeld_->push_back(zeppenfeld);
  invmasstagjetN_ = njets;

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
  deltaetatagjetInvMass_->push_back(invMass);
  deltaetatagjetDeltaEta_->push_back(deltaEta);
  deltaetatagjetZeppenfeld_->push_back(zeppenfeld);
  deltaetatagjetN_ = njets;


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

  zeptagjetInvMass_->push_back(invMass);
  zeptagjetDeltaEta_->push_back(deltaEta);
  zeptagjetZeppenfeld_->push_back(zeppenfeld);
  zeptagjetN_ = njets;
}

// --------------------------------------------------------------------
void  SimpleNtple::FillTrack(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //  std::cout << "[SimpleNtple::FillTrack]" << std::endl;

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

  //  std::cout << "[SimpleNtple::FillGenParticle]" << std::endl;

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


void SimpleNtple::FillGenJet(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //  std::cout << "[SimpleNtple::FillGenJet]" << std::endl;

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
void SimpleNtple::FillGenMet(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //  std::cout << "[SimpleNtple::FillGenMet]" << std::endl;

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

  //  std::cout << "[SimpleNtple::FillGenMet] DONE" << std::endl;
}


// ------------ method called once each job just before starting event loop  ------------
void SimpleNtple::beginJob(const edm::EventSetup& iSetup)
{

  //  std::cout << "[SimpleNtple::beginJob]" << std::endl;

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

  // vector of the TLorentz Vectors of electron
  eleLongLived_       = new std::vector<bool>;
  elePdgID_           = new std::vector<int>;
  eleCaloEnergy_      = new std::vector<double>;
  eleCaloEnergyError_ = new std::vector<double>;
  eleCaloEnergySig_   = new std::vector<double>;
  eleEnergy_          = new std::vector<double>;
  eleHadOverEm_       = new std::vector<double>;
  eleEt_              = new std::vector<double>;
  elePt_              = new std::vector<double>;
  eleIsoVal_          = new std::vector<double>;
  eleIsoSumPt_        = new std::vector<double>;
  eleIsoNtrack_       = new std::vector<int>;
  eleD0_              = new std::vector<double>;
  eleD0Error_         = new std::vector<double>;
  eleD0Sig_           = new std::vector<double>;
  eleDxy_             = new std::vector<double>;
  eleDxyError_        = new std::vector<double>;
  eleDxySig_          = new std::vector<double>;
  eleDz_              = new std::vector<double>;
  eleDzError_         = new std::vector<double>;
  eleDzSig_           = new std::vector<double>;
  eleEtaError_        = new std::vector<double>;
  elePhiError_        = new std::vector<double>;
  eleNormChi2_        = new std::vector<double>;
  eleQoverP_          = new std::vector<double>;
  eleQoverPError_     = new std::vector<double>;
  eleCharge_          = new std::vector<int>;
  eleP4_        = new TClonesArray ("TLorentzVector");
  elePrimVtxP3_ = new TClonesArray ("TVector3");
  mytree_->Branch("electronP4",       "TClonesArray",&eleP4_,       256000,0);
  mytree_->Branch("electronPrimVtxP3","TClonesArray",&elePrimVtxP3_,256000,0);
  mytree_->Branch("eleLongLived"       ,"std::vector<bool>",  &eleLongLived_      );
  mytree_->Branch("elePdgID"           ,"std::vector<int>",   &elePdgID_          );
  mytree_->Branch("eleCaloEnergy"      ,"std::vector<double>",&eleCaloEnergy_     );
  mytree_->Branch("eleCaloEnergyError" ,"std::vector<double>",&eleCaloEnergyError_);
  mytree_->Branch("eleCaloEnergySig"   ,"std::vector<double>",&eleCaloEnergySig_  );
  mytree_->Branch("eleEnergy"          ,"std::vector<double>",&eleEnergy_         );
  mytree_->Branch("eleHadOverEm"       ,"std::vector<double>",&eleHadOverEm_      );
  mytree_->Branch("eleEt"              ,"std::vector<double>",&eleEt_             );
  mytree_->Branch("elePt"              ,"std::vector<double>",&elePt_             );
  mytree_->Branch("eleIsoVal"          ,"std::vector<double>",&eleIsoVal_         );         
  mytree_->Branch("eleIsoSumPt"        ,"std::vector<double>",&eleIsoSumPt_       );
  mytree_->Branch("eleIsoNtrack"       ,"std::vector<int>",   &eleIsoNtrack_      );
  mytree_->Branch("eleD0"              ,"std::vector<double>",&eleD0_             );
  mytree_->Branch("eleD0Error"         ,"std::vector<double>",&eleD0Error_        );
  mytree_->Branch("eleD0Sig"           ,"std::vector<double>",&eleD0Sig_          );
  mytree_->Branch("eleDxy"             ,"std::vector<double>",&eleDxy_            );
  mytree_->Branch("eleDxyError"        ,"std::vector<double>",&eleDxyError_       );
  mytree_->Branch("eleDxySig"          ,"std::vector<double>",&eleDxySig_         );
  mytree_->Branch("eleDz"              ,"std::vector<double>",&eleDz_             );
  mytree_->Branch("eleDzError"         ,"std::vector<double>",&eleDzError_        );
  mytree_->Branch("eleDzSig"           ,"std::vector<double>",&eleDzSig_          );
  mytree_->Branch("eleEtaError"        ,"std::vector<double>",&eleEtaError_       );
  mytree_->Branch("elePhiError"        ,"std::vector<double>",&elePhiError_       );
  mytree_->Branch("eleNormChi2"        ,"std::vector<double>",&eleNormChi2_       );
  mytree_->Branch("eleQoverP"          ,"std::vector<double>",&eleQoverP_         );
  mytree_->Branch("eleQoverPError"     ,"std::vector<double>",&eleQoverPError_    );
  mytree_->Branch("eleCharge"          ,"std::vector<int>",   &eleCharge_         );

  // muon
  glbmuPromptTightFlag_ = new std::vector<bool>;
  glbmuP4_              = new TClonesArray ("TLorentzVector");
  glbmuPrimVtxP3_       = new TClonesArray ("TVector3");      
  glbmuCharge_          = new std::vector<int>;
  glbmuPdgID_           = new std::vector<int>;
  glbmuEmEnergy_        = new std::vector<double>;   
  glbmuEmS9Energy_      = new std::vector<double>; 
  glbmuHadEnergy_       = new std::vector<double>;  
  glbmuHadS9Energy_     = new std::vector<double>;
  glbmuHoEnergy_        = new std::vector<double>;   
  glbmuHoS9Energy_      = new std::vector<double>; 
  glbmuIso03emEt_       = new std::vector<double>;
  glbmuIso03hadEt_      = new std::vector<double>;
  glbmuIso03hoEt_       = new std::vector<double>;
  glbmuIso03nJets_      = new std::vector<int>;
  glbmuIso03nTracks_    = new std::vector<int>;
  glbmuIso03sumPt_      = new std::vector<double>;
  glbmuIso05emEt_       = new std::vector<double>;
  glbmuIso05hadEt_      = new std::vector<double>;
  glbmuIso05hoEt_       = new std::vector<double>;
  glbmuIso05nJets_      = new std::vector<int>;
  glbmuIso05nTracks_    = new std::vector<int>;
  glbmuIso05sumPt_      = new std::vector<double>;
  glbmuChi2_            = new std::vector<double>;
  glbmuNdof_            = new std::vector<double>;
  glbmud0_              = new std::vector<double>;
  glbmud0Err_           = new std::vector<double>;
  glbmudz_              = new std::vector<double>;
  glbmudzErr_           = new std::vector<double>;
  mytree_->Branch("glbmuonP4",        "TClonesArray",&glbmuP4_,        256000,0);
  mytree_->Branch("glbmuonPrimVtxP3", "TClonesArray",&glbmuPrimVtxP3_, 256000,0);
  mytree_->Branch("glbmuonPromptTightFlag","std::vector<bool>",  &glbmuPromptTightFlag_);
  mytree_->Branch("glbmuonCharge",         "std::vector<int>",   &glbmuCharge_         );
  mytree_->Branch("glbmuonPdgID",          "std::vector<int>",   &glbmuPdgID_          );
  mytree_->Branch("glbmuonEmEnergy",       "std::vector<double>",&glbmuEmEnergy_       );
  mytree_->Branch("glbmuonEmS9Energy",     "std::vector<double>",&glbmuEmS9Energy_     );
  mytree_->Branch("glbmuonHadEnergy",      "std::vector<double>",&glbmuHadEnergy_      );
  mytree_->Branch("glbmuonHadS9Energy",    "std::vector<double>",&glbmuHadS9Energy_    );
  mytree_->Branch("glbmuonHoEnergy",       "std::vector<double>",&glbmuHoEnergy_       );
  mytree_->Branch("glbmuonHoS9Energy",     "std::vector<double>",&glbmuHoS9Energy_     );
  mytree_->Branch("glbmuonIso03emEt",      "std::vector<double>",&glbmuIso03emEt_      );
  mytree_->Branch("glbmuonIso03hadEt",     "std::vector<double>",&glbmuIso03hadEt_     );
  mytree_->Branch("glbmuonIso03hoEt",      "std::vector<double>",&glbmuIso03hoEt_      );
  mytree_->Branch("glbmuonIso03nJets",     "std::vector<int>",   &glbmuIso03nJets_     );
  mytree_->Branch("glbmuonIso03nTracks",   "std::vector<int>",   &glbmuIso03nTracks_   );
  mytree_->Branch("glbmuonIso03sumPt",     "std::vector<double>",&glbmuIso03sumPt_     );
  mytree_->Branch("glbmuonIso05emEt",      "std::vector<double>",&glbmuIso05emEt_      );
  mytree_->Branch("glbmuonIso05hadEt",     "std::vector<double>",&glbmuIso05hadEt_     );
  mytree_->Branch("glbmuonIso05hoEt",      "std::vector<double>",&glbmuIso05hoEt_      );
  mytree_->Branch("glbmuonIso05nJets",     "std::vector<int>",   &glbmuIso05nJets_     );
  mytree_->Branch("glbmuonIso05nTracks",   "std::vector<int>",   &glbmuIso05nTracks_   );
  mytree_->Branch("glbmuonIso05sumPt",     "std::vector<double>",&glbmuIso05sumPt_     );
  mytree_->Branch("glbmuonChi2",           "std::vector<double>",&glbmuChi2_           );
  mytree_->Branch("glbmuonNdof",           "std::vector<double>",&glbmuNdof_           );
  mytree_->Branch("glbmuond0",             "std::vector<double>",&glbmud0_             );
  mytree_->Branch("glbmuond0Err",          "std::vector<double>",&glbmud0Err_          );
  mytree_->Branch("glbmuondz",             "std::vector<double>",&glbmudz_             );
  mytree_->Branch("glbmuondzErr",          "std::vector<double>",&glbmudzErr_          );

  // vector of the TLorentz Vectors of other jets with b tag
  btagjetP4_        = new TClonesArray ("TLorentzVector");
  btagjetPrimVtxP3_ = new TClonesArray ("TVector3");
  btagjetSecVtxP3_  = new TClonesArray ("TVector3");
  btagjetIP2d_             = new std::vector< std::vector<double> >;
  btagjetIP2dErr_          = new std::vector< std::vector<double> >;
  btagjetIP2dSig_          = new std::vector< std::vector<double> >;
  btagjetIP3d_             = new std::vector< std::vector<double> >;	 
  btagjetIP3dErr_          = new std::vector< std::vector<double> >;
  btagjetIP3dSig_          = new std::vector< std::vector<double> >;
  btagjetNtrack_           = new std::vector<int>;
  btagjetEmEnergyFraction_ = new std::vector<double>;
  btagjetChFrac_           = new std::vector<double>;
  btagjetCorEt_            = new std::vector<double>;
  btagjetCorPt_            = new std::vector<double>;
  btagjetCompoSVbTagDiscr_ = new std::vector<double>;
  btagjetHighEFFbTagDiscr_ = new std::vector<double>;
  btagjetHighPURbTagDiscr_ = new std::vector<double>;
  btagjetEtaetaMoment_     = new std::vector<double>;
  btagjetPhiphiMoment_     = new std::vector<double>;
  btagjetEtaphiMoment_     = new std::vector<double>;
  btagjetMaxEInEmTowers_   = new std::vector<double>;	
  btagjetMaxEInHadTowers_  = new std::vector<double>;
  mytree_->Branch ("btagjetP4",              "TClonesArray",       &btagjetP4_,        256000,0);
  mytree_->Branch ("btagjetPrimVtxP3",       "TClonesArray",       &btagjetPrimVtxP3_, 256000,0);
  mytree_->Branch ("btagjetSecVtxP3",        "TClonesArray",       &btagjetSecVtxP3_,  256000,0);
  mytree_->Branch ("btagjetIP2d",            "std::vector<std::vector<double> >",&btagjetIP2d_);
  mytree_->Branch ("btagjetIP2dErr",         "std::vector<std::vector<double> >",&btagjetIP2dErr_);
  mytree_->Branch ("btagjetIP2dSig",         "std::vector<std::vector<double> >",&btagjetIP2dSig_);
  mytree_->Branch ("btagjetIP3d",            "std::vector<std::vector<double> >",&btagjetIP3d_   );	 
  mytree_->Branch ("btagjetIP3dErr",         "std::vector<std::vector<double> >",&btagjetIP3dErr_);
  mytree_->Branch ("btagjetIP3dSig",         "std::vector<std::vector<double> >",&btagjetIP3dSig_);
  mytree_->Branch ("btagjetNtrack",          "std::vector<int>",   &btagjetNtrack_);
  mytree_->Branch ("btagjetEmFrac",          "std::vector<double>",&btagjetEmEnergyFraction_);
  mytree_->Branch ("btagjetChFrac",          "std::vector<double>",&btagjetChFrac_          );
  mytree_->Branch ("btagjetCorEt",           "std::vector<double>",&btagjetCorEt_           );
  mytree_->Branch ("btagjetCorPt",           "std::vector<double>",&btagjetCorPt_           );
  mytree_->Branch ("btagjetcompoSVbTagDiscr","std::vector<double>",&btagjetCompoSVbTagDiscr_);
  mytree_->Branch ("btagjethighEFFbTagDiscr","std::vector<double>",&btagjetHighEFFbTagDiscr_);
  mytree_->Branch ("btagjethighPURbTagDiscr","std::vector<double>",&btagjetHighPURbTagDiscr_);
  mytree_->Branch("btagjetEtaetaMoment",     "std::vector<double>",&btagjetEtaetaMoment_    );
  mytree_->Branch("btagjetPhiphiMoment",     "std::vector<double>",&btagjetPhiphiMoment_    );
  mytree_->Branch("btagjetEtaphiMoment",     "std::vector<double>",&btagjetEtaphiMoment_    );
  mytree_->Branch("btagjetMaxEInEmTowers",   "std::vector<double>",&btagjetMaxEInEmTowers_  );
  mytree_->Branch("btagjetMaxEInHadTowers",  "std::vector<double>",&btagjetMaxEInHadTowers_ );

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


  // vector of the TLorentz Vectors of met
  metP4_        = new TClonesArray ("TLorentzVector");
  metPrimVtxP3_ = new TClonesArray ("TVector3");
  metSig_   = new std::vector<double>;
  metSumEt_ = new std::vector<double>;
  mytree_->Branch ("metP4",       "TClonesArray", &metP4_,        256000,0);
  mytree_->Branch ("metPrimVtxP3","TClonesArray", &metPrimVtxP3_, 256000,0);
  mytree_->Branch ("metSig",  "std::vector<double>",&metSig_);
  mytree_->Branch ("metSumEt","std::vector<double>",&metSumEt_);


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
SimpleNtple::endJob() {

  //  std::cout << "[SimpleNtple::endJob]" << std::endl;

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

  delete eleP4_ ;
  delete eleLongLived_;
  delete elePdgID_;
  delete eleCaloEnergy_;
  delete eleCaloEnergyError_;
  delete eleCaloEnergySig_;
  delete eleEnergy_;
  delete eleHadOverEm_;
  delete eleEt_;
  delete elePt_;
  delete eleIsoVal_;
  delete eleIsoSumPt_;
  delete eleIsoNtrack_;
  delete eleD0_;
  delete eleD0Error_;
  delete eleD0Sig_;
  delete eleDxy_;
  delete eleDxyError_;
  delete eleDxySig_;
  delete eleDz_;
  delete eleDzError_;
  delete eleDzSig_;
  delete eleEtaError_;
  delete elePhiError_;
  delete eleNormChi2_;
  delete eleQoverP_;
  delete eleQoverPError_;
  delete eleCharge_;
  delete elePrimVtxP3_;

  delete glbmuPromptTightFlag_;
  delete glbmuP4_;
  delete glbmuPrimVtxP3_;
  delete glbmuCharge_;
  delete glbmuPdgID_;
  delete glbmuEmEnergy_;
  delete glbmuEmS9Energy_;
  delete glbmuHadEnergy_;
  delete glbmuHadS9Energy_;
  delete glbmuHoEnergy_;
  delete glbmuHoS9Energy_;
  delete glbmuIso03emEt_;
  delete glbmuIso03hadEt_;
  delete glbmuIso03hoEt_;
  delete glbmuIso03nJets_;
  delete glbmuIso03nTracks_;
  delete glbmuIso03sumPt_;
  delete glbmuIso05emEt_;
  delete glbmuIso05hadEt_;
  delete glbmuIso05hoEt_;
  delete glbmuIso05nJets_;
  delete glbmuIso05nTracks_;
  delete glbmuIso05sumPt_  ;
  delete glbmuChi2_;
  delete glbmuNdof_;
  delete glbmud0_;
  delete glbmud0Err_;
  delete glbmudz_;
  delete glbmudzErr_;

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

  delete btagjetP4_;
  delete btagjetPrimVtxP3_; 
  delete btagjetSecVtxP3_; 
  delete btagjetNtrack_;
  delete btagjetIP2d_;
  delete btagjetIP2dErr_;
  delete btagjetIP2dSig_;
  delete btagjetIP3d_;
  delete btagjetIP3dErr_;
  delete btagjetIP3dSig_;
  delete btagjetEmEnergyFraction_;
  delete btagjetChFrac_;
  delete btagjetCorEt_;
  delete btagjetCorPt_;
  delete btagjetCompoSVbTagDiscr_;
  delete btagjetHighEFFbTagDiscr_;
  delete btagjetHighPURbTagDiscr_;
  delete btagjetEtaetaMoment_;
  delete btagjetPhiphiMoment_;
  delete btagjetEtaphiMoment_;
  delete btagjetMaxEInEmTowers_;	
  delete btagjetMaxEInHadTowers_;

  delete metP4_ ;
  delete metPrimVtxP3_;
  delete metSig_ ;
  delete metSumEt_;

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
  
  //  std::cout << "[SimpleNtple::endJob]" << std::endl;
    
}


// --------------------------------------------------------------------




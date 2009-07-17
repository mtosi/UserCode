// -*- C++ -*-
//
// Package:    VBFHZZllbbBTagInfoAnalyzer
// Class:      VBFHZZllbbBTagInfoAnalyzer
// 
/**\class VBFHZZllbbBTagInfoAnalyzer VBFHZZllbbBTagInfoAnalyzer.cc HiggsAnalysis/VBFHiggsToZZto2l2b/src/VBFHZZllbbBTagInfoAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Mia TOSI
//         Created:  Mon Feb  2 17:31:44 CET 2009
// $Id: VBFHZZllbbBTagInfoAnalyzer.cc,v 1.2 2009/07/06 13:15:48 tosi Exp $
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

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbBTagInfoAnalyzer.h"

#include "PhysicsTools/UtilAlgos/interface/TFileDirectory.h"
#include "PhysicsTools/UtilAlgos/interface/TH1AddDirectorySentry.h" 

#include "PhysicsTools/Utilities/interface/deltaR.h"

#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/CorJetWithBTagDiscr.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/VBFHZZllbbUtils.h"

using namespace edm;
using namespace reco;
using namespace vbfhzz2l2b;
using namespace std;

VBFHZZllbbBTagInfoAnalyzer::VBFHZZllbbBTagInfoAnalyzer(const edm::ParameterSet& iConfig) :
  impactParameterTagInfosLabel_ ( iConfig.getParameter<edm::InputTag> ( "impactParameterTagInfosLabel" ) ),
  secondaryVertexTagInfosLabel_ ( iConfig.getParameter<edm::InputTag> ( "secondaryVertexTagInfosLabel" ) )
{
   //now do what ever initialization is needed
}

VBFHZZllbbBTagInfoAnalyzer::~VBFHZZllbbBTagInfoAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
VBFHZZllbbBTagInfoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<reco::MuonCollection> muonHandle;
  iEvent.getByLabel ("muons",muonHandle);
  for ( reco::MuonCollection::const_iterator muon_itr = muonHandle->begin();
	muon_itr != muonHandle->end(); ++muon_itr ) {
    if ( muon_itr->isGood(reco::Muon::GlobalMuonPromptTight) ) {
      muonpvX_ -> Fill( muon_itr->vertex().x());
      muonpvY_ -> Fill( muon_itr->vertex().y());
      muonpvZ_ -> Fill( muon_itr->vertex().z());
    }
  }


  edm::Handle<reco::VertexCollection> primaryVerticesWithBSHandle;
  iEvent.getByLabel("offlinePrimaryVerticesWithBS",primaryVerticesWithBSHandle);
  edm::Handle<reco::VertexCollection> primaryVerticesHandle;
  iEvent.getByLabel("offlinePrimaryVertices",primaryVerticesHandle);
  edm::Handle<reco::VertexCollection> pixelVerticesHandle;
  iEvent.getByLabel("pixelVertices",pixelVerticesHandle);

  for ( reco::VertexCollection::const_iterator pV_itr = primaryVerticesHandle->begin();
	pV_itr != primaryVerticesHandle->end();
	++pV_itr ) {
    pvX_ -> Fill(pV_itr -> x());
    pvY_ -> Fill(pV_itr -> y());
    pvZ_ -> Fill(pV_itr -> z());
  }
  for ( reco::VertexCollection::const_iterator pV_itr = primaryVerticesWithBSHandle->begin();
	pV_itr != primaryVerticesWithBSHandle->end();
	++pV_itr ) {
    BSpvX_ -> Fill(pV_itr -> x());
    BSpvY_ -> Fill(pV_itr -> y());
    BSpvZ_ -> Fill(pV_itr -> z());
  }
  for ( reco::VertexCollection::const_iterator pV_itr = pixelVerticesHandle->begin();
	pV_itr != pixelVerticesHandle->end();
	++pV_itr ) {
    PXpvX_ -> Fill(pV_itr -> x());
    PXpvY_ -> Fill(pV_itr -> y());
    PXpvZ_ -> Fill(pV_itr -> z());
  }



  edm::Handle<reco::SecondaryVertexTagInfoCollection> secondaryVtxTagInfosHandle;
  iEvent.getByLabel( secondaryVertexTagInfosLabel_, secondaryVtxTagInfosHandle );  

  const reco::TrackIPTagInfoRef & trkIPtagInfoRef = secondaryVtxTagInfosHandle -> at(0) . trackIPTagInfoRef ();
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
  double x = primVtxColl->x(); // x coordinate
  double y = primVtxColl->y(); // y coordinate
  double z = primVtxColl->z(); // z coordinate
  double xError = primVtxColl->xError(); //	error on x
  double yError = primVtxColl->yError(); //	error on y
  double zError = primVtxColl->zError(); //	error on z
  
  tagpvX_ -> Fill (x);
  tagpvY_ -> Fill (y);
  tagpvZ_ -> Fill (z);

  std::cout << "x: " << x << "+o-" << xError << std::endl;
  std::cout << "y: " << y << "+o-" << yError << std::endl;
  std::cout << "z: " << z << "+o-" << zError << std::endl;
  std::cout << "isFake: " << isFake << " <--> isValid: " << isValid << std::endl;
  std::cout << "chi2: " << chi2 << " ndof: " << ndof << " => normalizedChi2: " << normalizedChi2 << "(" << chi2/ndof << ")" << std::endl;
  std::cout << "tracksSize: " << tracksSize << std::endl;
  std::cout << "****************************************************************" << std::endl;
    

  std::cout << "*************************************************************" << std::endl;
  std::cout << "secondaryVtxTagInfosHandle" << std::endl;
  std::cout << "secondaryVtxTagInfosHandle->size(): " << secondaryVtxTagInfosHandle->size() << std::endl;
  for ( reco::SecondaryVertexTagInfoCollection::const_iterator secondaryVtxTagInfos_itr = secondaryVtxTagInfosHandle->begin();
	secondaryVtxTagInfos_itr != secondaryVtxTagInfosHandle->end(); ++secondaryVtxTagInfos_itr ) {

    // const GlobalVector & flightDirection = secondaryVtxTagInfos_itr -> flightDirection (unsigned int index) const
    // Measurement1D        flightDistance  = secondaryVtxTagInfos_itr -> flightDistance (unsigned int index, bool in2d=false) const
    RefToBase< Jet >        sVjetRef          = secondaryVtxTagInfos_itr -> jet();
    std::cout << "sVjet et: " << sVjetRef->et() <<std::endl;
    const JetTracksAssociationRef & jetTrackAssRef = secondaryVtxTagInfos_itr -> jtaRef();
    unsigned int   nSelTrks    = secondaryVtxTagInfos_itr -> nSelectedTracks ();
    TrackRefVector selTrkColl = secondaryVtxTagInfos_itr -> selectedTracks ();
    TrackRefVector trkColl    = secondaryVtxTagInfos_itr -> tracks (); // returns a list of tracks associated to the jet [belongs to BaseTagInfo]
    std::cout << "nSelTrks:           " << nSelTrks << std::endl;
    std::cout << "selTrkColl.size(): " << selTrkColl.size() << std::endl;
    std::cout << " trkColl.size():   " << trkColl.size() << std::endl;
    for ( unsigned int index = 0; index != nSelTrks; index++ ) {
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

      distanceTo1track_  -> Fill ( selTrkIPData.distanceToFirstTrack );
      distanceToJetAxis_ -> Fill ( selTrkIPData.distanceToJetAxis    );
      ip2d_              -> Fill ( selTrkIPData.ip2d.value()         );
      ip3d_              -> Fill ( selTrkIPData.ip3d.value()         );
      ip2dSig_           -> Fill ( selTrkIPData.ip2d.significance()  );
      ip3dSig_           -> Fill ( selTrkIPData.ip3d.significance()  );

    }
    const TrackIPTagInfoRef & trkIPtagInfoRef = secondaryVtxTagInfos_itr -> trackIPTagInfoRef ();

    std::cout << "---------------------------------------------------------------" << std::endl;
    std::cout << "impactParameterTagInfos" << std::endl;
    std::cout << "(trkIPtagInfoRef.product())->size(): " << (trkIPtagInfoRef.product())->size()<< std::endl;
    unsigned int trackIPTagInfos_index = 0;
    for( reco::TrackIPTagInfoCollection::const_iterator trackIPTagInfos_itr = (trkIPtagInfoRef.product())->begin(); 
	 trackIPTagInfos_itr !=  (trkIPtagInfoRef.product())->end(); ++trackIPTagInfos_itr, trackIPTagInfos_index++ ) {
      bool hasProb = trackIPTagInfos_itr -> hasProbabilities (); // check if probability information is globally available
                                                               // impact parameters in the collection
      if(!hasProb) std::cout << "hasProb is FALSE!!" << std::endl;

      // look @ TrackIPData
      const std::vector<reco::TrackIPTagInfo::TrackIPData> & ipColl = trackIPTagInfos_itr -> impactParameterData (); // vectors of TrackIPData orderd as the selectedTracks()
      std::cout << "ipColl.size(): " << ipColl.size() << std::endl;
      
      // SortCriteria
      //  - IP3DSig 	
      //  - Prob3D 	
      //  - IP2DSig 	
      //  - Prob2D 	
      //  - IP3DValue 	
      //  - IP2DValue 	
  
      std::vector< size_t > sortedIndicesIP3DSig = trackIPTagInfos_itr -> sortedIndexes (reco::TrackIPTagInfo::IP3DSig); // return the list of track index sorted by IP3DSig mode.
      for ( unsigned int index = 0; index != sortedIndicesIP3DSig.size(); index++ ) 
	std::cout << "sortedIndicesIP3DSig[" << index << "]: " << sortedIndicesIP3DSig[index] << std::endl;
      std::vector< size_t > sortedIndicesIP2DSig = trackIPTagInfos_itr -> sortedIndexes (reco::TrackIPTagInfo::IP2DSig); // return the list of track index sorted by IP2DSig mode.
      for ( unsigned int index = 0; index != sortedIndicesIP2DSig.size(); index++ ) 
	std::cout << "sortedIndicesIP2DSig[" << index << "]: " << sortedIndicesIP2DSig[index] << std::endl;
      std::vector< size_t > sortedIndicesIP3DSigWithCut = trackIPTagInfos_itr -> sortedIndexesWithCut (1., reco::TrackIPTagInfo::IP3DSig); // return the list of track index sorted by mode 
                                                                                                                                           // a cut can be specified to select only tracks with IP value or significance > cut 
      for ( unsigned int index = 0; index != sortedIndicesIP3DSigWithCut.size(); index++ ) 
      std::cout << "sortedIndicesIP3DSigWithCut[" << index << "]: " << sortedIndicesIP3DSigWithCut[index] << std::endl;
                                                                                                                   // or probability < cut (according to the specified mode)
      reco::TrackRefVector sortedTracks = trackIPTagInfos_itr -> sortedTracks (sortedIndicesIP3DSig);
      std::cout << "sortedTracks.size(): " << sortedTracks.size() << std::endl;

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

      // take the vector of TrackIPData orderd as the selectedTracks()  vector
      std::vector<reco::TrackIPTagInfo::TrackIPData>::const_iterator ipColl_itr = ipColl.begin();
      for ( reco::TrackRefVector::const_iterator selTrkColl_itr = selTrkColl.begin();
	    selTrkColl_itr != selTrkColl.end(); ++selTrkColl_itr, ++ipColl_itr ) {
	++selTrkNum;
	selTrkSumPt += (*selTrkColl_itr)->pt();

	std::cout << "ip2d.significance: " << ipColl_itr->ip2d.significance() << std::endl;
	std::cout << "ip3d.significance: " << ipColl_itr->ip3d.significance() << std::endl;
	// take tracks with minimum IP significance
	if ( ipColl_itr->ip3d.significance() > 1. ) {
	  selTrkColl_S1.push_back(*selTrkColl_itr);
	  selTrkSumPt_S1 += (*selTrkColl_itr)->pt();
	  ++selTrkNum_S1;
	}
	if ( ipColl_itr->ip3d.significance() > 2. ) {
	  selTrkColl_S2.push_back(*selTrkColl_itr);
	  selTrkSumPt_S2 += (*selTrkColl_itr)->pt();
	  ++selTrkNum_S2;
	}
	if ( ipColl_itr->ip3d.significance() > 3. ) {
	  selTrkColl_S3.push_back(*selTrkColl_itr);
	  selTrkSumPt_S3 += (*selTrkColl_itr)->pt();
	  ++selTrkNum_S3;
	}
      }
      
      // evaluate tag tracks invariant mass
      double bTagTkInvMass    = tracksInvariantMass( selTrkColl    );
      double bTagTkInvMass_S1 = tracksInvariantMass( selTrkColl_S1 );
      double bTagTkInvMass_S2 = tracksInvariantMass( selTrkColl_S2 );
      double bTagTkInvMass_S3 = tracksInvariantMass( selTrkColl_S3 );

      std::cout << " bTagTkInvMass   : " << bTagTkInvMass   
		<< " bTagTkInvMass_S1: " << bTagTkInvMass_S1
		<< " bTagTkInvMass_S2: " << bTagTkInvMass_S2
		<< " bTagTkInvMass_S3: " << bTagTkInvMass_S3 << std::endl;
      
      std::cout << " selTrkSumPt   : " << selTrkSumPt    
		<< " selTrkSumPt_S1: " << selTrkSumPt_S1 
		<< " selTrkSumPt_S2: " << selTrkSumPt_S2 
		<< " selTrkSumPt_S3: " << selTrkSumPt_S3 << std::endl;

      std::cout << " selTrkNum   : " << selTrkNum    
		<< " selTrkNum_S1: " << selTrkNum_S1 
		<< " selTrkNum_S2: " << selTrkNum_S2 
		<< " selTrkNum_S3: " << selTrkNum_S3 << std::endl;

      
      selTrkNum_     -> Fill ( selTrkNum        );
      selTrkNumS1_   -> Fill ( selTrkNum_S1     );
      selTrkNumS2_   -> Fill ( selTrkNum_S2     );
      selTrkNumS3_   -> Fill ( selTrkNum_S3     );
      selTrkSumPt_   -> Fill ( selTrkSumPt      );
      selTrkSumPtS1_ -> Fill ( selTrkSumPt_S1   );
      selTrkSumPtS2_ -> Fill ( selTrkSumPt_S2   );
      selTrkSumPtS3_ -> Fill ( selTrkSumPt_S3   );
      selTrkMass_    -> Fill ( bTagTkInvMass    );
      selTrkMassS1_  -> Fill ( bTagTkInvMass_S1 );
      selTrkMassS2_  -> Fill ( bTagTkInvMass_S2 );
      selTrkMassS3_  -> Fill ( bTagTkInvMass_S3 );      
      
      
    }
    std::cout << "---------------------------------------------------------------" << std::endl;

    unsigned int nSecVtxCand = secondaryVtxTagInfos_itr -> nVertexCandidates ();
    unsigned int nSecVtx     = secondaryVtxTagInfos_itr -> nVertices ();
    std::cout << "nSecVtxCand: " << nSecVtxCand << std::endl;
    std::cout << "nSecVtx:     " << nSecVtx << std::endl;
    for ( unsigned int index = 0; index != nSecVtxCand; index++ ) {
      unsigned int nSecVtxTrks = secondaryVtxTagInfos_itr -> nVertexTracks (index);
      std::cout << "nSecVtxTrks: " << nSecVtxTrks << std::endl;
      const Vertex   secVtx = secondaryVtxTagInfos_itr -> secondaryVertex (index);
      std::cout << "secondary vertex: " << std::endl;
      double x = secVtx.x(); // x coordinate
      double y = secVtx.y(); // y coordinate
      double z = secVtx.z(); // z coordinate
      double xError = secVtx.xError(); //	error on x
      double yError = secVtx.yError(); //	error on y
      double zError = secVtx.zError(); //	error on z
      std::cout << "x: " << x << "+o-" << xError << std::endl;
      std::cout << "y: " << y << "+o-" << yError << std::endl;
      std::cout << "z: " << z << "+o-" << zError << std::endl;

      tagsvX_ -> Fill (x);
      tagsvY_ -> Fill (y);
      tagsvZ_ -> Fill (z);

      TrackRefVector vtxTrkColl = secondaryVtxTagInfos_itr -> vertexTracks (index);
      std::cout << "vtxTrks.size(): " << vtxTrkColl.size() << std::endl;
      unsigned int trkIndex = secondaryVtxTagInfos_itr -> findTrack (vtxTrkColl[0]);
      std::cout << "trkIndex: " << trkIndex << " <--> 0" << std::endl;
    }
  }
  std::cout << "*************************************************************" << std::endl;



}


// ------------ method called once each job just before starting event loop  ------------
void 
VBFHZZllbbBTagInfoAnalyzer::beginJob(const edm::EventSetup&)
{

  TFileDirectory muoninfoSubDir = fs->mkdir( "muoninfo" );
  muonpvX_ = muoninfoSubDir.make<TH1D>("muonpvX","primary vertex x (cm)",100,-0.1,0.1);
  muonpvY_ = muoninfoSubDir.make<TH1D>("muonpvY","primary vertex y (cm)",100,-0.1,0.1);
  muonpvZ_ = muoninfoSubDir.make<TH1D>("muonpvZ","primary vertex z (cm)",100,-10.,10.);

  TFileDirectory pvinfoSubDir = fs->mkdir( "pvinfo" );
  pvX_ = pvinfoSubDir.make<TH1D>("pvX","primary vertex x (cm)",100,-0.1,0.1);
  pvY_ = pvinfoSubDir.make<TH1D>("pvY","primary vertex y (cm)",100,-0.1,0.1);
  pvZ_ = pvinfoSubDir.make<TH1D>("pvZ","primary vertex z (cm)",100,-10.,10.);
  BSpvX_ = pvinfoSubDir.make<TH1D>("BSpvX","primary vertex x (cm)",100,-0.1,0.1);
  BSpvY_ = pvinfoSubDir.make<TH1D>("BSpvY","primary vertex y (cm)",100,-0.1,0.1);
  BSpvZ_ = pvinfoSubDir.make<TH1D>("BSpvZ","primary vertex z (cm)",100,-10.,10.);
  PXpvX_ = pvinfoSubDir.make<TH1D>("PXpvX","primary vertex x (cm)",100,-0.1,0.1);
  PXpvY_ = pvinfoSubDir.make<TH1D>("PXpvY","primary vertex y (cm)",100,-0.1,0.1);
  PXpvZ_ = pvinfoSubDir.make<TH1D>("PXpvZ","primary vertex z (cm)",100,-10.,10.);

  TFileDirectory taginfoSubDir = fs->mkdir( "taginfo" );
  tagpvX_ = taginfoSubDir.make<TH1D>("tagpvX","primary vertex x (cm)",100,-0.1,0.1);
  tagpvY_ = taginfoSubDir.make<TH1D>("tagpvY","primary vertex y (cm)",100,-0.1,0.1);
  tagpvZ_ = taginfoSubDir.make<TH1D>("tagpvZ","primary vertex z (cm)",100,-10.,10.);
  distanceTo1track_  = taginfoSubDir.make<TH1D>("distanceTo1track","distance to first track",50,0.,1.);
  distanceToJetAxis_ = taginfoSubDir.make<TH1D>("distanceToJetAxis","distance to jet axis",100,-1.,1.);
  ip2d_    = taginfoSubDir.make<TH1D>("ip2d","2dim impact parameter",100,-1.,1.);
  ip3d_    = taginfoSubDir.make<TH1D>("ip3d","3dim impact parameter",100,-1.,1.);
  ip2dSig_ = taginfoSubDir.make<TH1D>("ip2dSig","2dim impact parameter significance",100,-10.,15.);
  ip3dSig_ = taginfoSubDir.make<TH1D>("ip3dSig","3dim impact parameter significance",100,-10.,15.);
  tagsvX_ = taginfoSubDir.make<TH1D>("tagsvX","secondary vertex x (cm)",100,-0.1,0.1);
  tagsvY_ = taginfoSubDir.make<TH1D>("tagsvY","secondary vertex y (cm)",100,-0.1,0.1);
  tagsvZ_ = taginfoSubDir.make<TH1D>("tagsvZ","secondary vertex z (cm)",100,-10.,10.);
  selTrkNum_   = taginfoSubDir.make<TH1D>("selTrkNum",  "number of selected tracks",                       10,0.,10.);
  selTrkNumS1_ = taginfoSubDir.make<TH1D>("selTrkNumS1","number of selected tracks w/ ip significance > 1",10,0.,10.);
  selTrkNumS2_ = taginfoSubDir.make<TH1D>("selTrkNumS2","number of selected tracks w/ ip significance > 2",10,0.,10.);
  selTrkNumS3_ = taginfoSubDir.make<TH1D>("selTrkNumS3","number of selected tracks w/ ip significance > 3",10,0.,10.);
  selTrkSumPt_   = taginfoSubDir.make<TH1D>("selTrkSumPt",  "Sum p_{T} of selected tracks",                       100,0.,100.);
  selTrkSumPtS1_ = taginfoSubDir.make<TH1D>("selTrkSumPtS1","Sum p_{T} of selected tracks w/ ip significance > 1",100,0.,100.);
  selTrkSumPtS2_ = taginfoSubDir.make<TH1D>("selTrkSumPtS2","Sum p_{T} of selected tracks w/ ip significance > 2",100,0.,100.);
  selTrkSumPtS3_ = taginfoSubDir.make<TH1D>("selTrkSumPtS3","Sum p_{T} of selected tracks w/ ip significance > 3",100,0.,100.);
  selTrkMass_   = taginfoSubDir.make<TH1D>("selTrkMass",  "invariant mass of selected tracks",                       250,0.,25.);
  selTrkMassS1_ = taginfoSubDir.make<TH1D>("selTrkMassS1","invariant mass of selected tracks w/ ip significance > 1",250,0.,25.);
  selTrkMassS2_ = taginfoSubDir.make<TH1D>("selTrkMassS2","invariant mass of selected tracks w/ ip significance > 2",250,0.,25.);
  selTrkMassS3_ = taginfoSubDir.make<TH1D>("selTrkMassS3","invariant mass of selected tracks w/ ip significance > 3",250,0.,25.);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
VBFHZZllbbBTagInfoAnalyzer::endJob() {
}

/** Taken from CMSSW/RecoTauTag/RecoTau/src/CaloRecoTauAlgorithm.cc
 * Evaluate jet tracks invariant mass
 */
double VBFHZZllbbBTagInfoAnalyzer::tracksInvariantMass( const reco::TrackRefVector & selectedTracks ) const {
  double chargedPI_mass_ = 0.13957018; // GeV/c^2

  int trksChargeSum = 0;
  // setting invariant mass of the tracks system
  TLorentzVector tracks4DVec(0.,0.,0.,0.);

  if( int(selectedTracks.size()) != 0 ) {
    int tracksCharge = 0;

    for ( int index=0; index < (int)selectedTracks.size(); ++index ) {
      tracksCharge += (selectedTracks)[index]->charge();
      TLorentzVector chargedPIcand_fromTrk_4DVec( (selectedTracks)[index]->momentum().x(),
						  (selectedTracks)[index]->momentum().y(),
						  (selectedTracks)[index]->momentum().z(),
						  sqrt(pow((double)(selectedTracks)[index]->momentum().r(),2) + pow(chargedPI_mass_,2)));
      
      tracks4DVec += chargedPIcand_fromTrk_4DVec;
    }
    trksChargeSum = tracksCharge;    
  }

  std::cout << "trksChargeSum: " << trksChargeSum << " tagMass: " << tracks4DVec.M() << std::endl;

  return tracks4DVec.M();
}

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
// $Id: VBFHZZllbbBTagInfoAnalyzer.cc,v 1.1 2009/05/14 10:52:00 tosi Exp $
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

  edm::Handle<reco::SecondaryVertexTagInfoCollection> secondaryVtxTagInfosHandle;
  iEvent.getByLabel( secondaryVertexTagInfosLabel_, secondaryVtxTagInfosHandle );  

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
      TrackRef selTrkRef = secondaryVtxTagInfos_itr -> track (index);
      const reco::SecondaryVertexTagInfo::TrackData trkData1  = secondaryVtxTagInfos_itr -> trackData (index);
      const reco::SecondaryVertexTagInfo::TrackData trkData2  = secondaryVtxTagInfos_itr -> trackData (selTrkRef);
      const reco::TrackIPTagInfo::TrackIPData trkIPData1 = secondaryVtxTagInfos_itr -> trackIPData (index);
      const reco::TrackIPTagInfo::TrackIPData trkIPData2 = secondaryVtxTagInfos_itr -> trackIPData (selTrkRef);
      std::cout << " trkIPData1 -> distanceToFirstTrack: "  << trkIPData1.distanceToFirstTrack << std::endl;
      std::cout << " trkIPData1 -> distanceToJetAxis:    "  << trkIPData1.distanceToJetAxis    << std::endl;
      std::cout << " trkIPData1 -> ip2d:		 "  << trkIPData1.ip2d.significance()  << std::endl;   
      std::cout << " trkIPData1 -> ip3d:                 "  << trkIPData1.ip3d.significance()  << std::endl;
    }
    const TrackIPTagInfoRef & trkIPtagInfoRef = secondaryVtxTagInfos_itr -> trackIPTagInfoRef ();

    std::cout << "---------------------------------------------------------------" << std::endl;
    std::cout << "impactParameterTagInfos" << std::endl;
    std::cout << "(trkIPtagInfoRef.product())->size(): " << (trkIPtagInfoRef.product())->size()<< std::endl;
    
  
    unsigned int trackIPTagInfos_index = 0;
    for( reco::TrackIPTagInfoCollection::const_iterator trackIPTagInfos_itr = (trkIPtagInfoRef.product())->begin(); 
	 trackIPTagInfos_itr !=  (trkIPtagInfoRef.product())->end(); ++trackIPTagInfos_itr, trackIPTagInfos_index++ ) {
      if ( trackIPTagInfos_index == 0 ) {
	edm::Ref< reco::VertexCollection > primVtxColl = trackIPTagInfos_itr->primaryVertex ();
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
	//      float  trackWeight = primVtxColl->trackWeight(const TrackRef &r); // returns the weight with which a Track has contributed to the vertex-fit.
	double x = primVtxColl->x(); // x coordinate
	double y = primVtxColl->y(); // y coordinate
	double z = primVtxColl->z(); // z coordinate
	double xError = primVtxColl->xError(); //	error on x
	double yError = primVtxColl->yError(); //	error on y
	double zError = primVtxColl->zError(); //	error on z
	
	std::cout << "x: " << x << "+o-" << xError << std::endl;
	std::cout << "y: " << y << "+o-" << yError << std::endl;
	std::cout << "z: " << z << "+o-" << zError << std::endl;
	std::cout << "isFake: " << isFake << " <--> isValid: " << isValid << std::endl;
	std::cout << "chi2: " << chi2 << " ndof: " << ndof << " => normalizedChi2: " << normalizedChi2 << "(" << chi2/ndof << ")" << std::endl;
	std::cout << "tracksSize: " << tracksSize << std::endl;
      }

      bool hasProb = trackIPTagInfos_itr -> hasProbabilities (); // check if probability information is globally available
                                                               // impact parameters in the collection
      if(!hasProb) std::cout << "hasProb is FALSE!!" << std::endl;

      // look @ TrackIPData
      const std::vector<reco::TrackIPTagInfo::TrackIPData> & ipColl = trackIPTagInfos_itr -> impactParameterData (); // vectors of TrackIPData orderd as the selectedTracks()
      std::cout << "ipColl.size(): " << ipColl.size() << std::endl;
      for ( unsigned int index = 0; index != ipColl.size(); index++ ) {
	
	const std::vector<float> & probabilities = trackIPTagInfos_itr -> probabilities (index);
	std::cout << "ipColl[index]: " << index << std::endl;
	std::cout << "probabilities.size(): " << probabilities.size() << std::endl;
	GlobalPoint 	closestToFirstTrk = ipColl[index].closestToFirstTrack;
	GlobalPoint 	closestToJet      = ipColl[index].closestToJetAxis;
	float 	        distanceToFirstTrk = ipColl[index].distanceToFirstTrack;
	float 	        distanceToJet      = ipColl[index].distanceToJetAxis;
	Measurement1D 	ip2d = ipColl[index].ip2d;
	Measurement1D 	ip3d = ipColl[index].ip3d;
	std::cout << " trkIPData1 -> distanceToFirstTrack: "  << ipColl[index].distanceToFirstTrack << std::endl;
	std::cout << " trkIPData1 -> distanceToJetAxis:    "  << ipColl[index].distanceToJetAxis    << std::endl;
	std::cout << " trkIPData1 -> ip2d:		 "  << ipColl[index].ip2d.significance()  << std::endl;   
	std::cout << " trkIPData1 -> ip3d:                 "  << ipColl[index].ip3d.significance()  << std::endl;
      }
      
      RefToBase< Jet > trkIPjetRef = trackIPTagInfos_itr -> jet();
      std::cout << "trkIPjet et: " << trkIPjetRef->et() <<std::endl;
      const JetTracksAssociationRef & jetTrackAssRef = trackIPTagInfos_itr -> jtaRef();
      const reco::TrackRefVector & selTrkColl = trackIPTagInfos_itr -> selectedTracks(); // return the vector of tracks for which the IP information is available 
                                                                                         // quality cuts are applied to reject fake tracks
      TrackRefVector trkColl    = trackIPTagInfos_itr -> tracks (); // returns a list of tracks associated to the jet [belongs to BaseTagInfo]
      std::cout << "selTrkColl.size(): " << selTrkColl.size() << std::endl;
      std::cout << " trkColl.size():   " << trkColl.size() << std::endl;
      
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

      int probNum_0 = trackIPTagInfos_itr->probabilities(0).size();
      int probNum_1 = trackIPTagInfos_itr->probabilities(1).size();

      std::cout << "probNum_0: " << probNum_0 << " <--> probNum_1: " << probNum_1 << std::endl;
      
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
      
    }
    std::cout << "---------------------------------------------------------------" << std::endl;

    unsigned int nVtxCand = secondaryVtxTagInfos_itr -> nVertexCandidates ();
    unsigned int nVtx     = secondaryVtxTagInfos_itr -> nVertices ();
    std::cout << "nVtxCand: " << nVtxCand << std::endl;
    std::cout << "nVtx:     " << nVtx << std::endl;
    for ( unsigned int index = 0; index != nVtxCand; index++ ) {
      unsigned int nVtxTrks = secondaryVtxTagInfos_itr -> nVertexTracks (index);
      std::cout << "nVtxTrks: " << nVtxTrks << std::endl;
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

  TFileDirectory bjetInfoSubDir = fs->mkdir( "bjetInfo" );
  TFileDirectory ipTagInfoSubDir     = bjetInfoSubDir.mkdir( "ipTagInfo" );
  TFileDirectory secVtxTagInfoSubDir = bjetInfoSubDir.mkdir( "secVtxTagInfo" );


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

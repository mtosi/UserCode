#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/VBFHZZllbbElectronTrackIsolationAlgos.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "Math/VectorUtil.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

// ------------------------------------------------------------------------------------------------


VBFHZZllbbElectronTrackIsolationAlgos::VBFHZZllbbElectronTrackIsolationAlgos (double coneRadius,
									      double vetoRadius,
									      double otherVetoRadius,
									      double ptMin,
									      double lipMax,
									      bool   useTrkQuality) :
  coneRadius_      (coneRadius) ,
  vetoRadius_      (vetoRadius) ,
  otherVetoRadius_ (otherVetoRadius) ,
  ptMin_           (ptMin) ,
  lipMax_          (lipMax) ,
  useTrkQuality_   (useTrkQuality)
{
  /*
  std::cout << "coneRadius_:      " << coneRadius_      << std::endl;
  std::cout << "vetoRadius_:      " << vetoRadius_      << std::endl;
  std::cout << "otherVetoRadius_: " << otherVetoRadius_ << std::endl;
  std::cout << "ptMin_:           " << ptMin_           << std::endl;
  std::cout << "lipMax_:          " << lipMax_          << std::endl;
  std::cout << "useTrkQuality_:   " << useTrkQuality_   << std::endl;
  */

}                                                        


// ------------------------------------------------------------------------------------------------


//VBFHZZllbbElectronTrackIsolationAlgos::~VBFHZZllbbElectronTrackIsolationAlgos () {}


// ------------------------------------------------------------------------------------------------


int 
VBFHZZllbbElectronTrackIsolationAlgos::numOfTrks (const edm::Handle<reco::PixelMatchGsfElectronCollection> & electronHandle,
						  const edm::Handle<reco::TrackCollection>                 & trackHandle,
						  const reco::PixelMatchGsfElectron                        & mainElectron) const 
{
  int trkCounter = 0 ;

  // take the electron track
  math::XYZVector tmpElectronMomentumAtVtx = mainElectron.trackMomentumAtVtx () ; 
  math::XYZVector tmpElectronPositionAtVtx = mainElectron.TrackPositionAtVtx () ; 

  // loop over track collection
  for (reco::TrackCollection::const_iterator track_itr = trackHandle->begin () ;
       track_itr != trackHandle->end ();  ++track_itr) {    

    // min pT threshold
    math::XYZVector tmpTrackMomentumAtVtx = track_itr->innerMomentum () ; 
    double tmpTrackMomentumAtVtx_pt  = sqrt (tmpTrackMomentumAtVtx.Perp2 ()) ;
    
    int nhits = track_itr->recHitsSize () ;
    if ( nhits < 5 ) continue;
    if ( tmpTrackMomentumAtVtx_pt < ptMin_ ) continue ;  
    
    // impact parameter threshold
    if (fabs( track_itr->vz () - tmpElectronPositionAtVtx.z() ) > lipMax_) continue ;
    
    if (useTrkQuality_ && !vbfhzz2l2b::testTrackerTrack (track_itr)) continue ;
    
    bool countTrack = true ;
    if (otherVetoRadius_ > 0.0001) { // to avoid useless caclulations 
      
      // loop over electrons
      for (reco::PixelMatchGsfElectronCollection::const_iterator ele_itr = electronHandle->begin () ; 
	   ele_itr != electronHandle->end () ; ++ele_itr) {
	
	if (&mainElectron == &*ele_itr) {
	  std::cout << "same electron" << std::endl;
	  continue ;
	}
	
	math::XYZVector eleMomentumAtVtx = ele_itr->trackMomentumAtVtx () ; 
	double eleDR = ROOT::Math::VectorUtil::DeltaR (tmpTrackMomentumAtVtx,eleMomentumAtVtx) ;
	if (eleDR < otherVetoRadius_) {
	  countTrack = false ;
	  break ;
	}
      } // loop over electrons
    } // if (otherVetoRadius_ > 0.0001)
    
    double dR = ROOT::Math::VectorUtil::DeltaR (tmpTrackMomentumAtVtx,tmpElectronMomentumAtVtx) ;
    if ( countTrack && 
	 (dR < coneRadius_ && dR >= vetoRadius_ ) )
      ++trkCounter ;
    
  } // loop over tracks
  
  return trkCounter ;

}                                              


// ------------------------------------------------------------------------------------------------


double 
VBFHZZllbbElectronTrackIsolationAlgos::sumPt (const edm::Handle<reco::PixelMatchGsfElectronCollection> & electronHandle,
					      const edm::Handle<reco::TrackCollection>                 & trackHandle,
					      const reco::PixelMatchGsfElectron                        & mainElectron) const
{
  double trkPtSum = 0 ;

  // take the electron track
  math::XYZVector tmpElectronMomentumAtVtx = mainElectron.trackMomentumAtVtx () ; 
  math::XYZVector tmpElectronPositionAtVtx = mainElectron.TrackPositionAtVtx () ; 

  // loop over track collection
  for (reco::TrackCollection::const_iterator track_itr = trackHandle->begin () ;
       track_itr != trackHandle->end () ; ++track_itr) {    

    // min pT threshold
    math::XYZVector tmpTrackMomentumAtVtx = track_itr->innerMomentum () ; 
    double tmpTrackMomentumAtVtx_pt  = sqrt (tmpTrackMomentumAtVtx.Perp2 ()) ;
    
    int nhits = track_itr->recHitsSize () ;
    if ( nhits < 5 ) continue;
    if ( tmpTrackMomentumAtVtx_pt < ptMin_ ) continue ;  
    
    // impact parameter threshold
    if (fabs( track_itr->vz () - tmpElectronPositionAtVtx.z() ) > lipMax_) continue ;

    if (useTrkQuality_ && !vbfhzz2l2b::testTrackerTrack (track_itr)) continue ;

    bool countTrack = true ;
    if (otherVetoRadius_ > 0.0001) { // to avoid useless caclulations 

      // loop over electrons
      for (reco::PixelMatchGsfElectronCollection::const_iterator ele_itr = electronHandle->begin () ; 
	   ele_itr != electronHandle->end () ; ++ele_itr) {

	if (&mainElectron == &*ele_itr) {
	  std::cout << "same electron" << std::endl;
	  continue ;
	}
	
/*
	edm::RefToBase<reco::PixelMatchGsfElectron> electronBaseRef = electronHandle->at (ele_itr - electronHandle->begin ()) ;
	reco::PixelMatchGsfElectronRef              electronRef     = electronBaseRef.castTo<reco::PixelMatchGsfElectronRef> () ;
	if (electronRef == mainElectron) continue ;
*/
	
	math::XYZVector eleMomentumAtVtx = ele_itr->trackMomentumAtVtx () ; 
	double eleDR = ROOT::Math::VectorUtil::DeltaR (tmpTrackMomentumAtVtx,eleMomentumAtVtx) ;
	if (eleDR < otherVetoRadius_) {
	  countTrack = false ;
	  break ;
	}
      } // loop over electrons
    } // if (otherVetoRadius_ > 0.0001)
    
    double dR = ROOT::Math::VectorUtil::DeltaR (tmpTrackMomentumAtVtx,tmpElectronMomentumAtVtx) ;
    if ( countTrack && 
	 (dR < coneRadius_ && dR >= vetoRadius_ ) )
      trkPtSum += tmpTrackMomentumAtVtx_pt ;

  } // loop over tracks
  
  return trkPtSum ;
}


// ------------------------------------------------------------------------------------------------


double
VBFHZZllbbElectronTrackIsolationAlgos::isoValue (const edm::Handle<reco::PixelMatchGsfElectronCollection> & electronHandle,
						 const edm::Handle<reco::TrackCollection>                 & trackHandle,
						 const reco::PixelMatchGsfElectron                        & mainElectron   ) const
{

  double isoVal = sumPt (electronHandle, trackHandle, mainElectron) ;
  isoVal /= mainElectron.pt () ;

  return isoVal ;
}


// ------------------------------------------------------------------------------------------------




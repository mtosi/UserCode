#ifndef VBFHZZllbbElectronTrackIsolationAlgos_h
#define VBFHZZllbbElectronTrackIsolationAlgos_h

// my includes
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackExtraBase.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/Handle.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/VBFHZZllbbUtils.h"

////CLHEP
//#include "Math/LorentzVector.h"

class VBFHZZllbbElectronTrackIsolationAlgos {

  public:
  
    //! constructor
    VBFHZZllbbElectronTrackIsolationAlgos (){}
    VBFHZZllbbElectronTrackIsolationAlgos (double coneRadius,
					   double vetoRadius,
					   double otherVetoRadius,
					   double ptMin,
					   double lipMax,
					   bool   useTrkQuality = true) ;
  
    //!destructor 
    ~VBFHZZllbbElectronTrackIsolationAlgos () {};
  
    int numOfTrks (const edm::Handle<reco::PixelMatchGsfElectronCollection> & electronHandle,
		   const edm::Handle<reco::TrackCollection>                 & trackHandle,
		   const reco::PixelMatchGsfElectron                        & mainElectron) const ;

    double sumPt (const edm::Handle<reco::PixelMatchGsfElectronCollection> & electronHandle,
		  const edm::Handle<reco::TrackCollection>                 & trackHandle,
		  const reco::PixelMatchGsfElectron                        & mainElectron) const ;

    double isoValue (const edm::Handle<reco::PixelMatchGsfElectronCollection> & electronHandle,
		     const edm::Handle<reco::TrackCollection>                 & trackHandle,
		     const reco::PixelMatchGsfElectron                        & mainElectron) const ;
  

  private :
  
    double coneRadius_ ;
    double vetoRadius_ ;
    double otherVetoRadius_ ;
    double ptMin_ ;
    double lipMax_ ;
    
    bool useTrkQuality_ ;

};

#endif

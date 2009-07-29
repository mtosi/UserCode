#ifndef CandidateTkIsolation_h
#define CandidateTkIsolation_h

//C++ includes
//#include <vector>
//#include <functional>

//CMSSW includes 
//#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"


class CandidateTkIsolation {
 public:
  
  //constructors
  CandidateTkIsolation () ;
  CandidateTkIsolation (const reco::Candidate * c, 
			const reco::CandidateCollection * c1coll ) ;
  CandidateTkIsolation (const reco::Candidate * c, int trackertrack ) ;
  CandidateTkIsolation (const reco::GsfElectron * electron, 
			const edm::View<reco::Track>* trackCollection,
			const edm::View<reco::GsfElectron>* electronCollection,
			int trackertrack); 

  CandidateTkIsolation (const reco::GsfElectron * electron, 
			const edm::View<reco::Track>* trackCollection,
			const edm::View<reco::GsfElectron>* electronCollection,
			const edm::View<reco::Muon>* muonCollection,
			int trackertrack); 

   //methods
  void setExtRadius (double extRadius) ;
  void setIntRadius (double intRadius) ;
  void setPtLow (double ptLow) ;
  void setLip (double lip) ;
  void setZmassinf ( double Zmassinf);
  void setZmasssup ( double Zmasssup);


  int getNumberTracks(const reco::CandidateCollection * c1coll ) const ;
  double getPtTracks (const reco::CandidateCollection * c1coll ) const ;
  double getPtTracks2(const reco::CandidateCollection * c1coll ) const ;
  std::vector<double> getPtTracksCorr(const reco::CandidateCollection * c1coll,const reco::CandidateCollection * c2coll ) const ;

  int getNumberTracks() const ;
  double getPtTracks () const ;
  double getPtTracks2() const ;
  std::vector<double> getPtTracksCorr() const ;

  bool testTrackerTrack(edm::View<reco::Track>::const_iterator iterTrack, const reco::Muon* muon);

  //destructor 
  ~CandidateTkIsolation() ;
  
 private:

  const reco::GsfElectron*  electron_ ;
  const edm::View<reco::Track> *trackCollection_ ;
  const edm::View<reco::GsfElectron> *electronCollection_ ;
  const edm::View<reco::Muon> *muonCollection_ ;

  double extRadius_ ;
  double intRadius_ ;
  double ptLow_ ;
  double lip_ ;
  float ZMASSPDG_;
  double Zmass_inf_;
  double Zmass_sup_;

  bool tracker_track_;
  bool track_quality_tt(double d0, double ed0,double dz, double edz, int nhits, double pt) const;
  bool track_quality(double pt, double dz) const;
    
};

#endif

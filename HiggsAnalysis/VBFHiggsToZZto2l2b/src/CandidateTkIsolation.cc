
//C++ includes
#include <vector>
#include <functional>

//ROOT includes
#include <Math/VectorUtil.h>

//CMSSW includes
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/CandidateTkIsolation.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/TrajectoryCleaning/interface/TrajectoryCleanerBySharedHits.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h"
#include "TrackingTools/PatternTools/interface/TrajectoryFitter.h"
#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

using namespace std ;
using namespace reco;
using namespace ROOT::Math::VectorUtil ;

CandidateTkIsolation::CandidateTkIsolation ()
{
}


CandidateTkIsolation::CandidateTkIsolation (const reco::Candidate * electron,
					    const reco::CandidateCollection * trackCollection)  
{
  extRadius_ = 0.25 ;
  intRadius_ = 0.015 ;
  ptLow_     = 1.5 ; 
  lip_       = 0.1 ; 
  Zmass_inf_ = 40;
  Zmass_sup_ = 10;
  ZMASSPDG_  = 91.1876;
}

CandidateTkIsolation::CandidateTkIsolation (const reco::Candidate *electron,  int trackertrack)
{

  extRadius_ = 0.25 ;
  intRadius_ = 0.015 ;
  ptLow_     = 1.5 ; 
  lip_       = 0.1 ; 
  Zmass_inf_ = 40;
  Zmass_sup_ = 10;
  ZMASSPDG_  = 91.1876;
  if      ( trackertrack==1 ) tracker_track_ = true;
  else if ( trackertrack==0 ) tracker_track_ = false;
  else {
    cout<< "TkIsolation::not a good value for tracker track, chosed trackertrack= false"<<endl;
    tracker_track_= false;
  }
}  


CandidateTkIsolation::CandidateTkIsolation (const reco::GsfElectron* electron, 
					    const edm::View<reco::Track>* trackCollection,
					    const edm::View<reco::GsfElectron>* electronCollection,
					    int trackertrack) : 
  electron_(electron) ,
  trackCollection_(trackCollection) ,
  electronCollection_(electronCollection) 
{
  extRadius_ = 0.25 ;
  intRadius_ = 0.015 ;
  ptLow_     = 1.5 ; 
  lip_       = 0.1 ; 
  Zmass_inf_ = 40;
  Zmass_sup_ = 10;
  ZMASSPDG_  = 91.1876;
  if      (trackertrack==1) tracker_track_= true;
  else if (trackertrack==0) tracker_track_= false;
  else {
    cout<< "TkIsolation::not a good value for tracker track, chosed trackertrack= false"<<endl;
    tracker_track_= false;
  }
}  

CandidateTkIsolation::CandidateTkIsolation (const reco::GsfElectron* electron, 
					    const edm::View<reco::Track>* trackCollection,
					    const edm::View<reco::GsfElectron>* electronCollection,
					    const edm::View<reco::Muon> *muonCollection,
					    int trackertrack) : 
  electron_(electron) ,
  trackCollection_(trackCollection) ,
  electronCollection_(electronCollection) ,
  muonCollection_(muonCollection)
{
  extRadius_ = 0.25 ;
  intRadius_ = 0.015 ;
  ptLow_     = 1.5 ; 
  lip_       = 0.1 ; 
  Zmass_inf_ = 40;
  Zmass_sup_ = 10;
  ZMASSPDG_  = 91.1876;
  if      (trackertrack==1) tracker_track_= true;
  else if (trackertrack==0) tracker_track_= false;
  else {
    cout<< "TkIsolation::not a good value for tracker track, chosed trackertrack= false"<<endl;
    tracker_track_= false;
  }
} 


CandidateTkIsolation::~CandidateTkIsolation ()
{
}

void CandidateTkIsolation::setExtRadius (double extRadius)
{
  extRadius_ = extRadius ;
}

void CandidateTkIsolation::setIntRadius (double intRadius)
{  
  intRadius_ = intRadius ;
}

void CandidateTkIsolation::setPtLow (double ptLow)
{  
  ptLow_ = ptLow ;
}

void CandidateTkIsolation::setLip (double lip)
{  
  lip_ = lip ;
}

void CandidateTkIsolation::setZmassinf ( double Zmassinf)
{
  Zmass_inf_ = Zmassinf;
}

void CandidateTkIsolation::setZmasssup ( double Zmasssup)
{
  Zmass_sup_ = Zmasssup;
}

int CandidateTkIsolation::getNumberTracks () const
{  
  //counter for the tracks in the isolation cone
  int dummyCounter = 0 ;    

  //Take the electron track
  if(electron_!=0) {
    
    reco::GsfTrackRef tmpTrack = electron_->gsfTrack() ;
    math::XYZVector tmpElectronMomentumAtVtx = (*tmpTrack).innerMomentum () ; 
    
    //    for (reco::CandidateCollection::const_iterator trackcand = trackCollection->begin(); 
    //	 trackcand != trackCollection->end(); ++trackcand ) {

    //const Candidate * cand;
    //cand = &* trackcand->masterClone();
    //const reco::Track * itrTr = dynamic_cast<const reco::Track *>(cand);
    for ( edm::View<reco::Track>::const_iterator itrTr  = trackCollection_->begin() ; 
      itrTr != trackCollection_->end(); ++itrTr ) {
      //math::XYZVector tmpTrackMomentumAtVtx = (*itrTr).innerMomentum () ;
      math::XYZVector tmpTrackMomentumAtVtx = itrTr->innerMomentum () ;
      double this_pt  = sqrt( tmpTrackMomentumAtVtx.Perp2 () );
      
      double d0 = fabs(itrTr->d0());
      double ed0 = itrTr->d0Error();
      double dz = itrTr->dz();
      double edz_tr = itrTr->dzError();
      double edz_ele= (*tmpTrack).dzError();
      double deltaz= fabs( dz - (*tmpTrack).dz());
      double edeltaz= sqrt(edz_tr*edz_tr + edz_ele*edz_ele);
      int nhits = itrTr->recHitsSize();
      if (tracker_track_ && !track_quality_tt(d0, ed0, deltaz, edeltaz, nhits,this_pt) ) continue;
      if (!tracker_track_ && !track_quality(this_pt,deltaz)) continue;
      
      double dr = DeltaR(tmpTrackMomentumAtVtx,tmpElectronMomentumAtVtx) ;
      if ( fabs(dr) < extRadius_ && 
	   fabs(dr) > intRadius_ )
	++dummyCounter ;  
    }//end loop over tracks                 
    
  }

  return dummyCounter ;
}

//double CandidateTkIsolation::getPtTracks (const reco::CandidateCollection * trackCollection) const
double CandidateTkIsolation::getPtTracks () const
{
  //dummy counter for the pT of tracks inside the cone
  double dummypT = 0 ;

  //Take the electron track
  if(electron_!=0) {
    
    reco::GsfTrackRef tmpTrack = electron_->gsfTrack() ;
    math::XYZVector tmpElectronMomentumAtVtx = (*tmpTrack).innerMomentum () ; 
    
    //    for (reco::CandidateCollection::const_iterator trackcand = trackCollection->begin(); 
    //	 trackcand != trackCollection->end(); ++trackcand ) {
      
    // const Candidate * cand;
    // cand = &* trackcand->masterClone();
    // const reco::Track * itrTr = dynamic_cast<const reco::Track *>(cand);
      
    for ( edm::View<reco::Track>::const_iterator itrTr  = (*trackCollection_).begin() ; 
     	  itrTr != (*trackCollection_).end()   ; 
      	  ++itrTr) 
      {
      math::XYZVector tmpTrackMomentumAtVtx = itrTr->innerMomentum () ; 
      double this_pt  = sqrt( tmpTrackMomentumAtVtx.Perp2 () );
      double d0 = fabs(itrTr->d0());
      double ed0 = itrTr->d0Error();
      double dz = itrTr->dz();
      double edz_tr = itrTr->dzError();
      double edz_ele= (*tmpTrack).dzError();
      double deltaz= fabs( dz - (*tmpTrack).dz());
      double edeltaz= sqrt(edz_tr*edz_tr + edz_ele*edz_ele);
      int nhits = itrTr->recHitsSize();
      if (tracker_track_ && !track_quality_tt(d0, ed0, deltaz, edeltaz, nhits,this_pt) ) continue;
      if (!tracker_track_ && !track_quality(this_pt, deltaz)) continue;
      
      double dr = DeltaR(tmpTrackMomentumAtVtx,tmpElectronMomentumAtVtx) ;
      if ( fabs(dr) < extRadius_ && 
	   fabs(dr) > intRadius_ )
	dummypT += this_pt ;
    }
    
  }
  return dummypT ;
}

//double CandidateTkIsolation::getPtTracks2 (const reco::CandidateCollection * trackCollection) const
double CandidateTkIsolation::getPtTracks2 () const
{
  //dummy counter for the pT of tracks inside the cone
  double dummypT = 0 ;
  
  if(electron_!=0) {
    //Take the electron track
    reco::GsfTrackRef tmpTrack = electron_->gsfTrack() ;
    math::XYZVector tmpElectronMomentumAtVtx = (*tmpTrack).innerMomentum () ; 
    
    //  for (reco::CandidateCollection::const_iterator trackcand = trackCollection->begin(); 
    //	 trackcand != trackCollection->end(); ++trackcand ) {
      
    //const Candidate * cand;
    //cand = &* trackcand->masterClone();
    //const reco::Track * itrTr = dynamic_cast<const reco::Track *>(cand);
      
    for ( edm::View<reco::Track>::const_iterator itrTr  = (*trackCollection_).begin() ; 
	  itrTr != (*trackCollection_).end()   ; 
	  ++itrTr) 
      {
      math::XYZVector tmpTrackMomentumAtVtx = itrTr->innerMomentum () ; 
      double this_pt  = sqrt( tmpTrackMomentumAtVtx.Perp2 () );
      double d0 = fabs(itrTr->d0());
      double ed0 = itrTr->d0Error();
      double dz = itrTr->dz();
      double edz_tr = itrTr->dzError();
      double edz_ele= (*tmpTrack).dzError();
      double deltaz= fabs( dz - (*tmpTrack).dz());
      double edeltaz= sqrt(edz_tr*edz_tr + edz_ele*edz_ele);
      int nhits = itrTr->recHitsSize();
      if (tracker_track_ && !track_quality_tt(d0, ed0, deltaz, edeltaz, nhits,this_pt) ) continue;
      if (!tracker_track_ && !track_quality(this_pt, deltaz)) continue;
      
      double dr = DeltaR(tmpTrackMomentumAtVtx,tmpElectronMomentumAtVtx) ;
      if ( fabs(dr) < extRadius_ && 
	   fabs(dr) > intRadius_ )
	dummypT += this_pt*this_pt ;
    }
    
  }
  return dummypT ;
}

//std::vector<double> CandidateTkIsolation::getPtTracksCorr (const reco::CandidateCollection * trackCollection,
//						      const reco::CandidateCollection * electronCollection) const
//{
std::vector<double> CandidateTkIsolation::getPtTracksCorr () const
{
  //dummy counter for the pT of tracks inside the cone
  //cout<<"getpttrackscor method"<<endl;
  vector<double> dummypTout;
  double dummypT=0,dummypT1=0,dummypT2=0;
  double tmpElPt1 = 0.,tmpElPt2 =0, drTrackMatch = 0.;
  //Take the electron track
  //  if(electron_!=0) {
  
  reco::GsfTrackRef tmpTrack = electron_->gsfTrack() ;
  math::XYZVector thisElectronMomentumAtVtx = (*tmpTrack).innerMomentum () ; 
  
  //    for (reco::CandidateCollection::const_iterator trackcand = trackCollection->begin(); 
  //	 trackcand != trackCollection->end(); ++trackcand ) {
  
  // const Candidate * cand;
  // cand = &* trackcand->masterClone();
  //const reco::Track * itrTr = dynamic_cast<const reco::Track *>(cand);

  //track

  for ( edm::View<reco::Track>::const_iterator itrTr  = (*trackCollection_).begin() ; 
	itrTr != (*trackCollection_).end()   ; 
	++itrTr) {
    math::XYZVector tmpTrackMomentumAtVtx = itrTr->innerMomentum () ; 
    double this_pt  = sqrt( tmpTrackMomentumAtVtx.Perp2 () );
    double d0 = fabs(itrTr->d0());
    double ed0 = itrTr->d0Error();
    double dz = itrTr->dz();
    double edz_tr = itrTr->dzError();
    double edz_ele= (*tmpTrack).dzError();
    double deltaz= fabs( dz - (*tmpTrack).dz());
    double edeltaz= sqrt(edz_tr*edz_tr + edz_ele*edz_ele);
    //if (deltaz!= fabs((*itrTr).vz()  - (*tmpTrack).vz()))	cout<<"dz= "<<deltaz<<" vz= "<<fabs((*itrTr).vz()  - (*tmpTrack).vz())<<endl;
    int nhits = itrTr->recHitsSize();
    if (tracker_track_ && !track_quality_tt(d0, ed0, deltaz, edeltaz, nhits,this_pt) ) continue;
    if (!tracker_track_ && !track_quality(this_pt, deltaz)) continue;
    
    
    double dr = DeltaR(tmpTrackMomentumAtVtx,thisElectronMomentumAtVtx) ;
    if ( fabs(dr) < extRadius_ && 
	 fabs(dr) > intRadius_ )
      dummypT += this_pt ;
  }

  if (dummypT==0) cout << "Electron isolated: no tracks around it.." << endl;
  

  for ( edm::View<reco::Muon>::const_iterator itrTr  = muonCollection_->begin() ; 
 	itrTr != muonCollection_->end()   ; 
 	++itrTr) {
    //cout << "Muon eta" << itrTr->p4().eta() << endl;
    math::XYZVector tmpTrackMomentumAtVtxMu(itrTr->p4().x(),itrTr->p4().y(),itrTr->p4().z()) ; 
    double thismuon_pt  = sqrt( tmpTrackMomentumAtVtxMu.Perp2 () );
    
//     //     double d0 = fabs((*itrTr).combinedMuon().d0());
//     // //     double ed0 = (*itrTr).d0Error();
//     // //     double dz = (*itrTr).dz();
//     // //     double edz_tr = (*itrTr).dzError();
//     // //     double edz_ele= (*tmpTrack).dzError();
//     // //     double deltaz= fabs( dz - (*tmpTrack).dz());
//     // //     double edeltaz= sqrt(edz_tr*edz_tr + edz_ele*edz_ele);
//     // //     //if (deltaz!= fabs((*itrTr).vz()  - (*tmpTrack).vz()))	cout<<"dz= "<<deltaz<<" vz= "<<fabs((*itrTr).vz()  - (*tmpTrack).vz())<<endl;
//     // //     int nhits = (*itrTr).recHitsSize();
//     // //     if (tracker_track_ && !track_quality_tt(d0, ed0, deltaz, edeltaz, nhits,this_pt) ) continue;
//     // //     if (!tracker_track_ && !track_quality(this_pt, deltaz)) continue;
    
    //double deltadr=sqrt((itrTr->p4().eta()-thisElectronMomentumAtVtx.eta())*(itrTr->p4().eta()-thisElectronMomentumAtVtx.eta())+(itrTr->p4().phi()-thisElectronMomentumAtVtx.phi())*(itrTr->p4().phi()-thisElectronMomentumAtVtx.phi()) );
    double deltadr = DeltaR(tmpTrackMomentumAtVtxMu,thisElectronMomentumAtVtx) ;
    if ( fabs(deltadr) < extRadius_ && 
 	 fabs(deltadr) > intRadius_ ){
      cout << "Found a muon in the isolation cone of electrons...subtracting the pT= " << thismuon_pt << endl;
      dummypT -= thismuon_pt ;
    }
  }
  

  // correct for close pairs from Z(Z*/gamma* decay)
  // idea: look for an electron candidate of oposite sign in the cone, substract its pT from the sum if 
  // invariant mass with original electron is greater then a threshold
  
  double invMassThr = 12. ;
  //    double tmpElPt1=0.,tmpElPt2=0.,drTrackMatch = 0.;
  
  int iloop=0,iloop2=0;
  //begin loop over electrons
  //    for (reco::CandidateCollection::const_iterator elecand = electronCollection->begin(); 
  // elecand != electronCollection->end(); ++elecand ) {
  
  //const Candidate * cand2;
  //cand2 = &* elecand->masterClone();
  //const reco::GsfElectron * ele = dynamic_cast<const reco::GsfElectron *>(cand2);
  
  for ( edm::View<reco::GsfElectron>::const_iterator ele  = electronCollection_->begin(); 
	ele != electronCollection_->end(); 
	++ele ) {
    //cout <<"electron number: "<<iloop<<endl;
    bool elveto=false;
    reco::GsfTrackRef tmpTrack = ele->gsfTrack() ;
    math::XYZVector tmpElectronMomentumAtVtx = (*tmpTrack).innerMomentum () ; 
    
    if (thisElectronMomentumAtVtx == tmpElectronMomentumAtVtx) {
      //   cout<<"same electron"<<endl;
      iloop++;
      continue ;
    }
    
    double dr = DeltaR(thisElectronMomentumAtVtx,tmpElectronMomentumAtVtx);
    //   cout<<"dr= "<<dr<<endl;
    if ( fabs(dr) < extRadius_  && 
	 fabs(dr) > intRadius_  ) { 
      // cout<<"inside cone"<<endl;
      math::XYZTLorentzVector thisElectronLorentzVector = electron_->p4() ;
      math::XYZTLorentzVector tmpElectronLorentzVector  = ele->p4() ;
      math::XYZTLorentzVector twoElectrons ;
      twoElectrons= thisElectronLorentzVector + tmpElectronLorentzVector ;
      if (ele->charge() != electron_->charge() && (twoElectrons.mass() > invMassThr)){
	elveto= true;
	//	cout<<"e vetoed Z*"<<endl;
      }
      
      float massmin=100;
      
      //	for (reco::CandidateCollection::const_iterator elecand = electronCollection->begin(); 
      //   elecand != electronCollection->end(); ++elecand ) {
      
      //const Candidate * cand3;
      //cand3 = &* elecand->masterClone();
      //const reco::GsfElectron * oele = dynamic_cast<const reco::GsfElectron *>(cand3);
      
      for(edm::View<reco::GsfElectron>::const_iterator oele  = electronCollection_->begin(); 
	  oele != electronCollection_->end(); ++oele ) {
	cout<<"checking electrons outside cone : "<<iloop2<<endl;
	math::XYZTLorentzVector otherElectronLorentzVector = oele->p4() ;
	if ((oele->gsfTrack()->innerMomentum() != tmpElectronMomentumAtVtx) 
	    && (oele->gsfTrack()->innerMomentum() != thisElectronMomentumAtVtx)){
	  if (oele->charge() != ele->charge()){
	    twoElectrons = tmpElectronLorentzVector + otherElectronLorentzVector;
	    float mass= twoElectrons.mass();
	    // cout<<"mass= "<<mass<<endl;
	    if (abs(mass-ZMASSPDG_) <abs( massmin))
	      massmin= mass-ZMASSPDG_;
	  }
	  if (oele->charge() != electron_->charge()){
	    twoElectrons = thisElectronLorentzVector + otherElectronLorentzVector;
	    float mass= twoElectrons.mass();
	    //cout<<"mass= "<<mass<<endl;
	    if (abs(mass-ZMASSPDG_) <abs( massmin))
	      massmin= mass-ZMASSPDG_;
	  }
	}
	iloop2++;
      }
      // cout<<"massmin= "<<massmin<<endl;
      //cout<<"zmassinf= "<<Zmass_inf_<<" Zmasssup= "<<Zmass_sup_<<endl;
      if (massmin> Zmass_inf_ && massmin <Zmass_sup_) elveto=true;
      
      // if (elveto)cout<<"this electron is vetoed"<<endl;
      //else cout<<"this electron is not vetoed"<<endl;
      // make a loop on found KF tracks and find if any of them matches this electron
      //	for (reco::CandidateCollection::const_iterator trackcand = trackCollection->begin(); 
      //   trackcand != trackCollection->end(); ++trackcand ) {
      
      //const Candidate * cand;
      //cand = &* trackcand->masterClone();
      //const reco::Track * itrTr = dynamic_cast<const reco::Track *>(cand);      

      for ( edm::View<reco::Track>::const_iterator itrTr  = (*trackCollection_).begin() ; 
	    itrTr != (*trackCollection_).end()   ; 
	    ++itrTr) 
	{
	  double tmpElPt1=0.,tmpElPt2=0.,drTrackMatch = 0.;
	  math::XYZVector thisTrackMomentumAtVtx = itrTr->innerMomentum () ;
	  drTrackMatch = DeltaR(tmpElectronMomentumAtVtx,thisTrackMomentumAtVtx);
	  
	  if ( drTrackMatch < intRadius_ ) 
	    {
	      tmpElPt1=sqrt( thisTrackMomentumAtVtx.Perp2 () );
	      // cout<<"all ele pt= "<<tmpElPt1<<endl;
	      if (elveto) {
		tmpElPt2=sqrt( thisTrackMomentumAtVtx.Perp2 () );
		//	cout<<"veto ele pt= "<<tmpElPt2<<endl;
	      }
	      break;	     
	    }
	  // else cout<<"no KF track match"<<endl;
	}//end loop KF tracks
    }//end if in internal cone
    iloop++;
  }//end loop over electrons
  
  //  }
  
  //cout<<"tmpelpt1 (all)= "<<tmpElPt1<<" tmpelpt2 (Z)= "<<tmpElPt2<<endl;
  dummypT1=dummypT;
  dummypT2= dummypT;
  if (tmpElPt1 != 0.) 
    if (tmpElPt1 <= dummypT) 
      dummypT1=  dummypT - tmpElPt1; 
  
  if (tmpElPt2 != 0.) 
    if (tmpElPt2 <= dummypT) 
      dummypT2=  dummypT - tmpElPt2;
  
  
  //else --->FIXME 
  
  dummypTout.push_back(dummypT1);
  dummypTout.push_back(dummypT2);

  return dummypTout;

}



bool CandidateTkIsolation::track_quality_tt(double d0, double ed0,double dz, double edz, int nhits, double pt) const
{

  //  cout<<"track quality d0 = "<< d0<< " err d0 = "<< ed0<<" dz= "<<dz<<" errdz = "<<edz<<" nb hits= "<<nhits<<" pt= "<<pt<<endl;
  if (pt < 1) return false;
  if ( nhits < 8 ) {
    if ( d0 > 0.04 ) return false;
    if ( dz > 0.50 ) return false;
    if ( d0 / ed0 > 7 ) return false;
    if ( dz / edz > 7 ) return false;
  }
  else if ( nhits < 10 ) {
    if ( d0 > 0.20 ) return false;
    if ( dz > 2.00 ) return false;
    if ( d0 / ed0 > 10. ) return false;
    if ( dz / edz > 10. ) return false;
  }
  else {
    if ( d0 > 1.00 ) return false;
    if ( dz > 5.00 ) return false;
  }

  //    cout<<"quality true"<<endl;
  return true ;
}

bool CandidateTkIsolation::track_quality (double pt , double dz) const
{
  if (pt < ptLow_) return false;
  if (dz > lip_) return false;
  return true;

}



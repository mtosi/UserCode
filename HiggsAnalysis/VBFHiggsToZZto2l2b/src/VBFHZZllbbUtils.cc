#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/VBFHZZllbbUtils.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include <iostream>
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"

namespace vbfhzz2l2b
{

  
  // calculate the distance between two lorentz vectors 
  // using DeltaR(eta, phi) or normal space angle(theta, phi)
  double distance(const math::XYZTLorentzVector& v1, const math::XYZTLorentzVector& v2, bool useDeltaR = true) {
    if(useDeltaR) return ROOT::Math::VectorUtil::DeltaR(v1, v2);
    return ROOT::Math::VectorUtil::Angle(v1, v2);
  }
  

  int bTaggerCode ( const std::string& bTagger ) {
    
    int code = -1;
    if ( bTagger == "HIGHEFF" )            code = HIGHEFF;
    else if ( bTagger == "HIGHPUR"       ) code = HIGHPUR;
    else if ( bTagger == "COMBSECVTX"    ) code = COMBSECVTX;
    else if ( bTagger == "COMBSECVTXMVA" ) code = COMBSECVTXMVA;
    else if ( bTagger == "SOFTMUON"      ) code = SOFTMUON;
    else if ( bTagger == "SOFTELECTRON"  ) code = SOFTELECTRON;
    else if ( bTagger == "JETPROB"       ) code = JETPROB;
    else
      std::cout << "[VBFHZZllbbUtils::bTaggerCode] --> WARNING: bTagger " << bTagger << " NOT IMPLEMENTED!" << std::endl;
    
    return code;
 }

  double resolution ( double & recValue, double & refValue ) {
    return (recValue-refValue)/refValue;
  }
  
  
  void setVertex (TVector3 &myvector, 
		  const double& x, const double& y, const double& z) {
    myvector.SetX (x);
    myvector.SetY (y);
    myvector.SetZ (z);
  }
  
  void setVertex (TVector3 &myvector, 
		  const math::XYZPoint & mom) {
    myvector.SetX (mom.X());
    myvector.SetY (mom.Y());
    myvector.SetZ (mom.Z());
  }

  void setVertex (TVector3 &myvector, 
		  const TVector3 & mom) {
    myvector.SetX (mom.X());
    myvector.SetY (mom.Y());
    myvector.SetZ (mom.Z());
  }

  void setVertex (TVector3 &myvector, 
		  const reco::Vertex & mom) {
    myvector.SetX (mom.x());
    myvector.SetY (mom.y());
    myvector.SetZ (mom.z());
  }

  void setMomentum (TLorentzVector & myvector, 
		    const reco::Candidate & gen) {
    myvector.SetPx (gen.px ()) ;
    myvector.SetPy (gen.py ()) ;
    myvector.SetPz (gen.pz ()) ;
    myvector.SetE (gen.energy ()) ;
  }

  void setMomentum (TLorentzVector & v1, 
		    const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& v2) {
    v1.SetPx ( v2.Px() );
    v1.SetPy ( v2.Py() );
    v1.SetPz ( v2.Pz() );
    v1.SetE  ( v2.E()  );
  }


  double tracksInvariantMass( const reco::TrackRefVector & selectedTracks ) {
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


  bool BhadronTable(int pdgcode) {
    bool isBhadron = false;
    int Bmeson[53] = {511,521,10511,10521,513,523,10513,10523,20513,20523,515,
		      525,531,10531,533,10533,20533,535,541,10541,543,10543,
		      20543,545,551,10551,100551,110551,200551,210551,553,10553,
		      20553,30553,100553,110553,120553,130553,200553,210553,
		      220553,300553,9000553,9010553,555,10555,20555,100555,
		      110555,120555,200555,557,10557};
    
    int Bbaryon[35] = {5122,5112,5212,5222,5114,5214,5224,5132,5232,5312,5322,
		       5314,5324,5332,5334,5142,5242,5412,5422,5414,5424,5342,
		       5432,5434,5442,5444,5512,5522,5514,5524,5532,5534,5542,
		       5544,5554};
    
    for(int i=0;i<53;i++){
      if(abs(pdgcode) == Bmeson[i]){
	isBhadron = true;
	return isBhadron;
      }
    }
    for(int i=0;i<35;i++){
      if(abs(pdgcode) == Bbaryon[i]){
	isBhadron = true;
	return isBhadron;
      }
    }
    return isBhadron;
  }

  
  double dzVtxMomentum (const math::XYZVector & vertex,
			const math::XYZVector & momentum) { 
    double invpt = 1. / sqrt ( momentum.Perp2 () ) ;
    return vertex.z () -  (vertex.x () * momentum.x () + vertex.y () * momentum.y ()) * invpt * (momentum.z () * invpt) ; 
  }


  bool testTrackerTrack (reco::TrackCollection::const_iterator & track_itr) {
    
    // extract track properties
    double d0    = track_itr->d0 () ; 
    double d0Err = track_itr->d0Error () ;
    double d0Sig = fabs(d0)/d0Err;
    double dz    = track_itr->dz () ;
    double dzErr = track_itr->dzError () ;
    double dzSig = fabs(dz)/dzErr;
    int    nhits = track_itr->recHitsSize () ;

    if ( nhits < 8 ) {
      if ( fabs(d0) >  0.04 ) return false;
      if ( fabs(dz) >  0.50 ) return false;
      if ( d0Sig    >  7.00 ) return false;
      if ( dzSig    > 10.00 ) return false;
    }
    else if ( nhits < 10 ) {
      if ( fabs(d0) >  0.20 ) return false;
      if ( fabs(dz) >  2.00 ) return false;
      if ( d0Sig    > 10.00 ) return false;
      if ( dzSig    > 10.00 ) return false;
    }
    else {
      if ( fabs(d0) > 1.00 ) return false;
      if ( fabs(dz) > 5.00 ) return false;
    }
    
    return true;
    
  }

  bool testTrackerTrack (reco::TrackCollection::const_iterator & track_itr, const reco::Muon* muon) {

  // Extract properties at Vertex
  float vz  = track_itr -> vz();
  float edz = track_itr -> dzError();
  float d0  = track_itr -> d0();
  float ed0 = track_itr -> d0Error();
  // Difference with lepton/track Z vertex
  float dz  = muon -> vertex().z() - vz;
  
  dz  = fabs(dz);   //impact parameter in the (r, z) plane
  edz = fabs(edz);  //error on dz 
  d0  = fabs(d0);   //impact parameter in the (r, phi) plane
  ed0 = fabs(ed0);  //error on d0   

  reco::Particle::LorentzVector track(track_itr -> px(), track_itr -> py(), track_itr -> pz(), track_itr -> p());
  float track_pt =  track.pt();
  int nhits = track_itr -> recHitsSize(); 

  if ( nhits < 8 ) {
    if ( track_pt < 1.00) return false;
    if ( d0       > 0.04) return false;
    if ( dz       > 0.50) return false;
    if ( d0 / ed0 > 7.00) return false;
    if ( dz / edz > 7.00) return false;
  }
  else if ( nhits < 10 ) {
    if ( track_pt <  1.00) return false;
    if ( d0       >  0.20) return false;
    if ( dz       >  2.00) return false;
    if ( d0 / ed0 > 10.00) return false;
    if ( dz / edz > 10.00) return false;
  }
  else {
    if ( track_pt < 1.) return false;
    if ( d0 > 1.00) return false;
    if ( dz > 5.00) return false;
  }
  return true;

}


}


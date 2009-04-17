/**\class VBFHZZllbbMuonIsolationProducer
 *
 * based on the code provided by A. Graziano, A. Drotzevsky
 * 
 * Compute isolation for cones around muon candidates
 */

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbMuonIsolationProducer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Muons
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h" 
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "RecoMuon/TrackingTools/interface/MuonPatternRecoDumper.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/VBFHZZllbbMuonAssociationMap.h"

// Tracker tracks
#include "DataFormats/TrackReco/interface/Track.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/CandidateTkIsolation.h"

// Electrons
#include <DataFormats/EgammaCandidates/interface/GsfElectron.h>
#include <DataFormats/EgammaCandidates/interface/GsfElectronFwd.h>
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

//
#include <memory>
#include <vector>

using namespace edm;
using namespace std;
using namespace reco;

// constructor
VBFHZZllbbMuonIsolationProducer::VBFHZZllbbMuonIsolationProducer(const edm::ParameterSet& pset) {
  
  theECALIsoDepositLabel_    = pset.getParameter<edm::InputTag>("ECALIsoDepositLabel");        
  theHCALIsoDepositLabel_    = pset.getParameter<edm::InputTag>("HCALIsoDepositLabel");        
  theHOCALIsoDepositLabel_   = pset.getParameter<edm::InputTag>("HOCALIsoDepositLabel");      
  theTrackerIsoDepositLabel_ = pset.getParameter<edm::InputTag>("TrackerIsoDepositLabel"); 
  muonLabel_     = pset.getParameter<edm::InputTag>("MuonLabel"    );
  electronLabel_ = pset.getParameter<edm::InputTag>("ElectronLabel");
  trackLabel_    = pset.getParameter<edm::InputTag>("TrackLabel"   );
  mainConeSize_ = pset.getParameter<double>("isolationCone"    );
  vetoConeSize_ = pset.getParameter<double>("isolationConeVeto");
  isoCut_       = pset.getParameter<double>("isolationCut"     );

  string alias;
  string iName = "MuonIsolation";

  produces<vector<double> >( alias = iName + "X" ).setBranchAlias( alias );
  produces<vector<double> >( alias = iName + "CalIso" ).setBranchAlias( alias );
  produces<vector<double> >( alias = iName + "ECalIso" ).setBranchAlias( alias );
  produces<vector<double> >( alias = iName + "HCalIso" ).setBranchAlias( alias );
  produces<vector<double> >( alias = iName + "SumpT" ).setBranchAlias( alias );
  produces<reco::MuonCollection>();
  produces<VBFHZZllbbMuonAssociationMap>();
}


// destructor
VBFHZZllbbMuonIsolationProducer::~VBFHZZllbbMuonIsolationProducer() {

}


void VBFHZZllbbMuonIsolationProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  auto_ptr<vector<double> > XAfterVetos     ( new vector<double> );
  auto_ptr<vector<double> > CalIsoAfterVetos( new vector<double> );
  auto_ptr<vector<double> > EcalAfterVetos  ( new vector<double> );
  auto_ptr<vector<double> > HcalAfterVetos  ( new vector<double> );
  auto_ptr<vector<double> > sumPtAfterVetos ( new vector<double> );
  auto_ptr<reco::MuonCollection> Isolmuon( new reco::MuonCollection );
  auto_ptr<VBFHZZllbbMuonAssociationMap> IsoMuMap( new VBFHZZllbbMuonAssociationMap() );
  VBFHZZllbbMuonAssociationMap::Filler filler( *IsoMuMap );

  // Get tracks
  Handle<edm::View<reco::Track> > trackHandle;
  iEvent.getByLabel(trackLabel_, trackHandle);

  //get electrons
  Handle<edm::View<reco::GsfElectron> > electronHandle;
  iEvent.getByLabel(electronLabel_, electronHandle);

  // Get muons used to build a map
  Handle<edm::View<reco::Muon> > muonHandle;
  iEvent.getByLabel(muonLabel_, muonHandle);

  vector<const reco::Muon*> muons;   

  edm::View<reco::Muon>::const_iterator muon_itr = muonHandle->begin();  
  unsigned int i=0;
  for( ; muon_itr != muonHandle->end(); ++muon_itr){
    muons.push_back(&(*muon_itr));
    i++;
  }

  //**************** MAPS ******************
  //to get the deposits inside an isolation cone
  
  Handle<reco::IsoDepositMap> ecalIsoHandle;    
  iEvent.getByLabel(theECALIsoDepositLabel_, ecalIsoHandle);

  Handle<reco::IsoDepositMap> hcalIsoHandle;
  iEvent.getByLabel(theHCALIsoDepositLabel_, hcalIsoHandle);
  
  Handle<reco::IsoDepositMap> hocalIsoHandle;
  iEvent.getByLabel(theHOCALIsoDepositLabel_, hocalIsoHandle);
  
  Handle<reco::IsoDepositMap> trackerIsoHandle;
  iEvent.getByLabel(theTrackerIsoDepositLabel_, trackerIsoHandle);

  //************** weights ******************

  double muon_muondeltaR = 0.;
  double depVeto03Trk = 0., depVeto03Ecal = 0., depVeto03Hcal = 0., depVeto03HOcal = 0., depVeto03CalIso = 0., depVeto03X = 0.;
  double TrkIso = 0., CalIso = 0., X1 = 0., X2 = 0., X3 = 0., X4 = 0., X5 = 0.;
  
  vector<double> sumPtOverMuPtAfterVetos, HOcalAfterVetos; 
  vector<double> X1Vec, X2Vec, X3Vec, X4Vec, X5Vec;

  unsigned int index=0;

  size_t n = muonHandle->size();
  std::vector<float> iso(n);

  for ( unsigned int h = 0 ; h < muonHandle->size() ; h++ ){ 

    cout << "Candidate muon found for tight isolation" << endl;
    
    //********* veto against the other muons in an isolation cone ************
    edm::Ref<edm::View<reco::Muon> > muonRef(muonHandle,h);
    
    reco::IsoDeposit depTracker((*trackerIsoHandle)[muonRef]); // get sumPt around h-th muon
    reco::IsoDeposit depEcal   ((*ecalIsoHandle)   [muonRef]); // get ECal dep around h-th muon
    reco::IsoDeposit depHcal   ((*hcalIsoHandle)   [muonRef]); // get HCal dep around h-th muon
    reco::IsoDeposit depHOcal =((*hocalIsoHandle)  [muonRef]); // get HoCal dep around h-th muon
    
    IsoDeposit::Vetos vetoTrkVector;                             //vector of vetos    
    reco::IsoDeposit::Direction dirMuon1 = depTracker.direction(); //direction of the muon you're isolating
    
    for(unsigned int k = 0; k < muonHandle->size(); k++) {    //inner loop over muons
      edm::Ref<edm::View<reco::Muon> > muonRefbis(muonHandle,k);

      if(k == h) continue; //skip same muons

      reco::IsoDeposit depTracker2((*trackerIsoHandle)[muonRefbis]); //get sumPt around h-th muon

      reco::IsoDeposit::Direction dirMuon2 = depTracker2.direction();                     
      
      muon_muondeltaR = dirMuon1.deltaR(dirMuon2);         //calculate delta_R
      
      if(muon_muondeltaR < mainConeSize_) {             //i.e. if another muon falls into the cone
	
	IsoDeposit::Veto vetoTrk(dirMuon2, vetoConeSize_);
	vetoTrkVector.push_back(vetoTrk);   //fill veto vector for tracker
      }
    }

    edm::View<reco::GsfElectron>::const_iterator electron_itr = electronHandle->begin();
    for ( ; electron_itr != electronHandle->end(); ++electron_itr) {

      reco::Particle::LorentzVector p4Trackelectron(electron_itr ->px(), electron_itr ->py(), electron_itr -> pz(), electron_itr -> p());      
      reco::IsoDeposit::Direction dirTrackelectron(p4Trackelectron.eta(), p4Trackelectron.phi());

      float muon_electrondeltaR = dirMuon1.deltaR(dirTrackelectron);
      
      if(muon_electrondeltaR < mainConeSize_ ) {
	
	IsoDeposit::Veto vetoTrk3(dirTrackelectron, vetoConeSize_);
	vetoTrkVector.push_back(vetoTrk3);   //fill veto vector for tracker	
      }
    }
    
    edm::View<reco::Track>::const_iterator track_itr = trackHandle->begin();
    for ( ; track_itr != trackHandle->end(); track_itr++) {    
      
      bool goodTrack = testTrackerTrack(track_itr, muons[h]);
      
      reco::Particle::LorentzVector p4Track(track_itr -> px(), track_itr -> py(), track_itr -> pz(), track_itr -> p());      
      reco::IsoDeposit::Direction dirTrack(p4Track.eta(), p4Track.phi());
      
      float muon_trackdeltaR = dirMuon1.deltaR(dirTrack);
      
      if(muon_trackdeltaR < mainConeSize_ && goodTrack == false) {
	
	IsoDeposit::Veto vetoTrk2(dirTrack, vetoConeSize_);
	vetoTrkVector.push_back(vetoTrk2);   //fill veto vector for tracker	
      }
    }


    //****** redefining all the isolation variables after vetos
    depVeto03Trk    = depTracker.depositWithin(mainConeSize_, vetoTrkVector, false); // sumPt
    depVeto03Ecal   = depEcal.depositWithin(   mainConeSize_, vetoTrkVector, false); // ECal
    depVeto03Hcal   = depHcal.depositWithin(   mainConeSize_, vetoTrkVector, false); // HCal
    depVeto03HOcal  = depHOcal.depositWithin(  mainConeSize_, vetoTrkVector, false); // HOCal
    depVeto03CalIso = 1.5 * depVeto03Ecal + depVeto03Hcal;                           // CalIso
    depVeto03X      = 2.0 * depVeto03Trk + depVeto03CalIso;                          // X
    
    TrkIso = depVeto03Trk;
    CalIso = depVeto03CalIso;
    
    X1 = TrkIso + CalIso;                                //linear1
    X2 = 1.5 * TrkIso + 0.5 * CalIso;                    //linear2
    X3 = TrkIso * CalIso;                                //hyperbolic
    X4 = sqrt( (TrkIso*TrkIso) + (CalIso*CalIso));       //circular
    X5 = sqrt( 4*(TrkIso*TrkIso) + (CalIso*CalIso) );    //elliptical
    
    //double ptmu = Muons[h] -> pt();   //pT of the muon the cone is drawn around
    double ptmu = muonHandle->at(h).pt();  
    if(ptmu != 0.) {                  //do not divide by zero!
      double sumptovermupt = depVeto03Trk/ptmu;                   //sumPt/pT_mu
      sumPtOverMuPtAfterVetos.push_back(sumptovermupt);
    }
    
    cout << "X variable:" << depVeto03X << endl;
    if ( depVeto03X < isoCut_ ){
      sumPtAfterVetos->push_back(depVeto03Trk);        //sumPt
      EcalAfterVetos->push_back(depVeto03Ecal);        //ECal
      HcalAfterVetos->push_back(depVeto03Hcal);        //HCal
      HOcalAfterVetos.push_back(depVeto03HOcal);       //HOCal
      CalIsoAfterVetos->push_back(depVeto03CalIso);    //CalIso
      XAfterVetos->push_back(depVeto03X);              //X
      X1Vec.push_back(X1);                             //X1
      X2Vec.push_back(X2);                             //X2
      X3Vec.push_back(X3);                             //X3
      X4Vec.push_back(X4);                             //X4
      X5Vec.push_back(X5);                             //X5
      
      Isolmuon->push_back(muonHandle->at(h));
      //edm::Ref<edm::View<reco::Muon> > mutrackref(muonHandle,index);
      //edm::RefToBase<reco::Muon> mutrackref = muonHandle->refAt(index);      
      //IsoMuMap->insert(mutrackref,float(depVeto03X));
      iso[index] = float(depVeto03X);
    }
    index++;
  }

  // filling map
  cout << "Size of iso=" << iso.size() << endl;
  filler.insert(muonHandle, iso.begin(), iso.end());
  filler.fill();  

  //
  // sort(XAfterVetos->begin(),XAfterVetos->end());
  //   sort(CalIsoAfterVetos->begin(),CalIsoAfterVetos->end());
  //   sort(EcalAfterVetos->begin(),EcalAfterVetos->end());
  //   sort(HcalAfterVetos->begin(),HcalAfterVetos->end());
  //   sort(sumPtAfterVetos->begin(),sumPtAfterVetos->end());
  

  const string & isoName = "MuonIsolation";
  iEvent.put( XAfterVetos,      isoName + "X" );
  iEvent.put( CalIsoAfterVetos, isoName + "CalIso" );
  iEvent.put( EcalAfterVetos,   isoName + "ECalIso" );  
  iEvent.put( HcalAfterVetos,   isoName + "HCalIso" );
  iEvent.put( sumPtAfterVetos,  isoName + "SumpT" );
  const string iName = "";
  iEvent.put( Isolmuon, iName );
  iEvent.put( IsoMuMap, iName );

}


bool VBFHZZllbbMuonIsolationProducer::testTrackerTrack(edm::View<reco::Track>::const_iterator iterTrack, const reco::Muon* muon) {

  /************* quality track requirements ************/  

  bool keep = true;
  // Extract properties at Vertex
  float vz  = iterTrack -> vz();
  float edz = iterTrack -> dzError();
  float d0  = iterTrack -> d0();
  float ed0 = iterTrack -> d0Error();
  // Difference with lepton/Z vertex
  float dz  = muon -> vertex().z() - vz;
  
  dz  = fabs(dz);                //impact parameter in the (r, z) plane
  edz = fabs(edz);               //error on dz 
  d0  = fabs(d0);                //impact parameter in the (r, phi) plane
  ed0 = fabs(ed0);               //error on d0   

  reco::Particle::LorentzVector track(iterTrack -> px(), iterTrack -> py(), iterTrack -> pz(), iterTrack -> p());
  float track_pt =  track.pt();
  int nhits = iterTrack -> recHitsSize(); 

  if ( nhits < 8 ) {
    if ( track_pt < 1.) return false;
    if ( d0 > 0.04) return false;
    if ( dz > 0.50) return false;
    if ( d0 / ed0 > 7.) return false;
    if ( dz / edz > 7.) return false;
  }
  else if ( nhits < 10 ) {
    if ( track_pt < 1.) return false;
    if ( d0 > 0.20) return false;
    if ( dz > 2.00) return false;
    if ( d0 / ed0 > 10.) return false;
    if ( dz / edz > 10.) return false;
  }
  else {
    if ( track_pt < 1.) return false;
    if ( d0 > 1.00) return false;
    if ( dz > 5.00) return false;
  }
  return keep;
}

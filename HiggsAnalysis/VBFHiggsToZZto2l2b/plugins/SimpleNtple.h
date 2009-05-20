// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "AnalysisDataFormats/Egamma/interface/ElectronID.h"
#include "AnalysisDataFormats/Egamma/interface/ElectronIDAssociation.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackExtraBase.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"


#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

#include "TParticle.h"

//
// class decleration
//

class SimpleNtple : public edm::EDAnalyzer {
 public:
  explicit SimpleNtple(const edm::ParameterSet&);
  ~SimpleNtple();

  typedef edm::View<reco::Track> trackCollection ;
  typedef math::XYZTLorentzVector LorentzVector ;
  
 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  void InitObjs();
  void FillEvent        (const edm::Event&, const edm::EventSetup&);
  void FillElectron     (const edm::Event&, const edm::EventSetup&);
  void FillMuon         (const edm::Event&, const edm::EventSetup&);
  void FillMet          (const edm::Event&, const edm::EventSetup&);
  void FillTagJet       (const edm::Event&, const edm::EventSetup&);
  void FillZhad         (const edm::Event&, const edm::EventSetup&);
  void FillZlep         (const edm::Event&, const edm::EventSetup&);
  void FillcorIC5CaloJetsWithBTag (const edm::Event&, const edm::EventSetup&); 
  void FillcorIC5PFJetsWithBTag   (const edm::Event&, const edm::EventSetup&);   
  void FillTrack        (const edm::Event&, const edm::EventSetup&);
  void FillGenParticle  (const edm::Event&, const edm::EventSetup&);
  void FillGenJet       (const edm::Event&, const edm::EventSetup&);
  void FillGenMet       (const edm::Event&, const edm::EventSetup&);

  void setVertex (TVector3 &, const TVector3 &);
  // ----------member data ---------------------------
  TTree *mytree_;

  // event
  int evtID_;
  int evtRun_, evtEvent_;
  std::vector<double> * tagjetInvMass_;    // depends on the tag jet definition
  std::vector<double> * tagjetDeltaEta_;
  std::vector<double> * tagjetZeppenfeld_;
  std::vector<double> * zjetInvMass_;    // depends on the z jet definition => btagger?
  std::vector<double> * zjetDeltaEta_;
  std::vector<double> * zjetZeppenfeld_;
  //electrons;
  int eleN_;
  TClonesArray        * eleP4_ ;
  std::vector<double> * eleEt_;
  std::vector<double> * elePt_;
  std::vector<double> * eleIsoSumPt_;
  std::vector<double> * eleIsoNtrack_;
  std::vector<double> * eleD0_;
  std::vector<double> * eleDxy_;
  std::vector<double> * eleDxyError_;
  std::vector<int>    * eleID_;
  TClonesArray        * eleVtxP3_;
  //muons
  int muN_;
  TClonesArray        * muP4_ ;
  std::vector<double> * muEt_;
  std::vector<double> * muPt_;
  std::vector<double> * muIsoSumPt_;
  std::vector<double> * muIsoNtrack_;
  std::vector<double> * muNormChi2_;
  std::vector<double> * muD0_;
  std::vector<double> * muDxy_;
  std::vector<double> * muDxyError_;
  std::vector<int>    * muID_;
  TClonesArray        * muVtxP3_;
  // tag jets
  int tagjetN_;
  int tagjetNtrack_;
  TClonesArray        * tagjetP4_;
  std::vector<double> * tagjetEmFrac_;
  std::vector<double> * tagjetChFrac_;
  std::vector<double> * tagjetCorEt_;
  std::vector<double> * tagjetCorPt_;
  TClonesArray        * tagjetVtxP3_;
  // other jets with b tag
  int btagjetN_;
  int btagjetNtrack_;
  TClonesArray        * btagjetP4_;
  std::vector<double> * btagjetEmFrac_;
  std::vector<double> * btagjetChFrac_;
  std::vector<double> * btagjetCorEt_;
  std::vector<double> * btagjetCorPt_;
  std::vector<double> * btagjetCompoSVbTagDiscr_;
  std::vector<double> * btagjetHighEFFbTagDiscr_;
  std::vector<double> * btagjetHighPURbTagDiscr_;
  TClonesArray        * btagjetVtxP3_;
 

  // met
  TClonesArray        * metP4_ ;
  std::vector<double> * metSig_ ;

  TClonesArray * trackP4_ ;
  TClonesArray * genparticleP4_ ;
  TClonesArray * genjetP4_;
  TClonesArray * genmetP4_;

  // utils  
  TLorentzVector myvector_ ;
  TVector3       myvertex_ ;

  // input parameter and tag label
  int           whichSim_;
  edm::InputTag trackLabel_;
  edm::InputTag muonLabel_;
  edm::InputTag electronLabel_;
  edm::InputTag metLabel_;
  edm::InputTag tagJetLabel_;
  std::string   corIC5CaloJetsWithBTagLabel_;
  std::string   corIC5PFJetsWithBTagLabel_;
  bool          corIC5PFJetsWithBTagFlag_;

  edm::InputTag genParticleLabel_;
  edm::InputTag genJetLabel_;
  edm::InputTag genMetLabel_;

};

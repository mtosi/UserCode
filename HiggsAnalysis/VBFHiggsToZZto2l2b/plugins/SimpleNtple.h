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

  void Init();
  void FillEle          (const edm::Event&, const edm::EventSetup&);
  void FillMu           (const edm::Event&, const edm::EventSetup&);
  void FillMet          (const edm::Event&, const edm::EventSetup&);
  void FillTagJet       (const edm::Event&, const edm::EventSetup&);
  void FillJet          (const edm::Event&, const edm::EventSetup&);
  void FillCaloJet      (const edm::Handle<reco::CaloJetCollection>&, TClonesArray &);
  void FillPFJet        (const edm::Handle<reco::PFJetCollection>&, TClonesArray &);
  void FillIC5CaloJets            (const edm::Event&, const edm::EventSetup&);	       
  void FillIC5PFJets              (const edm::Event&, const edm::EventSetup&);	       
  void FillcorIC5CaloJetsWithBTag (const edm::Event&, const edm::EventSetup&); 
  void FillcorIC5PFJetsWithBTag   (const edm::Event&, const edm::EventSetup&);   
  void FillSC5CaloJets            (const edm::Event&, const edm::EventSetup&);	       
  void FillSC5PFJets              (const edm::Event&, const edm::EventSetup&);	       
  void FillcorSC5CaloJetsWithBTag (const edm::Event&, const edm::EventSetup&); 
  void FillcorSC5PFJetsWithBTag   (const edm::Event&, const edm::EventSetup&);   

  void FillTracks       (const edm::Event&, const edm::EventSetup&);
  void FillKindEvent    (const edm::Event&, const edm::EventSetup&);
  void FillGenParticles (const edm::Event&, const edm::EventSetup&);
  void FillGenJet       (const edm::Event&, const edm::EventSetup&);
  void FillGenMet       (const edm::Event&, const edm::EventSetup&);

  void setMomentum (TLorentzVector &myvector, const LorentzVector & mom) ;
  
  // ----------member data ---------------------------
  TTree *mytree_;
  //electrons;
  int nEle;
  float IsolEleSumPt[30],IsolEleNTracks[30];
  int EleId[30];
  //muons
  int nMu;
  int IdEvent;
  float IsolMuSumPt[30],IsolMuNTracks[30];
  // tag jets
  float MinvTags;
  //other jets
  std::vector<double> * emFrac_;
  // other jets with b tag
  std::vector<double> * emFracWithBTag_;
  std::vector<double> * corEtWithBTag_;
  std::vector<double> * compoSVbTagDiscrWithBTag_;
  std::vector<double> * highEFFbTagDiscrWithBTag_;
  std::vector<double> * highPURbTagDiscrWithBTag_;
  std::vector<std::vector<double> > * discriminatorVecWithBTag_;
 

  TClonesArray * m_tagJets ;
  TClonesArray * m_otherJets ;
  TClonesArray * m_otherJets_IterativeCone5CaloJets ;
  TClonesArray * m_otherJets_IterativeCone5PFJets;
  TClonesArray * m_otherJets_corIterativeCone5CaloJetsWithBTag ;
  TClonesArray * m_otherJets_corIterativeCone5PFJetsWithBtag;
  TClonesArray * m_otherJets_SisCone5CaloJets ;
  TClonesArray * m_otherJets_SisCone5PFJets;
  TClonesArray * m_otherJets_corSisCone5CaloJetsWithBTag ;
  TClonesArray * m_otherJets_corSisCone5PFJetsWithBtag;            
  TClonesArray * m_electrons ;
  TClonesArray * m_muons ;
  TClonesArray * m_MET ;
  TClonesArray * m_tracks ;
  TClonesArray * m_genParticles ;
  TClonesArray * m_genJets;
  TClonesArray * m_genMet;
  
  TLorentzVector myvector ;
  
  int           whichSim_;
  edm::InputTag TracksTag_;
  edm::InputTag EleTag_;
  edm::InputTag IsolEleTag_;
  edm::InputTag MuTag_;
  edm::InputTag IsolMuTag_;
  edm::InputTag MetTag_;
  edm::InputTag TagJetTag_;
  edm::InputTag JetTag_;
  bool          bool_IterativeCone5CaloJetsTag_;
  bool          bool_IterativeCone5PFJetsTag_;
  bool          bool_corIterativeCone5CaloJetsWithBTag_;
  bool          bool_corIterativeCone5PFJetsWithBTag_;
  bool          bool_SisCone5CaloJetsTag_;
  bool          bool_SisCone5PFJetsTag_;
  bool          bool_corSisCone5CaloJetsWithBTag_;
  bool          bool_corSisCone5PFJetsWithBTag_;
  edm::InputTag IterativeCone5CaloJetsTag_;
  edm::InputTag IterativeCone5PFJetsTag_;
  std::string   corIterativeCone5CaloJetsWithBTag_;
  std::string   corIterativeCone5PFJetsWithBTag_;
  edm::InputTag SisCone5CaloJetsTag_;
  edm::InputTag SisCone5PFJetsTag_;
  std::string   corSisCone5CaloJetsWithBTag_;
  std::string   corSisCone5PFJetsWithBTag_;

  edm::InputTag MCtruthTag_;
  edm::InputTag genJetTag_;
  edm::InputTag genMetTag_;

};

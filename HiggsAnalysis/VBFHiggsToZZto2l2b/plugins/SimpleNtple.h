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

 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  void InitObjs();
  void FillEvent           (const edm::Event&, const edm::EventSetup&);
  void FillElectron        (const edm::Event&, const edm::EventSetup&);
  void FillMuon            (const edm::Event&, const edm::EventSetup&);
  void FillMet             (const edm::Event&, const edm::EventSetup&);
  void FillTagJet          (const edm::Event&, const edm::EventSetup&);
  void FillZhad            (const edm::Event&, const edm::EventSetup&);
  void FillZlep            (const edm::Event&, const edm::EventSetup&);
  void FillcorIC5CaloJetsWithBTag (const edm::Event&, const edm::EventSetup&); 
  void FillcorIC5PFJetsWithBTag   (const edm::Event&, const edm::EventSetup&);   
  void FillTrack           (const edm::Event&, const edm::EventSetup&);
  void FillGenParticle     (const edm::Event&, const edm::EventSetup&);
  void FillGenJet          (const edm::Event&, const edm::EventSetup&);
  void FillGenMet          (const edm::Event&, const edm::EventSetup&);

  void setVertex (TVector3 &, const TVector3 &);
  // ----------member data ---------------------------
  TTree *mytree_;

  // ***************************************
  // **************** event ****************
  int evtID_;
  int evtRun_, evtEvent_;
  int jetN_;
  int btagjetN_;
  int invmasstagjetN_;
  int deltaetatagjetN_;
  int zeptagjetN_;
  int muN_;
  int glbmuN_;
  int glbmuPromptTightN_;
  int eleN_;
  std::vector<double> * invmasstagjetInvMass_;    // depends on the tag jet definition
  std::vector<double> * invmasstagjetDeltaEta_;
  std::vector<double> * invmasstagjetZeppenfeld_;
  std::vector<double> * deltaetatagjetInvMass_;    // depends on the tag jet definition
  std::vector<double> * deltaetatagjetDeltaEta_;
  std::vector<double> * deltaetatagjetZeppenfeld_;
  std::vector<double> * zeptagjetInvMass_;    // depends on the tag jet definition
  std::vector<double> * zeptagjetDeltaEta_;
  std::vector<double> * zeptagjetZeppenfeld_;
  std::vector<double> * zjetInvMass_;    // depends on the z jet definition => btagger?
  std::vector<double> * zjetDeltaEta_;
  std::vector<double> * zjetZeppenfeld_;
  // ***************************************

  // utils  
  TLorentzVector myvector_ ;
  TVector3       myvertex_ ;
  // isolation algorithm
  
  // input parameter and tag label
  int           whichSim_;
  edm::InputTag vertexLabel_;
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

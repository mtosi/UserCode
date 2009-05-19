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

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/SimpleNtpleObj.h"

#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

#include "TParticle.h"

//
// class decleration
//

using namespace vbfhzz2l2b;
using namespace vbfhzz2l2b::SimpleNtpleObj;

class newSimpleNtple : public edm::EDAnalyzer {
 public:
  explicit newSimpleNtple(const edm::ParameterSet&);
  ~newSimpleNtple();

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
  void FillTracks       (const edm::Event&, const edm::EventSetup&);
  void FillGenParticles (const edm::Event&, const edm::EventSetup&);
  void FillGenJet       (const edm::Event&, const edm::EventSetup&);
  void FillGenMet       (const edm::Event&, const edm::EventSetup&);

  // ----------member data ---------------------------
  TTree *mytree_;
  int evtID_;

  // tag jets
  std::vector<double> * invMassTagJet_;
  std::vector<double> * deltaEtaTagJet_;
  std::vector<double> * zeppenfeldTagJet_;


  TLorentzVector myvector_ ;
  
  // SimpleNtpleObj
  TClonesArray *CloneEvt_;      vbfhzz2l2b::SimpleNtpleObj::EVT*      evt_;
  TClonesArray *CloneJet_;      vbfhzz2l2b::SimpleNtpleObj::JET*      jet_;
  TClonesArray *CloneMuon_;     vbfhzz2l2b::SimpleNtpleObj::MUON*     muon_;
  TClonesArray *CloneElectron_; vbfhzz2l2b::SimpleNtpleObj::ELECTRON* electron_;
  TClonesArray *CloneZhad_;     vbfhzz2l2b::SimpleNtpleObj::ZHAD*     Zhad_;
  // not yet implemented in SimpleNtpleObj
  TClonesArray * m_tagJets;
  TClonesArray * m_MET;
  TClonesArray * m_tracks;
  TClonesArray * m_genParticles;
  TClonesArray * m_genJets;
  TClonesArray * m_genMet;
  

  // input parameter and tag label
  int           whichSim_;
  edm::InputTag tracksLabel_;
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

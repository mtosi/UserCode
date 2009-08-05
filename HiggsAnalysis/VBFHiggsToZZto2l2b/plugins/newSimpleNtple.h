// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/SimpleNtpleObj.h"

#include "TTree.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TClonesArray.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/VBFHZZllbbElectronTrackIsolationAlgos.h"

//
// class decleration
//

class newSimpleNtple : public edm::EDAnalyzer {
 public:
  explicit newSimpleNtple(const edm::ParameterSet&);
  ~newSimpleNtple();

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

  // ----------member data ---------------------------
  TTree *mytree_;

  int evtID_;
  // SimpleNtpleObj
  TClonesArray *CloneEvt_;      EVT*      evt_;
  TClonesArray *CloneJet_;      JET*      jet_;
  TClonesArray *CloneMuon_;     MUON*     muon_;
  TClonesArray *CloneElectron_; ELECTRON* electron_;
  TClonesArray *CloneZhad_;     ZHAD*     Zhad_;

  // not yet implemented in SimpleNtpleObj
  // ***************************************************************
  // **************** tag jets w/ inv mass criteria ****************
  int invmasstagjetN_;
  std::vector<int>    * invmasstagjetNtrack_;
  TClonesArray        * invmasstagjetP4_;
  std::vector<double> * invmasstagjetEmEnergyFraction_;
  std::vector<double> * invmasstagjetChFrac_;
  std::vector<double> * invmasstagjetCorEt_;
  std::vector<double> * invmasstagjetCorPt_;
  std::vector<double> * invmasstagjetCompoSVbTagDiscr_;
  std::vector<double> * invmasstagjetHighEFFbTagDiscr_;
  std::vector<double> * invmasstagjetHighPURbTagDiscr_;
  TClonesArray        * invmasstagjetPrimVtxP3_;
  // ****************************************************************
  // ****************************************************************
  // **************** tag jets w/ delta eta criteria ****************
  int deltaetatagjetN_;
  TClonesArray        * deltaetatagjetP4_;
  std::vector<double> * deltaetatagjetEmEnergyFraction_;
  std::vector<double> * deltaetatagjetChFrac_;
  std::vector<double> * deltaetatagjetCorEt_;
  std::vector<double> * deltaetatagjetCorPt_;
  std::vector<double> * deltaetatagjetCompoSVbTagDiscr_;
  std::vector<double> * deltaetatagjetHighEFFbTagDiscr_;
  std::vector<double> * deltaetatagjetHighPURbTagDiscr_;
  TClonesArray        * deltaetatagjetPrimVtxP3_;
  // ****************************************************************
  // *****************************************************************
  // **************** tag jets w/ zeppenfeld criteria ****************
  int zeptagjetN_;
  TClonesArray        * zeptagjetP4_;
  std::vector<double> * zeptagjetEmEnergyFraction_;
  std::vector<double> * zeptagjetChFrac_;
  std::vector<double> * zeptagjetCorEt_;
  std::vector<double> * zeptagjetCorPt_;
  std::vector<double> * zeptagjetCompoSVbTagDiscr_;
  std::vector<double> * zeptagjetHighEFFbTagDiscr_;
  std::vector<double> * zeptagjetHighPURbTagDiscr_;
  TClonesArray        * zeptagjetPrimVtxP3_;
  // *****************************************************************
  // ****************************************
  // **************** tracks ****************
  TClonesArray * trackP4_ ;
  // ****************************************
  // **********************************************
  // **************** gen particle ****************
  TClonesArray                      * genparticleP4_ ;
  TClonesArray                      * genparticlePrimVtxP3_;
  std::vector<int>                  * genparticlePdgID_;
  std::vector<int>                  * genparticleStatus_;
  std::vector<int>                  * genparticleIndex_;
  std::vector<int>                  * genparticleMomN_;
  std::vector< std::pair<int,int> > * genparticleMomPdgID_;
  //  std::vector< std::pair<int,int> > * genparticleMomPdgIndex_;
  std::vector<int>                  * genparticleKidN_;
  std::vector< std::pair<int,int> > * genparticleKidPdgID_;
  //  std::vector< std::pair<int,int> > * genparticleKidPdgIndex_;
  // **********************************************
  // *****************************************
  // **************** gen jet ****************
  TClonesArray                    * genjetP4_;
  TClonesArray                    * genjetPrimVtxP3_;
  std::vector<int>                * genjetKidN_;
  std::vector< std::vector<int> > * genjetKidPdgID_;

  // *****************************************

  // *****************************************
  // **************** gen met ****************
  TClonesArray * genmetP4_;
  TClonesArray * genmetPrimVtxP3_;
  // *****************************************

  // utils  
  TLorentzVector myvector_ ;
  TVector3       myvertex_ ;
  // isolation algorithm
  VBFHZZllbbElectronTrackIsolationAlgos * eleTrkIsoAlgo_;
  
  // input parameter and tag label
  int           whichSim_;

  edm::InputTag vertexLabel_;
  edm::InputTag trackLabel_;
  edm::InputTag muonLabel_;
  edm::InputTag electronLabel_;
  bool          eleTrkIsoAlgoFlag_;
  edm::InputTag metLabel_;
  edm::InputTag tagJetLabel_;
  std::string   corIC5CaloJetsWithBTagLabel_;
  std::string   corIC5PFJetsWithBTagLabel_;
  bool          corIC5PFJetsWithBTagFlag_;

  edm::InputTag genParticleLabel_;
  edm::InputTag genJetLabel_;
  edm::InputTag genMetLabel_;
};

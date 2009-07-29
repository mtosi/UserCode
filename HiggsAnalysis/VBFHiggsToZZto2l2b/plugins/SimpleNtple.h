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

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/VBFHZZllbbElectronTrackIsolationAlgos.h"
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

  // *******************************************
  // **************** electrons ****************
  int eleN_;
  TClonesArray        * eleP4_ ;
  std::vector<bool>   * eleLongLived_;
  std::vector<int>    * elePdgID_;
  std::vector<double> * eleCaloEnergy_;
  std::vector<double> * eleCaloEnergyError_;
  std::vector<double> * eleCaloEnergySig_;
  std::vector<double> * eleEnergy_;
  std::vector<double> * eleHadOverEm_;
  std::vector<double> * eleEt_;
  std::vector<double> * elePt_;
  std::vector<double> * eleIsoVal_;
  std::vector<double> * eleIsoSumPt_;
  std::vector<int>    * eleIsoNtrack_;
  std::vector<double> * eleD0_;
  std::vector<double> * eleD0Error_;
  std::vector<double> * eleD0Sig_;
  std::vector<double> * eleDxy_;
  std::vector<double> * eleDxyError_;
  std::vector<double> * eleDxySig_;
  std::vector<double> * eleDz_;
  std::vector<double> * eleDzError_;
  std::vector<double> * eleDzSig_;
  std::vector<double> * eleEtaError_;
  std::vector<double> * elePhiError_;
  std::vector<double> * eleNormChi2_;
  std::vector<double> * eleQoverP_;
  std::vector<double> * eleQoverPError_;
  std::vector<int>    * eleCharge_;
  TClonesArray        * elePrimVtxP3_;
  // *******************************************

  // **********************************************
  // **************** global muons ****************
  int muN_;
  int glbmuN_;
  int glbmuPromptTightN_;
  std::vector<bool>   * glbmuPromptTightFlag_;
  TClonesArray        * glbmuP4_ ;
  TClonesArray        * glbmuPrimVtxP3_;
  std::vector<int>    * glbmuCharge_;
  std::vector<int>    * glbmuPdgID_;
  // calorimeter info
  std::vector<double> * glbmuEmEnergy_;   
  std::vector<double> * glbmuEmS9Energy_; 
  std::vector<double> * glbmuHadEnergy_;  
  std::vector<double> * glbmuHadS9Energy_;
  std::vector<double> * glbmuHoEnergy_;   
  std::vector<double> * glbmuHoS9Energy_; 
  // tracking info
  std::vector<double> * glbmuIso03emEt_;
  std::vector<double> * glbmuIso03hadEt_;
  std::vector<double> * glbmuIso03hoEt_;
  std::vector<int>    * glbmuIso03nJets_;
  std::vector<int>    * glbmuIso03nTracks_;
  std::vector<double> * glbmuIso03sumPt_;
  std::vector<double> * glbmuIso05emEt_;
  std::vector<double> * glbmuIso05hadEt_;
  std::vector<double> * glbmuIso05hoEt_;
  std::vector<int>    * glbmuIso05nJets_;
  std::vector<int>    * glbmuIso05nTracks_;
  std::vector<double> * glbmuIso05sumPt_;
  std::vector<double> * glbmuChi2_;
  std::vector<double> * glbmuNdof_;
  std::vector<double> * glbmud0_;
  std::vector<double> * glbmud0Err_;
  std::vector<double> * glbmudz_;
  std::vector<double> * glbmudzErr_;
  // **********************************************

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

  // *******************************************************
  // **************** other jets with b tag ****************
  int jetN_;
  int btagjetN_;
  std::vector<bool>   * btagjetFlag_;
  TClonesArray        * btagjetP4_;
  TClonesArray        * btagjetPrimVtxP3_;
  TClonesArray        * btagjetSecVtxP3_;
  std::vector< std::vector<double> > * btagjetIP2d_;
  std::vector< std::vector<double> > * btagjetIP2dErr_;
  std::vector< std::vector<double> > * btagjetIP2dSig_;
  std::vector< std::vector<double> > * btagjetIP3d_;
  std::vector< std::vector<double> > * btagjetIP3dErr_;
  std::vector< std::vector<double> > * btagjetIP3dSig_;
  std::vector<double> * btagjetChFrac_;
  std::vector<double> * btagjetCorEt_;
  std::vector<double> * btagjetCorPt_;
  std::vector<double> * btagjetCompoSVbTagDiscr_;
  std::vector<double> * btagjetHighEFFbTagDiscr_;
  std::vector<double> * btagjetHighPURbTagDiscr_;
  std::vector<double> * btagjetEtaetaMoment_;
  std::vector<double> * btagjetPhiphiMoment_;
  std::vector<double> * btagjetEtaphiMoment_;
  std::vector<double> * btagjetMaxDistance_;   // maximum distance from jet to constituent 
  std::vector<int>    * btagjetNconstituents_; // number of constituents 
  std::vector<double> * btagjetPhysicsEtaQuick_; // (float fZVertex)
  std::vector<double> * btagjetPhysicsEtaDetailed_; // (float fZVertex)
  // Constituents getJetConstituents () const;
  // std::vector<const reco::Candidate*> getJetConstituentsQuick_; // quick list of constituents
  std::vector<double> * btagjetMaxEInEmTowers_;	
  std::vector<double> * btagjetMaxEInHadTowers_;
  std::vector<double> * btagjetEmEnergyFraction_;
  std::vector<int>    * btagjetNtrack_;
  // ******************************************************* 

  // ************************************* 
  // **************** met ****************
  TClonesArray        * metP4_ ;
  TClonesArray        * metPrimVtxP3_;
  std::vector<double> * metSig_ ;
  std::vector<double> * metSumEt_;
  // ************************************* 

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

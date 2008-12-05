// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

// Classes to be accessed
// ----------------------
#include "AnalysisExamples/AnalysisObjects/interface/OfflineMEt.h"
#include "AnalysisExamples/AnalysisObjects/interface/OfflineJet.h"
#include "AnalysisExamples/AnalysisObjects/interface/GlobalMuon.h"
#include "AnalysisExamples/AnalysisObjects/interface/SimpleElectron.h"
#include "AnalysisExamples/AnalysisObjects/interface/SimpleTau.h"
#include "AnalysisExamples/AnalysisObjects/interface/MCParticle.h"

// Root includes
// -------------
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TFile.h"
#include <memory>

// TMVAtreeMaker
#include "AnalysisExamples/AnalysisClasses/interface/TMVAtreeWriter.h"

// Data includes
// -------------
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

// GenJets
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

// Associator for the jets
#include "AnalysisExamples/AnalysisClasses/interface/Associator.h"
#include "AnalysisExamples/AnalysisClasses/interface/AssociatorEt.h"

// Needed to use the classes trained with the TMVA
#include "TMVA/Reader.h"

//
// class declaration
//
class VBFHZZkinematicsAnalyzer : public edm::EDAnalyzer {

 public:
  explicit VBFHZZkinematicsAnalyzer(const edm::ParameterSet&);
  ~VBFHZZkinematicsAnalyzer();

  //
  // constants, enums and typedefs
  //

//
// static data member definitions
//

//
// constructors and destructor
//

 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // signal-like event TTree
  auto_ptr<TMVAtreeWriter<Float_t> > tmvaSignalTreeWriterPtr_;
  // combinatorial-like event TTree
  auto_ptr<TMVAtreeWriter<Float_t> > tmvaCombinTreeWriterPtr_;
  // event variables vectors [value and name]
  vector<Float_t> eventSignalVariablesVector_;
  vector<Float_t> eventCombinVariablesVector_;
  vector<TString> eventVariablesNamesVector_;

  unsigned nZ_;
  unsigned nZleptons_;
  unsigned nZjets_;
  unsigned nforwardjets_;
  unsigned njets_;

  int eventcounter_;
  int goodLeptonicEventCounter_;
  int twoJetsAboveEtCutEventsCounter_;
  int twoJetsAboveEtCut_DeltaR2016EventsCounter_;
  int twoJetsAboveEtCut_DeltaR2CutEventsCounter_;
  int fourJetsAboveEtCutEventsCounter_;

  int nbin_;

  math::XYZTLorentzVector null_XYZTLorentzVector_;
  GlobalMuon null_globalMuon_;
  OfflineJet null_offlinejet_;

  // ----------member data ---------------------------
  edm::ParameterSet conf_;
  edm::InputTag offlineJetLabel_;
  //  edm::InputTag offlineMEtLabel_;
  edm::InputTag globalMuonLabel_;
  edm::InputTag simpleElectronLabel_;
  //  edm::InputTag simpleTauLabel_;
  //  edm::InputTag combinedSVBJetTagsLabel_;
  std::string   MCParticleLabel_;      
  //  edm::InputTag simVtxLabel_;
  double        leptonPtCut_;	          
  double        jetEtCut_;
  double        jetPartonDeltaR2Cut_;
  double        jetLeptonDeltaRCut_;
  double        jetEMfracCut_;
  std::string   tmvaSignalSuffix_;
  std::string   tmvaCombinSuffix_;

  bool writeTMVA_;

  int variablesNumber_;

  // Array of variables for the Reader
  Float_t * variables_;

  // The Reader object
  TMVA::Reader * reader_;

  TFile* OutputFile;
  // Use a dynamic construction,
  // or the TFile problem will crash the job
  // when moving from one input file to another.
  // The histograms must be created after the TFile is opened.

  std::vector< std::pair<int,int> > vecTypeIndexPair;

  TH1D * eventsNumber_;
  TH1D * goodLeptonicEventsNumber_;
  TH1D * twoJetsAboveEtCutEventsNumber_;
  TH1D * twoJetsAboveEtCut_DeltaR2016EventsNumber_;
  TH1D * twoJetsAboveEtCut_DeltaR2CutEventsNumber_;
  TH1D * fourJetsAboveEtCutEventsNumber_;

  TH1D * ZpartonEta_;
  TH1D * ZpartonPt_;
  TH1D * ZpartonE_;
  TH1D * ZpartonsDeltaEta_;
  TH1D * ZpartonsDeltaR_;
  TH1D * ZpartonsMass_;
  TH1D * ZpartonsCollinearity_;
  TH1D * TagPartonEta_;
  TH1D * TagPartonPt_;
  TH1D * TagPartonE_;
  TH1D * TagPartonsDeltaEta_;
  TH1D * TagPartonsDeltaR_;
  TH1D * TagPartonsMass_;
  TH1D * TagPartonsCollinearity_;

  TH1D * hadronicZMass_;
  TH1D * hadronicZPt_;
  TH1D * hadronicZCollinearity_;

  TH1D * jetNumber_;
  TH1D * jetUncorrEt_;
  TH1D * jetEt_;
  TH1D * jetPhi_;
  TH1D * jetEta_;

  TH1D * muonJetMatchedNumber_;
  TH1D * muonsMatchedNumber_;

  TH1D * jetEta_aboveEtCut_;
  TH1D * jetNumber_aboveEtCut_;
  TH1D * jetsDeltaR2_aboveEtCut_;

  TH1D * overlapJetsDeltaEt_;
  TH1D * overlapJetsDeltaEta_;
  TH2D * overlapJetsDeltaEtVSDeltaR2_;
  TH2D * overlapJetsDeltaEtaVSDeltaR2_;
  TH2D * overlapJetsDeltaEtVSDeltaEta_;

  TH1D * jetParton_deltaR2_;
  TH1D * jetParton_deltaR2_zoom_;
  TH1D * jetParton_deltaR2max_;

  TH2D * jetParton_deltaEtVSdeltaR2_;
  TH2D * jetParton_deltaEtmeanVSdeltaR2mean_;
  TH2D * jetParton_deltaEVSdeltaR2_;
  TH2D * jetParton_deltaEtaVSdeltaR2_;
  TProfile * jetParton_deltaEtVSdeltaR2_profile_;
  TProfile * jetParton_deltaEtmeanVSdeltaR2mean_profile_;
  TProfile * jetParton_deltaEVSdeltaR2_profile_;
  TProfile * jetParton_deltaEtaVSdeltaR2_profile_;
  TH1D * jetParton_deltaEta_;
  TH1D * jetParton_deltaEt_;

  TH1D * hadronicZrecMass_;
  TH1D * hadronicZrecMassResolution_;
  TH1D * hadronicZrecPt_;
  TH1D * hadronicZrecPtResolution_;
  TH1D * hadronicZrecCollinearity_;
  TH1D * hadronicZrecCollinearityResolution_;
  TH2D * hadronicZrecDeltaPhiVSDeltaEta_;

  TH1D * muonNumber_;
  TH1D * muonEta_;	 
  TH1D * muonPhi_;	 
  TH1D * muonPt_;	 

  TH1D * electronNumber_;
  TH1D * electronEt_;	    
  TH1D * electronEta_;	    
  TH1D * electronPhi_;	    
  TH1D * electronPt_;	    

  TH1D *Bothok_05;
  TH1D *Bothok_02;
  TH1D *Bothok_05_30pc;
  TH1D *Bothok_02_30pc;
  TH1D *Bothok_05_20pc;
  TH1D *Bothok_02_20pc;
  TH1D *Bothok_05_10pc;
  TH1D *Bothok_02_10pc;
  TH1D *Allevents;
  TProfile *DeltaEtvsEt;

  TH1D * Dphiqq;
  TH1D * Detaqq;
  TH1D * Ptminqq;
  TH1D * Ptmaxqq;
  TH1D * Etaminqq;
  TH1D * Etamaxqq;
  TH1D * Dphihh;
  TH1D * Detahh;
  TH1D * Ptminhh;
  TH1D * Ptmaxhh;
  TH1D * Etaminhh;
  TH1D * Etamaxhh;
  TH1D * Ptqq;
  TH1D * Mqq;
  TH1D * Etaqq;
  TH1D * Pthh;
  TH1D * Mhh;
  TH1D * Etahh;
  TH1D * Ptll;
  TH1D * Mll;
  TH1D * Etall;
  TH1D * DphiTjetZjet;
  TH1D * DphiTjetZlep;
  TH1D * DphiminTZ;
  TH1D * DetaTjetZjet;
  TH1D * DetaTjetZlep;
  TH1D * DetaminTZ;
  TH1D * DphiZjetZlep;
  TH1D * DetaZjetZlep;
  TH1D * MassTjetZjet;
  TH1D * MassTjetZlep;
  TH1D * MassZjetZlep;
  TH1D * MassTZZ;
  TH1D * EtaTZZ;
  TH1D * PtTZZ;
 
  TH1D * F_Dphiqq;
  TH1D * F_Detaqq;
  TH1D * F_Ptminqq;
  TH1D * F_Ptmaxqq;
  TH1D * F_Etaminqq;
  TH1D * F_Etamaxqq;
  TH1D * F_Dphihh;
  TH1D * F_Detahh;
  TH1D * F_Ptminhh;
  TH1D * F_Ptmaxhh;
  TH1D * F_Etaminhh;
  TH1D * F_Etamaxhh;
  TH1D * F_Ptqq;
  TH1D * F_Mqq;
  TH1D * F_Etaqq;
  TH1D * F_Pthh;
  TH1D * F_Mhh;
  TH1D * F_Etahh;
  TH1D * F_Ptll;
  TH1D * F_Mll;
  TH1D * F_Etall;
  TH1D * F_DphiTjetZjet;
  TH1D * F_DphiTjetZlep;
  TH1D * F_DphiminTZ;
  TH1D * F_DetaTjetZjet;
  TH1D * F_DetaTjetZlep;
  TH1D * F_DetaminTZ;
  TH1D * F_DphiZjetZlep;
  TH1D * F_DetaZjetZlep;
  TH1D * F_MassTjetZjet;
  TH1D * F_MassTjetZlep;
  TH1D * F_MassZjetZlep;
  TH1D * F_MassTZZ;
  TH1D * F_EtaTZZ;
  TH1D * F_PtTZZ;
 


};

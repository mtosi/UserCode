// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

// Root includes
// -------------
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TFile.h"
#include <memory>

// TMVAtreeMaker
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/TMVAtreeWriter.h"

// GenJets
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"

// Needed to use the classes trained with the TMVA
#include "TMVA/Reader.h"

//
// class declaration
//
class TMVAntple : public edm::EDAnalyzer {

 public:
  explicit TMVAntple(const edm::ParameterSet&);
  ~TMVAntple();

 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // signal-like event TTree
  auto_ptr<TMVAtreeWriter<Float_t> > tmvaSignalTreeWriterPtr_;
  // combinatorial-like event TTree
  auto_ptr<TMVAtreeWriter<Float_t> > tmvaCombinTreeWriterPtr_;
  // signal-like w/ eta cut event TTree
  auto_ptr<TMVAtreeWriter<Float_t> > tmvaEtaCutTreeWriterPtr_;
  // signal-like w/ closest mjj to mZ cut event TTree
  auto_ptr<TMVAtreeWriter<Float_t> > tmvaMZcutTreeWriterPtr_;
  // event variables vectors [value and name]
  vector<TString> eventVariablesNamesVector_;
  vector<Float_t> eventSignalVariablesVector_;
  vector<Float_t> eventCombinVariablesVector_;
  vector<Float_t> eventEtaCutVariablesVector_;
  vector<Float_t> eventMZcutVariablesVector_;

  int eventcounter_;

  int lljjEventCounter_;
  int tmvaMatchingCounter_;
  int combinatorialCounter_;
  int etaCutMatchingCounter_;
  int mZMatchingCounter_;

  int nbin_;

  math::XYZTLorentzVector null_XYZTLorentzVector_;

  // ----------member data ---------------------------
  int           signal_;
  int           whichSim_;
  edm::InputTag electronLabel_;
  edm::InputTag muonLabel_;
  edm::InputTag metLabel_;
  edm::InputTag jetLabel_;
  std::string   corJetsWithBTagLabel_;
  edm::InputTag mcParticleLabel_;
  edm::InputTag genJetLabel_;
  edm::InputTag genMetLabel_;

  double leptonPtCut_;
  double jetEtCut_;
  double jetPartonDeltaR2Cut_;
  double jetLeptonDeltaRCut_;
  double jetEMfracCut_;
  double mHres_;
  double mH_;
  std::string tmvaSignalSuffix_;
  std::string tmvaCombinSuffix_;
  std::string tmvaEtaCutSuffix_;
  std::string tmvaMZcutSuffix_;

  int variablesNumber_;

  // Array of variables for the Reader
  Float_t * variables_;

  TFile* OutputFile;
  // Use a dynamic construction,
  // or the TFile problem will crash the job
  // when moving from one input file to another.
  // The histograms must be created after the TFile is opened.

  std::vector< std::pair<int,int> > vecTypeIndexPair;

  TH1D * eventsNumber_;
  TH1D * goodLeptonicEventsNumber_;
  TH1D * twoJetsAboveEtCutEventsNumber_;
  TH1D * fourJetsAboveEtCutEventsNumber_;
  TH1D * lljjEventNumber_;
  TH1D * tmvaMatchingEventNumber_;
  TH1D * combinatorialMatchEventNumber_;
  TH1D * etaCutMatchingEventNumber_;
  TH1D * mZMatchingEventNumber_;

};

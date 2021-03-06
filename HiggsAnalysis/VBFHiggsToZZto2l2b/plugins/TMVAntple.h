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
  auto_ptr<TMVAtreeWriter<Float_t> > tmvaTreeWriterPtr_;
  // event variables vectors [value and name]
  vector<TString> eventVariablesNamesVector_;
  vector<Float_t> eventVariablesVector_;

  int eventcounter_;
  int tmvaMatchingCounter_;

  int nbin_;

  math::XYZTLorentzVector null_XYZTLorentzVector_;

  // ----------member data ---------------------------
  int signal_;
  int whichSim_;
  edm::InputTag electronLabel_;
  edm::InputTag muonLabel_;
  edm::InputTag metLabel_;
  edm::InputTag mcParticleLabel_;
  std::string leptonicZLabel_;
  std::string hadronicZLabel_;
  std::string tagSystemLabel_;
  std::string tmvaSuffix_;

  int variablesNumber_;

  // Array of variables for the Reader
  Float_t * variables_;

  TFile* OutputFile;
  // Use a dynamic construction,
  // or the TFile problem will crash the job
  // when moving from one input file to another.
  // The histograms must be created after the TFile is opened.

  TH1D * eventsNumber_;
  TH1D * tmvaMatchingEventNumber_;

};

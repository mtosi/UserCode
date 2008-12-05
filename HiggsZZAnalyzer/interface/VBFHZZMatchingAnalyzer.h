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

// Associator for the jets
#include "AnalysisExamples/AnalysisClasses/interface/Associator.h"
#include "AnalysisExamples/AnalysisClasses/interface/AssociatorEt.h"


//
// class declaration
//
class VBFHZZMatchingAnalyzer : public edm::EDAnalyzer {

 public:
  explicit VBFHZZMatchingAnalyzer(const edm::ParameterSet&);
  ~VBFHZZMatchingAnalyzer();

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

  TH2D * ZpartonsDeltaRVSjetMass_;
  TH2D * ZpartonsEtResVSjetMass_;
  TH2D * ZpartonsDeltaRVSetRes_;
  TProfile * ZpartonsDeltaRVSjetMass_profile_;
  TProfile * ZpartonsEtResVSjetMass_profile_;
  TProfile * ZpartonsDeltaRVSetRes_profile_;
 
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
  TH1D * jetParton_deltaEtRes_;
  TH1D * jetParton_004deltaEtRes_;
  TH2D * jetParton_deltaRVSdeltaR_;
  TProfile * jetParton_deltaRVSdeltaR_profile_;

  TH1D * hadronicZrecMass_;
  TH1D * hadronicZrecMassResolution_;
  TH1D * hadronicZrecPt_;
  TH1D * hadronicZrecPtResolution_;
  TH1D * hadronicZrecCollinearity_;
  TH1D * hadronicZrecCollinearityResolution_;

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
  TH1D *Bothok_05_10pc;
  TH1D *Bothok_05_20pc;
  TH1D *Bothok_05_30pc;
  TH1D *Bothok_02_10pc;
  TH1D *Bothok_02_20pc;
  TH1D *Bothok_02_30pc;
  TH1D *Allevents;
};

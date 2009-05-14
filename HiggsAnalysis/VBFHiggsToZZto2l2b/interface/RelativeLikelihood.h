#ifndef RELATIVELIKELIHOOD_H
#define RELATIVELIKELIHOOD_H

/**
 * class RelativeLikelihood RelativeLikelihood.cc HiggsAnalysis/tHMEtplusJetsAnalyzer/src/RelativeLikelihood.cc
 * Package:    RelativeLikelihood
 * Class:      RelativeLikelihood
 *
 * Original Author: Marco De Mattia
 * Creation date: 25/8/2008
 * Mail: demattia@pd.infn.it
 *
 */

// EDM includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
// ------------

// For the level 1 trigger
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/L1Trig.h"

#include "HiggsAnalysis/ttHMEtplusJetsAnalyzer/interface/ttHdecaysCounter.h"
#include "HiggsAnalysis/ttHMEtplusJetsAnalyzer/interface/EventVariables.h"
#include "HiggsAnalysis/ttHMEtplusJetsAnalyzer/interface/QCDbTagMatrix.h"

// Classes to be accessed
// ----------------------
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/BaseJet.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/BaseMEt.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/OfflineMEt.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/OfflineJet.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/SimpleTrack.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/MCParticle.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/GlobalMuon.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/SimpleElectron.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/SimpleTau.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/Summary.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/DeltaR.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/SimpleJet.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2bs/interface/Particle.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/JetVertexAssociator.h"

#include "TFile.h"
#include "TH1D.h"

#include <string>

using namespace std;
using namespace vbfhzz2l2b;
using namespace edm;

// Class declaration
// -----------------
class RelativeLikelihood : public edm::EDAnalyzer {
public:
  explicit RelativeLikelihood(const edm::ParameterSet&);
  ~RelativeLikelihood();
  /// To check if a good muon is found in the event
  bool goodMuon( const GlobalMuonCollection & globalMuons, const OfflineJetCollection & caloJets );
  /// To check if a good electron is found in the event
  bool goodElectron( const SimpleElectronCollection & simpleElectrons, const OfflineJetCollection & caloJets );

protected:
  /// fills the histogram with the likelihood value for this event
  void evaluateLikelihood( const vector<double> & eventVariablesVector );
  virtual void beginJob(const edm::EventSetup&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  edm::ParameterSet conf_;

  edm::InputTag cenJetLabel_;
  edm::InputTag forJetLabel_;
  edm::InputTag tauJetLabel_;
  edm::InputTag l1MEtLabel_;
  edm::InputTag offlineJetLabel_;
  edm::InputTag offlineMEtLabel_;
  edm::InputTag MCParticleLabel_;
  edm::InputTag globalMuonLabel_;
  edm::InputTag simpleElectronLabel_;
  edm::InputTag simpleTauLabel_;
  edm::InputTag summaryLabel_;
  edm::InputTag vtxLabel_;
  bool withL1ForwardJets_;
  bool vtxAssoc_;
  string higgsFileName_;
  string hadronicTopFileName_;
  string qcdFileName_;
  string qcdHistoFileName_;
  double jetEtCut_;
  double jetEtaCut_;
  string countTTHdecaysFileName_;
  string countTTHdecays2tagsFileName_;
  string inputFileNameSignal_;
  string inputFileNameBackground_;
  string outputFileName_;
  string tmvaSuffix_;

  bool useTagMatrixForQCD_;

  int eventCounter_;

  // Declare as static so that only one exists, even if more
  // than one TDAna object is created
  // -------------------------------------------------------
  static L1Trig L1Trigger;

  int l1Eff_;

  TFile * inputFileSignal_;
  TFile * inputFileBackground_;
  TFile * outputFile_;
  TDirectory * outputDir_;

  // Histograms: the number depends on those written by ttHMEtplusJetsAnalyzer
  vector<TH1D> histogramVariableSignal_;
  vector<TH1D> histogramVariableBackground_;

  TH1D * relativeLikelihood_;

  ttHdecaysCounter * countTTHdecays_;
  ttHdecaysCounter * countTTHdecays2tags_;

  // Class to fill histograms on event variables
  EventVariables * eventVariables2Tags_;
  QCDbTagMatrix * qcdbTagMatrixMultiplier_;
 
  JetVertexAssociator * jetVertexAssociator_;

};


#endif // RELATIVELIKELIHOOD_CC

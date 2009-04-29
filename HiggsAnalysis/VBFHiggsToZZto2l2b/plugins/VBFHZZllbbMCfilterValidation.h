// -*- C++ -*-
//
// Package:    VBFHZZllbbMCfilterValidation
// Class:      VBFHZZllbbMCfilterValidation
// 
/**\class VBFHZZllbbMCfilterValidation VBFHZZllbbMCfilterValidation.cc HiggsAnalysis/VBFHiggsToZZto2l2b/src/VBFHZZllbbMCfilterValidation.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Mia TOSI
//         Created:  Tue Jan 20 15:48:58 CET 2009
// $Id: VBFHZZllbbMCfilterValidation.h,v 1.2 2009/04/29 11:19:34 tosi Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Root includes
// -------------
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TFile.h"

// For file output
// ---------------
#include <fstream>
#include <sstream>
#include <cmath>
#include <memory>
#include <string>
#include <iostream>
#include <iomanip>


// class decleration
//

class VBFHZZllbbMCfilterValidation : public edm::EDAnalyzer {
   public:
      explicit VBFHZZllbbMCfilterValidation(const edm::ParameterSet&);
      ~VBFHZZllbbMCfilterValidation();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

  int getEventID(const edm::Event&, const int &);

//
// constants, enums and typedefs
//

//
// static data member definitions
//

  // ----------member data ---------------------------
  unsigned int eventcounter_;
  unsigned int eventVBFcounter_;
  unsigned int eventGGFcounter_;

  unsigned int jetNumberCut_eventcounter_;
  unsigned int firstJetPtCut_eventcounter_;
  unsigned int secondJetPtCut_eventcounter_;
  unsigned int invMassCut_eventcounter_;
  unsigned int deltaEtaCut_eventcounter_;
  unsigned int leptonPtCut_eventcounter_;
  unsigned int jetNumberCut_eventVBFcounter_;
  unsigned int firstJetPtCut_eventVBFcounter_;
  unsigned int secondJetPtCut_eventVBFcounter_;
  unsigned int invMassCut_eventVBFcounter_;
  unsigned int deltaEtaCut_eventVBFcounter_;
  unsigned int leptonPtCut_eventVBFcounter_;
  unsigned int jetNumberCut_eventGGFcounter_;
  unsigned int firstJetPtCut_eventGGFcounter_;
  unsigned int secondJetPtCut_eventGGFcounter_;
  unsigned int invMassCut_eventGGFcounter_;
  unsigned int deltaEtaCut_eventGGFcounter_;
  unsigned int leptonPtCut_eventGGFcounter_;


  int           whichSim_;
  bool          signal_;
  edm::InputTag genJetLabel_;
  edm::InputTag genParticleLabel_;

  int    jetNumberCut_;
  double firstJetPtCut_;
  double secondJetPtCut_;
  double invMassCut_;
  double deltaEtaCut_;
  double leptonPtCut_;

  TH1D * eventsNumber_;
  TH1D * eventsIDNumber_;

  TH1D * VBFfirstJetPt_; 
  TH1D * VBFsecondJetPt_; 
  TH1D * VBFthirdJetPt_;
  TH1D * VBFfourthJetPt_;
  TH1D * VBFmaxDeltaEtaJetJet_;
  TH1D * VBFmaxInvMassJetJet_;
  TH1D * ggFfirstJetPt_; 
  TH1D * ggFsecondJetPt_; 
  TH1D * ggFthirdJetPt_;
  TH1D * ggFfourthJetPt_;
  TH1D * ggFmaxDeltaEtaJetJet_;
  TH1D * ggFmaxInvMassJetJet_;

  TH1D * firstElectronPt_;
  TH1D * secondElectronPt_;
  TH1D * firstMuonPt_;
  TH1D * secondMuonPt_;

};


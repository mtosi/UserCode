#ifndef StarterKit_test_PatAnalyzerKit_h
#define StarterKit_test_PatAnalyzerKit_h



// -*- C++ -*-
//// -*- C++ -*-
//
// Package:    PatAnalyzerKit
// Class:      PatAnalyzerKit
//
/**


*/
//-------------------------------------------------------------------------------------
//!\class PatAnalyzerKit PatAnalyzerKit.cc PhysicsTools/StarterKit/plugins/PatAnalyzerKit.cc
//!\brief PatAnalyzerKit is a generic ED analyzer with predefined physics histograms.
//!
//!  This is an ED analyzer which creates and fills a bunch of
//!  histograms of various physics quantities.  However, in order to
//!  make this also work in FWLite, most of the actual work is
//!  performed by another object called PhysicsHistograms, to which this
//!  ED analyzer delegates most of the work, except the interactions with
//!  EDM and Framework, like:
//!
//!  - obtaining parameters from edm::ParameterSet.  PhysicsHistograms receives
//!    commands with lists of sub-components to manipulate (usually only enable or
//!    disable).
//!
//!  - fetching collections from the event (the iteration over collections is done
//!    by PhysicsHistograms).
//!
//!  - putting single numbers back to the event -- this is how flat ntuples
//!    are made. (PhysicsHistograms provides several lists of PhysVarHisto* pointers:
//!       1. all PhysVarHistos
//!       2. those that need to be histogrammed
//!       3. those that need to be ntupled
//!       4. those that need to be filled in each events (i.e. "active" ones, which
//!          is a union of (2) and (3), and a subset of (1).
//!
//-------------------------------------------------------------------------------------
//
// Original Author:  Eric Vaandering, Salvatore Rappoccio
//         Created:  Wed Nov 28 15:31:57 CST 2007
// $Id: PatAnalyzerKit.h,v 1.3 2008/10/24 21:18:27 srappocc Exp $
//
// Revision History:
//       -  Sal Rappoccio, Fri Nov 30 12:49:44 CST 2007: Added other objects as first
//          stab of full analysis toolkit.
//       -  Sal Rappoccio, Mon Mar 03 12:00:00 CST 2008: Added photons and taus
//       -  Sal Rappoccio, Mon Apr 14 13:45:59 CDT 2008: Added CSA07 information.
//       -  Sal Rappoccio, Thu May 29 12:05:58 CDT 2008: Restructured SK, renamed this to PatAnalyzerKit
//-------------------------------------------------------------------------------------

// system include files
#include <memory>
#include <fstream>

// user include files
#include "PhysicsTools/StarterKit/interface/PatKitHelper.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "PhysicsTools/UtilAlgos/interface/TFileDirectory.h"

//
// class declaration
//

class PatAnalyzerKit : public edm::EDProducer  // public edm::EDAnalyzer
{
public:
  explicit PatAnalyzerKit(const edm::ParameterSet&);
  virtual ~PatAnalyzerKit();

protected:

  // beginJob
  virtual void beginJob(const edm::EventSetup&) ;
  // produce is where the ntuples are made
  virtual void produce( edm::Event &, const edm::EventSetup & );
  // endJob
  virtual void endJob() ;

  // The main sub-object which does the real work
  pat::PatKitHelper    helper_;

  // Verbosity
  int             verboseLevel_;

  // Physics objects handles
  edm::Handle<std::vector<pat::Muon> >     muonHandle_;
  edm::Handle<std::vector<pat::Electron> > electronHandle_;
  edm::Handle<std::vector<pat::Tau> >      tauHandle_;
  edm::Handle<std::vector<pat::Jet> >      jetHandle_;
  edm::Handle<std::vector<pat::MET> >      METHandle_;
  edm::Handle<std::vector<pat::Photon> >   photonHandle_;
  edm::Handle<std::vector<reco::RecoChargedCandidate> >   trackHandle_;
  edm::Handle<std::vector<reco::GenParticle> > genParticlesHandle_;

};



#endif

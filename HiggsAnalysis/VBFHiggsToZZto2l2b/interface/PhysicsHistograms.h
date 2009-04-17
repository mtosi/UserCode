#ifndef HiggsAnalysis_VBFHiggsToZZto2l2b_PhysicsHistograms_h
#define HiggsAnalysis_VBFHiggsToZZto2l2b_PhysicsHistograms_h

//------------------------------------------------------------------------
//!  \class PhysicsHistograms
//!  \brief Object to manage and fill various physics histograms
//!
//!  The order how the operations must be executed.
//!
//!  1. we first build our own default histogram groups (electrons, muons, etc)
//!
//!  2. the user-defined histogram groups are added by add*HistoGroup() methods.
//!
//!  3. configure starts:
//!     all PhysVarHisto pointers are collected in one big flat array for
//!     easy access and speedy processing.
//!
//!  4. various histograms are disabled.
//!
//!  5. various histograms are enabled.  configure ends.
//!
//!  At this point we're good to go and ready to see the events.
//------------------------------------------------------------------------



// system include files
#include <memory>
#include <fstream>


// user include files
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/HistoMuon.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/HistoElectron.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/HistoJet.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/HistoMET.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/HistoGenParticle.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "PhysicsTools/UtilAlgos/interface/TFileDirectory.h"

//
//--- Class declaration.
//


// Function to print out candidates
std::ostream & operator<<( std::ostream & out, const reco::Candidate & cand );

class PhysicsHistograms  {
public:

  struct KinAxisLimits {
  
    double pt1, pt2, m1, m2;

    KinAxisLimits( double apt1=0, double apt2=0, double am1=0, double am2=0 ) :
      pt1(apt1), pt2(apt2), m1(am1), m2(am2)
    {
    }
    ~KinAxisLimits() {}
  };


  explicit PhysicsHistograms ( KinAxisLimits const & muonAxis, 
			       KinAxisLimits const & electronAxis, 
			       KinAxisLimits const & jetAxis, 
			       KinAxisLimits const & METAxis, 
			       KinAxisLimits const & genParticlesAxis
			       );
  virtual ~PhysicsHistograms();

  //--- Standard methods used in the event processing, called either by ED analyzer
  //    (from the methods of the same name), or by the FWLite macro which does the
  //    event loop).
  //
  virtual void beginJob();  //!<  initialize before seeing any events
  virtual void   endJob();  //!<  do whatever is needed after seeing all events


  //--- Configuration.
  virtual void configure( std::string & histos_to_disable,   // comma separated list of names
			  std::string & histos_to_enable );  // comma separated list of names


  //--- Selection of a subset of PhysVarHistos.
  virtual void select( std::string  vars_to_select,   // comma separated list of names
		       std::vector< vbfhzz2l2b::PhysVarHisto * > & selectedHistos );

  //--- Clear cache vector for PhysVarHisto
  virtual void clearVec();

  //--- Specific actions for the event.
  // &&& Design note: we could have used overloaded fill() everywhere, but
  // &&&              the novices may find it confusing.

  //--- Process a whole collection of Muons...
  //
  inline void fillCollection( const reco::MuonCollection & coll, double w = 1.0 )
    { muonHistograms_->fillCollection(coll,w); }

  //--- ... or Electrons...
  //
  inline void fillCollection( const reco::GsfElectronCollection & coll, double w = 1.0 )
    { electronHistograms_->fillCollection(coll,w); }

  //--- ... or Jets...
  //
  inline void fillCollection( const reco::CaloJetCollection & coll, double w = 1.0 )
    { jetHistograms_->fillCollection(coll,w); }

  //--- ... or MET.
  //
  inline void fillCollection( const reco::CaloMETCollection & coll, double w = 1.0 )
    { metHistograms_->fillCollection(coll,w); }

  //--- ... or GenParticle.
  //
  inline void fillCollection( const std::vector<reco::GenParticle> & coll, double w = 1.0 )
    { genParticleHistograms_->fillCollection(coll,w); }



  // &&& Design note: again, let's be explicit.  This could be compressed into
  // &&&              fewer functions, but at the expense of more complicated
  // &&&              code under the hood, and also an interface which is a teeny
  // &&&              harder to master (and we are trying to avoid that; the
  // &&&              interface should be as dumb as possible).

  //--- Add one histo to muon group, or a whole group of muon histograms
  //
  inline void addMuonHisto ( vbfhzz2l2b::PhysVarHisto * h )
    { muonHistograms_->addHisto(h); }
  inline void addMuonHistoGroup( vbfhzz2l2b::HistoMuon * hgr )
    { muonHistograms_->addHistoGroup(hgr); }

  //--- Add one histo to electron group, or a whole group of electron histograms
  //
  inline void addElectronHisto ( vbfhzz2l2b::PhysVarHisto * h )
    { electronHistograms_->addHisto(h); }
  inline void addElectronHistoGroup( vbfhzz2l2b::HistoElectron * hgr )
    { electronHistograms_->addHistoGroup(hgr); }

  //--- Add one histo to jet group, or a whole group of jet histograms
  //
  inline void addJetHisto ( vbfhzz2l2b::PhysVarHisto * h )
    { jetHistograms_->addHisto(h); }
  inline void addJetHistoGroup( vbfhzz2l2b::HistoJet * hgr )
    { jetHistograms_->addHistoGroup(hgr); }

  //--- Add one histo to MET group, or a whole group of MET histograms
  //
  inline void addMetHisto ( vbfhzz2l2b::PhysVarHisto * h )
    { metHistograms_->addHisto(h); }
  inline void addMetHistoGroup( vbfhzz2l2b::HistoMET * hgr )
    { metHistograms_->addHistoGroup(hgr); }

  //--- Add one histo to genParticle group, or a whole group of genParticle histograms
  //
  inline void addGenParticleHisto ( vbfhzz2l2b::PhysVarHisto * h )
    { genParticleHistograms_->addHisto(h); }
  inline void addGenParticleHistoGroup( vbfhzz2l2b::HistoGenParticle * hgr )
    { genParticleHistograms_->addHistoGroup(hgr); }

  //--- Add one generic histo to list
  inline void addHisto( vbfhzz2l2b::PhysVarHisto * h )
    { allVarHistos_.push_back( h ); }



private:

  // Parameters for running
  std::string     outputTextName_;

  // Histogram server
  edm::Service<TFileService> fs;

  // Histogram objects that make "standard" plots for each object
  vbfhzz2l2b::HistoMuon        * muonHistograms_;
  vbfhzz2l2b::HistoElectron    * electronHistograms_;
  vbfhzz2l2b::HistoMET         * metHistograms_;
  vbfhzz2l2b::HistoJet         * jetHistograms_;
  vbfhzz2l2b::HistoGenParticle * genParticleHistograms_;

  //--- The summary of all PhysVarHistos.
  // &&& Is this still needed?
  std::vector< vbfhzz2l2b::PhysVarHisto* > allVarHistos_ ;
  std::vector< vbfhzz2l2b::PhysVarHisto* > enabledVarHistos_ ;

  //--- This is a nice feature but let's not worry about it for now. &&&
  ofstream        outputFile_;
};

#endif

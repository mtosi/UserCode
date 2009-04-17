//------------------------------------------------------------------------
// -*- C++ -*-
//
//! \class PhysicsHistograms PhysicsHistograms.cc Demo/TempAnaToolkit/src/PhysicsHistograms.cc
//!
//!  Description: Demonstration of a simple analysis toolkit for starter analyses
//!
//
// Original Author:  Petar Maksimovic
//         Created:  Christmas 2007
// $Id: PhysicsHistograms.cc,v 1.6 2008/10/24 21:18:29 srappocc Exp $
//
// Revision History:
//------------------------------------------------------------------------


#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/PhysicsHistograms.h"

#include <string>
#include <sstream>

//------------------------------------------------------------------------
//!  Create the objects that manage the histogram groups.
//------------------------------------------------------------------------
PhysicsHistograms::PhysicsHistograms( KinAxisLimits const & muonAxis, 
				      KinAxisLimits const & electronAxis, 
				      KinAxisLimits const & jetAxis, 
				      KinAxisLimits const & METAxis, 
				      KinAxisLimits const & genParticlesAxis
				      )
{
  //--- Initialize histogram objects
  std::cout << "PhysicsHistograms: Creating muon histograms" << std::endl;
  muonHistograms_     = 
    new vbfhzz2l2b::HistoMuon    ( "muon", "Muon", "mu", muonAxis.pt1, muonAxis.pt2, muonAxis.m1, muonAxis.m2 );
  std::cout << "PhysicsHistograms: Creating electron histograms" << std::endl;
  electronHistograms_ = 
    new vbfhzz2l2b::HistoElectron( "electron", "Electron", "e", electronAxis.pt1, electronAxis.pt2, electronAxis.m1, electronAxis.m2 );
  std::cout << "PhysicsHistograms: Creating jet histograms" << std::endl;
  jetHistograms_      =
    new vbfhzz2l2b::HistoJet     ( "jet", "Jet", "jet", jetAxis.pt1, jetAxis.pt2, jetAxis.m1, jetAxis.m2 );
  std::cout << "PhysicsHistograms: Creating met histograms" << std::endl;
  metHistograms_      = 
    new vbfhzz2l2b::HistoMET     ( "met", "MET", "met", METAxis.pt1, METAxis.pt2, METAxis.m1, METAxis.m2 );
  std::cout << "PhysicsHistograms: Creating genParticle histograms" << std::endl;
  genParticleHistograms_    = 
    new vbfhzz2l2b::HistoGenParticle   ( "genParticle", "GenParticle", "genParticle", genParticlesAxis.pt1, genParticlesAxis.pt2, genParticlesAxis.m1, genParticlesAxis.m2 );
}


//------------------------------------------------------------------------
//!  Destroy the objects that manage the histogram groups.
//!
//!  Note that the TH1's used by PhysVarHistos managed by these histo
//!  groups will *not* be deleted in the PhysVarHisto's destructor. So
//!  it's safe to delete both HistoGroups and PhysVarHistos.
//------------------------------------------------------------------------
PhysicsHistograms::~PhysicsHistograms()
{
  delete muonHistograms_     ;
  delete electronHistograms_ ;
  delete jetHistograms_      ;
  delete metHistograms_      ;
  delete genParticleHistograms_    ;

  for ( unsigned int i = 0; i < allVarHistos_.size(); ++i ) {
    if ( allVarHistos_[i] ) delete allVarHistos_[i];
  }

  outputFile_.close();
}



//------------------------------------------------------------------------
//!  Methods to configure (enable or disable) various PhysVarHistos one at the time,
//!  or in whole groups.
//------------------------------------------------------------------------
void
PhysicsHistograms::configure( std::string & histos_to_disable,   // comma separated list of names
			      std::string & histos_to_enable )   // comma separated list of names
{
  std::cout << "PhysicsHistograms:: configuring..."
	    << "\n   First disabling: " << histos_to_disable
	    << "\n   Then  enabling : " << histos_to_enable
	    << std::endl;


  //--- Pass this information to histogramGroups
  muonHistograms_        ->configure( histos_to_disable, histos_to_enable ) ;
  electronHistograms_    ->configure( histos_to_disable, histos_to_enable ) ;
  metHistograms_         ->configure( histos_to_disable, histos_to_enable ) ;
  jetHistograms_         ->configure( histos_to_disable, histos_to_enable ) ;
  genParticleHistograms_ ->configure( histos_to_disable, histos_to_enable ) ;

}


//------------------------------------------------------------------------
//!  Selection of a subset of PhysVarHistos.
//------------------------------------------------------------------------
void
PhysicsHistograms::select( std::string  vars_to_select,   // comma separated list of names
			   std::vector< vbfhzz2l2b::PhysVarHisto * > & selectedVars )
{
  std::cout << "PhysicsHistograms:: selecting the following variables:\n\t"
	    << vars_to_select
	    << std::endl;


  //--- Pass this information to histogramGroups
  muonHistograms_        ->select( vars_to_select, selectedVars ) ;
  electronHistograms_    ->select( vars_to_select, selectedVars ) ;
  metHistograms_         ->select( vars_to_select, selectedVars ) ;
  jetHistograms_         ->select( vars_to_select, selectedVars ) ;
  genParticleHistograms_ ->select( vars_to_select, selectedVars ) ;
  
  std::vector<vbfhzz2l2b::PhysVarHisto*>::iterator i = allVarHistos_.begin();
  std::vector<vbfhzz2l2b::PhysVarHisto*>::iterator end = allVarHistos_.end();
  std::string temp = "," + vars_to_select + ",";
  for ( ; i != end; ++i  ) {
    std::string test = "," + (*i)->name() + ",";
    std::cout << "testing " << test << std::endl;
    if ( temp.find( test ) != std::string::npos || temp == ",all," ) {
      std::cout << "FOUND!" << std::endl;
      selectedVars.push_back ( *i );
    }
  }

  std::cout << "PhysicsHistograms:: selected " << selectedVars.size()
	    << " variables." << std::endl;
}

void PhysicsHistograms::clearVec()
{
  muonHistograms_        ->clearVec() ;
  electronHistograms_    ->clearVec() ;
  metHistograms_         ->clearVec() ;
  jetHistograms_         ->clearVec() ;
  genParticleHistograms_ ->clearVec() ;
  for ( uint i = 0; i < allVarHistos_.size(); i++ ) {
    allVarHistos_[i]->clearVec();
  }
}

//------------------------------------------------------------------------
//!  Method called before seeing all events (in CMSSW) or before the
//!  event loop in FWLite.
//------------------------------------------------------------------------
void
PhysicsHistograms::beginJob()
{
  // Dummy for now
}


//------------------------------------------------------------------------
//!  Method called after seeing all events (in CMSSW) or after the
//!  event loop in FWLite.
//------------------------------------------------------------------------
void
PhysicsHistograms::endJob()
{
  // Dummy for now
}



//------------------------------------------------------------------------
//!  Method to print out reco::Candidates.
// &&& Design suggestion: this should go into HistoGroup<> instead.
//------------------------------------------------------------------------
std::ostream & operator<<( std::ostream & out, const reco::Candidate & cand )
{
  char buff[1000];
  sprintf( buff, "Pt, Eta, Phi, M = (%6.2f, %6.2f, %6.2f, %6.2f)",
           cand.pt(), cand.eta(), cand.phi(), cand.mass() );
  out << buff;
  return out;
}

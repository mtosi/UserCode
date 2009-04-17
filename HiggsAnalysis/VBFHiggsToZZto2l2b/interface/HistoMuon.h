#ifndef HiggsAnalysis_VBFHiggsToZZto2l2b_HistoMuon_h
#define HiggsAnalysis_VBFHiggsToZZto2l2b_HistoMuon_h

//------------------------------------------------------------
// Title: HistoMuon.h
// Purpose: To histogram Muons
//------------------------------------------------------------
//
// Interface:
//
//   HistoMuon ( TFile * file );
//   Description: Constructor.
//
//   void fill( TK::Muon * );
//   Description: Fill object. Will fill relevant muon variables
//
//   void write();
//   Description: Write object to file in question.
//
//   ~HistoMuon
//    Description: Destructor. Deallocates memory.
//
//------------------------------------------------------------
// This package's include files
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/HistoGroup.h"

// CMSSW include files
#include "DataFormats/MuonReco/interface/Muon.h"

// STL include files
#include <string>
#include <vector>

// ROOT include files
#include <TH1D.h>

namespace vbfhzz2l2b {

  using namespace reco;

  class HistoMuon : public HistoGroup<Muon> {

  public:
    HistoMuon(std::string dir = "muon", std::string group = "Muon", std::string pre ="mu",
	      double pt1=0., double pt2=200., double m1=0., double m2=200.,
	      TFileDirectory * parentDir=0);
    virtual ~HistoMuon() { } ;

    // fill a plain ol' muon:
    virtual void fill( const Muon *muon, uint iPart = 1, double weight = 1.0 );
    virtual void fill( const Muon &muon, uint iPart = 1, double weight = 1.0 ) { 
      fill(&muon, iPart,weight); 
    }

    // fill a muon that is a shallow clone, and take kinematics from 
    // shallow clone but detector plots from the muon itself
    virtual void fill( const reco::ShallowClonePtrCandidate *muon, uint iPart = 1, double weight = 1.0 );
    virtual void fill( const reco::ShallowClonePtrCandidate &muon, uint iPart = 1, double weight = 1.0 ) { 
      fill(&muon, iPart,weight); 
    }

    virtual void fillCollection( const std::vector<Muon> & coll, double weight = 1.0 );

    // Clear ntuple cache
    void clearVec();

  protected:
    PhysVarHisto * h_caloE_    ;

    //muon energy deposit analyzer
    PhysVarHisto * ecalDepEnergy_;
    PhysVarHisto * ecalS9DepEnergy_ ;
    PhysVarHisto * hcalDepEnergy_ ;
    PhysVarHisto * hcalS9DepEnergy_ ;
    PhysVarHisto * hoDepEnergy_ ;
    PhysVarHisto * hoS9DepEnergy_ ;

  };
}
#endif

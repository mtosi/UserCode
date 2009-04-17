#ifndef HiggsAnalysis_VBFHiggsToZZto2l2b_HistoComposite_h
#define HiggsAnalysis_VBFHiggsToZZto2l2b_HistoComposite_h

//------------------------------------------------------------
// Title: HistoComposite.h
// Purpose: To histogram Composites
//------------------------------------------------------------
//
// Interface:
//
//   HistoComposite ( TFile * file );
//   Description: Constructor.
//
//   void fill( TK::Composite * );
//   Description: Fill object. Will fill relevant muon variables
//
//   void write();
//   Description: Write object to file in question.
//
//   ~HistoComposite
//    Description: Destructor. Deallocates memory.
//
//------------------------------------------------------------
//
// Modification History:
//
//------------------------------------------------------------


// CMSSW include files
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/HistoMuon.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/HistoElectron.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/HistoJet.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/HistoMET.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/ShallowClonePtrCandidate.h"
#include "DataFormats/Common/interface/Ptr.h"

// STL include files
#include <string>

// ROOT include files
#include <TH1D.h>
#include <TFile.h>

namespace vbfhzz2l2b {

  using namespace reco;

  template <class T>
  class HistoMap {
  public:
    typedef std::string                              key_type;
    typedef T *                                      data_type;
    typedef std::map<key_type, data_type >           map_type;
    
    map_type map;
  };

  class HistoComposite : public HistoGroup<reco::CompositeCandidate> {

   public:


    typedef edm::Ptr<Muon>        MuonPtr;
    typedef edm::Ptr<GsfElectron> ElectronPtr;
    typedef edm::Ptr<CaloJet>     JetPtr;
    typedef edm::Ptr<CaloMET>     METPtr;
    
    
    HistoComposite(std::string dir, std::string candTitle, std::string candName,
		   double pt1=0, double pt2=200, double m1=0, double m2=200,
		   TFileDirectory * parentDir = 0 );
    virtual ~HistoComposite();

    // void fill( reco::CompositeCandidate * cand );
    void fill( const reco::CompositeCandidate * cand, double weight = 1.0 );
    void fill( const reco::CompositeCandidate & cand, double weight = 1.0 ) { return fill(&cand, weight); }


    void fill( const reco::ShallowClonePtrCandidate * pshallow, double weight = 1.0 );
    void fill( const reco::ShallowClonePtrCandidate & pshallow, double weight = 1.0 )
    { fill(&pshallow, weight); }

   protected:
    std::string       candName_;


    HistoMap<HistoMuon>        histoMuon_;
    HistoMap<HistoElectron>    histoElectron_;
    HistoMap<HistoJet>         histoJet_;
    HistoMap<HistoMET>         histoMET_;
    HistoMap<HistoComposite>   histoComposite_;
  };


}
#endif

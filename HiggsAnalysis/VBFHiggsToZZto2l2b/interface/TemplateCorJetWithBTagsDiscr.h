#ifndef TEMPLATECORJETWITHBTAGSDISCR_H
#define TEMPLATECORJETWITHBTAGSDISCR_H

#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"

//#include "DataFormats/Common/interface/Ref.h"
/**
 *
 * Used for offline jets. Includes b-tagging informations.
 *
 */

namespace vbfhzz2l2b {

  template<class T>
  class TemplateCorJetWithBTagsDiscr {
  public:
    TemplateCorJetWithBTagsDiscr( const T & JETREF,
			  const double & UNCORRET, 
			  const float & DISCRIMINATORHIGHEFF, 
			  const float & DISCRIMINATORHIGHPUR, 
			  const float & DISCRIMINATORCOMBSECVTX,   
			  const float & DISCRIMINATORCOMBSECVTXMVA,
			  const float & DISCRIMINATORSOFTMUON,     
			  const float & DISCRIMINATORSOFTELECTRON, 
			  const float & DISCRIMINATORJETPROB      
			  ) :
      jetRef_ ( JETREF ),
      uncorrEt_( UNCORRET ),
      discriminatorHighEff_      ( DISCRIMINATORHIGHEFF       ),
      discriminatorHighPur_      ( DISCRIMINATORHIGHPUR       ),
      discriminatorCombSecVtx_   ( DISCRIMINATORCOMBSECVTX    ),   
      discriminatorCombSecVtxMVA_( DISCRIMINATORCOMBSECVTXMVA ),
      discriminatorSoftMuon_     ( DISCRIMINATORSOFTMUON      ),     
      discriminatorSoftElectron_ ( DISCRIMINATORSOFTELECTRON  ), 
      discriminatorJetProb_      ( DISCRIMINATORJETPROB       )      
    {
    }
    // Default constructor, only needed for classes.h
    TemplateCorJetWithBTagsDiscr() {}

    T jetRef()                const { return jetRef_; }
    T* jet()                        { return &jetRef_; }
    double uncorrEt()         const { return uncorrEt_; }
    double emEnergyFraction() const { return emEnergyFraction_; }
    float discriminatorHighEff()       const { return discriminatorHighEff_;       }
    float discriminatorHighPur()       const { return discriminatorHighPur_;       }
    float discriminatorCombSecVtx()    const { return discriminatorCombSecVtx_;    }
    float discriminatorCombSecVtxMVA() const { return discriminatorCombSecVtxMVA_; }
    float discriminatorSoftMuon()      const { return discriminatorSoftMuon_;      }
    float discriminatorSoftElectron()  const { return discriminatorSoftElectron_;  }
    float discriminatorJetProb()       const { return discriminatorJetProb_;       }

    void setJetRef          ( const T      & JETREF           ) { jetRef_           = JETREF;           }
    void setUncorrEt        ( const double & UNCORRET         ) { uncorrEt_         = UNCORRET;         }
    void setEmEnergyFraction( const double & EMENERGYFRACTION ) { emEnergyFraction_ = EMENERGYFRACTION; }
    void setDiscriminatorHighEff      ( const float & DISCRIMINATORHIGHEFF       ) { discriminatorHighEff_       = DISCRIMINATORHIGHEFF;       }
    void setDiscriminatorHighPur      ( const float & DISCRIMINATORHIGHPUR       ) { discriminatorHighPur_       = DISCRIMINATORHIGHPUR;       }
    void setDiscriminatorCombSecVtx   ( const float & DISCRIMINATORCOMBSECVTX    ) { discriminatorCombSecVtx_    = DISCRIMINATORCOMBSECVTX;    }
    void setDiscriminatorCombSecVtxMVA( const float & DISCRIMINATORCOMBSECVTXMVA ) { discriminatorCombSecVtxMVA_ = DISCRIMINATORCOMBSECVTXMVA; }
    void setDiscriminatorSoftMuon     ( const float & DISCRIMINATORSOFTMUON	 ) { discriminatorSoftMuon_      = DISCRIMINATORSOFTMUON;      }
    void setDiscriminatorSoftElectron ( const float & DISCRIMINATORSOFTELECTRON  ) { discriminatorSoftElectron_  = DISCRIMINATORSOFTELECTRON;  }
    void setDiscriminatorJetProb      ( const float & DISCRIMINATORJETPROB       ) { discriminatorJetProb_       = DISCRIMINATORJETPROB;       }

  protected:
    T jetRef_;
    double uncorrEt_;
    double emEnergyFraction_;
    float discriminatorHighEff_;
    float discriminatorHighPur_;
    float discriminatorCombSecVtx_;
    float discriminatorCombSecVtxMVA_;
    float discriminatorSoftMuon_;
    float discriminatorSoftElectron_;
    float discriminatorJetProb_;

  };

  typedef std::vector<vbfhzz2l2b::TemplateCorJetWithBTagsDiscr<reco::CaloJet> > TemplateCorCaloJetWithBTagsDiscrCollection;
  typedef std::vector<vbfhzz2l2b::TemplateCorJetWithBTagsDiscr<reco::GenJet> >  TemplateCorGenJetWithBTagsDiscrCollection;
  typedef std::vector<vbfhzz2l2b::TemplateCorJetWithBTagsDiscr<reco::PFJet> >   TemplateCorPFJetWithBTagsDiscrCollection;

}

#endif // TEMPLATECORJETWITHBTAGSDISCR_H

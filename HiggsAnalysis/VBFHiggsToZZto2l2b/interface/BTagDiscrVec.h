#ifndef BTAGDISCRVEC_H
#define BTAGDISCRVEC_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/Common/interface/RefVector.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/VBFHZZllbbUtils.h"
#include <string>

namespace vbfhzz2l2b {

  class BTagDiscrVec {
    
  public:
    BTagDiscrVec( const std::vector<double> & DISCRIMINATORSVEC ) :
      discriminatorsVec_( DISCRIMINATORSVEC )
    {
    }
    
    // Default constructor, only needed for classes.h
    BTagDiscrVec() {}
    ~BTagDiscrVec() {}

    std::vector<double> discriminators() const { return discriminatorsVec_;                }
    double discriminatorHighEff()        const { 
      std::cout << "[BTagDiscrVec::discriminatorHighEff] discriminator: " << discriminatorsVec_[vbfhzz2l2b::HIGHEFF] << std::endl;
      return discriminatorsVec_[vbfhzz2l2b::HIGHEFF];       
    }
    double discriminatorHighPur()        const { return discriminatorsVec_[vbfhzz2l2b::HIGHPUR];       }
    double discriminatorCombSecVtx()     const { return discriminatorsVec_[vbfhzz2l2b::COMBSECVTX];    }
    double discriminatorCombSecVtxMVA()  const { return discriminatorsVec_[vbfhzz2l2b::COMBSECVTXMVA]; }
    double discriminatorSoftMuon()       const { return discriminatorsVec_[vbfhzz2l2b::SOFTMUON];      }
    double discriminatorSoftElectron()   const { return discriminatorsVec_[vbfhzz2l2b::SOFTELECTRON];  }
    double discriminatorJetProb()        const { return discriminatorsVec_[vbfhzz2l2b::JETPROB];       }
    double discriminator( const std::string& BTAGGER ) const { return discriminatorsVec_[vbfhzz2l2b::bTaggerCode(BTAGGER)]; }
   

    void setDiscriminators            ( const std::vector<double> DISCRIMINATORSVEC ) { discriminatorsVec_ = DISCRIMINATORSVEC;                         }
    void setDiscriminatorHighEff      ( const double & DISCRIMINATORHIGHEFF         ) { discriminatorsVec_[HIGHEFF]       = DISCRIMINATORHIGHEFF;       }
    void setDiscriminatorHighPur      ( const double & DISCRIMINATORHIGHPUR         ) { discriminatorsVec_[HIGHPUR]       = DISCRIMINATORHIGHPUR;       }
    void setDiscriminatorCombSecVtx   ( const double & DISCRIMINATORCOMBSECVTX      ) { discriminatorsVec_[COMBSECVTX]    = DISCRIMINATORCOMBSECVTX;    }
    void setDiscriminatorCombSecVtxMVA( const double & DISCRIMINATORCOMBSECVTXMVA   ) { discriminatorsVec_[COMBSECVTXMVA] = DISCRIMINATORCOMBSECVTXMVA; }
    void setDiscriminatorSoftMuon     ( const double & DISCRIMINATORSOFTMUON	    ) { discriminatorsVec_[SOFTMUON]      = DISCRIMINATORSOFTMUON;      }
    void setDiscriminatorSoftElectron ( const double & DISCRIMINATORSOFTELECTRON    ) { discriminatorsVec_[SOFTELECTRON]  = DISCRIMINATORSOFTELECTRON;  }
    void setDiscriminatorJetProb      ( const double & DISCRIMINATORJETPROB         ) { discriminatorsVec_[JETPROB]       = DISCRIMINATORJETPROB;       }
    void setDiscriminator( const double & DISCRIMINATOR, const std::string& BTAGGER ) { int index = bTaggerCode(BTAGGER);
    discriminatorsVec_[index] = DISCRIMINATOR;
    }

  protected:
    std::vector<double> discriminatorsVec_;
    
  };

  typedef std::vector<BTagDiscrVec> BTagDiscrVecCollection;
  typedef edm::Ref<BTagDiscrVecCollection> BTagDiscrVecRef;
  typedef edm::RefProd<BTagDiscrVecCollection> BTagDiscrVecRefProd;
  typedef edm::RefVector<BTagDiscrVecCollection> BTagDiscrVecRefVector;

}

#endif // BTAGDISCRVEC_H

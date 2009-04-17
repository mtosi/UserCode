#ifndef VBFHiggsToZZto2l2b_CorJetBTagDiscrAssociation_h
#define VBFHiggsToZZto2l2b_CorJetBTagDiscrAssociation_h

/** \class CorJetBTagDiscrAssociation
 *
 * \short Association between jets and b tag jet information
 ************************************************************/

#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/RefVector.h"

namespace fwlite {
  class Event;
}

namespace vbfhzz2l2b {

  typedef std::vector<double> DiscriminatorVector; 

  namespace CorJetBTagDiscrAssociation {
    class CorBTagDiscrData;

    typedef vbfhzz2l2b::CorJetBTagDiscrAssociation::CorBTagDiscrData Value;
    typedef std::vector<Value> Values;
    typedef edm::AssociationVector<reco::JetRefBaseProd, Values> Container;
    typedef Container::value_type value_type;
    typedef Container::transient_vector_type transient_vector_type;
    typedef edm::Ref <Container> Ref;
    typedef edm::RefProd <Container> RefProd;
    typedef edm::RefVector <Container> RefVector;


    /// jet corrected Et
    double jetCorEt (const Container&, const reco::JetBaseRef&);
    double jetCorEt (const Container&, const reco::Jet&);
    /// vector of b-tag discriminators associated to jet
    const DiscriminatorVector& discriminatorVec (const Container&, const reco::JetBaseRef&);
    const DiscriminatorVector& discriminatorVec (const Container&, const reco::Jet&);
    /// associate jet with value. 
    /// Returns false and associate nothing if jet is already associated
    bool setValue (Container&, const reco::JetBaseRef&, const CorBTagDiscrData&);
    bool setValue (Container*, const reco::JetBaseRef&, const CorBTagDiscrData&);
    /// get value for the association. 
    /// Throw exception if no association found
    const CorBTagDiscrData& getValue (const Container&, const reco::JetBaseRef&);
    const CorBTagDiscrData& getValue (const Container&, const reco::Jet&);
    /// fill list of all jets associated with values. 
    /// return # of jets in the list
    std::vector<reco::JetBaseRef > allJets (const Container&);
    /// check if jet is associated
    bool hasJet (const Container&, const reco::JetBaseRef&);
    bool hasJet (const Container&, const reco::Jet&);


    class CorBTagDiscrData {
    public:
      CorBTagDiscrData ();
      ~CorBTagDiscrData () {}
      double highEffDiscr_;
      double_t highPurDiscr_;
      double compoSVDiscr_;
      double jetProbDiscr_;
      DiscriminatorVector discrVec_;
      double corEt_;
      double corPt_;
    };
  }
}

#endif

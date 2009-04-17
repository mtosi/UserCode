#include "DataFormats/JetReco/src/JetAssociationTemplate.icc"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/CorJetBTagDiscrAssociation.h"

    /// jet corrected Et
double vbfhzz2l2b::CorJetBTagDiscrAssociation::jetCorEt (const Container& fContainer, const reco::JetBaseRef& fJet) {
  return getValue (fContainer, fJet).corEt_;
}
double vbfhzz2l2b::CorJetBTagDiscrAssociation::jetCorEt (const Container& fContainer, const reco::Jet& fJet) {
  return getValue (fContainer, fJet).corEt_;
}
/// vector of b-tag discriminators associated to jet
const vbfhzz2l2b::DiscriminatorVector& vbfhzz2l2b::CorJetBTagDiscrAssociation::discriminatorVec (const Container& fContainer, const reco::JetBaseRef& fJet) {
  return getValue (fContainer, fJet).discrVec_;
}
const vbfhzz2l2b::DiscriminatorVector& vbfhzz2l2b::CorJetBTagDiscrAssociation::discriminatorVec (const Container& fContainer, const reco::Jet& fJet) {
  return getValue (fContainer, fJet).discrVec_;
}
/// associate jet with value. 
/// Returns false and associate nothing if jet is already associated
bool vbfhzz2l2b::CorJetBTagDiscrAssociation::setValue (Container& fContainer, 
						       const reco::JetBaseRef& fJet, 
						       const vbfhzz2l2b::CorJetBTagDiscrAssociation::CorBTagDiscrData& fValue) {
  return JetAssociationTemplate::setValue (fContainer, fJet, fValue);
}
bool vbfhzz2l2b::CorJetBTagDiscrAssociation::setValue (Container* fContainer, 
						       const reco::JetBaseRef& fJet, 
						       const vbfhzz2l2b::CorJetBTagDiscrAssociation::CorBTagDiscrData& fValue) {
  return JetAssociationTemplate::setValue (fContainer, fJet, fValue);
}

const vbfhzz2l2b::CorJetBTagDiscrAssociation::CorBTagDiscrData& 
vbfhzz2l2b::CorJetBTagDiscrAssociation::getValue (const Container& fContainer, 
					const reco::Jet& fJet) {
  return JetAssociationTemplate::getValue<Container, Value> (fContainer, fJet);
}
/// get value for the association. 
/// Throw exception if no association found
const vbfhzz2l2b::CorJetBTagDiscrAssociation::CorBTagDiscrData& 
vbfhzz2l2b::CorJetBTagDiscrAssociation::getValue (const Container& fContainer, 
						      const reco::JetBaseRef& fJet) {
  return JetAssociationTemplate::getValue<Container, Value> (fContainer, fJet);
}

bool vbfhzz2l2b::CorJetBTagDiscrAssociation::hasJet (const Container& fContainer, 
						     const reco::JetBaseRef& fJet) {
  return JetAssociationTemplate::hasJet (fContainer, fJet);
}
bool vbfhzz2l2b::CorJetBTagDiscrAssociation::hasJet (const Container& fContainer, 
							 const reco::Jet& fJet) {
  return JetAssociationTemplate::hasJet (fContainer, fJet);
}
std::vector<reco::JetBaseRef > vbfhzz2l2b::CorJetBTagDiscrAssociation::allJets (const Container& fContainer) {
  return JetAssociationTemplate::allJets (fContainer);
}

vbfhzz2l2b::CorJetBTagDiscrAssociation::CorBTagDiscrData::CorBTagDiscrData () 
  : corEt_(0.), corPt_(0.)
{
  //  std::vector<double> null;
  //  DiscriminatorVector discrVec_(null);
}

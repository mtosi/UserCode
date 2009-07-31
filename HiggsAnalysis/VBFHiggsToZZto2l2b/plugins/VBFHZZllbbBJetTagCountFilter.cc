#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetFloatAssociation.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbBJetTagCountFilter.h"

#include <cmath>

VBFHZZllbbBJetTagCountFilter::VBFHZZllbbBJetTagCountFilter(const edm::ParameterSet& cfg)
{
    src_ = cfg.getParameter<edm::InputTag>("src");
    minDiscriminator_ = cfg.getParameter<double>("minDiscriminator");
    minJetEt_ = cfg.getParameter<double>("minJetEt");
    maxJetEta_ = cfg.getParameter<double>("maxJetEta");
    minNumber_ = cfg.getParameter<unsigned int>("minNumber");
}

// --------------------------------------------------------------------------------

VBFHZZllbbBJetTagCountFilter::~VBFHZZllbbBJetTagCountFilter()
{
}

// --------------------------------------------------------------------------------

bool
VBFHZZllbbBJetTagCountFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    unsigned int numberOfJets = 0;
    edm::Handle<reco::JetFloatAssociation::Container> jFA;
    iEvent.getByLabel(src_, jFA);
    std::vector<reco::JetBaseRef> jets = reco::JetFloatAssociation::allJets(*jFA);

    for( std::vector<reco::JetBaseRef>::const_iterator jet = jets.begin(); jet != jets.end(); ++jet )
    {
        if (
            (reco::JetFloatAssociation::getValue(*jFA, **jet) > minDiscriminator_) &&
            ((*jet)->et() > minJetEt_) &&
            (fabs((*jet)->eta()) < maxJetEta_)
            ) numberOfJets++;
    }
    if ( numberOfJets < minNumber_ ) return false;
    return true;
}

// --------------------------------------------------------------------------------


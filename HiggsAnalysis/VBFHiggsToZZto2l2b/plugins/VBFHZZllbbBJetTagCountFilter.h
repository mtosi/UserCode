#ifndef VBFHZZllbbBJetTagCountFilter_h
#define VBFHZZllbbBJetTagCountFilter_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

class VBFHZZllbbBJetTagCountFilter : public edm::EDFilter {
public:
    VBFHZZllbbBJetTagCountFilter(const edm::ParameterSet&);
    ~VBFHZZllbbBJetTagCountFilter();
private:
    virtual bool filter(edm::Event&, const edm::EventSetup&);
    edm::InputTag src_;
    double minDiscriminator_;
    double minJetEt_;
    double maxJetEta_;
    unsigned int minNumber_;
};

#endif

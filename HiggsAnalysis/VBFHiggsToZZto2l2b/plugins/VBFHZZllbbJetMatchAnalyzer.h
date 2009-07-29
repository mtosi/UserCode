// system include files
#include <memory>
#include <string>
#include <iostream>
#include <algorithm>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/CandMatchMapMany.h"

//  to access TFileService within a Module
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include <Math/VectorUtil.h>

using namespace std;
using namespace reco;
using namespace edm;
using namespace ROOT::Math::VectorUtil;

class VBFHZZllbbJetMatchAnalyzer : public edm::EDAnalyzer {
  public:
    explicit VBFHZZllbbJetMatchAnalyzer(const edm::ParameterSet&);
    ~VBFHZZllbbJetMatchAnalyzer() {}
     virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
     virtual void beginJob(const edm::EventSetup& iSetup) ;
     virtual void endJob() ;

  private:

    InputTag source_;
    InputTag matched_;
    InputTag matchedjetsOne_;   
    InputTag matchedjetsMany_;   

    Handle<CandidateCollection> source;
    Handle<CandidateCollection> matched;
    Handle<CandViewMatchMap>    matchedjetsOne;
    Handle<CandMatchMapMany>    matchedjetsMany;


    unsigned eventcounter_;

    edm::Service<TFileService> fs;
    TH1D * deltaR_;
    TH1D * deltaPt_;
    TH1D * resPt_; 
    TH2D * deltaPtVSdeltaR_;
    TH2D * resPtVSdeltaR_; 
};

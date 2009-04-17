#ifndef HiggsAnalysis_VBFHiggsToZZto2l2b_interface_VBFHelper_h
#define HiggsAnalysis_VBFHiggsToZZto2l2b_interface_VBFHelper_h




//-------------------------------------------------------------------------------------
//
// Original Author:  Salvatore Rappoccio
//         Created:  Mon Jul  7 10:37:27 CDT 2008
//-------------------------------------------------------------------------------------
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/PhysicsHistograms.h"

#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "PhysicsTools/UtilAlgos/interface/TFileDirectory.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/EDProducer.h"


namespace vbfhzz2l2b {

  class VBFHelper {

  public:

    VBFHelper(edm::ParameterSet const & parameters);
    ~VBFHelper();

    // Pull out a struct for the axis limits from the config file
    PhysicsHistograms::KinAxisLimits getAxisLimits( std::string name );
						    

    // Book histograms
    void bookHistos( edm::EDProducer * producer );

    // Get handles
    void getHandles( edm::Event  & event,
		     edm::Handle<reco::MuonCollection >           & muonHandle,	    
		     edm::Handle<reco::GsfElectronCollection >    & electronHandle, 
		     edm::Handle<reco::CaloJetCollection >        & caloJetHandle,  	
		     edm::Handle<reco::CaloMETCollection >        & caloMETHandle,  
		     edm::Handle<std::vector<reco::GenParticle> > & genParticlesHandle    
		     );

    
    // fill histograms
    void fillHistograms( edm::Event & event,
			 edm::Handle<reco::MuonCollection >           & muonHandle,
			 edm::Handle<reco::GsfElectronCollection >    & electronHandle,
			 edm::Handle<reco::CaloJetCollection >        & caloJetHandle,
			 edm::Handle<reco::CaloMETCollection >        & caloMETHandle,
			 edm::Handle<std::vector<reco::GenParticle> > & genParticlesHandle
			 );
    
    
    // Function to add ntuple variables to the EDProducer
    void addNtupleVar ( edm::EDProducer * prod, std::string name, std::string type );

    // Save ntuple variables to event evt
    void saveNtuple (  edm::Event & event,
		       const std::vector<vbfhzz2l2b::PhysVarHisto*> & ntvars);
    
    // Helper function template to write objects to event
    template <class T>
      void saveNtupleVar(  edm::Event & event,
			   std::string name, T value);

    // Helper function template to write vectors of objects to event
    template <class T>
      void saveNtupleVec(  edm::Event & event,
			   std::string name, const std::vector<T> & invec);

    
    // verbose switch
    int verboseLevel_;

    // Keep a version of the parameter set in question
    edm::ParameterSet         parameters_;

    // Here is where the histograms go
    PhysicsHistograms  *      physHistos_;

    // File service for histograms
    edm::Service<TFileService> fs_;
    
    // List of ntuple variables
    std::vector< vbfhzz2l2b::PhysVarHisto* > ntVars_ ;

    
    // run and event numbers
    vbfhzz2l2b::PhysVarHisto *  h_runNumber_;
    vbfhzz2l2b::PhysVarHisto *  h_eventNumber_;
    
  };

}


#endif

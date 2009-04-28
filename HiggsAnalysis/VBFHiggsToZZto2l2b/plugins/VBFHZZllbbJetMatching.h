// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

//  to access TFileService within a Module
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

// Root includes
// -------------
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TMath.h"

#include <cmath>

//
// class declaration
//
class VBFHZZllbbJetMatching : public edm::EDAnalyzer {

 public:
  explicit VBFHZZllbbJetMatching(const edm::ParameterSet&);
  ~VBFHZZllbbJetMatching();

 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  int eventcounter_;

  int nbin_;

  // ----------member data ---------------------------
  edm::InputTag electronLabel_;
  edm::InputTag muonLabel_;
  edm::InputTag metLabel_;
  edm::InputTag jetLabel_;
  std::string   corJetsWithBTagLabel_;
  edm::InputTag mcParticleLabel_;
  edm::InputTag genJetLabel_;
  edm::InputTag genMetLabel_;

  //  edm::Service<TFileService> fs;

  TH1D * eventsNumber_;

  TH1D * ZpartonEta_;	       
  TH1D * ZpartonPt_;	       
  TH1D * ZpartonEt_;	       
  TH1D * ZpartonE_;	       
  TH1D * ZpartonsDeltaEta_;    
  TH1D * ZpartonsDeltaR_;      
  TH1D * ZpartonsMass_;	       
  TH1D * ZpartonsPt_;	       
  TH1D * ZpartonsEt_;	       
  TH1D * ZpartonsE_;	       
  TH1D * ZpartonsCollinearity_;

  TH2D * ZpartonsDeltaRVSjetMass_;
  TH2D * ZpartonsEtResVSjetMass_;
  TH2D * ZpartonsDeltaRVSetRes_;
  TProfile * ZpartonsDeltaRVSjetMass_profile_;
  TProfile * ZpartonsEtResVSjetMass_profile_;
  TProfile * ZpartonsDeltaRVSetRes_profile_;
 
  TH1D * TAGpartonEta_;
  TH1D * TAGpartonPt_;
  TH1D * TAGpartonEt_;
  TH1D * TAGpartonE_;
  TH1D * TAGpartonsDeltaEta_;
  TH1D * TAGpartonsDeltaR_;
  TH1D * TAGpartonsMass_;
  TH1D * TAGpartonsPt_;
  TH1D * TAGpartonsEt_;
  TH1D * TAGpartonsE_;
  TH1D * TAGpartonsCollinearity_;

  TH1D * jetNumber_;
  TH1D * jetUncorrEt_;
  TH1D * jetCorrEt_;
  TH1D * jetUncorrPt_;
  TH1D * jetCorrPt_;
  TH1D * jetPhi_;
  TH1D * jetEta_;
  TH1D * jetMass_;
  TH1D * jetEMfrac_;
  TH1D * jetHIGHEFFdiscr_;   	  
  TH1D * jetHIGHPURdiscr_;   	  
  TH1D * jetCOMBSECVTXdiscr_;
  TH1D * jetJETPROBdiscr_;   

  TH2D * jetEMfracVSeta_;	      
  TH2D * jetEMfracVScorrEt_;	      
  TH2D * jetEMfracVScorrPt_;	      
  TH2D * jetEMfracVSuncorrEt_;	      
  TH2D * jetEMfracVSuncorrPt_;	      
  TH2D * jetEMfracVShighEFFdiscr_;    
  TH2D * jetEMfracVShighPURdiscr_;    
  TH2D * jetEMfracVScomboSECVTXdiscr_;
  TH2D * jetEMfracVSjetPROBdiscr_;    
  TProfile * jetEMfracVSeta_profile_;	      
  TProfile * jetEMfracVScorrEt_profile_;	      
  TProfile * jetEMfracVScorrPt_profile_;	      
  TProfile * jetEMfracVSuncorrEt_profile_;	      
  TProfile * jetEMfracVSuncorrPt_profile_;	      
  TProfile * jetEMfracVShighEFFdiscr_profile_;    
  TProfile * jetEMfracVShighPURdiscr_profile_;    
  TProfile * jetEMfracVScomboSECVTXdiscr_profile_;
  TProfile * jetEMfracVSjetPROBdiscr_profile_;    


  TH1D * muonJetMatchedNumber_;

  TH1D * jetParton_deltaR_;
  TH1D * jetParton_deltaR_zoom_;
  TH1D * jetParton_deltaRmax_;

  TH2D * jetParton_deltaEtVSdeltaR_;
  TH2D * jetParton_deltaEtmeanVSdeltaRmean_;
  TH2D * jetParton_deltaEVSdeltaR_;
  TH2D * jetParton_deltaEtaVSdeltaR_;
  TProfile * jetParton_deltaEtVSdeltaR_profile_;
  TProfile * jetParton_deltaEtmeanVSdeltaRmean_profile_;
  TProfile * jetParton_deltaEVSdeltaR_profile_;
  TProfile * jetParton_deltaEtaVSdeltaR_profile_;
  TH1D * jetParton_deltaEta_;
  TH1D * jetParton_deltaEt_;
  TH1D * jetParton_deltaEtRes_;
  TH1D * jetParton_004deltaEtRes_;
  TH2D * jetParton_deltaRVSdeltaR_;
  TProfile * jetParton_deltaRVSdeltaR_profile_;

  TH1D * ZjetEta_;
  TH1D * ZjetPt_;
  TH1D * ZjetEt_;
  TH1D * ZjetE_;
  TH1D * ZjetMass_;
  TH1D * ZjetsDeltaEta_;
  TH1D * ZjetsDeltaR_;
  TH1D * ZjetsMass_;
  TH1D * ZjetsMassResolution_;
  TH1D * ZjetsPt_;
  TH1D * ZjetsPtResolution_;
  TH1D * ZjetsE_;
  TH1D * ZjetsEResolution_;
  TH1D * ZjetsCollinearity_;
  TH1D * ZjetsCollinearityResolution_;

  TH1D * TAGjetEta_;
  TH1D * TAGjetPt_;
  TH1D * TAGjetEt_;
  TH1D * TAGjetE_;
  TH1D * TAGjetMass_;
  TH1D * TAGjetsDeltaEta_;
  TH1D * TAGjetsDeltaR_;
  TH1D * TAGjetsMass_;
  TH1D * TAGjetsMassResolution_;
  TH1D * TAGjetsPt_;
  TH1D * TAGjetsPtResolution_;
  TH1D * TAGjetsE_;
  TH1D * TAGjetsEResolution_;
  TH1D * TAGjetsCollinearity_;
  TH1D * TAGjetsCollinearityResolution_;

};

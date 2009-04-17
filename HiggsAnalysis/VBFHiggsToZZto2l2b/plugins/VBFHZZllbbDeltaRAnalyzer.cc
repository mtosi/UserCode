// -*- C++ -*-
//
// Package:    VBFHZZllbbDeltaRAnalyzer
// Class:      VBFHZZllbbDeltaRAnalyzer
// 
/**\class VBFHZZllbbDeltaRAnalyzer VBFHZZllbbDeltaRAnalyzer.cc HiggsAnalysis/VBFHiggsToZZto2l2b/src/VBFHZZllbbDeltaRAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Mia TOSI
//         Created:  Mon Feb  2 17:31:44 CET 2009
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbDeltaRAnalyzer.h"

#include "PhysicsTools/UtilAlgos/interface/TFileDirectory.h"
#include "PhysicsTools/UtilAlgos/interface/TH1AddDirectorySentry.h" 

#include "PhysicsTools/Utilities/interface/deltaR.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TH1F.h"

// class decleration
//

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
VBFHZZllbbDeltaRAnalyzer::VBFHZZllbbDeltaRAnalyzer(const edm::ParameterSet& iConfig) :
  genParticleLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "genParticleLabel" ) ),
  caloJetLabel_    ( iConfig.getUntrackedParameter<edm::InputTag>( "caloJetLabel"     ) ),
  muonLabel_       ( iConfig.getUntrackedParameter<edm::InputTag>( "muonLabel"        ) ),
  electronLabel_   ( iConfig.getUntrackedParameter<edm::InputTag>( "electronLabel"    ) ),
  bTagLabel_       ( iConfig.getUntrackedParameter<edm::InputTag>( "bTagLabel"        ) ),
  metLabel_        ( iConfig.getUntrackedParameter<edm::InputTag>( "metLabel"         ) ),
  ptMax_( iConfig.getUntrackedParameter<double>( "ptMax" ) ) 
{
   //now do what ever initialization is needed

}


VBFHZZllbbDeltaRAnalyzer::~VBFHZZllbbDeltaRAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
VBFHZZllbbDeltaRAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;

  edm::Handle<reco::CaloJetCollection> caloJetHandle;
  iEvent.getByLabel(caloJetLabel_,caloJetHandle);

  edm::Handle<reco::SecondaryVertexTagInfo> combinedSecondaryVertexTagInfosHandle;
  edm::Handle<reco::JetTagCollection> combinedSecondaryVertexHandle; 
  iEvent.getByLabel(bTagLabel_,combinedSecondaryVertexHandle);

  edm::Handle<reco::MuonCollection> muonHandle;
  iEvent.getByLabel(muonLabel_,muonHandle);

  edm::Handle<reco::GsfElectronCollection> electronHandle;
  iEvent.getByLabel(electronLabel_,electronHandle);

  edm::Handle<reco::CaloMETCollection> metHandle;
  iEvent.getByLabel(metLabel_,metHandle);

  edm::Handle<reco::VertexCollection> primaryVerticesWithBSHandle;
  iEvent.getByLabel("offlinePrimaryVerticesWithBS",primaryVerticesWithBSHandle);
  edm::Handle<reco::VertexCollection> primaryVerticesHandle;
  iEvent.getByLabel("offlinePrimaryVertices",primaryVerticesHandle);
  edm::Handle<reco::VertexCollection> pixelVerticesHandle;
  iEvent.getByLabel("pixelVertices",pixelVerticesHandle);

  std::cout << "primaryVerticesWithBSHandle->size(): " << primaryVerticesWithBSHandle->size() << std::endl;
  std::cout << "primaryVerticesHandle->size(): " << primaryVerticesHandle->size() << std::endl;
  std::cout << "pixelVerticesHandle->size(): " << pixelVerticesHandle->size() << std::endl;
  
  reco::VertexCollection::const_iterator primaryVertex_itr = primaryVerticesWithBSHandle->begin();
  double primaryVertexZ = primaryVertex_itr->z();

  reco::CaloJetCollection::const_iterator jet_itr = caloJetHandle->begin(); 
  for ( ; jet_itr != caloJetHandle->end(); ++jet_itr ) {
    h_jet_pt->    Fill( jet_itr->pt()  );
    h_jet_eta->   Fill( jet_itr->eta() );
    h_jet_phi->   Fill( jet_itr->phi() );
    double res = 0.;
    if (jet_itr->eta() != 0. ) res = (jet_itr->eta()-jet_itr->physicsEtaQuick(primaryVertexZ))/jet_itr->eta();
    h_jet_etaVSphysicsEtaResolution->Fill(res);
    h_jet_etaVSphysicsEta->Fill(jet_itr->eta(),jet_itr->physicsEtaQuick(primaryVertexZ));
    p_jet_etaVSphysicsEta->Fill(jet_itr->eta(),jet_itr->physicsEtaQuick(primaryVertexZ));
  }
  reco::MuonCollection::const_iterator muon_itr = muonHandle->begin(); 
  for ( ; muon_itr != muonHandle->end(); ++muon_itr ) {
    h_muon_pt ->Fill( muon_itr->pt()  );
    h_muon_eta->Fill( muon_itr->eta() );
    h_muon_phi->Fill( muon_itr->phi() );
  }

  reco::GsfElectronCollection::const_iterator electron_itr = electronHandle->begin();
  for ( ; electron_itr != electronHandle->end(); ++electron_itr ) {
    h_electron_pt-> Fill( electron_itr->pt()  );
    h_electron_eta->Fill( electron_itr->eta() );
    h_electron_phi->Fill( electron_itr->phi() );
  }

  const reco::CaloMETCollection * metCollection = metHandle.product();
  reco::CaloMETCollection::const_iterator met = metCollection->begin();
  for ( ; met != metCollection->end(); ++met ) {
    h_met_pt->Fill(  met->pt()  );
    h_met_eta->Fill( met->eta() );
    h_met_phi->Fill( met->phi() );
  }

  const reco::JetTagCollection & bTagCollection = *(combinedSecondaryVertexHandle.product());
  reco::JetTagCollection::const_iterator bTagCollection_itr = bTagCollection.begin();
  for ( ; bTagCollection_itr != bTagCollection.end(); ++bTagCollection_itr ) {
    if (bTagCollection_itr->second > 0) { 
      muon_itr = muonHandle->begin(); 
      for ( ; muon_itr != muonHandle->end(); ++muon_itr ) {
	double deltaR = DeltaR<reco::Muon, reco::Jet>()(*muon_itr, *bTagCollection_itr->first);
	h_bjetmuon_deltaR->Fill(deltaR);
      }
      electron_itr = electronHandle->begin(); 
      for ( ; electron_itr != electronHandle->end(); ++electron_itr ) {
	double deltaR = DeltaR<reco::GsfElectron, reco::Jet>()(*electron_itr, *bTagCollection_itr->first);
	h_bjetelectron_deltaR->Fill(deltaR);
      }
      met = metCollection->begin();
      for ( ; met != metCollection->end(); ++met ) {
	double deltaR = DeltaR<reco::MET, reco::Jet>()(*met, *bTagCollection_itr->first);
	h_bjetmet_deltaR->Fill(deltaR);
      }
    }
    h_jet_bDiscr->Fill(bTagCollection_itr->second);
  }
  
  jet_itr = caloJetHandle->begin(); 
  for ( ; jet_itr != caloJetHandle->end(); ++jet_itr ) {
    muon_itr = muonHandle->begin(); 
    for ( ; muon_itr != muonHandle->end(); ++muon_itr ) {
      double deltaR = DeltaR<reco::Muon, reco::CaloJet>()(*muon_itr, *jet_itr);
      h_jetmuon_deltaR->Fill(deltaR);
    }
    electron_itr = electronHandle->begin(); 
    for ( ; electron_itr != electronHandle->end(); ++electron_itr ) {
      double deltaR = DeltaR<reco::GsfElectron, reco::CaloJet>()(*electron_itr, *jet_itr);
      h_jetelectron_deltaR->Fill(deltaR);
    }
    met = metCollection->begin();
    for ( ; met != metCollection->end(); ++met ) {
      double deltaR = DeltaR<reco::MET, reco::CaloJet>()(*met, *jet_itr);
      h_jetmet_deltaR->Fill(deltaR);
    }
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
VBFHZZllbbDeltaRAnalyzer::beginJob(const edm::EventSetup&)
{
  TFileDirectory jetSubDir = fs->mkdir( "jet" );
  h_jet_pt     = jetSubDir.make<TH1D>( "jet_pt",     "jet p_{T}",                100,  0.  , 200.   );
  h_jet_eta    = jetSubDir.make<TH1D>( "jet_eta",    "jet #eta",                 100,-10.  ,  10.   );
  h_jet_phi    = jetSubDir.make<TH1D>( "jet_phi",    "jet #phi",                 100, -3.14,   3.14 );
  h_jet_bDiscr = jetSubDir.make<TH1D>( "jet_bDiscr", "jet b tag discriminator",  100,-10.,     2.   );
  h_jetmuon_deltaR     = jetSubDir.make<TH1D>( "jetmuon_deltaR",     "#Delta R_{j,#mu}",  50,  0.,  10. );
  h_jetelectron_deltaR = jetSubDir.make<TH1D>( "jetelectron_deltaR", "#Delta R_{j,e}",    50,  0.,  10. );
  h_jetmet_deltaR      = jetSubDir.make<TH1D>( "jetmet_deltaR",      "#Delta R_{j,met}",  50,  0.,  10. );

  h_jet_etaVSphysicsEtaResolution = jetSubDir.make<TH1D>( "jet_etaVSphysicsEtaResolution","jet #eta VS physicsEtaQuick resolution",40,-1.,1.);
  h_jet_etaVSphysicsEta = jetSubDir.make<TH2D>( "jet_etaVSphysicsEta","jet #eta VS physicsEtaQuick",100,-10.,10.,100,-10.,10.);
  p_jet_etaVSphysicsEta = jetSubDir.make<TProfile>( "jet_etaVSphysicsEta_profile","jet #eta VS physicsEtaQuick",100,-10.,10.,-10.,10.);

  TFileDirectory bjetSubDir = jetSubDir.mkdir( "bjet" );
  h_bjetmuon_deltaR     = bjetSubDir.make<TH1D>( "bjetmuon_deltaR",     "#Delta R_{j_{b},#mu}",  50,  0.,    10. );
  h_bjetelectron_deltaR = bjetSubDir.make<TH1D>( "bjetelectron_deltaR", "#Delta R_{j_{b},e}",    50,  0.,    10. );
  h_bjetmet_deltaR      = bjetSubDir.make<TH1D>( "bjetmet_deltaR",      "#Delta R_{j_{b},met}",  50,  0.,    10. );

  TFileDirectory muonSubDir = fs->mkdir( "muon" );
  h_muon_pt  = muonSubDir.make<TH1D>( "muon_pt",  "muon p_{T}", 100,  0.  , 200.   );
  h_muon_eta = muonSubDir.make<TH1D>( "muon_eta", "muon #eta",  100,-10.  ,  10.   );
  h_muon_phi = muonSubDir.make<TH1D>( "muon_phi", "muon #phi",  100, -3.14,   3.14 );
  
  TFileDirectory electronSubDir = fs->mkdir( "electron" );
  h_electron_pt  = electronSubDir.make<TH1D>( "electron_pt",  "electron p_{T}", 100,  0.  , 200.   );
  h_electron_eta = electronSubDir.make<TH1D>( "electron_eta", "electron #eta",  100,-10.  ,  10.   );
  h_electron_phi = electronSubDir.make<TH1D>( "electron_phi", "electron #phi",  100, -3.14,   3.14 );
  
  TFileDirectory metSubDir = fs->mkdir( "met" );
  h_met_pt  = metSubDir.make<TH1D>( "met_pt",  "met p_{T}", 100,  0.  , 200.   );
  h_met_eta = metSubDir.make<TH1D>( "met_eta", "met #eta",  100,-10.  ,  10.   );
  h_met_phi = metSubDir.make<TH1D>( "met_phi", "met #phi",  100, -3.14,   3.14 );
  
  //  TFileDirectory subSubDir = subDir.mkdir( "mySubSubDirectory" );
  //  TH1D * h_pt = subDir.make<TH1D>( "pt"  , "p_{t}", 100,  0., 100. );

}

// ------------ method called once each job just after ending the event loop  ------------
void 
VBFHZZllbbDeltaRAnalyzer::endJob() {
}

//define this as a plug-in
//DEFINE_FWK_MODULE(VBFHZZllbbDeltaRAnalyzer);

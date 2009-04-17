// -*- C++ -*-
//
// Package:    VBFHZZllbbRAWCORBTAGtesting
// Class:      VBFHZZllbbRAWCORBTAGtesting
// 
/**\class VBFHZZllbbRAWCORBTAGtesting VBFHZZllbbRAWCORBTAGtesting.cc HiggsAnalysis/VBFHiggsToZZto2l2b/src/VBFHZZllbbRAWCORBTAGtesting.cc

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

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbRAWCORBTAGtesting.h"

#include "PhysicsTools/UtilAlgos/interface/TFileDirectory.h"
#include "PhysicsTools/UtilAlgos/interface/TH1AddDirectorySentry.h" 

#include "PhysicsTools/Utilities/interface/deltaR.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

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
VBFHZZllbbRAWCORBTAGtesting::VBFHZZllbbRAWCORBTAGtesting(const edm::ParameterSet& iConfig) :
  caloJetLabel_    ( iConfig.getUntrackedParameter<edm::InputTag>( "jetLabel"      ) ),
  corrCaloJetLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "corrJetLabel"  ) ),
  muonLabel_       ( iConfig.getUntrackedParameter<edm::InputTag>( "muonLabel"     ) ),
  electronLabel_   ( iConfig.getUntrackedParameter<edm::InputTag>( "electronLabel" ) ),
  metLabel_        ( iConfig.getUntrackedParameter<edm::InputTag>( "metLabel"      ) ),
//  bTagLabel_       ( iConfig.getUntrackedParameter<edm::InputTag>( "bTagLabel"        ) ),
  bTagConfigLabel_ ( iConfig.getParameter< vector<edm::ParameterSet> >( "bTagConfig" ) ) ,
  jetCorrectionServiceLabel_ ( iConfig.getParameter<std::string>( "jetCorrectionService" ) )
{
   //now do what ever initialization is needed

}

VBFHZZllbbRAWCORBTAGtesting::~VBFHZZllbbRAWCORBTAGtesting()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
VBFHZZllbbRAWCORBTAGtesting::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;

  edm::Handle<reco::CaloJetCollection> caloJetHandle;
  iEvent.getByLabel(caloJetLabel_,caloJetHandle);

  edm::Handle<reco::CaloJetCollection> corrCaloJetHandle;
  iEvent.getByLabel(corrCaloJetLabel_,corrCaloJetHandle);

//  for (unsigned int iModule = 0; iModule != bTagConfigLabel_.size(); ++iModule) {
//    string dataFormatType = "JetTag";
//    InputTag bTagLabel = bTagConfigLabel_[iModule].getParameter<InputTag>("label");
//    bJetTagInputTags_.push_back( bTagLabel );
//  }

  edm::Handle<reco::MuonCollection> muonHandle;
  iEvent.getByLabel(muonLabel_,muonHandle);

  edm::Handle<reco::GsfElectronCollection> electronHandle;
  iEvent.getByLabel(electronLabel_,electronHandle);

  edm::Handle<reco::CaloMETCollection> metHandle;
  iEvent.getByLabel(metLabel_,metHandle);

  reco::MuonCollection::const_iterator muon_itr = muonHandle->begin(); 
  reco::GsfElectronCollection::const_iterator electron_itr = electronHandle->begin();
  const reco::CaloMETCollection * metCollection = metHandle.product();
  reco::CaloMETCollection::const_iterator met = metCollection->begin();

  const JetCorrector* corrector = JetCorrector::getJetCorrector (jetCorrectionServiceLabel_,iSetup);

  double corScale = 1.;

  //  edm::Handle<reco::SecondaryVertexTagInfo> combinedSecondaryVertexTagInfosHandle;
  //  edm::Handle<reco::JetTagCollection> combinedSecondaryVertexHandle; 
  //  iEvent.getByLabel(bTagLabel_,combinedSecondaryVertexHandle);
  for ( unsigned int bJetLabelIndex = 0; bJetLabelIndex != bJetTagInputTags_.size(); bJetLabelIndex++ ) {
    edm::Handle<reco::JetTagCollection> bTagHandle;
    iEvent.getByLabel(bJetTagInputTags_[bJetLabelIndex], bTagHandle);

    const reco::JetTagCollection & bTagCollection = *(bTagHandle.product());
    std::cout << "[VBFHZZllbbRAWCORBTAGtesting::analyze] Found " << bTagCollection.size() 
	      << " B candidates in collection " << bJetTagInputTags_[bJetLabelIndex] << std::endl;
    reco::JetTagCollection::const_iterator bTagCollection_itr = bTagCollection.begin();
    for ( ; bTagCollection_itr != bTagCollection.end(); ++bTagCollection_itr ) {
      std::cout << "bJetLabelIndex: " << bJetLabelIndex << " pt: " << (*bTagCollection_itr->first).pt() << std::endl;
      hvec_bjet_pt [bJetLabelIndex]->Fill((*bTagCollection_itr->first).pt() );
      hvec_bjet_eta[bJetLabelIndex]->Fill((*bTagCollection_itr->first).eta());
      hvec_bjet_phi[bJetLabelIndex]->Fill((*bTagCollection_itr->first).phi());
      
      muon_itr = muonHandle->begin(); 
      for ( ; muon_itr != muonHandle->end(); ++muon_itr ) {
	double deltaR = DeltaR<reco::Muon, reco::Jet>()(*muon_itr, *bTagCollection_itr->first);
	hvec_bjetmuon_deltaR[bJetLabelIndex]->Fill(deltaR);
      }
      electron_itr = electronHandle->begin(); 
      for ( ; electron_itr != electronHandle->end(); ++electron_itr ) {
	double deltaR = DeltaR<reco::GsfElectron, reco::Jet>()(*electron_itr, *bTagCollection_itr->first);
	hvec_bjetelectron_deltaR[bJetLabelIndex]->Fill(deltaR);
      }
      met = metCollection->begin();
      for ( ; met != metCollection->end(); ++met ) {
	double deltaR = DeltaR<reco::MET, reco::Jet>()(*met, *bTagCollection_itr->first);
	hvec_bjetmet_deltaR[bJetLabelIndex]->Fill(deltaR);
      }
      hvec_bjet_bDiscr[bJetLabelIndex]->Fill(bTagCollection_itr->second);
      
      const reco::CaloJet * referenceJet = dynamic_cast<const reco::CaloJet*>(&*bTagCollection_itr->first); 
      double uncorrEt = referenceJet->et();
      double emFrac   = referenceJet->emEnergyFraction();
      hvec_bjet_emFrac[bJetLabelIndex]->Fill(emFrac);
      
      corScale = corrector->correction( *referenceJet );
      std::cout << "corScale: " << corScale << std::endl;
      h_bcoronflyjet_pt ->Fill(corScale*(*bTagCollection_itr->first).pt() );
      
    }
  }
  std::cout << "---> done w/ bTagCollection" << std::endl;

//  reco::CaloJetCollection::const_iterator jet_itr = caloJetHandle->begin(); 
//  for ( ; jet_itr != caloJetHandle->end(); ++jet_itr ) {
//    h_jet_pt->    Fill( jet_itr->pt()               );
//    h_jet_eta->   Fill( jet_itr->eta()              );
//    h_jet_phi->   Fill( jet_itr->phi()              );
//    h_jet_emFrac->Fill( jet_itr->emEnergyFraction() );
//
//    muon_itr = muonHandle->begin(); 
//    for ( ; muon_itr != muonHandle->end(); ++muon_itr ) {
//      double deltaR = DeltaR<reco::Muon, reco::CaloJet>()(*muon_itr, *jet_itr);
//      h_jetmuon_deltaR->Fill(deltaR);
//    }
//    electron_itr = electronHandle->begin(); 
//    for ( ; electron_itr != electronHandle->end(); ++electron_itr ) {
//      double deltaR = DeltaR<reco::GsfElectron, reco::CaloJet>()(*electron_itr, *jet_itr);
//      h_jetelectron_deltaR->Fill(deltaR);
//    }
//    met = metCollection->begin();
//    for ( ; met != metCollection->end(); ++met ) {
//      double deltaR = DeltaR<reco::MET, reco::CaloJet>()(*met, *jet_itr);
//      h_jetmet_deltaR->Fill(deltaR);
//    }
//
//    corScale = corrector->correction( *jet_itr );
//    h_coronflyjet_pt ->Fill(corScale*(jet_itr->pt() ) );
//  }
//  std::cout << "---> done w/ CaloCollection" << std::endl;

//  reco::CaloJetCollection::const_iterator corJet_itr = corrCaloJetHandle->begin(); 
//  for ( ; corJet_itr != corrCaloJetHandle->end(); ++corJet_itr ) {
//    h_corjet_pt->  Fill( corJet_itr->pt()  );
//    h_corjet_eta-> Fill( corJet_itr->eta() );
//    h_corjet_phi-> Fill( corJet_itr->phi() );
//
//    muon_itr = muonHandle->begin(); 
//    for ( ; muon_itr != muonHandle->end(); ++muon_itr ) {
//      double deltaR = DeltaR<reco::Muon, reco::CaloJet>()(*muon_itr, *corJet_itr);
//      h_corjetmuon_deltaR->Fill(deltaR);
//    }
//    electron_itr = electronHandle->begin(); 
//    for ( ; electron_itr != electronHandle->end(); ++electron_itr ) {
//      double deltaR = DeltaR<reco::GsfElectron, reco::CaloJet>()(*electron_itr, *corJet_itr);
//      h_corjetelectron_deltaR->Fill(deltaR);
//    }
//    met = metCollection->begin();
//    for ( ; met != metCollection->end(); ++met ) {
//      double deltaR = DeltaR<reco::MET, reco::CaloJet>()(*met, *corJet_itr);
//      h_corjetmet_deltaR->Fill(deltaR);
//    }
//  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
VBFHZZllbbRAWCORBTAGtesting::beginJob(const edm::EventSetup&)
{
  TFileDirectory jetSubDir = fs->mkdir( "jet" );
  h_jet_pt     = jetSubDir.make<TH1D>( "jet_pt",     "jet p_{T}",                100,  0.  , 200.   );
  h_jet_eta    = jetSubDir.make<TH1D>( "jet_eta",    "jet #eta",                 100,-10.  ,  10.   );
  h_jet_phi    = jetSubDir.make<TH1D>( "jet_phi",    "jet #phi",                 100, -3.14,   3.14 );
  h_jet_emFrac = jetSubDir.make<TH1D>( "jet_emFrac", "jet EM fraction",          100, -1.,     1.5  );
  h_jetmuon_deltaR     = jetSubDir.make<TH1D>( "jetmuon_deltaR",     "#Delta R_{j,#mu}",  50,  0.,  10. );
  h_jetelectron_deltaR = jetSubDir.make<TH1D>( "jetelectron_deltaR", "#Delta R_{j,e}",    50,  0.,  10. );
  h_jetmet_deltaR      = jetSubDir.make<TH1D>( "jetmet_deltaR",      "#Delta R_{j,met}",  50,  0.,  10. );
  
  TFileDirectory corJetSubDir = fs->mkdir( "corjet" );
  h_corjet_pt     = corJetSubDir.make<TH1D>( "corjet_pt",     "corrected jet p_{T}",                100,  0.  , 200.   );
  h_corjet_eta    = corJetSubDir.make<TH1D>( "corjet_eta",    "corrected jet #eta",                 100,-10.  ,  10.   );
  h_corjet_phi    = corJetSubDir.make<TH1D>( "corjet_phi",    "corrected jet #phi",                 100, -3.14,   3.14 );
  h_corjet_emFrac = corJetSubDir.make<TH1D>( "corjet_emFrac", "corrected jet EM fraction",          100, -1.,     1.5  );
  h_corjetmuon_deltaR     = corJetSubDir.make<TH1D>( "corjetmuon_deltaR",     "#Delta R_{j_{cor},#mu}",  50,  0.,  10. );
  h_corjetelectron_deltaR = corJetSubDir.make<TH1D>( "corjetelectron_deltaR", "#Delta R_{j_{cor},e}",    50,  0.,  10. );
  h_corjetmet_deltaR      = corJetSubDir.make<TH1D>( "corjetmet_deltaR",      "#Delta R_{j_{cor},met}",  50,  0.,  10. );
  h_coronflyjet_pt     = corJetSubDir.make<TH1D>( "coronflyjet_pt",     "jet p_{T}",                100,  0.  , 200.   );
  h_bcoronflyjet_pt    = corJetSubDir.make<TH1D>( "bcoronflyjet_pt",     "jet p_{T}",                100,  0.  , 200.   );
  
  TFileDirectory bjetSubDir = fs->mkdir( "bjet" );

  for (unsigned int iModule = 0; iModule != bTagConfigLabel_.size(); ++iModule) {
    string dataFormatType = "JetTag";
    edm::InputTag bTagLabel = bTagConfigLabel_[iModule].getParameter<edm::InputTag>("label");
    bJetTagInputTags_.push_back( bTagLabel );

    TString label(bTagLabel.label());
    std::cout << "[VBFHZZllbbRAWCORBTAGtesting::beginJob] label: " << label << std::endl;
    TH1D * h_bjet_pt     = bjetSubDir.make<TH1D>( "bjet_"+label+"_pt",    label+" jet p_{T}", 100,  0.  , 200.   );
    TH1D * h_bjet_eta    = bjetSubDir.make<TH1D>( "bjet_"+label+"_eta",   label+" jet #eta",                 100,  -10.  ,  10.   );
    TH1D * h_bjet_phi    = bjetSubDir.make<TH1D>( "bjet_"+label+"_phi",   label+" jet #phi",                 100,   -3.14,   3.14 );
    TH1D * h_bjet_emFrac = bjetSubDir.make<TH1D>( "bjet_"+label+"_emFrac",label+" jet EM fraction",          100,   -1.,     1.5  );
    TH1D * h_bjet_bDiscr = bjetSubDir.make<TH1D>( "bjet_"+label+"_bDiscr",label+" jet b tag discriminator",  1000,-100.,     2.   );
    TH1D * h_bjetmuon_deltaR     = bjetSubDir.make<TH1D>( "bjet"+label+"muon_deltaR",    label+"#Delta R_{j_{b},#mu}",  50,  0.,    10. );
    TH1D * h_bjetelectron_deltaR = bjetSubDir.make<TH1D>( "bjet"+label+"electron_deltaR",label+"#Delta R_{j_{b},e}",    50,  0.,    10. );
    TH1D * h_bjetmet_deltaR      = bjetSubDir.make<TH1D>( "bjet"+label+"met_deltaR",     label+"#Delta R_{j_{b},met}",  50,  0.,    10. );

    hvec_bjet_pt.push_back(     h_bjet_pt     );
    hvec_bjet_eta.push_back(    h_bjet_eta    );
    hvec_bjet_phi.push_back(    h_bjet_phi    );
    hvec_bjet_emFrac.push_back( h_bjet_emFrac );
    hvec_bjet_bDiscr.push_back( h_bjet_bDiscr );
    hvec_bjetmuon_deltaR.push_back(     h_bjetmuon_deltaR     );
    hvec_bjetelectron_deltaR.push_back( h_bjetelectron_deltaR );
    hvec_bjetmet_deltaR.push_back(      h_bjetmet_deltaR      );

  }
  //  TFileDirectory subSubDir = subDir.mkdir( "mySubSubDirectory" );
  //  TH1D * h_pt = subDir.make<TH1D>( "pt"  , "p_{t}", 100,  0., 100. );

}

// ------------ method called once each job just after ending the event loop  ------------
void 
VBFHZZllbbRAWCORBTAGtesting::endJob() {
}

//define this as a plug-in
//DEFINE_FWK_MODULE(VBFHZZllbbRAWCORBTAGtesting);

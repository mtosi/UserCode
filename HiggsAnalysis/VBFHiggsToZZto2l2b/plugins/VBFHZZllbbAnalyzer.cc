// -*- C++ -*-
//
// Package:    VBFHZZllbbAnalyzer
// Class:      VBFHZZllbbAnalyzer
// 
/**\class VBFHZZllbbAnalyzer VBFHZZllbbAnalyzer.cc HiggsAnalysis/VBFHiggsToZZto2l2b/src/VBFHZZllbbAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Mia TOSI
//         Created:  Mon Feb  2 17:31:44 CET 2009
// $Id: VBFHZZllbbAnalyzer.cc,v 1.3 2009/07/29 09:57:29 tosi Exp $
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

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/RefToBase.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"


#include "PhysicsTools/UtilAlgos/interface/TFileDirectory.h"
#include "PhysicsTools/UtilAlgos/interface/TH1AddDirectorySentry.h" 

#include "PhysicsTools/Utilities/interface/deltaR.h"
#include "PhysicsTools/Utilities/interface/deltaPhi.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/JetReco/interface/GenJet.h"

#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"

// b-tagging
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"

// jet correction
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/CorJetWithBTagDiscr.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbAnalyzer.h"


// utilities
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/VBFHZZllbbUtils.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/ProcessIndex.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/PythiaParticleIndex.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/ParticlesMass.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/ParticlesCharge.h"

// For file output
// ---------------
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace edm;
using namespace reco;
using namespace vbfhzz2l2b;

// class decleration
//

//
// constants, enums and typedefs
//
const double PI_ = 3.141593;
//
// static data member definitions
//

//
// constructors and destructor
//
VBFHZZllbbAnalyzer::VBFHZZllbbAnalyzer(const edm::ParameterSet& iConfig) :
  whichSim_          ( iConfig.getParameter<int> ( "whichSim"  ) ), // 0:FastSim, 1:FullSim
  vertexLabel_       ( iConfig.getParameter<edm::InputTag> ( "vertexLabel"      ) ),
  trackLabel_        ( iConfig.getParameter<edm::InputTag> ( "trackLabel"       ) ),
  muonLabel_         ( iConfig.getParameter<edm::InputTag> ( "muonLabel"        ) ),
  electronLabel_     ( iConfig.getParameter<edm::InputTag> ( "electronLabel"    ) ),
  eleTrkIsoAlgoFlag_ ( iConfig.getParameter<bool>          ( "eleTrkIsoAlgoFlag") ),
  metLabel_          ( iConfig.getParameter<edm::InputTag> ( "metLabel"         ) ),
  caloJetLabel_      ( iConfig.getParameter<edm::InputTag> ( "caloJetLabel"     ) ),
  corIC5CaloJetsWithBTagLabel_ ( iConfig.getParameter<std::string> ( "corIC5CaloJetsWithBTagLabel" ) ),
  corIC5PFJetsWithBTagFlag_    ( iConfig.getParameter<bool>        ( "corIC5PFJetsWithBTagFlag"    ) ),
  genParticleLabel_  ( iConfig.getParameter<edm::InputTag> ( "genParticleLabel" ) ),
  genJetLabel_       ( iConfig.getParameter<edm::InputTag> ( "genJetLabel"      ) ),
  genMetLabel_       ( iConfig.getParameter<edm::InputTag> ( "genMetLabel"      ) ), 
  leptonPtCut_         ( iConfig.getParameter<double>("leptonPtCut"        ) ),
  jetEtCut_            ( iConfig.getParameter<double>("jetEtCut"           ) ),
  jetPartonDeltaR2Cut_ ( iConfig.getParameter<double>("jetPartonDeltaR2Cut") ),
  jetLeptonDeltaRCut_  ( iConfig.getParameter<double>("jetLeptonDeltaRCut" ) ),
  jetEMfracCut_        ( iConfig.getParameter<double>("jetEMfracCut"       ) )
{

  if ( corIC5PFJetsWithBTagFlag_ ) 
    corIC5PFJetsWithBTagLabel_ = iConfig.getParameter<std::string> ( "corIC5PFJetsWithBTagLabel" );
  

  if ( eleTrkIsoAlgoFlag_ )
    eleTrkIsoAlgo_ = new VBFHZZllbbElectronTrackIsolationAlgos(
	iConfig.getParameter         <double> ("coneRadius") ,
	iConfig.getParameter         <double> ("vetoRadius") ,
	iConfig.getParameter         <double> ("otherVetoRadius") ,
	iConfig.getParameter         <double> ("ptMin") ,
	iConfig.getParameter         <double> ("lipMax") ,
	iConfig.getUntrackedParameter<bool>   ("useTkQuality",true)
	);

  //now do what ever initialization is needed
  edm::Service<TFileService> fs ;
  //  mytree_  = fs->make <TTree>("VBFHZZhistos","VBFHZZhistos"); 


  // Now do what ever initialization is needed
  // -----------------------------------------

  //
  // constants, enums and typedefs
  nleptons_     = 2;
  njets_        = 3;
  nZjets_       = 2;
  nforwardjets_ = 2;

  eventcounter_             = 0;
  badeventcounter_          = 0;
  goodLeptonicEventCounter_ = 0;
  twoJetsAboveEtCutEventsCounter_ = 0;
  twoJetsAboveEtCut_DeltaR2016EventsCounter_ = 0;
  twoJetsAboveEtCut_DeltaR2CutEventsCounter_ = 0;
  fourJetsAboveEtCutEventsCounter_ = 0;

  nbin_ = 100;

  //  gROOT->Time();

}


VBFHZZllbbAnalyzer::~VBFHZZllbbAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  delete eleTrkIsoAlgo_;
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
VBFHZZllbbAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  eventcounter_++;
  if ( eventcounter_/1000 == float(eventcounter_)/100. ) {
    std::cout << "Event number " << eventcounter_ << std::endl;
  }

  // muons
  // -----
  edm::Handle<reco::MuonCollection> muonHandle;
  iEvent.getByLabel (muonLabel_,muonHandle);
  // electrons
  // ---------
  edm::Handle<reco::PixelMatchGsfElectronCollection> electronHandle ;
  iEvent.getByLabel (electronLabel_,electronHandle) ;
  // tracks
  // ------
  edm::Handle<reco::TrackCollection> trackHandle ;
  iEvent.getByLabel (trackLabel_, trackHandle) ;
  // Calo jets
  // ---------
  edm::Handle<reco::CaloJetCollection> caloJetHandle;
  iEvent.getByLabel (caloJetLabel_, caloJetHandle) ;
  // corIC5calo jet w/ b tag info
  // ----------------------------
  edm::Handle<vbfhzz2l2b::CorJetWithBTagDiscrCollection> corJetWithBTagHandle;
  iEvent.getByLabel(corIC5CaloJetsWithBTagLabel_,"corJetWithBTagDiscr",corJetWithBTagHandle);
  // secondaryVertexTagInfos
  // -----------------------
  edm::Handle<reco::SecondaryVertexTagInfoCollection> secondaryVtxTagInfosHandle;
  iEvent.getByLabel( "secondaryVertexTagInfos", secondaryVtxTagInfosHandle );  
  // calo met
  // --------
  edm::Handle<reco::CaloMETCollection> metCollectionHandle;
  iEvent.getByLabel (metLabel_ , metCollectionHandle);
  // genparticles
  // ------------
  edm::Handle<reco::GenParticleCollection> genParticleHandle; 
  iEvent.getByLabel (genParticleLabel_,genParticleHandle);

  if ( ( muonHandle->size() < nleptons_ && electronHandle->size() < nleptons_ ) ||
       corJetWithBTagHandle->size() < njets_ ) {
    badeventcounter_++;
    throw cms::Exception("NotFound") 
      << "BAD event:\n"
      << "# of muons:     " << muonHandle->size()           << " <--> " << nleptons_ << "\n"
      << "# of electrons: " << electronHandle->size()       << " <--> " << nleptons_ << "\n"
      << "# of jets:      " << corJetWithBTagHandle->size() << " <--> " << njets_    << "\n";
  }


  /////////////////////////////////// QUICK LOOP over LEPTONS ////////////////////////////
  // muons
  // -----
  h_muon_N->Fill(muonHandle->size());
  unsigned goodmuoncounter = 0;
  std::vector<reco::Muon> goodMuonVec;
  bool goodMuonPair = kFALSE;    
  for ( reco::MuonCollection::const_iterator muon_itr = muonHandle->begin();
	muon_itr != muonHandle->end(); ++muon_itr ) {
    h_muon_pt ->Fill( muon_itr->pt()  );
    h_muon_eta->Fill( muon_itr->eta() );
    h_muon_phi->Fill( muon_itr->phi() );

    if ( muon_itr->pt() < leptonPtCut_ ) continue;
    goodmuoncounter++;
    if (!goodMuonPair) 
      if ( muon_itr != (muonHandle->end()-1) )
	for ( reco::MuonCollection::const_iterator muon_itr2 = muon_itr+1;
	      muon_itr2 != muonHandle->end(); ++muon_itr2 ){
	  if ( muon_itr->charge()*muon_itr2->charge() < 0 ) {
	    goodMuonVec.push_back(*muon_itr);
	    goodMuonVec.push_back(*muon_itr2);
	    goodMuonPair = kTRUE;
	  } // if opposite charged muons pair
	} // loop over muon_itr2
  } // loop over muon_itr
  //  std::cout << "goodmuoncounter: " << goodmuoncounter 
  //	    << " <---> goodMuonVec: "  << goodMuonVec.size() << std::endl;

  if (goodMuonVec.size() != 0) {
    //    std::cout << "dimuonMass: " << (goodMuonVec[0].p4()+goodMuonVec[1].p4()).M() << std::endl;
    h_dimuon_mass -> Fill ((goodMuonVec[0].p4()+goodMuonVec[1].p4()).M()   );
    h_dimuon_pt   -> Fill ((goodMuonVec[0].p4()+goodMuonVec[1].p4()).pt()  );
    h_dimuon_eta  -> Fill ((goodMuonVec[0].p4()+goodMuonVec[1].p4()).eta() );
    h_dimuon_phi  -> Fill ((goodMuonVec[0].p4()+goodMuonVec[1].p4()).phi() );
    h_dimuon_deltaEta -> Fill (fabs(goodMuonVec[0].eta()-goodMuonVec[1].eta()));				 
    h_dimuon_deltaR   -> Fill (DeltaR  <reco::Muon, reco::Muon>()(goodMuonVec[0],goodMuonVec[1]));
    h_dimuon_deltaPhi -> Fill (DeltaPhi<reco::Muon, reco::Muon>()(goodMuonVec[0],goodMuonVec[1]));
    //PI_ - std::fabs( std::fabs( goodMuonVec[0].phi() - goodMuonVec[1].phi() ) - PI_ ) );
    h_dimuon_pz1pz2   -> Fill (goodMuonVec[0].pz()*goodMuonVec[1].pz());					 
    h_dimuon_eta1eta2 -> Fill (goodMuonVec[0].eta()*goodMuonVec[1].eta());                                       
    h_dimuon_massVSpt       -> Fill((goodMuonVec[0].p4()+goodMuonVec[1].p4()).M(),(goodMuonVec[0].p4()+goodMuonVec[1].p4()).pt() );
    h_dimuon_massVSeta      -> Fill((goodMuonVec[0].p4()+goodMuonVec[1].p4()).M(),(goodMuonVec[0].p4()+goodMuonVec[1].p4()).eta());
    h_dimuon_massVSphi      -> Fill((goodMuonVec[0].p4()+goodMuonVec[1].p4()).M(),(goodMuonVec[0].p4()+goodMuonVec[1].p4()).phi());
    h_dimuon_massVSdeltaEta -> Fill((goodMuonVec[0].p4()+goodMuonVec[1].p4()).M(),fabs(goodMuonVec[0].eta()-goodMuonVec[1].eta()) 				    );
    h_dimuon_massVSdeltaR   -> Fill((goodMuonVec[0].p4()+goodMuonVec[1].p4()).M(),DeltaR<reco::Muon, reco::Muon>()(goodMuonVec[0],goodMuonVec[1])		    );
    h_dimuon_massVSdeltaPhi -> Fill((goodMuonVec[0].p4()+goodMuonVec[1].p4()).M(),DeltaPhi<reco::Muon, reco::Muon>()(goodMuonVec[0],goodMuonVec[1]));
    //PI_ - std::fabs( std::fabs( goodMuonVec[0].phi() - goodMuonVec[1].phi() ) - PI_ ) );
    h_dimuon_massVSpz1pz2   -> Fill((goodMuonVec[0].p4()+goodMuonVec[1].p4()).M(),goodMuonVec[0].pz()*goodMuonVec[1].pz()					    );
    h_dimuon_massVSeta1eta2 -> Fill((goodMuonVec[0].p4()+goodMuonVec[1].p4()).M(),goodMuonVec[0].eta()*goodMuonVec[1].eta()                                         );
  }
  // electrons
  // ---------   
  h_electron_N->Fill(electronHandle->size());
  std::vector<reco::PixelMatchGsfElectron> goodElectronVec; 
  bool goodElectronPair = kFALSE;    
  for ( reco::PixelMatchGsfElectronCollection::const_iterator electron_itr = electronHandle->begin();
	electron_itr != electronHandle->end(); ++electron_itr ) {
    h_electron_pt-> Fill( electron_itr->pt()  );
    h_electron_eta->Fill( electron_itr->eta() );
    h_electron_phi->Fill( electron_itr->phi() );
    if ( electron_itr->pt() < leptonPtCut_ ) continue;
    if (!goodElectronPair) 
      if ( electron_itr != (electronHandle->end()-1) )
	for ( reco::PixelMatchGsfElectronCollection::const_iterator electron_itr2 = electron_itr+1;
	      electron_itr2 != electronHandle->end(); ++electron_itr2 ){
	  if ( electron_itr->charge()*electron_itr2->charge() < 0 ) {
	    goodElectronVec.push_back(*electron_itr);
	    goodElectronVec.push_back(*electron_itr2);
	    goodElectronPair = kTRUE;
	  } // if opposite charged electrons pair
	} // loop over electron_itr2
  } // loop over electron_itr
  //  std::cout << "goodElectronVec: "  << goodElectronVec.size() << std::endl;
  //  if (goodElectronVec.size() != 0)
  //    std::cout << "dielectronMass: " << (goodElectronVec[0].p4()+goodElectronVec[1].p4()).M() << std::endl;

  /////////////////////////////////// QUICK LOOP over JETS ////////////////////////////
  std::vector<vbfhzz2l2b::CorJetWithBTagDiscrCollection> goodJetVec;
  h_jet_N->Fill(corJetWithBTagHandle->size());
  bool goodJetPair = kFALSE;    
  unsigned bjetcounter = 0;
  std::vector<reco::JetBaseRef> jets = vbfhzz2l2b::CorJetBTagDiscrAssociation::allJets(*corJetWithBTagHandle);  
  for ( std::vector<reco::JetBaseRef>::const_iterator jet_itr = jets.begin();
	jet_itr != jets.end(); ++jet_itr ) {
    std::vector<double> discrVec = (*corJetWithBTagHandle)[*jet_itr].discrVec_;
    if ( discrVec[vbfhzz2l2b::COMBSECVTX] > 0.8 ) bjetcounter++;

    for ( reco::MuonCollection::const_iterator muon_itr = muonHandle->begin();
	  muon_itr != muonHandle->end(); ++muon_itr ) {
      double deltaR = DeltaR<reco::Muon, reco::Jet>()(*muon_itr, **jet_itr);
      h_jetmuon_deltaR->Fill(deltaR);
      if (discrVec[vbfhzz2l2b::COMBSECVTX] > 0.8) h_bjetmuon_deltaR->Fill(deltaR); 
    }
    for ( reco::PixelMatchGsfElectronCollection::const_iterator electron_itr = electronHandle->begin(); 
	  electron_itr != electronHandle->end(); ++electron_itr ) {
      double deltaR = DeltaR<reco::GsfElectron, reco::Jet>()(*electron_itr, **jet_itr);
      h_jetelectron_deltaR->Fill(deltaR);
      if (discrVec[vbfhzz2l2b::COMBSECVTX] > 0.8) h_bjetelectron_deltaR->Fill(deltaR); 
    }
    const reco::CaloMETCollection * metCollection = metCollectionHandle.product();
    for ( reco::CaloMETCollection::const_iterator met = metCollection->begin();
	  met != metCollection->end(); ++met ) {
      double deltaR = DeltaR<reco::CaloMET, reco::Jet>()(*met, **jet_itr);
      h_jetmet_deltaR->Fill(deltaR);
      if (discrVec[vbfhzz2l2b::COMBSECVTX] > 0.8) h_bjetmet_deltaR->Fill(deltaR); 
    }

    double corScale = (*corJetWithBTagHandle)[*jet_itr].corEt_/(*jet_itr)->et();
    h_jet_bDiscr->Fill( discrVec[vbfhzz2l2b::COMBSECVTX]        );
    h_jet_pt->    Fill( (*jet_itr)->pt()                        );
    h_jet_coret-> Fill( (*corJetWithBTagHandle)[*jet_itr].corEt_);
    h_jet_et->    Fill( (*jet_itr)->et()                        );
    h_jet_corpt-> Fill( (*jet_itr)->pt()*corScale               );
    h_jet_eta->   Fill( (*jet_itr)->eta()                       );
    h_jet_phi->   Fill( (*jet_itr)->phi()                       );
  }
  h_bjet_N->Fill(bjetcounter);

  const reco::CaloMETCollection * metCollection = metCollectionHandle.product();
  for ( reco::CaloMETCollection::const_iterator met = metCollection->begin();
	met != metCollection->end(); ++met ) {
    h_met_pt->Fill(  met->pt()  );
    h_met_eta->Fill( met->eta() );
    h_met_phi->Fill( met->phi() );
  }

  if ( goodMuonVec.size() < nleptons_ && goodElectronVec.size() < nleptons_ )
    throw cms::Exception("NotFound")
      << "goodMuonVec:     " << goodMuonVec.size()     << " <--> 2\n"
      << "goodElectronVec: " << goodElectronVec.size() << " <--> 2\n";
  // count number of events w/ at least 2 leptons
  // of the same flavour
  // opposite charge
  goodLeptonicEventCounter_++;




}


// ------------ method called once each job just before starting event loop  ------------
void 
VBFHZZllbbAnalyzer::beginJob(const edm::EventSetup&)
{
  TFileDirectory jetSubDir = fs->mkdir( "jet" );
  h_jet_N      = jetSubDir.make<TH1D>( "jet_N",      "number of jets",           nbin_,  0.  , 100.   );
  h_jet_et     = jetSubDir.make<TH1D>( "jet_et",     "jet E_{T}",                nbin_,  0.  , 200.   );
  h_jet_coret  = jetSubDir.make<TH1D>( "jet_coret",  "jet cor E_{T}",            nbin_,  0.  , 200.   );
  h_jet_pt     = jetSubDir.make<TH1D>( "jet_pt",     "jet p_{T}",                nbin_,  0.  , 200.   );
  h_jet_corpt  = jetSubDir.make<TH1D>( "jet_corpt",  "jet cor p_{T}",            nbin_,  0.  , 200.   );
  h_jet_eta    = jetSubDir.make<TH1D>( "jet_eta",    "jet #eta",                 nbin_,-10.  ,  10.   );
  h_jet_phi    = jetSubDir.make<TH1D>( "jet_phi",    "jet #phi",                 nbin_, -3.14,   3.14 );
  h_jet_bDiscr = jetSubDir.make<TH1D>( "jet_bDiscr", "jet b tag discriminator",  nbin_,-10.,     2.   );
  h_jetmuon_deltaR     = jetSubDir.make<TH1D>( "jetmuon_deltaR",     "#Delta R_{j,#mu}",nbin_,  0.,  10. );
  h_jetelectron_deltaR = jetSubDir.make<TH1D>( "jetelectron_deltaR", "#Delta R_{j,e}",  nbin_,  0.,  10. );
  h_jetmet_deltaR      = jetSubDir.make<TH1D>( "jetmet_deltaR",      "#Delta R_{j,met}",nbin_,  0.,  10. );
  
  TFileDirectory bjetSubDir = jetSubDir.mkdir( "bjet" );
  h_bjet_N              = bjetSubDir.make<TH1D>( "bjet_N",              "number of b-tagged jets", nbin_,  0.  , 100. );
  h_bjetmuon_deltaR     = bjetSubDir.make<TH1D>( "bjetmuon_deltaR",     "#Delta R_{j_{b},#mu}",    nbin_,  0.,    10. );
  h_bjetelectron_deltaR = bjetSubDir.make<TH1D>( "bjetelectron_deltaR", "#Delta R_{j_{b},e}",      nbin_,  0.,    10. );
  h_bjetmet_deltaR      = bjetSubDir.make<TH1D>( "bjetmet_deltaR",      "#Delta R_{j_{b},met}",    nbin_,  0.,    10. );

  TFileDirectory muonSubDir = fs->mkdir( "muon" );
  h_muon_N   = muonSubDir.make<TH1D>( "muon_N",   "number of muons", nbin_,  0.  , 100.   );
  h_muon_pt  = muonSubDir.make<TH1D>( "muon_pt",  "muon p_{T}",      nbin_,  0.  , 200.   );
  h_muon_eta = muonSubDir.make<TH1D>( "muon_eta", "muon #eta",       nbin_,-10.  ,  10.   );
  h_muon_phi = muonSubDir.make<TH1D>( "muon_phi", "muon #phi",       nbin_, -3.14,   3.14 );

  TFileDirectory dimuonSubDir = muonSubDir.mkdir( "dimuon" );
  h_dimuon_mass = dimuonSubDir.make<TH1D>( "dimuon_mass","di-muon invariant mass",nbin_,  0.  ,100.  );
  h_dimuon_pt   = dimuonSubDir.make<TH1D>( "dimuon_pt",  "di-muon p_{T}",         nbin_,  0.  ,100.  );
  h_dimuon_eta  = dimuonSubDir.make<TH1D>( "dimuon_eta", "di-muon #eta",          nbin_,-10.  , 10.  );
  h_dimuon_phi  = dimuonSubDir.make<TH1D>( "dimuon_phi", "di-muon #phi",          nbin_, -3.14,  3.14);
  h_dimuon_deltaEta = dimuonSubDir.make<TH1D>( "dimuon_deltaEta","di-muon #Delta#eta",          nbin_,     0.,   10.);
  h_dimuon_deltaR   = dimuonSubDir.make<TH1D>( "dimuon_deltaR",  "di-muon #DeltaR",             nbin_,     0.,   10.);
  h_dimuon_deltaPhi = dimuonSubDir.make<TH1D>( "dimuon_deltaPhi","di-muon #Delta#phi",          nbin_,     0.,   10.);
  h_dimuon_pz1pz2   = dimuonSubDir.make<TH1D>( "dimuon_pz1pz2",  "di-muon pz_{1}*pz_{2}",    20*nbin_, -4000., 4000.);
  h_dimuon_eta1eta2 = dimuonSubDir.make<TH1D>( "dimuon_eta1eta2","di-muon #eta_{1}*#eta_{2}",   nbin_,    -9.,    9.);
  h_dimuon_massVSpt       = dimuonSubDir.make<TH2D>( "dimuon_massVSpt",      "di-muon inv mass VS p_{T}",            nbin_,0.,100.,   nbin_,     0.  ,  100.  );
  h_dimuon_massVSeta      = dimuonSubDir.make<TH2D>( "dimuon_massVSeta",     "di-muon inv mass VS #eta",             nbin_,0.,100.,   nbin_,   -10.  ,   10.  );
  h_dimuon_massVSphi      = dimuonSubDir.make<TH2D>( "dimuon_massVSphi",     "di-muon inv mass VS #phi",             nbin_,0.,100.,   nbin_,    -3.14,    3.14);
  h_dimuon_massVSdeltaEta = dimuonSubDir.make<TH2D>( "dimuon_massVSdeltaEta","di-muon inv mass VS #Delta#eta",       nbin_,0.,100.,   nbin_,     0.  ,   10.  );
  h_dimuon_massVSdeltaR   = dimuonSubDir.make<TH2D>( "dimuon_massVSdeltaR",  "di-muon inv mass VS #DeltaR",          nbin_,0.,100.,   nbin_,     0.  ,   10.  );
  h_dimuon_massVSdeltaPhi = dimuonSubDir.make<TH2D>( "dimuon_massVSdeltaPhi","di-muon inv mass VS #Delta#phi",       nbin_,0.,100.,   nbin_,     0.  ,   10.  );
  h_dimuon_massVSpz1pz2   = dimuonSubDir.make<TH2D>( "dimuon_massVSpz1pz2",  "di-muon inv mass VS pz_{1}*pz_{2}",    nbin_,0.,100.,20*nbin_, -4000.  , 4000.  );
  h_dimuon_massVSeta1eta2 = dimuonSubDir.make<TH2D>( "dimuon_massVSeta1eta2","di-muon inv mass VS #eta_{1}*#eta_{2}",nbin_,0.,100.,   nbin_,    -9.  ,    9.  );

  TFileDirectory electronSubDir = fs->mkdir( "electron" );
  h_electron_N   = electronSubDir.make<TH1D>( "electron_N",   "number of electrons", nbin_,  0.  , 100.   );
  h_electron_pt  = electronSubDir.make<TH1D>( "electron_pt",  "electron p_{T}",      nbin_,  0.  , 200.   );
  h_electron_eta = electronSubDir.make<TH1D>( "electron_eta", "electron #eta",       nbin_,-10.  ,  10.   );
  h_electron_phi = electronSubDir.make<TH1D>( "electron_phi", "electron #phi",       nbin_, -3.14,   3.14 );
  
  TFileDirectory metSubDir = fs->mkdir( "met" );
  h_met_N   = metSubDir.make<TH1D>( "met_N",   "number of mets", nbin_,  0.  , 100.   );
  h_met_pt  = metSubDir.make<TH1D>( "met_pt",  "met p_{T}",      nbin_,  0.  , 200.   );
  h_met_eta = metSubDir.make<TH1D>( "met_eta", "met #eta",       nbin_,-10.  ,  10.   );
  h_met_phi = metSubDir.make<TH1D>( "met_phi", "met #phi",       nbin_, -3.14,   3.14 );
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
VBFHZZllbbAnalyzer::endJob() {
}


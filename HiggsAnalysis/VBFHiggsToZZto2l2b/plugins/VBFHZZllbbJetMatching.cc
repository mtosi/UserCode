//
// Original Author:  Mia Tosi
//         Created:  Fri Feb 22 17:56:22 CET 2008
// $Id: VBFHZZllbbJetMatching.cc,v 1.2 2009/04/28 16:45:01 tosi Exp $
//
//


// system include files
//#include <memory>

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbJetMatching.h"

// Classes to be accessed
// ----------------------
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

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

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/VBFHZZllbbUtils.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/PythiaParticleIndex.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/ParticlesMass.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/ParticlesCharge.h"

#include "PhysicsTools/Utilities/interface/deltaR.h"
#include "PhysicsTools/Utilities/interface/deltaPhi.h"
//#include "DataFormats/Math/interface/deltaR.h"
//#include "DataFormats/Math/interface/deltaPhi.h"

// For file output
// ---------------
#include <fstream>
#include <sstream>
#include <cmath>
#include <memory>
#include <string>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace edm;
using namespace reco;
using namespace vbfhzz2l2b;


enum { FASTSIM = 0,
       FULLSIM = 1
};
enum { NZ        = 2,
       NZLEPTONS = 2,
       NZJETS    = 2,
       NTAGJETS  = 2,
       NJETS     = 4,
       NBIN      = 100
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

// constructors
VBFHZZllbbJetMatching::VBFHZZllbbJetMatching(const edm::ParameterSet& iConfig) :
  electronLabel_        ( iConfig.getParameter<edm::InputTag> ( "electronLabel"        ) ),
  muonLabel_            ( iConfig.getParameter<edm::InputTag> ( "muonLabel"            ) ),
  metLabel_             ( iConfig.getParameter<edm::InputTag> ( "metLabel"             ) ),
  jetLabel_             ( iConfig.getParameter<edm::InputTag> ( "jetLabel"             ) ),
  corJetsWithBTagLabel_ ( iConfig.getParameter<std::string>   ( "corJetsWithBTagLabel" ) ),
  mcParticleLabel_      ( iConfig.getParameter<edm::InputTag> ( "mcParticleLabel"      ) ),
  genJetLabel_          ( iConfig.getParameter<edm::InputTag> ( "genJetLabel"          ) ),
  genMetLabel_          ( iConfig.getParameter<edm::InputTag> ( "genMetLabel"          ) )
{
  // Now do what ever initialization is needed
  // -----------------------------------------

  //
  // constants, enums and typedefs
  eventcounter_ = 0;

  std::cout << "[VBFHZZllbbJetMatching::VBFHZZllbbJetMatching]" << std::endl;

  //  gROOT->Time();

}


// destructor
VBFHZZllbbJetMatching::~VBFHZZllbbJetMatching()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  std::cout << "[VBFHZZllbbJetMatching::~VBFHZZllbbJetMatching]" << std::endl;

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
VBFHZZllbbJetMatching::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::cout << "[VBFHZZllbbJetMatching::analyze]" << std::endl;

  eventcounter_++;
  if ( eventcounter_/100 == float(eventcounter_)/100. ) {
    std::cout << "Event number " << eventcounter_ << std::endl;
  }

  // muons
  // -----
  edm::Handle<reco::MuonCollection> muonHandle;
  iEvent.getByLabel (muonLabel_,muonHandle);
  // electrons
  // ---------
  edm::Handle<edm::View<reco::PixelMatchGsfElectron> > electronHandle ;
  iEvent.getByLabel (electronLabel_,electronHandle) ;
  // calo jets
  // ---------
  edm::Handle<reco::CaloJetCollection> jetHandle ;
  iEvent.getByLabel (jetLabel_, jetHandle) ;
  // corrected jet w/ b tagging discriminator
  // ----------------------------------------
  edm::Handle<vbfhzz2l2b::CorJetWithBTagDiscrCollection> corJetWithBTagHandle;
  iEvent.getByLabel(corJetsWithBTagLabel_,"corJetWithBTagDiscr",corJetWithBTagHandle);
  // MCParticle
  // ----------
  edm::Handle<reco::GenParticleCollection> mcParticleHandle; 
  iEvent.getByLabel (mcParticleLabel_,mcParticleHandle);
  // met
  // ---
  edm::Handle<reco::CaloMETCollection> metHandle;
  iEvent.getByLabel (metLabel_ , metHandle);
  const CaloMETCollection *caloMETcoll = metHandle.product();
  const CaloMET *caloMET = &(caloMETcoll->front());


  edm::Handle<reco::VertexCollection> primaryVerticesWithBSHandle;
  iEvent.getByLabel("offlinePrimaryVerticesWithBS",primaryVerticesWithBSHandle);
  edm::Handle<reco::VertexCollection> primaryVerticesHandle;
  iEvent.getByLabel("offlinePrimaryVertices",primaryVerticesHandle);
  edm::Handle<reco::VertexCollection> pixelVerticesHandle;
  iEvent.getByLabel("pixelVertices",pixelVerticesHandle);

  std::cout << "primaryVerticesWithBSHandle->size(): " << primaryVerticesWithBSHandle->size() << std::endl;
  std::cout << "primaryVerticesHandle->size(): " << primaryVerticesHandle->size() << std::endl;
  std::cout << "pixelVerticesHandle->size(): " << pixelVerticesHandle->size() << std::endl;
  

  unsigned muonCollectionSize       = muonHandle->size();
  unsigned electronCollectionSize   = electronHandle->size();
  unsigned caloJetCollectionSize    = jetHandle->size();
  unsigned metCollectionSize        = metHandle->size();
  unsigned mcParticleCollectionSize = mcParticleHandle->size();

  std::cout << "mcParticleCollectionSize: " << mcParticleCollectionSize << std::endl;
  
  /////////////////////////////////// HEPG analysis ////////////////////////////
  // Z counters 
  unsigned int Zparticles = 0;
  unsigned int hadronicZparticles = 0;
  unsigned int leptonicZparticles = 0;
  unsigned int ZintoEcounter   = 0; 
  unsigned int ZintoMUcounter  = 0;
  unsigned int ZintoTAUcounter = 0;
  unsigned int ZintoBQUARKcounter = 0;

  std::vector< const Candidate* > hadronicZ;
  std::vector< const Candidate* > leptonicZ;
  std::vector< const Candidate* > H;
  std::vector< const Candidate* > tagSystem;
  
  const Candidate * Hmother1 = 0;
  const Candidate * Hmother2 = 0;
  const Candidate * Zmother1 = 0;
  const Candidate * Zmother2 = 0;

  GenParticleCollection::const_iterator MCparticle = mcParticleHandle->begin();
  MCparticle = mcParticleHandle->begin();
  for ( ; MCparticle != mcParticleHandle->end(); ++MCparticle ) {    
    int tmpParticleId = MCparticle->pdgId();
    if ( fabs(tmpParticleId) > pythiaH_ ) continue;
    int tmpMotherId = 0;
    if( MCparticle->mother() != 0 ) tmpMotherId = MCparticle->mother()->pdgId();
    
    // H boson
    if ( tmpParticleId == pythiaH_ && tmpParticleId != tmpMotherId ) {
      int tmpNumberOfMothers = MCparticle->numberOfMothers();
      if ( tmpNumberOfMothers != 2 ) {
	 std::cout << "--> Higgs number of mother " << tmpNumberOfMothers << "(";
	 for ( int motherIndex = 0; motherIndex != tmpNumberOfMothers; motherIndex++ )
	   std::cout << "ID: " << MCparticle->mother(motherIndex)->pdgId();
	 std::cout << ")" << std::endl;
      }
      else {
	Hmother1 = MCparticle->mother(0);
	Hmother2 = MCparticle->mother(1);
	//	std::cout << "Hmother1: " << Hmother1->pdgId() << std::endl;
	//	std::cout << "Hmother2: " << Hmother2->pdgId() << std::endl;
      }
    }
    // H decay
    if ( tmpMotherId == pythiaH_ && tmpParticleId != tmpMotherId ) {
      std::cout << "status: " << MCparticle->status() << std::endl;
      if ( tmpParticleId != pythiaZ_ ) 
	std::cout << "WARNING!! --> Higgs decaying into " << tmpParticleId << std::endl;
    }
    
    if (tagSystem.size() == 0) {
      int tmpNumberOfDaughters = MCparticle->numberOfDaughters();
      if ( fabs(tmpParticleId) <= pythiat_ && tmpNumberOfDaughters == 2 ) {    
	for (int daughterIndex = 0; daughterIndex != tmpNumberOfDaughters; daughterIndex++ ) {
	  std::cout << "sister " << daughterIndex << " ID: " << MCparticle->daughter(daughterIndex)->pdgId() << std::endl;
	  if ( fabs(MCparticle->daughter(daughterIndex)->pdgId()) <= pythiat_ )
	    tagSystem.push_back(MCparticle->daughter(daughterIndex));
	}
      }
      if (tagSystem.size() != 0) {
	std::cout << "tagSystem BEFORE: " << tagSystem.size() << " ";
	if (tagSystem.size() != 2) tagSystem.clear();
	else 
	  std::cout << "ID1: " << tagSystem[0]->pdgId() << " ID2: " << tagSystem[1]->pdgId() << std::endl;
	std::cout << " --> AFTER: " << tagSystem.size() << std::endl;
      }
    }
    
    // Z boson
    if ( tmpParticleId == pythiaZ_ && tmpMotherId != tmpParticleId ) {
      std::cout << "status: " << MCparticle->status() << std::endl;
      int tmpNumberOfMothers = MCparticle->numberOfMothers();
      if ( tmpNumberOfMothers == 1 ) {
	if ( tmpMotherId != pythiaH_ ) {
	  std::cout << "found " << Zparticles+1 << " Z not from H decay, from " << tmpMotherId << std::endl;
	  if (Zparticles == 0) Zmother1 = MCparticle->mother(0);
	  if (Zparticles == 1) Zmother2 = MCparticle->mother(0);
	}
      } else if ( tmpNumberOfMothers == 2 ) {
	std::cout << "found Z not from H decay, from " << MCparticle->mother(0)->pdgId() 
		  << " and " << MCparticle->mother(1)->pdgId() << std::endl;
	Zmother1 = MCparticle->mother(0);
	Zmother2 = MCparticle->mother(1);
      }
      Zparticles++; 

      int tmpNumberOfDaughters = MCparticle->numberOfDaughters();
      //      std::cout << "Z numberOfDaughters: " << tmpNumberOfDaughters << std::endl;
      for ( int daughterIndex = 0; daughterIndex != tmpNumberOfDaughters; daughterIndex++ ) {
	int tmpZdaughterId = MCparticle->daughter(daughterIndex)->pdgId();
	// hadronically decay Z
	if ( fabs(tmpZdaughterId) <= pythiat_ ) {
	  //	  std::cout << "found hadronically decaying Z into " << tmpZdaughterId << std::endl;
	  hadronicZparticles++;
	  hadronicZ.push_back(&*(MCparticle->daughter(daughterIndex)));
	  if ( fabs(tmpZdaughterId) == pythiab_ ) ZintoBQUARKcounter++;
	}
	// leptonically decay Z
	else if ( fabs(tmpZdaughterId) >= pythiae_ && fabs(tmpZdaughterId) <= pythiatau_ ) {
	  //	  std::cout << "found leptonically decaying Z into " << tmpZdaughterId << std::endl;
	  leptonicZparticles++;
	  if ( fabs(tmpZdaughterId) == pythiae_   ) ZintoEcounter++; 
	  if ( fabs(tmpZdaughterId) == pythiamu_  ) ZintoMUcounter++;
	  if ( fabs(tmpZdaughterId) == pythiatau_ ) ZintoTAUcounter++;
	  leptonicZ.push_back(&*(MCparticle->daughter(daughterIndex)));
	}
      }
    }
   } // end loop over MCparticle
   
  unsigned int heavyQcounter = 0; 
  unsigned int lightQcounter = 0;
  MCparticle = mcParticleHandle->begin();
  for ( ; MCparticle != mcParticleHandle->end(); ++MCparticle ) {    
    int tmpParticleId = MCparticle->pdgId();
    if ( fabs(tmpParticleId) > pythiat_ && tmpParticleId != pythiagluon_ ) continue;
    
    if ( Hmother1 != 0 || Hmother2 != 0 ) {
      if ( MCparticle->mother() == Hmother1 || MCparticle->mother() == Hmother2 ) {
	std::cout << "tag quark: " << tmpParticleId << std::endl;
	tagSystem.push_back(&*MCparticle);
      }
    }
    
    if ( Zmother1 != 0 || Zmother2 != 0 ) {
      if ( MCparticle->mother() == Zmother1 || MCparticle->mother() == Zmother2 ) {
	std::cout << "VQQ quark: " << tmpParticleId << std::endl;
	tagSystem.push_back(&*MCparticle);
      }
    }
  } // end loop over MCparticle
    
  if (tagSystem.size() == 2 ) {
    std::cout << "***************************************" << std::endl;
    std::cout << "pythiac: " << pythiac_ << std::endl;
    std::cout << "TAG SYSTEM: ";
    for ( unsigned int tagIndex = 0; tagIndex != tagSystem.size(); tagIndex++ ) {
      int tmpParticleId = tagSystem[tagIndex]->pdgId();
      std::cout << " " << " ID" << tagIndex+1 << ": " << tmpParticleId;
      if (tmpParticleId == pythiagluon_) ;
      else {
	if ( fabs(tmpParticleId) >= pythiac_ ) heavyQcounter++;
	else lightQcounter++;
      }
    }
    std::cout << "" << std::endl;
    std::cout << "***************************************" << std::endl;
  }

  std::cout << "heavyQcounter: " << heavyQcounter << std::endl;
  std::cout << "lightQcounter: " << lightQcounter << std::endl;

  ////////////////////////////////////// PLOT HEPG DISTRIBUTION ////////////////////////////
  std::vector< const Candidate* > partonsCandidates;
  if ( hadronicZ.size() == NZJETS ) {
    for ( std::vector< const Candidate* >::const_iterator hadronicZ_itr = hadronicZ.begin();
	  hadronicZ_itr != hadronicZ.end(); ++hadronicZ_itr) {
      double tmp_ZpartonEta = (*hadronicZ_itr)->eta();
      double tmp_ZpartonPt  = (*hadronicZ_itr)->pt();
      double tmp_ZpartonEt  = (*hadronicZ_itr)->et();
      double tmp_ZpartonE   = (*hadronicZ_itr)->energy();
      ZpartonEta_->Fill( tmp_ZpartonEta );
      ZpartonPt_ ->Fill( tmp_ZpartonPt  );
      ZpartonEt_ ->Fill( tmp_ZpartonEt  );
      ZpartonE_  ->Fill( tmp_ZpartonE   );

      partonsCandidates.push_back(*hadronicZ_itr);
    }
    double tmp_ZpartonsDeltaEta     = fabs(hadronicZ[0]->eta()-hadronicZ[1]->eta());
    double tmp_ZpartonsDeltaPhi     = deltaPhi(hadronicZ[0]->phi(),hadronicZ[1]->phi());
    double tmp_ZpartonsDeltaR       = deltaR(hadronicZ[0]->eta(),hadronicZ[0]->phi(),hadronicZ[1]->eta(),hadronicZ[1]->phi());
    double tmp_ZpartonsMass         = (hadronicZ[0]->p4()+(hadronicZ[1]->p4())).mass();
    double tmp_ZpartonsPt           = (hadronicZ[0]->p4()+(hadronicZ[1]->p4())).pt();
    //    double tmp_ZpartonsEt           = (hadronicZ[0]->p4()+(hadronicZ[1]->p4())).et();
    double tmp_ZpartonsEt           = 0.;
    double tmp_ZpartonsE            = (hadronicZ[0]->p4()+(hadronicZ[1]->p4())).e();
    double tmp_ZpartonsCollinearity = hadronicZ[0]->pz()*hadronicZ[1]->pz();
    ZpartonsDeltaEta_    ->Fill( tmp_ZpartonsDeltaEta     );
    ZpartonsDeltaR_      ->Fill( tmp_ZpartonsDeltaR       );
    ZpartonsMass_        ->Fill( tmp_ZpartonsMass         );
    ZpartonsPt_          ->Fill( tmp_ZpartonsPt           );
    ZpartonsEt_          ->Fill( tmp_ZpartonsEt           );
    ZpartonsE_           ->Fill( tmp_ZpartonsE            );
    ZpartonsCollinearity_->Fill( tmp_ZpartonsCollinearity );

  } 
  else
    std::cout << "hadronicZ.size: " << hadronicZ.size() << std::endl;

  if ( tagSystem.size() == NTAGJETS ) {
    std::vector< const Candidate* >::const_iterator tagSystem_itr = tagSystem.begin();
    for ( ; tagSystem_itr != tagSystem.end(); ++tagSystem_itr) {
      double tmp_TAGpartonEta = (*tagSystem_itr)->eta();
      double tmp_TAGpartonPt  = (*tagSystem_itr)->pt();
      double tmp_TAGpartonEt  = (*tagSystem_itr)->et();
      double tmp_TAGpartonE   = (*tagSystem_itr)->energy();
      TAGpartonEta_->Fill( tmp_TAGpartonEta );
      TAGpartonPt_ ->Fill( tmp_TAGpartonPt  );
      TAGpartonEt_ ->Fill( tmp_TAGpartonEt  );
      TAGpartonE_  ->Fill( tmp_TAGpartonE   );

      partonsCandidates.push_back(*tagSystem_itr);
    }
    double tmp_TAGpartonsDeltaEta     = fabs(tagSystem[0]->eta()-tagSystem[1]->eta());
    double tmp_TAGpartonsDeltaPhi     = deltaPhi(tagSystem[0]->phi(),tagSystem[1]->phi());
    double tmp_TAGpartonsDeltaR       = deltaR(tagSystem[0]->eta(),tagSystem[0]->phi(),tagSystem[1]->eta(),tagSystem[1]->phi());
    double tmp_TAGpartonsMass         = (tagSystem[0]->p4()+(tagSystem[1]->p4())).mass();
    double tmp_TAGpartonsPt           = (tagSystem[0]->p4()+(tagSystem[1]->p4())).pt();
    //    double tmp_TAGpartonsEt           = (tagSystem[0]->p4()+(tagSystem[1]->p4())).et();
    double tmp_TAGpartonsEt           = 0.;
    double tmp_TAGpartonsE            = (tagSystem[0]->p4()+(tagSystem[1]->p4())).e();
    double tmp_TAGpartonsCollinearity = tagSystem[0]->pz()*tagSystem[1]->pz();
    TAGpartonsDeltaEta_    ->Fill( tmp_TAGpartonsDeltaEta     );
    TAGpartonsDeltaR_      ->Fill( tmp_TAGpartonsDeltaR       );
    TAGpartonsMass_        ->Fill( tmp_TAGpartonsMass         );
    TAGpartonsPt_          ->Fill( tmp_TAGpartonsPt           );
    TAGpartonsEt_          ->Fill( tmp_TAGpartonsEt           );
    TAGpartonsE_           ->Fill( tmp_TAGpartonsE            );
    TAGpartonsCollinearity_->Fill( tmp_TAGpartonsCollinearity );
  }
  else
    std::cout << "tagSystem.size: " << tagSystem.size() << std::endl;
  
  //////////////////////////////////////////////////////////////////////////////////////////
  
  // corrected jet w/ b tagging discriminator
  // ----------------------------------------
  int jetIndex = 0;    
  std::vector<reco::JetBaseRef> jets = vbfhzz2l2b::CorJetBTagDiscrAssociation::allJets(*corJetWithBTagHandle);
  jetNumber_->Fill(jets.size());
  for ( std::vector<reco::JetBaseRef>::const_iterator jet_itr = jets.begin(); 
	jet_itr != jets.end();
	++jet_itr, jetIndex++ ) {
    std::cout << "jetIndex: " << jetIndex << std::endl;
    
    std::vector<double> discrVec = (*corJetWithBTagHandle)[*jet_itr].discrVec_;
    double tmp_jetCorrPt   = (*corJetWithBTagHandle)[*jet_itr].corPt_;
    double tmp_jetCorrEt   = (*corJetWithBTagHandle)[*jet_itr].corEt_;
    double tmp_jetUncorrPt = (*jet_itr)->pt();
    double tmp_jetUncorrEt = (*jet_itr)->et();
    double tmp_jetEta      = (*jet_itr)->eta();
    double tmp_jetPhi      = (*jet_itr)->phi();

    int tagger = -1;
    tagger = bTaggerCode("HIGHEFF"); double tmp_jetHIGHEFFdiscr = discrVec[tagger];
    std::cout << "HIGHEFF:    " << (*corJetWithBTagHandle)[*jet_itr].highEffDiscr_ 
	      << " <--> discrVec[" << tagger << "]: " << discrVec[tagger] << std::endl;
    tagger = bTaggerCode("HIGHPUR"); double tmp_jetHIGHPURdiscr = discrVec[tagger];	  
    std::cout << "HIGHPUR:    " << (*corJetWithBTagHandle)[*jet_itr].highPurDiscr_ 
	      << " <--> discrVec[" << tagger << "]: " << discrVec[tagger] << std::endl;
    tagger = bTaggerCode("COMBSECVTX"); double tmp_jetCOMBSECVTXdiscr = discrVec[tagger];
    std::cout << "COMBSECVTX: " << (*corJetWithBTagHandle)[*jet_itr].compoSVDiscr_
	      << " <--> discrVec[" << tagger << "]: " << discrVec[tagger] << std::endl;
    tagger = bTaggerCode("JETPROB"); double tmp_jetJETPROBdiscr = discrVec[tagger];   
    std::cout << "JETPROB:    " << (*corJetWithBTagHandle)[*jet_itr].jetProbDiscr_
	      << " <--> discrVec[" << tagger << "]: " << discrVec[tagger] << std::endl;

    double tmp_jetEMfrac = (dynamic_cast<const reco::CaloJet*>(&**jet_itr))->emEnergyFraction();
    std::cout << "emFrac: " << tmp_jetEMfrac << std::endl;

    jetUncorrEt_        ->Fill( tmp_jetUncorrEt        );
    jetCorrEt_          ->Fill( tmp_jetCorrEt          );
    jetUncorrPt_        ->Fill( tmp_jetUncorrPt        );
    jetCorrPt_          ->Fill( tmp_jetCorrPt          );
    jetPhi_             ->Fill( tmp_jetPhi             );
    jetEta_             ->Fill( tmp_jetEta             );
    jetEMfrac_          ->Fill( tmp_jetEMfrac          );
    jetHIGHEFFdiscr_    ->Fill( tmp_jetHIGHEFFdiscr    );	  
    jetHIGHPURdiscr_    ->Fill( tmp_jetHIGHPURdiscr    );	  
    jetCOMBSECVTXdiscr_ ->Fill( tmp_jetCOMBSECVTXdiscr );
    jetJETPROBdiscr_    ->Fill( tmp_jetJETPROBdiscr    );   

    jetEMfracVSeta_             ->Fill( tmp_jetEMfrac , tmp_jetEta             );	      
    jetEMfracVScorrEt_          ->Fill( tmp_jetEMfrac , tmp_jetCorrEt          );	      
    jetEMfracVScorrPt_          ->Fill( tmp_jetEMfrac , tmp_jetCorrPt          );	      
    jetEMfracVSuncorrEt_        ->Fill( tmp_jetEMfrac , tmp_jetUncorrEt        );	      
    jetEMfracVSuncorrPt_        ->Fill( tmp_jetEMfrac , tmp_jetUncorrPt        );	      
    jetEMfracVShighEFFdiscr_    ->Fill( tmp_jetEMfrac , tmp_jetHIGHEFFdiscr    );    
    jetEMfracVShighPURdiscr_    ->Fill( tmp_jetEMfrac , tmp_jetHIGHPURdiscr    );    
    jetEMfracVScomboSECVTXdiscr_->Fill( tmp_jetEMfrac , tmp_jetCOMBSECVTXdiscr );
    jetEMfracVSjetPROBdiscr_    ->Fill( tmp_jetEMfrac , tmp_jetJETPROBdiscr    );    

  }

    std::vector< reco::JetBaseRef > ZjetVec;
    std::vector< double >     ZjetDeltaRVec;
    std::vector<int>          ZjetInd;
    for ( unsigned int index = 0; index < NZJETS; index++ ) ZjetDeltaRVec.push_back(99.);

    std::vector< reco::JetBaseRef > TjetVec;
    std::vector< double >   TjetDeltaRVec;
    std::vector<int>        TjetInd;
    for ( unsigned int index = 0; index < NTAGJETS; index++) TjetDeltaRVec.push_back(99.);
    
    std::vector<std::vector<double> > closestJetDeltaRVec;
    std::vector<std::vector<int> >    closestJetIndexVec;
    unsigned int partonsCandidatesNumber = partonsCandidates.size();
    for ( unsigned index = 0; index < partonsCandidatesNumber; index++ ) {
      std::vector<double> null_closestJetDeltaRVec;
      std::vector<int>    null_closestJetIndexVec;
      for ( unsigned index = 0; index < partonsCandidatesNumber; index++ ) {
	null_closestJetDeltaRVec.push_back(99.);
	null_closestJetIndexVec.push_back(-1);
      }
      closestJetDeltaRVec.push_back(null_closestJetDeltaRVec);
      closestJetIndexVec.push_back(null_closestJetIndexVec);
    }
    
    std::vector< const Candidate* >::const_iterator partonsCandidates_itr = partonsCandidates.begin();
    unsigned int partonIndex = 0;
    for ( ; partonsCandidates_itr != partonsCandidates.end(); ++partonsCandidates_itr,
	    partonIndex++ ) {
      double partonPt  = partonsCandidates[partonIndex]->pt();
      double partonEta = partonsCandidates[partonIndex]->eta();
      double partonPhi = partonsCandidates[partonIndex]->phi();
      //	      std::cout << "parton[" << partonIndex << "]: pt: " << partonPt
      //			<< " eta: " << partonEta
      //			<< " phi: " << partonPhi << std::endl;

      unsigned int jetIndex = 0;
      for ( std::vector<reco::JetBaseRef>::const_iterator jet_itr = jets.begin(); 
	    jet_itr != jets.end();
	    ++jet_itr, jetIndex++ ) {
	
	double jetEt  = (*jet_itr)->et();
	double jetEta = (*jet_itr)->eta();
	double jetPhi = (*jet_itr)->phi();
	//		std::cout << "jet[" << jetIndex << "]: et: " << jetEt
	//			  << " eta: " << jetEta
	//			  << " phi: " << jetPhi << std::endl;
	
	double jetParton_deltaEta = jetEta - partonEta;
	double jetParton_deltaPhi = deltaPhi(jetPhi,partonPhi);
	double jetParton_deltaR   = deltaR(jetEta,jetPhi,partonEta,partonPhi);      
      //		std::cout << "deltaR: " << deltaR << std::endl;
      
      //		if ( deltaR <= jetPartonDeltaRCut_ ) {
      if ( jetParton_deltaR <= closestJetDeltaRVec[partonIndex][3] ) {
	if ( jetParton_deltaR <= closestJetDeltaRVec[partonIndex][2] ) {
	  if ( jetParton_deltaR <= closestJetDeltaRVec[partonIndex][1] ) {
	    if ( jetParton_deltaR <= closestJetDeltaRVec[partonIndex][0] ) {
	      closestJetDeltaRVec[partonIndex][3] = closestJetDeltaRVec[partonIndex][2];
	      closestJetDeltaRVec[partonIndex][2] = closestJetDeltaRVec[partonIndex][1];
	      closestJetDeltaRVec[partonIndex][1] = closestJetDeltaRVec[partonIndex][0];
	      closestJetDeltaRVec[partonIndex][0] = jetParton_deltaR;
	      closestJetIndexVec[partonIndex][3] = closestJetIndexVec[partonIndex][2];
	      closestJetIndexVec[partonIndex][2] = closestJetIndexVec[partonIndex][1];
	      closestJetIndexVec[partonIndex][1] = closestJetIndexVec[partonIndex][0];
	      closestJetIndexVec[partonIndex][0] = jetIndex;
	    } else {
	      closestJetDeltaRVec[partonIndex][3] = closestJetDeltaRVec[partonIndex][2];
	      closestJetDeltaRVec[partonIndex][2] = closestJetDeltaRVec[partonIndex][1];
	      closestJetDeltaRVec[partonIndex][1] = jetParton_deltaR;
	      closestJetIndexVec[partonIndex][3] = closestJetIndexVec[partonIndex][2];
	      closestJetIndexVec[partonIndex][2] = closestJetIndexVec[partonIndex][1];
	      closestJetIndexVec[partonIndex][1] = jetIndex;
	    }
	  } else {
	    closestJetDeltaRVec[partonIndex][3] = closestJetDeltaRVec[partonIndex][2];
	    closestJetDeltaRVec[partonIndex][2] = jetParton_deltaR;
	    closestJetIndexVec[partonIndex][3] = closestJetIndexVec[partonIndex][2];
	    closestJetIndexVec[partonIndex][2] = jetIndex;
	  }
	} else {
	  closestJetDeltaRVec[partonIndex][3] = jetParton_deltaR;
	  closestJetIndexVec[partonIndex][3] = jetIndex;
	}
      }
    } // end loop over jet
    for (int index = 0; index < NJETS; index++ ) {
      unsigned int jetMatchedIndex = closestJetIndexVec[partonIndex][index];
      double etRes = (jets[jetMatchedIndex]->et()-partonPt)/partonPt;
      double jetMass = jets[jetMatchedIndex]->p4().M();
      std::cout << "closestJetDeltaRVec[" << partonIndex << "][" << index << ": " << jetMatchedIndex 
		<< " etRes (%): " << etRes*100.  
		<< " jetMass: " << jetMass
		<< std::endl;      
    }
  } // end loop over partons
  
  unsigned int nAss = 0;
  std::vector<bool> partonAss;
  for ( unsigned int partonIndex = 0; partonIndex < partonsCandidatesNumber; partonIndex++ ) 
    partonAss.push_back(kFALSE);
  
  for ( unsigned int assIndex = 0; assIndex < partonsCandidatesNumber; assIndex++ ) {
    //	      std::cout << "assIndex: " << assIndex << std::endl;
    double minDeltaR = 99.;
    int minPartonIndex = -1;
    // find the best DeltaR matching to find the best parton association in the collection
    for ( unsigned int partonIndex = 0; partonIndex < partonsCandidatesNumber; partonIndex++ ) {
      //		std::cout << "partonIndex: " << partonIndex << " partonAss: " << partonAss[partonIndex] << std::endl;
      if ( !partonAss[partonIndex] && closestJetDeltaRVec[partonIndex][0] <= minDeltaR ) {
	minDeltaR = closestJetDeltaRVec[partonIndex][0];
	minPartonIndex = partonIndex;
      }
    } // end loop over partons [not matched yet]
    
    // parton association
    partonAss[minPartonIndex] = kTRUE;
    nAss++;
    // save the matched jet into the proper vector
    unsigned int jetIndex = closestJetIndexVec[minPartonIndex][0];
    //	      std::cout << "minPartonIndex: " << minPartonIndex 
    //			<< " to jetIndex: " << jetIndex
    //			<< " => nAss: " << nAss << std::endl;
    if ( minPartonIndex == 0 || minPartonIndex == 1 ) {
      ZjetVec.push_back(jets[jetIndex]);
      ZjetInd.push_back(jetIndex);
    }
    else if ( minPartonIndex == 2 || minPartonIndex == 3 ) {
      TjetVec.push_back(jets[jetIndex]);
      TjetInd.push_back(jetIndex);
    }
    
    // in case of "non-biunivocity" pop-up jet associations belong to worst DeltaR parton
    for ( unsigned int partonIndex = 0; partonIndex < partonsCandidatesNumber; partonIndex++ ) {
      if ( partonAss[partonIndex] ) continue;
      //		std::cout << "partonIndex: " << partonIndex 
      //			  << " partonAss: " << partonAss[partonIndex] << std::endl;
      //		std::cout << "closestJetIndex: " << closestJetIndexVec[partonIndex][0] << std::endl;
      for ( unsigned int iAssIndex = 0; iAssIndex < partonsCandidatesNumber; iAssIndex++ ) {
	if ( closestJetIndexVec[partonIndex][iAssIndex] != closestJetIndexVec[minPartonIndex][0] ) continue;
	//		  std::cout << "iAssIndex: " << iAssIndex << std::endl;
	for ( unsigned int jAssIndex = iAssIndex+1; jAssIndex < partonsCandidatesNumber; jAssIndex++ ) {
	  //		    std::cout << "jAssIndex: " << jAssIndex << std::endl;
	  closestJetDeltaRVec[partonIndex][jAssIndex-1] = closestJetDeltaRVec[partonIndex][jAssIndex];
	  closestJetIndexVec[partonIndex][jAssIndex-1] = closestJetIndexVec[partonIndex][jAssIndex];
	}
      }
      //		std::cout << "closestJetIndex: " << closestJetIndexVec[partonIndex][0] << std::endl;
    } // end loop over partons [not matched yet and w/ the same jet]
    
  } // end loop over association index
  
  std::cout << "************************************" << std::endl;
  std::cout << "ZjetVec: " << ZjetVec.size() << std::endl;
  std::cout << "TjetVec: " << TjetVec.size() << std::endl;
  
  for ( unsigned index = 0; index < NZJETS; index++ ) {
    ZjetEta_  ->Fill(ZjetVec[index]->eta());
    ZjetPt_   ->Fill(ZjetVec[index]->pt());
    ZjetEt_   ->Fill(ZjetVec[index]->et());
    //    ZjetE_    ->Fill(ZjetVec[index]->e());
    ZjetMass_ ->Fill(ZjetVec[index]->mass());
  }
  ZjetsDeltaEta_ ->Fill(fabs(ZjetVec[0]->eta()-ZjetVec[1]->eta()));
  ZjetsDeltaR_ 	 ->Fill(deltaR(ZjetVec[0]->eta(),ZjetVec[0]->phi(),ZjetVec[1]->eta(),ZjetVec[1]->phi()));
  double recMass = (ZjetVec[0]->p4()+(ZjetVec[1]->p4())).mass();
  double mcMass = (hadronicZ[0]->p4()+(hadronicZ[1]->p4())).mass();
  ZjetsMass_           ->Fill(recMass);
  ZjetsMassResolution_ ->Fill(resolution(recMass,mcMass));
  double recPt = (ZjetVec[0]->p4()+(ZjetVec[1]->p4())).pt();
  double mcPt = (hadronicZ[0]->p4()+(hadronicZ[1]->p4())).pt();
  ZjetsPt_           ->Fill(recPt);
  ZjetsPtResolution_ ->Fill(resolution(recPt,mcPt));
  double recE = (ZjetVec[0]->p4()+(ZjetVec[1]->p4())).e();
  double mcE = (hadronicZ[0]->p4()+(hadronicZ[1]->p4())).e();
  ZjetsE_           ->Fill(recE);
  ZjetsEResolution_ ->Fill(resolution(recE,mcE));
  double recCollinearity = ZjetVec[0]->pz()*ZjetVec[1]->pz();
  double mcCollinearity = hadronicZ[0]->pz()*hadronicZ[1]->pz();
  ZjetsCollinearity_           ->Fill(recCollinearity);
  ZjetsCollinearityResolution_ ->Fill(resolution(recCollinearity,mcCollinearity));
  
    

}

// ------------ method called once each job just before starting event loop  ------------
void 
VBFHZZllbbJetMatching::beginJob(const edm::EventSetup&)
{

  std::cout << "VBFHZZllbbJetMatching::beginJob]" << std::endl;

  int nbin = 100;

  std::cout << "nbin: " << nbin << std::endl;

  // histo from HEPG information
  // ---------------------------
  edm::Service<TFileService> fs;

  eventsNumber_ = fs->make<TH1D>("eventsNumber","total number of events",1,0.,1.);

  std::cout << "eventsNumber" << std::endl;

  TFileDirectory mcSubDir = fs->mkdir( "HEPGinfo" );
  std::cout << "HEPGinfo" << std::endl;

  TFileDirectory mcZSubDir = mcSubDir.mkdir( "Zpartons" );

  ZpartonEta_ = mcZSubDir.make<TH1D>("ZpartonEta","partons from Z #eta", nbin, -8.,   8.);
  ZpartonPt_  = mcZSubDir.make<TH1D>("ZpartonPt", "partons from Z p_{T}",nbin,  0., 500.);
  ZpartonEt_  = mcZSubDir.make<TH1D>("ZpartonEt", "partons from Z E_{T}",nbin,  0., 500.);
  ZpartonE_   = mcZSubDir.make<TH1D>("ZpartonE",  "partons from Z E",    nbin,-10.,7000.);
  ZpartonsDeltaEta_     = mcZSubDir.make<TH1D>("ZpartonsDeltaEta",    "#Delta#eta between partons from Z",nbin,    0.,    10.);
  ZpartonsDeltaR_       = mcZSubDir.make<TH1D>("ZpartonsDeltaR",      "#DeltaR between partons from Z",   nbin,    0.,    10.);
  ZpartonsMass_         = mcZSubDir.make<TH1D>("ZpartonsMass",        "invariant mass of partons from Z", nbin,    0.,   120.);
  ZpartonsPt_           = mcZSubDir.make<TH1D>("ZpartonsPt",          "invariant p_{T} of partons from Z",nbin,    0.,   500.);
  ZpartonsEt_           = mcZSubDir.make<TH1D>("ZpartonsEt",          "invariant E_{T} of partons from Z",nbin,    0.,   500.);
  ZpartonsE_            = mcZSubDir.make<TH1D>("ZpartonsE",           "invariant E of partons from Z",    nbin,    0.,   500.);
  ZpartonsCollinearity_ = mcZSubDir.make<TH1D>("ZpartonsCollinearity","collinearity of partons from Z",   nbin,-7000.,  7000.);


//  ZpartonsDeltaRVSjetMass_ = mcZSubDir.make<TH2D>("ZpartonsDeltaRVSjetMass","jet mass VS #DeltaR between the 2 partons from Z",nbin,0.,0.8,nbin,0.,150.);
//  ZpartonsEtResVSjetMass_  = mcZSubDir.make<TH2D>("ZpartonsEtResVSjetMass", "jet mass VS E_{T} resolution",nbin,0.,1.,nbin,0.,150.);
//  ZpartonsDeltaRVSetRes_   = mcZSubDir.make<TH2D>("ZpartonsDeltaRVSetRes",  "E_{T} resolution VS #DeltaR between the 2 partons from Z",nbin,0.,0.8,nbin,0.,1.);
//  ZpartonsDeltaRVSjetMass_profile_ = mcZSubDir.make<TProfile>("ZpartonsDeltaRVSjetMass_profile","jet mass VS #DeltaR between the 2 partons from Z",nbin,0.,0.8,0.,150.);
//  ZpartonsEtResVSjetMass_profile_  = mcZSubDir.make<TProfile>("ZpartonsEtResVSjetMass_profile", "jet mass VS E_{T} resolution",nbin,0.,1.,0.,150.);
//  ZpartonsDeltaRVSetRes_profile_   = mcZSubDir.make<TProfile>("ZpartonsDeltaRVSetRes_profile",  "E_{T} resolution VS #DeltaR between the 2 partons from Z",nbin,0.,0.8,0.,1.);
//
  TFileDirectory mcTSubDir = mcSubDir.mkdir( "tagSystem" );
  TAGpartonEta_ = mcTSubDir.make<TH1D>("TAGpartonEta","tag partons #eta", nbin, -8.,   8.);
  TAGpartonPt_  = mcTSubDir.make<TH1D>("TAGpartonPt", "tag partons p_{T}",nbin,  0., 500.);
  TAGpartonEt_  = mcTSubDir.make<TH1D>("TAGpartonEt", "tag partons E_{T}",nbin,-10.,7000.);
  TAGpartonE_   = mcTSubDir.make<TH1D>("TAGpartonE",  "tag partons E",    nbin,-10.,7000.);
  TAGpartonsDeltaEta_     = mcTSubDir.make<TH1D>("TAGpartonsDeltaEta",    "#Delta#eta between tag partons",       nbin,    0.,   10.);
  TAGpartonsDeltaR_       = mcTSubDir.make<TH1D>("TAGpartonsDeltaR",      "#DeltaR between tag partons",          nbin,    0.,   10.);
  TAGpartonsMass_         = mcTSubDir.make<TH1D>("TAGpartonsMass",        "invariant mass of tag partons",        nbin,    0.,10000.);
  TAGpartonsPt_           = mcTSubDir.make<TH1D>("TAGpartonsPt",          "invariant p_{T} of tag partons system",nbin,    0., 7000.);
  TAGpartonsEt_           = mcTSubDir.make<TH1D>("TAGpartonsEt",          "invariant E_{T} of tag partons system",nbin,    0., 7000.);
  TAGpartonsE_            = mcTSubDir.make<TH1D>("TAGpartonsE",           "invariant E of tag partons system",    nbin,    0., 7000.);
  TAGpartonsCollinearity_ = mcTSubDir.make<TH1D>("TAGpartonsCollinearity","collinearity of tag partons system",   nbin,-7000., 7000.);

  TFileDirectory jetSubDir = fs->mkdir( "jet" );
  jetNumber_          = jetSubDir.make<TH1D>("jetNumber",         "number of jets",              nbin,  0. ,100. );
  jetUncorrEt_        = jetSubDir.make<TH1D>("jetUncorrEt",       "jet uncorrected E_{T}",       nbin,  0. ,500. );
  jetCorrEt_          = jetSubDir.make<TH1D>("jetCorrEt",         "jet correctedE_{T}",          nbin,  0. ,500. );
  jetUncorrPt_        = jetSubDir.make<TH1D>("jetUncorrPt",       "jet uncorrected p_{T}",       nbin,  0. ,500. );
  jetCorrPt_          = jetSubDir.make<TH1D>("jetCorrPt",         "jet corrected p_{T}",         nbin,  0. ,500. );
  jetPhi_             = jetSubDir.make<TH1D>("jetPhi",            "jet #Phi",                    nbin, -TMath::Pi(),TMath::Pi() );
  jetEta_             = jetSubDir.make<TH1D>("jetEta",            "jet #eta",                    nbin, -8.,   8. );
  jetMass_            = jetSubDir.make<TH1D>("jetMass",           "jet mass",                    nbin,  0. ,100. );
  jetEMfrac_          = jetSubDir.make<TH1D>("jetEMfrac",         "jet EM fraction",             nbin, -0.5,  1.5);
  jetHIGHEFFdiscr_    = jetSubDir.make<TH1D>("jetHIGHEFFdiscr",   "jet HIGHEFF discriminator",   nbin, -100., 10.);   	  
  jetHIGHPURdiscr_    = jetSubDir.make<TH1D>("jetHIGHPURdiscr",   "jet HIGHPUR discriminator",   nbin, -100., 10.);   	  
  jetCOMBSECVTXdiscr_ = jetSubDir.make<TH1D>("jetCOMBSECVTXdiscr","jet COMBSECVTX discriminator",nbin, -100., 10.);
  jetJETPROBdiscr_    = jetSubDir.make<TH1D>("jetJETPROBdiscr",   "jet JETPROB discriminator",   nbin, -100., 10.);   

  jetEMfracVSeta_             = jetSubDir.make<TH2D>("jetEMfracVSeta"             ,"jet EM fraction VS #eta",                    nbin, -0.5, 1.5, nbin,   -8.,   8. );	      
  jetEMfracVScorrEt_          = jetSubDir.make<TH2D>("jetEMfracVScorrEt"          ,"jet EM fraction VS corrected E_{T}",         nbin, -0.5, 1.5, nbin,    0., 500. );	      
  jetEMfracVScorrPt_          = jetSubDir.make<TH2D>("jetEMfracVScorrPt"          ,"jet EM fraction VS corrected p_{T}",         nbin, -0.5, 1.5, nbin,    0., 500. );	      
  jetEMfracVSuncorrEt_        = jetSubDir.make<TH2D>("jetEMfracVSuncorrEt"        ,"jet EM fraction VS uncorrected E_{T}",       nbin, -0.5, 1.5, nbin,    0., 500. );	      
  jetEMfracVSuncorrPt_        = jetSubDir.make<TH2D>("jetEMfracVSuncorrPt"        ,"jet EM fraction VS uncorrected p_{T}",       nbin, -0.5, 1.5, nbin,    0., 500. );	      
  jetEMfracVShighEFFdiscr_    = jetSubDir.make<TH2D>("jetEMfracVShighEFFdiscr"    ,"jet EM fraction VS HIGHEFF discriminator",   nbin, -0.5, 1.5, nbin, -100.,   10.);    
  jetEMfracVShighPURdiscr_    = jetSubDir.make<TH2D>("jetEMfracVShighPURdiscr"    ,"jet EM fraction VS HIGHPUR discriminator",   nbin, -0.5, 1.5, nbin, -100.,   10.);    
  jetEMfracVScomboSECVTXdiscr_= jetSubDir.make<TH2D>("jetEMfracVScomboSECVTXdiscr","jet EM fraction VS COMBSECVTX discriminator",nbin, -0.5, 1.5, nbin, -100.,   10.);
  jetEMfracVSjetPROBdiscr_    = jetSubDir.make<TH2D>("jetEMfracVSjetPROBdiscr"    ,"jet EM fraction VS JETPROB discriminator",   nbin, -0.5, 1.5, nbin, -100.,   10.);    
  jetEMfracVSeta_profile_             = jetSubDir.make<TProfile>("jetEMfracVSeta_profile"             ,"jet EM fraction VS #eta",                    nbin, -0.5, 10.,   -8.,   8. );	      
  jetEMfracVScorrEt_profile_          = jetSubDir.make<TProfile>("jetEMfracVScorrEt_profile"          ,"jet EM fraction VS corrected E_{T}",         nbin, -0.5, 10.,    0., 500. );	      
  jetEMfracVScorrPt_profile_          = jetSubDir.make<TProfile>("jetEMfracVScorrPt_profile"          ,"jet EM fraction VS corrected p_{T}",         nbin, -0.5, 10.,    0., 500. );	      
  jetEMfracVSuncorrEt_profile_        = jetSubDir.make<TProfile>("jetEMfracVSuncorrEt_profile"        ,"jet EM fraction VS uncorrected E_{T}",       nbin, -0.5, 10.,    0., 500. );	      
  jetEMfracVSuncorrPt_profile_        = jetSubDir.make<TProfile>("jetEMfracVSuncorrPt_profile"        ,"jet EM fraction VS uncorrected p_{T}",       nbin, -0.5, 10.,    0., 500. );	      
  jetEMfracVShighEFFdiscr_profile_    = jetSubDir.make<TProfile>("jetEMfracVShighEFFdiscr_profile"    ,"jet EM fraction VS HIGHEFF discriminator",   nbin, -0.5, 10., -100.,   1.5);    
  jetEMfracVShighPURdiscr_profile_    = jetSubDir.make<TProfile>("jetEMfracVShighPURdiscr_profile"    ,"jet EM fraction VS HIGHPUR discriminator",   nbin, -0.5, 10., -100.,   1.5);    
  jetEMfracVScomboSECVTXdiscr_profile_= jetSubDir.make<TProfile>("jetEMfracVScomboSECVTXdiscr_profile","jet EM fraction VS COMBSECVTX discriminator",nbin, -0.5, 10., -100.,   1.5);
  jetEMfracVSjetPROBdiscr_profile_    = jetSubDir.make<TProfile>("jetEMfracVSjetPROBdiscr_profile"    ,"jet EM fraction VS JETPROB discriminator",   nbin, -0.5, 10., -100.,   1.5);    


 muonJetMatchedNumber_ = jetSubDir.make<TH1D>("muonJetMatchedNumber","number of jets w/ at least 1 muon in a distance #DeltaR<0.5", 10,0.,10.0);

 jetParton_deltaR_      = jetSubDir.make<TH1D>("jetParton_deltaR",     "#DeltaR^{2} between jet and partons from Z",    nbin,0.,4.0);
 jetParton_deltaR_zoom_ = jetSubDir.make<TH1D>("jetParton_deltaR_zoom","#DeltaR^{2} between jet and partons from Z",    nbin,0.,0.1);
 jetParton_deltaRmax_   = jetSubDir.make<TH1D>("jetParton_deltaRmax",  "max #DeltaR^{2} between jet and partons from Z",nbin,0.,4.0);

 jetParton_deltaEtVSdeltaR_         = jetSubDir.make<TH2D>("jetParton_deltaEtVSdeltaR",        "#DeltaE_{T} VS #DeltaR^{2} between partons and matched jets",          nbin,0.,0.04,nbin,-50.,50.);
 jetParton_deltaEtmeanVSdeltaRmean_ = jetSubDir.make<TH2D>("jetParton_deltaEtmeanVSdeltaRmean","mean #DeltaE_{T} VS mean #DeltaR^{2} between partons and matched jets",nbin,0.,0.04,nbin,  0.,50.);
 jetParton_deltaEVSdeltaR_   = jetSubDir.make<TH2D>("jetParton_deltaEVSdeltaR",  "#DeltaE VS #DeltaR^{2} between partons and matched jets",    nbin,0.,0.04,nbin,-50.,50.);
 jetParton_deltaEtaVSdeltaR_ = jetSubDir.make<TH2D>("jetParton_deltaEtaVSdeltaR","#Delta#eta VS #DeltaR^{2} between partons and matched jets", nbin,0.,0.04,nbin, -5., 5.);
 jetParton_deltaEtVSdeltaR_profile_         = jetSubDir.make<TProfile>("jetParton_deltaEtVSdeltaR_profile",        "#DeltaE_{T} VS #DeltaR^{2} between partons and matched jets",          nbin,0.,0.04,-50.,50.);
 jetParton_deltaEtmeanVSdeltaRmean_profile_ = jetSubDir.make<TProfile>("jetParton_deltaEtmeanVSdeltaRmean_profile","mean #DeltaE_{T} VS mean #DeltaR^{2} between partons and matched jets",nbin,0.,0.04,  0.,50.);
 jetParton_deltaEVSdeltaR_profile_   = jetSubDir.make<TProfile>("jetParton_deltaEVSdeltaR_profile",  "#DeltaE VS #DeltaR^{2} between partons and matched jets",    nbin,0.,0.04,-50.,50.);
 jetParton_deltaEtaVSdeltaR_profile_ = jetSubDir.make<TProfile>("jetParton_deltaEtaVSdeltaR_profile","#Delta#eta VS #DeltaR^{2} between partons and matched jets", nbin,0.,0.04, -5., 5.);
 jetParton_deltaEta_ = jetSubDir.make<TH1D>("jetParton_deltaEta","#Delta#eta between partons and matched jets", nbin, -5., 5.);
 jetParton_deltaEt_  = jetSubDir.make<TH1D>("jetParton_deltaEt", "#DeltaE_{T} between partons and matched jets",nbin,-50.,50.);
 jetParton_deltaEtRes_  = jetSubDir.make<TH1D>("jetParton_deltaEtRes", "E_{T} resolution between partons and matched jets",nbin,-1.,1.);
 jetParton_004deltaEtRes_ = jetSubDir.make<TH1D>("jetParton_004deltaEtRes", "E_{T} resolution between partons and matched jets (#DeltaR<=0.2)",nbin,-1.,1.);

 jetParton_deltaRVSdeltaR_ = jetSubDir.make<TH2D>("jetParton_deltaRVSdeltaR","#DeltaR between the 2jets VS #DeltaR between the 2 partons",nbin,0.,5.,nbin,0.,5.);
 jetParton_deltaRVSdeltaR_profile_ = jetSubDir.make<TProfile>("jetParton_deltaRVSdeltaR_profile","#DeltaR between the 2jets VS #DeltaR between the 2 partons",nbin,0.,5.,0.,5.);

 TFileDirectory jetZSubDir = fs->mkdir( "Zjet" );
 ZjetEta_  = jetZSubDir.make<TH1D>("ZjetEta", "#eta of jets matched to Z partons", nbin,-8.,8.);
 ZjetPt_   = jetZSubDir.make<TH1D>("ZjetPt",  "p_{T} of jets matched to Z partons",nbin,0.,500.);
 ZjetEt_   = jetZSubDir.make<TH1D>("ZjetEt",  "E_{T} of jets matched to Z partons",nbin,0.,500.);
 ZjetE_    = jetZSubDir.make<TH1D>("ZjetE",   "E of jets matched to Z partons",    nbin,0.,500.);
 ZjetMass_ = jetZSubDir.make<TH1D>("ZjetMass","mass of jets matched to Z partons", nbin,0.,100.); 
 ZjetsDeltaEta_ = jetZSubDir.make<TH1D>("ZjetsDeltaEta","#Delta#eta between jets matched to Z partons", nbin,0.,10.);
 ZjetsDeltaR_ = jetZSubDir.make<TH1D>("ZjetsDeltaR","#DeltaR between jets matched to Z partons", nbin,0.,10.);
 ZjetsMass_           = jetZSubDir.make<TH1D>("ZjetsMass",          "reconstructed Z mass from jets matched to Z partons",           nbin, 0.,150.);
 ZjetsMassResolution_ = jetZSubDir.make<TH1D>("ZjetsMassResolution","reconstructed Z mass resolution from jets matched to Z partons",nbin,-1.,  1.);
 ZjetsPt_           = jetZSubDir.make<TH1D>("ZjetsPt",          "reconstructed Z p_{T} from jets matched to Z partons",           nbin, 0.,500.);
 ZjetsPtResolution_ = jetZSubDir.make<TH1D>("ZjetsPtResolution","reconstructed Z p_{T} resolution from jets matched to Z partons",nbin,-1.,  1.);
 ZjetsE_           = jetZSubDir.make<TH1D>("ZjetsE",          "reconstructed Z E from jets matched to Z partons",           nbin, 0.,500.);
 ZjetsEResolution_ = jetZSubDir.make<TH1D>("ZjetsEResolution","reconstructed Z E resolution from jets matched to Z partons",nbin,-1.,  1.);
 ZjetsCollinearity_           = jetZSubDir.make<TH1D>("ZjetsCollinearity",          "reconstructed collinearity between jets matched to Z partons",        nbin, -6000.,6000.);
 ZjetsCollinearityResolution_ = jetZSubDir.make<TH1D>("ZjetsCollinearityResolution","reconstructed collinearity resolution from jets matched to Z partons",nbin,    -1.,   1.);

 TFileDirectory jetTAGSubDir = fs->mkdir( "TAGjet" );
 TAGjetEta_  = jetTAGSubDir.make<TH1D>("TAGjetEta", "#eta of jets matched to TAG partons", nbin,-8.,8.);
 TAGjetPt_   = jetTAGSubDir.make<TH1D>("TAGjetPt",  "p_{T} of jets matched to TAG partons",nbin,0.,500.);
 TAGjetEt_   = jetTAGSubDir.make<TH1D>("TAGjetEt",  "E_{T} of jets matched to TAG partons",nbin,0.,500.);
 TAGjetE_    = jetTAGSubDir.make<TH1D>("TAGjetE",   "E of jets matched to TAG partons",    nbin,0.,500.);
 TAGjetMass_ = jetTAGSubDir.make<TH1D>("TAGjetMass","mass of jets matched to TAG partons", nbin,0.,100.); 
 TAGjetsDeltaEta_ = jetTAGSubDir.make<TH1D>("TAGjetsDeltaEta","#Delta#eta between jets matched to TAG partons", nbin,0.,10.);
 TAGjetsDeltaR_ = jetTAGSubDir.make<TH1D>("TAGjetsDeltaR","#DeltaR between jets matched to TAG partons", nbin,0.,10.);
 TAGjetsMass_           = jetTAGSubDir.make<TH1D>("TAGjetsMass",          "reconstructed TAG system mass from jets matched to TAG partons",           nbin, 0.,150.);
 TAGjetsMassResolution_ = jetTAGSubDir.make<TH1D>("TAGjetsMassResolution","reconstructed TAG system mass resolution from jets matched to TAG partons",nbin,-1.,  1.);
 TAGjetsPt_           = jetTAGSubDir.make<TH1D>("TAGjetsPt",          "reconstructed TAG system p_{T} from jets matched to TAG partons",           nbin, 0.,500.);
 TAGjetsPtResolution_ = jetTAGSubDir.make<TH1D>("TAGjetsPtResolution","reconstructed TAG system p_{T} resolution from jets matched to TAG partons",nbin,-1.,  1.);
 TAGjetsE_           = jetTAGSubDir.make<TH1D>("TAGjetsE",          "reconstructed TAG system E from jets matched to TAG partons",           nbin, 0.,500.);
 TAGjetsEResolution_ = jetTAGSubDir.make<TH1D>("TAGjetsEResolution","reconstructed TAG system E resolution from jets matched to TAG partons",nbin,-1.,  1.);
 TAGjetsCollinearity_           = jetTAGSubDir.make<TH1D>("TAGjetsCollinearity",          "reconstructed collinearity between jets matched to TAG partons",        nbin, -6000.,6000.);
 TAGjetsCollinearityResolution_ = jetTAGSubDir.make<TH1D>("TAGjetsCollinearityResolution","reconstructed collinearity resolution from jets matched to TAG partons",nbin,    -1.,   1.);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
VBFHZZllbbJetMatching::endJob() {

  std::cout << "VBFHZZllbbJetMatching::endJob" << std::endl;

  // Fill histograms
  // ---------------
  eventsNumber_->SetBinContent(1,eventcounter_);

}

//define this as a plug-in
DEFINE_FWK_MODULE(VBFHZZllbbJetMatching);


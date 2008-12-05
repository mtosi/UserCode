//
// Original Author:  Mia Tosi
//         Created:  Fri Feb 22 17:56:22 CET 2008
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

#include "AnalysisExamples/HiggsZZAnalyzer/interface/VBFHZZanalyzer.h"

#include "AnalysisExamples/AnalysisClasses/interface/DeltaPhi.h"
#include "AnalysisExamples/AnalysisClasses/interface/DeltaR.h"
#include "AnalysisExamples/AnalysisClasses/interface/DiParticleMass.h"

#include "AnalysisExamples/AnalysisObjects/interface/PythiaParticleIndex.h"
#include "AnalysisExamples/AnalysisObjects/interface/ParticlesMass.h"
#include "AnalysisExamples/AnalysisObjects/interface/ParticlesCharge.h"

// For file output
// ---------------
#include <fstream>
#include <sstream>
#include <cmath>
#include <memory>
#include <string>
#include <iostream>
#include <iomanip>


VBFHZZanalyzer::VBFHZZanalyzer(const edm::ParameterSet& iConfig) :
  conf_( iConfig ),
  offlineJetLabel_(         iConfig.getUntrackedParameter<edm::InputTag>( "OfflineJets"                  ) ),
  //  offlineMEtLabel_(         iConfig.getUntrackedParameter<edm::InputTag>( "OfflineMEt"                   ) ),
  globalMuonLabel_(         iConfig.getUntrackedParameter<edm::InputTag>( "GlobalMuons"                  ) ),
  simpleElectronLabel_(     iConfig.getUntrackedParameter<edm::InputTag>( "SimpleElectrons"              ) ),
  //  simpleTauLabel_(          iConfig.getUntrackedParameter<edm::InputTag>( "SimpleTaus"                   ) ),
  //  combinedSVBJetTagsLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "combinedSVBJetTags"           ) ),
  //  simVtxLabel_(             iConfig.getUntrackedParameter<edm::InputTag>( "SimVtx"                       ) ),
  MCParticleLabel_(         iConfig.getUntrackedParameter<std::string>(   "MCParticles"                  ) ),
  leptonPtCut_(	       iConfig.getUntrackedParameter<double>("leptonPtCut"        ) ),
  jetEtCut_(           iConfig.getUntrackedParameter<double>("jetEtCut"           ) ),
  jetPartonDeltaR2Cut_(iConfig.getUntrackedParameter<double>("jetPartonDeltaR2Cut") ),
  jetLeptonDeltaRCut_( iConfig.getUntrackedParameter<double>("jetLeptonDeltaRCut" ) ),
  jetEMfracCut_(       iConfig.getUntrackedParameter<double>("jetEMfracCut"       ) )
{
  // Now do what ever initialization is needed
  // -----------------------------------------

  //
  // constants, enums and typedefs
  nZ_           = 2;
  nZleptons_    = 2;
  nZjets_       = 2;
  nforwardjets_ = 2;
  njets_        = 4;

  eventcounter_             = 0;
  goodLeptonicEventCounter_ = 0;
  twoJetsAboveEtCutEventsCounter_;
  twoJetsAboveEtCut_DeltaR2016EventsCounter_;
  twoJetsAboveEtCut_DeltaR2CutEventsCounter_;
  fourJetsAboveEtCutEventsCounter_;

  nbin_         = 100;

  null_XYZTLorentzVector_ = math::XYZTLorentzVector(0.,0.,0.,0.);
  null_globalMuon_ = GlobalMuon( 0.,0.,0.,0,
				 null_XYZTLorentzVector_,
				 math::XYZPoint(0.,0.,0.),
				 0.,0.,0.,0.,0.,0.,0.,0,0,0.,0,0.,0.);
  null_offlinejet_ = OfflineJet( 0.,0.,0.,0.,0.,
				 null_XYZTLorentzVector_,
				 math::XYZPoint(0.,0.,0.),
				 0.,0.,0.,0, 0., 0.,0, 0., 0.,0, 0., 0. );
  gROOT->Time();

}


VBFHZZanalyzer::~VBFHZZanalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

  // write out the file
  // ------------------
  OutputFile->Write();

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
VBFHZZanalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
  using namespace std;
  using namespace anaobj;

  eventcounter_++;
  if ( eventcounter_/100 == float(eventcounter_)/100. ) {
    std::cout << "Event number " << eventcounter_ << std::endl;
  }

  // Global muons
  // ------------
  edm::Handle < GlobalMuonCollection > globalMuons;
  iEvent.getByLabel( globalMuonLabel_, globalMuons );
  // SimpleElectrons
  // ---------------
  edm::Handle < SimpleElectronCollection > simpleElectrons;
  iEvent.getByLabel( simpleElectronLabel_, simpleElectrons);
  // Calorimeter jets
  // ----------------
  edm::Handle<OfflineJetCollection> offlineJets;
  iEvent.getByLabel( offlineJetLabel_, offlineJets );
  // MCParticle
  // ----------
  edm::Handle < CandidateCollection > MCparticles;
  iEvent.getByLabel( MCParticleLabel_, MCparticles );

  unsigned globalMuonsSize     = globalMuons->size();
  unsigned simpleElectronsSize = simpleElectrons->size();
  unsigned offlineJetsSize     = offlineJets->size();
  if ( ( globalMuonsSize >= nZleptons_ || simpleElectronsSize >= nZleptons_ ) &&
       offlineJetsSize >= njets_ ) {
    // count number of events w/ at least 2 leptons and 4 jets
    // as we expect in H->ZZ->lljj via VBF

    /////////////////////////////////// QUICK LOOP over LEPTONS ////////////////////////////
    // Global muons
    // ------------
    std::vector<GlobalMuon> goodMuonVec;
    muonNumber_->Fill(globalMuons->size());
    bool goodMuonPair = kFALSE;
    GlobalMuonCollection::const_iterator globalMuon_itr = globalMuons->begin();
    for ( ; globalMuon_itr != globalMuons->end(); ++globalMuon_itr ) {
      double tmpMuonPt  = globalMuon_itr->pt();
      double tmpMuonEta = globalMuon_itr->eta();
      double tmpMuonPhi = globalMuon_itr->phi();
      if ( tmpMuonPt >= leptonPtCut_ ) {
	if (!goodMuonPair) {
	  if ( globalMuon_itr != (globalMuons->end()-1) ) {
	    GlobalMuonCollection::const_iterator globalMuon_itr2 = globalMuon_itr+1;
	    for ( ; globalMuon_itr2 != globalMuons->end(); ++globalMuon_itr2 ){
	      if ( globalMuon_itr->charge()*globalMuon_itr2->charge() < 0 ) {
		goodMuonVec.push_back(*globalMuon_itr);
		goodMuonVec.push_back(*globalMuon_itr2);
		goodMuonPair = kTRUE;
	      }
	    }
	  }
	}
      }
    }
    //    std::cout << "goodMuonVec: "  << goodMuonVec.size() << std::endl;
    if (goodMuonVec.size() != 0) std::cout << "dimuonMass: " << (goodMuonVec[0].p4()+goodMuonVec[1].p4()).M() << std::endl;
    
    // SimpleElectrons
    // ---------------
    std::vector<SimpleElectron> goodElectronVec;
    electronNumber_->Fill(simpleElectrons->size());
    SimpleElectronCollection::const_iterator simpleElectron_itr = simpleElectrons->begin();
    for ( ; simpleElectron_itr != simpleElectrons->end(); ++simpleElectron_itr ) {
      double tmpElectronEt  = simpleElectron_itr->et();
      double tmpElectronEta = simpleElectron_itr->eta();
      double tmpElectronPhi = simpleElectron_itr->phi();
      double tmpElectronPt  = simpleElectron_itr->pt();
      if ( tmpElectronPt >= leptonPtCut_ ) {
      }
    }
  
    if ( goodMuonVec.size() >= nZleptons_ || goodElectronVec.size() >= nZleptons_ ) {

      goodLeptonicEventCounter_++;

      /////////////////////////////////// HEPG analysis ////////////////////////////
      std::vector< const Candidate* > Z1Candidates;
      std::vector< const Candidate* > Z2Candidates;
      std::vector< const Candidate* > TagCandidates;
      bool tmpHmother = kFALSE;
      const Candidate* tmpZMother = 0;
      // pre-loop to find Higgs mothers
      CandidateCollection::const_iterator MCparticle = MCparticles->begin();
      int MCparticleCounter = 0;
      for ( ; MCparticle != MCparticles->end(); ++MCparticle, MCparticleCounter++ ) {    
	int tmpParticleId = MCparticle->pdgId();
	if ( fabs(tmpParticleId) > pythiaH_ ) continue;
	int tmpMotherId = 0;
	if( MCparticle->mother() != 0 ) tmpMotherId = MCparticle->mother()->pdgId();
	// candidates from Zs
	if ( tmpMotherId == pythiaZ_ && tmpParticleId != tmpMotherId ) {
	  if ( tmpZMother == 0 ) {
	    Z1Candidates.push_back(&*MCparticle);
	    tmpZMother = MCparticle->mother();
	  } else {
	    if ( MCparticle->mother() == tmpZMother ) {
	      Z1Candidates.push_back(&*MCparticle);
	    } else {
	      Z2Candidates.push_back(&*MCparticle);
	    }
	  }
	}
	// candidates for tag quarks
	int tmpNumberOfDaughters = MCparticle->numberOfDaughters();
	for (int daughterIndex = 0; daughterIndex != tmpNumberOfDaughters; daughterIndex++ ) {
	  if ( MCparticle->daughter(daughterIndex)->pdgId() == pythiaH_ ) tmpHmother = kTRUE;
	}
	if ( tmpHmother ) {
	  for (int daughterIndex = 0; daughterIndex != tmpNumberOfDaughters; daughterIndex++ ) {
	    if ( MCparticle->daughter(daughterIndex)->pdgId() != pythiaH_ ) {
	      if ( MCparticle->daughter(daughterIndex)->numberOfMothers() == 2 ) {
		if (TagCandidates.size() < 2 ) TagCandidates.push_back(MCparticle->daughter(daughterIndex));
		else {
		  if (MCparticle->daughter(daughterIndex)->pdgId() != TagCandidates[0]->pdgId() &&
		      MCparticle->daughter(daughterIndex)->pdgId() != TagCandidates[1]->pdgId() ) 
		    std::cout << "WRONG tag partons matching" << std::endl;
		}
	      }
	    }
	  }
	  tmpHmother = kFALSE;
	}
      }
      //      std::cout << "candidates frrm Z1: "      << Z1Candidates.size()  << std::endl;
      //      std::cout << "candidates frrm Z2: "      << Z2Candidates.size()  << std::endl;
      //      std::cout << "candidates for tag jets: " << TagCandidates.size() << std::endl;


      // manage info of particles from Z decay
      // -------------------------------------
      std::vector< const Candidate* > hadronicZcandidates;
      std::vector< const Candidate* > leptonicZcandidates;
      Particle * hadronicZ = 0;
      Particle * leptonicZ = 0;
      int    Z1particlesId = 0;
      int    Z2particlesId = 0;
      double Z1mass = 0.;
      double Z2mass = 0.;
      if ( Z1Candidates.size() == 2 && Z2Candidates.size() == 2 ) {
	Z1particlesId = fabs(int(Z1Candidates[0]->pdgId()));
	Z2particlesId = fabs(int(Z2Candidates[0]->pdgId()));
	Particle * Z1particle = new Particle(ZCharge_,
					     Z1Candidates[0]->p4()+(Z1Candidates[1]->p4()),
					     Z1Candidates[0]->vertex(),
					     pythiaZ_,0,true);
	Particle * Z2particle = new Particle(ZCharge_,
					     Z2Candidates[0]->p4()+(Z2Candidates[1]->p4()),
					     Z2Candidates[0]->vertex(),
					     pythiaZ_,0,true);
	Z1mass = Z1particle->mass();
	Z2mass = Z2particle->mass();
	if ( fabs(Z1particlesId) <= pythiat_ && 
	     (fabs(Z2particlesId) == pythiae_ || fabs(Z2particlesId) == pythiamu_ ) ) {
	  hadronicZ = Z1particle;
	  leptonicZ = Z2particle;
	  hadronicZcandidates.push_back(Z1Candidates[0]);
	  hadronicZcandidates.push_back(Z1Candidates[1]);
	}
	else if ( fabs(Z2particlesId) <= pythiat_ && 
		  (fabs(Z1particlesId) == pythiae_ || fabs(Z1particlesId) == pythiamu_ ) ) {
	  hadronicZ = Z2particle;
	  leptonicZ = Z1particle;
	  hadronicZcandidates.push_back(Z2Candidates[0]);
	  hadronicZcandidates.push_back(Z2Candidates[1]);
	}
      }

      //      std::cout << "hadronicZcandidates: " << hadronicZcandidates.size() << std::endl;

      ////////////////////////////////////// PLOT HEPG DISTRIBUTION ////////////////////////////
      if ( hadronicZcandidates.size() != 0 ) {
	for ( std::vector< const Candidate* >::const_iterator hadronicZcandidates_itr = hadronicZcandidates.begin();
	      hadronicZcandidates_itr != hadronicZcandidates.end(); ++hadronicZcandidates_itr) {
	  double tmp_ZpartonEta = (*hadronicZcandidates_itr)->eta();
	  double tmp_ZpartonPt  = (*hadronicZcandidates_itr)->pt();
	  double tmp_ZpartonE   = (*hadronicZcandidates_itr)->energy();
	  ZpartonEta_->Fill(tmp_ZpartonEta);
	  ZpartonPt_->Fill( tmp_ZpartonPt );
	  ZpartonE_->Fill(  tmp_ZpartonE  );
	}
	double tmp_ZpartonsDeltaEta     = fabs(hadronicZcandidates[0]->eta()-hadronicZcandidates[1]->eta());
	double tmp_ZpartonsMass         = (hadronicZcandidates[0]->p4()+(hadronicZcandidates[1]->p4())).mass();
	double tmp_ZpartonsCollinearity = hadronicZcandidates[0]->pz()*hadronicZcandidates[1]->pz();
	ZpartonsDeltaEta_->Fill(    tmp_ZpartonsDeltaEta    );
	ZpartonsMass_->Fill(        tmp_ZpartonsMass        );
	ZpartonsCollinearity_->Fill(tmp_ZpartonsCollinearity);
      }
      if ( TagCandidates.size() != 0 ) {
	std::vector< const Candidate* >::const_iterator TagCandidates_itr = TagCandidates.begin();
	for ( ; TagCandidates_itr != TagCandidates.end(); ++TagCandidates_itr) {
	  double tmp_TagPartonEta = (*TagCandidates_itr)->eta();
	  double tmp_TagPartonPt  = (*TagCandidates_itr)->pt();
	  double tmp_TagPartonE   = (*TagCandidates_itr)->energy();
	  TagPartonEta_->Fill(tmp_TagPartonEta);
	  TagPartonPt_->Fill( tmp_TagPartonPt );
	  TagPartonE_->Fill( tmp_TagPartonE  );
	}
	double tmp_TagPartonsDeltaEta     = fabs(TagCandidates[0]->eta()-TagCandidates[1]->eta());
	double tmp_TagPartonsMass         = (TagCandidates[0]->p4()+(TagCandidates[1]->p4())).mass();
	double tmp_TagPartonsCollinearity = TagCandidates[0]->pz()*TagCandidates[1]->pz();
	TagPartonsDeltaEta_->Fill(    tmp_TagPartonsDeltaEta    );
	TagPartonsMass_->Fill(        tmp_TagPartonsMass        );
	TagPartonsCollinearity_->Fill(tmp_TagPartonsCollinearity);
      }

      if ( hadronicZcandidates.size() != 0 ) {
	//    /////////////////////////////////// JET-MUON MATCHING ////////////////////////////
	if ( globalMuonsSize != 0 ) {
	  int muonsMatchedCounter = 0;
	  GlobalMuonCollection::const_iterator globalMuon_itr = globalMuons->begin();
	  for ( ; globalMuon_itr != globalMuons->end(); ++globalMuon_itr ) {
	    int muonJetMatchedCounter = 0;
	    double tmp_muonEta = globalMuon_itr->eta();
	    double tmp_muonPhi = globalMuon_itr->eta();
	    OfflineJetCollection::const_iterator offlinejet = offlineJets->begin();
	    for ( ; offlinejet != offlineJets->end(); ++offlinejet ) {
	      double tmp_jetEta = offlinejet->eta();
	      double tmp_jetPhi = offlinejet->phi();
	      double tmp_jetMuon_deltaEta = tmp_jetEta-tmp_muonEta;
	      double tmp_jetMuon_deltaPhi = DeltaPhi(tmp_jetPhi,tmp_muonPhi);
	      double tmp_jetMuon_deltaR2  = pow(tmp_jetMuon_deltaEta,2)+pow(tmp_jetMuon_deltaPhi,2);
	      if ( tmp_jetMuon_deltaR2 < 0.25 ) muonJetMatchedCounter++;
	    }
	    if (muonJetMatchedCounter != 0 ) muonsMatchedCounter++;
	    muonJetMatchedNumber_->Fill(muonJetMatchedCounter);
	    muonsMatchedNumber_->Fill(double(muonsMatchedCounter/globalMuonsSize));
	  }
	}

	//    /////////////////////////////////// JET-PARTON MATCHING ////////////////////////////
	// Calorimeter jets
	// ----------------
	jetNumber_->Fill(offlineJets->size());


	// tom ideas
	// ---------
	// Et cut loop
	// -----------
	for ( int iet=0; iet<41; iet++ ) {
	  double etcut = 10.+(double)iet;

	  std::vector< OfflineJet > aboveEtCutJetVec;
	  OfflineJetCollection::const_iterator offlinejet = offlineJets->begin();
	  for ( ; offlinejet != offlineJets->end(); ++offlinejet ) {
	    double tmpOfflinejetEt  = offlinejet->et();
	    double tmpOfflinejetEta = offlinejet->eta();
	    if ( tmpOfflinejetEt >= etcut ) {
	      aboveEtCutJetVec.push_back(*offlinejet);
	      if ( etcut == jetEtCut_ ) jetEta_aboveEtCut_->Fill(tmpOfflinejetEta);
	    }
	  }
	  if ( etcut == jetEtCut_ ) {
	    jetNumber_aboveEtCut_->Fill(aboveEtCutJetVec.size());
	    if ( aboveEtCutJetVec.size() >= nZjets_ ) twoJetsAboveEtCutEventsCounter_++;
	    if ( aboveEtCutJetVec.size() >= njets_ )  fourJetsAboveEtCutEventsCounter_++;
	  }
	  if ( aboveEtCutJetVec.size() < nZjets_ ) continue;
	
	  std::vector< std::vector<OfflineJet> > partonClosestJetsVec;
	  std::vector< std::vector<double> >     partonClosestJetsDeltaR2Vec;
	  std::vector< const Candidate* >::const_iterator hadronicZcandidates_itr = hadronicZcandidates.begin();
	  for ( ; hadronicZcandidates_itr != hadronicZcandidates.end(); ++hadronicZcandidates_itr) {
	    std::vector<OfflineJet> null_matchedJetsVec;
	    std::vector<double>     null_deltaR2Vec;
	    for ( int index = 0; index != hadronicZcandidates.size(); index++ ) {
	      null_matchedJetsVec.push_back(null_offlinejet_);
	      null_deltaR2Vec.push_back(99.);
	    }
	    partonClosestJetsVec.push_back(null_matchedJetsVec);
	    partonClosestJetsDeltaR2Vec.push_back(null_deltaR2Vec);
	  }

	  std::cout << "----------><----------" << std::endl;
	  hadronicZcandidates_itr = hadronicZcandidates.begin();
	  int partonIndex = 0;
	  for ( ; hadronicZcandidates_itr != hadronicZcandidates.end(); ++hadronicZcandidates_itr, 
		  partonIndex++ ) {
	    std::cout << "parton[" << partonIndex << "]: pt: " << hadronicZcandidates[partonIndex]->pt()
		      << " eta: " << hadronicZcandidates[partonIndex]->eta()
		      << " phi: " << hadronicZcandidates[partonIndex]->phi() << std::endl;
	    OfflineJetCollection::const_iterator offlinejet = offlineJets->begin();
	    int jetIndex = 0;
	    for ( ; offlinejet != offlineJets->end(); ++offlinejet, jetIndex++ ) {
	      double tmpOfflinejetEt = offlinejet->et();
	      double tmpOfflinejetEta = offlinejet->eta();
	      double tmp_jetParton_deltaR2 = -99.;
	      //	    if ( tmpOfflinejetEt >= jetEtCut_ ) {
	      if ( tmpOfflinejetEt >= etcut ) {
		std::cout << "jet[" << jetIndex << "]: et: " << offlinejet->et()
			  << " eta: " << offlinejet->eta()
			  << " phi: " << offlinejet->phi() << std::endl;
		double tmpOfflinejetPhi = offlinejet->phi();
		double tmp_jetParton_deltaEta = tmpOfflinejetEta-(*hadronicZcandidates_itr)->eta();
		double tmp_jetParton_deltaPhi = DeltaPhi(tmpOfflinejetPhi,(*hadronicZcandidates_itr)->phi());
		tmp_jetParton_deltaR2 = pow(tmp_jetParton_deltaEta,2)+pow(tmp_jetParton_deltaPhi,2);
		if (tmp_jetParton_deltaR2 <= partonClosestJetsDeltaR2Vec[partonIndex][1] ) {
		  if (tmp_jetParton_deltaR2 <= partonClosestJetsDeltaR2Vec[partonIndex][0] ) {
		    partonClosestJetsDeltaR2Vec[partonIndex][1] = partonClosestJetsDeltaR2Vec[partonIndex][0];
		    partonClosestJetsVec[partonIndex][1] = partonClosestJetsVec[partonIndex][0];
		    partonClosestJetsDeltaR2Vec[partonIndex][0] = tmp_jetParton_deltaR2;
		    partonClosestJetsVec[partonIndex][0] = *offlinejet;
		  } else {
		    partonClosestJetsDeltaR2Vec[partonIndex][1] = tmp_jetParton_deltaR2;
		    partonClosestJetsVec[partonIndex][1] = *offlinejet;
		  }
		}
	      }
	      jetEt_->Fill( tmpOfflinejetEt );
	      jetEta_->Fill(tmpOfflinejetEta);
	    } // end loop over jets
	  } // end loop over partons from Z

	  std::vector< OfflineJet > ZjetVec;
	  std::vector< double >     ZjetDeltaR2Vec;
	  if (partonClosestJetsVec[0][0].et() == partonClosestJetsVec[1][0].et()) {
	    std::cout << "*************************" << std::endl;
	    std::cout << "partons matched to the same jet!" << std::endl;
	    std::cout << "BEFORE: partonClosestJetsDeltaR2Vec[0][0]: " << partonClosestJetsDeltaR2Vec[0][0] << std::endl;
	    std::cout << "BEFORE: partonClosestJetsDeltaR2Vec[0][1]: " << partonClosestJetsDeltaR2Vec[0][1] << std::endl;
	    std::cout << "BEFORE: partonClosestJetsDeltaR2Vec[1][0]: " << partonClosestJetsDeltaR2Vec[1][0] << std::endl;
	    std::cout << "BEFORE: partonClosestJetsDeltaR2Vec[1][1]: " << partonClosestJetsDeltaR2Vec[1][1] << std::endl;
	    std::cout << "*************************" << std::endl;
	    if ( partonClosestJetsDeltaR2Vec[0][0] <= partonClosestJetsDeltaR2Vec[1][0] ) {
	      ZjetVec.push_back(partonClosestJetsVec[0][0]);
	      ZjetVec.push_back(partonClosestJetsVec[1][1]);
	      ZjetDeltaR2Vec.push_back(partonClosestJetsDeltaR2Vec[0][0]);
	      ZjetDeltaR2Vec.push_back(partonClosestJetsDeltaR2Vec[1][1]);
	    } else {
	      ZjetVec.push_back(partonClosestJetsVec[0][1]);
	      ZjetVec.push_back(partonClosestJetsVec[1][0]);
	      ZjetDeltaR2Vec.push_back(partonClosestJetsDeltaR2Vec[0][1]);
	      ZjetDeltaR2Vec.push_back(partonClosestJetsDeltaR2Vec[1][0]);
	    }
	    std::cout << "AFTER: ZjetDeltaR2Vec[0]: " << ZjetDeltaR2Vec[0] << std::endl;
	    std::cout << "AFTER: ZjetDeltaR2Vec[1]: " << ZjetDeltaR2Vec[1] << std::endl;
	  } else {
	    int partonIndex = 0;
	    for ( ; partonIndex != hadronicZcandidates.size(); partonIndex++ ) {
	      ZjetVec.push_back(partonClosestJetsVec[partonIndex][0]);
	      ZjetDeltaR2Vec.push_back(partonClosestJetsDeltaR2Vec[partonIndex][0]);
	    }
	  }

	  if ( etcut == jetEtCut_ ) {
	    if ( ZjetDeltaR2Vec.size() != 0 ) {
	      if ( ZjetDeltaR2Vec.size() != 2 ) std::cout << "ZjetDeltaR2Vec.size() != 2" << std::endl;
	      for ( int partonIndex = 0; partonIndex != ZjetDeltaR2Vec.size(); partonIndex++ ) {
		//	    std::cout << "ZjetDeltaR2Vec[" << partonIndex << "]: " << ZjetDeltaR2Vec[partonIndex] << std::endl;
		jetParton_deltaR2_->Fill(ZjetDeltaR2Vec[partonIndex]);
		if (ZjetDeltaR2Vec[partonIndex] <= 1.) jetParton_deltaR2_zoom_->Fill(ZjetDeltaR2Vec[partonIndex]);
		double deltaEt  = hadronicZcandidates[partonIndex]->et()-ZjetVec[partonIndex].et();
		double deltaE   = hadronicZcandidates[partonIndex]->energy()-ZjetVec[partonIndex].e();
		double deltaEta = hadronicZcandidates[partonIndex]->eta()-ZjetVec[partonIndex].eta();
		double deltaR2  = ZjetDeltaR2Vec[partonIndex];
		jetParton_deltaEtVSdeltaR2_->Fill( deltaR2,deltaEt );
		jetParton_deltaEVSdeltaR2_->Fill(  deltaR2,deltaE  );
		jetParton_deltaEtaVSdeltaR2_->Fill(deltaR2,deltaEta);
		jetParton_deltaEtVSdeltaR2_profile_->Fill( deltaR2,deltaEt );
		jetParton_deltaEVSdeltaR2_profile_->Fill(  deltaR2,deltaE  );
		jetParton_deltaEtaVSdeltaR2_profile_->Fill(deltaR2,deltaEta);
		jetParton_deltaEta_->Fill(deltaEta);
		jetParton_deltaEt_->Fill(deltaEt);
		if (hadronicZcandidates[partonIndex]->et()!=0.)jetParton_deltaEtRes_->Fill(deltaEt/hadronicZcandidates[partonIndex]->et());
		if (ZjetDeltaR2Vec[partonIndex] <= jetPartonDeltaR2Cut_) {
		  if (hadronicZcandidates[partonIndex]->et()!=0.) jetParton_004deltaEtRes_->Fill(deltaEt/hadronicZcandidates[partonIndex]->et());
		}
	      }
	      double deltaR2_max = ZjetDeltaR2Vec[0];
	      if (ZjetDeltaR2Vec[0] < ZjetDeltaR2Vec[1]) deltaR2_max = ZjetDeltaR2Vec[1];
	      double deltaR2_mean = TMath::Sqrt(pow(ZjetDeltaR2Vec[0],2)+pow(ZjetDeltaR2Vec[1],2)); 
	      double deltaEt_mean = TMath::Sqrt(pow((ZjetVec[0].et()-hadronicZcandidates[0]->et()),2)+
						pow((ZjetVec[1].et()-hadronicZcandidates[1]->et()),2)); 
	      if (ZjetDeltaR2Vec[0]<=0.16 && ZjetDeltaR2Vec[1]<=0.16) twoJetsAboveEtCut_DeltaR2016EventsCounter_++;
	      if (ZjetDeltaR2Vec[0]<=jetPartonDeltaR2Cut_ && ZjetDeltaR2Vec[1]<=jetPartonDeltaR2Cut_) twoJetsAboveEtCut_DeltaR2016EventsCounter_++;
	      jetParton_deltaR2max_->Fill(deltaR2_max);
	      jetParton_deltaEtmeanVSdeltaR2mean_->Fill(deltaR2_mean,deltaEt_mean);
	      jetParton_deltaEtmeanVSdeltaR2mean_profile_->Fill(deltaR2_mean,deltaEt_mean);

	      Particle * hadronicZrec = 0;
	      if (ZjetDeltaR2Vec[0] != 99. && ZjetDeltaR2Vec[1] != 99.) {
		std::vector< double > partonMatchedJetScaleVec;
		for ( int partonIndex = 0; partonIndex != ZjetDeltaR2Vec.size(); partonIndex++ ) {
		  double scale = 1.;
		  if (ZjetVec[partonIndex].uncorrEt() != 0.) scale = ZjetVec[partonIndex].et()/ZjetVec[partonIndex].uncorrEt();
		  partonMatchedJetScaleVec.push_back(scale);
		}
		hadronicZrec = new Particle(ZCharge_,
					    ZjetVec[0].p4()*partonMatchedJetScaleVec[0]
					    +(ZjetVec[1].p4()*partonMatchedJetScaleVec[1]),
					    ZjetVec[0].vertex(),
					    pythiaZ_,0,true);
		double hadronicZrecMass = hadronicZrec->mass();
		double hadronicZMass    = hadronicZ->mass();
		double hadronicZrecMassResolution = -99.;
		if ( hadronicZMass != 0. ) hadronicZrecMassResolution = (hadronicZrecMass-hadronicZMass)/hadronicZMass;
		//	    std::cout << "hadronicZrecMass: " << hadronicZrecMass << std::endl;
		double hadronicZrecPt = hadronicZrec->pt();
		double hadronicZPt    = hadronicZ->pt();
		double hadronicZrecPtResolution = -99.;
		if ( hadronicZPt != 0. ) hadronicZrecPtResolution = (hadronicZrecPt-hadronicZPt)/hadronicZPt;
		double hadronicZrecCollinearity = (ZjetVec[0].p4()).pz()*partonMatchedJetScaleVec[0]*(ZjetVec[1].p4()).pz()*partonMatchedJetScaleVec[1];
		double hadronicZCollinearity    = hadronicZcandidates[0]->pz()*hadronicZcandidates[1]->pz();
		double hadronicZrecCollinearityResolution = -99.;
		if ( hadronicZCollinearity != 0. ) hadronicZrecCollinearityResolution = (hadronicZrecCollinearity-hadronicZCollinearity)/hadronicZCollinearity;
		double hadronicZDeltaR = TMath::Sqrt(pow(DeltaPhi(hadronicZcandidates[0]->phi(),hadronicZcandidates[1]->phi()),2)+pow(hadronicZcandidates[0]->eta()-hadronicZcandidates[1]->eta(),2));
		double hadronicZrecDeltaR = TMath::Sqrt(pow(DeltaPhi(ZjetVec[0].phi(),ZjetVec[1].phi()),2)+pow(ZjetVec[0].eta()-ZjetVec[1].eta(),2));
		hadronicZrecMass_->Fill(hadronicZrecMass);
		hadronicZMass_->Fill(hadronicZMass);
		hadronicZrecMassResolution_->Fill(hadronicZrecMassResolution);
		hadronicZrecPt_->Fill(hadronicZrecPt);
		hadronicZPt_->Fill(hadronicZPt);
		hadronicZrecPtResolution_->Fill(hadronicZrecPtResolution);
		hadronicZrecCollinearity_->Fill(hadronicZrecCollinearity);
		hadronicZCollinearity_->Fill(hadronicZCollinearity);
		hadronicZrecCollinearityResolution_->Fill(hadronicZrecCollinearityResolution);
		jetParton_deltaRVSdeltaR_->Fill(hadronicZDeltaR,hadronicZrecDeltaR);
		jetParton_deltaRVSdeltaR_profile_->Fill(hadronicZDeltaR,hadronicZrecDeltaR);
	      }
	      delete hadronicZrec;
	    }
	  }

	  // let's take stock:
	  // -----------------
	  if ( aboveEtCutJetVec.size() < njets_ ) continue;
	  Allevents->Fill(etcut);
	  if ( ZjetDeltaR2Vec[0]<0.25 && ZjetDeltaR2Vec[1]<0.25 ) Bothok_05->Fill(etcut);
	  if ( ZjetDeltaR2Vec[0]<0.04 && ZjetDeltaR2Vec[1]<0.04 ) Bothok_02->Fill(etcut);
	  if ( ZjetDeltaR2Vec[0]<0.25 && ZjetDeltaR2Vec[1]<0.25 &&
	       fabs(hadronicZcandidates[0]->et()-ZjetVec[0].et())/hadronicZcandidates[0]->et() < 0.3 &&
	       fabs(hadronicZcandidates[1]->et()-ZjetVec[1].et())/hadronicZcandidates[1]->et() < 0.3
	       ) Bothok_05_30pc->Fill(etcut);
	  if ( ZjetDeltaR2Vec[0]<0.04 && ZjetDeltaR2Vec[1]<0.04  &&
	       fabs(hadronicZcandidates[0]->et()-ZjetVec[0].et())/hadronicZcandidates[0]->et() < 0.3 &&
	       fabs(hadronicZcandidates[1]->et()-ZjetVec[1].et())/hadronicZcandidates[1]->et() < 0.3
	       ) Bothok_02_30pc->Fill(etcut);
	  if ( ZjetDeltaR2Vec[0]<0.25 && ZjetDeltaR2Vec[1]<0.25 &&
	       fabs(hadronicZcandidates[0]->et()-ZjetVec[0].et())/hadronicZcandidates[0]->et() < 0.2 &&
	       fabs(hadronicZcandidates[1]->et()-ZjetVec[1].et())/hadronicZcandidates[1]->et() < 0.2
	       ) Bothok_05_20pc->Fill(etcut);
	  if ( ZjetDeltaR2Vec[0]<0.04 && ZjetDeltaR2Vec[1]<0.04  &&
	       fabs(hadronicZcandidates[0]->et()-ZjetVec[0].et())/hadronicZcandidates[0]->et() < 0.2 &&
	       fabs(hadronicZcandidates[1]->et()-ZjetVec[1].et())/hadronicZcandidates[1]->et() < 0.2
	       ) Bothok_02_20pc->Fill(etcut);
	  if ( ZjetDeltaR2Vec[0]<0.25 && ZjetDeltaR2Vec[1]<0.25 &&
	       fabs(hadronicZcandidates[0]->et()-ZjetVec[0].et())/hadronicZcandidates[0]->et() < 0.1 &&
	       fabs(hadronicZcandidates[1]->et()-ZjetVec[1].et())/hadronicZcandidates[1]->et() < 0.1
	       ) Bothok_05_10pc->Fill(etcut);
	  if ( ZjetDeltaR2Vec[0]<0.04 && ZjetDeltaR2Vec[1]<0.04  &&
	       fabs(hadronicZcandidates[0]->et()-ZjetVec[0].et())/hadronicZcandidates[0]->et() < 0.1 &&
	       fabs(hadronicZcandidates[1]->et()-ZjetVec[1].et())/hadronicZcandidates[1]->et() < 0.1
	       ) Bothok_02_10pc->Fill(etcut);
	} // end loop on iet
	
      }
      delete hadronicZ;
      delete leptonicZ;
    }
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
VBFHZZanalyzer::beginJob(const edm::EventSetup&)
{

  vecTypeIndexPair.push_back(std::pair<int,int>(pythiab_, 0));
  vecTypeIndexPair.push_back(std::pair<int,int>(pythiae_, 1));
  vecTypeIndexPair.push_back(std::pair<int,int>(pythiamu_,2));

  // File for output histograms
  // --------------------------
  OutputFile = new TFile((conf_.getUntrackedParameter<std::string>("OutputName")).c_str() ,
			 "RECREATE","HiggsZZOutput");
  // The file must be opened first, 
  // so that becomes the default position for all the histograms
  OutputFile->cd();
  // White background for the canvases
  gROOT->SetStyle("Plain");
  eventsNumber_ = new TH1D("eventsNumber","total number of events",1,0.,1.);
  goodLeptonicEventsNumber_                  = new TH1D("goodLeptonicEventsNumber",                 "number of events w/ good leptonic part",                                         1,0.,1.);
  twoJetsAboveEtCutEventsNumber_             = new TH1D("twoJetsAboveEtCutEventsNumber",            "number of events w/ >= 2 jets above E_{T} cut",                                  1,0.,1.);
  twoJetsAboveEtCut_DeltaR2016EventsNumber_  = new TH1D("twoJetsAboveEtCut_DeltaR2016EventsNumber", "number of events w/ >= 2 jets above E_{T} cut and #DeltaR^{2}_{qj}<=0.16",       1,0.,1.);
  twoJetsAboveEtCut_DeltaR2CutEventsNumber_  = new TH1D("twoJetsAboveEtCut_DeltaR2CutEventsNumber", "number of events w/ >= 2 jets above E_{T} cut and #DeltaR^{2}_{qj}<=#DeltaR cut",1,0.,1.);
  fourJetsAboveEtCutEventsNumber_            = new TH1D("fourJetsAboveEtCutEventsNumber",           "number of events w/ >= 4 jets above E_{T} cut",                                  1,0.,1.);


  // histo from HEPG information
  // ---------------------------
  ZpartonEta_ = new TH1D("ZpartonEta","#eta of partons from Z", nbin_, -8.,   8.);
  ZpartonPt_  = new TH1D("ZpartonPt", "p_{T} of partons from Z",nbin_,  0., 500.);
  ZpartonE_   = new TH1D("ZpartonE",  "E of partons from Z",    nbin_,-10.,7000.);
  ZpartonsDeltaEta_     = new TH1D("ZpartonsDeltaEta",    "#Delta#eta between partons from Z",nbin_,    0.,   10.);
  ZpartonsMass_         = new TH1D("ZpartonsMass",        "invariant mass of partons from Z", nbin_,    0., 10000.);
  ZpartonsCollinearity_ = new TH1D("ZpartonsCollinearity","collinearity of partons from Z",   nbin_,-7000., 7000.);

  TagPartonEta_ = new TH1D("TagPartonEta","#eta of tag partons", nbin_, -8.,   8.);
  TagPartonPt_  = new TH1D("TagPartonPt", "p_{T} of tag partons",nbin_,  0., 500.);
  TagPartonE_   = new TH1D("TagPartonE",  "E of tag partons",    nbin_,-10.,7000.);
  TagPartonsDeltaEta_     = new TH1D("TagPartonsDeltaEta",    "#Delta#eta between tag partons",nbin_,    0.,   10.);
  TagPartonsMass_         = new TH1D("TagPartonsMass",        "invariant mass of tag partons", nbin_,    0.,10000.);
  TagPartonsCollinearity_ = new TH1D("TagPartonsCollinearity","collinearity of tag partons",   nbin_,-7000., 7000.);

  hadronicZMass_         = new TH1D("hadronicZMass",        "Z mass from partons",                nbin_,    0., 150.);
  hadronicZPt_           = new TH1D("hadronicZPt",          "Z p_{T} from partons",               nbin_,    0., 500.);
  hadronicZCollinearity_ = new TH1D("hadronicZCollinearity","collinearity between partons from Z",nbin_,-6000.,6000.);

  jetNumber_        = new TH1D("jetNumber",  "number of jets",         nbin_,   0.,100.);
  jetUncorrEt_      = new TH1D("jetUncorrEt","offline jet uncor E_{T}",nbin_,   0.,500.);
  jetEt_            = new TH1D("jetEt",      "offline jet E_{T}",      nbin_,   0.,500.);
  jetPhi_           = new TH1D("jetPhi",     "offline jet #Phi",       nbin_, -TMath::Pi(),TMath::Pi() );
  jetEta_           = new TH1D("jetEta",     "offline jet #eta",       nbin_,  -8., 8.);

  jetEta_aboveEtCut_      = new TH1D("jetEta_aboveEtCut",      "#eta of offline jet w/ E_{T} >= E_{T}^{CUT}",      nbin_,-8.,  8.);
  jetNumber_aboveEtCut_   = new TH1D("jetNumber_aboveEtCut" ,  "number of jets w/ E_{T} >= E_{T}^{CUT}",           nbin_, 0.,100.);
  jetsDeltaR2_aboveEtCut_ = new TH1D("jetsDeltaR2_aboveEtCut ","#DeltaR^{2} between jets w/ E_{T} >= E_{T}^{CUT}", nbin_,0.,  10.);

  muonJetMatchedNumber_ = new TH1D("muonJetMatchedNumber","number of jets w/ at least 1 muon in a distance #DeltaR<0.5", 10,0.,10.0);
  muonsMatchedNumber_   = new TH1D("muonsMatchedNumber",  "fraction of muons in a distance #DeltaR<0.5 to a jet",       110,0., 1.1);

  jetParton_deltaR2_      = new TH1D("jetParton_deltaR2",     "#DeltaR^{2} between jet and partons from Z",    nbin_,0.,4.0);
  jetParton_deltaR2_zoom_ = new TH1D("jetParton_deltaR2_zoom","#DeltaR^{2} between jet and partons from Z",    nbin_,0.,0.1);
  jetParton_deltaR2max_   = new TH1D("jetParton_deltaR2max",  "max #DeltaR^{2} between jet and partons from Z",nbin_,0.,4.0);

  jetParton_deltaEtVSdeltaR2_         = new TH2D("jetParton_deltaEtVSdeltaR2",        "#DeltaE_{T} VS #DeltaR^{2} between partons and matched jets",          nbin_,0.,0.04,nbin_,-50.,50.);
  jetParton_deltaEtmeanVSdeltaR2mean_ = new TH2D("jetParton_deltaEtmeanVSdeltaR2mean","mean #DeltaE_{T} VS mean #DeltaR^{2} between partons and matched jets",nbin_,0.,0.04,nbin_,  0.,50.);
  jetParton_deltaEVSdeltaR2_   = new TH2D("jetParton_deltaEVSdeltaR2",  "#DeltaE VS #DeltaR^{2} between partons and matched jets",    nbin_,0.,0.04,nbin_,-50.,50.);
  jetParton_deltaEtaVSdeltaR2_ = new TH2D("jetParton_deltaEtaVSdeltaR2","#Delta#eta VS #DeltaR^{2} between partons and matched jets", nbin_,0.,0.04,nbin_, -5., 5.);
  jetParton_deltaEtVSdeltaR2_profile_         = new TProfile("jetParton_deltaEtVSdeltaR2_profile",        "#DeltaE_{T} VS #DeltaR^{2} between partons and matched jets",          nbin_,0.,0.04,-50.,50.);
  jetParton_deltaEtmeanVSdeltaR2mean_profile_ = new TProfile("jetParton_deltaEtmeanVSdeltaR2mean_profile","mean #DeltaE_{T} VS mean #DeltaR^{2} between partons and matched jets",nbin_,0.,0.04,  0.,50.);
  jetParton_deltaEVSdeltaR2_profile_   = new TProfile("jetParton_deltaEVSdeltaR2_profile",  "#DeltaE VS #DeltaR^{2} between partons and matched jets",    nbin_,0.,0.04,-50.,50.);
  jetParton_deltaEtaVSdeltaR2_profile_ = new TProfile("jetParton_deltaEtaVSdeltaR2_profile","#Delta#eta VS #DeltaR^{2} between partons and matched jets", nbin_,0.,0.04, -5., 5.);
  jetParton_deltaEta_ = new TH1D("jetParton_deltaEta","#Delta#eta between partons and matched jets", nbin_, -5., 5.);
  jetParton_deltaEt_  = new TH1D("jetParton_deltaEt", "#DeltaE_{T} between partons and matched jets",nbin_,-50.,50.);
  jetParton_deltaEtRes_  = new TH1D("jetParton_deltaEtRes", "E_{T} resolution between partons and matched jets",nbin_,-1.,1.);
  jetParton_004deltaEtRes_ = new TH1D("jetParton_004deltaEtRes", "E_{T} resolution between partons and matched jets (#DeltaR<=0.2)",nbin_,-1.,1.);

  jetParton_deltaRVSdeltaR_ = new TH2D("jetParton_deltaRVSdeltaR","#DeltaR between the 2jets VS #DeltaR between the 2 partons",nbin_,0.,5.,nbin_,0.,5.);
  jetParton_deltaRVSdeltaR_profile_ = new TProfile("jetParton_deltaRVSdeltaR_profile","#DeltaR between the 2jets VS #DeltaR between the 2 partons",nbin_,0.,5.,0.,5.);

  hadronicZrecMass_           = new TH1D("hadronicZrecMass",          "reconstructed Z mass from jets matched to Z partons",           nbin_, 0.,150.);
  hadronicZrecMassResolution_ = new TH1D("hadronicZrecMassResolution","reconstructed Z mass resolution from jets matched to Z partons",nbin_,-1.,  1.);
  hadronicZrecPt_           = new TH1D("hadronicZrecPt",          "reconstructed Z p_{T} from jets matched to Z partons",           nbin_, 0.,500.);
  hadronicZrecPtResolution_ = new TH1D("hadronicZrecPtResolution","reconstructed Z p_{T} resolution from jets matched to Z partons",nbin_,-1.,  1.);
  hadronicZrecCollinearity_           = new TH1D("hadronicZrecCollinearity",          "reconstructed collinearity between jets matched to Z partons",        nbin_, -6000.,6000.);
  hadronicZrecCollinearityResolution_ = new TH1D("hadronicZrecCollinearityResolution","reconstructed collinearity resolution from jets matched to Z partons",nbin_,    -1.,   1.);

  muonNumber_     = new TH1D("muonNumber","number of muons",nbin_,   0.,         10.);
  muonEta_        = new TH1D("muonEta",   "muon #eta",      nbin_, -10.,         10.);
  muonPhi_        = new TH1D("muonPhi",   "muon #Phi",      nbin_, -TMath::Pi(),TMath::Pi() );
  muonPt_         = new TH1D("muonPt",    "muon p_{T}",     nbin_,   0.,       1000.);

  electronNumber_ = new TH1D("electronNumber","number of electrons",nbin_,   0.,         10.);
  electronEt_     = new TH1D("electronEt",    "electron E_{T}",     nbin_,   0.,       1000.);
  electronEta_    = new TH1D("electronEta",   "electron #eta",      nbin_, -10.,         10.);
  electronPhi_    = new TH1D("electronPhi",   "electron #Phi",      nbin_, -TMath::Pi(),TMath::Pi() );
  electronPt_     = new TH1D("electronPt",    "electron p_{T}",     nbin_,   0.,       1000.);

  // histo with fractions
  // ---------------------------
  Bothok_05 = new TH1D ("Bothok_05", "Evt frac. with both tag jets found within 0.5 rads", 41, 9.5, 50.5 );
  Bothok_02 = new TH1D ("Bothok_02", "Evt frac. with both tag jets found within 0.2 rads", 41, 9.5, 50.5 );
  Bothok_05_30pc = new TH1D ("Bothok_05_30pc", "Evt frac. with both tag jets found within 0.5 rads and <30% Et res.", 41, 9.5, 50.5 );
  Bothok_02_30pc = new TH1D ("Bothok_02_30pc", "Evt frac. with both tag jets found within 0.2 rads and <30% Et res.", 41, 9.5, 50.5 );
  Bothok_05_20pc = new TH1D ("Bothok_05_20pc", "Evt frac. with both tag jets found within 0.5 rads and <20% Et res.", 41, 9.5, 50.5 );
  Bothok_02_20pc = new TH1D ("Bothok_02_20pc", "Evt frac. with both tag jets found within 0.2 rads and <20% Et res.", 41, 9.5, 50.5 );
  Bothok_05_10pc = new TH1D ("Bothok_05_10pc", "Evt frac. with both tag jets found within 0.5 rads and <10% Et res.", 41, 9.5, 50.5 );
  Bothok_02_10pc = new TH1D ("Bothok_02_10pc", "Evt frac. with both tag jets found within 0.2 rads and <10% Et res.", 41, 9.5, 50.5 );
  Allevents = new TH1D ("Allevents", "All Evts.", 41, 9.5, 50.5 );
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
VBFHZZanalyzer::endJob() {

  OutputFile->cd();

  // Fill histograms
  // ---------------
  eventsNumber_->SetBinContent(1,eventcounter_);
  eventsNumber_->Write();

  goodLeptonicEventsNumber_->SetBinContent(1,goodLeptonicEventCounter_);
  goodLeptonicEventsNumber_->Write();

  twoJetsAboveEtCutEventsNumber_->SetBinContent(1,twoJetsAboveEtCutEventsCounter_);
  twoJetsAboveEtCutEventsNumber_->Write();
  twoJetsAboveEtCut_DeltaR2016EventsNumber_->SetBinContent(1,twoJetsAboveEtCut_DeltaR2016EventsCounter_);
  twoJetsAboveEtCut_DeltaR2016EventsNumber_->Write();
  twoJetsAboveEtCut_DeltaR2CutEventsNumber_->SetBinContent(1,twoJetsAboveEtCut_DeltaR2CutEventsCounter_);
  twoJetsAboveEtCut_DeltaR2CutEventsNumber_->Write();
  fourJetsAboveEtCutEventsNumber_->SetBinContent(1,fourJetsAboveEtCutEventsCounter_);
  fourJetsAboveEtCutEventsNumber_->Write();

  ZpartonEta_->Write();
  ZpartonPt_->Write();
  ZpartonE_->Write();
  ZpartonsDeltaEta_->Write();
  ZpartonsMass_->Write();
  ZpartonsCollinearity_->Write();

  TagPartonEta_->Write();
  TagPartonPt_->Write();
  TagPartonE_->Write();
  TagPartonsDeltaEta_->Write();
  TagPartonsMass_->Write();
  TagPartonsCollinearity_->Write();

  hadronicZMass_->Write();
  hadronicZPt_->Write();
  hadronicZCollinearity_->Write();

  jetNumber_->Write();
  jetUncorrEt_->Write();
  jetEt_->Write();
  jetPhi_->Write();
  jetEta_->Write();

  muonJetMatchedNumber_->Write();
  muonsMatchedNumber_->Write();

  jetEta_aboveEtCut_->Write();
  jetNumber_aboveEtCut_->Write();
  jetsDeltaR2_aboveEtCut_->Write();

  jetParton_deltaR2_->Write();
  jetParton_deltaR2_zoom_->Write();
  jetParton_deltaR2max_->Write();

  jetParton_deltaEtVSdeltaR2_->Write();
  jetParton_deltaEtmeanVSdeltaR2mean_->Write();
  jetParton_deltaEVSdeltaR2_->Write();
  jetParton_deltaEtaVSdeltaR2_->Write();
  jetParton_deltaEtVSdeltaR2_profile_->Write();
  jetParton_deltaEtmeanVSdeltaR2mean_profile_->Write();
  jetParton_deltaEVSdeltaR2_profile_->Write();
  jetParton_deltaEtaVSdeltaR2_profile_->Write();
  jetParton_deltaEta_->Write();
  jetParton_deltaEt_->Write();
  jetParton_deltaEtRes_->Write();
  jetParton_004deltaEtRes_->Write();
  jetParton_deltaRVSdeltaR_->Write();
  jetParton_deltaRVSdeltaR_profile_->Write();

  hadronicZrecMass_->Write();
  hadronicZrecMassResolution_->Write();
  hadronicZrecPt_->Write();
  hadronicZrecPtResolution_->Write();
  hadronicZrecCollinearity_->Write();
  hadronicZrecCollinearityResolution_->Write();

  muonNumber_->Write();
  muonEta_->Write();       
  muonPhi_->Write();       
  muonPt_->Write();        

  electronNumber_->Write();
  electronEt_->Write();	    
  electronEta_->Write();	    
  electronPhi_->Write();	    
  electronPt_->Write();	    

  Bothok_05->Write();
  Bothok_02->Write();
  Bothok_05_30pc->Write();
  Bothok_02_30pc->Write();
  Bothok_05_20pc->Write();
  Bothok_02_20pc->Write();
  Bothok_05_10pc->Write();
  Bothok_02_10pc->Write();
  Allevents->Write();
}
//define this as a plug-in
DEFINE_FWK_MODULE(VBFHZZanalyzer);

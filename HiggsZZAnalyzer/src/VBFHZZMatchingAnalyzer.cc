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

#include "AnalysisExamples/HiggsZZAnalyzer/interface/VBFHZZMatchingAnalyzer.h"

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


VBFHZZMatchingAnalyzer::VBFHZZMatchingAnalyzer(const edm::ParameterSet& iConfig) :
  conf_( iConfig ),
  offlineJetLabel_(    iConfig.getUntrackedParameter<edm::InputTag>( "OfflineJets"    ) ),
  globalMuonLabel_(    iConfig.getUntrackedParameter<edm::InputTag>( "GlobalMuons"    ) ),
  simpleElectronLabel_(iConfig.getUntrackedParameter<edm::InputTag>( "SimpleElectrons") ),
  MCParticleLabel_(    iConfig.getUntrackedParameter<std::string>(   "MCParticles"    ) ),
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
  twoJetsAboveEtCutEventsCounter_            = 0;
  twoJetsAboveEtCut_DeltaR2016EventsCounter_ = 0;
  twoJetsAboveEtCut_DeltaR2CutEventsCounter_ = 0;
  fourJetsAboveEtCutEventsCounter_           = 0;

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


VBFHZZMatchingAnalyzer::~VBFHZZMatchingAnalyzer()
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
VBFHZZMatchingAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
	leptonicZcandidates.push_back(Z2Candidates[0]);
	leptonicZcandidates.push_back(Z2Candidates[1]);
      }
      else if ( fabs(Z2particlesId) <= pythiat_ && 
		(fabs(Z1particlesId) == pythiae_ || fabs(Z1particlesId) == pythiamu_ ) ) {
	hadronicZ = Z2particle;
	leptonicZ = Z1particle;
	hadronicZcandidates.push_back(Z2Candidates[0]);
	hadronicZcandidates.push_back(Z2Candidates[1]);
	leptonicZcandidates.push_back(Z1Candidates[0]);
	leptonicZcandidates.push_back(Z1Candidates[1]);
      }
    }

    //      if (hadronicZcandidates.size() == 0 || leptonicZcandidates.size() == 0) 
    //	std::cout << "Z1CandidatesID: " << Z1Candidates[0]->pdgId() 
    //		  << " <--> Z2CandidatesID: " << Z2Candidates[0]->pdgId() << std::endl;

    ////////////////////////////////////// PLOT HEPG DISTRIBUTION ////////////////////////////
    std::vector< const Candidate* > partonsCandidates;
    if ( hadronicZcandidates.size() != 0 ) {
      for ( std::vector< const Candidate* >::const_iterator hadronicZcandidates_itr = hadronicZcandidates.begin();
	    hadronicZcandidates_itr != hadronicZcandidates.end(); ++hadronicZcandidates_itr) {
	double tmp_ZpartonEta = (*hadronicZcandidates_itr)->eta();
	double tmp_ZpartonPt  = (*hadronicZcandidates_itr)->pt();
	double tmp_ZpartonE   = (*hadronicZcandidates_itr)->energy();
	ZpartonEta_->Fill(tmp_ZpartonEta);
	ZpartonPt_->Fill( tmp_ZpartonPt );
	ZpartonE_->Fill(  tmp_ZpartonE  );
	partonsCandidates.push_back(*hadronicZcandidates_itr);
      }
      double tmp_ZpartonsDeltaEta     = fabs(hadronicZcandidates[0]->eta()-hadronicZcandidates[1]->eta());
      double tmp_ZpartonsDeltaPhi     = DeltaPhi(hadronicZcandidates[0]->phi(),hadronicZcandidates[1]->phi());
      double tmp_ZpartonsDeltaR       = TMath::Sqrt(pow(tmp_ZpartonsDeltaEta,2)+pow(tmp_ZpartonsDeltaPhi,2));
      double tmp_ZpartonsMass         = (hadronicZcandidates[0]->p4()+(hadronicZcandidates[1]->p4())).mass();
      double tmp_ZpartonsCollinearity = hadronicZcandidates[0]->pz()*hadronicZcandidates[1]->pz();
      ZpartonsDeltaEta_->Fill(    tmp_ZpartonsDeltaEta    );
      ZpartonsDeltaR_->Fill(      tmp_ZpartonsDeltaR      );
      ZpartonsMass_->Fill(        tmp_ZpartonsMass        );
      ZpartonsCollinearity_->Fill(tmp_ZpartonsCollinearity);
      if ( tmp_ZpartonsDeltaR <= 0.8 ) {
	OfflineJetCollection::const_iterator offlinejet = offlineJets->begin();
	for ( ; offlinejet != offlineJets->end(); ++offlinejet ) {
	  double tmp_jetEt  = offlinejet->et();
	  double tmp_jetEta = offlinejet->eta();
	  double tmp_jetPhi = offlinejet->phi();
	  double deltaEta = tmp_jetEta-hadronicZcandidates[0]->eta();
	  double deltaPhi = DeltaPhi(tmp_jetPhi,hadronicZcandidates[0]->phi());
	  double deltaR2  = pow(deltaPhi,2)+pow(deltaEta,2);
	  double etRes = fabs((tmp_jetEt-(hadronicZcandidates[0]->p4()+(hadronicZcandidates[1]->p4())).pt())/(hadronicZcandidates[0]->p4()+(hadronicZcandidates[1]->p4())).pt());

	  if ( deltaR2 <= 0.7 ) {
	    ZpartonsDeltaRVSjetMass_->Fill(tmp_ZpartonsDeltaR,offlinejet->p4().M());
	    ZpartonsDeltaRVSjetMass_profile_->Fill(tmp_ZpartonsDeltaR,offlinejet->p4().M());
	  }
	  if ( etRes <= 0.9 ) {
	    ZpartonsEtResVSjetMass_->Fill(etRes,offlinejet->p4().M());
	    ZpartonsEtResVSjetMass_profile_->Fill(etRes,offlinejet->p4().M());
	  }
	  if ( etRes <= 0.9 || deltaR2 <= 0.7 ) {
	    ZpartonsDeltaRVSetRes_->Fill(deltaR2,etRes);
	    ZpartonsDeltaRVSetRes_profile_->Fill(deltaR2,etRes);
	  }
	}
      }
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
	partonsCandidates.push_back(*TagCandidates_itr);
      }
      double tmp_TagPartonsDeltaEta     = fabs(TagCandidates[0]->eta()-TagCandidates[1]->eta());
      double tmp_TagPartonsDeltaPhi     = DeltaPhi(TagCandidates[0]->phi(),TagCandidates[1]->phi());
      double tmp_TagPartonsDeltaR       = TMath::Sqrt(pow(tmp_TagPartonsDeltaEta,2)+pow(tmp_TagPartonsDeltaPhi,2));
      double tmp_TagPartonsMass         = (TagCandidates[0]->p4()+(TagCandidates[1]->p4())).mass();
      double tmp_TagPartonsCollinearity = TagCandidates[0]->pz()*TagCandidates[1]->pz();
      TagPartonsDeltaEta_->Fill(    tmp_TagPartonsDeltaEta    );
      TagPartonsDeltaR_->Fill(      tmp_TagPartonsDeltaR      );
      TagPartonsMass_->Fill(        tmp_TagPartonsMass        );
      TagPartonsCollinearity_->Fill(tmp_TagPartonsCollinearity);
    }

    //      if (partonsCandidates.size()!=4) std::cout << "partonsCandidates.size: "    << partonsCandidates.size() 
    //						 << " hadronicZcandidates.size: " << hadronicZcandidates.size() 
    //						 << " TagCandidates.size: "       << TagCandidates.size()       << std::endl;


    /////////////////////////////////// QUICK LOOP over LEPTONS ////////////////////////////
    // Global muons
    // ------------
    std::vector< GlobalMuon > goodMuonVec;
    std::vector< std::pair<GlobalMuon,GlobalMuon> > goodMuonsPairVec;
    muonNumber_->Fill(globalMuons->size());
    GlobalMuonCollection::const_iterator globalMuon_itr = globalMuons->begin();
    for ( ; globalMuon_itr != globalMuons->end(); ++globalMuon_itr ) {
      double tmpMuonPt  = globalMuon_itr->pt();
      double tmpMuonEta = globalMuon_itr->eta();
      double tmpMuonPhi = globalMuon_itr->phi();
      if ( tmpMuonPt >= leptonPtCut_ ) {
	goodMuonVec.push_back(*globalMuon_itr);
	if ( globalMuon_itr != (globalMuons->end()-1) ) {
	  GlobalMuonCollection::const_iterator globalMuon_itr2 = globalMuon_itr+1;
	  for ( ; globalMuon_itr2 != globalMuons->end(); ++globalMuon_itr2 ){
	    if ( globalMuon_itr->charge()*globalMuon_itr2->charge() < 0 ) {
	      std::pair< GlobalMuon, GlobalMuon > tmpMuonsPair(*globalMuon_itr,*globalMuon_itr2);
	      goodMuonsPairVec.push_back(tmpMuonsPair);
	    }
	  }
	}
      }
      muonEta_->Fill(tmpMuonEta);
      muonPhi_->Fill(tmpMuonPhi);
      muonPt_->Fill( tmpMuonPt );
    }
    //    if ( goodMuonsPairVec.size() >= 1 ) {
    //      std::cout << "goodMuonsPairVec: " << goodMuonsPairVec.size() << std::endl;
    //      std::cout << "goodMuonVec: " << goodMuonVec.size() << std::endl;
    //    }

    // SimpleElectrons
    // ---------------
    std::vector< SimpleElectron > goodElectronVec;
    std::vector< std::pair<SimpleElectron,SimpleElectron> > goodElectronsPairVec;
    electronNumber_->Fill(simpleElectrons->size());
    SimpleElectronCollection::const_iterator simpleElectron_itr = simpleElectrons->begin();
    for ( ; simpleElectron_itr != simpleElectrons->end(); ++simpleElectron_itr ) {
      double tmpElectronEt  = simpleElectron_itr->et();
      double tmpElectronEta = simpleElectron_itr->eta();
      double tmpElectronPhi = simpleElectron_itr->phi();
      double tmpElectronPt  = simpleElectron_itr->pt();
      if ( tmpElectronPt >= leptonPtCut_ ) {
	goodElectronVec.push_back(*simpleElectron_itr);
	if ( simpleElectron_itr != (simpleElectrons->end()-1) ) {
	  SimpleElectronCollection::const_iterator simpleElectron_itr2 = simpleElectron_itr+1;
	  for ( ; simpleElectron_itr2 != simpleElectrons->end(); ++simpleElectron_itr2 ){
	    if ( simpleElectron_itr->charge()*simpleElectron_itr2->charge() < 0 ) {
	      std::pair< SimpleElectron, SimpleElectron> tmpElectronsPair(*simpleElectron_itr,*simpleElectron_itr2);
	      goodElectronsPairVec.push_back(tmpElectronsPair);
	    }
	  }
	}
      }
      electronEt_->Fill( tmpElectronEt );
      electronEta_->Fill(tmpElectronEta);
      electronPhi_->Fill(tmpElectronPhi);
      electronPt_->Fill( tmpElectronPt );
    }
    //    if ( goodElectronsPairVec.size() >= 1 ) {
    //      std::cout << "goodElectronsPairVec: " << goodElectronsPairVec.size() << std::endl;
    //      std::cout << "goodElectronVec: " << goodElectronVec.size() << std::endl;
    //    }
  
    std::vector<math::XYZTLorentzVector> ZlepVec;
    if ( goodMuonVec.size() >= nZleptons_ || goodElectronVec.size() >= nZleptons_ ) {

      //    /////////////////////////////////// LEPTON MATCHING ////////////////////////////
      if ( leptonicZcandidates.size() != 0 ) {
	int leptonId = fabs(leptonicZcandidates[0]->pdgId());
	if ( leptonId == pythiamu_ || leptonId == pythiae_ ) {
	  std::vector<int> ZlepChargeVec;

	  std::vector< const Candidate* >::const_iterator leptonicZcandidates_itr = leptonicZcandidates.begin();
	  for ( ; leptonicZcandidates_itr != leptonicZcandidates.end(); ++leptonicZcandidates_itr ) {
	    ZlepVec.push_back(null_XYZTLorentzVector_);
	    ZlepChargeVec.push_back(99);
	  }
	  leptonicZcandidates_itr = leptonicZcandidates.begin();
	  unsigned int leptonicZcandidateIndex = 0;

	  for ( ; leptonicZcandidates_itr != leptonicZcandidates.end(); ++leptonicZcandidates_itr, leptonicZcandidateIndex++ ) {
	    double tmp_leptonicZcandidateEta = (*leptonicZcandidates_itr)->eta();
	    double tmp_leptonicZcandidatePhi = (*leptonicZcandidates_itr)->phi();
	    double leptonDeltaR2 = 99.;
	    if ( leptonId == pythiamu_) {
	      std::vector< GlobalMuon >::const_iterator goodMuonVec_itr = goodMuonVec.begin();
	      for ( ; goodMuonVec_itr != goodMuonVec.end(); ++goodMuonVec_itr ) {
		double tmp_leptonEta    = goodMuonVec_itr->eta();
		double tmp_leptonPhi    = goodMuonVec_itr->phi();
		double tmp_leptonCharge = goodMuonVec_itr->charge();
		double tmp_deltaEta = tmp_leptonicZcandidateEta - tmp_leptonEta;
		double tmp_deltaPhi = DeltaPhi(tmp_leptonicZcandidatePhi,tmp_leptonPhi);
		double tmp_deltaR2 = pow(tmp_deltaEta,2)+pow(tmp_deltaPhi,2);
		if ( tmp_deltaR2 <= leptonDeltaR2 ) {
		  if ( ( leptonicZcandidateIndex == goodMuonVec.size()-1 ) &&
		       ( ZlepChargeVec[leptonicZcandidateIndex]*tmp_leptonCharge < 0 ) ) {
		    ZlepVec[leptonicZcandidateIndex] = goodMuonVec_itr->p4();
		    leptonDeltaR2 = tmp_deltaR2;
		    ZlepChargeVec[leptonicZcandidateIndex] = tmp_leptonCharge;
		  }
		}
	      }
	    } else if (leptonId == pythiae_ ) {
	      std::vector< SimpleElectron >::const_iterator goodElectronVec_itr = goodElectronVec.begin();
	      for ( ; goodElectronVec_itr != goodElectronVec.end(); ++goodElectronVec_itr ) {
		double tmp_leptonEta    = goodElectronVec_itr->eta();
		double tmp_leptonPhi    = goodElectronVec_itr->phi();
		double tmp_leptonCharge = goodElectronVec_itr->charge();
		double tmp_deltaEta = tmp_leptonicZcandidateEta - tmp_leptonEta;
		double tmp_deltaPhi = DeltaPhi(tmp_leptonicZcandidatePhi,tmp_leptonPhi);
		double tmp_deltaR2 = pow(tmp_deltaEta,2)+pow(tmp_deltaPhi,2);
		if ( tmp_deltaR2 <= leptonDeltaR2 ) {
		  if ( ( leptonicZcandidateIndex == goodElectronVec.size()-1 ) &&
		       ( ZlepChargeVec[leptonicZcandidateIndex]*tmp_leptonCharge < 0 ) ) {
		    ZlepVec[leptonicZcandidateIndex] = goodElectronVec_itr->p4();
		    leptonDeltaR2 = tmp_deltaR2;
		    ZlepChargeVec[leptonicZcandidateIndex] = tmp_leptonCharge;
		  }
		}
	      }
	    }
	  }
	  //	  std::cout << "ZlepVec: " << ZlepVec.size() << std::endl;
	}
      }
    }
    if ( ZlepVec.size() >= nZleptons_ ) {
      goodLeptonicEventCounter_++;

      if ( partonsCandidates.size() >= njets_ ) {

	//    /////////////////////////////////// JET-PARTON MATCHING ////////////////////////////
	// Calorimeter jets
	// ----------------
	jetNumber_->Fill(offlineJets->size());

	std::vector< OfflineJet > aboveEtCutJetVec;
	OfflineJetCollection::const_iterator offlinejet = offlineJets->begin();
	for ( ; offlinejet != offlineJets->end(); ++offlinejet ) {
	  double tmpOfflinejetEt  = offlinejet->et();
	  double tmpOfflinejetEta = offlinejet->eta();
	  //	  if ( tmpOfflinejetEt >= jetEtCut_ ) {
	    aboveEtCutJetVec.push_back(*offlinejet);
	    jetEta_aboveEtCut_->Fill(tmpOfflinejetEta);
	    //	  }
	}
	jetNumber_aboveEtCut_->Fill(aboveEtCutJetVec.size());
	if ( aboveEtCutJetVec.size() >= nZjets_ ) twoJetsAboveEtCutEventsCounter_++;
	if ( aboveEtCutJetVec.size() >= njets_ )  fourJetsAboveEtCutEventsCounter_++;
	if ( aboveEtCutJetVec.size() >= nZjets_ ) { 
	  if ( aboveEtCutJetVec.size() >= njets_ ) {
	
	    std::vector< OfflineJet > ZjetVec;
	    std::vector< double >     ZjetDeltaR2Vec;
	    std::vector<int>          ZjetInd;
	    unsigned int hadronicZcandidatesNumber = hadronicZcandidates.size();
	    for ( unsigned int index = 0; index < hadronicZcandidatesNumber; index++) {
	      //	      ZjetVec.push_back(null_offlinejet_);
	      ZjetDeltaR2Vec.push_back(99.);
	      ///	      ZjetInd.push_back(-1);
	    }
	    std::vector<OfflineJet> TjetVec;
	    std::vector< double >   TjetDeltaR2Vec;
	    std::vector<int>        TjetInd;
	    unsigned int TagCandidatesNumber = TagCandidates.size();
	    for ( unsigned int index = 0; index < TagCandidatesNumber; index++) {
	      //	      TjetVec.push_back(null_offlinejet_);
	      TjetDeltaR2Vec.push_back(99.);
	      //	      TjetInd.push_back(-1);
	    }
	    

	    std::vector<std::vector<double> > closestJetDeltaR2Vec;
	    std::vector<std::vector<int> >    closestJetIndexVec;
	    unsigned int partonsCandidatesNumber = partonsCandidates.size();
	    for ( unsigned index = 0; index < partonsCandidatesNumber; index++ ) {
	      std::vector<double> null_closestJetDeltaR2Vec;
	      std::vector<int>    null_closestJetIndexVec;
	      for ( unsigned index = 0; index < partonsCandidatesNumber; index++ ) {
		null_closestJetDeltaR2Vec.push_back(99.);
		null_closestJetIndexVec.push_back(-1);
	      }
	      closestJetDeltaR2Vec.push_back(null_closestJetDeltaR2Vec);
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
	      std::vector< OfflineJet >::const_iterator aboveEtCutJetVec_itr = aboveEtCutJetVec.begin();
	      unsigned int jetIndex = 0;
	      for ( ; aboveEtCutJetVec_itr != aboveEtCutJetVec.end(); ++aboveEtCutJetVec_itr, 
		      jetIndex++ ) {
		double jetEt  = aboveEtCutJetVec_itr->et();
		double jetEta = aboveEtCutJetVec_itr->eta();
		double jetPhi = aboveEtCutJetVec_itr->phi();
		//		std::cout << "jet[" << jetIndex << "]: et: " << jetEt
		//			  << " eta: " << jetEta
		//			  << " phi: " << jetPhi << std::endl;
		
		double deltaEta = jetEta - partonEta;
		double deltaPhi = DeltaPhi(jetPhi,partonPhi);
		double deltaR2 = pow(deltaEta,2)+pow(deltaPhi,2);

		//		std::cout << "deltaR2: " << deltaR2 << std::endl;

		//		if ( deltaR2 <= jetPartonDeltaR2Cut_ ) {
		  if ( deltaR2 <= closestJetDeltaR2Vec[partonIndex][3] ) {
		    if ( deltaR2 <= closestJetDeltaR2Vec[partonIndex][2] ) {
		      if ( deltaR2 <= closestJetDeltaR2Vec[partonIndex][1] ) {
			if ( deltaR2 <= closestJetDeltaR2Vec[partonIndex][0] ) {
			  closestJetDeltaR2Vec[partonIndex][3] = closestJetDeltaR2Vec[partonIndex][2];
			  closestJetDeltaR2Vec[partonIndex][2] = closestJetDeltaR2Vec[partonIndex][1];
			  closestJetDeltaR2Vec[partonIndex][1] = closestJetDeltaR2Vec[partonIndex][0];
			  closestJetDeltaR2Vec[partonIndex][0] = deltaR2;
			  closestJetIndexVec[partonIndex][3] = closestJetIndexVec[partonIndex][2];
			  closestJetIndexVec[partonIndex][2] = closestJetIndexVec[partonIndex][1];
			  closestJetIndexVec[partonIndex][1] = closestJetIndexVec[partonIndex][0];
			  closestJetIndexVec[partonIndex][0] = jetIndex;
			} else {
			  closestJetDeltaR2Vec[partonIndex][3] = closestJetDeltaR2Vec[partonIndex][2];
			  closestJetDeltaR2Vec[partonIndex][2] = closestJetDeltaR2Vec[partonIndex][1];
			  closestJetDeltaR2Vec[partonIndex][1] = deltaR2;
			  closestJetIndexVec[partonIndex][3] = closestJetIndexVec[partonIndex][2];
			  closestJetIndexVec[partonIndex][2] = closestJetIndexVec[partonIndex][1];
			  closestJetIndexVec[partonIndex][1] = jetIndex;
			}
		      } else {
			closestJetDeltaR2Vec[partonIndex][3] = closestJetDeltaR2Vec[partonIndex][2];
			closestJetDeltaR2Vec[partonIndex][2] = deltaR2;
			closestJetIndexVec[partonIndex][3] = closestJetIndexVec[partonIndex][2];
			closestJetIndexVec[partonIndex][2] = jetIndex;
		      }
		    } else {
		      closestJetDeltaR2Vec[partonIndex][3] = deltaR2;
		      closestJetIndexVec[partonIndex][3] = jetIndex;
		    }
		  }
		  //		}
	      } // end loop over jet above et cut
	      for (int index = 0; index < 4; index++ ) {
		int jetMatchedIndex = closestJetIndexVec[partonIndex][index];
		double etRes = (aboveEtCutJetVec[jetMatchedIndex].et()-partonPt)/partonPt;
		double jetMass = aboveEtCutJetVec[jetMatchedIndex].p4().M();
		std::cout << "closestJetDeltaR2Vec[" << partonIndex << "][" << index << ": " << jetMatchedIndex 
			  << " etRes (%): " << etRes*100.  << std::endl;
		
	      }
	    } // end loop over partons

	    unsigned int nAss = 0;
	    std::vector<bool> partonAss;
	    for ( unsigned int partonIndex = 0; partonIndex < partonsCandidatesNumber; partonIndex++ ) 
	      partonAss.push_back(kFALSE);

	    for ( unsigned int assIndex = 0; assIndex < partonsCandidatesNumber; assIndex++ ) {
	      //	      std::cout << "assIndex: " << assIndex << std::endl;
	      double minDeltaR2 = 99.;
	      int minPartonIndex = -1;
	      // find the best DeltaR2 matching to find the best parton association in the collection
	      for ( unsigned int partonIndex = 0; partonIndex < partonsCandidatesNumber; partonIndex++ ) {
		//		std::cout << "partonIndex: " << partonIndex << " partonAss: " << partonAss[partonIndex] << std::endl;
		if ( !partonAss[partonIndex] && closestJetDeltaR2Vec[partonIndex][0] <= minDeltaR2 ) {
		  minDeltaR2 = closestJetDeltaR2Vec[partonIndex][0];
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
		ZjetVec.push_back(aboveEtCutJetVec[jetIndex]);
		ZjetInd.push_back(jetIndex);
	      }
	      else if ( minPartonIndex == 2 || minPartonIndex == 3 ) {
		TjetVec.push_back(aboveEtCutJetVec[jetIndex]);
		TjetInd.push_back(jetIndex);
	      }

	      // in case of "non-biunivocity" pop-up jet associations belong to worst DeltaR2 parton
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
		    closestJetDeltaR2Vec[partonIndex][jAssIndex-1] = closestJetDeltaR2Vec[partonIndex][jAssIndex];
		    closestJetIndexVec[partonIndex][jAssIndex-1] = closestJetIndexVec[partonIndex][jAssIndex];
		  }
		}
		//		std::cout << "closestJetIndex: " << closestJetIndexVec[partonIndex][0] << std::endl;
	      } // end loop over partons [not matched yet and w/ the same jet]
	      
	    } // end loop over association index

//	    std::cout << "************************************" << std::endl;
//	    std::cout << "ZjetVec: " << ZjetVec.size() << std::endl;
//	    std::cout << "TjetVec: " << TjetVec.size() << std::endl;

	  } else std::cout << "WARNING!: less than 4 jets above et cut" << std::endl;
	} else std::cout << "WARNING!: less than 2 jets above et cut" << std::endl;
      }
    } 
    delete hadronicZ;
    delete leptonicZ;

  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
VBFHZZMatchingAnalyzer::beginJob(const edm::EventSetup&)
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
  ZpartonsDeltaR_       = new TH1D("ZpartonsDeltaR",      "#DeltaR between partons from Z",   nbin_,    0.,   10.);
  ZpartonsMass_         = new TH1D("ZpartonsMass",        "invariant mass of partons from Z", nbin_,    0., 10000.);
  ZpartonsCollinearity_ = new TH1D("ZpartonsCollinearity","collinearity of partons from Z",   nbin_,-7000., 7000.);

  ZpartonsDeltaRVSjetMass_ = new TH2D("ZpartonsDeltaRVSjetMass","jet mass VS #DelraR between the 2 partons from Z",nbin_,0.,0.8,nbin_,0.,150.);
  ZpartonsEtResVSjetMass_  = new TH2D("ZpartonsEtResVSjetMass", "jet mass VS E_{T} resolution",nbin_,0.,1.,nbin_,0.,150.);
  ZpartonsDeltaRVSetRes_   = new TH2D("ZpartonsDeltaRVSetRes",  "E_{T} resolution VS #DelraR between the 2 partons from Z",nbin_,0.,0.8,nbin_,0.,1.);
  ZpartonsDeltaRVSjetMass_profile_ = new TProfile("ZpartonsDeltaRVSjetMass_profile","jet mass VS #DelraR between the 2 partons from Z",nbin_,0.,0.8,0.,150.);
  ZpartonsEtResVSjetMass_profile_  = new TProfile("ZpartonsEtResVSjetMass_profile", "jet mass VS E_{T} resolution",nbin_,0.,1.,0.,150.);
  ZpartonsDeltaRVSetRes_profile_   = new TProfile("ZpartonsDeltaRVSetRes_profile",  "E_{T} resolution VS #DelraR between the 2 partons from Z",nbin_,0.,0.8,0.,1.);


  TagPartonEta_ = new TH1D("TagPartonEta","#eta of tag partons", nbin_, -8.,   8.);
  TagPartonPt_  = new TH1D("TagPartonPt", "p_{T} of tag partons",nbin_,  0., 500.);
  TagPartonE_   = new TH1D("TagPartonE",  "E of tag partons",    nbin_,-10.,7000.);
  TagPartonsDeltaEta_     = new TH1D("TagPartonsDeltaEta",    "#Delta#eta between tag partons",nbin_,    0.,   10.);
  TagPartonsDeltaR_       = new TH1D("TagPartonsDeltaR",      "#DeltaR between tag partons",   nbin_,    0.,   10.);
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
VBFHZZMatchingAnalyzer::endJob() {

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
  ZpartonsDeltaR_->Write();
  ZpartonsMass_->Write();
  ZpartonsCollinearity_->Write();

  ZpartonsDeltaRVSjetMass_->Write();
  ZpartonsEtResVSjetMass_->Write();
  ZpartonsDeltaRVSetRes_->Write();
  ZpartonsDeltaRVSjetMass_profile_->Write();
  ZpartonsEtResVSjetMass_profile_->Write();
  ZpartonsDeltaRVSetRes_profile_->Write();

  TagPartonEta_->Write();
  TagPartonPt_->Write();
  TagPartonE_->Write();
  TagPartonsDeltaEta_->Write();
  TagPartonsDeltaR_->Write();
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
DEFINE_FWK_MODULE(VBFHZZMatchingAnalyzer);

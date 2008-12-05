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

#include "AnalysisExamples/HiggsZZAnalyzer/interface/VBFHZZkinematicsAnalyzer.h"

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

// ROOT version magic to support TMVA interface changes in newer ROOT
#include <RVersion.h>

#include <TMVA/DataSet.h>
#include <TMVA/Types.h>
#include <TMVA/MethodBase.h>
#include <TMVA/Methods.h>

VBFHZZkinematicsAnalyzer::VBFHZZkinematicsAnalyzer(const edm::ParameterSet& iConfig) :
  conf_( iConfig ),
  offlineJetLabel_(         iConfig.getUntrackedParameter<edm::InputTag>( "OfflineJets"  ) ),
  //  offlineMEtLabel_(         iConfig.getUntrackedParameter<edm::InputTag>( "OfflineMEt" ) ),
  globalMuonLabel_(         iConfig.getUntrackedParameter<edm::InputTag>( "GlobalMuons"  ) ),
  simpleElectronLabel_(     iConfig.getUntrackedParameter<edm::InputTag>( "SimpleElectrons" ) ),
  //  simpleTauLabel_(          iConfig.getUntrackedParameter<edm::InputTag>( "SimpleTaus" ) ),
  //  combinedSVBJetTagsLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "combinedSVBJetTags" ) ),
  //  simVtxLabel_(             iConfig.getUntrackedParameter<edm::InputTag>( "SimVtx" ) ),
  MCParticleLabel_(         iConfig.getUntrackedParameter<std::string>(   "MCParticles" ) ),
  leptonPtCut_(	       iConfig.getUntrackedParameter<double>("leptonPtCut"        ) ),
  jetEtCut_(           iConfig.getUntrackedParameter<double>("jetEtCut"           ) ),
  jetPartonDeltaR2Cut_(iConfig.getUntrackedParameter<double>("jetPartonDeltaR2Cut") ),
  jetLeptonDeltaRCut_( iConfig.getUntrackedParameter<double>("jetLeptonDeltaRCut" ) ),
  jetEMfracCut_(       iConfig.getUntrackedParameter<double>("jetEMfracCut"       ) ),
  tmvaSignalSuffix_(   iConfig.getUntrackedParameter<std::string>("tmvaSignalSuffix") ),
  tmvaCombinSuffix_(   iConfig.getUntrackedParameter<std::string>("tmvaCombinSuffix") ),
  writeTMVA_( iConfig.getUntrackedParameter<bool>("writeTMVA") )
{
  // Now do what ever initialization is needed
  // -----------------------------------------

  // Total number of variables passed to the TMVA
  variablesNumber_ = 48;

  // fill the list of event variables for the TMVA
  // tag jets variables
  eventVariablesNamesVector_.push_back("etq1"    );
  eventVariablesNamesVector_.push_back("etq2"    );
  eventVariablesNamesVector_.push_back("etaq1"   );
  eventVariablesNamesVector_.push_back("etaq2"   );
  eventVariablesNamesVector_.push_back("phiq1"   );
  eventVariablesNamesVector_.push_back("phiq2"   );
  eventVariablesNamesVector_.push_back("dphiqq"  );
  eventVariablesNamesVector_.push_back("detaqq"  );
  eventVariablesNamesVector_.push_back("ptminqq" );
  eventVariablesNamesVector_.push_back("ptmaxqq" );
  eventVariablesNamesVector_.push_back("etaminqq");
  eventVariablesNamesVector_.push_back("etamaxqq");
  // Z jets variables
  eventVariablesNamesVector_.push_back("eth1"    );
  eventVariablesNamesVector_.push_back("eth2"    );
  eventVariablesNamesVector_.push_back("etah1"   );
  eventVariablesNamesVector_.push_back("etah2"   );
  eventVariablesNamesVector_.push_back("phih1"   );
  eventVariablesNamesVector_.push_back("phih2"   );
  eventVariablesNamesVector_.push_back("dphihh"  );
  eventVariablesNamesVector_.push_back("detahh"  );
  eventVariablesNamesVector_.push_back("ptminhh" );
  eventVariablesNamesVector_.push_back("ptmaxhh" );
  eventVariablesNamesVector_.push_back("etaminhh");
  eventVariablesNamesVector_.push_back("etamaxhh");
  // tag jets system variables
  eventVariablesNamesVector_.push_back("ptqq" );
  eventVariablesNamesVector_.push_back("mqq"  );
  eventVariablesNamesVector_.push_back("etaqq");
  // Z jets system variables
  eventVariablesNamesVector_.push_back("pthh" );
  eventVariablesNamesVector_.push_back("mhh"  );
  eventVariablesNamesVector_.push_back("etahh");
  // Z leptons system variables
  eventVariablesNamesVector_.push_back("ptll" );
  eventVariablesNamesVector_.push_back("mll"  );
  eventVariablesNamesVector_.push_back("etall");
  // 2-system variables
  eventVariablesNamesVector_.push_back("dphiTjetZjet");
  eventVariablesNamesVector_.push_back("dphiTjetZlep");
  eventVariablesNamesVector_.push_back("dphiminTZ"   );
  eventVariablesNamesVector_.push_back("detaTjetZjet");
  eventVariablesNamesVector_.push_back("detaTjetZlep");
  eventVariablesNamesVector_.push_back("detaminTZ"   );
  eventVariablesNamesVector_.push_back("dphiZjetZlep");
  eventVariablesNamesVector_.push_back("detaZjetZlep");
  eventVariablesNamesVector_.push_back("massTjetZjet");
  eventVariablesNamesVector_.push_back("massTjetZlep");
  eventVariablesNamesVector_.push_back("massZjetZlep");
  // 3-system variables
  eventVariablesNamesVector_.push_back("massTZZ"     );
  eventVariablesNamesVector_.push_back("etaTZZ"      );
  eventVariablesNamesVector_.push_back("ptTZZ"       );
  eventVariablesNamesVector_.push_back("eventNumber" );

  // reset the auto-pointers
  // signal-like event TTree
  tmvaSignalTreeWriterPtr_.reset(new TMVAtreeWriter<Float_t>(eventVariablesNamesVector_,tmvaSignalSuffix_));
  // combinatorial-like event TTree
  tmvaCombinTreeWriterPtr_.reset(new TMVAtreeWriter<Float_t>(eventVariablesNamesVector_,tmvaCombinSuffix_));


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

  if( !writeTMVA_ ) {

    // Array of variables for the Reader
    variables_ = new Float_t[variablesNumber_];

    // Create the Reader object
    // ------------------------
    reader_ = new TMVA::Reader("!Color");    

    // create a set of variables and declare them to the reader
    // - the variable names must corresponds in name and type to 
    // those given in the weight file(s) that you use

    // all minus the last one, which is the event number

    for( int iVar = 0; iVar < variablesNumber_-1; ++iVar ) {
      reader_->AddVariable( eventVariablesNamesVector_[iVar], &(variables_[iVar]) );
    }

    // book the MVA methods
    //
    string dir    = "AnalysisExamples/HiggsZZAnalyzer/data/";
    string prefix = "TMVAnalysis";
    string bdtName = "_BDT.weights.txt";

    edm::FileInPath fileWithFullPath(dir+prefix+bdtName);

    cout << "booking BDT method" << endl;
    reader_->BookMVA( "BDT method", fileWithFullPath.fullPath() );
    cout << "BDT method booked" << endl;
  }

  gROOT->Time();

}


VBFHZZkinematicsAnalyzer::~VBFHZZkinematicsAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

  if( !writeTMVA_ ) {
    delete[] variables_;
    delete reader_;
  }

  // write out the file
  // ------------------
  OutputFile->Write();

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
VBFHZZkinematicsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
      Z1particlesId = int(fabs(Z1Candidates[0]->pdgId()));
      Z2particlesId = int(fabs(Z2Candidates[0]->pdgId()));
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
	int leptonId = int(fabs(leptonicZcandidates[0]->pdgId()));
	if ( leptonId == pythiamu_ || leptonId == pythiae_ ) {

	  std::vector< const Candidate* >::const_iterator leptonicZcandidates_itr = leptonicZcandidates.begin();
	  for ( ; leptonicZcandidates_itr != leptonicZcandidates.end(); ++leptonicZcandidates_itr ) 
	    ZlepVec.push_back(null_XYZTLorentzVector_);

	  leptonicZcandidates_itr = leptonicZcandidates.begin();
	  unsigned int leptonicZcandidateIndex = 0;

	  for ( ; leptonicZcandidates_itr != leptonicZcandidates.end(); ++leptonicZcandidates_itr, leptonicZcandidateIndex++ ) {
	    double tmp_leptonicZcandidateEta = (*leptonicZcandidates_itr)->eta();
	    double tmp_leptonicZcandidatePhi = (*leptonicZcandidates_itr)->phi();
	    double tmp_leptonicZcandidateCharge = (*leptonicZcandidates_itr)->charge();
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
		if (  tmp_leptonicZcandidateCharge*tmp_leptonCharge > 0  ) {
		  if ( tmp_deltaR2 <= leptonDeltaR2 ) {
		    ZlepVec[leptonicZcandidateIndex] = goodMuonVec_itr->p4();
		    leptonDeltaR2 = tmp_deltaR2;
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
		if (  tmp_leptonicZcandidateCharge*tmp_leptonCharge > 0  ) {
		  if ( tmp_deltaR2 <= leptonDeltaR2 ) {
		    ZlepVec[leptonicZcandidateIndex] = goodElectronVec_itr->p4();
		    leptonDeltaR2 = tmp_deltaR2;
		  }
		}
	      }
	    }
	  }
	  //	  std::cout << "ZlepVec: " << ZlepVec.size() << std::endl;
	  //	  if ( ZlepVec[0].pt() == 0. || ZlepVec[1].pt() == 0 )
	  //	    std::cout << "WARNING: ZlepVec[0]->pt: " << ZlepVec[0].pt() << " <--> ZlepVec[1]->pt: " << ZlepVec[1].pt() << std::endl;
	}
      }
    }
    if ( ZlepVec.size() >= nZleptons_ && ZlepVec[0].pt() != 0. && ZlepVec[1].pt() != 0. ) {
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
	  if ( tmpOfflinejetEt >= jetEtCut_ ) {
	    aboveEtCutJetVec.push_back(*offlinejet);
	    jetEta_aboveEtCut_->Fill(tmpOfflinejetEta);
	  }
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
	    for ( unsigned int index = 0; index < hadronicZcandidatesNumber; index++) 
	      ZjetDeltaR2Vec.push_back(99.);
	    
	    std::vector<OfflineJet> TjetVec;
	    std::vector< double >   TjetDeltaR2Vec;
	    std::vector<int>        TjetInd;
	    unsigned int TagCandidatesNumber = TagCandidates.size();
	    for ( unsigned int index = 0; index < TagCandidatesNumber; index++) 
	      TjetDeltaR2Vec.push_back(99.);
	    

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
		//		std::cout << "closestJetDeltaR2Vec[" << partonIndex << "][" << index << "]: " << jetMatchedIndex 
		//			  << " --> deltaR2: " << closestJetDeltaR2Vec[partonIndex][index] << std::endl;
		
	      }
	    } // end loop over partons
	  


	    std::vector<bool> partonAss;
	    for ( unsigned int partonIndex = 0; partonIndex < partonsCandidatesNumber; partonIndex++ ) 
	      partonAss.push_back(kFALSE);

	    for ( unsigned int assIndex = 0; assIndex < partonsCandidatesNumber; assIndex++ ) {
	      //	      std::cout << "assIndex: " << assIndex << std::endl;
	      double minDeltaR2 = jetPartonDeltaR2Cut_;
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
	      if ( minPartonIndex >= 0 ) {
		partonAss[minPartonIndex] = kTRUE;
		// save the matched jet into the proper vector
		unsigned int jetIndex = closestJetIndexVec[minPartonIndex][0];
		//	      std::cout << "minPartonIndex: " << minPartonIndex 
		//			<< " to jetIndex: " << jetIndex
		//			<< " => nAss: " << nAss << std::endl;

		if ( minPartonIndex == 0 || minPartonIndex == 1 ) {
		  ZjetVec.push_back(aboveEtCutJetVec[jetIndex]);
		  ZjetInd.push_back(jetIndex);
		  //		  std::cout << "Z jetIndex: " <<  jetIndex << std::endl;
		}
		else if ( minPartonIndex == 2 || minPartonIndex == 3 ) {
		  TjetVec.push_back(aboveEtCutJetVec[jetIndex]);
		  TjetInd.push_back(jetIndex);
		  //		  std::cout << "T jetIndex: " <<  jetIndex << std::endl;
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
	      }
	    } // end loop over association index

	    //	    std::cout << "************************************" << std::endl;
	    //	    std::cout << "ZjetVec: " << ZjetVec.size() << std::endl;
	    //	    std::cout << "TjetVec: " << TjetVec.size() << std::endl;
	    //	    std::cout << "ZlepVec: " << ZlepVec.size() << std::endl;


	    // ====================================================================================================
	    if (ZjetVec.size()==2 && TjetVec.size() == 2 ){
	      //	    if (partonAss[0] && partonAss[1] && partonAss[2] && partonAss[3] ) {

	      Particle * Zjet = new Particle(0,ZjetVec[0].p4()+ZjetVec[1].p4(), ZjetVec[0].vertex(),     pythiaZ_,0,true);
	      Particle * Zlep = new Particle(0,ZlepVec[0]+ZlepVec[1],           math::XYZPoint(0.,0.,0.),pythiaZ_,0,true);
	      Particle * Tjet = new Particle(0,TjetVec[0].p4()+TjetVec[1].p4(), TjetVec[0].vertex(),     0,       0,true);

	      if ( ZjetVec[0].eta() == 0. ||
		   ZjetVec[1].eta() == 0. ||
		   TjetVec[0].eta() == 0. ||
		   TjetVec[1].eta() == 0. ||
		   ZjetVec[0].et() < jetEtCut_ ||
		   ZjetVec[1].et() < jetEtCut_ ||
		   TjetVec[0].et() < jetEtCut_ ||
		   TjetVec[1].et() < jetEtCut_ ) {
//		std::cout << "WARNING!!!" << std::endl;
//		std::cout << "ZjetVec[0].et: " << ZjetVec[0].et() << " eta: " << ZjetVec[0].eta() << " phi: " << ZjetVec[0].phi() << std::endl;
//		std::cout << "ZjetVec[1].et: " << ZjetVec[1].et() << " eta: " << ZjetVec[1].eta() << " phi: " << ZjetVec[1].phi() << std::endl;
//		std::cout << "ZlepVec[0].pt: " << ZlepVec[0].pt() << " eta: " << ZlepVec[0].eta() << " phi: " << ZlepVec[0].phi() << std::endl;
//		std::cout << "ZlepVec[1].pt: " << ZlepVec[1].pt() << " eta: " << ZlepVec[1].eta() << " phi: " << ZlepVec[1].phi() << std::endl;
//		std::cout << "TjetVec[0].et: " << TjetVec[0].et() << " eta: " << TjetVec[0].eta() << " phi: " << TjetVec[0].phi() << std::endl;
//		std::cout << "TjetVec[1].et: " << TjetVec[1].et() << " eta: " << TjetVec[1].eta() << " phi: " << TjetVec[1].phi() << std::endl;
		//std::cout << "Zjet.et: " << Zjet->pt() << std::endl;
		//std::cout << "Zlep.et: " << Zlep->pt() << std::endl;
		//std::cout << "Tjet.et: " << Tjet->pt() << std::endl;
		//std::cout << "Indices: " << TjetInd[0] << " " << TjetInd[1] << " " << ZjetInd[0] << " " << ZjetInd[1] << std::endl;
	      }

	      if ( ZlepVec[0].pt() != 0. && ZlepVec[1].pt() !=0. ) {

		// Study kinematics of true ZZqq system
		// ------------------------------------
		// qq: tag jets
		// ------------
		double dphiqq = DeltaPhi(TjetVec[0].phi(),TjetVec[1].phi());
		double detaqq = fabs(TjetVec[0].eta()-TjetVec[1].eta());
		double ptminqq = TjetVec[0].et();
		double ptmaxqq = TjetVec[1].et();
		if ( TjetVec[1].et()<ptminqq ) {
		  ptminqq = TjetVec[1].et();
		  ptmaxqq = TjetVec[0].et();
		}
		double etaminqq = fabs(TjetVec[0].eta());
		double etamaxqq = fabs(TjetVec[1].eta());
		if ( fabs(TjetVec[1].eta())<etaminqq ) {
		  etaminqq = fabs(TjetVec[1].eta());
		  etamaxqq = fabs(TjetVec[0].eta());
		}
		Dphiqq->Fill(dphiqq);
		Detaqq->Fill(detaqq);
		Ptminqq->Fill(ptminqq);
		Ptmaxqq->Fill(ptmaxqq);
		Etaminqq->Fill(etaminqq);
		Etamaxqq->Fill(etamaxqq);

		eventSignalVariablesVector_.push_back(TjetVec[0].et());
		eventSignalVariablesVector_.push_back(TjetVec[1].et());
		eventSignalVariablesVector_.push_back(TjetVec[0].eta());
		eventSignalVariablesVector_.push_back(TjetVec[1].eta());
		eventSignalVariablesVector_.push_back(TjetVec[0].phi());
		eventSignalVariablesVector_.push_back(TjetVec[1].phi());
		eventSignalVariablesVector_.push_back(dphiqq);
		eventSignalVariablesVector_.push_back(detaqq);
		eventSignalVariablesVector_.push_back(ptminqq);
		eventSignalVariablesVector_.push_back(ptmaxqq);
		eventSignalVariablesVector_.push_back(etaminqq);
		eventSignalVariablesVector_.push_back(etamaxqq);

		// hh: Z jets
		// ----------
		double dphihh = DeltaPhi(ZjetVec[0].phi(),ZjetVec[1].phi());
		double detahh = fabs(ZjetVec[0].eta()-ZjetVec[1].eta());
		double ptminhh = ZjetVec[0].et();
		double ptmaxhh = ZjetVec[1].et();
		if ( ZjetVec[1].et()<ptminhh ) {
		  ptminhh = ZjetVec[1].et();
		  ptmaxhh = ZjetVec[0].et();
		}
		double etaminhh = fabs(ZjetVec[0].eta());
		double etamaxhh = fabs(ZjetVec[1].eta());
		if ( fabs(ZjetVec[1].eta())<etaminhh ) {
		  etaminhh = fabs(ZjetVec[1].eta());
		  etamaxhh = fabs(ZjetVec[0].eta());
		}
		Dphihh->Fill(dphihh);
		Detahh->Fill(detahh);
		Ptminhh->Fill(ptminhh);
		Ptmaxhh->Fill(ptmaxhh);
		Etaminhh->Fill(etaminhh);
		Etamaxhh->Fill(etamaxhh);
		hadronicZrecDeltaPhiVSDeltaEta_->Fill(detahh,dphihh);

		eventSignalVariablesVector_.push_back(ZjetVec[0].et());
		eventSignalVariablesVector_.push_back(ZjetVec[1].et());
		eventSignalVariablesVector_.push_back(ZjetVec[0].eta());
		eventSignalVariablesVector_.push_back(ZjetVec[1].eta());
		eventSignalVariablesVector_.push_back(ZjetVec[0].phi());
		eventSignalVariablesVector_.push_back(ZjetVec[1].phi());
		eventSignalVariablesVector_.push_back(dphihh);
		eventSignalVariablesVector_.push_back(detahh);
		eventSignalVariablesVector_.push_back(ptminhh);
		eventSignalVariablesVector_.push_back(ptmaxhh);
		eventSignalVariablesVector_.push_back(etaminhh);
		eventSignalVariablesVector_.push_back(etamaxhh);

		// qq system
		// ---------
		double ptqq = Tjet->pt();
		double mqq = Tjet->mass();
		double etaqq = fabs(Tjet->eta());
		// hh system
		// ---------
		double pthh = Zjet->pt();
		double mhh = Zjet->mass();
		double etahh = fabs(Zjet->eta());
		// ll system
		// ---------
		double ptll = Zlep->pt();
		double mll = Zlep->mass();
		double etall = fabs(Zlep->eta());
		//std::cout << "Systems Pt:   " << ptqq  << " " << pthh  << " " << ptll  << std::endl;
		//std::cout << "Systems eta:  " << etaqq << " " << etahh << " " << etall << std::endl;
		//std::cout << "Systems mass: " << mqq   << " " << mhh   << " " << mll   << std::endl;
		Ptqq->Fill(ptqq);
		Mqq->Fill(mqq);
		Etaqq->Fill(etaqq);
		Pthh->Fill(pthh);
		Mhh->Fill(mhh);
		Etahh->Fill(etahh);
		Ptll->Fill(ptll);
		Mll->Fill(mll);
		Etall->Fill(etall);

		eventSignalVariablesVector_.push_back(ptqq);
		eventSignalVariablesVector_.push_back(mqq);
		eventSignalVariablesVector_.push_back(etaqq);
		eventSignalVariablesVector_.push_back(pthh);
		eventSignalVariablesVector_.push_back(mhh);
		eventSignalVariablesVector_.push_back(etahh);
		eventSignalVariablesVector_.push_back(ptll);
		eventSignalVariablesVector_.push_back(mll);
		eventSignalVariablesVector_.push_back(etall);

		// 2-particle vars
		// ---------------
		double dphiTjetZjet = DeltaPhi(Tjet,Zjet);
		double dphiTjetZlep = DeltaPhi(Tjet,Zlep);
		double dphiminTZ = dphiTjetZjet;
		if ( dphiminTZ>dphiTjetZlep ) dphiminTZ = dphiTjetZlep;
		double detaTjetZjet = fabs(Tjet->eta()-Zjet->eta());
		double detaTjetZlep = fabs(Tjet->eta()-Zlep->eta());
		double detaminTZ = detaTjetZjet;
		if ( detaminTZ>detaTjetZlep ) detaminTZ = detaTjetZlep;
		double dphiZjetZlep = DeltaPhi(Zjet,Zlep);
		double detaZjetZlep = fabs(Zjet->eta()-Zlep->eta());
		double massTjetZjet = (Tjet->p4()+Zjet->p4()).mass();
		double massTjetZlep = (Tjet->p4()+Zlep->p4()).mass();
		double massZjetZlep = (Zjet->p4()+Zlep->p4()).mass();
		DphiTjetZjet->Fill(dphiTjetZjet);
		DphiTjetZlep->Fill(dphiTjetZlep);
		DphiminTZ->Fill(dphiminTZ);
		DetaTjetZjet->Fill(detaTjetZjet);
		DetaTjetZlep->Fill(detaTjetZlep);
		DetaminTZ->Fill(detaminTZ);
		DphiZjetZlep->Fill(dphiZjetZlep);
		DetaZjetZlep->Fill(detaZjetZlep);
		MassTjetZjet->Fill(massTjetZjet);
		MassTjetZlep->Fill(massTjetZlep);
		MassZjetZlep->Fill(massZjetZlep);

		eventSignalVariablesVector_.push_back(dphiTjetZjet);
		eventSignalVariablesVector_.push_back(dphiTjetZlep);
		eventSignalVariablesVector_.push_back(dphiminTZ);
		eventSignalVariablesVector_.push_back(detaTjetZjet);
		eventSignalVariablesVector_.push_back(detaTjetZlep);
		eventSignalVariablesVector_.push_back(detaminTZ);
		eventSignalVariablesVector_.push_back(dphiZjetZlep);
		eventSignalVariablesVector_.push_back(detaZjetZlep);
		eventSignalVariablesVector_.push_back(massTjetZjet);
		eventSignalVariablesVector_.push_back(massTjetZlep);
		eventSignalVariablesVector_.push_back(massZjetZlep);

		// 3-particle vars
		// ---------------
		double massTZZ = (Tjet->p4()+Zjet->p4()+Zlep->p4()).mass();
		double etaTZZ = fabs((Tjet->p4()+Zjet->p4()+Zlep->p4()).eta());
		double ptTZZ = (Tjet->p4()+Zjet->p4()+Zlep->p4()).pt();
		MassTZZ->Fill(massTZZ);
		EtaTZZ->Fill(etaTZZ);
		PtTZZ->Fill(ptTZZ);
	
		eventSignalVariablesVector_.push_back(massTZZ);
		eventSignalVariablesVector_.push_back(etaTZZ);
		eventSignalVariablesVector_.push_back(ptTZZ);

                // Store the event number
                eventSignalVariablesVector_.push_back(eventcounter_);

                // Variable to hold the value of the discriminant for the true combination
                double discriminantSignal = 0.;

                if( writeTMVA_ ) {
                  // Fill the tree for the TMVA
                  tmvaSignalTreeWriterPtr_->fill(eventSignalVariablesVector_);
                }
                else {
                  // Fill the values in the array used by the reader
                  vector<Float_t>::const_iterator varIt = eventSignalVariablesVector_.begin();
                  int iVarSignal = 0;
                  for( ; varIt != eventSignalVariablesVector_.end(); ++varIt, ++iVarSignal ) {
                    variables_[iVarSignal] = *varIt;
                  }
                  discriminantSignal = reader_->EvaluateMVA( "BDT method" );
                }
                eventSignalVariablesVector_.clear();

                // Vector to hold the weights of all the combinatorics
                vector<double> discriminantBackground;

		//std::cout << "Total system mass: " << massTZZ << " Eta: " << etaTZZ << " Pt: " << ptTZZ << std::endl;

		// Loop on true combination and all false ones, and derive kinematics of configuration
		// -----------------------------------------------------------------------------------
		//	    delete Tjet;
		//	    delete Zjet;
		int i1=0;
		int i2=0;
		int i3=0;
		int i4=0;
		std::vector<OfflineJet>::const_iterator jet1 = aboveEtCutJetVec.begin();
		for ( ; jet1 != aboveEtCutJetVec.end()-1; ++jet1, i1++ ) {
		  std::vector<OfflineJet>::const_iterator jet2 = jet1+1;
		  for ( ; jet2 != aboveEtCutJetVec.end(); ++jet2, i2++) {
		    if ( jet2==jet1 ) continue;
		    std::vector<OfflineJet>::const_iterator jet3 = aboveEtCutJetVec.begin();
		    for ( ; jet3 != aboveEtCutJetVec.end()-1; ++jet3, i3++) {
		      if ( jet3==jet1 || jet3==jet2 ) continue;
		      std::vector<OfflineJet>::const_iterator jet4 = jet3+1;
		      for ( ; jet4 != aboveEtCutJetVec.end(); ++jet4, i4++) {
			if ( jet4==jet1 || jet4==jet2 || jet4==jet3 ) continue;
			// Remove true combo and homologues
			// --------------------------------
			if (i1==TjetInd[0] && i2==TjetInd[1] && i3==ZjetInd[0] && i4==ZjetInd[1]) continue;
			if (i1==TjetInd[1] && i2==TjetInd[0] && i3==ZjetInd[0] && i4==ZjetInd[1]) continue;
			if (i1==TjetInd[0] && i2==TjetInd[1] && i3==ZjetInd[1] && i4==ZjetInd[0]) continue;
			if (i1==TjetInd[1] && i2==TjetInd[0] && i3==ZjetInd[1] && i4==ZjetInd[0]) continue;
			// Ok, here we have a sensible quadruplet
			// --------------------------------------
			//std::cout << "Combination " << i1 << i2 << i3 << i4 << std::endl;
			Particle *Tjet = new Particle (0,jet1->p4()+jet2->p4(), jet1->vertex(),0,0,true);
			Particle *Zjet = new Particle (0,jet3->p4()+jet4->p4(), jet3->vertex(),0,0,true);

			// Study kinematics of random ZZqq systems
			// ---------------------------------------
			// qq: tag jets
			// ------------
			double f_dphiqq = DeltaPhi(jet1->phi(),jet2->phi());
			double f_detaqq = fabs(jet1->eta()-jet2->eta());
			double f_ptminqq = jet1->et();
			double f_ptmaxqq = jet2->et();
			if ( jet2->et()<f_ptminqq ) {
			  f_ptminqq = jet2->et();
			  f_ptmaxqq = jet1->et();
			}
			double f_etaminqq = fabs(jet1->eta());
			double f_etamaxqq = fabs(jet2->eta());
			if ( fabs(jet2->eta())<f_etaminqq ) {
			  f_etaminqq = fabs(jet2->eta());
			  f_etamaxqq = fabs(jet1->eta());
			}
			F_Dphiqq->Fill(f_dphiqq);
			F_Detaqq->Fill(f_detaqq);
			F_Ptminqq->Fill(f_ptminqq);
			F_Ptmaxqq->Fill(f_ptmaxqq);
			F_Etaminqq->Fill(f_etaminqq);
			F_Etamaxqq->Fill(f_etamaxqq);
			eventCombinVariablesVector_.push_back(jet1->et());
			eventCombinVariablesVector_.push_back(jet2->et());
			eventCombinVariablesVector_.push_back(jet1->eta());
			eventCombinVariablesVector_.push_back(jet2->eta());
			eventCombinVariablesVector_.push_back(jet1->phi());
			eventCombinVariablesVector_.push_back(jet2->phi());
			eventCombinVariablesVector_.push_back(f_dphiqq);
			eventCombinVariablesVector_.push_back(f_detaqq);
			eventCombinVariablesVector_.push_back(f_ptminqq);
			eventCombinVariablesVector_.push_back(f_ptmaxqq);
			eventCombinVariablesVector_.push_back(f_etaminqq);
			eventCombinVariablesVector_.push_back(f_etamaxqq);

			// hh: Z jets
			// ----------
			double f_dphihh = DeltaPhi(jet3->phi(),jet4->phi());
			double f_detahh = fabs(jet3->eta()-jet4->eta());
			double f_ptminhh = jet3->et();
			double f_ptmaxhh = jet4->et();
			if ( jet4->et()<f_ptminhh ) {
			  f_ptminhh = jet4->et();
			  f_ptmaxhh = jet3->et();
			}
			double f_etaminhh = fabs(jet3->eta());
			double f_etamaxhh = fabs(jet4->eta());
			if ( fabs(jet4->eta())<f_etaminhh ) {
			  f_etaminhh = fabs(jet4->eta());
			  f_etamaxhh = fabs(jet3->eta());
			}
			F_Dphihh->Fill(f_dphihh);
			F_Detahh->Fill(f_detahh);
			F_Ptminhh->Fill(f_ptminhh);
			F_Ptmaxhh->Fill(f_ptmaxhh);
			F_Etaminhh->Fill(f_etaminhh);
			F_Etamaxhh->Fill(f_etamaxhh);

			eventCombinVariablesVector_.push_back(jet3->et());
			eventCombinVariablesVector_.push_back(jet4->et());
			eventCombinVariablesVector_.push_back(jet3->eta());
			eventCombinVariablesVector_.push_back(jet4->eta());
			eventCombinVariablesVector_.push_back(jet3->phi());
			eventCombinVariablesVector_.push_back(jet4->phi());
			eventCombinVariablesVector_.push_back(f_dphihh);
			eventCombinVariablesVector_.push_back(f_detahh);
			eventCombinVariablesVector_.push_back(f_ptminhh);
			eventCombinVariablesVector_.push_back(f_ptmaxhh);
			eventCombinVariablesVector_.push_back(f_etaminhh);
			eventCombinVariablesVector_.push_back(f_etamaxhh);

			// qq system
			// ---------
			double f_ptqq = (Tjet->p4()).pt();
			double f_mqq = (Tjet->p4()).mass();
			double f_etaqq = fabs((Tjet->p4()).eta());
			// hh system
			// ---------
			double f_pthh = (Zjet->p4()).pt();
			double f_mhh = (Zjet->p4()).mass();
			double f_etahh = fabs((Zjet->p4()).eta());
			// ll system
			// ---------
			double f_ptll = (Zlep->p4()).pt();
			double f_mll = Zlep->mass();
			double f_etall = fabs(Zlep->eta());
			F_Ptqq->Fill(f_ptqq);
			F_Mqq->Fill(f_mqq);
			F_Etaqq->Fill(f_etaqq);
			F_Pthh->Fill(f_pthh);
			F_Mhh->Fill(f_mhh);
			F_Etahh->Fill(f_etahh);
			F_Ptll->Fill(f_ptll);
			F_Mll->Fill(f_mll);
			F_Etall->Fill(f_etall);

			eventCombinVariablesVector_.push_back(f_ptqq);
			eventCombinVariablesVector_.push_back(f_mqq);
			eventCombinVariablesVector_.push_back(f_etaqq);
			eventCombinVariablesVector_.push_back(f_pthh);
			eventCombinVariablesVector_.push_back(f_mhh);
			eventCombinVariablesVector_.push_back(f_etahh);
			eventCombinVariablesVector_.push_back(f_ptll);
			eventCombinVariablesVector_.push_back(f_mll);
			eventCombinVariablesVector_.push_back(f_etall);
			// 2-particle vars
			// ---------------
			double f_dphiTjetZjet = DeltaPhi((Tjet->p4()).phi(),(Zjet->p4()).phi());
			double f_dphiTjetZlep = DeltaPhi((Tjet->p4()).phi(),Zlep->phi());
			double f_dphiminTZ = f_dphiTjetZjet;
			if ( f_dphiminTZ>f_dphiTjetZlep ) f_dphiminTZ = f_dphiTjetZlep;
			double f_detaTjetZjet = fabs((Tjet->p4()).eta()-(Zjet->p4()).eta());
			double f_detaTjetZlep = fabs((Tjet->p4()).eta()-Zlep->eta());
			double f_detaminTZ = f_detaTjetZjet;
			if ( f_detaminTZ>f_detaTjetZlep ) f_detaminTZ = f_detaTjetZlep;
			double f_dphiZjetZlep = DeltaPhi((Zjet->p4()).phi(),Zlep->phi());
			double f_detaZjetZlep = fabs((Zjet->p4()).eta()-Zlep->eta());
			double f_massTjetZjet = (Tjet->p4()+Zjet->p4()).mass();
			double f_massTjetZlep = (Tjet->p4()+Zlep->p4()).mass();
			double f_massZjetZlep = (Zjet->p4()+Zlep->p4()).mass();
			F_DphiTjetZjet->Fill(f_dphiTjetZjet);
			F_DphiTjetZlep->Fill(f_dphiTjetZlep);
			F_DphiminTZ->Fill(f_dphiminTZ);
			F_DetaTjetZjet->Fill(f_detaTjetZjet);
			F_DetaTjetZlep->Fill(f_detaTjetZlep);
			F_DetaminTZ->Fill(f_detaminTZ);
			F_DphiZjetZlep->Fill(f_dphiZjetZlep);
			F_DetaZjetZlep->Fill(f_detaZjetZlep);
			F_MassTjetZjet->Fill(f_massTjetZjet);
			F_MassTjetZlep->Fill(f_massTjetZlep);
			F_MassZjetZlep->Fill(f_massZjetZlep);

			eventCombinVariablesVector_.push_back(f_dphiTjetZjet);
			eventCombinVariablesVector_.push_back(f_dphiTjetZlep);
			eventCombinVariablesVector_.push_back(f_dphiminTZ);
			eventCombinVariablesVector_.push_back(f_detaTjetZjet);
			eventCombinVariablesVector_.push_back(f_detaTjetZlep);
			eventCombinVariablesVector_.push_back(f_detaminTZ);
			eventCombinVariablesVector_.push_back(f_dphiZjetZlep);
			eventCombinVariablesVector_.push_back(f_detaZjetZlep);
			eventCombinVariablesVector_.push_back(f_massTjetZjet);
			eventCombinVariablesVector_.push_back(f_massTjetZlep);
			eventCombinVariablesVector_.push_back(f_massZjetZlep);
			// 3-particle vars
			// ---------------
			double f_massTZZ = (Tjet->p4()+Zjet->p4()+Zlep->p4()).mass();
			double f_etaTZZ = fabs((Tjet->p4()+Zjet->p4()+Zlep->p4()).eta());
			double f_ptTZZ = (Tjet->p4()+Zjet->p4()+Zlep->p4()).pt();
			F_MassTZZ->Fill(f_massTZZ);
			F_EtaTZZ->Fill(f_etaTZZ);
			F_PtTZZ->Fill(f_ptTZZ);

			eventCombinVariablesVector_.push_back(f_massTZZ);
			eventCombinVariablesVector_.push_back(f_etaTZZ);
			eventCombinVariablesVector_.push_back(f_ptTZZ);
		    
                        // Store the event number
                        eventCombinVariablesVector_.push_back(eventcounter_);

                        if( writeTMVA_ ) {
                          // Fill the tree for the TMVA
                          tmvaCombinTreeWriterPtr_->fill(eventCombinVariablesVector_);
                        }
                        else {
                          // Fill the values in the array used by the reader
                          vector<Float_t>::const_iterator varIt = eventSignalVariablesVector_.begin();
                          int iVarBackground = 0;
                          for( ; varIt != eventCombinVariablesVector_.end(); ++varIt, ++iVarBackground ) {
                            variables_[iVarBackground] = *varIt;
                          }
                          discriminantBackground.push_back( reader_->EvaluateMVA( "BDT method" ) );
                        }
                        eventCombinVariablesVector_.clear();
		      }
		    }
		  }
		}

                if( !writeTMVA_ ) {
                  // After the loop on all combinations sort the vector of background discriminants
                  sort( discriminantBackground.rbegin(), discriminantBackground.rend() );
                  int i = 0;
                  for( vector<double>::const_iterator discr = discriminantBackground.begin(); discr != discriminantBackground.end(); ++discr, ++i ) {
                    cout << "discr["<<i<<"] = " << *discr << endl;
                  }
                }

		// =============================================================================================
	      }

	      // =============================================================================================

	    }
	  } else std::cout << "WARNING!: less than 4 jets above et cut" << std::endl;
	} else std::cout << "WARNING!: less than 2 jets above et cut" << std::endl;
      } // end if partonsCandidates.size() >= njets_
    } // end if ZlepVec.size() >= nZleptons_
    delete hadronicZ;
    delete leptonicZ;
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
VBFHZZkinematicsAnalyzer::beginJob(const edm::EventSetup&)
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


  ZpartonEta_ = new TH1D("ZpartonEta","#eta of partons from Z", nbin_, -8.,   8.);
  ZpartonPt_  = new TH1D("ZpartonPt", "p_{T} of partons from Z",nbin_,  0., 500.);
  ZpartonE_   = new TH1D("ZpartonE",  "E of partons from Z",    nbin_,-10.,7000.);
  ZpartonsDeltaEta_     = new TH1D("ZpartonsDeltaEta",    "#Delta#eta between partons from Z",nbin_,    0.,   10.);
  ZpartonsDeltaR_       = new TH1D("ZpartonsDeltaR",      "#DeltaR between partons from Z",   nbin_,    0.,   10.);
  ZpartonsMass_         = new TH1D("ZpartonsMass",        "invariant mass of partons from Z", nbin_,    0., 10000.);
  ZpartonsCollinearity_ = new TH1D("ZpartonsCollinearity","collinearity of partons from Z",   nbin_,-7000., 7000.);

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

  overlapJetsDeltaEt_  = new TH1D("overlapJetsDeltaEt", "#DeltaE_{T} between 2 overlapping jets",nbin_,-100.,100.);
  overlapJetsDeltaEta_ = new TH1D("overlapJetsDeltaEta","#Delta#eta between 2 overlapping jets", nbin_,  -8.,  8.);
  overlapJetsDeltaEtVSDeltaR2_  = new TH2D("overlapJetsDeltaEtVSDeltaR2", "#DeltaE_{T} VS #DeltaR^{2} between 2 overlapping jets",nbin_, 0.,10.,nbin_,-100.,100.);
  overlapJetsDeltaEtaVSDeltaR2_ = new TH2D("overlapJetsDeltaEtaVSDeltaR2","#Delta#eta VS #DeltaR^{2} between 2 overlapping jets", nbin_, 0.,10.,nbin_,  -8.,  8.);
  overlapJetsDeltaEtVSDeltaEta_ = new TH2D("overlapJetsDeltaEtVSDeltaEta","#DeltaE_{T} VS #Delta#eta between 2 overlapping jets", nbin_,-8., 8.,nbin_,-100.,100.);

  muonJetMatchedNumber_ = new TH1D("muonJetMatchedNumber","number of jets w/ at least 1 muon in a distance #DeltaR<0.5", 10,0.,10.0);
  muonsMatchedNumber_   = new TH1D("muonsMatchedNumber",  "fraction of muons in a distance #DeltaR<0.5 to a jet",       110,0., 1.1);

  jetParton_deltaR2_      = new TH1D("jetParton_deltaR2",     "#DeltaR^{2} between jet and partons from Z",    nbin_,0.,4.0);
  jetParton_deltaR2_zoom_ = new TH1D("jetParton_deltaR2_zoom","#DeltaR^{2} between jet and partons from Z",    nbin_,0.,0.1);
  jetParton_deltaR2max_   = new TH1D("jetParton_deltaR2max",  "max #DeltaR^{2} between jet and partons from Z",nbin_,0.,4.0);

  jetParton_deltaEtVSdeltaR2_         = new TH2D("jetParton_deltaEtVSdeltaR2",        "#DeltaE_{T} VS #DeltaR^{2} between partons and matched jets",          nbin_,0.,0.25,nbin_,-50.,50.);
  jetParton_deltaEtmeanVSdeltaR2mean_ = new TH2D("jetParton_deltaEtmeanVSdeltaR2mean","mean #DeltaE_{T} VS mean #DeltaR^{2} between partons and matched jets",nbin_,0.,0.25,nbin_,  0.,50.);
  jetParton_deltaEVSdeltaR2_   = new TH2D("jetParton_deltaEVSdeltaR2",  "#DeltaE VS #DeltaR^{2} between partons and matched jets",    nbin_,0.,0.25,nbin_,-50.,50.);
  jetParton_deltaEtaVSdeltaR2_ = new TH2D("jetParton_deltaEtaVSdeltaR2","#Delta#eta VS #DeltaR^{2} between partons and matched jets", nbin_,0.,0.25,nbin_, -5., 5.);
  jetParton_deltaEtVSdeltaR2_profile_         = new TProfile("jetParton_deltaEtVSdeltaR2_profile",        "#DeltaE_{T} VS #DeltaR^{2} between partons and matched jets",          nbin_,0.,0.25,-50.,50.);
  jetParton_deltaEtmeanVSdeltaR2mean_profile_ = new TProfile("jetParton_deltaEtmeanVSdeltaR2mean_profile","mean #DeltaE_{T} VS mean #DeltaR^{2} between partons and matched jets",nbin_,0.,0.25,  0.,50.);
  jetParton_deltaEVSdeltaR2_profile_   = new TProfile("jetParton_deltaEVSdeltaR2_profile",  "#DeltaE VS #DeltaR^{2} between partons and matched jets",    nbin_,0.,0.25,-50.,50.);
  jetParton_deltaEtaVSdeltaR2_profile_ = new TProfile("jetParton_deltaEtaVSdeltaR2_profile","#Delta#eta VS #DeltaR^{2} between partons and matched jets", nbin_,0.,0.25, -5., 5.);
  jetParton_deltaEta_ = new TH1D("jetParton_deltaEta","#Delta#eta between partons and matched jets", nbin_, -5., 5.);
  jetParton_deltaEt_  = new TH1D("jetParton_deltaEt", "#DeltaE_{T} between partons and matched jets",nbin_,-50.,50.);

  hadronicZrecMass_           = new TH1D("hadronicZrecMass",          "reconstructed Z mass from jets matched to Z partons",           nbin_, 0.,150.);
  hadronicZrecMassResolution_ = new TH1D("hadronicZrecMassResolution","reconstructed Z mass resolution from jets matched to Z partons",nbin_,-1.,  1.);
  hadronicZrecPt_           = new TH1D("hadronicZrecPt",          "reconstructed Z p_{T} from jets matched to Z partons",           nbin_, 0.,500.);
  hadronicZrecPtResolution_ = new TH1D("hadronicZrecPtResolution","reconstructed Z p_{T} resolution from jets matched to Z partons",nbin_,-1.,  1.);
  hadronicZrecCollinearity_           = new TH1D("hadronicZrecCollinearity",          "reconstructed collinearity between jets matched to Z partons",        nbin_, -6000.,6000.);
  hadronicZrecCollinearityResolution_ = new TH1D("hadronicZrecCollinearityResolution","reconstructed collinearity resolution from jets matched to Z partons",nbin_,    -1.,   1.);
  hadronicZrecDeltaPhiVSDeltaEta_ = new TH2D("hadronicZrecDeltaPhiVSDeltaEta","",nbin_,0.,5.,nbin_,0.,TMath::Pi());

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
  Bothok_05 = new TH1D ("Bothok_05", "Evt frac. with both tag jets found within 0.5 rads", 41, -0.05, 4.05 );
  Bothok_02 = new TH1D ("Bothok_02", "Evt frac. with both tag jets found within 0.2 rads", 41, -0.05, 4.05 );
  Bothok_05_30pc = new TH1D ("Bothok_05_30pc", "Evt frac. with both tag jets found within 0.5 rads and <30% Et res.", 41, -0.05, 4.05 );
  Bothok_02_30pc = new TH1D ("Bothok_02_30pc", "Evt frac. with both tag jets found within 0.2 rads and <30% Et res.", 41, -0.05, 4.05 );
  Bothok_05_20pc = new TH1D ("Bothok_05_20pc", "Evt frac. with both tag jets found within 0.5 rads and <20% Et res.", 41, -0.05, 4.05 );
  Bothok_02_20pc = new TH1D ("Bothok_02_20pc", "Evt frac. with both tag jets found within 0.2 rads and <20% Et res.", 41, -0.05, 4.05 );
  Bothok_05_10pc = new TH1D ("Bothok_05_10pc", "Evt frac. with both tag jets found within 0.5 rads and <10% Et res.", 41, -0.05, 4.05 );
  Bothok_02_10pc = new TH1D ("Bothok_02_10pc", "Evt frac. with both tag jets found within 0.2 rads and <10% Et res.", 41, -0.05, 4.05 );
  Allevents = new TH1D ("Allevents", "All Evts.", 41, -0.05, 4.05 );
  DeltaEtvsEt = new TProfile("DeltaEtvsEt", "Et offset versus Parton Pt", 20, 0., 200.,-50.,50. );

  // =========================================================

  Dphiqq = new TH1D ("Dphiqq","Delta Phi tag jets", nbin_, 0., TMath::Pi());
  Detaqq = new TH1D ("Detaqq","Delta Eta tag jets", nbin_, 0., 10.0);
  Ptminqq = new TH1D ("Ptminqq","Pt of lowest Pt tag jet", nbin_, 0., 200.);
  Ptmaxqq = new TH1D ("Ptmaxqq","Pt of highest Pt tag jet", nbin_, 0., 500.);
  Etaminqq = new TH1D ("Etaminqq","Min |eta| of tag jet", nbin_, 0., 5.);
  Etamaxqq = new TH1D ("Etamaxqq","Max |eta| of tag jet", nbin_, 0., 5.);
  Dphihh = new TH1D ("Dphihh","Delta Phi Z jets", nbin_, 0., TMath::Pi());
  Detahh = new TH1D ("Detahh","Delta Eta Z jets", nbin_, 0., 5.0);
  Ptminhh = new TH1D ("Ptminhh","Pt of lowest Pt Z jet", nbin_, 0., 200.);
  Ptmaxhh = new TH1D ("Ptmaxhh","Pt of highest Pt Z jet", nbin_, 0., 500.);
  Etaminhh = new TH1D ("Etaminhh","Min |eta| of Z jet", nbin_, 0., 5.);
  Etamaxhh = new TH1D ("Etamaxhh","Max |eta| of Z jet", nbin_, 0., 5.);
  Ptqq = new TH1D ("Ptqq","Pt of tag jet system", nbin_, 0., 500.);
  Mqq = new TH1D ("Mqq","Mass of tag jet system", nbin_, 0., 2000.);
  Etaqq = new TH1D ("Etaqq","|Eta} of tag jet system", nbin_, 0., 5.);
  Pthh = new TH1D ("Pthh","Pt of had Z system", nbin_, 0., 500.);
  Mhh = new TH1D ("Mhh","Mass of had Z system", nbin_, 0., 2000.);
  Etahh = new TH1D ("Etahh","|Eta} of had Z system", nbin_, 0., 5.);  
  Ptll = new TH1D ("Ptll","Pt of lep Z system", nbin_, 0., 500.);
  Mll = new TH1D ("Mll","Mass of lep Z system", nbin_, 0., 2000.);
  Etall = new TH1D ("Etall","|Eta} of lep Z system", nbin_, 0., 5.);  
  DphiTjetZjet = new TH1D("DphiTjetZjet","Delta phi between T and had Z syst", nbin_, 0., TMath::Pi());
  DphiTjetZlep = new TH1D("DphiTjetZlep","Delta phi between T and lep Z syst", nbin_, 0., TMath::Pi());
  DphiminTZ = new TH1D("DphiminTZ","Min DP between T and one Z system", nbin_, 0., TMath::Pi());
  DetaTjetZjet = new TH1D("DetaTjetZjet","Delta eta between T and had Z syst", nbin_, 0., 5.);
  DetaTjetZlep = new TH1D("DetaTjetZlep","Delta eta between T and lep Z syst", nbin_, 0., 5.);
  DetaminTZ = new TH1D("DetaminTZ","Delta eta min between T and one Z system", nbin_, 0., 5.);
  DphiZjetZlep = new TH1D("DphiZjetZlep","Delta phi between had and lep Z syst", nbin_, 0., TMath::Pi());
  DetaZjetZlep = new TH1D("DetaZjetZlep","Delta eta between had and lep Z syst", nbin_, 0., 5.);
  MassTjetZjet = new TH1D("MassTjetZjet","Mass of T and had Z syst", nbin_, 0., 2000.);
  MassTjetZlep = new TH1D("MassTjetZlep","Mass of T and lep Z syst", nbin_, 0., 2000.);
  MassZjetZlep = new TH1D("MassZjetZlep","Mass of ZZ system", nbin_, 0., 2000.);
  MassTZZ = new TH1D("MassTZZ","Mass of TZZ system", nbin_, 0., 2000.);
  EtaTZZ = new TH1D("EtaTZZ","Eta of TZZ system", nbin_, 0., 5.);
  PtTZZ = new TH1D("PtTZZ","Pt of TZZ system", nbin_, 0., 500.);

  F_Dphiqq = new TH1D ("F_Dphiqq","Delta Phi tag jets", nbin_, 0., TMath::Pi());
  F_Detaqq = new TH1D ("F_Detaqq","Delta Eta tag jets", nbin_, 0., 10.0);
  F_Ptminqq = new TH1D ("F_Ptminqq","Pt of lowest Pt tag jet", nbin_, 0., 200.);
  F_Ptmaxqq = new TH1D ("F_Ptmaxqq","Pt of highest Pt tag jet", nbin_, 0., 500.);
  F_Etaminqq = new TH1D ("F_Etaminqq","Min |eta| of tag jet", nbin_, 0., 5.);
  F_Etamaxqq = new TH1D ("F_Etamaxqq","Max |eta| of tag jet", nbin_, 0., 5.);
  F_Dphihh = new TH1D ("F_Dphihh","Delta Phi Z jets", nbin_, 0., TMath::Pi());
  F_Detahh = new TH1D ("F_Detahh","Delta Eta Z jets", nbin_, 0., 5.0);
  F_Ptminhh = new TH1D ("F_Ptminhh","Pt of lowest Pt Z jet", nbin_, 0., 200.);
  F_Ptmaxhh = new TH1D ("F_Ptmaxhh","Pt of highest Pt Z jet", nbin_, 0., 500.);
  F_Etaminhh = new TH1D ("F_Etaminhh","Min |eta| of Z jet", nbin_, 0., 5.);
  F_Etamaxhh = new TH1D ("F_Etamaxhh","Max |eta| of Z jet", nbin_, 0., 5.);
  F_Ptqq = new TH1D ("F_Ptqq","Pt of tag jet system", nbin_, 0., 500.);
  F_Mqq = new TH1D ("F_Mqq","Mass of tag jet system", nbin_, 0., 2000.);
  F_Etaqq = new TH1D ("F_Etaqq","|Eta} of tag jet system", nbin_, 0., 5.);
  F_Pthh = new TH1D ("F_Pthh","Pt of had Z system", nbin_, 0., 500.);
  F_Mhh = new TH1D ("F_Mhh","Mass of had Z system", nbin_, 0., 2000.);
  F_Etahh = new TH1D ("F_Etahh","|Eta} of had Z system", nbin_, 0., 5.);  
  F_Ptll = new TH1D ("F_Ptll","Pt of lep Z system", nbin_, 0., 500.);
  F_Mll = new TH1D ("F_Mll","Mass of lep Z system", nbin_, 0., 2000.);
  F_Etall = new TH1D ("F_Etall","|Eta} of lep Z system", nbin_, 0., 5.);  
  F_DphiTjetZjet = new TH1D("F_DphiTjetZjet","Delta phi between T and had Z syst", nbin_, 0., TMath::Pi());
  F_DphiTjetZlep = new TH1D("F_DphiTjetZlep","Delta phi between T and lep Z syst", nbin_, 0., TMath::Pi());
  F_DphiminTZ = new TH1D("F_DphiminTZ","Min DP between T and one Z system", nbin_, 0., TMath::Pi());
  F_DetaTjetZjet = new TH1D("F_DetaTjetZjet","Delta eta between T and had Z syst", nbin_, 0., 5.);
  F_DetaTjetZlep = new TH1D("F_DetaTjetZlep","Delta eta between T and lep Z syst", nbin_, 0., 5.);
  F_DetaminTZ = new TH1D("F_DetaminTZ","Delta eta min between T and one Z system", nbin_, 0., 5.);
  F_DphiZjetZlep = new TH1D("F_DphiZjetZlep","Delta phi between had and lep Z syst", nbin_, 0., TMath::Pi());
  F_DetaZjetZlep = new TH1D("F_DetaZjetZlep","Delta eta between had and lep Z syst", nbin_, 0., 5.);
  F_MassTjetZjet = new TH1D("F_MassTjetZjet","Mass of T and had Z syst", nbin_, 0., 2000.);
  F_MassTjetZlep = new TH1D("F_MassTjetZlep","Mass of T and lep Z syst", nbin_, 0., 2000.);
  F_MassZjetZlep = new TH1D("F_MassZjetZlep","Mass of ZZ system", nbin_, 0., 2000.);
  F_MassTZZ = new TH1D("F_MassTZZ","Mass of TZZ system", nbin_, 0., 2000.);
  F_EtaTZZ = new TH1D("F_EtaTZZ","Eta of TZZ system", nbin_, 0., 5.);
  F_PtTZZ = new TH1D("F_PtTZZ","Pt of TZZ system", nbin_, 0., 500.);

  // =========================================================

}

// ------------ method called once each job just after ending the event loop  ------------
void 
VBFHZZkinematicsAnalyzer::endJob() {

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

  overlapJetsDeltaEt_->Write();
  overlapJetsDeltaEta_->Write();
  overlapJetsDeltaEtVSDeltaR2_->Write();
  overlapJetsDeltaEtaVSDeltaR2_->Write();
  overlapJetsDeltaEtVSDeltaEta_->Write();

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

  hadronicZrecMass_->Write();
  hadronicZrecMassResolution_->Write();
  hadronicZrecPt_->Write();
  hadronicZrecPtResolution_->Write();
  hadronicZrecCollinearity_->Write();
  hadronicZrecCollinearityResolution_->Write();
  hadronicZrecDeltaPhiVSDeltaEta_->Write();

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
  DeltaEtvsEt->Write();

  // ====================================================
  Dphiqq->Write();
  Detaqq->Write();
  Ptminqq->Write();
  Ptmaxqq->Write();
  Etaminqq->Write();
  Etamaxqq->Write();
  Dphihh->Write();
  Detahh->Write();
  Ptminhh->Write();
  Ptmaxhh->Write();
  Etaminhh->Write();
  Etamaxhh->Write();
  Ptqq->Write();
  Mqq->Write();
  Etaqq->Write();
  Pthh->Write();
  Mhh->Write();
  Etahh->Write();
  Ptll->Write();
  Mll->Write();
  Etall->Write();
  DphiTjetZjet->Write();
  DphiTjetZlep->Write();
  DphiminTZ->Write();
  DetaTjetZjet->Write();
  DetaTjetZlep->Write();
  DetaminTZ->Write();
  MassTjetZjet->Write();
  MassTjetZlep->Write();
  MassZjetZlep->Write();
  MassTZZ->Write();
  EtaTZZ->Write();
  PtTZZ->Write();
 
  F_Dphiqq->Write();
  F_Detaqq->Write();
  F_Ptminqq->Write();
  F_Ptmaxqq->Write();
  F_Etaminqq->Write();
  F_Etamaxqq->Write();
  F_Dphihh->Write();
  F_Detahh->Write();
  F_Ptminhh->Write();
  F_Ptmaxhh->Write();
  F_Etaminhh->Write();
  F_Etamaxhh->Write();
  F_Ptqq->Write();
  F_Mqq->Write();
  F_Etaqq->Write();
  F_Pthh->Write();
  F_Mhh->Write();
  F_Etahh->Write();
  F_Ptll->Write();
  F_Mll->Write();
  F_Etall->Write();
  F_DphiTjetZjet->Write();
  F_DphiTjetZlep->Write();
  F_DphiminTZ->Write();
  F_DetaTjetZjet->Write();
  F_DetaTjetZlep->Write();
  F_DetaminTZ->Write();
  F_MassTjetZjet->Write();
  F_MassTjetZlep->Write();
  F_MassZjetZlep->Write();
  F_MassTZZ->Write();
  F_EtaTZZ->Write();
  F_PtTZZ->Write();
 
  // ====================================================

}
//define this as a plug-in
DEFINE_FWK_MODULE(VBFHZZkinematicsAnalyzer);

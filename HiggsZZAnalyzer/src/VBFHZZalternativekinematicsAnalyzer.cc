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

#include "AnalysisExamples/HiggsZZAnalyzer/interface/VBFHZZalternativekinematicsAnalyzer.h"

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

VBFHZZalternativekinematicsAnalyzer::VBFHZZalternativekinematicsAnalyzer(const edm::ParameterSet& iConfig) :
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
  tmvaEtaCutSuffix_(   iConfig.getUntrackedParameter<std::string>("tmvaEtaCutSuffix") ),
  tmvaMZcutSuffix_(   iConfig.getUntrackedParameter<std::string>("tmvaMZcutSuffix") ),
  mHres_(             iConfig.getUntrackedParameter<double>("mHres") ),
  mH_(                iConfig.getUntrackedParameter<double>("mH") ),
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
  // signal-like w/ eta cut event TTree
  tmvaEtaCutTreeWriterPtr_.reset(new TMVAtreeWriter<Float_t>(eventVariablesNamesVector_,tmvaEtaCutSuffix_));
  // signal-like w/ closest mjj to mZ cut event TTree
  tmvaMZcutTreeWriterPtr_.reset(new TMVAtreeWriter<Float_t>(eventVariablesNamesVector_,tmvaMZcutSuffix_));


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

  lljjEventCounter_      = 0;
  tmvaMatchingCounter_   = 0;
  combinatorialCounter_  = 0;
  etaCutMatchingCounter_ = 0;
  mZMatchingCounter_     = 0;

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


VBFHZZalternativekinematicsAnalyzer::~VBFHZZalternativekinematicsAnalyzer()
{
   std::cout << "endJob" << std::endl;

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
VBFHZZalternativekinematicsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

    if (hadronicZcandidates.size()==2 && leptonicZcandidates.size()==2) lljjEventCounter_++;

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
		  double deltaR2 = closestJetDeltaR2Vec[minPartonIndex][0];
		  double deltaEt  = partonsCandidates[minPartonIndex]->et()-aboveEtCutJetVec[jetIndex].et();
		  double deltaEta = partonsCandidates[minPartonIndex]->eta()-aboveEtCutJetVec[jetIndex].eta();

		  deltaR2ZjetsZpartons_->Fill(deltaR2);
		  jetParton_deltaEt_->Fill(deltaEt);
		  jetParton_deltaEta_->Fill(deltaEta);
		  jetParton_deltaEtVSdeltaR2_->Fill( deltaR2,deltaEt );
		  jetParton_deltaEtaVSdeltaR2_->Fill(deltaR2,deltaEta);
		  jetParton_deltaEtVSdeltaR2_profile_->Fill( deltaR2,deltaEt );
		  jetParton_deltaEtaVSdeltaR2_profile_->Fill(deltaR2,deltaEta);
		  if (partonsCandidates[minPartonIndex]->et()!=0.) jetParton_deltaEtRes_->Fill(deltaEt/partonsCandidates[minPartonIndex]->et());
		  //		  std::cout << "Z jetIndex: " <<  jetIndex << std::endl;
		}
		else if ( minPartonIndex == 2 || minPartonIndex == 3 ) {
		  TjetVec.push_back(aboveEtCutJetVec[jetIndex]);
		  TjetInd.push_back(jetIndex);
		  deltaR2TjetsTpartons_->Fill(closestJetDeltaR2Vec[minPartonIndex][0]);
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


	      double hadronicZrecMass = Zjet->mass();
	      double hadronicZMass    = hadronicZ->mass();
	      double hadronicZrecMassResolution = -99.;
	      if ( hadronicZMass != 0. ) hadronicZrecMassResolution = (hadronicZrecMass-hadronicZMass)/hadronicZMass;
	      double hadronicZrecPt = Zjet->pt();
	      double hadronicZPt    = hadronicZ->pt();
	      double hadronicZrecPtResolution = -99.;
	      if ( hadronicZPt != 0. ) hadronicZrecPtResolution = (hadronicZrecPt-hadronicZPt)/hadronicZPt;
	      double hadronicZDeltaR = TMath::Sqrt(pow(DeltaPhi(hadronicZcandidates[0]->phi(),hadronicZcandidates[1]->phi()),2)+pow(hadronicZcandidates[0]->eta()-hadronicZcandidates[1]->eta(),2));
	      double hadronicZrecDeltaR = TMath::Sqrt(pow(DeltaPhi(ZjetVec[0].phi(),ZjetVec[1].phi()),2)+pow(ZjetVec[0].eta()-ZjetVec[1].eta(),2));
	      hadronicZrecMass_->Fill(hadronicZrecMass);
	      hadronicZMass_->Fill(hadronicZMass);
	      hadronicZrecMassResolution_->Fill(hadronicZrecMassResolution);
	      hadronicZrecPt_->Fill(hadronicZrecPt);
	      hadronicZPt_->Fill(hadronicZPt);
	      hadronicZrecPtResolution_->Fill(hadronicZrecPtResolution);

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

		tmvaMatchingCounter_++;

		// Study alternativekinematics of true ZZqq system
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

		// Loop on true combination and all false ones, and derive alternativekinematics of configuration
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

			// Study alternativekinematics of random ZZqq systems
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

	      

		// count combination due to jet multiplicity
		jet1 = aboveEtCutJetVec.begin();
		for ( ; jet1 != aboveEtCutJetVec.end()-1; ++jet1, i1++ ) {
		  std::vector<OfflineJet>::const_iterator jet2 = jet1+1;
		  for ( ; jet2 != aboveEtCutJetVec.end(); ++jet2, i2++) {
		    std::vector<OfflineJet>::const_iterator jet3 = aboveEtCutJetVec.begin();
		    for ( ; jet3 != aboveEtCutJetVec.end()-1; ++jet3, i3++) {
		      if ( jet3==jet1 || jet3==jet2 ) continue;
		      std::vector<OfflineJet>::const_iterator jet4 = jet3+1;
		      for ( ; jet4 != aboveEtCutJetVec.end(); ++jet4, i4++) {
			if ( jet4==jet1 || jet4==jet2 ) continue;
			combinatorialCounter_++;
		      }
		    }
		  }
		}
	      }

	      // =============================================================================================
	      // derive alternativekinematics of configuration: eta cut
	      // ------------------------------------------------------

	      // eT sorting
	      // ----------
	      std::vector<OfflineJet> etSortedJetVec;
	      std::vector<int>        etSortedJetIndexVec;
	      std::vector<bool> stored;
	      
	      for ( unsigned int times = 0; times < aboveEtCutJetVec.size(); times++ ) 
		stored.push_back(kFALSE);
	      
	      for ( unsigned int times = 0; times < aboveEtCutJetVec.size(); times++ ) {
		double maxEt = 0.;
		unsigned int maxIndex = 0;
		std::vector<OfflineJet>::const_iterator aboveEtCutJetVec_itr = aboveEtCutJetVec.begin(); 
		unsigned int index = 0;
		for ( ; aboveEtCutJetVec_itr != aboveEtCutJetVec.end(); ++aboveEtCutJetVec_itr, index++ ) {
		  if ( stored[index] ) continue;
		  if ( aboveEtCutJetVec_itr->et() >= maxEt ) {
		    maxEt = aboveEtCutJetVec_itr->et();
		    maxIndex = index;
		  }
		}
		stored[maxIndex] = kTRUE;
		etSortedJetVec.push_back(aboveEtCutJetVec[maxIndex]);
		etSortedJetIndexVec.push_back(maxIndex);
	      }
	      
	      //	      for ( unsigned int times = 0; times < aboveEtCutJetVec.size(); times++ ) {
	      //		std::cout << "aboveEtCutJetVec[" << times << "]: " << aboveEtCutJetVec[times].et()
	      //			  << " --> etSortedJetVec[" << times << "]: " << etSortedJetVec[times].et() << std::endl;
	      //		std::cout << "etSortedJetIndexVec[" << times << "]: " << etSortedJetIndexVec[times] << " -->eta: " << etSortedJetVec[times].eta() << std::endl;
	      //	      }
	      
	      std::vector<OfflineJet> centralJetsVector;
	      std::vector<int>        centralJetsIndexVector;
	      std::vector<OfflineJet> forwardJetsVector;
	      std::vector<int>        forwardJetsIndexVector;
	      std::vector<OfflineJet>::const_iterator etSortedJetVec_itr      = etSortedJetVec.begin();
	      std::vector<int>::const_iterator        etSortedJetIndexVec_itr = etSortedJetIndexVec.begin();
	      for ( ; etSortedJetVec_itr != etSortedJetVec.end(); ++etSortedJetVec_itr, ++etSortedJetIndexVec_itr ) {
		if ( fabs( etSortedJetVec_itr->eta() ) <= 2. ) {
		  centralJetsVector.push_back(*etSortedJetVec_itr);
		  centralJetsIndexVector.push_back(*etSortedJetIndexVec_itr);
		}
		else {
		  forwardJetsVector.push_back(*etSortedJetVec_itr);
		  forwardJetsIndexVector.push_back(*etSortedJetIndexVec_itr);
		}
	      }

	      if ( centralJetsVector.size() >= 2 && forwardJetsVector.size() >= 2 ) {

		//		for (unsigned int index = 0; index < centralJetsVector.size(); index++ ) 
		//		  std::cout << "centralJetsIndexVector[" << index << "]: " << centralJetsIndexVector[index] << ": et: " << centralJetsVector[index].et() << std::endl;
		//		for (unsigned int index = 0 ; index < forwardJetsIndexVector.size(); index++ )
		//		  std::cout << "forwardJetsIndexVector[" << index << "]: " << forwardJetsIndexVector[index] << ": et: " << forwardJetsVector[index].et() << std::endl;		

		std::vector<OfflineJet> TjetsVector;
		std::vector<int>        TjetsIndexVector;
		TjetsVector.push_back(forwardJetsVector[0]);
		TjetsIndexVector.push_back(forwardJetsIndexVector[0]);
		double eta1 = forwardJetsVector[0].eta();
		std::vector<OfflineJet>::const_iterator forwardJetsVector_itr2      = forwardJetsVector.begin()+1;
		std::vector<int>::const_iterator        forwardJetsIndexVector_itr2 = forwardJetsIndexVector.begin()+1;
		for ( ; forwardJetsVector_itr2 != forwardJetsVector.end(); ++forwardJetsVector_itr2, ++forwardJetsIndexVector_itr2 ){
		  if ( eta1*forwardJetsVector_itr2->eta() < 0 ) {
		    TjetsVector.push_back(*forwardJetsVector_itr2);
		    TjetsIndexVector.push_back(*forwardJetsIndexVector_itr2);
		    break;
		  }		      
		}
		//		for (unsigned int index = 0; index < TjetsIndexVector.size(); index++ )
		//		  std::cout << "TjetsIndexVector[" << index << "]: " << TjetsIndexVector[index] << std::endl;

		if ( TjetsVector.size() == 2 ) {
		  std::vector<OfflineJet> ZjetsVector;
		  std::vector<int>        ZjetsIndexVector;
		  for (unsigned int index = 0; index < 2; index++) {
		    ZjetsVector.push_back(centralJetsVector[index]);
		    ZjetsIndexVector.push_back(centralJetsIndexVector[index]);
		  }

		  if ( ZjetsVector.size() == 2 ) {
		    if (
			( (TjetsIndexVector[0]==TjetInd[0] && TjetsIndexVector[1]==TjetInd[1] ) ||
			  (TjetsIndexVector[1]==TjetInd[0] && TjetsIndexVector[0]==TjetInd[1] ) ) &&
			( (ZjetsIndexVector[0]==ZjetInd[0] && ZjetsIndexVector[1]==ZjetInd[1] ) ||
			  (ZjetsIndexVector[1]==ZjetInd[0] && ZjetsIndexVector[0]==ZjetInd[1] ) )
			) etaCutMatchingCounter_++;
		  
		    Particle *Tjet =  new Particle (0,TjetsVector[0].p4()+TjetsVector[1].p4(), TjetsVector[0].vertex(),0,0,true);
		    Particle *Zjet =  new Particle (0,ZjetsVector[0].p4()+ZjetsVector[1].p4(), ZjetsVector[0].vertex(),0,0,true);
		  
		    // Study alternativekinematics of ZZqq systems
		    // -------------------------------------------
		    // qq: tag jets
		    // ------------
		    double etaCut_dphiqq = DeltaPhi(TjetsVector[0].phi(),TjetsVector[1].phi());
		    double etaCut_detaqq = fabs(TjetsVector[0].eta()-TjetsVector[1].eta());
		    double etaCut_ptminqq = TjetsVector[0].et();
		    double etaCut_ptmaxqq = TjetsVector[1].et();
		    if ( TjetsVector[1].et()<etaCut_ptminqq ) {
		      etaCut_ptminqq = TjetsVector[1].et();
		      etaCut_ptmaxqq = TjetsVector[0].et();
		    }
		    double etaCut_etaminqq = fabs(TjetsVector[0].eta());
		    double etaCut_etamaxqq = fabs(TjetsVector[1].eta());
		    if ( fabs(TjetsVector[1].eta())<etaCut_etaminqq ) {
		      etaCut_etaminqq = fabs(TjetsVector[1].eta());
		      etaCut_etamaxqq = fabs(TjetsVector[0].eta());
		    }
		    ETACUT_Dphiqq->Fill(etaCut_dphiqq);
		    ETACUT_Detaqq->Fill(etaCut_detaqq);
		    ETACUT_Ptminqq->Fill(etaCut_ptminqq);
		    ETACUT_Ptmaxqq->Fill(etaCut_ptmaxqq);
		    ETACUT_Etaminqq->Fill(etaCut_etaminqq);
		    ETACUT_Etamaxqq->Fill(etaCut_etamaxqq);
		    eventEtaCutVariablesVector_.push_back(TjetsVector[0].et());
		    eventEtaCutVariablesVector_.push_back(TjetsVector[1].et());
		    eventEtaCutVariablesVector_.push_back(TjetsVector[0].eta());
		    eventEtaCutVariablesVector_.push_back(TjetsVector[1].eta());
		    eventEtaCutVariablesVector_.push_back(TjetsVector[0].phi());
		    eventEtaCutVariablesVector_.push_back(TjetsVector[1].phi());
		    eventEtaCutVariablesVector_.push_back(etaCut_dphiqq);
		    eventEtaCutVariablesVector_.push_back(etaCut_detaqq);
		    eventEtaCutVariablesVector_.push_back(etaCut_ptminqq);
		    eventEtaCutVariablesVector_.push_back(etaCut_ptmaxqq);
		    eventEtaCutVariablesVector_.push_back(etaCut_etaminqq);
		    eventEtaCutVariablesVector_.push_back(etaCut_etamaxqq);
		    
		    // hh: Z jets
		    // ----------
		    double etaCut_dphihh = DeltaPhi(ZjetsVector[0].phi(),ZjetsVector[1].phi());
		    double etaCut_detahh = fabs(ZjetsVector[0].eta()-ZjetsVector[1].eta());
		    double etaCut_ptminhh = ZjetsVector[0].et();
		    double etaCut_ptmaxhh = ZjetsVector[1].et();
		    if ( ZjetsVector[1].et()<etaCut_ptminhh ) {
		      etaCut_ptminhh = ZjetsVector[1].et();
		      etaCut_ptmaxhh = ZjetsVector[0].et();
		    }
		    double etaCut_etaminhh = fabs(ZjetsVector[0].eta());
		    double etaCut_etamaxhh = fabs(ZjetsVector[1].eta());
		    if ( fabs(ZjetsVector[1].eta())<etaCut_etaminhh ) {
		      etaCut_etaminhh = fabs(ZjetsVector[1].eta());
		      etaCut_etamaxhh = fabs(ZjetsVector[0].eta());
		    }
		    ETACUT_Dphihh->Fill(etaCut_dphihh);
		    ETACUT_Detahh->Fill(etaCut_detahh);
		    ETACUT_Ptminhh->Fill(etaCut_ptminhh);
		    ETACUT_Ptmaxhh->Fill(etaCut_ptmaxhh);
		    ETACUT_Etaminhh->Fill(etaCut_etaminhh);
		    ETACUT_Etamaxhh->Fill(etaCut_etamaxhh);
		    
		    eventEtaCutVariablesVector_.push_back(ZjetsVector[0].et());
		    eventEtaCutVariablesVector_.push_back(ZjetsVector[1].et());
		    eventEtaCutVariablesVector_.push_back(ZjetsVector[0].eta());
		    eventEtaCutVariablesVector_.push_back(ZjetsVector[1].eta());
		    eventEtaCutVariablesVector_.push_back(ZjetsVector[0].phi());
		    eventEtaCutVariablesVector_.push_back(ZjetsVector[1].phi());
		    eventEtaCutVariablesVector_.push_back(etaCut_dphihh);
		    eventEtaCutVariablesVector_.push_back(etaCut_detahh);
		    eventEtaCutVariablesVector_.push_back(etaCut_ptminhh);
		    eventEtaCutVariablesVector_.push_back(etaCut_ptmaxhh);
		    eventEtaCutVariablesVector_.push_back(etaCut_etaminhh);
		    eventEtaCutVariablesVector_.push_back(etaCut_etamaxhh);
		    
		    // qq system
		    // ---------
		    double etaCut_ptqq = (Tjet->p4()).pt();
		    double etaCut_mqq = (Tjet->p4()).mass();
		    double etaCut_etaqq = fabs((Tjet->p4()).eta());
		    // hh systems
		    // ---------
		    double etaCut_pthh = (Zjet->p4()).pt();
		    double etaCut_mhh = (Zjet->p4()).mass();
		    double etaCut_etahh = fabs((Zjet->p4()).eta());
		    // ll systems
		    // ---------
		    double etaCut_ptll = (Zlep->p4()).pt();
		    double etaCut_mll = Zlep->mass();
		    double etaCut_etall = fabs(Zlep->eta());
		    ETACUT_Ptqq->Fill(etaCut_ptqq);
		    ETACUT_Mqq->Fill(etaCut_mqq);
		    ETACUT_Etaqq->Fill(etaCut_etaqq);
		    ETACUT_Pthh->Fill(etaCut_pthh);
		    ETACUT_Mhh->Fill(etaCut_mhh);
		    ETACUT_Etahh->Fill(etaCut_etahh);
		    ETACUT_Ptll->Fill(etaCut_ptll);
		    ETACUT_Mll->Fill(etaCut_mll);
		    ETACUT_Etall->Fill(etaCut_etall);
		  
		    eventEtaCutVariablesVector_.push_back(etaCut_ptqq);
		    eventEtaCutVariablesVector_.push_back(etaCut_mqq);
		    eventEtaCutVariablesVector_.push_back(etaCut_etaqq);
		    eventEtaCutVariablesVector_.push_back(etaCut_pthh);
		    eventEtaCutVariablesVector_.push_back(etaCut_mhh);
		    eventEtaCutVariablesVector_.push_back(etaCut_etahh);
		    eventEtaCutVariablesVector_.push_back(etaCut_ptll);
		    eventEtaCutVariablesVector_.push_back(etaCut_mll);
		    eventEtaCutVariablesVector_.push_back(etaCut_etall);
		    // 2-particle vars
		    // ---------------
		    double etaCut_dphiTjetZjet = DeltaPhi((Tjet->p4()).phi(),(Zjet->p4()).phi());
		    double etaCut_dphiTjetZlep = DeltaPhi((Tjet->p4()).phi(),Zlep->phi());
		    double etaCut_dphiminTZ = etaCut_dphiTjetZjet;
		    if ( etaCut_dphiminTZ>etaCut_dphiTjetZlep ) etaCut_dphiminTZ = etaCut_dphiTjetZlep;
		    double etaCut_detaTjetZjet = fabs((Tjet->p4()).eta()-(Zjet->p4()).eta());
		    double etaCut_detaTjetZlep = fabs((Tjet->p4()).eta()-Zlep->eta());
		    double etaCut_detaminTZ = etaCut_detaTjetZjet;
		    if ( etaCut_detaminTZ>etaCut_detaTjetZlep ) etaCut_detaminTZ = etaCut_detaTjetZlep;
		    double etaCut_dphiZjetZlep = DeltaPhi((Zjet->p4()).phi(),Zlep->phi());
		    double etaCut_detaZjetZlep = fabs((Zjet->p4()).eta()-Zlep->eta());
		    double etaCut_massTjetZjet = (Tjet->p4()+Zjet->p4()).mass();
		    double etaCut_massTjetZlep = (Tjet->p4()+Zlep->p4()).mass();
		    double etaCut_massZjetZlep = (Zjet->p4()+Zlep->p4()).mass();
		    ETACUT_DphiTjetZjet->Fill(etaCut_dphiTjetZjet);
		    ETACUT_DphiTjetZlep->Fill(etaCut_dphiTjetZlep);
		    ETACUT_DphiminTZ->Fill(etaCut_dphiminTZ);
		    ETACUT_DetaTjetZjet->Fill(etaCut_detaTjetZjet);
		    ETACUT_DetaTjetZlep->Fill(etaCut_detaTjetZlep);
		    ETACUT_DetaminTZ->Fill(etaCut_detaminTZ);
		    ETACUT_DphiZjetZlep->Fill(etaCut_dphiZjetZlep);
		    ETACUT_DetaZjetZlep->Fill(etaCut_detaZjetZlep);
		    ETACUT_MassTjetZjet->Fill(etaCut_massTjetZjet);
		    ETACUT_MassTjetZlep->Fill(etaCut_massTjetZlep);
		    ETACUT_MassZjetZlep->Fill(etaCut_massZjetZlep);
		    
		    eventEtaCutVariablesVector_.push_back(etaCut_dphiTjetZjet);
		    eventEtaCutVariablesVector_.push_back(etaCut_dphiTjetZlep);
		    eventEtaCutVariablesVector_.push_back(etaCut_dphiminTZ);
		    eventEtaCutVariablesVector_.push_back(etaCut_detaTjetZjet);
		    eventEtaCutVariablesVector_.push_back(etaCut_detaTjetZlep);
		    eventEtaCutVariablesVector_.push_back(etaCut_detaminTZ);
		    eventEtaCutVariablesVector_.push_back(etaCut_dphiZjetZlep);
		    eventEtaCutVariablesVector_.push_back(etaCut_detaZjetZlep);
		    eventEtaCutVariablesVector_.push_back(etaCut_massTjetZjet);
		    eventEtaCutVariablesVector_.push_back(etaCut_massTjetZlep);
		    eventEtaCutVariablesVector_.push_back(etaCut_massZjetZlep);
		    // 3-particle vars
		    // ---------------
		    double etaCut_massTZZ = (Tjet->p4()+Zjet->p4()+Zlep->p4()).mass();
		    double etaCut_etaTZZ = fabs((Tjet->p4()+Zjet->p4()+Zlep->p4()).eta());
		    double etaCut_ptTZZ = (Tjet->p4()+Zjet->p4()+Zlep->p4()).pt();
		    ETACUT_MassTZZ->Fill(etaCut_massTZZ);
		    ETACUT_EtaTZZ->Fill(etaCut_etaTZZ);
		    ETACUT_PtTZZ->Fill(etaCut_ptTZZ);
		    
		    eventEtaCutVariablesVector_.push_back(etaCut_massTZZ);
		    eventEtaCutVariablesVector_.push_back(etaCut_etaTZZ);
		    eventEtaCutVariablesVector_.push_back(etaCut_ptTZZ);
		  
		    // Stored the event number
		    eventEtaCutVariablesVector_.push_back(eventcounter_);
		    
		    if( writeTMVA_ ) {
		      // Fill the tree for the TMVA
		      tmvaEtaCutTreeWriterPtr_->fill(eventEtaCutVariablesVector_);
		    }
		    /*
		      else {
		      // Fill the values in the array used by the reader
		      vector<Float_t>::const_iterator varIt = eventSignalVariablesVector_.begin();
		      int iVarBackground = 0;
		      for( ; varIt != eventEtaCutVariablesVector_.end(); ++varIt, ++iVarBackground ) {
		      variables_[iVarBackground] = *varIt;
		      }
		      discriminantBackground.push_back( reader_->EvaluateMVA( "BDT method" ) );
		      }
		    */
		    eventEtaCutVariablesVector_.clear();
		  }
	      
		  ZjetsVector.clear();
		  ZjetsIndexVector.clear();
		  // derive alternativekinematics of configuration: closest mjj to mZ
		  // ----------------------------------------------------------------
		  double mZres = 20.;
		  double mHres = mHres_;
		  unsigned int ZjetIndex1 = 0;
		  unsigned int ZjetIndex2 = 1;
		  double tmp_mZres = mZres;
		  double tmp_chi2 = 1000.;
		  std::vector<OfflineJet>::const_iterator centralJetsVector_itr1      = centralJetsVector.begin();
		  std::vector<int>::const_iterator        centralJetsIndexVector_itr1 = centralJetsIndexVector.begin();
		  unsigned int index1 = 0;
		  for ( ; centralJetsVector_itr1 != centralJetsVector.end(); ++centralJetsVector_itr1, ++centralJetsIndexVector_itr1, index1++ ) {
		    std::vector<OfflineJet>::const_iterator centralJetsVector_itr2      = centralJetsVector_itr1+1;
		    std::vector<int>::const_iterator        centralJetsIndexVector_itr2 = centralJetsIndexVector_itr1+1;
		    unsigned int index2 = index1+1;
		    for ( ; centralJetsVector_itr2 != centralJetsVector.end(); ++centralJetsVector_itr2, ++centralJetsIndexVector_itr2, index2++ ) {
		      double jjInvMass = (centralJetsVector_itr1->p4()+centralJetsVector_itr2->p4()).mass();
		      //		      std::cout << "jjInvMass[" << *centralJetsIndexVector_itr1 << ";" << *centralJetsIndexVector_itr2 << "]: " << jjInvMass << std::endl;
		      if ( fabs(jjInvMass-ZMass_) <= tmp_mZres ) {
			tmp_mZres = fabs(jjInvMass-ZMass_);
			ZjetIndex1 = index1;
			ZjetIndex2 = index2;
		      }
		      double jjInvMassChi2 = pow((jjInvMass-ZMass_)/mZres,2);
		      double lljjIinvMass = (centralJetsVector_itr1->p4()+centralJetsVector_itr2->p4()+Zlep->p4()).mass();
		      double lljjIinvMassChi2 = pow((lljjIinvMass-mH_)/mHres,2);
		      if ( jjInvMassChi2+lljjIinvMassChi2 < tmp_chi2 ) {
			tmp_chi2 = jjInvMassChi2+lljjIinvMassChi2;
			ZjetIndex1 = index1;
			ZjetIndex2 = index2;
			//			std::cout << "index1: " << *centralJetsIndexVector_itr1 << " <--> index2: " << *centralJetsIndexVector_itr2 << std::endl;
		      }
		    }
		  }
		
		  ZjetsVector.push_back(centralJetsVector[ZjetIndex1]);
		  ZjetsVector.push_back(centralJetsVector[ZjetIndex2]);
		  ZjetsIndexVector.push_back(centralJetsIndexVector[ZjetIndex1]);
		  ZjetsIndexVector.push_back(centralJetsIndexVector[ZjetIndex2]);
		    
//		  for (unsigned int index = 0; index < 2; index++) {
//		    std::cout << "ZjetsIndexVector[" << index << "]: " << ZjetsIndexVector[index] 
//			      << " <--> ZjetInd[" << index << "]: " << ZjetInd[index] << std::endl;
//		    std::cout << "TjetsIndexVector[" << index << "]: " << TjetsIndexVector[index]
//			      << " <--> TjetInd[" << index << "]: " << TjetInd[index] << std::endl;
//		  }
		
		  if ( ZjetsVector.size() == 2 ) {
		    if (
			( (TjetsIndexVector[0]==TjetInd[0] && TjetsIndexVector[1]==TjetInd[1] ) ||
			  (TjetsIndexVector[1]==TjetInd[0] && TjetsIndexVector[0]==TjetInd[1] ) ) &&
			( (ZjetsIndexVector[0]==ZjetInd[0] && ZjetsIndexVector[1]==ZjetInd[1] ) ||
			  (ZjetsIndexVector[1]==ZjetInd[0] && ZjetsIndexVector[0]==ZjetInd[1] ) ) ) mZMatchingCounter_++;

		    Particle *Tjet =  new Particle (0,TjetsVector[0].p4()+TjetsVector[1].p4(), TjetsVector[0].vertex(),0,0,true);
		    Particle *Zjet =  new Particle (0,ZjetsVector[0].p4()+ZjetsVector[1].p4(), ZjetsVector[0].vertex(),0,0,true);
		      


		    // Study alternativekinematics of random ZZqq systems
		    // ---------------------------------------
		    // qq: tag jets
		    // ------------
		    double mZcut_dphiqq = DeltaPhi(TjetsVector[0].phi(),TjetsVector[1].phi());
		    double mZcut_detaqq = fabs(TjetsVector[0].eta()-TjetsVector[1].eta());
		    double mZcut_ptminqq = TjetsVector[0].et();
		    double mZcut_ptmaxqq = TjetsVector[1].et();
		    if ( TjetsVector[1].et()<mZcut_ptminqq ) {
		      mZcut_ptminqq = TjetsVector[1].et();
		      mZcut_ptmaxqq = TjetsVector[2].et();
		    }
		    double mZcut_etaminqq = fabs(TjetsVector[0].eta());
		    double mZcut_etamaxqq = fabs(TjetsVector[1].eta());
		    if ( fabs(TjetsVector[1].eta())<mZcut_etaminqq ) {
		      mZcut_etaminqq = fabs(TjetsVector[1].eta());
		      mZcut_etamaxqq = fabs(TjetsVector[0].eta());
		    }
		    MZCUT_Dphiqq->Fill(mZcut_dphiqq);
		    MZCUT_Detaqq->Fill(mZcut_detaqq);
		    MZCUT_Ptminqq->Fill(mZcut_ptminqq);
		    MZCUT_Ptmaxqq->Fill(mZcut_ptmaxqq);
		    MZCUT_Etaminqq->Fill(mZcut_etaminqq);
		    MZCUT_Etamaxqq->Fill(mZcut_etamaxqq);
		    eventMZcutVariablesVector_.push_back(TjetsVector[0].et());
		    eventMZcutVariablesVector_.push_back(TjetsVector[1].et());
		    eventMZcutVariablesVector_.push_back(TjetsVector[0].eta());
		    eventMZcutVariablesVector_.push_back(TjetsVector[1].eta());
		    eventMZcutVariablesVector_.push_back(TjetsVector[0].phi());
		    eventMZcutVariablesVector_.push_back(TjetsVector[1].phi());
		    eventMZcutVariablesVector_.push_back(mZcut_dphiqq);
		    eventMZcutVariablesVector_.push_back(mZcut_detaqq);
		    eventMZcutVariablesVector_.push_back(mZcut_ptminqq);
		    eventMZcutVariablesVector_.push_back(mZcut_ptmaxqq);
		    eventMZcutVariablesVector_.push_back(mZcut_etaminqq);
		    eventMZcutVariablesVector_.push_back(mZcut_etamaxqq);

		    // hh: Z jets
		    // ----------
		    double mZcut_dphihh = DeltaPhi(ZjetsVector[0].phi(),ZjetsVector[1].phi());
		    double mZcut_detahh = fabs(ZjetsVector[0].eta()-ZjetsVector[1].eta());
		    double mZcut_ptminhh = ZjetsVector[0].et();
		    double mZcut_ptmaxhh = ZjetsVector[1].et();
		    if ( ZjetsVector[1].et()<mZcut_ptminhh ) {
		      mZcut_ptminhh = ZjetsVector[1].et();
		      mZcut_ptmaxhh = ZjetsVector[0].et();
		    }
		    double mZcut_etaminhh = fabs(ZjetsVector[0].eta());
		    double mZcut_etamaxhh = fabs(ZjetsVector[1].eta());
		    if ( fabs(ZjetsVector[1].eta())<mZcut_etaminhh ) {
		      mZcut_etaminhh = fabs(ZjetsVector[1].eta());
		      mZcut_etamaxhh = fabs(ZjetsVector[0].eta());
		    }
		    MZCUT_Dphihh->Fill(mZcut_dphihh);
		    MZCUT_Detahh->Fill(mZcut_detahh);
		    MZCUT_Ptminhh->Fill(mZcut_ptminhh);
		    MZCUT_Ptmaxhh->Fill(mZcut_ptmaxhh);
		    MZCUT_Etaminhh->Fill(mZcut_etaminhh);
		    MZCUT_Etamaxhh->Fill(mZcut_etamaxhh);

		    eventMZcutVariablesVector_.push_back(ZjetsVector[0].et());
		    eventMZcutVariablesVector_.push_back(ZjetsVector[1].et());
		    eventMZcutVariablesVector_.push_back(ZjetsVector[0].eta());
		    eventMZcutVariablesVector_.push_back(ZjetsVector[1].eta());
		    eventMZcutVariablesVector_.push_back(ZjetsVector[0].phi());
		    eventMZcutVariablesVector_.push_back(ZjetsVector[1].phi());
		    eventMZcutVariablesVector_.push_back(mZcut_dphihh);
		    eventMZcutVariablesVector_.push_back(mZcut_detahh);
		    eventMZcutVariablesVector_.push_back(mZcut_ptminhh);
		    eventMZcutVariablesVector_.push_back(mZcut_ptmaxhh);
		    eventMZcutVariablesVector_.push_back(mZcut_etaminhh);
		    eventMZcutVariablesVector_.push_back(mZcut_etamaxhh);

		    // qq system
		    // ---------
		    double mZcut_ptqq = (Tjet->p4()).pt();
		    double mZcut_mqq = (Tjet->p4()).mass();
		    double mZcut_etaqq = fabs((Tjet->p4()).eta());
		    // hh system
		    // ---------
		    double mZcut_pthh = (Zjet->p4()).pt();
		    double mZcut_mhh = (Zjet->p4()).mass();
		    double mZcut_etahh = fabs((Zjet->p4()).eta());
		    // ll system
		    // ---------
		    double mZcut_ptll = (Zlep->p4()).pt();
		    double mZcut_mll = Zlep->mass();
		    double mZcut_etall = fabs(Zlep->eta());
		    MZCUT_Ptqq->Fill(mZcut_ptqq);
		    MZCUT_Mqq->Fill(mZcut_mqq);
		    MZCUT_Etaqq->Fill(mZcut_etaqq);
		    MZCUT_Pthh->Fill(mZcut_pthh);
		    MZCUT_Mhh->Fill(mZcut_mhh);
		    MZCUT_Etahh->Fill(mZcut_etahh);
		    MZCUT_Ptll->Fill(mZcut_ptll);
		    MZCUT_Mll->Fill(mZcut_mll);
		    MZCUT_Etall->Fill(mZcut_etall);

		    eventMZcutVariablesVector_.push_back(mZcut_ptqq);
		    eventMZcutVariablesVector_.push_back(mZcut_mqq);
		    eventMZcutVariablesVector_.push_back(mZcut_etaqq);
		    eventMZcutVariablesVector_.push_back(mZcut_pthh);
		    eventMZcutVariablesVector_.push_back(mZcut_mhh);
		    eventMZcutVariablesVector_.push_back(mZcut_etahh);
		    eventMZcutVariablesVector_.push_back(mZcut_ptll);
		    eventMZcutVariablesVector_.push_back(mZcut_mll);
		    eventMZcutVariablesVector_.push_back(mZcut_etall);
		    // 2-particle vars
		    // ---------------
		    double mZcut_dphiTjetZjet = DeltaPhi((Tjet->p4()).phi(),(Zjet->p4()).phi());
		    double mZcut_dphiTjetZlep = DeltaPhi((Tjet->p4()).phi(),Zlep->phi());
		    double mZcut_dphiminTZ = mZcut_dphiTjetZjet;
		    if ( mZcut_dphiminTZ>mZcut_dphiTjetZlep ) mZcut_dphiminTZ = mZcut_dphiTjetZlep;
		    double mZcut_detaTjetZjet = fabs((Tjet->p4()).eta()-(Zjet->p4()).eta());
		    double mZcut_detaTjetZlep = fabs((Tjet->p4()).eta()-Zlep->eta());
		    double mZcut_detaminTZ = mZcut_detaTjetZjet;
		    if ( mZcut_detaminTZ>mZcut_detaTjetZlep ) mZcut_detaminTZ = mZcut_detaTjetZlep;
		    double mZcut_dphiZjetZlep = DeltaPhi((Zjet->p4()).phi(),Zlep->phi());
		    double mZcut_detaZjetZlep = fabs((Zjet->p4()).eta()-Zlep->eta());
		    double mZcut_massTjetZjet = (Tjet->p4()+Zjet->p4()).mass();
		    double mZcut_massTjetZlep = (Tjet->p4()+Zlep->p4()).mass();
		    double mZcut_massZjetZlep = (Zjet->p4()+Zlep->p4()).mass();
		    MZCUT_DphiTjetZjet->Fill(mZcut_dphiTjetZjet);
		    MZCUT_DphiTjetZlep->Fill(mZcut_dphiTjetZlep);
		    MZCUT_DphiminTZ->Fill(mZcut_dphiminTZ);
		    MZCUT_DetaTjetZjet->Fill(mZcut_detaTjetZjet);
		    MZCUT_DetaTjetZlep->Fill(mZcut_detaTjetZlep);
		    MZCUT_DetaminTZ->Fill(mZcut_detaminTZ);
		    MZCUT_DphiZjetZlep->Fill(mZcut_dphiZjetZlep);
		    MZCUT_DetaZjetZlep->Fill(mZcut_detaZjetZlep);
		    MZCUT_MassTjetZjet->Fill(mZcut_massTjetZjet);
		    MZCUT_MassTjetZlep->Fill(mZcut_massTjetZlep);
		    MZCUT_MassZjetZlep->Fill(mZcut_massZjetZlep);

		    eventMZcutVariablesVector_.push_back(mZcut_dphiTjetZjet);
		    eventMZcutVariablesVector_.push_back(mZcut_dphiTjetZlep);
		    eventMZcutVariablesVector_.push_back(mZcut_dphiminTZ);
		    eventMZcutVariablesVector_.push_back(mZcut_detaTjetZjet);
		    eventMZcutVariablesVector_.push_back(mZcut_detaTjetZlep);
		    eventMZcutVariablesVector_.push_back(mZcut_detaminTZ);
		    eventMZcutVariablesVector_.push_back(mZcut_dphiZjetZlep);
		    eventMZcutVariablesVector_.push_back(mZcut_detaZjetZlep);
		    eventMZcutVariablesVector_.push_back(mZcut_massTjetZjet);
		    eventMZcutVariablesVector_.push_back(mZcut_massTjetZlep);
		    eventMZcutVariablesVector_.push_back(mZcut_massZjetZlep);
		    // 3-particle vars
		    // ---------------
		    double mZcut_massTZZ = (Tjet->p4()+Zjet->p4()+Zlep->p4()).mass();
		    double mZcut_etaTZZ = fabs((Tjet->p4()+Zjet->p4()+Zlep->p4()).eta());
		    double mZcut_ptTZZ = (Tjet->p4()+Zjet->p4()+Zlep->p4()).pt();
		    MZCUT_MassTZZ->Fill(mZcut_massTZZ);
		    MZCUT_EtaTZZ->Fill(mZcut_etaTZZ);
		    MZCUT_PtTZZ->Fill(mZcut_ptTZZ);

		    eventMZcutVariablesVector_.push_back(mZcut_massTZZ);
		    eventMZcutVariablesVector_.push_back(mZcut_etaTZZ);
		    eventMZcutVariablesVector_.push_back(mZcut_ptTZZ);
		    
		    // Store the event number
		    eventMZcutVariablesVector_.push_back(eventcounter_);

		    if( writeTMVA_ ) {
		      // Fill the tree for the TMVA
		      tmvaMZcutTreeWriterPtr_->fill(eventMZcutVariablesVector_);
		    }
		    /*
		    else {
		      // Fill the values in the array used by the reader
		      vector<Float_t>::const_iterator varIt = eventSignalVariablesVector_.begin();
		      int iVarBackground = 0;
		      for( ; varIt != eventMZcutVariablesVector_.end(); ++varIt, ++iVarBackground ) {
		      variables_[iVarBackground] = *varIt;
		      }
		      discriminantBackground.push_back( reader_->EvaluateMVA( "BDT method" ) );
		    }
		    */
		    eventMZcutVariablesVector_.clear();
		  }
		}
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


  if ( eventcounter_/5000 == float(eventcounter_)/5000. ) {
    std::cout << "************************************************************************" << std::endl;
    std::cout << "lljjEventCounter: "      << lljjEventCounter_      << " => " << double(lljjEventCounter_     )/double(eventcounter_        )*100. << "%" << std::endl;
    std::cout << "tmvaMatchingCounter: "   << tmvaMatchingCounter_   << " => " << double(tmvaMatchingCounter_  )/double(lljjEventCounter_    )*100. << "%" << std::endl;
    std::cout << "combinatorialCounter: "  << combinatorialCounter_  << " => " << double(tmvaMatchingCounter_  )/double(combinatorialCounter_)*100. << "%" << std::endl;
    std::cout << "etaCutMatchingCounter: " << etaCutMatchingCounter_ << " => " << double(etaCutMatchingCounter_)/double(tmvaMatchingCounter_ )*100. << "%" << std::endl;
    std::cout << "mZMatchingCounter: "     << mZMatchingCounter_     << " => " << double(mZMatchingCounter_    )/double(tmvaMatchingCounter_ )*100. << "%" << std::endl; 
    std::cout << "************************************************************************" << std::endl;
  }

  
}


// ------------ method called once each job just before starting event loop  ------------
void 
VBFHZZalternativekinematicsAnalyzer::beginJob(const edm::EventSetup&)
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

  lljjEventNumber_               = new TH1D("lljjEventNumber",              "",1,0.,1.);
  tmvaMatchingEventNumber_       = new TH1D("tmvaMatchingEventNumber",      "",1,0.,1.);
  combinatorialMatchEventNumber_ = new TH1D("combinatorialMatchEventNumber","",1,0.,1.);
  etaCutMatchingEventNumber_     = new TH1D("etaCutMatchingEventNumber",    "",1,0.,1.);
  mZMatchingEventNumber_         = new TH1D("mZMatchingEventNumber",        "",1,0.,1.);

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
  jetParton_deltaEtaVSdeltaR2_ = new TH2D("jetParton_deltaEtaVSdeltaR2","#Delta#eta VS #DeltaR^{2} between partons and matched jets", nbin_,0.,0.25,nbin_, -5., 5.);
  jetParton_deltaEtVSdeltaR2_profile_         = new TProfile("jetParton_deltaEtVSdeltaR2_profile",        "#DeltaE_{T} VS #DeltaR^{2} between partons and matched jets",          nbin_,0.,0.25,-50.,50.);
  jetParton_deltaEtaVSdeltaR2_profile_ = new TProfile("jetParton_deltaEtaVSdeltaR2_profile","#Delta#eta VS #DeltaR^{2} between partons and matched jets", nbin_,0.,0.25, -5., 5.);
  jetParton_deltaEta_ = new TH1D("jetParton_deltaEta","#Delta#eta between partons and matched jets", nbin_, -5., 5.);
  jetParton_deltaEt_  = new TH1D("jetParton_deltaEt", "#DeltaE_{T} between partons and matched jets",nbin_,-50.,50.);
  jetParton_deltaEtRes_  = new TH1D("jetParton_deltaEtRes", "E_{T} resolution between partons and matched jets",nbin_,-1.,1.);

  deltaR2ZjetsZpartons_ = new TH1D("deltaR2ZjetsZpartons","",nbin_,0.,0.5);
  deltaR2TjetsTpartons_ = new TH1D("deltaR2TjetsTpartons","",nbin_,0.,0.5);

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
  Mhh = new TH1D ("Mhh","Mass of had Z system", nbin_, 0., 150.);
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
  MassZjetZlep = new TH1D("MassZjetZlep","Mass of ZZ system", nbin_, 0., 1000.);
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
  F_Mhh = new TH1D ("F_Mhh","Mass of had Z system", nbin_, 0., 150.);
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
  F_MassZjetZlep = new TH1D("F_MassZjetZlep","Mass of ZZ system", nbin_, 0., 1000.);
  F_MassTZZ = new TH1D("F_MassTZZ","Mass of TZZ system", nbin_, 0., 2000.);
  F_EtaTZZ = new TH1D("F_EtaTZZ","Eta of TZZ system", nbin_, 0., 5.);
  F_PtTZZ = new TH1D("F_PtTZZ","Pt of TZZ system", nbin_, 0., 500.);

  ETACUT_Dphiqq = new TH1D ("ETACUT_Dphiqq","Delta Phi tag jets", nbin_, 0., TMath::Pi());
  ETACUT_Detaqq = new TH1D ("ETACUT_Detaqq","Delta Eta tag jets", nbin_, 0., 10.0);
  ETACUT_Ptminqq = new TH1D ("ETACUT_Ptminqq","Pt of lowest Pt tag jet", nbin_, 0., 200.);
  ETACUT_Ptmaxqq = new TH1D ("ETACUT_Ptmaxqq","Pt of highest Pt tag jet", nbin_, 0., 500.);
  ETACUT_Etaminqq = new TH1D ("ETACUT_Etaminqq","Min |eta| of tag jet", nbin_, 0., 5.);
  ETACUT_Etamaxqq = new TH1D ("ETACUT_Etamaxqq","Max |eta| of tag jet", nbin_, 0., 5.);
  ETACUT_Dphihh = new TH1D ("ETACUT_Dphihh","Delta Phi Z jets", nbin_, 0., TMath::Pi());
  ETACUT_Detahh = new TH1D ("ETACUT_Detahh","Delta Eta Z jets", nbin_, 0., 5.0);
  ETACUT_Ptminhh = new TH1D ("ETACUT_Ptminhh","Pt of lowest Pt Z jet", nbin_, 0., 200.);
  ETACUT_Ptmaxhh = new TH1D ("ETACUT_Ptmaxhh","Pt of highest Pt Z jet", nbin_, 0., 500.);
  ETACUT_Etaminhh = new TH1D ("ETACUT_Etaminhh","Min |eta| of Z jet", nbin_, 0., 5.);
  ETACUT_Etamaxhh = new TH1D ("ETACUT_Etamaxhh","Max |eta| of Z jet", nbin_, 0., 5.);
  ETACUT_Ptqq = new TH1D ("ETACUT_Ptqq","Pt of tag jet system", nbin_, 0., 500.);
  ETACUT_Mqq = new TH1D ("ETACUT_Mqq","Mass of tag jet system", nbin_, 0., 2000.);
  ETACUT_Etaqq = new TH1D ("ETACUT_Etaqq","|Eta} of tag jet system", nbin_, 0., 5.);
  ETACUT_Pthh = new TH1D ("ETACUT_Pthh","Pt of had Z system", nbin_, 0., 500.);
  ETACUT_Mhh = new TH1D ("ETACUT_Mhh","Mass of had Z system", nbin_, 0., 150.);
  ETACUT_Etahh = new TH1D ("ETACUT_Etahh","|Eta} of had Z system", nbin_, 0., 5.);  
  ETACUT_Ptll = new TH1D ("ETACUT_Ptll","Pt of lep Z system", nbin_, 0., 500.);
  ETACUT_Mll = new TH1D ("ETACUT_Mll","Mass of lep Z system", nbin_, 0., 2000.);
  ETACUT_Etall = new TH1D ("ETACUT_Etall","|Eta} of lep Z system", nbin_, 0., 5.);  
  ETACUT_DphiTjetZjet = new TH1D("ETACUT_DphiTjetZjet","Delta phi between T and had Z syst", nbin_, 0., TMath::Pi());
  ETACUT_DphiTjetZlep = new TH1D("ETACUT_DphiTjetZlep","Delta phi between T and lep Z syst", nbin_, 0., TMath::Pi());
  ETACUT_DphiminTZ = new TH1D("ETACUT_DphiminTZ","Min DP between T and one Z system", nbin_, 0., TMath::Pi());
  ETACUT_DetaTjetZjet = new TH1D("ETACUT_DetaTjetZjet","Delta eta between T and had Z syst", nbin_, 0., 5.);
  ETACUT_DetaTjetZlep = new TH1D("ETACUT_DetaTjetZlep","Delta eta between T and lep Z syst", nbin_, 0., 5.);
  ETACUT_DetaminTZ = new TH1D("ETACUT_DetaminTZ","Delta eta min between T and one Z system", nbin_, 0., 5.);
  ETACUT_DphiZjetZlep = new TH1D("ETACUT_DphiZjetZlep","Delta phi between had and lep Z syst", nbin_, 0., TMath::Pi());
  ETACUT_DetaZjetZlep = new TH1D("ETACUT_DetaZjetZlep","Delta eta between had and lep Z syst", nbin_, 0., 5.);
  ETACUT_MassTjetZjet = new TH1D("ETACUT_MassTjetZjet","Mass of T and had Z syst", nbin_, 0., 2000.);
  ETACUT_MassTjetZlep = new TH1D("ETACUT_MassTjetZlep","Mass of T and lep Z syst", nbin_, 0., 2000.);
  ETACUT_MassZjetZlep = new TH1D("ETACUT_MassZjetZlep","Mass of ZZ system", nbin_, 0., 1000.);
  ETACUT_MassTZZ = new TH1D("ETACUT_MassTZZ","Mass of TZZ system", nbin_, 0., 2000.);
  ETACUT_EtaTZZ = new TH1D("ETACUT_EtaTZZ","Eta of TZZ system", nbin_, 0., 5.);
  ETACUT_PtTZZ = new TH1D("ETACUT_PtTZZ","Pt of TZZ system", nbin_, 0., 500.);

  MZCUT_Dphiqq = new TH1D ("MZCUT_Dphiqq","Delta Phi tag jets", nbin_, 0., TMath::Pi());
  MZCUT_Detaqq = new TH1D ("MZCUT_Detaqq","Delta Eta tag jets", nbin_, 0., 10.0);
  MZCUT_Ptminqq = new TH1D ("MZCUT_Ptminqq","Pt of lowest Pt tag jet", nbin_, 0., 200.);
  MZCUT_Ptmaxqq = new TH1D ("MZCUT_Ptmaxqq","Pt of highest Pt tag jet", nbin_, 0., 500.);
  MZCUT_Etaminqq = new TH1D ("MZCUT_Etaminqq","Min |eta| of tag jet", nbin_, 0., 5.);
  MZCUT_Etamaxqq = new TH1D ("MZCUT_Etamaxqq","Max |eta| of tag jet", nbin_, 0., 5.);
  MZCUT_Dphihh = new TH1D ("MZCUT_Dphihh","Delta Phi Z jets", nbin_, 0., TMath::Pi());
  MZCUT_Detahh = new TH1D ("MZCUT_Detahh","Delta Eta Z jets", nbin_, 0., 5.0);
  MZCUT_Ptminhh = new TH1D ("MZCUT_Ptminhh","Pt of lowest Pt Z jet", nbin_, 0., 200.);
  MZCUT_Ptmaxhh = new TH1D ("MZCUT_Ptmaxhh","Pt of highest Pt Z jet", nbin_, 0., 500.);
  MZCUT_Etaminhh = new TH1D ("MZCUT_Etaminhh","Min |eta| of Z jet", nbin_, 0., 5.);
  MZCUT_Etamaxhh = new TH1D ("MZCUT_Etamaxhh","Max |eta| of Z jet", nbin_, 0., 5.);
  MZCUT_Ptqq = new TH1D ("MZCUT_Ptqq","Pt of tag jet system", nbin_, 0., 500.);
  MZCUT_Mqq = new TH1D ("MZCUT_Mqq","Mass of tag jet system", nbin_, 0., 2000.);
  MZCUT_Etaqq = new TH1D ("MZCUT_Etaqq","|Eta} of tag jet system", nbin_, 0., 5.);
  MZCUT_Pthh = new TH1D ("MZCUT_Pthh","Pt of had Z system", nbin_, 0., 500.);
  MZCUT_Mhh = new TH1D ("MZCUT_Mhh","Mass of had Z system", nbin_, 0., 150.);
  MZCUT_Etahh = new TH1D ("MZCUT_Etahh","|Eta} of had Z system", nbin_, 0., 5.);  
  MZCUT_Ptll = new TH1D ("MZCUT_Ptll","Pt of lep Z system", nbin_, 0., 500.);
  MZCUT_Mll = new TH1D ("MZCUT_Mll","Mass of lep Z system", nbin_, 0., 2000.);
  MZCUT_Etall = new TH1D ("MZCUT_Etall","|Eta} of lep Z system", nbin_, 0., 5.);  
  MZCUT_DphiTjetZjet = new TH1D("MZCUT_DphiTjetZjet","Delta phi between T and had Z syst", nbin_, 0., TMath::Pi());
  MZCUT_DphiTjetZlep = new TH1D("MZCUT_DphiTjetZlep","Delta phi between T and lep Z syst", nbin_, 0., TMath::Pi());
  MZCUT_DphiminTZ = new TH1D("MZCUT_DphiminTZ","Min DP between T and one Z system", nbin_, 0., TMath::Pi());
  MZCUT_DetaTjetZjet = new TH1D("MZCUT_DetaTjetZjet","Delta eta between T and had Z syst", nbin_, 0., 5.);
  MZCUT_DetaTjetZlep = new TH1D("MZCUT_DetaTjetZlep","Delta eta between T and lep Z syst", nbin_, 0., 5.);
  MZCUT_DetaminTZ = new TH1D("MZCUT_DetaminTZ","Delta eta min between T and one Z system", nbin_, 0., 5.);
  MZCUT_DphiZjetZlep = new TH1D("MZCUT_DphiZjetZlep","Delta phi between had and lep Z syst", nbin_, 0., TMath::Pi());
  MZCUT_DetaZjetZlep = new TH1D("MZCUT_DetaZjetZlep","Delta eta between had and lep Z syst", nbin_, 0., 5.);
  MZCUT_MassTjetZjet = new TH1D("MZCUT_MassTjetZjet","Mass of T and had Z syst", nbin_, 0., 2000.);
  MZCUT_MassTjetZlep = new TH1D("MZCUT_MassTjetZlep","Mass of T and lep Z syst", nbin_, 0., 2000.);
  MZCUT_MassZjetZlep = new TH1D("MZCUT_MassZjetZlep","Mass of ZZ system", nbin_, 0., 1000.);
  MZCUT_MassTZZ = new TH1D("MZCUT_MassTZZ","Mass of TZZ system", nbin_, 0., 2000.);
  MZCUT_EtaTZZ = new TH1D("MZCUT_EtaTZZ","Eta of TZZ system", nbin_, 0., 5.);
  MZCUT_PtTZZ = new TH1D("MZCUT_PtTZZ","Pt of TZZ system", nbin_, 0., 500.);

  // =========================================================

}

// ------------ method called once each job just after ending the event loop  ------------
void 
VBFHZZalternativekinematicsAnalyzer::endJob() {

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


  lljjEventNumber_->SetBinContent(1,double(lljjEventCounter_));
  tmvaMatchingEventNumber_->SetBinContent(1,double(tmvaMatchingCounter_));
  combinatorialMatchEventNumber_->SetBinContent(1,double(combinatorialCounter_));
  etaCutMatchingEventNumber_->SetBinContent(1,double(etaCutMatchingCounter_));
  mZMatchingEventNumber_->SetBinContent(1,double(mZMatchingCounter_));

  lljjEventNumber_->Write();
  tmvaMatchingEventNumber_->Write();
  combinatorialMatchEventNumber_->Write();
  etaCutMatchingEventNumber_->Write();
  mZMatchingEventNumber_->Write();

  std::cout << "lljjEventCounter: "      << lljjEventCounter_      << " => " << double(lljjEventCounter_     )/double(eventcounter_        )*100. << "%" << std::endl;
  std::cout << "tmvaMatchingCounter: "   << tmvaMatchingCounter_   << " => " << double(tmvaMatchingCounter_  )/double(lljjEventCounter_    )*100. << "%" << std::endl;
  std::cout << "combinatorialCounter: "  << combinatorialCounter_  << " => " << double(tmvaMatchingCounter_  )/double(combinatorialCounter_)*100. << "%" << std::endl;
  std::cout << "etaCutMatchingCounter: " << etaCutMatchingCounter_ << " => " << double(etaCutMatchingCounter_)/double(tmvaMatchingCounter_ )*100. << "%" << std::endl;
  std::cout << "mZMatchingCounter: "     << mZMatchingCounter_     << " => " << double(mZMatchingCounter_    )/double(tmvaMatchingCounter_ )*100. << "%" << std::endl; 


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
  jetParton_deltaEtaVSdeltaR2_->Write();
  jetParton_deltaEtVSdeltaR2_profile_->Write();
  jetParton_deltaEtaVSdeltaR2_profile_->Write();
  jetParton_deltaEta_->Write();
  jetParton_deltaEt_->Write();

  deltaR2ZjetsZpartons_->Write();
  deltaR2TjetsTpartons_->Write();

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
 
  ETACUT_Dphiqq->Write();
  ETACUT_Detaqq->Write();
  ETACUT_Ptminqq->Write();
  ETACUT_Ptmaxqq->Write();
  ETACUT_Etaminqq->Write();
  ETACUT_Etamaxqq->Write();
  ETACUT_Dphihh->Write();
  ETACUT_Detahh->Write();
  ETACUT_Ptminhh->Write();
  ETACUT_Ptmaxhh->Write();
  ETACUT_Etaminhh->Write();
  ETACUT_Etamaxhh->Write();
  ETACUT_Ptqq->Write();
  ETACUT_Mqq->Write();
  ETACUT_Etaqq->Write();
  ETACUT_Pthh->Write();
  ETACUT_Mhh->Write();
  ETACUT_Etahh->Write();
  ETACUT_Ptll->Write();
  ETACUT_Mll->Write();
  ETACUT_Etall->Write();
  ETACUT_DphiTjetZjet->Write();
  ETACUT_DphiTjetZlep->Write();
  ETACUT_DphiminTZ->Write();
  ETACUT_DetaTjetZjet->Write();
  ETACUT_DetaTjetZlep->Write();
  ETACUT_DetaminTZ->Write();
  ETACUT_MassTjetZjet->Write();
  ETACUT_MassTjetZlep->Write();
  ETACUT_MassZjetZlep->Write();
  ETACUT_MassTZZ->Write();
  ETACUT_EtaTZZ->Write();
  ETACUT_PtTZZ->Write();
 
  MZCUT_Dphiqq->Write();
  MZCUT_Detaqq->Write();
  MZCUT_Ptminqq->Write();
  MZCUT_Ptmaxqq->Write();
  MZCUT_Etaminqq->Write();
  MZCUT_Etamaxqq->Write();
  MZCUT_Dphihh->Write();
  MZCUT_Detahh->Write();
  MZCUT_Ptminhh->Write();
  MZCUT_Ptmaxhh->Write();
  MZCUT_Etaminhh->Write();
  MZCUT_Etamaxhh->Write();
  MZCUT_Ptqq->Write();
  MZCUT_Mqq->Write();
  MZCUT_Etaqq->Write();
  MZCUT_Pthh->Write();
  MZCUT_Mhh->Write();
  MZCUT_Etahh->Write();
  MZCUT_Ptll->Write();
  MZCUT_Mll->Write();
  MZCUT_Etall->Write();
  MZCUT_DphiTjetZjet->Write();
  MZCUT_DphiTjetZlep->Write();
  MZCUT_DphiminTZ->Write();
  MZCUT_DetaTjetZjet->Write();
  MZCUT_DetaTjetZlep->Write();
  MZCUT_DetaminTZ->Write();
  MZCUT_MassTjetZjet->Write();
  MZCUT_MassTjetZlep->Write();
  MZCUT_MassZjetZlep->Write();
  MZCUT_MassTZZ->Write();
  MZCUT_EtaTZZ->Write();
  MZCUT_PtTZZ->Write();
 
  // ====================================================

}
//define this as a plug-in
DEFINE_FWK_MODULE(VBFHZZalternativekinematicsAnalyzer);

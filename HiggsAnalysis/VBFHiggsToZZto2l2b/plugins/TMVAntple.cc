//
// Original Author:  Mia Tosi
//         Created:  Fri Feb 22 17:56:22 CET 2008
// $Id: VBFHZZalternativekinematicsAnalyzer.cc,v 1.1 2008/12/05 17:22:51 tosi Exp $
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

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/TMVAntple.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

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

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/PythiaParticleIndex.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/ParticlesMass.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/ParticlesCharge.h"

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
       NMINJETS  = 4
};


TMVAntple::TMVAntple(const edm::ParameterSet& iConfig) :
  signal_               ( iConfig.getParameter<int>           ( "signal"               ) ), // 0:Signal,  1:Background
  whichSim_             ( iConfig.getParameter<int>           ( "whichSim"             ) ), // 0:FastSim, 1:FullSim
  electronLabel_        ( iConfig.getParameter<edm::InputTag> ( "electronLabel"        ) ),
  muonLabel_            ( iConfig.getParameter<edm::InputTag> ( "muonLabel"            ) ),
  metLabel_             ( iConfig.getParameter<edm::InputTag> ( "metLabel"             ) ),
  jetLabel_             ( iConfig.getParameter<edm::InputTag> ( "jetLabel"             ) ),
  corJetsWithBTagLabel_ ( iConfig.getParameter<std::string>   ( "corJetsWithBTagLabel" ) ),
  mcParticleLabel_      ( iConfig.getParameter<edm::InputTag> ( "mcParticleLabel"      ) ),
  genJetLabel_          ( iConfig.getParameter<edm::InputTag> ( "genJetLabel"          ) ),
  genMetLabel_          ( iConfig.getParameter<edm::InputTag> ( "genMetLabel"          ) ),

  leptonPtCut_         ( iConfig.getUntrackedParameter<double> ("leptonPtCut"        ) ),
  jetEtCut_            ( iConfig.getUntrackedParameter<double> ("jetEtCut"           ) ),
  jetPartonDeltaR2Cut_ ( iConfig.getUntrackedParameter<double> ("jetPartonDeltaR2Cut") ),
  jetLeptonDeltaRCut_  ( iConfig.getUntrackedParameter<double> ("jetLeptonDeltaRCut" ) ),
  jetEMfracCut_        ( iConfig.getUntrackedParameter<double> ("jetEMfracCut"       ) ),
  mHres_               ( iConfig.getUntrackedParameter<double> ("mHres"              ) ),
  mH_                  ( iConfig.getUntrackedParameter<double> ("mH"                 ) ),
  tmvaSignalSuffix_ ( iConfig.getUntrackedParameter<std::string> ("tmvaSignalSuffix") ),
  tmvaCombinSuffix_ ( iConfig.getUntrackedParameter<std::string> ("tmvaCombinSuffix") ),
  tmvaEtaCutSuffix_ ( iConfig.getUntrackedParameter<std::string> ("tmvaEtaCutSuffix") ),
  tmvaMZcutSuffix_  ( iConfig.getUntrackedParameter<std::string> ("tmvaMZcutSuffix" ) )

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
  eventcounter_             = 0;

  lljjEventCounter_      = 0;
  tmvaMatchingCounter_   = 0;
  combinatorialCounter_  = 0;
  etaCutMatchingCounter_ = 0;
  mZMatchingCounter_     = 0;

  nbin_         = 100;

  null_XYZTLorentzVector_ = math::XYZTLorentzVector(0.,0.,0.,0.);

  edm::Service<TFileService> fs ;

  gROOT->Time();

}


TMVAntple::~TMVAntple()
{
  std::cout << "endJob" << std::endl;

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
TMVAntple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
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


  unsigned muonCollectionSize       = muonHandle->size();
  unsigned electronCollectionSize   = electronHandle->size();
  unsigned jetCollectionSize        = jetHandle->size();
  unsigned mcParticleCollectionSize = mcParticleHandle->size();

  if ( ( muonCollectionSize >= NZLEPTONS || electronCollectionSize >= NZLEPTONS ) &&
       jetCollectionSize >= NMINJETS ) {

    /////////////////////////////////// HEPG analysis ////////////////////////////
  // Z counters 
    unsigned int Zparticles = 0;
    unsigned int hadronicZparticles = 0;
    unsigned int leptonicZparticles = 0;
    unsigned int ZintoEcounter   = 0; 
    unsigned int ZintoMUcounter  = 0;
    unsigned int ZintoTAUcounter = 0;
    unsigned int ZintoBQUARKcounter = 0;
    
    std::vector< const Candidate * > tagSystem;
    std::vector< const Candidate * > H;
    std::vector< const Candidate * > hadronicZ;
    std::vector< const Candidate * > leptonicZ;
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
    
//    /////////////////////////////////// QUICK LOOP over LEPTONS ////////////////////////////
//    // Global muons
//    // ------------
//    std::vector< GlobalMuon > goodMuonVec;
//    std::vector< std::pair<GlobalMuon,GlobalMuon> > goodMuonsPairVec;
//    GlobalMuonCollection::const_iterator globalMuon_itr = globalMuons->begin();
//    for ( ; globalMuon_itr != globalMuons->end(); ++globalMuon_itr ) {
//      double tmpMuonPt  = globalMuon_itr->pt();
//      double tmpMuonEta = globalMuon_itr->eta();
//      double tmpMuonPhi = globalMuon_itr->phi();
//      if ( tmpMuonPt >= leptonPtCut_ ) {
//	goodMuonVec.push_back(*globalMuon_itr);
//	if ( globalMuon_itr != (globalMuons->end()-1) ) {
//	  GlobalMuonCollection::const_iterator globalMuon_itr2 = globalMuon_itr+1;
//	  for ( ; globalMuon_itr2 != globalMuons->end(); ++globalMuon_itr2 ){
//	    if ( globalMuon_itr->charge()*globalMuon_itr2->charge() < 0 ) {
//	      std::pair< GlobalMuon, GlobalMuon > tmpMuonsPair(*globalMuon_itr,*globalMuon_itr2);
//	      goodMuonsPairVec.push_back(tmpMuonsPair);
//	    }
//	  }
//	}
//      }
//    }
//    //    if ( goodMuonsPairVec.size() >= 1 ) {
//    //      std::cout << "goodMuonsPairVec: " << goodMuonsPairVec.size() << std::endl;
//    //      std::cout << "goodMuonVec: " << goodMuonVec.size() << std::endl;
//    //    }
//
//    // SimpleElectrons
//    // ---------------
//    std::vector< SimpleElectron > goodElectronVec;
//    std::vector< std::pair<SimpleElectron,SimpleElectron> > goodElectronsPairVec;
//    SimpleElectronCollection::const_iterator simpleElectron_itr = simpleElectrons->begin();
//    for ( ; simpleElectron_itr != simpleElectrons->end(); ++simpleElectron_itr ) {
//      double tmpElectronEt  = simpleElectron_itr->et();
//      double tmpElectronEta = simpleElectron_itr->eta();
//      double tmpElectronPhi = simpleElectron_itr->phi();
//      double tmpElectronPt  = simpleElectron_itr->pt();
//      if ( tmpElectronPt >= leptonPtCut_ ) {
//	goodElectronVec.push_back(*simpleElectron_itr);
//	if ( simpleElectron_itr != (simpleElectrons->end()-1) ) {
//	  SimpleElectronCollection::const_iterator simpleElectron_itr2 = simpleElectron_itr+1;
//	  for ( ; simpleElectron_itr2 != simpleElectrons->end(); ++simpleElectron_itr2 ){
//	    if ( simpleElectron_itr->charge()*simpleElectron_itr2->charge() < 0 ) {
//	      std::pair< SimpleElectron, SimpleElectron> tmpElectronsPair(*simpleElectron_itr,*simpleElectron_itr2);
//	      goodElectronsPairVec.push_back(tmpElectronsPair);
//	    }
//	  }
//	}
//      }
//    }
//    //    if ( goodElectronsPairVec.size() >= 1 ) {
//    //      std::cout << "goodElectronsPairVec: " << goodElectronsPairVec.size() << std::endl;
//    //      std::cout << "goodElectronVec: " << goodElectronVec.size() << std::endl;
//    //    }
//  
//    std::vector<math::XYZTLorentzVector> ZlepVec;
//    if ( goodMuonVec.size() >= NZLEPTONS || goodElectronVec.size() >= NZLEPTONS ) {
//
//      //    /////////////////////////////////// LEPTON MATCHING ////////////////////////////
//      if ( leptonicZcandidates.size() != 0 ) {
//	int leptonId = int(fabs(leptonicZcandidates[0]->pdgId()));
//	if ( leptonId == pythiamu_ || leptonId == pythiae_ ) {
//
//	  std::vector< const Candidate* >::const_iterator leptonicZcandidates_itr = leptonicZcandidates.begin();
//	  for ( ; leptonicZcandidates_itr != leptonicZcandidates.end(); ++leptonicZcandidates_itr ) 
//	    ZlepVec.push_back(null_XYZTLorentzVector_);
//
//	  leptonicZcandidates_itr = leptonicZcandidates.begin();
//	  unsigned int leptonicZcandidateIndex = 0;
//
//	  for ( ; leptonicZcandidates_itr != leptonicZcandidates.end(); ++leptonicZcandidates_itr, leptonicZcandidateIndex++ ) {
//	    double tmp_leptonicZcandidateEta = (*leptonicZcandidates_itr)->eta();
//	    double tmp_leptonicZcandidatePhi = (*leptonicZcandidates_itr)->phi();
//	    double tmp_leptonicZcandidateCharge = (*leptonicZcandidates_itr)->charge();
//	    double leptonDeltaR2 = 99.;
//	    if ( leptonId == pythiamu_) {
//	      std::vector< GlobalMuon >::const_iterator goodMuonVec_itr = goodMuonVec.begin();
//	      for ( ; goodMuonVec_itr != goodMuonVec.end(); ++goodMuonVec_itr ) {
//		double tmp_leptonEta    = goodMuonVec_itr->eta();
//		double tmp_leptonPhi    = goodMuonVec_itr->phi();
//		double tmp_leptonCharge = goodMuonVec_itr->charge();
//		double tmp_deltaEta = tmp_leptonicZcandidateEta - tmp_leptonEta;
//		double tmp_deltaPhi = DeltaPhi(tmp_leptonicZcandidatePhi,tmp_leptonPhi);
//		double tmp_deltaR2 = pow(tmp_deltaEta,2)+pow(tmp_deltaPhi,2);
//		if (  tmp_leptonicZcandidateCharge*tmp_leptonCharge > 0  ) {
//		  if ( tmp_deltaR2 <= leptonDeltaR2 ) {
//		    ZlepVec[leptonicZcandidateIndex] = goodMuonVec_itr->p4();
//		    leptonDeltaR2 = tmp_deltaR2;
//		  }
//		}
//	      }
//	    } else if (leptonId == pythiae_ ) {
//	      std::vector< SimpleElectron >::const_iterator goodElectronVec_itr = goodElectronVec.begin();
//	      for ( ; goodElectronVec_itr != goodElectronVec.end(); ++goodElectronVec_itr ) {
//		double tmp_leptonEta    = goodElectronVec_itr->eta();
//		double tmp_leptonPhi    = goodElectronVec_itr->phi();
//		double tmp_leptonCharge = goodElectronVec_itr->charge();
//		double tmp_deltaEta = tmp_leptonicZcandidateEta - tmp_leptonEta;
//		double tmp_deltaPhi = DeltaPhi(tmp_leptonicZcandidatePhi,tmp_leptonPhi);
//		double tmp_deltaR2 = pow(tmp_deltaEta,2)+pow(tmp_deltaPhi,2);
//		if (  tmp_leptonicZcandidateCharge*tmp_leptonCharge > 0  ) {
//		  if ( tmp_deltaR2 <= leptonDeltaR2 ) {
//		    ZlepVec[leptonicZcandidateIndex] = goodElectronVec_itr->p4();
//		    leptonDeltaR2 = tmp_deltaR2;
//		  }
//		}
//	      }
//	    }
//	  }
//	  //	  std::cout << "ZlepVec: " << ZlepVec.size() << std::endl;
//	  //	  if ( ZlepVec[0].pt() == 0. || ZlepVec[1].pt() == 0 )
//	  //	    std::cout << "WARNING: ZlepVec[0]->pt: " << ZlepVec[0].pt() << " <--> ZlepVec[1]->pt: " << ZlepVec[1].pt() << std::endl;
//	}
//      }
//    }
//    if ( ZlepVec.size() >= NZLEPTONS && ZlepVec[0].pt() != 0. && ZlepVec[1].pt() != 0. ) {
//
//      if ( partonsCandidates.size() >= NMINJETS ) {
//
//	//    /////////////////////////////////// JET-PARTON MATCHING ////////////////////////////
//	// Calorimeter jets
//	// ----------------
//	std::vector< OfflineJet > aboveEtCutJetVec;
//	OfflineJetCollection::const_iterator offlinejet = offlineJets->begin();
//	for ( ; offlinejet != offlineJets->end(); ++offlinejet ) {
//	  double tmpOfflinejetEt  = offlinejet->et();
//	  double tmpOfflinejetEta = offlinejet->eta();
//	  if ( tmpOfflinejetEt >= jetEtCut_ ) {
//	    aboveEtCutJetVec.push_back(*offlinejet);
//	  }
//	}
//	if ( aboveEtCutJetVec.size() >= NZJETS ) { 
//	  if ( aboveEtCutJetVec.size() >= njets_ ) {
//	
//	    std::vector< OfflineJet > ZjetVec;
//	    std::vector< double >     ZjetDeltaR2Vec;
//	    std::vector<int>          ZjetInd;
//	    unsigned int hadronicZcandidatesNumber = hadronicZcandidates.size();
//	    for ( unsigned int index = 0; index < hadronicZcandidatesNumber; index++) 
//	      ZjetDeltaR2Vec.push_back(99.);
//	    
//	    std::vector<OfflineJet> TjetVec;
//	    std::vector< double >   TjetDeltaR2Vec;
//	    std::vector<int>        TjetInd;
//	    unsigned int TagCandidatesNumber = TagCandidates.size();
//	    for ( unsigned int index = 0; index < TagCandidatesNumber; index++) 
//	      TjetDeltaR2Vec.push_back(99.);
//	    
//
//	    std::vector<std::vector<double> > closestJetDeltaR2Vec;
//	    std::vector<std::vector<int> >    closestJetIndexVec;
//	    unsigned int partonsCandidatesNumber = partonsCandidates.size();
//	    for ( unsigned index = 0; index < partonsCandidatesNumber; index++ ) {
//	      std::vector<double> null_closestJetDeltaR2Vec;
//	      std::vector<int>    null_closestJetIndexVec;
//	      for ( unsigned index = 0; index < partonsCandidatesNumber; index++ ) {
//		null_closestJetDeltaR2Vec.push_back(99.);
//		null_closestJetIndexVec.push_back(-1);
//	      }
//	      closestJetDeltaR2Vec.push_back(null_closestJetDeltaR2Vec);
//	      closestJetIndexVec.push_back(null_closestJetIndexVec);
//	    }
//
//	    std::vector< const Candidate* >::const_iterator partonsCandidates_itr = partonsCandidates.begin();
//	    unsigned int partonIndex = 0;
//	    for ( ; partonsCandidates_itr != partonsCandidates.end(); ++partonsCandidates_itr,
//		    partonIndex++ ) {
//	      double partonPt  = partonsCandidates[partonIndex]->pt();
//	      double partonEta = partonsCandidates[partonIndex]->eta();
//	      double partonPhi = partonsCandidates[partonIndex]->phi();
//	      //	      std::cout << "parton[" << partonIndex << "]: pt: " << partonPt
//	      //			<< " eta: " << partonEta
//	      //			<< " phi: " << partonPhi << std::endl;
//	      std::vector< OfflineJet >::const_iterator aboveEtCutJetVec_itr = aboveEtCutJetVec.begin();
//	      unsigned int jetIndex = 0;
//	      for ( ; aboveEtCutJetVec_itr != aboveEtCutJetVec.end(); ++aboveEtCutJetVec_itr, 
//		      jetIndex++ ) {
//		double jetEt  = aboveEtCutJetVec_itr->et();
//		double jetEta = aboveEtCutJetVec_itr->eta();
//		double jetPhi = aboveEtCutJetVec_itr->phi();
//		//		std::cout << "jet[" << jetIndex << "]: et: " << jetEt
//		//			  << " eta: " << jetEta
//		//			  << " phi: " << jetPhi << std::endl;
//		
//		double deltaEta = jetEta - partonEta;
//		double deltaPhi = DeltaPhi(jetPhi,partonPhi);
//		double deltaR2 = pow(deltaEta,2)+pow(deltaPhi,2);
//		//		std::cout << "deltaR2: " << deltaR2 << std::endl;
//
//		//		if ( deltaR2 <= jetPartonDeltaR2Cut_ ) {
//
//		if ( deltaR2 <= closestJetDeltaR2Vec[partonIndex][3] ) {
//		  if ( deltaR2 <= closestJetDeltaR2Vec[partonIndex][2] ) {
//		    if ( deltaR2 <= closestJetDeltaR2Vec[partonIndex][1] ) {
//		      if ( deltaR2 <= closestJetDeltaR2Vec[partonIndex][0] ) {
//			closestJetDeltaR2Vec[partonIndex][3] = closestJetDeltaR2Vec[partonIndex][2];
//			closestJetDeltaR2Vec[partonIndex][2] = closestJetDeltaR2Vec[partonIndex][1];
//			closestJetDeltaR2Vec[partonIndex][1] = closestJetDeltaR2Vec[partonIndex][0];
//			closestJetDeltaR2Vec[partonIndex][0] = deltaR2;
//			closestJetIndexVec[partonIndex][3] = closestJetIndexVec[partonIndex][2];
//			closestJetIndexVec[partonIndex][2] = closestJetIndexVec[partonIndex][1];
//			closestJetIndexVec[partonIndex][1] = closestJetIndexVec[partonIndex][0];
//			closestJetIndexVec[partonIndex][0] = jetIndex;
//			  
//		      } else {
//			closestJetDeltaR2Vec[partonIndex][3] = closestJetDeltaR2Vec[partonIndex][2];
//			closestJetDeltaR2Vec[partonIndex][2] = closestJetDeltaR2Vec[partonIndex][1];
//			closestJetDeltaR2Vec[partonIndex][1] = deltaR2;
//			closestJetIndexVec[partonIndex][3] = closestJetIndexVec[partonIndex][2];
//			closestJetIndexVec[partonIndex][2] = closestJetIndexVec[partonIndex][1];
//			closestJetIndexVec[partonIndex][1] = jetIndex;
//		      }
//		    } else {
//		      closestJetDeltaR2Vec[partonIndex][3] = closestJetDeltaR2Vec[partonIndex][2];
//		      closestJetDeltaR2Vec[partonIndex][2] = deltaR2;
//		      closestJetIndexVec[partonIndex][3] = closestJetIndexVec[partonIndex][2];
//		      closestJetIndexVec[partonIndex][2] = jetIndex;
//		    }
//		  } else {
//		    closestJetDeltaR2Vec[partonIndex][3] = deltaR2;
//		    closestJetIndexVec[partonIndex][3] = jetIndex;
//		  }
//		}
//		//		}
//	      } // end loop over jet above et cut
//	      for (int index = 0; index < 4; index++ ) {
//		int jetMatchedIndex = closestJetIndexVec[partonIndex][index];
//		double etRes = (aboveEtCutJetVec[jetMatchedIndex].et()-partonPt)/partonPt;
//		double jetMass = aboveEtCutJetVec[jetMatchedIndex].p4().M();
//		//		std::cout << "closestJetDeltaR2Vec[" << partonIndex << "][" << index << "]: " << jetMatchedIndex 
//		//			  << " --> deltaR2: " << closestJetDeltaR2Vec[partonIndex][index] << std::endl;
//		
//	      }
//	    } // end loop over partons
//	  
//
//
//	    std::vector<bool> partonAss;
//	    for ( unsigned int partonIndex = 0; partonIndex < partonsCandidatesNumber; partonIndex++ ) 
//	      partonAss.push_back(kFALSE);
//
//	    for ( unsigned int assIndex = 0; assIndex < partonsCandidatesNumber; assIndex++ ) {
//	      //	      std::cout << "assIndex: " << assIndex << std::endl;
//	      double minDeltaR2 = jetPartonDeltaR2Cut_;
//	      int minPartonIndex = -1;
//	      // find the best DeltaR2 matching to find the best parton association in the collection
//	      for ( unsigned int partonIndex = 0; partonIndex < partonsCandidatesNumber; partonIndex++ ) {
//		//		std::cout << "partonIndex: " << partonIndex << " partonAss: " << partonAss[partonIndex] << std::endl;
//		if ( !partonAss[partonIndex] && closestJetDeltaR2Vec[partonIndex][0] <= minDeltaR2 ) {
//		  minDeltaR2 = closestJetDeltaR2Vec[partonIndex][0];
//		  minPartonIndex = partonIndex;
//		}
//	      } // end loop over partons [not matched yet]
//
//	      // parton association
//	      if ( minPartonIndex >= 0 ) {
//		partonAss[minPartonIndex] = kTRUE;
//		// save the matched jet into the proper vector
//		unsigned int jetIndex = closestJetIndexVec[minPartonIndex][0];
//		//	      std::cout << "minPartonIndex: " << minPartonIndex 
//		//			<< " to jetIndex: " << jetIndex
//		//			<< " => nAss: " << nAss << std::endl;
//
//		if ( minPartonIndex == 0 || minPartonIndex == 1 ) {
//		  ZjetVec.push_back(aboveEtCutJetVec[jetIndex]);
//		  ZjetInd.push_back(jetIndex);
//		  double deltaR2 = closestJetDeltaR2Vec[minPartonIndex][0];
//		  double deltaEt  = partonsCandidates[minPartonIndex]->et()-aboveEtCutJetVec[jetIndex].et();
//		  double deltaEta = partonsCandidates[minPartonIndex]->eta()-aboveEtCutJetVec[jetIndex].eta();
//
//		  //		  std::cout << "Z jetIndex: " <<  jetIndex << std::endl;
//		}
//		else if ( minPartonIndex == 2 || minPartonIndex == 3 ) {
//		  TjetVec.push_back(aboveEtCutJetVec[jetIndex]);
//		  TjetInd.push_back(jetIndex);
//		  //		  std::cout << "T jetIndex: " <<  jetIndex << std::endl;
//		}
//		// in case of "non-biunivocity" pop-up jet associations belong to worst DeltaR2 parton
//		for ( unsigned int partonIndex = 0; partonIndex < partonsCandidatesNumber; partonIndex++ ) {
//		  if ( partonAss[partonIndex] ) continue;
//		  //		std::cout << "partonIndex: " << partonIndex 
//		  //			  << " partonAss: " << partonAss[partonIndex] << std::endl;
//		  //		std::cout << "closestJetIndex: " << closestJetIndexVec[partonIndex][0] << std::endl;
//		  for ( unsigned int iAssIndex = 0; iAssIndex < partonsCandidatesNumber; iAssIndex++ ) {
//		    if ( closestJetIndexVec[partonIndex][iAssIndex] != closestJetIndexVec[minPartonIndex][0] ) continue;
//		    //		  std::cout << "iAssIndex: " << iAssIndex << std::endl;
//		    for ( unsigned int jAssIndex = iAssIndex+1; jAssIndex < partonsCandidatesNumber; jAssIndex++ ) {
//		      //		    std::cout << "jAssIndex: " << jAssIndex << std::endl;
//		      closestJetDeltaR2Vec[partonIndex][jAssIndex-1] = closestJetDeltaR2Vec[partonIndex][jAssIndex];
//		      closestJetIndexVec[partonIndex][jAssIndex-1] = closestJetIndexVec[partonIndex][jAssIndex];
//		    }
//		  }
//		  //		std::cout << "closestJetIndex: " << closestJetIndexVec[partonIndex][0] << std::endl;
//		} // end loop over partons [not matched yet and w/ the same jet]
//	      }
//	    } // end loop over association index
//
//	    //	    std::cout << "************************************" << std::endl;
//	    //	    std::cout << "ZjetVec: " << ZjetVec.size() << std::endl;
//	    //	    std::cout << "TjetVec: " << TjetVec.size() << std::endl;
//	    //	    std::cout << "ZlepVec: " << ZlepVec.size() << std::endl;
//
//
//	    // ====================================================================================================
//	    if (ZjetVec.size()==2 && TjetVec.size() == 2 ){
//	      //	    if (partonAss[0] && partonAss[1] && partonAss[2] && partonAss[3] ) {
//	      
//
//	      Particle * Zjet = new Particle(0,ZjetVec[0].p4()+ZjetVec[1].p4(), ZjetVec[0].vertex(),     pythiaZ_,0,true);
//	      Particle * Zlep = new Particle(0,ZlepVec[0]+ZlepVec[1],           math::XYZPoint(0.,0.,0.),pythiaZ_,0,true);
//	      Particle * Tjet = new Particle(0,TjetVec[0].p4()+TjetVec[1].p4(), TjetVec[0].vertex(),     0,       0,true);
//
//
//	      double hadronicZrecMass = Zjet->mass();
//	      double hadronicZMass    = hadronicZ->mass();
//	      double hadronicZrecMassResolution = -99.;
//	      if ( hadronicZMass != 0. ) hadronicZrecMassResolution = (hadronicZrecMass-hadronicZMass)/hadronicZMass;
//	      double hadronicZrecPt = Zjet->pt();
//	      double hadronicZPt    = hadronicZ->pt();
//	      double hadronicZrecPtResolution = -99.;
//	      if ( hadronicZPt != 0. ) hadronicZrecPtResolution = (hadronicZrecPt-hadronicZPt)/hadronicZPt;
//	      double hadronicZDeltaR = TMath::Sqrt(pow(DeltaPhi(hadronicZcandidates[0]->phi(),hadronicZcandidates[1]->phi()),2)+pow(hadronicZcandidates[0]->eta()-hadronicZcandidates[1]->eta(),2));
//	      double hadronicZrecDeltaR = TMath::Sqrt(pow(DeltaPhi(ZjetVec[0].phi(),ZjetVec[1].phi()),2)+pow(ZjetVec[0].eta()-ZjetVec[1].eta(),2));
//
//	      if ( ZjetVec[0].eta() == 0. ||
//		   ZjetVec[1].eta() == 0. ||
//		   TjetVec[0].eta() == 0. ||
//		   TjetVec[1].eta() == 0. ||
//		   ZjetVec[0].et() < jetEtCut_ ||
//		   ZjetVec[1].et() < jetEtCut_ ||
//		   TjetVec[0].et() < jetEtCut_ ||
//		   TjetVec[1].et() < jetEtCut_ ) {
//		//		std::cout << "WARNING!!!" << std::endl;
//		//		std::cout << "ZjetVec[0].et: " << ZjetVec[0].et() << " eta: " << ZjetVec[0].eta() << " phi: " << ZjetVec[0].phi() << std::endl;
//		//		std::cout << "ZjetVec[1].et: " << ZjetVec[1].et() << " eta: " << ZjetVec[1].eta() << " phi: " << ZjetVec[1].phi() << std::endl;
//		//		std::cout << "ZlepVec[0].pt: " << ZlepVec[0].pt() << " eta: " << ZlepVec[0].eta() << " phi: " << ZlepVec[0].phi() << std::endl;
//		//		std::cout << "ZlepVec[1].pt: " << ZlepVec[1].pt() << " eta: " << ZlepVec[1].eta() << " phi: " << ZlepVec[1].phi() << std::endl;
//		//		std::cout << "TjetVec[0].et: " << TjetVec[0].et() << " eta: " << TjetVec[0].eta() << " phi: " << TjetVec[0].phi() << std::endl;
//		//		std::cout << "TjetVec[1].et: " << TjetVec[1].et() << " eta: " << TjetVec[1].eta() << " phi: " << TjetVec[1].phi() << std::endl;
//		//std::cout << "Zjet.et: " << Zjet->pt() << std::endl;
//		//std::cout << "Zlep.et: " << Zlep->pt() << std::endl;
//		//std::cout << "Tjet.et: " << Tjet->pt() << std::endl;
//		//std::cout << "Indices: " << TjetInd[0] << " " << TjetInd[1] << " " << ZjetInd[0] << " " << ZjetInd[1] << std::endl;
//	      }
//
//	      if ( ZlepVec[0].pt() != 0. && ZlepVec[1].pt() !=0. ) {
//
//		tmvaMatchingCounter_++;
//
//		// Study alternativekinematics of true ZZqq system
//		// ------------------------------------
//		// qq: tag jets
//		// ------------
//		double dphiqq = DeltaPhi(TjetVec[0].phi(),TjetVec[1].phi());
//		double detaqq = fabs(TjetVec[0].eta()-TjetVec[1].eta());
//		double ptminqq = TjetVec[0].et();
//		double ptmaxqq = TjetVec[1].et();
//		if ( TjetVec[1].et()<ptminqq ) {
//		  ptminqq = TjetVec[1].et();
//		  ptmaxqq = TjetVec[0].et();
//		}
//		double etaminqq = fabs(TjetVec[0].eta());
//		double etamaxqq = fabs(TjetVec[1].eta());
//		if ( fabs(TjetVec[1].eta())<etaminqq ) {
//		  etaminqq = fabs(TjetVec[1].eta());
//		  etamaxqq = fabs(TjetVec[0].eta());
//		}
//		eventSignalVariablesVector_.push_back(TjetVec[0].et());
//		eventSignalVariablesVector_.push_back(TjetVec[1].et());
//		eventSignalVariablesVector_.push_back(TjetVec[0].eta());
//		eventSignalVariablesVector_.push_back(TjetVec[1].eta());
//		eventSignalVariablesVector_.push_back(TjetVec[0].phi());
//		eventSignalVariablesVector_.push_back(TjetVec[1].phi());
//		eventSignalVariablesVector_.push_back(dphiqq);
//		eventSignalVariablesVector_.push_back(detaqq);
//		eventSignalVariablesVector_.push_back(ptminqq);
//		eventSignalVariablesVector_.push_back(ptmaxqq);
//		eventSignalVariablesVector_.push_back(etaminqq);
//		eventSignalVariablesVector_.push_back(etamaxqq);
//
//		// hh: Z jets
//		// ----------
//		double dphihh = DeltaPhi(ZjetVec[0].phi(),ZjetVec[1].phi());
//		double detahh = fabs(ZjetVec[0].eta()-ZjetVec[1].eta());
//		double ptminhh = ZjetVec[0].et();
//		double ptmaxhh = ZjetVec[1].et();
//		if ( ZjetVec[1].et()<ptminhh ) {
//		  ptminhh = ZjetVec[1].et();
//		  ptmaxhh = ZjetVec[0].et();
//		}
//		double etaminhh = fabs(ZjetVec[0].eta());
//		double etamaxhh = fabs(ZjetVec[1].eta());
//		if ( fabs(ZjetVec[1].eta())<etaminhh ) {
//		  etaminhh = fabs(ZjetVec[1].eta());
//		  etamaxhh = fabs(ZjetVec[0].eta());
//		}
//		eventSignalVariablesVector_.push_back(ZjetVec[0].et());
//		eventSignalVariablesVector_.push_back(ZjetVec[1].et());
//		eventSignalVariablesVector_.push_back(ZjetVec[0].eta());
//		eventSignalVariablesVector_.push_back(ZjetVec[1].eta());
//		eventSignalVariablesVector_.push_back(ZjetVec[0].phi());
//		eventSignalVariablesVector_.push_back(ZjetVec[1].phi());
//		eventSignalVariablesVector_.push_back(dphihh);
//		eventSignalVariablesVector_.push_back(detahh);
//		eventSignalVariablesVector_.push_back(ptminhh);
//		eventSignalVariablesVector_.push_back(ptmaxhh);
//		eventSignalVariablesVector_.push_back(etaminhh);
//		eventSignalVariablesVector_.push_back(etamaxhh);
//
//		// qq system
//		// ---------
//		double ptqq = Tjet->pt();
//		double mqq = Tjet->mass();
//		double etaqq = fabs(Tjet->eta());
//		// hh system
//		// ---------
//		double pthh = Zjet->pt();
//		double mhh = Zjet->mass();
//		double etahh = fabs(Zjet->eta());
//		// ll system
//		// ---------
//		double ptll = Zlep->pt();
//		double mll = Zlep->mass();
//		double etall = fabs(Zlep->eta());
//		//std::cout << "Systems Pt:   " << ptqq  << " " << pthh  << " " << ptll  << std::endl;
//		//std::cout << "Systems eta:  " << etaqq << " " << etahh << " " << etall << std::endl;
//		//std::cout << "Systems mass: " << mqq   << " " << mhh   << " " << mll   << std::endl;
//		eventSignalVariablesVector_.push_back(ptqq);
//		eventSignalVariablesVector_.push_back(mqq);
//		eventSignalVariablesVector_.push_back(etaqq);
//		eventSignalVariablesVector_.push_back(pthh);
//		eventSignalVariablesVector_.push_back(mhh);
//		eventSignalVariablesVector_.push_back(etahh);
//		eventSignalVariablesVector_.push_back(ptll);
//		eventSignalVariablesVector_.push_back(mll);
//		eventSignalVariablesVector_.push_back(etall);
//
//		// 2-particle vars
//		// ---------------
//		double dphiTjetZjet = DeltaPhi(Tjet,Zjet);
//		double dphiTjetZlep = DeltaPhi(Tjet,Zlep);
//		double dphiminTZ = dphiTjetZjet;
//		if ( dphiminTZ>dphiTjetZlep ) dphiminTZ = dphiTjetZlep;
//		double detaTjetZjet = fabs(Tjet->eta()-Zjet->eta());
//		double detaTjetZlep = fabs(Tjet->eta()-Zlep->eta());
//		double detaminTZ = detaTjetZjet;
//		if ( detaminTZ>detaTjetZlep ) detaminTZ = detaTjetZlep;
//		double dphiZjetZlep = DeltaPhi(Zjet,Zlep);
//		double detaZjetZlep = fabs(Zjet->eta()-Zlep->eta());
//		double massTjetZjet = (Tjet->p4()+Zjet->p4()).mass();
//		double massTjetZlep = (Tjet->p4()+Zlep->p4()).mass();
//		double massZjetZlep = (Zjet->p4()+Zlep->p4()).mass();
//		eventSignalVariablesVector_.push_back(dphiTjetZjet);
//		eventSignalVariablesVector_.push_back(dphiTjetZlep);
//		eventSignalVariablesVector_.push_back(dphiminTZ);
//		eventSignalVariablesVector_.push_back(detaTjetZjet);
//		eventSignalVariablesVector_.push_back(detaTjetZlep);
//		eventSignalVariablesVector_.push_back(detaminTZ);
//		eventSignalVariablesVector_.push_back(dphiZjetZlep);
//		eventSignalVariablesVector_.push_back(detaZjetZlep);
//		eventSignalVariablesVector_.push_back(massTjetZjet);
//		eventSignalVariablesVector_.push_back(massTjetZlep);
//		eventSignalVariablesVector_.push_back(massZjetZlep);
//
//		// 3-particle vars
//		// ---------------
//		double massTZZ = (Tjet->p4()+Zjet->p4()+Zlep->p4()).mass();
//		double etaTZZ = fabs((Tjet->p4()+Zjet->p4()+Zlep->p4()).eta());
//		double ptTZZ = (Tjet->p4()+Zjet->p4()+Zlep->p4()).pt();
//		eventSignalVariablesVector_.push_back(massTZZ);
//		eventSignalVariablesVector_.push_back(etaTZZ);
//		eventSignalVariablesVector_.push_back(ptTZZ);
//
//                // Store the event number
//                eventSignalVariablesVector_.push_back(eventcounter_);
//
//                // Variable to hold the value of the discriminant for the true combination
//                double discriminantSignal = 0.;
//
//                if( writeTMVA_ ) {
//                  // Fill the tree for the TMVA
//                  tmvaSignalTreeWriterPtr_->fill(eventSignalVariablesVector_);
//                }
//                else {
//                  // Fill the values in the array used by the reader
//                  vector<Float_t>::const_iterator varIt = eventSignalVariablesVector_.begin();
//                  int iVarSignal = 0;
//                  for( ; varIt != eventSignalVariablesVector_.end(); ++varIt, ++iVarSignal ) {
//                    variables_[iVarSignal] = *varIt;
//                  }
//                  discriminantSignal = reader_->EvaluateMVA( "BDT method" );
//                }
//                eventSignalVariablesVector_.clear();
//
//                // Vector to hold the weights of all the combinatorics
//                vector<double> discriminantBackground;
//
//		//std::cout << "Total system mass: " << massTZZ << " Eta: " << etaTZZ << " Pt: " << ptTZZ << std::endl;
//
//		// Loop on true combination and all false ones, and derive alternativekinematics of configuration
//		// -----------------------------------------------------------------------------------
//		//	    delete Tjet;
//		//	    delete Zjet;
//		int i1=0;
//		int i2=0;
//		int i3=0;
//		int i4=0;
//		std::vector<OfflineJet>::const_iterator jet1 = aboveEtCutJetVec.begin();
//		for ( ; jet1 != aboveEtCutJetVec.end()-1; ++jet1, i1++ ) {
//		  std::vector<OfflineJet>::const_iterator jet2 = jet1+1;
//		  for ( ; jet2 != aboveEtCutJetVec.end(); ++jet2, i2++) {
//		    if ( jet2==jet1 ) continue;
//		    std::vector<OfflineJet>::const_iterator jet3 = aboveEtCutJetVec.begin();
//		    for ( ; jet3 != aboveEtCutJetVec.end()-1; ++jet3, i3++) {
//		      if ( jet3==jet1 || jet3==jet2 ) continue;
//		      std::vector<OfflineJet>::const_iterator jet4 = jet3+1;
//		      for ( ; jet4 != aboveEtCutJetVec.end(); ++jet4, i4++) {
//			if ( jet4==jet1 || jet4==jet2 || jet4==jet3 ) continue;
//			// Remove true combo and homologues
//			// --------------------------------
//			if (i1==TjetInd[0] && i2==TjetInd[1] && i3==ZjetInd[0] && i4==ZjetInd[1]) continue;
//			if (i1==TjetInd[1] && i2==TjetInd[0] && i3==ZjetInd[0] && i4==ZjetInd[1]) continue;
//			if (i1==TjetInd[0] && i2==TjetInd[1] && i3==ZjetInd[1] && i4==ZjetInd[0]) continue;
//			if (i1==TjetInd[1] && i2==TjetInd[0] && i3==ZjetInd[1] && i4==ZjetInd[0]) continue;
//			// Ok, here we have a sensible quadruplet
//			// --------------------------------------
//			//std::cout << "Combination " << i1 << i2 << i3 << i4 << std::endl;
//			Particle *Tjet = new Particle (0,jet1->p4()+jet2->p4(), jet1->vertex(),0,0,true);
//			Particle *Zjet = new Particle (0,jet3->p4()+jet4->p4(), jet3->vertex(),0,0,true);
//
//			// Study alternativekinematics of random ZZqq systems
//			// ---------------------------------------
//			// qq: tag jets
//			// ------------
//			double f_dphiqq = DeltaPhi(jet1->phi(),jet2->phi());
//			double f_detaqq = fabs(jet1->eta()-jet2->eta());
//			double f_ptminqq = jet1->et();
//			double f_ptmaxqq = jet2->et();
//			if ( jet2->et()<f_ptminqq ) {
//			  f_ptminqq = jet2->et();
//			  f_ptmaxqq = jet1->et();
//			}
//			double f_etaminqq = fabs(jet1->eta());
//			double f_etamaxqq = fabs(jet2->eta());
//			if ( fabs(jet2->eta())<f_etaminqq ) {
//			  f_etaminqq = fabs(jet2->eta());
//			  f_etamaxqq = fabs(jet1->eta());
//			}
//			eventCombinVariablesVector_.push_back(jet1->et());
//			eventCombinVariablesVector_.push_back(jet2->et());
//			eventCombinVariablesVector_.push_back(jet1->eta());
//			eventCombinVariablesVector_.push_back(jet2->eta());
//			eventCombinVariablesVector_.push_back(jet1->phi());
//			eventCombinVariablesVector_.push_back(jet2->phi());
//			eventCombinVariablesVector_.push_back(f_dphiqq);
//			eventCombinVariablesVector_.push_back(f_detaqq);
//			eventCombinVariablesVector_.push_back(f_ptminqq);
//			eventCombinVariablesVector_.push_back(f_ptmaxqq);
//			eventCombinVariablesVector_.push_back(f_etaminqq);
//			eventCombinVariablesVector_.push_back(f_etamaxqq);
//
//			// hh: Z jets
//			// ----------
//			double f_dphihh = DeltaPhi(jet3->phi(),jet4->phi());
//			double f_detahh = fabs(jet3->eta()-jet4->eta());
//			double f_ptminhh = jet3->et();
//			double f_ptmaxhh = jet4->et();
//			if ( jet4->et()<f_ptminhh ) {
//			  f_ptminhh = jet4->et();
//			  f_ptmaxhh = jet3->et();
//			}
//			double f_etaminhh = fabs(jet3->eta());
//			double f_etamaxhh = fabs(jet4->eta());
//			if ( fabs(jet4->eta())<f_etaminhh ) {
//			  f_etaminhh = fabs(jet4->eta());
//			  f_etamaxhh = fabs(jet3->eta());
//			}
//			eventCombinVariablesVector_.push_back(jet3->et());
//			eventCombinVariablesVector_.push_back(jet4->et());
//			eventCombinVariablesVector_.push_back(jet3->eta());
//			eventCombinVariablesVector_.push_back(jet4->eta());
//			eventCombinVariablesVector_.push_back(jet3->phi());
//			eventCombinVariablesVector_.push_back(jet4->phi());
//			eventCombinVariablesVector_.push_back(f_dphihh);
//			eventCombinVariablesVector_.push_back(f_detahh);
//			eventCombinVariablesVector_.push_back(f_ptminhh);
//			eventCombinVariablesVector_.push_back(f_ptmaxhh);
//			eventCombinVariablesVector_.push_back(f_etaminhh);
//			eventCombinVariablesVector_.push_back(f_etamaxhh);
//
//			// qq system
//			// ---------
//			double f_ptqq = (Tjet->p4()).pt();
//			double f_mqq = (Tjet->p4()).mass();
//			double f_etaqq = fabs((Tjet->p4()).eta());
//			// hh system
//			// ---------
//			double f_pthh = (Zjet->p4()).pt();
//			double f_mhh = (Zjet->p4()).mass();
//			double f_etahh = fabs((Zjet->p4()).eta());
//			// ll system
//			// ---------
//			double f_ptll = (Zlep->p4()).pt();
//			double f_mll = Zlep->mass();
//			double f_etall = fabs(Zlep->eta());
//			eventCombinVariablesVector_.push_back(f_ptqq);
//			eventCombinVariablesVector_.push_back(f_mqq);
//			eventCombinVariablesVector_.push_back(f_etaqq);
//			eventCombinVariablesVector_.push_back(f_pthh);
//			eventCombinVariablesVector_.push_back(f_mhh);
//			eventCombinVariablesVector_.push_back(f_etahh);
//			eventCombinVariablesVector_.push_back(f_ptll);
//			eventCombinVariablesVector_.push_back(f_mll);
//			eventCombinVariablesVector_.push_back(f_etall);
//			// 2-particle vars
//			// ---------------
//			double f_dphiTjetZjet = DeltaPhi((Tjet->p4()).phi(),(Zjet->p4()).phi());
//			double f_dphiTjetZlep = DeltaPhi((Tjet->p4()).phi(),Zlep->phi());
//			double f_dphiminTZ = f_dphiTjetZjet;
//			if ( f_dphiminTZ>f_dphiTjetZlep ) f_dphiminTZ = f_dphiTjetZlep;
//			double f_detaTjetZjet = fabs((Tjet->p4()).eta()-(Zjet->p4()).eta());
//			double f_detaTjetZlep = fabs((Tjet->p4()).eta()-Zlep->eta());
//			double f_detaminTZ = f_detaTjetZjet;
//			if ( f_detaminTZ>f_detaTjetZlep ) f_detaminTZ = f_detaTjetZlep;
//			double f_dphiZjetZlep = DeltaPhi((Zjet->p4()).phi(),Zlep->phi());
//			double f_detaZjetZlep = fabs((Zjet->p4()).eta()-Zlep->eta());
//			double f_massTjetZjet = (Tjet->p4()+Zjet->p4()).mass();
//			double f_massTjetZlep = (Tjet->p4()+Zlep->p4()).mass();
//			double f_massZjetZlep = (Zjet->p4()+Zlep->p4()).mass();
//			eventCombinVariablesVector_.push_back(f_dphiTjetZjet);
//			eventCombinVariablesVector_.push_back(f_dphiTjetZlep);
//			eventCombinVariablesVector_.push_back(f_dphiminTZ);
//			eventCombinVariablesVector_.push_back(f_detaTjetZjet);
//			eventCombinVariablesVector_.push_back(f_detaTjetZlep);
//			eventCombinVariablesVector_.push_back(f_detaminTZ);
//			eventCombinVariablesVector_.push_back(f_dphiZjetZlep);
//			eventCombinVariablesVector_.push_back(f_detaZjetZlep);
//			eventCombinVariablesVector_.push_back(f_massTjetZjet);
//			eventCombinVariablesVector_.push_back(f_massTjetZlep);
//			eventCombinVariablesVector_.push_back(f_massZjetZlep);
//			// 3-particle vars
//			// ---------------
//			double f_massTZZ = (Tjet->p4()+Zjet->p4()+Zlep->p4()).mass();
//			double f_etaTZZ = fabs((Tjet->p4()+Zjet->p4()+Zlep->p4()).eta());
//			double f_ptTZZ = (Tjet->p4()+Zjet->p4()+Zlep->p4()).pt();
//			eventCombinVariablesVector_.push_back(f_massTZZ);
//			eventCombinVariablesVector_.push_back(f_etaTZZ);
//			eventCombinVariablesVector_.push_back(f_ptTZZ);
//		    
//                        // Store the event number
//                        eventCombinVariablesVector_.push_back(eventcounter_);
//
//                        if( writeTMVA_ ) {
//                          // Fill the tree for the TMVA
//                          tmvaCombinTreeWriterPtr_->fill(eventCombinVariablesVector_);
//                        }
//                        else {
//                          // Fill the values in the array used by the reader
//                          vector<Float_t>::const_iterator varIt = eventSignalVariablesVector_.begin();
//                          int iVarBackground = 0;
//                          for( ; varIt != eventCombinVariablesVector_.end(); ++varIt, ++iVarBackground ) {
//                            variables_[iVarBackground] = *varIt;
//                          }
//                          discriminantBackground.push_back( reader_->EvaluateMVA( "BDT method" ) );
//                        }
//                        eventCombinVariablesVector_.clear();
//		      }
//		    }
//		  }
//		}
//	      
//                if( !writeTMVA_ ) {
//                  // After the loop on all combinations sort the vector of background discriminants
//                  sort( discriminantBackground.rbegin(), discriminantBackground.rend() );
//                  int i = 0;
//                  for( vector<double>::const_iterator discr = discriminantBackground.begin(); discr != discriminantBackground.end(); ++discr, ++i ) {
//                    cout << "discr["<<i<<"] = " << *discr << endl;
//                  }
//                }
//
//	      
//
//		// count combination due to jet multiplicity
//		jet1 = aboveEtCutJetVec.begin();
//		for ( ; jet1 != aboveEtCutJetVec.end()-1; ++jet1, i1++ ) {
//		  std::vector<OfflineJet>::const_iterator jet2 = jet1+1;
//		  for ( ; jet2 != aboveEtCutJetVec.end(); ++jet2, i2++) {
//		    std::vector<OfflineJet>::const_iterator jet3 = aboveEtCutJetVec.begin();
//		    for ( ; jet3 != aboveEtCutJetVec.end()-1; ++jet3, i3++) {
//		      if ( jet3==jet1 || jet3==jet2 ) continue;
//		      std::vector<OfflineJet>::const_iterator jet4 = jet3+1;
//		      for ( ; jet4 != aboveEtCutJetVec.end(); ++jet4, i4++) {
//			if ( jet4==jet1 || jet4==jet2 ) continue;
//			combinatorialCounter_++;
//		      }
//		    }
//		  }
//		}
//	      }
//
//	      // =============================================================================================
//	      // derive alternativekinematics of configuration: eta cut
//	      // ------------------------------------------------------
//
//	      // eT sorting
//	      // ----------
//	      std::vector<OfflineJet> etSortedJetVec;
//	      std::vector<int>        etSortedJetIndexVec;
//	      std::vector<bool> stored;
//	      
//	      for ( unsigned int times = 0; times < aboveEtCutJetVec.size(); times++ ) 
//		stored.push_back(kFALSE);
//	      
//	      for ( unsigned int times = 0; times < aboveEtCutJetVec.size(); times++ ) {
//		double maxEt = 0.;
//		unsigned int maxIndex = 0;
//		std::vector<OfflineJet>::const_iterator aboveEtCutJetVec_itr = aboveEtCutJetVec.begin(); 
//		unsigned int index = 0;
//		for ( ; aboveEtCutJetVec_itr != aboveEtCutJetVec.end(); ++aboveEtCutJetVec_itr, index++ ) {
//		  if ( stored[index] ) continue;
//		  if ( aboveEtCutJetVec_itr->et() >= maxEt ) {
//		    maxEt = aboveEtCutJetVec_itr->et();
//		    maxIndex = index;
//		  }
//		}
//		stored[maxIndex] = kTRUE;
//		etSortedJetVec.push_back(aboveEtCutJetVec[maxIndex]);
//		etSortedJetIndexVec.push_back(maxIndex);
//	      }
//	      
//	      //	      for ( unsigned int times = 0; times < aboveEtCutJetVec.size(); times++ ) {
//	      //		std::cout << "aboveEtCutJetVec[" << times << "]: " << aboveEtCutJetVec[times].et()
//	      //			  << " --> etSortedJetVec[" << times << "]: " << etSortedJetVec[times].et() << std::endl;
//	      //		std::cout << "etSortedJetIndexVec[" << times << "]: " << etSortedJetIndexVec[times] << " -->eta: " << etSortedJetVec[times].eta() << std::endl;
//	      //	      }
//	      
//	      std::vector<OfflineJet> centralJetsVector;
//	      std::vector<int>        centralJetsIndexVector;
//	      std::vector<OfflineJet> forwardJetsVector;
//	      std::vector<int>        forwardJetsIndexVector;
//	      std::vector<OfflineJet>::const_iterator etSortedJetVec_itr      = etSortedJetVec.begin();
//	      std::vector<int>::const_iterator        etSortedJetIndexVec_itr = etSortedJetIndexVec.begin();
//	      for ( ; etSortedJetVec_itr != etSortedJetVec.end(); ++etSortedJetVec_itr, ++etSortedJetIndexVec_itr ) {
//		if ( fabs( etSortedJetVec_itr->eta() ) <= 2. ) {
//		  centralJetsVector.push_back(*etSortedJetVec_itr);
//		  centralJetsIndexVector.push_back(*etSortedJetIndexVec_itr);
//		}
//		else {
//		  forwardJetsVector.push_back(*etSortedJetVec_itr);
//		  forwardJetsIndexVector.push_back(*etSortedJetIndexVec_itr);
//		}
//	      }
//
//	      if ( centralJetsVector.size() >= 2 && forwardJetsVector.size() >= 2 ) {
//
//		//		for (unsigned int index = 0; index < centralJetsVector.size(); index++ ) 
//		//		  std::cout << "centralJetsIndexVector[" << index << "]: " << centralJetsIndexVector[index] << ": et: " << centralJetsVector[index].et() << std::endl;
//		//		for (unsigned int index = 0 ; index < forwardJetsIndexVector.size(); index++ )
//		//		  std::cout << "forwardJetsIndexVector[" << index << "]: " << forwardJetsIndexVector[index] << ": et: " << forwardJetsVector[index].et() << std::endl;		
//
//		std::vector<OfflineJet> TjetsVector;
//		std::vector<int>        TjetsIndexVector;
//		TjetsVector.push_back(forwardJetsVector[0]);
//		TjetsIndexVector.push_back(forwardJetsIndexVector[0]);
//		double eta1 = forwardJetsVector[0].eta();
//		std::vector<OfflineJet>::const_iterator forwardJetsVector_itr2      = forwardJetsVector.begin()+1;
//		std::vector<int>::const_iterator        forwardJetsIndexVector_itr2 = forwardJetsIndexVector.begin()+1;
//		for ( ; forwardJetsVector_itr2 != forwardJetsVector.end(); ++forwardJetsVector_itr2, ++forwardJetsIndexVector_itr2 ){
//		  if ( eta1*forwardJetsVector_itr2->eta() < 0 ) {
//		    TjetsVector.push_back(*forwardJetsVector_itr2);
//		    TjetsIndexVector.push_back(*forwardJetsIndexVector_itr2);
//		    break;
//		  }		      
//		}
//		//		for (unsigned int index = 0; index < TjetsIndexVector.size(); index++ )
//		//		  std::cout << "TjetsIndexVector[" << index << "]: " << TjetsIndexVector[index] << std::endl;
//
//		if ( TjetsVector.size() == 2 ) {
//		  std::vector<OfflineJet> ZjetsVector;
//		  std::vector<int>        ZjetsIndexVector;
//		  for (unsigned int index = 0; index < 2; index++) {
//		    ZjetsVector.push_back(centralJetsVector[index]);
//		    ZjetsIndexVector.push_back(centralJetsIndexVector[index]);
//		  }
//
//		  if ( ZjetsVector.size() == 2 ) {
//		    if (
//			( (TjetsIndexVector[0]==TjetInd[0] && TjetsIndexVector[1]==TjetInd[1] ) ||
//			  (TjetsIndexVector[1]==TjetInd[0] && TjetsIndexVector[0]==TjetInd[1] ) ) &&
//			( (ZjetsIndexVector[0]==ZjetInd[0] && ZjetsIndexVector[1]==ZjetInd[1] ) ||
//			  (ZjetsIndexVector[1]==ZjetInd[0] && ZjetsIndexVector[0]==ZjetInd[1] ) )
//			) etaCutMatchingCounter_++;
//		  
//		    Particle *Tjet =  new Particle (0,TjetsVector[0].p4()+TjetsVector[1].p4(), TjetsVector[0].vertex(),0,0,true);
//		    Particle *Zjet =  new Particle (0,ZjetsVector[0].p4()+ZjetsVector[1].p4(), ZjetsVector[0].vertex(),0,0,true);
//		  
//		    // Study alternativekinematics of ZZqq systems
//		    // -------------------------------------------
//		    // qq: tag jets
//		    // ------------
//		    double etaCut_dphiqq = DeltaPhi(TjetsVector[0].phi(),TjetsVector[1].phi());
//		    double etaCut_detaqq = fabs(TjetsVector[0].eta()-TjetsVector[1].eta());
//		    double etaCut_ptminqq = TjetsVector[0].et();
//		    double etaCut_ptmaxqq = TjetsVector[1].et();
//		    if ( TjetsVector[1].et()<etaCut_ptminqq ) {
//		      etaCut_ptminqq = TjetsVector[1].et();
//		      etaCut_ptmaxqq = TjetsVector[0].et();
//		    }
//		    double etaCut_etaminqq = fabs(TjetsVector[0].eta());
//		    double etaCut_etamaxqq = fabs(TjetsVector[1].eta());
//		    if ( fabs(TjetsVector[1].eta())<etaCut_etaminqq ) {
//		      etaCut_etaminqq = fabs(TjetsVector[1].eta());
//		      etaCut_etamaxqq = fabs(TjetsVector[0].eta());
//		    }
//		    eventEtaCutVariablesVector_.push_back(TjetsVector[0].et());
//		    eventEtaCutVariablesVector_.push_back(TjetsVector[1].et());
//		    eventEtaCutVariablesVector_.push_back(TjetsVector[0].eta());
//		    eventEtaCutVariablesVector_.push_back(TjetsVector[1].eta());
//		    eventEtaCutVariablesVector_.push_back(TjetsVector[0].phi());
//		    eventEtaCutVariablesVector_.push_back(TjetsVector[1].phi());
//		    eventEtaCutVariablesVector_.push_back(etaCut_dphiqq);
//		    eventEtaCutVariablesVector_.push_back(etaCut_detaqq);
//		    eventEtaCutVariablesVector_.push_back(etaCut_ptminqq);
//		    eventEtaCutVariablesVector_.push_back(etaCut_ptmaxqq);
//		    eventEtaCutVariablesVector_.push_back(etaCut_etaminqq);
//		    eventEtaCutVariablesVector_.push_back(etaCut_etamaxqq);
//		    
//		    // hh: Z jets
//		    // ----------
//		    double etaCut_dphihh = DeltaPhi(ZjetsVector[0].phi(),ZjetsVector[1].phi());
//		    double etaCut_detahh = fabs(ZjetsVector[0].eta()-ZjetsVector[1].eta());
//		    double etaCut_ptminhh = ZjetsVector[0].et();
//		    double etaCut_ptmaxhh = ZjetsVector[1].et();
//		    if ( ZjetsVector[1].et()<etaCut_ptminhh ) {
//		      etaCut_ptminhh = ZjetsVector[1].et();
//		      etaCut_ptmaxhh = ZjetsVector[0].et();
//		    }
//		    double etaCut_etaminhh = fabs(ZjetsVector[0].eta());
//		    double etaCut_etamaxhh = fabs(ZjetsVector[1].eta());
//		    if ( fabs(ZjetsVector[1].eta())<etaCut_etaminhh ) {
//		      etaCut_etaminhh = fabs(ZjetsVector[1].eta());
//		      etaCut_etamaxhh = fabs(ZjetsVector[0].eta());
//		    }
//		    eventEtaCutVariablesVector_.push_back(ZjetsVector[0].et());
//		    eventEtaCutVariablesVector_.push_back(ZjetsVector[1].et());
//		    eventEtaCutVariablesVector_.push_back(ZjetsVector[0].eta());
//		    eventEtaCutVariablesVector_.push_back(ZjetsVector[1].eta());
//		    eventEtaCutVariablesVector_.push_back(ZjetsVector[0].phi());
//		    eventEtaCutVariablesVector_.push_back(ZjetsVector[1].phi());
//		    eventEtaCutVariablesVector_.push_back(etaCut_dphihh);
//		    eventEtaCutVariablesVector_.push_back(etaCut_detahh);
//		    eventEtaCutVariablesVector_.push_back(etaCut_ptminhh);
//		    eventEtaCutVariablesVector_.push_back(etaCut_ptmaxhh);
//		    eventEtaCutVariablesVector_.push_back(etaCut_etaminhh);
//		    eventEtaCutVariablesVector_.push_back(etaCut_etamaxhh);
//		    
//		    // qq system
//		    // ---------
//		    double etaCut_ptqq = (Tjet->p4()).pt();
//		    double etaCut_mqq = (Tjet->p4()).mass();
//		    double etaCut_etaqq = fabs((Tjet->p4()).eta());
//		    // hh systems
//		    // ---------
//		    double etaCut_pthh = (Zjet->p4()).pt();
//		    double etaCut_mhh = (Zjet->p4()).mass();
//		    double etaCut_etahh = fabs((Zjet->p4()).eta());
//		    // ll systems
//		    // ---------
//		    double etaCut_ptll = (Zlep->p4()).pt();
//		    double etaCut_mll = Zlep->mass();
//		    double etaCut_etall = fabs(Zlep->eta());
//		    eventEtaCutVariablesVector_.push_back(etaCut_ptqq);
//		    eventEtaCutVariablesVector_.push_back(etaCut_mqq);
//		    eventEtaCutVariablesVector_.push_back(etaCut_etaqq);
//		    eventEtaCutVariablesVector_.push_back(etaCut_pthh);
//		    eventEtaCutVariablesVector_.push_back(etaCut_mhh);
//		    eventEtaCutVariablesVector_.push_back(etaCut_etahh);
//		    eventEtaCutVariablesVector_.push_back(etaCut_ptll);
//		    eventEtaCutVariablesVector_.push_back(etaCut_mll);
//		    eventEtaCutVariablesVector_.push_back(etaCut_etall);
//		    // 2-particle vars
//		    // ---------------
//		    double etaCut_dphiTjetZjet = DeltaPhi((Tjet->p4()).phi(),(Zjet->p4()).phi());
//		    double etaCut_dphiTjetZlep = DeltaPhi((Tjet->p4()).phi(),Zlep->phi());
//		    double etaCut_dphiminTZ = etaCut_dphiTjetZjet;
//		    if ( etaCut_dphiminTZ>etaCut_dphiTjetZlep ) etaCut_dphiminTZ = etaCut_dphiTjetZlep;
//		    double etaCut_detaTjetZjet = fabs((Tjet->p4()).eta()-(Zjet->p4()).eta());
//		    double etaCut_detaTjetZlep = fabs((Tjet->p4()).eta()-Zlep->eta());
//		    double etaCut_detaminTZ = etaCut_detaTjetZjet;
//		    if ( etaCut_detaminTZ>etaCut_detaTjetZlep ) etaCut_detaminTZ = etaCut_detaTjetZlep;
//		    double etaCut_dphiZjetZlep = DeltaPhi((Zjet->p4()).phi(),Zlep->phi());
//		    double etaCut_detaZjetZlep = fabs((Zjet->p4()).eta()-Zlep->eta());
//		    double etaCut_massTjetZjet = (Tjet->p4()+Zjet->p4()).mass();
//		    double etaCut_massTjetZlep = (Tjet->p4()+Zlep->p4()).mass();
//		    double etaCut_massZjetZlep = (Zjet->p4()+Zlep->p4()).mass();
//		    eventEtaCutVariablesVector_.push_back(etaCut_dphiTjetZjet);
//		    eventEtaCutVariablesVector_.push_back(etaCut_dphiTjetZlep);
//		    eventEtaCutVariablesVector_.push_back(etaCut_dphiminTZ);
//		    eventEtaCutVariablesVector_.push_back(etaCut_detaTjetZjet);
//		    eventEtaCutVariablesVector_.push_back(etaCut_detaTjetZlep);
//		    eventEtaCutVariablesVector_.push_back(etaCut_detaminTZ);
//		    eventEtaCutVariablesVector_.push_back(etaCut_dphiZjetZlep);
//		    eventEtaCutVariablesVector_.push_back(etaCut_detaZjetZlep);
//		    eventEtaCutVariablesVector_.push_back(etaCut_massTjetZjet);
//		    eventEtaCutVariablesVector_.push_back(etaCut_massTjetZlep);
//		    eventEtaCutVariablesVector_.push_back(etaCut_massZjetZlep);
//		    // 3-particle vars
//		    // ---------------
//		    double etaCut_massTZZ = (Tjet->p4()+Zjet->p4()+Zlep->p4()).mass();
//		    double etaCut_etaTZZ = fabs((Tjet->p4()+Zjet->p4()+Zlep->p4()).eta());
//		    double etaCut_ptTZZ = (Tjet->p4()+Zjet->p4()+Zlep->p4()).pt();
//		    eventEtaCutVariablesVector_.push_back(etaCut_massTZZ);
//		    eventEtaCutVariablesVector_.push_back(etaCut_etaTZZ);
//		    eventEtaCutVariablesVector_.push_back(etaCut_ptTZZ);
//		  
//		    // Stored the event number
//		    eventEtaCutVariablesVector_.push_back(eventcounter_);
//		    
//		    if( writeTMVA_ ) {
//		      // Fill the tree for the TMVA
//		      tmvaEtaCutTreeWriterPtr_->fill(eventEtaCutVariablesVector_);
//		    }
//		    /*
//		      else {
//		      // Fill the values in the array used by the reader
//		      vector<Float_t>::const_iterator varIt = eventSignalVariablesVector_.begin();
//		      int iVarBackground = 0;
//		      for( ; varIt != eventEtaCutVariablesVector_.end(); ++varIt, ++iVarBackground ) {
//		      variables_[iVarBackground] = *varIt;
//		      }
//		      discriminantBackground.push_back( reader_->EvaluateMVA( "BDT method" ) );
//		      }
//		    */
//		    eventEtaCutVariablesVector_.clear();
//		  }
//	      
//		  ZjetsVector.clear();
//		  ZjetsIndexVector.clear();
//		  // derive alternativekinematics of configuration: closest mjj to mZ
//		  // ----------------------------------------------------------------
//		  double mZres = 20.;
//		  double mHres = mHres_;
//		  unsigned int ZjetIndex1 = 0;
//		  unsigned int ZjetIndex2 = 1;
//		  double tmp_mZres = mZres;
//		  double tmp_chi2 = 1000.;
//		  std::vector<OfflineJet>::const_iterator centralJetsVector_itr1      = centralJetsVector.begin();
//		  std::vector<int>::const_iterator        centralJetsIndexVector_itr1 = centralJetsIndexVector.begin();
//		  unsigned int index1 = 0;
//		  for ( ; centralJetsVector_itr1 != centralJetsVector.end(); ++centralJetsVector_itr1, ++centralJetsIndexVector_itr1, index1++ ) {
//		    std::vector<OfflineJet>::const_iterator centralJetsVector_itr2      = centralJetsVector_itr1+1;
//		    std::vector<int>::const_iterator        centralJetsIndexVector_itr2 = centralJetsIndexVector_itr1+1;
//		    unsigned int index2 = index1+1;
//		    for ( ; centralJetsVector_itr2 != centralJetsVector.end(); ++centralJetsVector_itr2, ++centralJetsIndexVector_itr2, index2++ ) {
//		      double jjInvMass = (centralJetsVector_itr1->p4()+centralJetsVector_itr2->p4()).mass();
//		      //		      std::cout << "jjInvMass[" << *centralJetsIndexVector_itr1 << ";" << *centralJetsIndexVector_itr2 << "]: " << jjInvMass << std::endl;
//		      if ( fabs(jjInvMass-ZMass_) <= tmp_mZres ) {
//			tmp_mZres = fabs(jjInvMass-ZMass_);
//			ZjetIndex1 = index1;
//			ZjetIndex2 = index2;
//		      }
//		      double jjInvMassChi2 = pow((jjInvMass-ZMass_)/mZres,2);
//		      double lljjIinvMass = (centralJetsVector_itr1->p4()+centralJetsVector_itr2->p4()+Zlep->p4()).mass();
//		      double lljjIinvMassChi2 = pow((lljjIinvMass-mH_)/mHres,2);
//		      if ( jjInvMassChi2+lljjIinvMassChi2 < tmp_chi2 ) {
//			tmp_chi2 = jjInvMassChi2+lljjIinvMassChi2;
//			ZjetIndex1 = index1;
//			ZjetIndex2 = index2;
//			//			std::cout << "index1: " << *centralJetsIndexVector_itr1 << " <--> index2: " << *centralJetsIndexVector_itr2 << std::endl;
//		      }
//		    }
//		  }
//		
//		  ZjetsVector.push_back(centralJetsVector[ZjetIndex1]);
//		  ZjetsVector.push_back(centralJetsVector[ZjetIndex2]);
//		  ZjetsIndexVector.push_back(centralJetsIndexVector[ZjetIndex1]);
//		  ZjetsIndexVector.push_back(centralJetsIndexVector[ZjetIndex2]);
//		    
//		  //		  for (unsigned int index = 0; index < 2; index++) {
//		  //		    std::cout << "ZjetsIndexVector[" << index << "]: " << ZjetsIndexVector[index] 
//		  //			      << " <--> ZjetInd[" << index << "]: " << ZjetInd[index] << std::endl;
//		  //		    std::cout << "TjetsIndexVector[" << index << "]: " << TjetsIndexVector[index]
//		  //			      << " <--> TjetInd[" << index << "]: " << TjetInd[index] << std::endl;
//		  //		  }
//		
//		  if ( ZjetsVector.size() == 2 ) {
//		    if (
//			( (TjetsIndexVector[0]==TjetInd[0] && TjetsIndexVector[1]==TjetInd[1] ) ||
//			  (TjetsIndexVector[1]==TjetInd[0] && TjetsIndexVector[0]==TjetInd[1] ) ) &&
//			( (ZjetsIndexVector[0]==ZjetInd[0] && ZjetsIndexVector[1]==ZjetInd[1] ) ||
//			  (ZjetsIndexVector[1]==ZjetInd[0] && ZjetsIndexVector[0]==ZjetInd[1] ) ) ) mZMatchingCounter_++;
//
//		    Particle *Tjet =  new Particle (0,TjetsVector[0].p4()+TjetsVector[1].p4(), TjetsVector[0].vertex(),0,0,true);
//		    Particle *Zjet =  new Particle (0,ZjetsVector[0].p4()+ZjetsVector[1].p4(), ZjetsVector[0].vertex(),0,0,true);
//		      
//
//
//		    // Study alternativekinematics of random ZZqq systems
//		    // ---------------------------------------
//		    // qq: tag jets
//		    // ------------
//		    double mZcut_dphiqq = DeltaPhi(TjetsVector[0].phi(),TjetsVector[1].phi());
//		    double mZcut_detaqq = fabs(TjetsVector[0].eta()-TjetsVector[1].eta());
//		    double mZcut_ptminqq = TjetsVector[0].et();
//		    double mZcut_ptmaxqq = TjetsVector[1].et();
//		    if ( TjetsVector[1].et()<mZcut_ptminqq ) {
//		      mZcut_ptminqq = TjetsVector[1].et();
//		      mZcut_ptmaxqq = TjetsVector[2].et();
//		    }
//		    double mZcut_etaminqq = fabs(TjetsVector[0].eta());
//		    double mZcut_etamaxqq = fabs(TjetsVector[1].eta());
//		    if ( fabs(TjetsVector[1].eta())<mZcut_etaminqq ) {
//		      mZcut_etaminqq = fabs(TjetsVector[1].eta());
//		      mZcut_etamaxqq = fabs(TjetsVector[0].eta());
//		    }
//		    eventMZcutVariablesVector_.push_back(TjetsVector[0].et());
//		    eventMZcutVariablesVector_.push_back(TjetsVector[1].et());
//		    eventMZcutVariablesVector_.push_back(TjetsVector[0].eta());
//		    eventMZcutVariablesVector_.push_back(TjetsVector[1].eta());
//		    eventMZcutVariablesVector_.push_back(TjetsVector[0].phi());
//		    eventMZcutVariablesVector_.push_back(TjetsVector[1].phi());
//		    eventMZcutVariablesVector_.push_back(mZcut_dphiqq);
//		    eventMZcutVariablesVector_.push_back(mZcut_detaqq);
//		    eventMZcutVariablesVector_.push_back(mZcut_ptminqq);
//		    eventMZcutVariablesVector_.push_back(mZcut_ptmaxqq);
//		    eventMZcutVariablesVector_.push_back(mZcut_etaminqq);
//		    eventMZcutVariablesVector_.push_back(mZcut_etamaxqq);
//
//		    // hh: Z jets
//		    // ----------
//		    double mZcut_dphihh = DeltaPhi(ZjetsVector[0].phi(),ZjetsVector[1].phi());
//		    double mZcut_detahh = fabs(ZjetsVector[0].eta()-ZjetsVector[1].eta());
//		    double mZcut_ptminhh = ZjetsVector[0].et();
//		    double mZcut_ptmaxhh = ZjetsVector[1].et();
//		    if ( ZjetsVector[1].et()<mZcut_ptminhh ) {
//		      mZcut_ptminhh = ZjetsVector[1].et();
//		      mZcut_ptmaxhh = ZjetsVector[0].et();
//		    }
//		    double mZcut_etaminhh = fabs(ZjetsVector[0].eta());
//		    double mZcut_etamaxhh = fabs(ZjetsVector[1].eta());
//		    if ( fabs(ZjetsVector[1].eta())<mZcut_etaminhh ) {
//		      mZcut_etaminhh = fabs(ZjetsVector[1].eta());
//		      mZcut_etamaxhh = fabs(ZjetsVector[0].eta());
//		    }
//		    eventMZcutVariablesVector_.push_back(ZjetsVector[0].et());
//		    eventMZcutVariablesVector_.push_back(ZjetsVector[1].et());
//		    eventMZcutVariablesVector_.push_back(ZjetsVector[0].eta());
//		    eventMZcutVariablesVector_.push_back(ZjetsVector[1].eta());
//		    eventMZcutVariablesVector_.push_back(ZjetsVector[0].phi());
//		    eventMZcutVariablesVector_.push_back(ZjetsVector[1].phi());
//		    eventMZcutVariablesVector_.push_back(mZcut_dphihh);
//		    eventMZcutVariablesVector_.push_back(mZcut_detahh);
//		    eventMZcutVariablesVector_.push_back(mZcut_ptminhh);
//		    eventMZcutVariablesVector_.push_back(mZcut_ptmaxhh);
//		    eventMZcutVariablesVector_.push_back(mZcut_etaminhh);
//		    eventMZcutVariablesVector_.push_back(mZcut_etamaxhh);
//
//		    // qq system
//		    // ---------
//		    double mZcut_ptqq = (Tjet->p4()).pt();
//		    double mZcut_mqq = (Tjet->p4()).mass();
//		    double mZcut_etaqq = fabs((Tjet->p4()).eta());
//		    // hh system
//		    // ---------
//		    double mZcut_pthh = (Zjet->p4()).pt();
//		    double mZcut_mhh = (Zjet->p4()).mass();
//		    double mZcut_etahh = fabs((Zjet->p4()).eta());
//		    // ll system
//		    // ---------
//		    double mZcut_ptll = (Zlep->p4()).pt();
//		    double mZcut_mll = Zlep->mass();
//		    double mZcut_etall = fabs(Zlep->eta());
//		    eventMZcutVariablesVector_.push_back(mZcut_ptqq);
//		    eventMZcutVariablesVector_.push_back(mZcut_mqq);
//		    eventMZcutVariablesVector_.push_back(mZcut_etaqq);
//		    eventMZcutVariablesVector_.push_back(mZcut_pthh);
//		    eventMZcutVariablesVector_.push_back(mZcut_mhh);
//		    eventMZcutVariablesVector_.push_back(mZcut_etahh);
//		    eventMZcutVariablesVector_.push_back(mZcut_ptll);
//		    eventMZcutVariablesVector_.push_back(mZcut_mll);
//		    eventMZcutVariablesVector_.push_back(mZcut_etall);
//		    // 2-particle vars
//		    // ---------------
//		    double mZcut_dphiTjetZjet = DeltaPhi((Tjet->p4()).phi(),(Zjet->p4()).phi());
//		    double mZcut_dphiTjetZlep = DeltaPhi((Tjet->p4()).phi(),Zlep->phi());
//		    double mZcut_dphiminTZ = mZcut_dphiTjetZjet;
//		    if ( mZcut_dphiminTZ>mZcut_dphiTjetZlep ) mZcut_dphiminTZ = mZcut_dphiTjetZlep;
//		    double mZcut_detaTjetZjet = fabs((Tjet->p4()).eta()-(Zjet->p4()).eta());
//		    double mZcut_detaTjetZlep = fabs((Tjet->p4()).eta()-Zlep->eta());
//		    double mZcut_detaminTZ = mZcut_detaTjetZjet;
//		    if ( mZcut_detaminTZ>mZcut_detaTjetZlep ) mZcut_detaminTZ = mZcut_detaTjetZlep;
//		    double mZcut_dphiZjetZlep = DeltaPhi((Zjet->p4()).phi(),Zlep->phi());
//		    double mZcut_detaZjetZlep = fabs((Zjet->p4()).eta()-Zlep->eta());
//		    double mZcut_massTjetZjet = (Tjet->p4()+Zjet->p4()).mass();
//		    double mZcut_massTjetZlep = (Tjet->p4()+Zlep->p4()).mass();
//		    double mZcut_massZjetZlep = (Zjet->p4()+Zlep->p4()).mass();
//		    eventMZcutVariablesVector_.push_back(mZcut_dphiTjetZjet);
//		    eventMZcutVariablesVector_.push_back(mZcut_dphiTjetZlep);
//		    eventMZcutVariablesVector_.push_back(mZcut_dphiminTZ);
//		    eventMZcutVariablesVector_.push_back(mZcut_detaTjetZjet);
//		    eventMZcutVariablesVector_.push_back(mZcut_detaTjetZlep);
//		    eventMZcutVariablesVector_.push_back(mZcut_detaminTZ);
//		    eventMZcutVariablesVector_.push_back(mZcut_dphiZjetZlep);
//		    eventMZcutVariablesVector_.push_back(mZcut_detaZjetZlep);
//		    eventMZcutVariablesVector_.push_back(mZcut_massTjetZjet);
//		    eventMZcutVariablesVector_.push_back(mZcut_massTjetZlep);
//		    eventMZcutVariablesVector_.push_back(mZcut_massZjetZlep);
//		    // 3-particle vars
//		    // ---------------
//		    double mZcut_massTZZ = (Tjet->p4()+Zjet->p4()+Zlep->p4()).mass();
//		    double mZcut_etaTZZ = fabs((Tjet->p4()+Zjet->p4()+Zlep->p4()).eta());
//		    double mZcut_ptTZZ = (Tjet->p4()+Zjet->p4()+Zlep->p4()).pt();
//		    eventMZcutVariablesVector_.push_back(mZcut_massTZZ);
//		    eventMZcutVariablesVector_.push_back(mZcut_etaTZZ);
//		    eventMZcutVariablesVector_.push_back(mZcut_ptTZZ);
//		    
//		    // Store the event number
//		    eventMZcutVariablesVector_.push_back(eventcounter_);
//
//		    if( writeTMVA_ ) {
//		      // Fill the tree for the TMVA
//		      tmvaMZcutTreeWriterPtr_->fill(eventMZcutVariablesVector_);
//		    }
//		    /*
//		      else {
//		      // Fill the values in the array used by the reader
//		      vector<Float_t>::const_iterator varIt = eventSignalVariablesVector_.begin();
//		      int iVarBackground = 0;
//		      for( ; varIt != eventMZcutVariablesVector_.end(); ++varIt, ++iVarBackground ) {
//		      variables_[iVarBackground] = *varIt;
//		      }
//		      discriminantBackground.push_back( reader_->EvaluateMVA( "BDT method" ) );
//		      }
//		    */
//		    eventMZcutVariablesVector_.clear();
//		  }
//		}
//	      }
//
//	      // =============================================================================================
//		
//	    }
//	  } else std::cout << "WARNING!: less than 4 jets above et cut" << std::endl;
//	} else std::cout << "WARNING!: less than 2 jets above et cut" << std::endl;
//      } // end if partonsCandidates.size() >= njets_
//    } // end if ZlepVec.size() >= NZLEPTONS
    delete Hmother1;
    delete Hmother2;
    delete Zmother1;
    delete Zmother2;
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
TMVAntple::beginJob(const edm::EventSetup&)
{

  vecTypeIndexPair.push_back(std::pair<int,int>(pythiab_, 0));
  vecTypeIndexPair.push_back(std::pair<int,int>(pythiae_, 1));
  vecTypeIndexPair.push_back(std::pair<int,int>(pythiamu_,2));

  // File for output histograms
  // --------------------------
  // White background for the canvases
  gROOT->SetStyle("Plain");
  eventsNumber_ = new TH1D("eventsNumber","total number of events",1,0.,1.);

  lljjEventNumber_               = new TH1D("lljjEventNumber",              "",1,0.,1.);
  tmvaMatchingEventNumber_       = new TH1D("tmvaMatchingEventNumber",      "",1,0.,1.);
  combinatorialMatchEventNumber_ = new TH1D("combinatorialMatchEventNumber","",1,0.,1.);
  etaCutMatchingEventNumber_     = new TH1D("etaCutMatchingEventNumber",    "",1,0.,1.);
  mZMatchingEventNumber_         = new TH1D("mZMatchingEventNumber",        "",1,0.,1.);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
TMVAntple::endJob() {

  OutputFile->cd();

  // Fill histograms
  // ---------------
  eventsNumber_->SetBinContent(1,eventcounter_);
  eventsNumber_->Write();

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

}
//define this as a plug-in
DEFINE_FWK_MODULE(TMVAntple);

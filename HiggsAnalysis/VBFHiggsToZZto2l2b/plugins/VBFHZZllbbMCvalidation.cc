// -*- C++ -*-
//
// Package:    VBFHZZllbbMCvalidation
// Class:      VBFHZZllbbMCvalidation
// 
/**\class VBFHZZllbbMCvalidation VBFHZZllbbMCvalidation.cc HiggsAnalysis/VBFHiggsToZZto2l2bs/VQQmadgraphValidationAnalyzer/src/VBFHZZllbbMCvalidation.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Mia TOSI
//         Created:  Tue Jan 20 15:48:58 CET 2009
// $Id: VBFHZZllbbMCvalidation.cc,v 1.2 2009/04/28 16:45:01 tosi Exp $
//
//


// system include files
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbMCvalidation.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"

// constructors and destructor
//
VBFHZZllbbMCvalidation::VBFHZZllbbMCvalidation(const edm::ParameterSet& iConfig) :
  conf_( iConfig ),
  MCParticleLabel_(iConfig.getUntrackedParameter<std::string>("MCParticles") ),
  signal_(iConfig.getUntrackedParameter<bool>("signal") )
{
   //now do what ever initialization is needed
  eventcounter_ = 0;
  eventVBFcounter_ = 0;
  eventZcounter_  = 0;
  eventZZcounter_ = 0;
  eventHadronicZcounter_ = 0;
  eventLeptonicZcounter_ = 0;
  eventZintoEcounter_   = 0;
  eventZintoMUcounter_  = 0;
  eventZintoTAUcounter_ = 0;
  eventHeavyQcounter_ = 0;
  eventLightQcounter_ = 0;

  //  gROOT->Time();

}


VBFHZZllbbMCvalidation::~VBFHZZllbbMCvalidation()
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
VBFHZZllbbMCvalidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
  using namespace vbfhzz2l2b;
  using namespace std;

  eventcounter_++;
  if ( eventcounter_/100 == float(eventcounter_)/100. ) {
    std::cout << "Event number " << eventcounter_ << std::endl;
  }

  // MCParticle
  // ----------
  edm::Handle < GenParticleCollection > MCparticles;
  iEvent.getByLabel( MCParticleLabel_, MCparticles );

  unsigned int MCparticlesSize = MCparticles->size();
  
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
  std::vector< const Candidate * > Z1;
  std::vector< const Candidate * > Z2;
  const Candidate * Hmother1 = 0;
  const Candidate * Hmother2 = 0;
  const Candidate * Zmother1 = 0;
  const Candidate * Zmother2 = 0;

  GenParticleCollection::const_iterator MCparticle = MCparticles->begin();
  MCparticle = MCparticles->begin();
  for ( ; MCparticle != MCparticles->end(); ++MCparticle ) {    
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
      Zrapidity_ ->Fill(MCparticle->rapidity());
      Zeta_      ->Fill(MCparticle->eta()     );
      Zpt_       ->Fill(MCparticle->pt()      );
      Zphi_      ->Fill(MCparticle->phi()     );
      Zmass_     ->Fill(MCparticle->mass()    );


      int tmpNumberOfDaughters = MCparticle->numberOfDaughters();
      //      std::cout << "Z numberOfDaughters: " << tmpNumberOfDaughters << std::endl;
      for ( int daughterIndex = 0; daughterIndex != tmpNumberOfDaughters; daughterIndex++ ) {
	int tmpZdaughterId = MCparticle->daughter(daughterIndex)->pdgId();
	// hadronically decay Z
	if ( fabs(tmpZdaughterId) <= pythiat_ ) {
	  //	  std::cout << "found hadronically decaying Z into " << tmpZdaughterId << std::endl;
	  hadronicZparticles++;
	  Z1.push_back(&*(MCparticle->daughter(daughterIndex)));
	  if ( fabs(tmpZdaughterId) == pythiab_ ) ZintoBQUARKcounter++;
	  ZquarkFlavour_->Fill(tmpZdaughterId);
	}
	// leptonically decay Z
	else if ( fabs(tmpZdaughterId) >= pythiae_ && fabs(tmpZdaughterId) <= pythiatau_ ) {
	  //	  std::cout << "found leptonically decaying Z into " << tmpZdaughterId << std::endl;
	  leptonicZparticles++;
	  if ( fabs(tmpZdaughterId) == pythiae_   ) ZintoEcounter++; 
	  if ( fabs(tmpZdaughterId) == pythiamu_  ) ZintoMUcounter++;
	  if ( fabs(tmpZdaughterId) == pythiatau_ ) ZintoTAUcounter++;
	  Z2.push_back(&*(MCparticle->daughter(daughterIndex)));
	  ZleptonFlavour_->Fill(tmpZdaughterId);
	}
      }
    }
   } // end loop over MCparticle
   
  unsigned int heavyQcounter = 0; 
  unsigned int lightQcounter = 0;
  MCparticle = MCparticles->begin();
  for ( ; MCparticle != MCparticles->end(); ++MCparticle ) {    
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
      if (tmpParticleId == pythiagluon_) tagQuarkFlavour_->Fill(0.);
      else {
	tagQuarkFlavour_->Fill(tmpParticleId);
	if ( fabs(tmpParticleId) >= pythiac_ ) heavyQcounter++;
	else lightQcounter++;
      }
    }
    std::cout << "" << std::endl;
    std::cout << "***************************************" << std::endl;
  }

  std::cout << "heavyQcounter: " << heavyQcounter << std::endl;
  std::cout << "lightQcounter: " << lightQcounter << std::endl;

  if ( tagSystem.size() == 2 ) eventVBFcounter_++;
  if ( Zparticles >= 1 ) eventZcounter_++;
  if ( Zparticles >= 2 ) eventZZcounter_++;
  if ( Zparticles >= 1 && hadronicZparticles >= 1 ) eventHadronicZcounter_++;
  if ( Zparticles >= 1 && leptonicZparticles >= 1 ) eventLeptonicZcounter_++;
  if ( ZintoEcounter   >= 1 ) eventZintoEcounter_++;
  if ( ZintoMUcounter  >= 1 ) eventZintoMUcounter_++;
  if ( ZintoTAUcounter >= 1 ) eventZintoTAUcounter_++;
  if ( ZintoBQUARKcounter >= 1 ) eventZintoBQUARKcounter_++;
  if ( heavyQcounter >= 1 ) eventHeavyQcounter_++;
  if ( heavyQcounter == 2 ) eventHeavyQcounter_++;
  if ( lightQcounter >= 1 ) eventLightQcounter_++;
  if ( lightQcounter == 2 ) eventLightQcounter_++;

  std::cout << "eventVBFcounter: " << eventVBFcounter_ << std::endl;
  std::cout << "eventHeavyQcounter: " << eventHeavyQcounter_ << std::endl;
  std::cout << "eventLightQcounter: " << eventLightQcounter_ << std::endl;

  std::cout << "number of Z: " << Zparticles << std::endl;
  std::cout << "hadronically decaying: " << hadronicZparticles << " <--> leptonically decaying: " << leptonicZparticles << std::endl;
  std::cout << "Z1.size(): " << Z1.size() << " <--> Z2.size(): " << Z2.size() << std::endl;
  std::cout << "tagSystem.size(): " << tagSystem.size() << std::endl;

  if (Z2.size() == 2) {
    leptonicZrapidity_->Fill((Z2[0]->p4()+Z2[1]->p4()).Rapidity());
    leptonicZeta_     ->Fill((Z2[0]->p4()+Z2[1]->p4()).eta()     );
    leptonicZpt_      ->Fill((Z2[0]->p4()+Z2[1]->p4()).pt()      );
    leptonicZphi_     ->Fill((Z2[0]->p4()+Z2[1]->p4()).phi()     );
    leptonicZmass_    ->Fill((Z2[0]->p4()+Z2[1]->p4()).mass()    );
  }
  if (Z1.size() == 2) {
    hadronicZrapidity_->Fill((Z1[0]->p4()+Z1[1]->p4()).Rapidity());
    hadronicZeta_     ->Fill((Z1[0]->p4()+Z1[1]->p4()).eta()     );
    hadronicZpt_      ->Fill((Z1[0]->p4()+Z1[1]->p4()).pt()      );
    hadronicZphi_     ->Fill((Z1[0]->p4()+Z1[1]->p4()).phi()     );
    hadronicZmass_    ->Fill((Z1[0]->p4()+Z1[1]->p4()).mass()    );

    zQUARKrapidity_->Fill(Z1[0]->rapidity());
    zQUARKeta_     ->Fill(Z1[0]->eta()     );
    zQUARKpt_      ->Fill(Z1[0]->pt()      );
    zQUARKphi_     ->Fill(Z1[0]->phi()     );
    zQUARKrapidity_->Fill(Z1[1]->rapidity());
    zQUARKeta_     ->Fill(Z1[1]->eta()     );
    zQUARKpt_      ->Fill(Z1[1]->pt()      );
    zQUARKphi_     ->Fill(Z1[1]->phi()     );

    zQUARKetaVSphi_->Fill(Z1[0]->eta(),Z1[0]->phi());
    zQUARKetaVSphi_->Fill(Z1[1]->eta(),Z1[1]->phi());
    zQUARKetaVSeta_->Fill(Z1[0]->eta(),Z1[1]->eta());
  }
  if (tagSystem.size() == 2) {
    tagSYSTEMrapidity_->Fill((tagSystem[0]->p4()+tagSystem[1]->p4()).Rapidity());
    tagSYSTEMeta_     ->Fill((tagSystem[0]->p4()+tagSystem[1]->p4()).eta()     );
    tagSYSTEMpt_      ->Fill((tagSystem[0]->p4()+tagSystem[1]->p4()).pt()      );
    tagSYSTEMphi_     ->Fill((tagSystem[0]->p4()+tagSystem[1]->p4()).phi()     );
    tagSYSTEMmass_    ->Fill((tagSystem[0]->p4()+tagSystem[1]->p4()).mass()    );
    
    tagQUARKrapidity_->Fill(tagSystem[0]->rapidity());
    tagQUARKeta_     ->Fill(tagSystem[0]->eta()     );
    tagQUARKpt_      ->Fill(tagSystem[0]->pt()      );
    tagQUARKphi_     ->Fill(tagSystem[0]->phi()     );
    tagQUARKrapidity_->Fill(tagSystem[1]->rapidity());
    tagQUARKeta_     ->Fill(tagSystem[1]->eta()     );
    tagQUARKpt_      ->Fill(tagSystem[1]->pt()      );
    tagQUARKphi_     ->Fill(tagSystem[1]->phi()     );
    
    tagQUARKetaVSphi_->Fill(tagSystem[0]->eta(),tagSystem[0]->phi());
    tagQUARKetaVSphi_->Fill(tagSystem[1]->eta(),tagSystem[1]->phi());
    tagQUARKetaVSeta_->Fill(tagSystem[0]->eta(),tagSystem[1]->eta());
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
VBFHZZllbbMCvalidation::beginJob(const edm::EventSetup&)
{
  // File for output histograms
  // --------------------------
  OutputFile = new TFile((conf_.getUntrackedParameter<std::string>("OutputName")).c_str() ,
			 "RECREATE","VQQ-madgraphFall08_IDEAL_GEN-SIM-RECOOutput");
  // The file must be opened first, 
  // so that becomes the default position for all the histograms
  OutputFile->cd();
  // White background for the canvases
  //  gROOT->SetStyle("Plain");
  eventsNumber_ = new TH1D("eventsNumber","total number of events",1,0.,1.);

  tagQuarkFlavour_ = new TH1D("tagQuarkFlavour","number of different tag partons flavour",13,-6.,7.);
  tagQuarkFlavour_profile_ = new TProfile("tagQuarkFlavour_profile","number of different tag partons flavour",13,-6.,7.,0.,5000.);
  ZleptonFlavour_ = new TH1D("ZleptonFlavour","number of different leptons flavour",31,-15.,16.);
  ZleptonFlavour_profile_ = new TProfile("ZleptonFlavour_profile","number of different leptons flavour",31,-15.,16.,0.,5000.);
  ZquarkFlavour_ = new TH1D("ZquarkFlavour","number of different quark flavour",13,-6.,7.);
  ZquarkFlavour_profile_ = new TProfile("ZquarkFlavour_profile","number of different quark flavour",13,-6.,7.,0.,5000.);

  int nbin_ = 10;
  Zrapidity_ = new TH1D("Zrapidity","Z rapidity",nbin_,-10., 10.);
  Zeta_      = new TH1D("Zeta",     "Z #eta",    nbin_,-10., 10.);
  Zpt_       = new TH1D("Zpt",      "Z p_{T}",   nbin_,  0.,400.);
  Zphi_      = new TH1D("Zphi",     "Z #phi",    nbin_,-3.5, 3.5);
  Zmass_     = new TH1D("Zmass",    "Z mass",    nbin_,  0.,160.);
  leptonicZrapidity_ = new TH1D("leptonicZrapidity","leptonic Z rapidity",nbin_,-10., 10.);
  leptonicZeta_      = new TH1D("leptonicZeta",     "leptonic Z #eta",    nbin_,-10., 10.);
  leptonicZpt_       = new TH1D("leptonicZpt",      "leptonic Z p_{T}",   nbin_,  0.,400.);
  leptonicZphi_      = new TH1D("leptonicZphi",     "leptonic Z #phi",    nbin_,-3.5, 3.5);
  leptonicZmass_     = new TH1D("leptonicZmass",    "leptonic Z mass",    nbin_,  0.,160.);
  hadronicZrapidity_ = new TH1D("hadronicZrapidity","hadronic Z rapidity",nbin_,-10., 10.);
  hadronicZeta_      = new TH1D("hadronicZeta",     "hadronic Z #eta",    nbin_,-10., 10.);
  hadronicZpt_       = new TH1D("hadronicZpt",      "hadronic Z p_{T}",   nbin_,  0.,400.);
  hadronicZphi_      = new TH1D("hadronicZphi",     "hadronic Z #phi",    nbin_,-3.5, 3.5);
  hadronicZmass_     = new TH1D("hadronicZmass",    "hadronic Z mass",    nbin_,  0.,160.);

  zQUARKrapidity_  = new TH1D("zQUARKrapidity", "Z quark rapidity", nbin_,-10.,  10.);
  zQUARKeta_       = new TH1D("zQUARKeta",      "Z quark #eta",     nbin_,-10.,  10.);
  zQUARKpt_        = new TH1D("zQUARKpt",       "Z quark p_{T}",    nbin_,  0., 400.);
  zQUARKphi_       = new TH1D("zQUARKphi",      "Z quark #phi",     nbin_,-3.5,  3.5);

  tagSYSTEMrapidity_ = new TH1D("tagSYSTEMrapidity","tag system rapidity",nbin_,-10.,  10.);
  tagSYSTEMeta_      = new TH1D("tagSYSTEMeta",     "tag system #eta",    nbin_,-10.,  10.);
  tagSYSTEMpt_       = new TH1D("tagSYSTEMpt",      "tag system p_{T}",   nbin_,  0., 400.);
  tagSYSTEMphi_      = new TH1D("tagSYSTEMphi",     "tag system #phi",    nbin_,-3.5,  3.5);
  tagSYSTEMmass_     = new TH1D("tagSYSTEMmass",    "tag system mass",    nbin_,  0.,6000.);
  tagQUARKrapidity_  = new TH1D("tagQUARKrapidity", "tag quark rapidity", nbin_,-10.,  10.);
  tagQUARKeta_       = new TH1D("tagQUARKeta",      "tag quark #eta",     nbin_,-10.,  10.);
  tagQUARKpt_        = new TH1D("tagQUARKpt",       "tag quark p_{T}",    nbin_,  0., 400.);
  tagQUARKphi_       = new TH1D("tagQUARKphi",      "tag quark #phi",     nbin_,-3.5,  3.5);

  tagQUARKetaVSphi_ = new TH2D("tagQUARKetaVSphi","tag quark #eta-#phi plane",    nbin_,-10.,10.,nbin_,-3.5,3.5);
  tagQUARKetaVSeta_ = new TH2D("tagQUARKetaVSeta","tag quark #eta-#eta",          nbin_,-10.,10.,nbin_,-10.,10.);
  zQUARKetaVSphi_   = new TH2D("zQUARKetaVSphi",  "quark from Z #eta-#phi plane", nbin_,-10.,10.,nbin_,-3.5,3.5);
  zQUARKetaVSeta_   = new TH2D("zQUARKetaVSeta",  "quark from Z #eta-#eta",       nbin_,-10.,10.,nbin_,-10.,10.);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
VBFHZZllbbMCvalidation::endJob() {

  OutputFile->cd();

  std::cout << "number of VBF events: "             << eventVBFcounter_ << std::endl;
  std::cout << "number of events w/ at least 1 Z: " << eventZcounter_  << std::endl;
  std::cout << "number of events w/ 2 Z: "          << eventZZcounter_ << std::endl;
  std::cout << "number of leptonically Z decays: "  << eventLeptonicZcounter_ << std::endl;
  std::cout << "number of hadronically Z decays: "  << eventHadronicZcounter_ << std::endl;
  std::cout << "number of light quarks produced w/ V: " << eventLightQcounter_ << std::endl;
  std::cout << "number of heavy quarks produced w/ V: " << eventHeavyQcounter_ << std::endl;

  eventsNumber_->SetBinContent(1,eventcounter_);
  eventsNumber_->Write();

  tagQuarkFlavour_->Write();
  int nbin = tagQuarkFlavour_->GetNbinsX();
  for ( int ibin = 1; ibin <= nbin; ibin++ ) {
    double ibinEntries = tagQuarkFlavour_->GetBinContent(ibin);
    tagQuarkFlavour_profile_->Fill(double(ibin)-7.,ibinEntries);
  }
  tagQuarkFlavour_profile_->Write();

  ZleptonFlavour_->Write();
  nbin = ZleptonFlavour_->GetNbinsX();
  for ( int ibin = 1; ibin <= nbin; ibin++ ) {
    double ibinEntries = ZleptonFlavour_->GetBinContent(ibin);
    ZleptonFlavour_profile_->Fill(double(ibin)-16.,ibinEntries);
  }
  ZleptonFlavour_profile_->Write();
  ZquarkFlavour_->Write();
  nbin = ZquarkFlavour_->GetNbinsX();
  for ( int ibin = 1; ibin <= nbin; ibin++ ) {
    double ibinEntries = ZquarkFlavour_->GetBinContent(ibin);
    ZquarkFlavour_profile_->Fill(double(ibin)-7.,ibinEntries);
  }
  ZquarkFlavour_profile_->Write();


  Zrapidity_->Write();
  Zeta_     ->Write();
  Zpt_      ->Write();
  Zphi_     ->Write();
  Zmass_    ->Write();

  leptonicZrapidity_->Write();
  leptonicZeta_     ->Write();
  leptonicZpt_      ->Write();
  leptonicZphi_     ->Write();
  leptonicZmass_    ->Write();

  hadronicZrapidity_->Write();
  hadronicZeta_     ->Write();
  hadronicZpt_      ->Write();
  hadronicZphi_     ->Write();
  hadronicZmass_    ->Write();

  zQUARKrapidity_ ->Write();
  zQUARKeta_      ->Write();
  zQUARKpt_       ->Write();
  zQUARKphi_      ->Write();

  tagSYSTEMrapidity_->Write();
  tagSYSTEMeta_     ->Write();
  tagSYSTEMpt_      ->Write();
  tagSYSTEMphi_     ->Write();
  tagSYSTEMmass_    ->Write();
  tagQUARKrapidity_ ->Write();
  tagQUARKeta_      ->Write();
  tagQUARKpt_       ->Write();
  tagQUARKphi_      ->Write();

  tagQUARKetaVSphi_->Write();
  tagQUARKetaVSeta_->Write();
  zQUARKetaVSphi_  ->Write();
  zQUARKetaVSeta_  ->Write();


}

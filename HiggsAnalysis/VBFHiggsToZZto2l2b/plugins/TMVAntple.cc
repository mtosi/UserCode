//
// Original Author:  Mia Tosi
//         Created:  Fri Feb 22 17:56:22 CET 2008
// $Id: TMVAntple.cc,v 1.2 2009/04/28 16:45:01 tosi Exp $
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

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

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
  signal_          ( iConfig.getParameter<int>           ( "signal"           ) ), // 1:Signal,  0:Background
  whichSim_        ( iConfig.getParameter<int>           ( "whichSim"         ) ), // 0:FastSim, 1:FullSim
  electronLabel_   ( iConfig.getParameter<edm::InputTag> ( "electronLabel"    ) ),
  muonLabel_       ( iConfig.getParameter<edm::InputTag> ( "muonLabel"        ) ),
  metLabel_        ( iConfig.getParameter<edm::InputTag> ( "metLabel"         ) ),
  mcParticleLabel_ ( iConfig.getParameter<edm::InputTag> ( "mcParticleLabel"  ) ),
  leptonicZLabel_  ( iConfig.getParameter<std::string>   ( "leptonicZLabel"   ) ),
  hadronicZLabel_  ( iConfig.getParameter<std::string>   ( "hadronicZLabel"   ) ),
  tagSystemLabel_  ( iConfig.getParameter<std::string>   ( "tagSystemLabel"   ) ),
  tmvaSuffix_( iConfig.getParameter<std::string>   ( "tmvaSuffix" ) )
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
  tmvaTreeWriterPtr_.reset(new TMVAtreeWriter<Float_t>(eventVariablesNamesVector_,tmvaSuffix_));

  // constants, enums and typedefs
  eventcounter_             = 0;

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

  // reconstructed leptonicZ
  edm::Handle<reco::MuonCollection> leptonicZHandle;
  iEvent.getByLabel(leptonicZLabel_,"corJetWithBTagDiscr",leptonicZHandle);
  // reconstructed hadronicZ
  edm::Handle<reco::CaloJetCollection> hadronicZHandle;
  iEvent.getByLabel(hadronicZLabel_,"corJetWithBTagDiscr",hadronicZHandle);
  // reconstructed tagSystem
  edm::Handle<reco::CaloJetCollection> tagSystemHandle;
  iEvent.getByLabel(tagSystemLabel_,"corJetWithBTagDiscr",tagSystemHandle);

  unsigned muonCollectionSize       = muonHandle->size();
  unsigned electronCollectionSize   = electronHandle->size();
  unsigned mcParticleCollectionSize = mcParticleHandle->size();
  unsigned leptonicZCollectionSize  = leptonicZHandle->size();
  unsigned hadronicZCollectionSize  = hadronicZHandle->size();
  unsigned tagSystemCollectionSize  = tagSystemHandle->size();

  const reco::MuonCollection & lepZmuonCol = *(leptonicZHandle.product());
  const reco::CaloJetCollection & hadZjetCol = *(hadronicZHandle.product());
  const reco::CaloJetCollection & tagSjetCol = *(tagSystemHandle.product());

  if ( leptonicZCollectionSize >= NZLEPTONS && hadronicZCollectionSize >= NZJETS && tagSystemCollectionSize >= NTAGJETS ) {

    tmvaMatchingCounter_++;
    
    Particle * hadZ = new Particle(0,hadZjetCol[0].p4()+hadZjetCol[1].p4(),hadZjetCol[0].vertex(),    pythiaZ_,0,true);
    Particle * lepZ = new Particle(0,lepZmuonCol[0].p4()+lepZmuonCol[1].p4(),math::XYZPoint(0.,0.,0.),pythiaZ_,0,true);
    Particle * tagS = new Particle(0,tagSjetCol[0].p4()+tagSjetCol[1].p4(),tagSjetCol[0].vertex(),           0,0,true);
    
    
    // qq: tag jets
    // ------------
    double dphiqq = deltaPhi(tagSjetCol[0].phi(),tagSjetCol[1].phi());
    double detaqq = fabs(tagSjetCol[0].eta()-tagSjetCol[1].eta());
    double drqq   = deltaR(tagSjetCol[0].eta(),tagSjetCol[0].phi(),tagSjetCol[1].eta(),tagSjetCol[1].phi());
    double ptminqq = tagSjetCol[0].et();
    double ptmaxqq = tagSjetCol[1].et();
    if ( tagSjetCol[1].et() < ptminqq ) {
      ptminqq = ptmaxqq;
      ptmaxqq = tagSjetCol[0].et();
    }
    double etaminqq = fabs(tagSjetCol[0].eta());
    double etamaxqq = fabs(tagSjetCol[1].eta());
    if ( fabs(tagSjetCol[1].eta()) < etaminqq ) {
      etaminqq = etamaxqq;
      etamaxqq = fabs(tagSjetCol[0].eta());
    }
    eventVariablesVector_.push_back(tagSjetCol[0].et());
    eventVariablesVector_.push_back(tagSjetCol[1].et());
    eventVariablesVector_.push_back(tagSjetCol[0].eta());
    eventVariablesVector_.push_back(tagSjetCol[1].eta());
    eventVariablesVector_.push_back(tagSjetCol[0].phi());
    eventVariablesVector_.push_back(tagSjetCol[1].phi());
    eventVariablesVector_.push_back(dphiqq);
    eventVariablesVector_.push_back(detaqq);
    eventVariablesVector_.push_back(ptminqq);
    eventVariablesVector_.push_back(ptmaxqq);
    eventVariablesVector_.push_back(etaminqq);
    eventVariablesVector_.push_back(etamaxqq);
    
    // hh: Z jets
    // ----------
    double dphihh = deltaPhi(hadZjetCol[0].phi(),hadZjetCol[1].phi());
    double detahh = fabs(hadZjetCol[0].eta()-hadZjetCol[1].eta());
    double drhh   = deltaR(hadZjetCol[0].eta(),hadZjetCol[0].phi(),hadZjetCol[1].eta(),hadZjetCol[1].phi());
    double ptminhh = hadZjetCol[0].et();
    double ptmaxhh = hadZjetCol[1].et();
    if ( hadZjetCol[1].et() < ptminhh ) {
      ptminhh = ptmaxhh;
      ptmaxhh = hadZjetCol[0].et();
    }
    double etaminhh = fabs(hadZjetCol[0].eta());
    double etamaxhh = fabs(hadZjetCol[1].eta());
    if ( fabs(hadZjetCol[1].eta()) < etaminhh ) {
      etaminhh = etamaxhh;
      etamaxhh = fabs(hadZjetCol[0].eta());
    }
    eventVariablesVector_.push_back(hadZjetCol[0].et());
    eventVariablesVector_.push_back(hadZjetCol[1].et());
    eventVariablesVector_.push_back(hadZjetCol[0].eta());
    eventVariablesVector_.push_back(hadZjetCol[1].eta());
    eventVariablesVector_.push_back(hadZjetCol[0].phi());
    eventVariablesVector_.push_back(hadZjetCol[1].phi());
    eventVariablesVector_.push_back(dphihh);
    eventVariablesVector_.push_back(detahh);
    eventVariablesVector_.push_back(ptminhh);
    eventVariablesVector_.push_back(ptmaxhh);
    eventVariablesVector_.push_back(etaminhh);
    eventVariablesVector_.push_back(etamaxhh);
    
    // qq system
    // ---------
    double ptqq  = tagS->pt();
    double mqq   = tagS->mass();
    double etaqq = fabs(tagS->eta());
    eventVariablesVector_.push_back(ptqq);
    eventVariablesVector_.push_back(mqq);
    eventVariablesVector_.push_back(etaqq);

    // hh system
    // ---------
    double pthh  = hadZ->pt();
    double mhh   = hadZ->mass();
    double etahh = fabs(hadZ->eta());
    eventVariablesVector_.push_back(pthh);
    eventVariablesVector_.push_back(mhh);
    eventVariablesVector_.push_back(etahh);

    // ll system
    // ---------
    double ptll  = lepZ->pt();
    double mll   = lepZ->mass();
    double etall = fabs(lepZ->eta());
    eventVariablesVector_.push_back(ptll);
    eventVariablesVector_.push_back(mll);
    eventVariablesVector_.push_back(etall);

    std::cout << "Systems Pt:   " << ptqq  << " " << pthh  << " " << ptll  << std::endl;
    std::cout << "Systems eta:  " << etaqq << " " << etahh << " " << etall << std::endl;
    std::cout << "Systems mass: " << mqq   << " " << mhh   << " " << mll   << std::endl;
    
    // 2-particle vars
    // ---------------
    double dphiTZjet = deltaPhi(tagS->phi(),hadZ->phi());
    double dphiTZlep = deltaPhi(tagS->phi(),lepZ->phi());
    double dphiminTZ = dphiTZjet;
    if ( dphiminTZ > dphiTZlep ) dphiminTZ = dphiTZlep;
    double detaTZjet = fabs(tagS->eta()-hadZ->eta());
    double detaTZlep = fabs(tagS->eta()-lepZ->eta());
    double detaminTZ = detaTZjet;
    if ( detaminTZ > detaTZlep ) detaminTZ = detaTZlep;
    double dphiZjetZlep = deltaPhi(hadZ->phi(),lepZ->phi());
    double detaZjetZlep = fabs(hadZ->eta()-lepZ->eta());
    double massTZjet    = (tagS->p4()+(hadZ->p4())).mass();
    double massTZlep    = (tagS->p4()+(lepZ->p4())).mass();
    double massZjetZlep = (hadZ->p4()+(lepZ->p4())).mass();
    eventVariablesVector_.push_back(dphiTZjet);
    eventVariablesVector_.push_back(dphiTZlep);
    eventVariablesVector_.push_back(dphiminTZ);
    eventVariablesVector_.push_back(detaTZjet);
    eventVariablesVector_.push_back(detaTZlep);
    eventVariablesVector_.push_back(detaminTZ);
    eventVariablesVector_.push_back(dphiZjetZlep);
    eventVariablesVector_.push_back(detaZjetZlep);
    eventVariablesVector_.push_back(massTZjet);
    eventVariablesVector_.push_back(massTZlep);
    eventVariablesVector_.push_back(massZjetZlep);
    
    // 3-particle vars
    // ---------------
    double massTZZ = (tagS->p4()+(hadZ->p4())+(lepZ->p4())).mass();
    double ptTZZ   = (tagS->p4()+(hadZ->p4())+(lepZ->p4())).pt();
    double etaTZZ  = fabs((tagS->p4()+(hadZ->p4())+(lepZ->p4())).eta());
    eventVariablesVector_.push_back(massTZZ);
    eventVariablesVector_.push_back(etaTZZ);
    eventVariablesVector_.push_back(ptTZZ);
    
    // Store the event number
    eventVariablesVector_.push_back(eventcounter_);
    
    // Fill the tree for the TMVA
    tmvaTreeWriterPtr_->fill(eventVariablesVector_);

    eventVariablesVector_.clear();

    delete tagS;
    delete hadZ;
    delete lepZ;
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
TMVAntple::beginJob(const edm::EventSetup&)
{

  // File for output histograms
  // --------------------------
  // White background for the canvases
  gROOT->SetStyle("Plain");
  eventsNumber_ = new TH1D("eventsNumber","total number of events",1,0.,1.);

  tmvaMatchingEventNumber_       = new TH1D("tmvaMatchingEventNumber",      "",1,0.,1.);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TMVAntple::endJob() {

  OutputFile->cd();

  // Fill histograms
  // ---------------
  eventsNumber_->SetBinContent(1,eventcounter_);
  eventsNumber_->Write();

  tmvaMatchingEventNumber_->SetBinContent(1,double(tmvaMatchingCounter_));

  tmvaMatchingEventNumber_->Write();

  std::cout << "tmvaMatchingCounter: " << tmvaMatchingCounter_   
	    << " => " << double(tmvaMatchingCounter_  )/double(eventcounter_ )*100. << "%" 
	    << std::endl;

}
//define this as a plug-in
DEFINE_FWK_MODULE(TMVAntple);

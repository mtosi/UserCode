// system include files
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbMCfilterValidation.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/RefToBase.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/JetReco/interface/GenJet.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/VBFHZZllbbUtils.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/PythiaParticleIndex.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/ParticlesMass.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/ParticlesCharge.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/ProcessIndex.h"


enum { FASTSIM = 0,
       FULLSIM = 1
};

typedef reco::GenJetCollection::const_iterator genJetItr;

// constructors and destructor
//
VBFHZZllbbMCfilterValidation::VBFHZZllbbMCfilterValidation(const edm::ParameterSet& iConfig) :
  whichSim_      ( iConfig.getParameter<int>           ( "whichSim"      ) ), // 0:FastSim, 1:FullSim
  signal_        ( iConfig.getParameter<bool>          ( "signal"        ) ),
  genJetLabel_   ( iConfig.getParameter<edm::InputTag> ( "genJetLabel"   ) ),
  muonLabel_     ( iConfig.getParameter<edm::InputTag> ( "muonLabel"     ) ),
  electronLabel_ ( iConfig.getParameter<edm::InputTag> ( "electronLabel" ) ),
  jetNumberCut_   ( iConfig.getParameter<int>    ( "jetNumberCut"   ) ),
  firstJetPtCut_  ( iConfig.getParameter<double> ( "firstJetPtCut"  ) ),
  secondJetPtCut_ ( iConfig.getParameter<double> ( "secondJetPtCut" ) ),
  invMassCut_	  ( iConfig.getParameter<double> ( "invMassCut"     ) ),
  deltaEtaCut_	  ( iConfig.getParameter<double> ( "deltaEtaCut"    ) ),
  leptonPtCut_    ( iConfig.getParameter<double> ( "leptonPtCut"    ) )


{
   //now do what ever initialization is needed
  eventcounter_ = 0;
  eventVBFcounter_ = 0;
  eventGGFcounter_ = 0;

  //  gROOT->Time();

}


VBFHZZllbbMCfilterValidation::~VBFHZZllbbMCfilterValidation()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
VBFHZZllbbMCfilterValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
  using namespace vbfhzz2l2b;
  using namespace std;

  eventcounter_++;
  if ( eventcounter_/100 == float(eventcounter_)/100. ) {
    std::cout << "Event number " << eventcounter_ << std::endl;
  }

  int IDevent = getEventID(iEvent, whichSim_);
  if ( IDevent == HWWFusion_ || IDevent == HZZFusion_ ) eventVBFcounter_++;
  else if ( IDevent == HggFusion_ ) eventGGFcounter_++;
  

  edm::Handle< reco::GenJetCollection > genJetsHandle ;
  iEvent.getByLabel( genJetLabel_, genJetsHandle ) ;
  std::cout << "genJetsHandle->size(): " << genJetsHandle->size() << std::endl;  

  // pt sorting
  GenJetCollection genJetSortedColl=*genJetsHandle;
  std::sort(genJetSortedColl.begin(),genJetSortedColl.end(),PtGreater());

  // looking for the highest invariant mass jets pair
  std::pair<genJetItr,genJetItr> maxInvMassPair = 
    vbfhzz2l2b::findPair_maxInvMass_ptMinCut<genJetItr>(genJetsHandle->begin(), genJetsHandle->end(),
							firstJetPtCut_, secondJetPtCut_);
  double maxInvMass = ( (maxInvMassPair.first)->p4() + ((maxInvMassPair.second)->p4()) ).M();
  std::cout << "[VBFHZZllbbMCfilterValidation::analyze] max inv mass: " << maxInvMass << std::endl;

  // looking for the highest delta eta jets pair
  std::pair<genJetItr,genJetItr> maxDeltaEtaPair = 
    vbfhzz2l2b::findPair_maxDeltaEta_ptMinCut<genJetItr>(genJetsHandle->begin(), genJetsHandle->end(),
							 firstJetPtCut_, secondJetPtCut_);
  double maxDeltaEta = fabs( (maxDeltaEtaPair.first)->eta() - (maxDeltaEtaPair.second)->eta() );
  std::cout << "[VBFHZZllbbMCfilterValidation::analyze] max delta eta: " << maxDeltaEta << std::endl;

  
  // fill proper histograms
  if ( IDevent == HWWFusion_ || IDevent == HZZFusion_ ) {
    VBFfirstJetPt_        -> Fill(genJetSortedColl[0].pt());    
    VBFsecondJetPt_       -> Fill(genJetSortedColl[1].pt());
    VBFthirdJetPt_        -> Fill(genJetSortedColl[2].pt());
    VBFfourthJetPt_       -> Fill(genJetSortedColl[3].pt());
    VBFmaxInvMassJetJet_  -> Fill(maxInvMass);
    VBFmaxDeltaEtaJetJet_ -> Fill(maxDeltaEta);
  } else if ( IDevent == HggFusion_ ) {
    ggFfirstJetPt_        -> Fill(genJetSortedColl[0].pt());    
    ggFsecondJetPt_       -> Fill(genJetSortedColl[1].pt());
    ggFthirdJetPt_        -> Fill(genJetSortedColl[2].pt());
    ggFfourthJetPt_       -> Fill(genJetSortedColl[3].pt());
    ggFmaxInvMassJetJet_  -> Fill(maxInvMass);
    ggFmaxDeltaEtaJetJet_ -> Fill(maxDeltaEta);
  }
 
}



// ------------ method called once each job just before starting event loop  ------------
void 
VBFHZZllbbMCfilterValidation::beginJob(const edm::EventSetup&)
{
  // file for output histograms
  edm::Service<TFileService> fs;

  eventsNumber_    = fs->make<TH1D>("eventsNumber",   "total number of events",    1,0.,1.);
  eventsVBFNumber_ = fs->make<TH1D>("eventsVBFNumber","total number of VBF events",1,0.,1.);
  eventsGGFNumber_ = fs->make<TH1D>("eventsGGFNumber","total number of ggF events",1,0.,1.);

  int nbin = 100;

  TFileDirectory VBFSubDir = fs->mkdir( "VBF" );
  VBFfirstJetPt_        = VBFSubDir.make<TH1D>("VBFfirstJetPt",       "1^{st} gen-jet p_{T} distribution",nbin,0., 500.); 
  VBFsecondJetPt_       = VBFSubDir.make<TH1D>("VBFsecondJetPt",      "2^{nd} gen-jet p_{T} distribution",nbin,0., 500.); 
  VBFthirdJetPt_        = VBFSubDir.make<TH1D>("VBFthirdJetPt",       "3^{rd} gen-jet p_{T} distribution",nbin,0., 500.); 
  VBFfourthJetPt_       = VBFSubDir.make<TH1D>("VBFfourthJetPt",      "4^{th} gen-jet p_{T} distribution",nbin,0., 500.); 
  VBFmaxDeltaEtaJetJet_ = VBFSubDir.make<TH1D>("VBFmaxDeltaEtaJetJet","maximum #Delta#eta between 2 jets",nbin,0.,  10.);
  VBFmaxInvMassJetJet_  = VBFSubDir.make<TH1D>("VBFmaxInvMassJetJet", "maximum 2 jets invariant mass",    nbin,0.,4000.);

  TFileDirectory ggFSubDir = fs->mkdir( "ggF" );
  ggFfirstJetPt_        = ggFSubDir.make<TH1D>("ggFfirstJetPt",       "1^{st} gen-jet p_{T} distribution",nbin,0., 500.); 
  ggFsecondJetPt_       = ggFSubDir.make<TH1D>("ggFsecondJetPt",      "2^{nd} gen-jet p_{T} distribution",nbin,0., 500.); 
  ggFthirdJetPt_        = ggFSubDir.make<TH1D>("ggFthirdJetPt",       "3^{rd} gen-jet p_{T} distribution",nbin,0., 500.); 
  ggFfourthJetPt_       = ggFSubDir.make<TH1D>("ggFfourthJetPt",      "4^{th} gen-jet p_{T} distribution",nbin,0., 500.); 
  ggFmaxDeltaEtaJetJet_ = ggFSubDir.make<TH1D>("ggFmaxDeltaEtaJetJet","maximum #Delta#eta between 2 jets",nbin,0.,  10.);
  ggFmaxInvMassJetJet_  = ggFSubDir.make<TH1D>("ggFmaxInvMassJetJet", "maximum 2 jets invariant mass",    nbin,0.,4000.);


  firstElectronPt_  = fs->make<TH1D>("firstElectronPt", "1^{st} electron p_{T} distribution",nbin,0.,200.);
  secondElectronPt_ = fs->make<TH1D>("secondElectronPt","2^{nd} electron p_{T} distribution",nbin,0.,200.);
  firstMuonPt_      = fs->make<TH1D>("firstMuonPt",     "1^{st} muon p_{T} distribution",    nbin,0.,200.);
  secondMuonPt_     = fs->make<TH1D>("secondMuonPt",    "2^{nd} muon p_{T} distribution",    nbin,0.,200.);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
VBFHZZllbbMCfilterValidation::endJob() {

  std::cout << "number of events:     " << eventcounter_    << std::endl;
  std::cout << "number of VBF events: " << eventVBFcounter_ << std::endl;
  std::cout << "number of ggF events: " << eventGGFcounter_ << std::endl;

  eventsNumber_    -> SetBinContent(1,eventcounter_   );
  eventsVBFNumber_ -> SetBinContent(1,eventVBFcounter_);
  eventsGGFNumber_ -> SetBinContent(1,eventGGFcounter_);

}

int
VBFHZZllbbMCfilterValidation::getEventID (const edm::Event & iEvent, const int & whichSim_ ) {

  int eventID = -1;
  if ( whichSim_ == FULLSIM ) {
    edm::Handle<edm::HepMCProduct> evtMC;
    iEvent.getByLabel("source", evtMC);
    
    const HepMC::GenEvent * mcEv = evtMC->GetEvent();
    eventID = mcEv->signal_process_id();
  }
  else if ( whichSim_ == FASTSIM ) {
    edm::Handle<int> genProcessID;
    iEvent.getByLabel( "genEventProcID", genProcessID );
    
    eventID = *genProcessID;
  }
  else
    std::cout << "--> WARNING: simulation not specificied!!" << std::endl;
  return eventID;
}    



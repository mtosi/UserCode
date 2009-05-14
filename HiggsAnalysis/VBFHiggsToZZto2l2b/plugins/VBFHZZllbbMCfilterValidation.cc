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
  whichSim_ ( iConfig.getParameter<int>  ( "whichSim" ) ), // 0:FastSim, 1:FullSim
  signal_   ( iConfig.getParameter<bool> ( "signal"   ) ),
  genJetLabel_      ( iConfig.getParameter<edm::InputTag> ( "genJetLabel"      ) ),
  genParticleLabel_ ( iConfig.getParameter<edm::InputTag> ( "genParticleLabel" ) ),
  jetNumberCut_       ( iConfig.getParameter<int>    ( "jetNumberCut"       ) ),
  firstJetPtCut_      ( iConfig.getParameter<double> ( "firstJetPtCut"      ) ),
  secondJetPtCut_     ( iConfig.getParameter<double> ( "secondJetPtCut"     ) ),
  jetpairInvMassCut_  ( iConfig.getParameter<double> ( "jetpairInvMassCut"  ) ),
  jetpairDeltaEtaCut_ ( iConfig.getParameter<double> ( "jetpairDeltaEtaCut" ) ),
  leptonEtaCut_   ( iConfig.getParameter<double> ( "leptonEtaCut" ) ),
  leptonPtCut_    ( iConfig.getParameter<double> ( "leptonPtCut"  ) )


{
   //now do what ever initialization is needed
  eventcounter_ = 0;
  eventVBFcounter_ = 0;
  eventGGFcounter_ = 0;

  jetNumberCut_eventcounter_   = 0;
  firstJetPtCut_eventcounter_  = 0;
  secondJetPtCut_eventcounter_ = 0;
  jetpairInvMassCut_eventcounter_     = 0;
  jetpairDeltaEtaCut_eventcounter_    = 0;
  leptonEtaCut_eventcounter_   = 0;
  leptonPtCut_eventcounter_    = 0;
  jetNumberCut_eventVBFcounter_   = 0;
  firstJetPtCut_eventVBFcounter_  = 0;
  secondJetPtCut_eventVBFcounter_ = 0;
  jetpairInvMassCut_eventVBFcounter_     = 0;
  jetpairDeltaEtaCut_eventVBFcounter_    = 0;
  leptonEtaCut_eventVBFcounter_   = 0;
  leptonPtCut_eventVBFcounter_    = 0;
  jetNumberCut_eventGGFcounter_   = 0;
  firstJetPtCut_eventGGFcounter_  = 0;
  secondJetPtCut_eventGGFcounter_ = 0;
  jetpairInvMassCut_eventGGFcounter_     = 0;
  jetpairDeltaEtaCut_eventGGFcounter_    = 0;
  leptonEtaCut_eventGGFcounter_   = 0;
  leptonPtCut_eventGGFcounter_    = 0;
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

  // pt sorting
  GenJetCollection genJetSortedColl=*genJetsHandle;
  std::sort(genJetSortedColl.begin(),genJetSortedColl.end(),PtGreater());

  // looking for the highest invariant mass jets pair
  std::pair<genJetItr,genJetItr> maxInvMassPair = 
    vbfhzz2l2b::findPair_maxInvMass_ptMinCut<genJetItr>(genJetsHandle->begin(), genJetsHandle->end(),
							firstJetPtCut_, secondJetPtCut_);
  double maxInvMass = 0.;
  if (maxInvMassPair.first != maxInvMassPair.second) 
    maxInvMass = ( (maxInvMassPair.first)->p4() + ((maxInvMassPair.second)->p4()) ).M();

  // looking for the highest delta eta jets pair
  std::pair<genJetItr,genJetItr> maxDeltaEtaPair = 
    vbfhzz2l2b::findPair_maxDeltaEta_ptMinCut<genJetItr>(genJetsHandle->begin(), genJetsHandle->end(),
							 firstJetPtCut_, secondJetPtCut_);
  double maxDeltaEta = 0.;
  if(maxDeltaEtaPair.first != maxDeltaEtaPair.second) 
    maxDeltaEta = fabs( (maxDeltaEtaPair.first)->eta() - (maxDeltaEtaPair.second)->eta() );

  edm::Handle < reco::GenParticleCollection > genParticlesHandle;
  iEvent.getByLabel( genParticleLabel_, genParticlesHandle ) ;

  std::vector<reco::GenParticle> electronsVec;
  std::vector<reco::GenParticle> muonsVec;
  std::vector<reco::GenParticle> leptonsVec;

  for ( reco::GenParticleCollection::const_iterator particle_itr = genParticlesHandle->begin();
	particle_itr != genParticlesHandle->end(); ++particle_itr ) {
    if ( fabs(particle_itr->pdgId()) == pythiae_  ) {
      electronsVec.push_back (*particle_itr);
      leptonsVec.push_back   (*particle_itr);
    }
    if ( fabs(particle_itr->pdgId()) == pythiamu_ ) {
      muonsVec.push_back   (*particle_itr);
      leptonsVec.push_back (*particle_itr);
    }
  }

  // pt sorting
  std::sort(electronsVec.begin(),electronsVec.end(),PtGreater());
  std::sort(muonsVec.begin(),muonsVec.end(),PtGreater());
  std::sort(leptonsVec.begin(),leptonsVec.end(),PtGreater());

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
  if (electronsVec.size() >= 1) firstElectronPt_  -> Fill(electronsVec[0].pt()); 
  if (electronsVec.size() >= 2) secondElectronPt_ -> Fill(electronsVec[1].pt()); 
  if (muonsVec.size() >= 1)     firstMuonPt_      -> Fill(muonsVec[0].pt()); 
  if (muonsVec.size() >= 2)     secondMuonPt_     -> Fill(muonsVec[1].pt()); 



  // filter efficiency
  if ( genJetSortedColl[0].pt() >= firstJetPtCut_  ) firstJetPtCut_eventcounter_++;
  if ( genJetSortedColl[1].pt() >= secondJetPtCut_ ) secondJetPtCut_eventcounter_++;
  if ( maxInvMass >= jetpairInvMassCut_            ) jetpairInvMassCut_eventcounter_++;
  if ( maxDeltaEta >= jetpairDeltaEtaCut_          ) jetpairDeltaEtaCut_eventcounter_++;
  /*
  if ( electronsVec.size() >= 1 ) 
    if ( electronsVec[0].pt() >= leptonPtCut_      ) leptonPtCut_eventcounter_++;
  if ( muonsVec.size() >= 1 ) 
    if ( muonsVec[0].pt() >= leptonPtCut_          ) leptonPtCut_eventcounter_++;
  */
  if ( leptonsVec.size() >= 1 )
    if ( leptonsVec[0].pt() >= leptonPtCut_        ) leptonPtCut_eventcounter_++;
  if ( IDevent == HWWFusion_ || IDevent == HZZFusion_ ) {
    if ( genJetSortedColl[0].pt() >= firstJetPtCut_  ) firstJetPtCut_eventVBFcounter_++;
    if ( genJetSortedColl[1].pt() >= secondJetPtCut_ ) secondJetPtCut_eventVBFcounter_++;
    if ( maxInvMass >= jetpairInvMassCut_            ) jetpairInvMassCut_eventVBFcounter_++;
    if ( maxDeltaEta >= jetpairDeltaEtaCut_          ) jetpairDeltaEtaCut_eventVBFcounter_++;
    /*
    if ( electronsVec.size() >= 1 ) 
      if ( electronsVec[0].pt() >= leptonPtCut_      ) leptonPtCut_eventVBFcounter_++;
    if ( muonsVec.size() >= 1 )
      if ( muonsVec[0].pt() >= leptonPtCut_          ) leptonPtCut_eventVBFcounter_++;
    */
    if ( leptonsVec.size() >= 1 )
      if ( leptonsVec[0].pt() >= leptonPtCut_        ) leptonPtCut_eventVBFcounter_++;
  } else if ( IDevent == HggFusion_ ) {
    if ( genJetSortedColl[0].pt() >= firstJetPtCut_  ) firstJetPtCut_eventGGFcounter_++;
    if ( genJetSortedColl[1].pt() >= secondJetPtCut_ ) secondJetPtCut_eventGGFcounter_++;
    if ( maxInvMass >= jetpairInvMassCut_            ) jetpairInvMassCut_eventGGFcounter_++;
    if ( maxDeltaEta >= jetpairDeltaEtaCut_          ) jetpairDeltaEtaCut_eventGGFcounter_++;
    /*
    if ( electronsVec.size() >= 1 ) 
      if ( electronsVec[0].pt() >= leptonPtCut_      ) leptonPtCut_eventGGFcounter_++;
    if ( muonsVec.size() >= 1 ) 
      if ( muonsVec[0].pt() >= leptonPtCut_          ) leptonPtCut_eventGGFcounter_++;
    */
    if ( leptonsVec.size() >= 1 )
      if ( leptonsVec[0].pt() >= leptonPtCut_        ) leptonPtCut_eventGGFcounter_++;

  }

}



// ------------ method called once each job just before starting event loop  ------------
void 
VBFHZZllbbMCfilterValidation::beginJob(const edm::EventSetup&)
{
  // file for output histograms
  edm::Service<TFileService> fs;

  eventsNumber_   = fs->make<TH1D>("eventsNumber",   "total number of events",    1,0.,1.);
  eventsIDNumber_ = fs->make<TH1D>("eventsVBFNumber","total number of VBF events",2,0.,2.);

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

  eventsNumber_   -> SetBinContent(1,eventcounter_   );
  eventsIDNumber_ -> SetBinContent(1,eventVBFcounter_);
  eventsIDNumber_ -> SetBinContent(2,eventGGFcounter_);

  std::cout << "************************************************************" << std::endl;
  std::cout << "number of events read: " << eventcounter_    << std::endl;
  std::cout << "number of VBF events:  " << eventVBFcounter_ << " (" << double(eventVBFcounter_)/double(eventcounter_)*100. << "%)" << std::endl;
  std::cout << "number of ggF events:  " << eventGGFcounter_ << " (" << double(eventGGFcounter_)/double(eventcounter_)*100. << "%)" << std::endl;
  std::cout << "*** Efficiency for the various subsamples *** " <<  std::endl;
  std::cout << "---> jet number cut" << std::endl;
  std::cout << "     w.r.t. all events: " << double(jetNumberCut_eventcounter_)/double(eventcounter_)*100. << "%" << std::endl; 
  std::cout << "     w.r.t. VBF events: " << double(jetNumberCut_eventVBFcounter_)/double(eventVBFcounter_)*100. << "%" << std::endl; 
  std::cout << "     w.r.t. ggF events: " << double(jetNumberCut_eventGGFcounter_)/double(eventGGFcounter_)*100. << "%" << std::endl; 
  std::cout << "---> first jet pt cut" << std::endl;
  std::cout << "     w.r.t. all events: " << double(firstJetPtCut_eventcounter_)/double(eventcounter_)*100. << "%" << std::endl; 
  std::cout << "     w.r.t. VBF events: " << double(firstJetPtCut_eventVBFcounter_)/double(eventVBFcounter_)*100. << "%" << std::endl; 
  std::cout << "     w.r.t. ggF events: " << double(firstJetPtCut_eventGGFcounter_)/double(eventGGFcounter_)*100. << "%" << std::endl; 
  std::cout << "---> second jet pt cut" << std::endl;
  std::cout << "     w.r.t. all events: " << double(secondJetPtCut_eventcounter_)/double(eventcounter_)*100. << "%" << std::endl; 
  std::cout << "     w.r.t. VBF events: " << double(secondJetPtCut_eventVBFcounter_)/double(eventVBFcounter_)*100. << "%" << std::endl; 
  std::cout << "     w.r.t. ggF events: " << double(secondJetPtCut_eventGGFcounter_)/double(eventGGFcounter_)*100. << "%" << std::endl; 
  std::cout << "---> invariant mass cut" << std::endl;
  std::cout << "     w.r.t. all events: " << double(jetpairInvMassCut_eventcounter_)/double(eventcounter_)*100. << "%" << std::endl; 
  std::cout << "     w.r.t. VBF events: " << double(jetpairInvMassCut_eventVBFcounter_)/double(eventVBFcounter_)*100. << "%" << std::endl; 
  std::cout << "     w.r.t. ggF events: " << double(jetpairInvMassCut_eventGGFcounter_)/double(eventGGFcounter_)*100. << "%" << std::endl; 
  std::cout << "---> delta eta cut" << std::endl;
  std::cout << "     w.r.t. all events: " << double(jetpairDeltaEtaCut_eventcounter_)/double(eventcounter_)*100. << "%" << std::endl; 
  std::cout << "     w.r.t. VBF events: " << double(jetpairDeltaEtaCut_eventVBFcounter_)/double(eventVBFcounter_)*100. << "%" << std::endl; 
  std::cout << "     w.r.t. ggF events: " << double(jetpairDeltaEtaCut_eventGGFcounter_)/double(eventGGFcounter_)*100. << "%" << std::endl; 
  std::cout << "---> lepton pt cut" << std::endl;
  std::cout << "     w.r.t. all events: " << double(leptonPtCut_eventcounter_)/double(eventcounter_)*100. << "%" << std::endl; 
  std::cout << "     w.r.t. VBF events: " << double(leptonPtCut_eventVBFcounter_)/double(eventVBFcounter_)*100. << "%" << std::endl; 
  std::cout << "     w.r.t. ggF events: " << double(leptonPtCut_eventGGFcounter_)/double(eventGGFcounter_)*100. << "%" << std::endl; 
  std::cout << "************************************************************" << std::endl;
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



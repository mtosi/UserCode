#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbDisplay.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"


#include <TMarker.h>


VBFHZZllbbDisplay::VBFHZZllbbDisplay (const edm::ParameterSet& iConfig) :
  genParticleLabel_( iConfig.getUntrackedParameter<edm::InputTag>( "genParticleLabel" ) ),
  jetLabel_        ( iConfig.getUntrackedParameter<edm::InputTag>( "jetLabel"         ) ),
  muonLabel_       ( iConfig.getUntrackedParameter<edm::InputTag>( "muonLabel"        ) ),
  electronLabel_   ( iConfig.getUntrackedParameter<edm::InputTag>( "electronLabel"    ) ),
  tagJetLabel_     ( iConfig.getUntrackedParameter<edm::InputTag>( "tagJetLabel"      ) ),
  metLabel_        ( iConfig.getUntrackedParameter<edm::InputTag>( "metLabel"         ) )
{}


// --------------------------------------------------------------------


VBFHZZllbbDisplay::~VBFHZZllbbDisplay ()
{
 
   delete canvas_ ;
   delete histo2D_ ;
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


// --------------------------------------------------------------------


void
VBFHZZllbbDisplay::analyze (const edm::Event& iEvent, 
                             const edm::EventSetup& iSetup)
{

  using namespace reco;
  using namespace std;
  using namespace anaobj;

  eventcounter_++;
  /*
  //get tag jets
  edm::Handle<reco::RecoChargedCandidateCollection> tagJetHandle ;
  iEvent.getByLabel (tagJetLabel_, tagJetHandle) ;
//  fillgraph (tagJet_eta_, tagJet_phi_, tagJet_graph_, tagJetHandle, 24, 1) ;

  tagJet_eta_ = new double [tagJetHandle->size ()] ;
  tagJet_phi_ = new double [tagJetHandle->size ()] ;

  for (int it = 0 ; it < tagJetHandle->size () ; ++it)
    {
      tagJet_eta_[it] = (*tagJetHandle)[it].eta () ;
      tagJet_phi_[it] = (*tagJetHandle)[it].phi () ;
    }
  tagJet_graph_ = new TGraph (tagJetHandle->size (), tagJet_eta_, tagJet_phi_) ;
  tagJet_graph_->SetMarkerStyle (24) ;
  tagJet_graph_->SetMarkerColor (4) ;
  */

  // get the other jets
  edm::Handle<reco::CaloJetCollection> jetHandle ;
  iEvent.getByLabel (jetLabel_, jetHandle) ;

  jet_eta_ = new double [jetHandle->size ()] ;
  jet_phi_ = new double [jetHandle->size ()] ;

  int index = 0;
  reco::CaloJetCollection::const_iterator jet_itr = jetHandle->begin(); 
  for ( ; jet_itr != jetHandle->end(); ++jet_itr, index++ ) {
      jet_eta_[index] = jet_itr->eta () ;
      jet_phi_[index] = jet_itr->phi () ;
  }
  jet_graph_ = new TGraph (jetHandle->size (), jet_eta_, jet_phi_) ;
  jet_graph_->SetMarkerStyle (25) ; // empty square
  jet_graph_->SetMarkerColor (1) ;

  // get the electrons collection
  edm::Handle<reco::GsfElectronCollection> electronHandle;
  iEvent.getByLabel(electronLabel_,electronHandle);

  electron_eta_ = new double [electronHandle->size ()] ;
  electron_phi_ = new double [electronHandle->size ()] ;

  index = 0;
  reco::GsfElectronCollection::const_iterator electron_itr = electronHandle->begin();
  for ( ; electron_itr != electronHandle->end(); ++electron_itr, index++ ) {
    electron_eta_[index] = electron_itr->eta () ;
    electron_phi_[index] = electron_itr->phi () ;
  }
  electron_graph_ = new TGraph (electronHandle->size (), electron_eta_, electron_phi_) ;
  electron_graph_->SetMarkerStyle (30) ; // empty star
  electron_graph_->SetMarkerColor (3) ;

  // get the muons collection
  edm::Handle<reco::MuonCollection> muonHandle;
  iEvent.getByLabel(muonLabel_,muonHandle);

  muon_eta_ = new double [muonHandle->size ()] ;
  muon_phi_ = new double [muonHandle->size ()] ;

  index = 0;
  reco::MuonCollection::const_iterator muon_itr = muonHandle->begin(); 
  for ( ; muon_itr != muonHandle->end(); ++muon_itr, index++ ) {
    muon_eta_[index] = muon_itr->eta () ;
    muon_phi_[index] = muon_itr->phi () ;
  }
  muon_graph_ = new TGraph (muonHandle->size (), muon_eta_, muon_phi_) ;
  muon_graph_->SetMarkerStyle (26) ; // empty triangle
  muon_graph_->SetMarkerColor (2) ;

  // get the calo MET
  edm::Handle<reco::CaloMETCollection> metHandle;
  iEvent.getByLabel (metLabel_, metHandle);
  const reco::CaloMETCollection * metCollection = metHandle.product();

  MET_eta_ = new double [1] ;
  MET_phi_ = new double [1] ;

  reco::CaloMETCollection::const_iterator met = metCollection->begin();
  for ( ; met != metCollection->end(); ++met ) {
    MET_eta_[0] = met->eta () ;
    MET_phi_[0] = met->phi () ;
  }
  MET_graph_ = new TGraph (1, MET_eta_, MET_phi_) ;
  MET_graph_->SetMarkerStyle (24) ; // empty circle 
  MET_graph_->SetMarkerColor (8) ;

  // MC particles
  TMarker * mcpoints[6] = {0,0,0,0,0,0};
  edm::Handle< GenParticleCollection > genParticleHandle ; 
  iEvent.getByLabel (genParticleLabel_, genParticleHandle) ;
  
  // loop over generated particles  
  int num = 0 ;
  int particleIndex = 0;
  for (GenParticleCollection::const_iterator MCparticle = genParticleHandle->begin () ; 
       MCparticle != genParticleHandle->end () ;  ++MCparticle, particleIndex++)  {
    int tmpParticleID = MCparticle->pdgId () ;
    int tmpParticleSTATUS = MCparticle->status () ;
    if (abs (tmpParticleID) == pythiaZ_ && tmpParticleSTATUS == 3) {
      for (unsigned int daughterIndex = 0; daughterIndex < MCparticle->numberOfDaughters () ; ++daughterIndex) {
	int tmpDaughterID = MCparticle->daughter (daughterIndex)->pdgId () ;    
	if (abs (tmpDaughterID) == pythiae_) {
	  mcpoints[num] = new TMarker (MCparticle->daughter (daughterIndex)->eta (),
				       MCparticle->daughter (daughterIndex)->phi (), 21) ;
	  mcpoints[num]->SetMarkerSize (2) ;
	  mcpoints[num]->SetMarkerColor (4) ;
	  ++num ;
	}
	else if (abs (tmpDaughterID) == pythiamu_) {
	  mcpoints[num] = new TMarker (MCparticle->daughter (daughterIndex)->eta (),
				       MCparticle->daughter (daughterIndex)->phi (), 21) ;
	  mcpoints[num]->SetMarkerSize (2) ;
	  mcpoints[num]->SetMarkerColor (4) ;
	  ++num ;
	}
	else if (abs (tmpDaughterID) <= pythiat_) {
	  mcpoints[num] = new TMarker (MCparticle->daughter (daughterIndex)->eta (),
				       MCparticle->daughter (daughterIndex)->phi (), 20) ;
	  mcpoints[num]->SetMarkerSize (2) ;
	  mcpoints[num]->SetMarkerColor (6) ;
	  ++num ;
	}
      }
    } // Z
    
    if ( particleIndex == 7 || particleIndex == 8 ){
      mcpoints[num] = new TMarker (MCparticle->eta (),
				   MCparticle->phi (), 20) ;
      mcpoints[num]->SetMarkerSize (2) ;
      mcpoints[num]->SetMarkerColor (2) ;
      ++num ;
    }

    if ((tmpParticleID == pythiaH_) && (tmpParticleSTATUS == 3)) {
      if (MCparticle->numberOfMothers () < 1) continue ;
      const Candidate * interact0 = MCparticle->mother (0) ;
      if (interact0->numberOfDaughters () < 2) continue ;
      
      mcpoints[num] = new TMarker (interact0->daughter (0)->eta (),
				   interact0->daughter (0)->phi (), 20) ;
      mcpoints[num]->SetMarkerSize (3) ;
      mcpoints[num]->SetMarkerColor (2) ;
      ++num ;
      mcpoints[num] = new TMarker (interact0->daughter (1)->eta (),
				   interact0->daughter (1)->phi (), 20) ;
      mcpoints[num]->SetMarkerSize (3) ;
      mcpoints[num]->SetMarkerColor (21) ;
      ++num ;
    } // look for jet tags
  } // loop over generated particles
  
// -------------------------------------------------

  canvas_->cd () ;
  histo2D_->Draw () ;
  for (int i = 0 ; i < 6 ; ++i) if(mcpoints[i]!=0) mcpoints[i]->Draw () ;
  //  tagJet_graph_->Draw ("P") ;
  if(jet_graph_->GetN() != 0)      jet_graph_->Draw ("P") ;
  if(electron_graph_->GetN() != 0) electron_graph_->Draw ("P") ;
  if(muon_graph_->GetN() != 0)     muon_graph_->Draw ("P") ;
  if(MET_graph_->GetN() != 0)      MET_graph_->Draw ("P") ;
  legend_ = new TLegend (0.9, 0.85, 1., 1.) ;
  char labelName[10] ;

  //  sprintf (labelName,"TAG %d",tagJetHandle->size ()) ;
  //  legend_->AddEntry (tagJet_graph_,labelName,"P") ;
  sprintf (labelName,"JET %d",jetHandle->size ()) ;   
  legend_->AddEntry (jet_graph_,labelName,"P") ;
  sprintf (labelName,"ELECTRON %d",electronHandle->size ()) ;    
  legend_->AddEntry (electron_graph_,labelName,"P") ;
  sprintf (labelName,"MUON %d",muonHandle->size ()) ;   
  legend_->AddEntry (muon_graph_,labelName,"P") ;
  legend_->AddEntry (MET_graph_,"MET","P") ;
  legend_->Draw () ;
  char name[30] ;
  sprintf (name,"plot/event_%d_%d.gif",eventcounter_, iEvent.id ().run ()) ;
  canvas_->Print (name,"gif") ;
  
  delete legend_ ;
  for (int i = 0 ; i < 6 ; ++i) delete mcpoints[i] ;
//  delete [] tagJet_eta_ ;
//  delete [] tagJet_phi_ ;
//  delete tagJet_graph_ ;

  delete [] jet_eta_ ;
  delete [] jet_phi_ ;
  delete jet_graph_ ;

  delete [] electron_eta_ ;
  delete [] electron_phi_ ;
  delete electron_graph_ ;

  delete [] muon_eta_ ;
  delete [] muon_phi_ ;
  delete muon_graph_ ;

  delete [] MET_eta_ ;
  delete [] MET_phi_ ;
  delete MET_graph_ ;
}


// --------------------------------------------------------------------


void 
VBFHZZllbbDisplay::beginJob (const edm::EventSetup&)
{
  eventcounter_ = 0;

//  gROOT->SetStyle ("Plain") ;
  canvas_ = new TCanvas () ;
  canvas_->SetGrid () ;
  histo2D_ = new TH2F ("bkg", "",10,-5,5,10,-3.14,3.14) ;
  histo2D_->SetStats (0) ;
  histo2D_->GetXaxis ()->SetTitle ("#eta") ;
  histo2D_->GetYaxis ()->SetTitle ("#phi") ;

}


// --------------------------------------------------------------------


void 
VBFHZZllbbDisplay::endJob () 
{
}

//define this as a plug-in
//DEFINE_FWK_MODULE(VBFHZZllbbDisplay);

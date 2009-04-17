/* \class VBFHZZllbbCommonPreselection
 *
 *
 * H->ZZ->llbb analysis preselection:
 * skim input
 * electron and muon selection
 * m_ll, m4l constraints
 * loose isolation on electrons and muons
 */


// system include files
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbCommonPreselection.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include <DataFormats/EgammaCandidates/interface/GsfElectron.h>
#include <DataFormats/EgammaCandidates/interface/GsfElectronFwd.h>
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <memory>

// namespaces
using namespace edm;
using namespace std;
using namespace reco;

struct SortCandByDecreasingPt {
  bool operator()( const Candidate &c1, const Candidate &c2) const {
    return c1.pt() > c2.pt();
  }
};

// Constructor
VBFHZZllbbCommonPreselection::VBFHZZllbbCommonPreselection(const edm::ParameterSet& pset) {

  // Decay Channel
  decaychannel = pset.getParameter<std::string>("decaychannel");
  hLabel_      = pset.getParameter<edm::InputTag>("HLabel");
  leptonLabel_ = pset.getParameter<edm::InputTag>("LeptonLabel");
  zToLLLabel_  = pset.getParameter<edm::InputTag>("ZllLabel");
  leptonLooseIsolLabel_     = pset.getParameter<edm::InputTag>("leptonLooseIsolLabel");

  edm::ParameterSet cutsConf = pset.getParameter<edm::ParameterSet>("cuts");
  nLeptonCut_             = cutsConf.getParameter<int>   ( "nLeptonCut"            );
  llMassCut_              = cutsConf.getParameter<double>( "llMassCut"             );
  fourBodyMassCut_        = cutsConf.getParameter<double>( "fourBodyMassCut"       );
  numberOfLeptonCombsCut_ = cutsConf.getParameter<int>   ( "numberOfLeptonCombsCut");
  numberOf2l2bCombsCut_   = cutsConf.getParameter<int>   ( "numberOf2l2bCombsCut"  );
  nLooseLeptonCut_        = cutsConf.getParameter<int>   ( "nLooseLeptonCut"       );

  cout << "Starting preselection for channel " << decaychannel << endl;

  string alias;
  if (decaychannel=="2e2b"){
    produces<bool> (alias = decaychannel + "PreselAtleast2Ele"  ).setBranchAlias( alias );
    produces<bool> (alias = decaychannel + "PreselAtleast1ZEE"  ).setBranchAlias( alias );
    produces<bool> (alias = decaychannel + "PreselLoose2IsolEle").setBranchAlias( alias );
  }
  else if (decaychannel=="2mu2b"){
    produces<bool> (alias = decaychannel + "PreselAtleast2Mu"   ).setBranchAlias( alias );
    produces<bool> (alias = decaychannel + "PreselAtleast1ZMuMu").setBranchAlias( alias );
    produces<bool> (alias = decaychannel + "PreselLoose2IsolMu" ).setBranchAlias( alias );
  }
  produces<bool> (alias = decaychannel + "PreselAtleast1H"      ).setBranchAlias( alias );
  produces<bool> (alias = decaychannel + "Presel"               ).setBranchAlias( alias );
}


// Destructor
VBFHZZllbbCommonPreselection::~VBFHZZllbbCommonPreselection() {

}


// Filter event (event preselection)
void VBFHZZllbbCommonPreselection::produce(edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  //2e2b
  auto_ptr<bool> atleast2Ele   ( new bool );
  auto_ptr<bool> atleast1ZEE   ( new bool );
  auto_ptr<bool> Loose2IsoEle  ( new bool );
  
  *atleast2Ele   = false;
  *atleast1ZEE   = false;
  *Loose2IsoEle  = false;
  
  //2mu2b
  auto_ptr<bool> atleast2Mu    ( new bool);
  auto_ptr<bool> atleast1ZMuMu ( new bool );
  auto_ptr<bool> Loose2IsoMu   ( new bool );

  *atleast2Mu    = false;
  *atleast1ZMuMu = false;
  *Loose2IsoMu   = false;
    
  auto_ptr<bool> atleast1H ( new bool );
  *atleast1H     = false;

  // Preselection flag
  auto_ptr<bool> boolPresel ( new bool );
  *boolPresel    = false;

   // Selected Electrons
  int posElectronCounter=0,negElectronCounter=0;
  if ( decaychannel=="2e2b" ) {
    edm::Handle<edm::View<GsfElectron> > electronsHandle;
    iEvent.getByLabel(leptonLabel_, electronsHandle);
    edm::View<GsfElectron>::const_iterator electron_itr = electronsHandle->begin();
    for ( ; electron_itr != electronsHandle->end(); ++electron_itr ) {
      std::cout << "Electron w/ pt= " <<  electron_itr->pt() 
		<< " and eta" << electron_itr->eta() 
		<< " p=" <<  electron_itr->p() << std::endl;
      if      ( electron_itr->charge() > 0 ) posElectronCounter++;
      else if ( electron_itr->charge() < 0 ) negElectronCounter++;
    }
  }

  // Selected Muons
  int posMuonCounter=0,negMuonCounter=0;
  if (decaychannel=="2mu2b" ) {
    Handle<edm::View<Muon> > muonsHandle;
    iEvent.getByLabel(leptonLabel_, muonsHandle);
    edm::View<Muon>::const_iterator muon_itr = muonsHandle->begin();
    for ( ; muon_itr != muonsHandle->end(); ++muon_itr ) {
      std::cout << "Muon with pt= " <<  muon_itr->pt() 
		<< " and eta" << muon_itr->eta() 
		<< " p=" <<  muon_itr->p() << std::endl;
      if      ( muon_itr->charge() > 0 ) posMuonCounter++;
      else if ( muon_itr->charge() < 0 ) negMuonCounter++;
    }
  }
        
  /// channel conditions
  int nElectron = 0;
  int nMuon     = 0;
  if (decaychannel=="2e2b") {
    if ( posElectronCounter>=1 && negElectronCounter>=1 ) {
      nElectron=posElectronCounter+negElectronCounter;
    }
  }
  if (decaychannel=="2mu2b") {
    if ( posMuonCounter>=1 && negMuonCounter>=1 ) {
      nMuon=posMuonCounter+negMuonCounter;
    }
  }

  // Pairs of LL
  int nZEEcounter   = 0;
  int nZMuMucounter = 0;
  Handle<CompositeCandidateCollection> zLLCandidates;
  iEvent.getByLabel(zToLLLabel_, zLLCandidates);    
  CompositeCandidateCollection::const_iterator zIter = zLLCandidates->begin();
  for ( ; zIter!= zLLCandidates->end(); ++zIter ) {
    if ( zIter->p4().mass() > llMassCut_ ){ 
      if (decaychannel=="2e2b"  ) nZEEcounter++;
      if (decaychannel=="2e2mu" ) nZMuMucounter++;
    }  
  }

  // 4 body combinations
  Handle<CompositeCandidateCollection> higgsCandidates;
  iEvent.getByLabel(hLabel_, higgsCandidates);
  int nHiggsCounter=0;
  CompositeCandidateCollection::const_iterator hIter=higgsCandidates->begin();
  for ( ; hIter!= higgsCandidates->end(); ++hIter ) {
    if ( hIter->p4().mass() > fourBodyMassCut_ ) nHiggsCounter++;
  }    

  // Loose isolation for electrons
  int nLooseIsolElectronCounter = 0;
  if (decaychannel=="2e2mu" ) {
    Handle<edm::View<GsfElectron> > looseIsoElectrons;
    iEvent.getByLabel(leptonLooseIsolLabel_, looseIsoElectrons);
    edm::View<GsfElectron>::const_iterator electron_itr = looseIsoElectrons->begin();
    for ( ; electron_itr != looseIsoElectrons->end(); ++electron_itr ) {
      nLooseIsolElectronCounter++;
    }
    cout <<"Number of loose isolated electrons= " << nLooseIsolElectronCounter << endl;
  }


  // Loose isolation for muons
  int nLooseIsolMuonCounter = 0;
  if (decaychannel=="2e2mu" ) {
    Handle<edm::View<Muon> > looseIsoMuons;
    iEvent.getByLabel(leptonLooseIsolLabel_.label(), looseIsoMuons);
    edm::View<Muon>::const_iterator muon_itr = looseIsoMuons->begin();
    for ( ; muon_itr != looseIsoMuons->end(); ++muon_itr ) {
      nLooseIsolMuonCounter++;
    }
    cout <<"Number of loose isolated muons= " << nLooseIsolMuonCounter << endl;
  } 
  
  // 4Leptons channel
  if (decaychannel=="2e2b"){
    if (nElectron                 >= nLeptonCut_             ) *atleast2Ele  = true;
    if (nZEEcounter               >= numberOfLeptonCombsCut_ ) *atleast1ZEE  = true;
    if (nLooseIsolElectronCounter >= nLooseLeptonCut_        ) *Loose2IsoEle = true;
    iEvent.put(atleast2Ele,  decaychannel + "PreselAtleast2Ele"  );
    iEvent.put(atleast1ZEE,  decaychannel + "PreselAtleast1ZEE"  );
    iEvent.put(Loose2IsoEle, decaychannel + "PreselLoose2IsolEle");

    if ( (nElectron                 >= nLeptonCut_             ) &&
	 (nZEEcounter               >= numberOfLeptonCombsCut_ ) && 
	 (nHiggsCounter             >= numberOf2l2bCombsCut_   ) &&
         (nLooseIsolElectronCounter >= nLooseLeptonCut_        )    ) *boolPresel=true;
  }
  else if (decaychannel=="2mu2b"){
    if (nMuon                  >= nLeptonCut_             ) *atleast2Mu    = true;
    if (nZMuMucounter          >= numberOfLeptonCombsCut_ ) *atleast1ZMuMu = true;
    if (nLooseIsolMuonCounter  >= nLooseLeptonCut_        ) *Loose2IsoMu   = true;
    iEvent.put(atleast2Mu,    decaychannel + "PreselAtleast2Mu"   ); 
    iEvent.put(atleast1ZMuMu, decaychannel + "PreselAtleast1ZMuMu");
    iEvent.put(Loose2IsoMu,   decaychannel + "PreselLoose2IsolMu" );   

    if ( (nMuon                 >= nLeptonCut_             ) &&
	 (nZMuMucounter         >= numberOfLeptonCombsCut_ ) && 
	 (nHiggsCounter         >= numberOf2l2bCombsCut_   ) &&
	 (nLooseIsolMuonCounter >= nLooseLeptonCut_        )    ) *boolPresel=true;
  }
         
  if (nHiggsCounter >= numberOf2l2bCombsCut_) *atleast1H = true;
  iEvent.put(atleast1H,  decaychannel + "PreselAtleast1H");
  iEvent.put(boolPresel, decaychannel + "Presel"         );


}

void VBFHZZllbbCommonPreselection::beginJob(const edm::EventSetup& iSetup) {
  cout << "Starting preselection" << endl;
}

void VBFHZZllbbCommonPreselection::endJob() {
  cout << "Create preselection variables" << endl;
}



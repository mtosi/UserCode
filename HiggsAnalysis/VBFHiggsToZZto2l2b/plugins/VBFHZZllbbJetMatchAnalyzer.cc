#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbJetMatchAnalyzer.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"

using namespace std;
using namespace reco;
using namespace edm;
using namespace ROOT::Math::VectorUtil;

VBFHZZllbbJetMatchAnalyzer::VBFHZZllbbJetMatchAnalyzer(const edm::ParameterSet& iConfig) :

  source_          ( iConfig.getParameter<InputTag> ("src")          ),
  matched_         ( iConfig.getParameter<InputTag> ("matched")      ),
  matchedjetsOne_  ( iConfig.getParameter<InputTag> ("matchMapOne")  ),
  matchedjetsMany_ ( iConfig.getParameter<InputTag> ("matchMapMany") )
{

  edm::Service<TFileService> fs ;

}

void VBFHZZllbbJetMatchAnalyzer::beginJob( const edm::EventSetup& iSetup)
{
  eventcounter_             = 0;

  deltaR_          = fs->make<TH1D> ( "deltaR",          "#DeltaR between matched jets",                      100,   0.,  5. );
  deltaPt_         = fs->make<TH1D> ( "deltaPt",         "#Deltap_{T} between matched jets",                  100, -20., 20. );
  resPt_           = fs->make<TH1D> ( "resPt",           "p_{T} resolution between matched jets",             100,  -1.,  1. );
  deltaPtVSdeltaR_ = fs->make<TH2D> ( "deltaPtVSdeltaR", "#Deltap_{T} VS #DeltaR between matched jets",       100,   0.,  5. , 100, -20., 20. );
  resPtVSdeltaR_   = fs->make<TH2D> ( "resPtVSdeltaR",    "p_{T} resolution VS #DeltaR between matched jets", 100,   0.,  5. , 100,  -1.,  1. );

}

void VBFHZZllbbJetMatchAnalyzer::endJob()
{      

}

void VBFHZZllbbJetMatchAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  cout << "***********************************************" << std::endl;
  eventcounter_++;
  if ( eventcounter_/1000 == float(eventcounter_)/100. ) {
    std::cout << "Event number " << eventcounter_ << std::endl;
  }

  //  cout << "[GenJetTest] analysing event " << iEvent.id() << endl;
  
  try {
    iEvent.getByLabel (source_,source);
    iEvent.getByLabel (matched_,matched);
    iEvent.getByLabel (matchedjetsOne_ , matchedjetsOne );
    iEvent.getByLabel (matchedjetsMany_, matchedjetsMany);
  } catch(std::exception& ce) {
    cerr << "[GenJetTest] caught std::exception " << ce.what() << endl;
    return;
  }
  
  //
  // Printout for OneToOne matcher
  //
  double dR    = -1.;
  double dPt   = -99.;
  double resPt = -99.; 
  cout << "**********************" << endl;
  cout << "* OneToOne Printout  *" << endl;
  cout << "**********************" << endl;
  for( CandViewMatchMap::const_iterator f  = matchedjetsOne->begin();
                                        f != matchedjetsOne->end();
                                        f++) {

      const Candidate *sourceRef = &*(f->key);
      const Candidate *matchRef  = &*(f->val);
      dR= DeltaR( sourceRef->p4() , matchRef->p4() );
      dPt = sourceRef->pt() - matchRef->pt();
      resPt = (sourceRef->pt() - matchRef->pt())/sourceRef->pt();
      printf("[GenJetTest] (pt,eta,phi) source = %6.2f %5.2f %5.2f matched = %6.2f %5.2f %5.2f dR=%5.3f\n",
	     sourceRef->et(),
	     sourceRef->eta(),
	     sourceRef->phi(), 
	     matchRef->et(), 
	     matchRef->eta(),
	     matchRef->phi(), 
	     dR);
      deltaR_          -> Fill ( dR        );
      deltaPt_         -> Fill ( dPt       );
      resPt_           -> Fill ( resPt     );
      deltaPtVSdeltaR_ -> Fill ( dR, dPt   );
      resPtVSdeltaR_   -> Fill ( dR, resPt );
  }
  //
  // Printout for OneToMany matcher
  //
  cout << "**********************" << endl;
  cout << "* OneToMany Printout *" << endl;
  cout << "**********************" << endl;
  for( CandMatchMapMany::const_iterator f  = matchedjetsMany->begin();
                                        f != matchedjetsMany->end();
                                        f++) {
    const Candidate *sourceRef = &*(f->key);

    printf("[GenJetTest] (pt,eta,phi) source = %6.2f %5.2f %5.2f\n",
	   sourceRef->et(),
	   sourceRef->eta(),
	   sourceRef->phi()  );

    const vector< pair<CandidateRef, double> > vectMatched = f->val;
    vector< pair<CandidateRef, double> >::const_iterator matchIT;
    for ( matchIT = vectMatched.begin(); matchIT != vectMatched.end(); matchIT++) {
      const Candidate *matchedRef = &*( (*matchIT).first );
      double deltaR = (*matchIT).second;
      printf("             (pt,eta,phi) matched = %6.2f %5.2f %5.2f - dR=%5.2f\n",
	     matchedRef->et(),
	     matchedRef->eta(),
	     matchedRef->phi(),
	     deltaR);
      
    }
  }
  cout << "***********************************************" << std::endl;

}


#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbJetMatchAnalyzer.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"

using namespace std;
using namespace reco;
using namespace edm;
using namespace ROOT::Math::VectorUtil;

VBFHZZllbbJetMatchAnalyzer::VBFHZZllbbJetMatchAnalyzer(const edm::ParameterSet& iConfig) :

  source_          ( iConfig.getParameter<InputTag> ("source")         ),
  matched_         ( iConfig.getParameter<InputTag> ("matched")        ),
  matchedjetsOne_  ( iConfig.getParameter<InputTag> ("matchMapOne")    ),
  matchedjetsMany_ ( iConfig.getParameter<InputTag> ("matchMapMany")   ),
  histomakerFlag_  ( iConfig.getParameter<bool>     ("histomakerFlag") ),
  dRcut_           ( iConfig.getParameter<double>   ("dRcut")          )
{

  if (histomakerFlag_)
    edm::Service<TFileService> fs ;

}

void VBFHZZllbbJetMatchAnalyzer::beginJob( const edm::EventSetup& iSetup)
{

  eventcounter_              = 0;
  eventmatchedcounter_       = 0;
  eventmatchedcounter_dRcut_ = 0;

  if (histomakerFlag_) {
    deltaR_          = fs->make<TH1D> ( "deltaR",          "#DeltaR between matched jets",                            100,   0.,  5. );
    totalLenght_     = fs->make<TH1D> ( "totalLenght",     "Total lenght",                                            100,   0.,  5. );
    deltaPt_         = fs->make<TH1D> ( "deltaPt",         "#Deltap_{T} between matched jets",                        100, -20., 20. );
    resPt_           = fs->make<TH1D> ( "resPt",           "p_{T} resolution between matched jets",                   100,  -1.,  1. );
    deltaPtVSdeltaR_ = fs->make<TH2D> ( "deltaPtVSdeltaR", "#Deltap_{T} VS #DeltaR between matched jets",             100,   0.,  5. , 100, -20., 20. );
    resPtVSdeltaR_   = fs->make<TH2D> ( "resPtVSdeltaR",    "p_{T} resolution VS #DeltaR between matched jets",       100,   0.,  5. , 100,  -1.,  1. );
    deltaPt_dRcut_   = fs->make<TH1D> ( "deltaPt_dRcut",    "#Deltap_{T} between matched jets w/ #DeltaR<=dRcut",     100, -20., 20. );
    resPt_dRcut_     = fs->make<TH1D> ( "resPt_dRcut",      "p_{T} resolution between matched jets w/ #DeltaR<=dRcut",100,  -1.,  1. );
  }
}

void VBFHZZllbbJetMatchAnalyzer::endJob()
{      
  cout << "--------------------------------------------------" << std::endl;
  std::cout << "number of matched jets:              " << eventmatchedcounter_       << std::endl;
  std::cout << "number of matched jets w/ dR<=dRcut: " << eventmatchedcounter_dRcut_ << std::endl;
  std::cout << " --> efficiency: " << double(eventmatchedcounter_dRcut_)/double(eventmatchedcounter_);
  cout << "--------------------------------------------------" << std::endl;
}

void VBFHZZllbbJetMatchAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

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
  unsigned matchedcounter       = 0;
  unsigned matchedcounter_dRcut = 0;
  double totalLenght = 0.;

  cout << "**********************" << endl;
  cout << "* OneToOne Printout  *" << endl;
  cout << "**********************" << endl;
  for( CandViewMatchMap::const_iterator f  = matchedjetsOne->begin();
                                        f != matchedjetsOne->end();
                                        f++) {
    matchedcounter++;
    eventmatchedcounter_++;
    const Candidate *sourceRef = &*(f->key);
    const Candidate *matchRef  = &*(f->val);
    double dR= DeltaR( sourceRef->p4() , matchRef->p4() );
    double dPt = sourceRef->pt() - matchRef->pt();
    double resPt = (sourceRef->pt() - matchRef->pt())/sourceRef->pt();
    totalLenght += dR;
    printf("[GenJetTest] (pt,eta,phi) source = %6.2f %5.2f %5.2f matched = %6.2f %5.2f %5.2f dR=%5.3f\n",
	   sourceRef->et(),
	   sourceRef->eta(),
	   sourceRef->phi(), 
	   matchRef->et(), 
	   matchRef->eta(),
	   matchRef->phi(), 
	   dR);
    if (histomakerFlag_) {
      deltaR_          -> Fill ( dR        );
      deltaPt_         -> Fill ( dPt       );
      resPt_           -> Fill ( resPt     );
      deltaPtVSdeltaR_ -> Fill ( dR, dPt   );
      resPtVSdeltaR_   -> Fill ( dR, resPt );
    }      
    if (dR<=dRcut_) {
      matchedcounter_dRcut++;
      eventmatchedcounter_dRcut_++;
      if (histomakerFlag_) {
	deltaPt_dRcut_ -> Fill ( dPt   );
	resPt_dRcut_   -> Fill ( resPt );
      }
    }
  }

  totalLenght_ -> Fill ( totalLenght );

  std::cout << "number of matched jets:              " << matchedcounter       << std::endl;
  std::cout << "number of matched jets w/ dR<=dRcut: " << matchedcounter_dRcut << std::endl;
  std::cout << " --> efficiency: " << double(matchedcounter_dRcut)/double(matchedcounter);

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
      printf("             (pt,eta,phi) matched = %6.2f %5.2f %5.2f - deltaR=%5.2f\n",
	     matchedRef->et(),
	     matchedRef->eta(),
	     matchedRef->phi(),
	     deltaR);
      
    }
  }
}


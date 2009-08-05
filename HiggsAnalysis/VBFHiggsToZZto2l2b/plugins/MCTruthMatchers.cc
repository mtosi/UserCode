#include "PhysicsTools/UtilAlgos/interface/PhysObjectMatcher.h"
#include "PhysicsTools/UtilAlgos/interface/MCMatchSelector.h"
#include "PhysicsTools/UtilAlgos/interface/MatchByDRDPt.h"
#include "PhysicsTools/UtilAlgos/interface/MatchLessByDPt.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

// Match by deltaR and deltaPt, ranking by deltaR (default)
typedef reco::PhysObjectMatcher<
  reco::CandidateView,
  reco::GenParticleCollection,
  reco::MCMatchSelector<reco::CandidateView::value_type,
			reco::GenParticleCollection::value_type>,
  reco::MatchByDRDPt<reco::CandidateView::value_type,
		     reco::GenParticleCollection::value_type>
> MCMatcher;

// Alternative: match by deltaR and deltaPt, ranking by deltaPt
typedef reco::PhysObjectMatcher<
  reco::CandidateView,
  reco::GenParticleCollection,
  reco::MCMatchSelector<reco::CandidateView::value_type,
			reco::GenParticleCollection::value_type>,
  reco::MatchByDRDPt<reco::CandidateView::value_type,
		     reco::GenParticleCollection::value_type>,
  reco::MatchLessByDPt<reco::CandidateView,
		       reco::GenParticleCollection>
> MCMatcherByPt;

// JET Match by deltaR, ranking by deltaR (default)
typedef reco::PhysObjectMatcher<
  reco::CandidateView,
  reco::GenJetCollection,
  reco::MCMatchSelector<reco::CandidateView::value_type,
			reco::GenJetCollection::value_type>,
  reco::MatchByDR<reco::CandidateView::value_type,
		  reco::CandidateView::value_type>
> GenJetMatcher;

// Calo JET Match by deltaR, ranking by deltaR (default)
typedef reco::PhysObjectMatcher<
  reco::CandidateView,
  reco::CaloJetCollection,
  reco::MCMatchSelector<reco::CandidateView::value_type,
			reco::CaloJetCollection::value_type>,
  reco::MatchByDR<reco::CandidateView::value_type,
		  reco::CandidateView::value_type>
> CaloJetMatcher;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MCMatcher);
DEFINE_FWK_MODULE(MCMatcherByPt);
DEFINE_FWK_MODULE(GenJetMatcher);
DEFINE_FWK_MODULE(CaloJetMatcher);

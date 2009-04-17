#ifndef VBFHZZllbbUTILS_H
#define VBFHZZllbbUTILS_H

// system include files
#include <memory>
#include <string>

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"

namespace vbfhzz2l2b
{

  enum { HIGHEFF       = 1,
	 HIGHPUR       = 2,
	 COMBSECVTX    = 0,
	 COMBSECVTXMVA = 4,
	 SOFTMUON      = 5,
	 SOFTELECTRON  = 6,
	 JETPROB       = 3 };
  
  
  int bTaggerCode ( const std::string& bTagger ) {
    
    int code = -1;
    if ( bTagger == "HIGHEFF" )            code = HIGHEFF;
    else if ( bTagger == "HIGHPUR"       ) code = HIGHPUR;
    else if ( bTagger == "COMBSECVTX"    ) code = COMBSECVTX;
    else if ( bTagger == "COMBSECVTXMVA" ) code = COMBSECVTXMVA;
    else if ( bTagger == "SOFTMUON"      ) code = SOFTMUON;
    else if ( bTagger == "SOFTELECTRON"  ) code = SOFTELECTRON;
    else if ( bTagger == "JETPROB"       ) code = JETPROB;
    else
      std::cout << "[VBFHZZllbbUtils::bTaggerCode] --> WARNING: bTagger " << bTagger << " NOT IMPLEMENTED!" << std::endl;
    
    return code;
 }
  
  


typedef reco::CaloJetCollection::const_iterator caloJetItr ;

template <typename T>
std::pair<T,T> findCaloJetsPair_maxInvMass (T &begin, T &end,
					    double jetPtMin,
					    double jetEtaMax) {

  std::pair<T,T> jetsPair (begin,begin) ;
  double maxInvMass = 0. ;

  // first loop over jets
  for ( T firstJet = begin; firstJet != end; ++firstJet ) {
    if (firstJet->pt() < jetPtMin || 
	fabs( firstJet->eta() ) > jetEtaMax) continue ;

      math::XYZTLorentzVector jetsSumP4 = firstJet->p4 () ;

      // second loop over jets
      for ( T secondJet = firstJet + 1; secondJet != end; ++secondJet ) {
	if (secondJet->pt() < jetPtMin || 
	    fabs( secondJet->eta () ) > jetEtaMax ) continue ;
	
	jetsSumP4 += secondJet->p4 ();

	if (jetsSumP4.M () > maxInvMass) {
	    maxInvMass = jetsSumP4.M () ;
	    jetsPair.first  = firstJet ;
	    jetsPair.second = secondJet ;
	}
      } // end second loop over jets
  } // end first loop over jets
  
  return jetsPair ;


}

template <typename T>
std::pair<T,T> findCaloJetsPair_maxPt (T &begin, T &end,
				       double jetPtMin, 
				       double jetEtaMax) {
  std::pair<caloJetItr,caloJetItr> jetsPair (begin, begin) ; 

  double maxPt1 = 0.;
  double maxPt2 = 0.;

  caloJetItr jet1;
  caloJetItr jet2;

  // loop over jets
  for ( caloJetItr jet = begin; jet != end; ++jet ) {
    if (jet->pt() < jetPtMin || 
	fabs( jet->eta() ) > jetEtaMax) continue ;
    
    if (jet->p4().Pt() > maxPt1) {
      jet2 = jet1;
      jet1 = jet;
      maxPt2 = maxPt1;
      maxPt1 = jet->p4().Pt() ;
    } 
    else if (jet->p4().Pt() > maxPt2) {
      jet2 = jet;
      maxPt2 = jet->p4().Pt() ;
    }
  }
  
  jetsPair.first = jet1;
  jetsPair.second = jet2 ;
  
  return jetsPair;
  
}

template <typename T1, typename T2>
std::pair<T1,T1> findCaloJetsPair_sortT (T1 &begin, T1 &end) {

  std::vector<std::pair<T1,T1> > jetsPairsVector;
  
    // first loop over jets
    for ( T1 firstJet = begin; firstJet != end ; ++firstJet ) {
  
      // second loop over jets
      for ( T1 secondJet = firstJet + 1 ; secondJet != end; ++secondJet ) {
	
	jetsPairsVector.push_back ( make_pair (firstJet, secondJet)) ;
      } // end second loop over jets
    } // end first loop over jetsPairsVector
    
    sort ( jetsPairsVector.begin () , jetsPairsVector.end () , T2 () ) ;
    return jetsPairsVector.back () ;
}
  
 
template <typename T>
struct ptSorting {
  typedef T first_argument_type;
  typedef T second_argument_type;
  bool operator() ( const T &t1, 
		    const T &t2 ) const { 
    return t1.p4().Pt() > t2.p4().Pt(); 
  }
};

template <typename T>
struct maxInvMassSorting {
  bool operator() ( const std::pair<T,T> & couple1,
		    const std::pair<T,T> & couple2 ) const {
    double invMass1 = (couple1.first->p4()+couple1.second->p4()).M();
    double invMass2 = (couple2.first->p4()+couple2.second->p4()).M();
    return ( invMass1 < invMass2 );
  }
};

template <typename T>
struct maxPtSorting { 
  bool operator() ( const std::pair<T,T> & couple1,
		    const std::pair<T,T> & couple2 ) const {
    return ( couple1.first->p4().Pt()  < couple2.first->p4().Pt() &&
	     couple1.second->p4().Pt() < couple2.second->p4().Pt() );
  }
};

template <typename T>
struct maxSumPtSorting {
  bool operator() ( const std::pair<T,T> & couple1,
		    const std::pair<T,T> & couple2 ) {
    double sumPt1 = couple1.first->p4().Pt() + couple1.second->p4().Pt();
    double sumPt2 = couple2.first->p4().Pt() + couple2.second->p4().Pt();
    return ( sumPt1 < sumPt2 );
  }
};



}

#endif // VBFHZZLLBBUTILS_H



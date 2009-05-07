#ifndef VBFHZZllbbUTILS_H
#define VBFHZZllbbUTILS_H

// system include files

#include "TLorentzVector.h"

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
  
  // calculate the distance between two lorentz vectors 
  // using DeltaR(eta, phi) or normal space angle(theta, phi)
  double distance(const math::XYZTLorentzVector&, const math::XYZTLorentzVector&, bool);

  int bTaggerCode ( const std::string& );

  double resolution ( double &, double & );

  bool BhadronTable(int pdgcode);

  void setMomentum (TLorentzVector &, 
		    const reco::Candidate &);


  // *********************************************************
  // *** find couple of objs w/ the highest invariant mass ***
  // *********************************************************
  template <typename T>
  std::pair<T,T> findPair_maxInvMass (const T begin, const T end) {
    std::pair<T,T> objPair (begin,begin) ;
    double maxInvMass = 0. ;
    // first loop over objects
    for ( T firstObj = begin; firstObj != end; ++firstObj ) {
        
      // second loop over objs
      for ( T secondObj = firstObj + 1; secondObj != end; ++secondObj ) {
	
	math::XYZTLorentzVector objsSumP4 = firstObj->p4 () + secondObj->p4 ();
	
	if ( objsSumP4.M () > maxInvMass ) {
	  maxInvMass = objsSumP4.M () ;
	  objPair.first  = firstObj ;
	  objPair.second = secondObj ;
	}
      } // end second loop over objs
    } // end first loop over objs  
    return objPair ;
  }
  // ------------------------------------------------------------


  // ****************************************************
  // *** find couple of objs which satisfy pt min cut ***
  // *** w/ the highest invariant mass                ***
  // ****************************************************
  template <typename T>
  std::pair<T,T> findPair_maxInvMass_ptMinCut (const T begin, const T end,
					       double ptMin1, double ptMin2 = -1.) {

    std::pair<T,T> objPair(begin,begin) ;
    double maxInvMass = 0. ;
    // first loop over objects
    for ( T firstObj = begin; firstObj != end; ++firstObj ) {
      if (firstObj->pt() < ptMin1) continue ;
      
      // second loop over objs
      for ( T secondObj = firstObj + 1; secondObj != end; ++secondObj ) {
	if (secondObj->pt() < ptMin2) continue ;
	
	math::XYZTLorentzVector objsSumP4 = firstObj->p4 () + secondObj->p4 ();

	if ( objsSumP4.M () > maxInvMass ) {
	  maxInvMass = objsSumP4.M () ;
	  objPair.first  = firstObj ;
	  objPair.second = secondObj ;
	}
      } // end second loop over objs
    } // end first loop over objs  
    return objPair ;
  }
  // ------------------------------------------------------------

  template <typename T>
  std::pair<T,T> findPair_maxInvMass_ptMinCut_etaMaxCut (const T begin, const T end,
							 double ptMin1, double ptMin2 = -1.,
							 double etaMax1, double etaMax2 = 10.) {
    std::pair<T,T> objPair (begin,begin) ;
    double maxInvMass = 0. ;
    // first loop over objects
    for ( T firstObj = begin; firstObj != end; ++firstObj ) {
      if (firstObj->pt() < ptMin1 || 
	  fabs( firstObj->eta() ) > etaMax1) continue ;
          
      // second loop over objs
      for ( T secondObj = firstObj + 1; secondObj != end; ++secondObj ) {
	if (secondObj->pt() < ptMin2 || 
	    fabs( secondObj->eta () ) > etaMax2 ) continue ;
	
	math::XYZTLorentzVector objsSumP4 = firstObj->p4 () + secondObj->p4 ();
	
	if ( objsSumP4.M () > maxInvMass ) {
	  maxInvMass = objsSumP4.M () ;
	  objPair.first  = firstObj ;
	  objPair.second = secondObj ;
	}
      } // end second loop over objs
    } // end first loop over objs  
    return objPair ;
  }

  template <typename T>
    std::pair<T,T> findPair_maxPt (T &begin, T &end) {
    
    std::pair<T,T> objPair (begin, begin) ; 
    
    double maxPt1 = 0.;
    double maxPt2 = 0.;

    T obj1, obj2;

    // loop over objs
    for ( T obj = begin; obj != end; ++obj ) {
    
      if (obj->p4().Pt() > maxPt1) {
	obj2 = obj1;
	obj1 = obj;
	maxPt2 = maxPt1;
	maxPt1 = obj->p4().Pt() ;
      } 
      else if (obj->p4().Pt() > maxPt2) {
	obj2 = obj;
	maxPt2 = obj->p4().Pt() ;
      }
    }
  
    objPair.first  = obj1;
    objPair.second = obj2 ;
  
    return objPair;
  }

  template <typename T>
    std::pair<T,T> findPair_maxPt_ptMinCut (T &begin, T &end,
					    double ptMin) {

    std::pair<T,T> objPair (begin, begin) ; 
    
    double maxPt1 = 0.;
    double maxPt2 = 0.;

    T obj1, obj2;

    // loop over objs
    for ( T obj = begin; obj != end; ++obj ) {
      if (obj->pt() < ptMin) continue ;
    
      if (obj->p4().Pt() > maxPt1) {
	obj2 = obj1;
	obj1 = obj;
	maxPt2 = maxPt1;
	maxPt1 = obj->p4().Pt() ;
      } 
      else if (obj->p4().Pt() > maxPt2) {
	obj2 = obj;
	maxPt2 = obj->p4().Pt() ;
      }
    }
    
    objPair.first  = obj1;
    objPair.second = obj2 ;
    
    return objPair;
    
  }

  template <typename T>
    std::pair<T,T> findPair_maxPt_ptMinCut_etaMaxCut (T &begin, T &end,
						      double ptMin, double etaMax) {
    
    std::pair<T,T> objPair (begin, begin) ; 
    
    double maxPt1 = 0.;
    double maxPt2 = 0.;
    
    T obj1, obj2;
    
    // loop over objs
    for ( T obj = begin; obj != end; ++obj ) {
      if (obj->pt() < ptMin || 
	  fabs( obj->eta() ) > etaMax) continue ;
      
      if (obj->p4().Pt() > maxPt1) {
	obj2 = obj1;
	obj1 = obj;
	maxPt2 = maxPt1;
	maxPt1 = obj->p4().Pt() ;
      } 
      else if (obj->p4().Pt() > maxPt2) {
	obj2 = obj;
	maxPt2 = obj->p4().Pt() ;
      }
    }
    objPair.first  = obj1;
    objPair.second = obj2 ;
    
    return objPair;
  }
  
 template<class T>
   std::pair<T,T> findPair_maxInvMass_oppositeEta (T & begin, T & end ) {

   std::pair<T,T> objPair (begin,begin) ;
   double maxInvMass = 0. ;
   // first loop over jets
   for (T firstObj = begin ; firstObj != end ; ++firstObj ) {

      // second loop over objs
      for (T secondObj = firstObj + 1 ; secondObj != end ; ++secondObj ) {
	
	if (firstObj->eta ()*secondObj->eta () > 0) continue ;
	
	math::XYZTLorentzVector objsSumP4 = firstObj->p4 () + secondObj->p4 () ;

	if (objsSumP4.M () > maxInvMass) {
	  maxInvMass = objsSumP4.M () ;
	  objPair.first  = firstObj ;
	  objPair.second = secondObj ;
	}
      } // second loop over objs
   } // first loop over objs
   
   return objPair ;
 }

 template<class T>
   std::pair<T,T> findPair_maxInvMass_oppositeEta_ptMinCut (T & begin, T & end,
						 double ptMin ) {

   std::pair<T,T> objPair (begin,begin) ;
   double maxInvMass = 0. ;
   // first loop over jets
   for (T firstObj = begin ; firstObj != end ; ++firstObj ) {

     if (firstObj->pt () < ptMin) continue ;

      // second loop over objs
      for (T secondObj = firstObj + 1 ; secondObj != end ; ++secondObj ) {
	
	if (secondObj->pt () < ptMin) continue ;
	
	if (firstObj->eta ()*secondObj->eta () > 0) continue ;
	
	math::XYZTLorentzVector objsSumP4 = firstObj->p4 () + secondObj->p4 () ;

	if (objsSumP4.M () > maxInvMass) {
	  maxInvMass = objsSumP4.M () ;
	  objPair.first  = firstObj ;
	  objPair.second = secondObj ;
	}
      } // second loop over objs
   } // first loop over objs
   
   return objPair ;
 }

 template<class T>
   std::pair<T,T> findPair_maxInvMass_oppositeEta_ptMinCut_etaMaxCut (T & begin, T & end,
							   double ptMin, double etaMax) {

   std::pair<T,T> objPair (begin,begin) ;
   double maxInvMass = 0. ;
   // first loop over jets
   for (T firstObj = begin ; firstObj != end ; ++firstObj ) {

      if (firstObj->pt () < ptMin || 
          fabs (firstObj->eta ()) > etaMax) continue ;

      // second loop over objs
      for (T secondObj = firstObj + 1 ; secondObj != end ; ++secondObj ) {
	
	if (secondObj->pt () < ptMin || 
	    fabs (secondObj->eta ()) > etaMax) continue ;
	
	if (firstObj->eta ()*secondObj->eta () > 0) continue ;
	
	math::XYZTLorentzVector objsSumP4 = firstObj->p4 () + secondObj->p4 () ;

	if (objsSumP4.M () > maxInvMass) {
	  maxInvMass = objsSumP4.M () ;
	  objPair.first  = firstObj ;
	  objPair.second = secondObj ;
	}
      } // second loop over objs
   } // first loop over objs
   
   return objPair ;
 }


  // ****************************************************
  // *** find couple of objs w/ the highest delta eta ***
  // ****************************************************
  template <typename T>
  std::pair<T,T> findPair_maxDeltaEta (const T begin, const T end) {

    std::pair<T,T> objPair (begin,begin) ;
    double maxDeltaEta = 0. ;
    // first loop over objects
    for ( T firstObj = begin; firstObj != end; ++firstObj ) {
      
      // second loop over objs
      for ( T secondObj = firstObj + 1; secondObj != end; ++secondObj ) {
	
	double deltaEta = fabs(firstObj->eta() - secondObj->eta());
	
	if ( deltaEta > maxDeltaEta ) {
	  maxDeltaEta = deltaEta ;
	  objPair.first  = firstObj ;
	  objPair.second = secondObj ;
	}
      } // end second loop over objs
    } // end first loop over objs  
    return objPair ;
  }
  // ------------------------------------------------------------

  // ****************************************************
  // *** find couple of objs which satisfy pt min cut ***
  // *** w/ the highest delta eta                     ***
  // ****************************************************
  template <typename T>
  std::pair<T,T> findPair_maxDeltaEta_ptMinCut (const T begin, const T end,
						double ptMin1, double ptMin2 = -1.) {

    std::pair<T,T> objPair (begin,begin) ;
    double maxDeltaEta = 0. ;
    // first loop over objects
    for ( T firstObj = begin; firstObj != end; ++firstObj ) {
      if (firstObj->pt() < ptMin1) continue ;
      
      // second loop over objs
      for ( T secondObj = firstObj + 1; secondObj != end; ++secondObj ) {
	if (secondObj->pt() < ptMin2) continue ;
	
	double deltaEta = fabs(firstObj->eta() - secondObj->eta());
	if ( deltaEta > maxDeltaEta ) {
	  maxDeltaEta = deltaEta ;
	  objPair.first  = firstObj ;
	  objPair.second = secondObj ;
	}
      } // end second loop over objs
    } // end first loop over objs  
    return objPair ;
  }
  // ------------------------------------------------------------


  // *****************************************
  // *** find couple of objs which satisfy ***
  // *** - pt min cut                      ***
  // *** - eta max cut                     ***
  // *** w/ the highest delta eta          ***
  // *****************************************
  template <typename T>
  std::pair<T,T> findPair_maxDeltaEta_ptMinCut_etaMaxCut (const T begin, const T end,
							  double ptMin1, double ptMin2 = -1.,
							  double etaMax1, double etaMax2 = 10.) {

    std::pair<T,T> objPair (begin,begin) ;
    double maxDeltaEta = 0. ;
    // first loop over objects
    for ( T firstObj = begin; firstObj != end; ++firstObj ) {
      if (firstObj->pt() < ptMin1 || 
	  fabs( firstObj->eta() ) > etaMax1) continue ;
      
      // second loop over objs
      for ( T secondObj = firstObj + 1; secondObj != end; ++secondObj ) {
	if (secondObj->pt() < ptMin2 || 
	    fabs( secondObj->eta () ) > etaMax2 ) continue ;
	
	double deltaEta = fabs(firstObj->eta() - secondObj->eta());
	
	if ( deltaEta > maxDeltaEta ) {
	  maxDeltaEta = deltaEta ;
	  objPair.first  = firstObj ;
	  objPair.second = secondObj ;
	}
      } // end second loop over objs
    } // end first loop over objs  
    return objPair ;
  }
  // ------------------------------------------------------------


// --------------------------------------------------------------------

template <class T>
  std::pair<T,T> findObjPair_maxPt_oppositeEta_ptMinCut_etaMaxCut (T & begin, T & end) {

  std::pair<T,T> objPair (begin, begin) ; 

  double ptMax1 = 0;
  double ptMax2 = 0;

  T obj1, obj2;

  for (T obj = begin ; obj != end ; ++obj) {

    if (obj->p4().Pt() > ptMax1) {
      obj2 = obj1;
      obj1 = obj;
      ptMax2 = ptMax1;
      ptMax1 = obj->p4().Pt() ;
    } 
    else if ( (obj->p4().Pt() > ptMax2) && 
	      (obj1->p4().Eta() * obj->p4().Eta() < 0) ) {
      obj2 = obj;
      ptMax2 = obj->p4().Pt() ;
    }
  }
  
  objPair.first  = obj1;
  objPair.second = obj2 ;
  
  return objPair ;
  
}

template <class T>
std::pair<T,T> findObjPair_maxPt_oppositeEta_ptMinCut_etaMaxCut (T & begin, T & end,
								 double ptMin) {

  std::pair<T,T> objPair (begin, begin) ; 

  double ptMax1 = 0;
  double ptMax2 = 0;

  T obj1, obj2;

  for (T obj = begin ; obj != end ; ++obj) {

    if (obj->pt () < ptMin) continue ;

    if (obj->p4().Pt() > ptMax1) {
      obj2 = obj1;
      obj1 = obj;
      ptMax2 = ptMax1;
      ptMax1 = obj->p4().Pt() ;
    } 
    else if ( (obj->p4().Pt() > ptMax2) && 
	      (obj1->p4().Eta() * obj->p4().Eta() < 0) ) {
      obj2 = obj;
      ptMax2 = obj->p4().Pt() ;
    }
  }
  
  objPair.first  = obj1;
  objPair.second = obj2 ;
  
  return objPair ;
  
}

template <class T>
std::pair<T,T> findObjPair_maxPt_oppositeEta_ptMinCut_etaMaxCut (T & begin, T & end,
								 double ptMin, double etaMax) {

  std::pair<T,T> objPair (begin, begin) ; 

  double ptMax1 = 0;
  double ptMax2 = 0;

  T obj1, obj2;

  for (T obj = begin ; obj != end ; ++obj) {

    if (obj->pt () < ptMin || 
	fabs (obj->eta ()) > etaMax) continue ;

    if (obj->p4().Pt() > ptMax1) {
      obj2 = obj1;
      obj1 = obj;
      ptMax2 = ptMax1;
      ptMax1 = obj->p4().Pt() ;
    } 
    else if ( (obj->p4().Pt() > ptMax2) && 
	      (obj1->p4().Eta() * obj->p4().Eta() < 0) ) {
      obj2 = obj;
      ptMax2 = obj->p4().Pt() ;
    }
  }
  
  objPair.first  = obj1;
  objPair.second = obj2 ;
  
  return objPair ;
  
 }

  class PtGreater {
  public:
    template <typename T> bool operator () (const T& i, const T& j) {
      return (i.pt() > j.pt());
    }
  };
  

template <typename T>
struct maxPtSorting {
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
struct maxPairPtSorting { 
  bool operator() ( const std::pair<T,T> & couple1,
		    const std::pair<T,T> & couple2 ) const {
    return ( couple1.first->p4().Pt()  < couple2.first->p4().Pt() &&
	     couple1.second->p4().Pt() < couple2.second->p4().Pt() );
  }
};

template <typename T>
struct maxPairSumPtSorting {
  bool operator() ( const std::pair<T,T> & couple1,
		    const std::pair<T,T> & couple2 ) {
    double sumPt1 = couple1.first->p4().Pt() + couple1.second->p4().Pt();
    double sumPt2 = couple2.first->p4().Pt() + couple2.second->p4().Pt();
    return ( sumPt1 < sumPt2 );
  }
};

}

#endif // VBFHZZLLBBUTILS_H



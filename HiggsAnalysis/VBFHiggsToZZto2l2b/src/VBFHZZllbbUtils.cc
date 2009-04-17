#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/VBFHZZllbbUtils.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include <iostream>

#include "TLorentzVector.h"

namespace vbfhzz2l2b
{

void setMomentum (TLorentzVector & myvector, 
		  const reco::Candidate & gen) {
  myvector.SetPx (gen.px ()) ;
  myvector.SetPy (gen.py ()) ;
  myvector.SetPz (gen.pz ()) ;
  myvector.SetE (gen.energy ()) ;
}


std::pair<caloJetItr,caloJetItr>	
findTagJets (caloJetItr begin, caloJetItr end,
             double jetPtMin, double jetEtaMax) 
{

  std::pair<caloJetItr,caloJetItr> tagJets (begin,begin) ;
  double maxInvMass = 0. ;

  //PG find the tagging jets

  //PG first loop over jets
  for (caloJetItr firstJet = begin ; 
       firstJet != end ; 
       ++firstJet ) 
    {
      if (firstJet->pt () < jetPtMin || 
          fabs (firstJet->eta ()) > jetEtaMax) continue ;
      math::XYZTLorentzVector firstLV = firstJet->p4 () ;
      //PG second loop over jets
      for (caloJetItr secondJet = firstJet + 1 ; 
           secondJet != end ; 
           ++secondJet ) 
        {
          if (secondJet->pt () < jetPtMin || 
              fabs (secondJet->eta ()) > jetEtaMax) continue ;
	  //if (firstJet->eta ()*secondJet->eta () > 0) continue ; // only tags in opposite emispheres are considered
          math::XYZTLorentzVector sumLV = secondJet->p4 () + firstLV ;
          if (sumLV.M () > maxInvMass)
            {
              maxInvMass = sumLV.M () ;
              tagJets.first = firstJet ;
              tagJets.second = secondJet ;
            }
        } //PG second loop over jets
    } //PG first loop over jets

  return tagJets ;
}


// --------------------------------------------------------------------


std::pair<caloJetItr,caloJetItr>
findMaxPtJetsPair (caloJetItr begin, caloJetItr end,
             double jetPtMin, double jetEtaMax)
{
  std::pair<caloJetItr,caloJetItr> tagJets (begin, begin) ; 

  double ptMax1 = 0;
  double ptMax2 = 0;

  caloJetItr myJet1;
  caloJetItr myJet2;
  for (caloJetItr Jet = begin ;
       Jet != end ;
       ++Jet)
    {
      if (Jet->p4().Pt() > ptMax1)
	{
	  myJet1 = Jet;
	  myJet2 = myJet1;
	  ptMax2 = ptMax1;
	  ptMax1 = Jet->p4().Pt() ;

        } 
      // else if ((Jet->p4().Pt() > ptMax2) && (myJet1->p4().Eta() *Jet->p4().Eta() < 0))
      else if (Jet->p4().Pt() > ptMax2) 
	{
         myJet2 = Jet;
         ptMax2 = Jet->p4().Pt() ;
	}
    }

  tagJets.first = myJet1;
  tagJets.second = myJet2 ;

  return tagJets ;

    }


// --------------------------------------------------------------------------------


double 
deltaPhi (double phi1, double phi2)
{

  double deltaphi=fabs(phi1-phi2);
  if (deltaphi > 6.283185308) deltaphi -= 6.283185308;
  if (deltaphi > 3.141592654) deltaphi = 6.283185308-deltaphi;
  return deltaphi;
}


}


#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/HistoComposite.h"

#include <iostream>

using namespace std;
using namespace reco;
using namespace vbfhzz2l2b;

HistoComposite::
HistoComposite( std::string dir, std::string candTitle, std::string candName,
		double pt1, double pt2, double m1, double m2, TFileDirectory * parentDir )
  :
  HistoGroup<reco::CompositeCandidate>( dir, candTitle, candName, pt1, pt2, m1, m2, parentDir ),
  candName_(candName)
{
}

HistoComposite::~HistoComposite()
{
  if ( histoMuon_.map.size() > 0 ) {
    HistoMap<HistoMuon>::map_type::iterator begin = histoMuon_.map.begin();
    HistoMap<HistoMuon>::map_type::iterator end = histoMuon_.map.end();
    for ( HistoMap<HistoMuon>::map_type::iterator i = begin;
	  i != end; ++i ) if ( i->second ) delete i->second ;
  }

  if ( histoElectron_.map.size() > 0 ) {
    HistoMap<HistoElectron>::map_type::iterator begin = histoElectron_.map.begin();
    HistoMap<HistoElectron>::map_type::iterator end = histoElectron_.map.end();
    for ( HistoMap<HistoElectron>::map_type::iterator i = begin;
	  i != end; ++i ) if ( i->second ) delete i->second ;
  }

  if ( histoJet_.map.size() > 0 ) {
    HistoMap<HistoJet>::map_type::iterator begin = histoJet_.map.begin();
    HistoMap<HistoJet>::map_type::iterator end = histoJet_.map.end();
    for ( HistoMap<HistoJet>::map_type::iterator i = begin;
	  i != end; ++i ) if ( i->second ) delete i->second ;
  }

  if ( histoMET_.map.size() > 0 ) {
    HistoMap<HistoMET>::map_type::iterator begin = histoMET_.map.begin();
    HistoMap<HistoMET>::map_type::iterator end = histoMET_.map.end();
    for ( HistoMap<HistoMET>::map_type::iterator i = begin;
	  i != end; ++i ) if ( i->second ) delete i->second ;
  }

  if ( histoComposite_.map.size() > 0 ) {
    HistoMap<HistoComposite>::map_type::iterator begin = histoComposite_.map.begin();
    HistoMap<HistoComposite>::map_type::iterator end = histoComposite_.map.end();
    for ( HistoMap<HistoComposite>::map_type::iterator i = begin;
	  i != end; ++i ) if ( i->second ) delete i->second ;
  }
}

void HistoComposite::fill( const reco::CompositeCandidate * cand, double weight )
{


  // Fill 4-vector information for candidate
  HistoGroup<reco::CompositeCandidate>::fill( cand, 1, weight );

  const vector<string> & roles = cand->roles();

  if ( roles.size() != cand->numberOfDaughters() ) {
    cout << "HistoComposite::fill: Error: Nroles should match Ndaughters" << endl;
    return;
  }


  // Now fill information for daughters
  for (unsigned int i = 0; i < cand->numberOfDaughters(); ++i ) {
//      cout << "-------------processing component " << i << endl;
    const reco::Candidate * c = cand->daughter(i);
    string role = roles[i];

//     cout << "Role = " << roles[i] << endl;
//     cout << "pdgid = " << c->pdgId() << endl;
//     cout << "pt = " << c->pt() << endl;

    // Figure out what the candidate is based on type
    const Muon                     * pcmuon      = dynamic_cast<const Muon*>                    ( c );
    const GsfElectron              * pcelectron  = dynamic_cast<const GsfElectron*>             ( c );
    const CaloJet                  * pcjet       = dynamic_cast<const CaloJet*>                 ( c );
    const CaloMET                  * pcmet       = dynamic_cast<const CaloMET*>                 ( c );
    const reco::CompositeCandidate * pccomposite = dynamic_cast<const reco::CompositeCandidate*>( c );

    // The pointers might be in shallow clones, so check for that too
    const reco::ShallowClonePtrCandidate * pshallow = dynamic_cast<const reco::ShallowClonePtrCandidate *>(c);

    if ( pcmuon == 0 && c->hasMasterClonePtr() )  pcmuon = dynamic_cast<const Muon*>( &*(c->masterClonePtr()) );
    if ( pcelectron == 0 && c->hasMasterClonePtr() )  pcelectron = dynamic_cast<const GsfElectron*>( &*(c->masterClonePtr()) );
    if ( pcjet == 0 && c->hasMasterClonePtr() )  pcjet = dynamic_cast<const CaloJet*>( &*(c->masterClonePtr()) );
    if ( pcmet == 0 && c->hasMasterClonePtr() )  pcmet = dynamic_cast<const CaloMET*>( &*(c->masterClonePtr()) );

    if ( pccomposite == 0 && c->hasMasterClonePtr() )  pccomposite = dynamic_cast<const reco::CompositeCandidate*>( &*(c->masterClonePtr()) );

    // ------------------------------------------------------
    // Fill histograms if the candidate is a muon
    // ------------------------------------------------------
    if      ( pcmuon != 0 ) {
//         cout << "Filling muon" << endl;
       // Here is where we do not yet have a histogram for this role
      if ( histoMuon_.map.find( role ) == histoMuon_.map.end() ) {
	histoMuon_.map[role] = new HistoMuon( role, role, role,
					      pt1_, pt2_, m1_, m2_,
					      currDir_ ) ;
      }
      // Here is if the candidate is a shallow clone, we need to
      // fill kinematic information from the shallow clone and 
      // detector information from the base object
      if ( c->hasMasterClonePtr() ) {
	histoMuon_.map[role]    ->fill( pshallow, 1, weight );
      }
      // Here is if the candidate is a straightforward pointer
      else {
	histoMuon_.map[role]    ->fill( pcmuon, 1, weight );
      }
    }

    // ------------------------------------------------------
    // Fill histograms if the candidate is a electron
    // ------------------------------------------------------
    if      ( pcelectron != 0 ) {
//         cout << "Filling electron" << endl;
       // Here is where we do not yet have a histogram for this role
      if ( histoElectron_.map.find( role ) == histoElectron_.map.end() ) {
	histoElectron_.map[role] = new HistoElectron( role, role, role,
					      pt1_, pt2_, m1_, m2_,
					      currDir_  ) ;
      }
      // Here is if the candidate is a shallow clone, we need to
      // fill kinematic information from the shallow clone and 
      // detector information from the base object
      if ( c->hasMasterClonePtr() ) {
	histoElectron_.map[role]    ->fill( pshallow, 1, weight );
      }
      // Here is if the candidate is a straightforward pointer
      else {
	histoElectron_.map[role]    ->fill( pcelectron, 1, weight );
      }
    }


    // ------------------------------------------------------
    // Fill histograms if the candidate is a jet
    // ------------------------------------------------------
    if      ( pcjet != 0 ) {
//         cout << "Filling jet" << endl;
       // Here is where we do not yet have a histogram for this role
      if ( histoJet_.map.find( role ) == histoJet_.map.end() ) {
	histoJet_.map[role] = new HistoJet( role, role, role,
					      pt1_, pt2_, m1_, m2_,
					      currDir_  ) ;
      }
      // Here is if the candidate is a shallow clone, we need to
      // fill kinematic information from the shallow clone and 
      // detector information from the base object
      if ( c->hasMasterClonePtr() ) {
	histoJet_.map[role]    ->fill( pshallow, 1, weight );
      }
      // Here is if the candidate is a straightforward pointer
      else {
	histoJet_.map[role]    ->fill( pcjet, 1, weight );
      }
    }


    // ------------------------------------------------------
    // Fill histograms if the candidate is a met
    // ------------------------------------------------------
    if      ( pcmet != 0 ) {
//         cout << "Filling met" << endl;
       // Here is where we do not yet have a histogram for this role
      if ( histoMET_.map.find( role ) == histoMET_.map.end() ) {
	histoMET_.map[role] = new HistoMET( role, role, role,
					      pt1_, pt2_, m1_, m2_,
					      currDir_  ) ;
      }
      // Here is if the candidate is a shallow clone, we need to
      // fill kinematic information from the shallow clone and 
      // detector information from the base object
      if ( c->hasMasterClonePtr() ) {
	histoMET_.map[role]    ->fill( pshallow, 1, weight );
      }
      // Here is if the candidate is a straightforward pointer
      else {
	histoMET_.map[role]    ->fill( pcmet, 1, weight );
      }
    }



    // ------------------------------------------------------
    // Fill histograms if the candidate is a composite
    // ------------------------------------------------------
    if      ( pccomposite != 0 ) {
//       cout << "Filling composite with role " << role << endl;
       // Here is where we do not yet have a histogram for this role
      if ( histoComposite_.map.find( role ) == histoComposite_.map.end() ) {
	histoComposite_.map[role] = new HistoComposite( role, role, role,
							pt1_, pt2_, m1_, m2_,
							currDir_ ) ;
      }
      // Here is if the candidate is a shallow clone, we need to
      // fill kinematic information from the shallow clone and 
      // detector information from the base object
      if ( c->hasMasterClonePtr() ) {
	histoComposite_.map[role]    ->fill( pshallow, weight );
      }
      // Here is if the candidate is a straightforward pointer
      else {
	histoComposite_.map[role]    ->fill( pccomposite, weight );
      }
    }


  }
}




void HistoComposite::fill( const reco::ShallowClonePtrCandidate * pshallow, double weight )
{

  const reco::CompositeCandidate * cand = dynamic_cast<reco::CompositeCandidate const *>(pshallow);

  if ( cand == 0 ) {
    cout << "Error! Was passed a shallow clone that is not at heart a composite candidate" << endl;
    return;
  }

  // Fill 4-vector information for candidate
  HistoGroup<reco::CompositeCandidate>::fill( pshallow, 1, weight );

  const vector<string> & roles = cand->roles();

  if ( roles.size() != cand->numberOfDaughters() ) {
    cout << "HistoComposite::fill: Error: Nroles should match Ndaughters" << endl;
    return;
  }


  // Now fill information for daughters
  for (unsigned int i = 0; i < cand->numberOfDaughters(); ++i ) {
//      cout << "-------------processing component " << i << endl;
    const reco::Candidate * c = cand->daughter(i);
    string role = roles[i];

//     cout << "Role = " << roles[i] << endl;
//     cout << "pdgid = " << c->pdgId() << endl;
//     cout << "pt = " << c->pt() << endl;

    // Figure out what the candidate is based on type
    const Muon                     * pcmuon      = dynamic_cast<const Muon*>                    ( c );
    const GsfElectron              * pcelectron  = dynamic_cast<const GsfElectron*>             ( c );
    const CaloJet                  * pcjet       = dynamic_cast<const CaloJet*>                 ( c );
    const CaloMET                  * pcmet       = dynamic_cast<const CaloMET*>                 ( c );
    const reco::CompositeCandidate * pccomposite = dynamic_cast<const reco::CompositeCandidate*>( c );

    // The pointers might be in shallow clones, so check for that too
    const reco::ShallowClonePtrCandidate * pshallow_da = dynamic_cast<const reco::ShallowClonePtrCandidate *>(c);

    if ( pcmuon == 0 && c->hasMasterClonePtr() )  pcmuon = dynamic_cast<const Muon*>( &*(c->masterClonePtr()) );
    if ( pcelectron == 0 && c->hasMasterClonePtr() )  pcelectron = dynamic_cast<const GsfElectron*>( &*(c->masterClonePtr()) );

    if ( pcjet == 0 && c->hasMasterClonePtr() )  pcjet = dynamic_cast<const CaloJet*>( &*(c->masterClonePtr()) );
    if ( pcmet == 0 && c->hasMasterClonePtr() )  pcmet = dynamic_cast<const CaloMET*>( &*(c->masterClonePtr()) );
    if ( pccomposite == 0 && c->hasMasterClonePtr() )  pccomposite = dynamic_cast<const reco::CompositeCandidate*>( &*(c->masterClonePtr()) );

    // ------------------------------------------------------
    // Fill histograms if the candidate is a muon
    // ------------------------------------------------------
    if      ( pcmuon != 0 ) {
//         cout << "Filling muon" << endl;
       // Here is where we do not yet have a histogram for this role
      if ( histoMuon_.map.find( role ) == histoMuon_.map.end() ) {
	histoMuon_.map[role] = new HistoMuon( role, role, role,
					      pt1_, pt2_, m1_, m2_,
					      currDir_ ) ;
      }
      // Here is if the candidate is a shallow clone, we need to
      // fill kinematic information from the shallow clone and 
      // detector information from the base object
      if ( c->hasMasterClonePtr() ) {
	histoMuon_.map[role]    ->fill( pshallow_da, 1, weight );
      }
      // Here is if the candidate is a straightforward pointer
      else {
	histoMuon_.map[role]    ->fill( pcmuon, 1, weight );
      }
    }

    // ------------------------------------------------------
    // Fill histograms if the candidate is a electron
    // ------------------------------------------------------
    if      ( pcelectron != 0 ) {
//         cout << "Filling electron" << endl;
       // Here is where we do not yet have a histogram for this role
      if ( histoElectron_.map.find( role ) == histoElectron_.map.end() ) {
	histoElectron_.map[role] = new HistoElectron( role, role, role,
					      pt1_, pt2_, m1_, m2_,
					      currDir_  ) ;
      }
      // Here is if the candidate is a shallow clone, we need to
      // fill kinematic information from the shallow clone and 
      // detector information from the base object
      if ( c->hasMasterClonePtr() ) {
	histoElectron_.map[role]    ->fill( pshallow_da, 1, weight );
      }
      // Here is if the candidate is a straightforward pointer
      else {
	histoElectron_.map[role]    ->fill( pcelectron, 1, weight );
      }
    }


    // ------------------------------------------------------
    // Fill histograms if the candidate is a jet
    // ------------------------------------------------------
    if      ( pcjet != 0 ) {
//         cout << "Filling jet" << endl;
       // Here is where we do not yet have a histogram for this role
      if ( histoJet_.map.find( role ) == histoJet_.map.end() ) {
	histoJet_.map[role] = new HistoJet( role, role, role,
					      pt1_, pt2_, m1_, m2_,
					      currDir_  ) ;
      }
      // Here is if the candidate is a shallow clone, we need to
      // fill kinematic information from the shallow clone and 
      // detector information from the base object
      if ( c->hasMasterClonePtr() ) {
	histoJet_.map[role]    ->fill( pshallow_da, 1, weight );
      }
      // Here is if the candidate is a straightforward pointer
      else {
	histoJet_.map[role]    ->fill( pcjet, 1, weight );
      }
    }


    // ------------------------------------------------------
    // Fill histograms if the candidate is a met
    // ------------------------------------------------------
    if      ( pcmet != 0 ) {
//         cout << "Filling met" << endl;
       // Here is where we do not yet have a histogram for this role
      if ( histoMET_.map.find( role ) == histoMET_.map.end() ) {
	histoMET_.map[role] = new HistoMET( role, role, role,
					      pt1_, pt2_, m1_, m2_,
					      currDir_  ) ;
      }
      // Here is if the candidate is a shallow clone, we need to
      // fill kinematic information from the shallow clone and 
      // detector information from the base object
      if ( c->hasMasterClonePtr() ) {
	histoMET_.map[role]    ->fill( pshallow_da, 1, weight );
      }
      // Here is if the candidate is a straightforward pointer
      else {
	histoMET_.map[role]    ->fill( pcmet, 1, weight );
      }
    }



    // ------------------------------------------------------
    // Fill histograms if the candidate is a composite
    // ------------------------------------------------------
    if      ( pccomposite != 0 ) {
//       cout << "Filling composite with role " << role << endl;
       // Here is where we do not yet have a histogram for this role
      if ( histoComposite_.map.find( role ) == histoComposite_.map.end() ) {
	histoComposite_.map[role] = new HistoComposite( role, role, role,
							pt1_, pt2_, m1_, m2_,
							currDir_ ) ;
      }
      // Here is if the candidate is a shallow clone, we need to
      // fill kinematic information from the shallow clone and 
      // detector information from the base object
      if ( c->hasMasterClonePtr() ) {
	histoComposite_.map[role]    ->fill( pshallow_da, weight );
      }
      // Here is if the candidate is a straightforward pointer
      else {
	histoComposite_.map[role]    ->fill( pccomposite, weight );
      }
    }


  }
}

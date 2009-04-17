#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/HistoJet.h"


#include <iostream>

using vbfhzz2l2b::HistoJet;
using namespace std;

// Constructor:

HistoJet::HistoJet(  std::string dir, std::string group,std::string pre,
		     double pt1, double pt2, double m1, double m2,
		     TFileDirectory * parentDir ) 
  : HistoGroup<CaloJet>( dir, group, pre, pt1, pt2, m1, m2, parentDir)
{
  // book relevant jet histograms

  addHisto( mMultiplicity_         =
            new PhysVarHisto( pre + "Multiplicity", "Multiplicity", 100, 0, 100, currDir_, "", "vD" )
            );

  addHisto( mEta_         =
            new PhysVarHisto( pre + "Eta", "Eta",100, 0, 100, currDir_, "", "vD" )
            );


  addHisto( mPhi_         =
            new PhysVarHisto( pre + "Phi", "Phi", 100, 0, 100,currDir_, "", "vD" )
            );


  addHisto( mE_         =
            new PhysVarHisto( pre + "E",   "E", 100, 0, 100,currDir_, "", "vD" )
            );
  addHisto( mP_         =
            new PhysVarHisto( pre + "P",  "P", 100, 0, 100,currDir_, "", "vD" )
            );

  addHisto( mPt_         =
            new PhysVarHisto( pre + "Pt",  "Pt", 100, 0, 100,currDir_, "", "vD" )
            );

  addHisto( mPt_1_         =
            new PhysVarHisto( pre + "Pt1", "Pt1", 100, 0, 100,currDir_, "", "vD" )
            );
  addHisto( mPt_2_         =
            new PhysVarHisto( pre + "Pt2", "Pt2", 100, 0, 300,currDir_, "", "vD" )
            );
  addHisto( mPt_3_         =
            new PhysVarHisto( pre + "Pt3", "Pt3", 100, 0, 5000,currDir_, "", "vD" )
            );

  addHisto( mMass_         =
            new PhysVarHisto( pre + "Mass", "Mass", 100, 0, 25,currDir_, "", "vD" )
            );

  //  addHisto( mNTracks_         =
  //	    new PhysVarHisto( pre + "NTracks", "Jet N_{TRK}", 51, -0.5, 50.5, currDir_, "", "vD" )
  //	    );

  addHisto( mConstituents_         =
            new PhysVarHisto( pre + "Constituents", "# of Constituents", 100, 0, 30,currDir_, "", "vD" )
            );
 // CaloJet specific
  addHisto( mHadEnergyInHO_         =
            new PhysVarHisto( pre + "HadEnergyInHO", "HadEnergyInHO", 100, 0, 10, currDir_, "", "vD" )
            );

  addHisto( mHadEnergyInHB_         =
            new PhysVarHisto( pre + "HadEnergyInHB", "HadEnergyInHB", 100, 0, 50, currDir_, "", "vD" )
            );
  addHisto( mHadEnergyInHF_         =
            new PhysVarHisto( pre + "HadEnergyInHF", "HadEnergyInHF", 100, 0, 50, currDir_, "", "vD" )
            );
  addHisto( mHadEnergyInHE_         =
            new PhysVarHisto( pre + "HadEnergyInHE", "HadEnergyInHE", 100, 0, 100, currDir_, "", "vD" )
            );
  addHisto( mEmEnergyInEB_         =
            new PhysVarHisto( pre + "EmEnergyInEB", "EmEnergyInEB", 100, 0, 10, currDir_, "", "vD" )
            );
  addHisto( mEmEnergyInEE_         =
            new PhysVarHisto( pre + "EmEnergyInEE", "EmEnergyInEE", 100, 0, 50, currDir_, "", "vD" )
            );
  addHisto( mEmEnergyInHF_         =
            new PhysVarHisto( pre + "EmEnergyInHF", "EmEnergyInHF", 120, -20, 100, currDir_, "", "vD" )
            );
  addHisto( mHadEnergyFraction_         =
            new PhysVarHisto( pre + "HadEnergyFraction", "HAD Energy Fraction", 120, -0.1, 1.1, currDir_, "", "vD" )
            );
  addHisto( mEmEnergyFraction_         =
            new PhysVarHisto( pre + "EmEnergyFraction", "EM Energy Fraction", 120, -0.1, 1.1, currDir_, "", "vD" )
            );
  //  addHisto( mChargeFraction_         =
  //            new PhysVarHisto( pre + "ChargeFraction", "Fraction of charged tracks pt",500,0,5, currDir_, "", "vD" )
  //            );

}

HistoJet::~HistoJet()
{
}


void HistoJet::fill( const CaloJet * jet, uint iJet, double weight )
{

  // First fill common 4-vector histograms
  HistoGroup<CaloJet>::fill( jet, iJet, weight );

  // fill relevant jet histograms
  jetME_->fill(1,   iJet,   weight);
  if (mEta_)          mEta_->fill (        jet->eta(),                     iJet,   weight);
  if (mPhi_)          mPhi_->fill (        jet->phi(),                     iJet,   weight);
  if (mE_)            mE_->fill   (        jet->energy(),                  iJet,   weight);
  if (mP_)            mP_->fill   (        jet->p(),                       iJet,   weight);
  if (mPt_)           mPt_->fill  (        jet->pt(),                      iJet,   weight);
  if (mPt_1_)         mPt_1_->fill(        jet->pt(),                      iJet,   weight);
  if (mPt_2_)         mPt_2_->fill(        jet->pt(),                      iJet,   weight);
  if (mPt_3_)         mPt_3_->fill(        jet->pt(),                      iJet,   weight);
  if (mMass_)         mMass_->fill(        jet->mass(),                    iJet,   weight);
  //  if (mNTracks_)      mNTracks_->fill(     jet->associatedTracks().size(), iJet, weight );
  if (mConstituents_) mConstituents_->fill(jet->nConstituents(),           iJet,   weight);

  if (mHadEnergyInHO_)     mHadEnergyInHO_->fill(         jet->hadEnergyInHO(),          iJet,   weight);
  if (mHadEnergyInHB_)     mHadEnergyInHB_->fill(         jet->hadEnergyInHB(),          iJet,   weight);
  if (mHadEnergyInHF_)     mHadEnergyInHF_->fill(         jet->hadEnergyInHF(),          iJet,   weight);
  if (mHadEnergyInHE_)     mHadEnergyInHE_->fill(         jet->hadEnergyInHE(),          iJet,   weight);
  if (mEmEnergyInEB_)      mEmEnergyInEB_->fill(          jet->emEnergyInEB(),           iJet,   weight);
  if (mEmEnergyInEE_)      mEmEnergyInEE_->fill(          jet->emEnergyInEE(),           iJet,   weight);
  if (mEmEnergyInHF_)      mEmEnergyInHF_->fill(          jet->emEnergyInHF(),           iJet,   weight);
  if (mHadEnergyFraction_) mHadEnergyFraction_->fill(     jet->energyFractionHadronic(), iJet,   weight);
  if (mEmEnergyFraction_)  mEmEnergyFraction_->fill(      jet->emEnergyFraction(),       iJet,   weight);
  //  if (mChargeFraction_)   mChargeFraction_->fill( (JetTracksAssociation::tracksP4(*jetTracks,*i_caljet)).pt()/pt, iJet, weight);

}

void HistoJet::fill( const reco::ShallowClonePtrCandidate * pshallow, uint iJet, double weight )
{

  // Get the underlying object that the shallow clone represents
  const CaloJet * jet = dynamic_cast<const CaloJet*>(&*(pshallow->masterClonePtr()));

  if ( jet == 0 ) {
    cout << "Error! Was passed a shallow clone that is not at heart a jet" << endl;
    return;
  }

  // First fill common 4-vector histograms from shallow clone
  HistoGroup<CaloJet>::fill( pshallow, iJet, weight);

  // fill relevant jet histograms
  jetME_->fill(1,   iJet,   weight);

  if (mEta_)          mEta_->fill(         jet->eta(),                     iJet,   weight);
  if (mPhi_)          mPhi_->fill(         jet->phi(),                     iJet,   weight);
  if (mE_)            mE_->fill(           jet->energy(),                  iJet,   weight);
  if (mP_)            mP_->fill(           jet->p(),                       iJet,   weight);
  if (mPt_)           mPt_->fill(          jet->pt(),                      iJet,   weight);
  if (mPt_1_)         mPt_1_->fill(        jet->pt(),                      iJet,   weight);
  if (mPt_2_)         mPt_2_->fill(        jet->pt(),                      iJet,   weight);
  if (mPt_3_)         mPt_3_->fill(        jet->pt(),                      iJet,   weight);
  if (mMass_)         mMass_->fill(        jet->mass(),                    iJet,   weight);
  //  if (mNTracks_)      mNTracks_->fill(     jet->associatedTracks().size(), iJet, weight );
  if (mConstituents_) mConstituents_->fill(jet->nConstituents(),           iJet,   weight);

}


void HistoJet::fillCollection( const std::vector<CaloJet> & coll, double weight ) 
{
 
  h_size_->fill( coll.size(), 1, weight );     //! Save the size of the collection.

  if (mMultiplicity_) mMultiplicity_->fill( coll.size(), 1, weight);

  std::vector<CaloJet>::const_iterator
    iobj = coll.begin(),
    iend = coll.end();

  uint i = 1;              //! Fortran-style indexing
  for ( ; iobj != iend; ++iobj, ++i ) {
    fill( &*iobj, i, weight);      //! &*iobj dereferences to the pointer to a PHYS_OBJ*
  } 
}

void HistoJet::clearVec()
{
  HistoGroup<CaloJet>::clearVec();

  jetME_->clearVec( );

  // Generic Jet Parameters
  mMultiplicity_->clearVec( );
  mEta_->clearVec( );
  mPhi_->clearVec( );
  mE_->clearVec( );
  mP_->clearVec( );
  mPt_->clearVec( );
  mPt_1_->clearVec( );
  mPt_2_->clearVec( );
  mPt_3_->clearVec( );
  mMass_->clearVec( );
  //  mNTracks_->clearVec( );
  mConstituents_->clearVec( );

}

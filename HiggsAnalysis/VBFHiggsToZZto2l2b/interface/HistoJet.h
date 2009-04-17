#ifndef HiggsAnalysis_VBFHiggsToZZto2l2b_HistoJet_h
#define HiggsAnalysis_VBFHiggsToZZto2l2b_HistoJet_h

//------------------------------------------------------------
// Title: HistoJet.h
// Purpose: To histogram Jets
//------------------------------------------------------------
//
// Interface:
//
//   HistoJet ( TFile * file );
//   Description: Constructor.
//
//   void fill( TK::Jet * );
//   Description: Fill object. Will fill relevant jet variables
//
//   void write();
//   Description: Write object to file in question.
//
//   ~HistoJet
//    Description: Destructor. Deallocates memory.
//
//------------------------------------------------------------
//
// Modification History:
//
//------------------------------------------------------------


// CMSSW include files
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/HistoGroup.h"

// STL include files
#include <string>

// ROOT include files
#include <TH1D.h>
#include <TFile.h>

namespace vbfhzz2l2b {

  using namespace reco;

  class HistoJet : public HistoGroup<CaloJet> {

  public:
    HistoJet( std::string dir = "jet",std::string group = "Jet",
	      std::string pre="jet",
	      double pt1=0, double pt2=200, double m1=0, double m2=200,
	      TFileDirectory * parentDir=0 );
    virtual ~HistoJet();


    // fill a plain ol' jet:
    virtual void fill( const CaloJet *jet, uint iPart = 1, double weight = 1.0 );
    virtual void fill( const CaloJet &jet, uint iPart = 1, double weight = 1.0 ) { fill(&jet, iPart,weight); }

    // fill a jet that is a shallow clone, and take kinematics from 
    // shallow clone but detector plots from the jet itself
    virtual void fill( const reco::ShallowClonePtrCandidate *jet, uint iPart = 1, double weight = 1.0 );
    virtual void fill( const reco::ShallowClonePtrCandidate &jet, uint iPart = 1, double weight = 1.0 )
    { fill(&jet, iPart, weight); }

    virtual void fillCollection( const std::vector<CaloJet> & coll, double weight = 1.0 );

    // Clear ntuple cache
    void clearVec();
  protected:
    PhysVarHisto *    jetME_;

   // Generic Jet Parameters
    PhysVarHisto *    mMultiplicity_;
    PhysVarHisto *    mEta_;
    PhysVarHisto *    mPhi_;
    PhysVarHisto *    mE_;
    PhysVarHisto *    mP_;
    PhysVarHisto *    mPt_;
    PhysVarHisto *    mPt_1_;
    PhysVarHisto *    mPt_2_;
    PhysVarHisto *    mPt_3_;
    PhysVarHisto *    mMass_;
    //    PhysVarHisto *    mNTracks_;
    PhysVarHisto *    mConstituents_;

 // CaloJet specific
    PhysVarHisto *    mHadEnergyInHO_;
    PhysVarHisto *    mHadEnergyInHB_;
    PhysVarHisto *    mHadEnergyInHF_;
    PhysVarHisto *    mHadEnergyInHE_;
    PhysVarHisto *    mEmEnergyInEB_;
    PhysVarHisto *    mEmEnergyInEE_;
    PhysVarHisto *    mEmEnergyInHF_;
    PhysVarHisto *    mHadEnergyFraction_;
    PhysVarHisto *    mEmEnergyFraction_;
    PhysVarHisto *    mChargeFraction_;

//  PhysVarHisto * mEBfractionVsEta_;
//  PhysVarHisto * mEEfractionVsEta_;
//  PhysVarHisto * mHBfractionVsEta_;
//  PhysVarHisto * mHOfractionVsEta_;
//  PhysVarHisto * mHEfractionVsEta_;
//  PhysVarHisto * mHFfractionVsEta_; 
//  PhysVarHisto * mCaloEnergyVsEta_;
//  PhysVarHisto * memEnergyVsEta_;
//  PhysVarHisto * mhadEnergyVsEta_;

// generator distributions
    PhysVarHisto * mGenJetMulti_;
    PhysVarHisto * mPtHat_;
    PhysVarHisto * mGenPt_;
    PhysVarHisto * mGenEta_;
    PhysVarHisto * mGenPhi_;
    PhysVarHisto * mGenInvMassLeading_;
    PhysVarHisto * mdR_;
    
//    PhysVarHisto * mGenEnergyVsEta_;
//    PhysVarHisto * mrespVsPtBarrel_;
//    PhysVarHisto * mCaloErespVsEta_;
//    PhysVarHisto * memErespVsEta_;
//    PhysVarHisto * mhadErespVsEta_;


  };

}
#endif

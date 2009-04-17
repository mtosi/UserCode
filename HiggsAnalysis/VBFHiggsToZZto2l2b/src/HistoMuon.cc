#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/HistoMuon.h"

#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"



#include <iostream>
#include <sstream>

using vbfhzz2l2b::HistoMuon;
using namespace std;

// Constructor:
HistoMuon::HistoMuon(std::string dir, std::string group,std::string pre,
		     double pt1, double pt2, double m1, double m2,
		     TFileDirectory * parentDir)
  : HistoGroup<Muon>( dir, group, pre, pt1, pt2, m1, m2, parentDir)
{
  std::string histname = "";
  histname = "TrackIso";

  histname = "CaloE";
  addHisto( h_caloE_ =
            new PhysVarHisto( pre + histname, "Muon Calorimeter Energy", 50, 0, 50, currDir_, "", "vD" )
            );

  //muon energy deposit analyzer
  histname = "ecalDepositedEnergy";
  addHisto( ecalDepEnergy_ =
            new PhysVarHisto( pre + histname, histname, 50, 0., 20., currDir_, "", "vD" )
            );

  histname = "ecalS9DepositedEnergy";
  addHisto( ecalS9DepEnergy_ =
            new PhysVarHisto( pre + histname, histname, 80, 0. ,40., currDir_, "", "vD" )
            );
  histname = "hadDepositedEnergy";
  addHisto( hcalDepEnergy_ =
            new PhysVarHisto( pre + histname, histname, 65, -0.5, 64.5, currDir_, "", "vD" )
            );

  histname = "hadS9DepositedEnergy";
  addHisto( hcalS9DepEnergy_ =
            new PhysVarHisto( pre + histname, histname, 80, 0. ,40., currDir_, "", "vD" )
            );

  histname = "hoDepositedEnergy";
  addHisto( hoDepEnergy_ =
            new PhysVarHisto( pre + histname, histname, 50, 0. ,20. , currDir_, "", "vD" )
            );

  histname = "hoS9DepositedEnergy";
  addHisto( hoS9DepEnergy_ =
            new PhysVarHisto( pre + histname, histname, 50, 0. ,20. , currDir_, "", "vD" )
            );
}


// fill a plain ol' muon
void HistoMuon::fill( const Muon *muon, uint iMu, double weight )
{

  // First fill common 4-vector histograms

  HistoGroup<Muon>::fill( muon, iMu, weight);

  // fill relevant muon histograms
  reco::MuonEnergy muEnergy = muon->calEnergy();
  h_caloE_->fill( muEnergy.em+muEnergy.had+muEnergy.ho, iMu , weight);

///////////////////////////////////////

  // get all the mu energy deposits
  ecalDepEnergy_->fill(muEnergy.em, iMu, weight);

  hcalDepEnergy_->fill(muEnergy.had, iMu, weight);

  hoDepEnergy_->fill(muEnergy.ho, iMu, weight);

  ecalS9DepEnergy_->fill(muEnergy.emS9, iMu, weight);

  hcalS9DepEnergy_->fill(muEnergy.hadS9, iMu, weight);

  hoS9DepEnergy_->fill(muEnergy.hoS9, iMu, weight);

}


// fill a muon that is a shallow clone, and take kinematics from 
// shallow clone but detector plots from the muon itself
void HistoMuon::fill( const reco::ShallowClonePtrCandidate *pshallow, uint iMu, double weight )
{

  // Get the underlying object that the shallow clone represents
  const Muon * muon = dynamic_cast<const Muon*>(&*(pshallow->masterClonePtr()));

  if ( muon == 0 ) {
    cout << "Error! Was passed a shallow clone that is not at heart a muon" << endl;
    return;
  }

  

  // First fill common 4-vector histograms from shallow clone

  HistoGroup<Muon>::fill( pshallow, iMu, weight);

  // fill relevant muon histograms from muon
  reco::MuonEnergy muEnergy = muon->calEnergy();

  h_caloE_->fill( muEnergy.em+muEnergy.had+muEnergy.ho, iMu , weight);

///////////////////////////////////////

  // get all the mu energy deposits
  ecalDepEnergy_->fill(muEnergy.em, iMu, weight);

  hcalDepEnergy_->fill(muEnergy.had, iMu, weight);

  hoDepEnergy_->fill(muEnergy.ho, iMu, weight);

  ecalS9DepEnergy_->fill(muEnergy.emS9, iMu, weight);

  hcalS9DepEnergy_->fill(muEnergy.hadS9, iMu, weight);

  hoS9DepEnergy_->fill(muEnergy.hoS9, iMu, weight);

}

void HistoMuon::fillCollection( const std::vector<Muon> & coll, double weight )
{

  h_size_->fill( coll.size(), 1, weight );     //! Save the size of the collection.

  std::vector<Muon>::const_iterator
    iobj = coll.begin(),
    iend = coll.end();

  uint i = 1;              //! Fortran-style indexing
  for ( ; iobj != iend; ++iobj, ++i ) {
    fill( &*iobj, i, weight);      //! &*iobj dereferences to the pointer to a PHYS_OBJ*
  }
}


void HistoMuon::clearVec()
{
  HistoGroup<Muon>::clearVec();

  h_caloE_->clearVec();
  
//muon energy deposit analyzer
  ecalDepEnergy_->clearVec();
  ecalS9DepEnergy_->clearVec();
  hcalDepEnergy_->clearVec();
  hcalS9DepEnergy_->clearVec();
  hoDepEnergy_->clearVec();
  hoS9DepEnergy_->clearVec();

}

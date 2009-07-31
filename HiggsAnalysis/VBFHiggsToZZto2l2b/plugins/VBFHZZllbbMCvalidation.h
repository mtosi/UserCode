// -*- C++ -*-
//
// Package:    VBFHZZllbbMCvalidation
// Class:      VBFHZZllbbMCvalidation
// 
/**\class VBFHZZllbbMCvalidation VBFHZZllbbMCvalidation.cc HiggsAnalysis/VBFHiggsToZZto2l2b/src/VBFHZZllbbMCvalidation.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Mia TOSI
//         Created:  Tue Jan 20 15:48:58 CET 2009
// $Id: VBFHZZllbbMCvalidation.h,v 1.3 2009/05/14 16:52:18 tosi Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/PythiaParticleIndex.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/ParticlesMass.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/ParticlesCharge.h"

// Root includes
// -------------
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TFile.h"

// For file output
// ---------------
#include <fstream>
#include <sstream>
#include <cmath>
#include <memory>
#include <string>
#include <iostream>
#include <iomanip>


// class decleration
//

class VBFHZZllbbMCvalidation : public edm::EDAnalyzer {
   public:
      explicit VBFHZZllbbMCvalidation(const edm::ParameterSet&);
      ~VBFHZZllbbMCvalidation();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

//
// constants, enums and typedefs
//

//
// static data member definitions
//

  // ----------member data ---------------------------
  unsigned int eventcounter_;
  unsigned int eventVBFcounter_;
  unsigned int eventZcounter_;
  unsigned int eventZZcounter_;
  unsigned int eventHadronicZcounter_;
  unsigned int eventLeptonicZcounter_;
  unsigned int eventZintoEcounter_;
  unsigned int eventZintoMUcounter_;
  unsigned int eventZintoTAUcounter_;
  unsigned int eventZintoBQUARKcounter_;
  unsigned int eventHeavyQcounter_;
  unsigned int eventLightQcounter_;

  edm::ParameterSet conf_;
  std::string MCParticleLabel_;
  int         signal_;

  TFile* OutputFile;

  TH1D * eventsNumber_;
  
  TH1D     * tagQuarkFlavour_;
  TProfile * tagQuarkFlavour_profile_;
  TH1D     * ZleptonFlavour_;
  TProfile * ZleptonFlavour_profile_;
  TH1D     * ZquarkFlavour_;
  TProfile * ZquarkFlavour_profile_;

  TH1D * Zrapidity_;
  TH1D * Zeta_;
  TH1D * Zpt_;
  TH1D * Zphi_;
  TH1D * Zmass_;
  TH1D * leptonicZrapidity_;
  TH1D * leptonicZeta_;
  TH1D * leptonicZpt_;
  TH1D * leptonicZphi_;
  TH1D * leptonicZmass_;
  TH1D * hadronicZrapidity_;
  TH1D * hadronicZeta_;
  TH1D * hadronicZpt_;
  TH1D * hadronicZphi_;
  TH1D * hadronicZmass_;

  TH1D * zQUARKrapidity_;
  TH1D * zQUARKeta_;
  TH1D * zQUARKpt_;
  TH1D * zQUARKphi_;

  TH1D * tagSYSTEMrapidity_;
  TH1D * tagSYSTEMeta_;
  TH1D * tagSYSTEMpt_;
  TH1D * tagSYSTEMphi_;
  TH1D * tagSYSTEMmass_;
  TH1D * tagQUARKrapidity_;
  TH1D * tagQUARKeta_;
  TH1D * tagQUARKpt_;
  TH1D * tagQUARKphi_;

  TH2D * tagQUARKetaVSphi_;
  TH2D * tagQUARKetaVSeta_;
  TH2D * zQUARKetaVSphi_;
  TH2D * zQUARKetaVSeta_;
};

//define this as a plug-in
//DEFINE_FWK_MODULE(VBFHZZllbbMCvalidation);

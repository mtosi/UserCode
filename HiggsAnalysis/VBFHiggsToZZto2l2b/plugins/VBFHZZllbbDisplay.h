// -*- C++ -*-
//
// Package:    VBFHZZllbbDisplay
// Class:      VBFHZZllbbDisplay
// 
//class VBFHZZllbbDisplay VBFHZZllbbDisplay.cc HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbDisplay.cc
//
//

#ifndef VBFHZZllbbDisplay_H
#define VBFHZZllbbDisplay_H

// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include "AnalysisExamples/AnalysisObjects/interface/PythiaParticleIndex.h"
#include "AnalysisExamples/AnalysisObjects/interface/ParticlesMass.h"
#include "AnalysisExamples/AnalysisObjects/interface/ParticlesCharge.h"


//root include
#include <TCanvas.h>
#include <TLegend.h>
#include <TH2.h>
#include <TTree.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TVector3.h>


class VBFHZZllbbDisplay : public edm::EDAnalyzer {

  public:
    
      explicit VBFHZZllbbDisplay (const edm::ParameterSet&);
      ~VBFHZZllbbDisplay ();

   private:

      virtual void beginJob (const edm::EventSetup&);
      virtual void analyze (const edm::Event&, const edm::EventSetup&);
      virtual void endJob ();

   private:

      edm::InputTag genParticleLabel_;
      edm::InputTag jetLabel_;
      edm::InputTag muonLabel_;
      edm::InputTag electronLabel_;
      edm::InputTag tagJetLabel_;
      edm::InputTag metLabel_;


      unsigned int eventcounter_;

      // common graphic objects
      TCanvas * canvas_ ;
      TH2F * histo2D_ ;
      TLegend * legend_ ;

      double * tagJet_eta_;
      double * tagJet_phi_;
      TGraph * tagJet_graph_;

      double * jet_eta_;
      double * jet_phi_;
      TGraph * jet_graph_;

      double * electron_eta_;
      double * electron_phi_;
      TGraph * electron_graph_;

      double * muon_eta_;
      double * muon_phi_;
      TGraph * muon_graph_;

      double * MET_eta_;
      double * MET_phi_;
      TGraph * MET_graph_;

      TH2F * MCparticle_;

  template <class handle>
  void fillgraph (double * eta, double * phi, TGraph * map,
		  handle & inputColl,
		  int marker, int color)
  {
    if (inputColl->size () == 0) return ;
    eta = new double [inputColl->size ()] ;
    phi = new double [inputColl->size ()] ;
    
    for (int it = 0 ; it < inputColl->size () ; ++it)
      {
	eta[it] = (*inputColl)[it].eta () ;
	phi[it] = (*inputColl)[it].phi () ;
           }
    map = new TGraph (inputColl->size (), eta, phi) ;
    map->SetMarkerStyle (marker) ;
    map->SetMarkerColor (color) ;
  }                     
};

#endif

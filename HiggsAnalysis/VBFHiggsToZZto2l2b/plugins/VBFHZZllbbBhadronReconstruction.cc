// -*- C++ -*-
//
// Package:    VBFHZZllbbBhadronReconstruction
// Class:      VBFHZZllbbBhadronReconstruction
// 
/**\class VBFHZZllbbBhadronReconstruction VBFHZZllbbBhadronReconstruction.cc HiggsAnalysis/VBFHiggsToZZto2l2b/src/VBFHZZllbbBhadronReconstruction.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Mia TOSI
//         Created:  Mon Feb  2 17:31:44 CET 2009
// $Id$
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

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbBhadronReconstruction.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/CorJetWithBTagDiscr.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace vbfhzz2l2b;

// class decleration
//

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
VBFHZZllbbBhadronReconstruction::VBFHZZllbbBhadronReconstruction(const edm::ParameterSet& iConfig) :
  corJetWithBTagLabel_( iConfig.getUntrackedParameter<std::string>("corJetWithBTagLabel") )
{
   //now do what ever initialization is needed
  MB_  = 5.28;       // B hadron mass (GeV/c^2)
  MB2_ = 5.28*5.28;  // B hadron mass (GeV/c^2)

}

VBFHZZllbbBhadronReconstruction::~VBFHZZllbbBhadronReconstruction()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
VBFHZZllbbBhadronReconstruction::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<vbfhzz2l2b::CorJetWithBTagDiscrCollection> corJetWithBTagHandle;
  //  edm::Handle<vbfhzz2l2b::CorJetWithBTagsDiscrCollection> corJetWithBTagHandle;
  iEvent.getByLabel(corJetWithBTagLabel_,"corJetWithBTagDiscr",corJetWithBTagHandle);
  std::cout << "corJetWithBTagHandle->size(): " << corJetWithBTagHandle->size() << std::endl;

  int jetIndex = 0;    
  std::vector<reco::JetBaseRef> jets = vbfhzz2l2b::CorJetBTagDiscrAssociation::allJets(*corJetWithBTagHandle);
  for ( std::vector<reco::JetBaseRef>::const_iterator jet = jets.begin(); 
	jet != jets.end();
	++jet, jetIndex++ ) {
    std::cout << "jetIndex: " << jetIndex << std::endl;
    
    double corrEt = (*corJetWithBTagHandle)[*jet].corEt_;
    double uncorrEt = (*jet)->et();
    std::cout << "uncorrEt: " << uncorrEt << std::endl;
    std::cout << "corrEt: " << corrEt << std::endl;
    double emFrac = (dynamic_cast<const reco::CaloJet*>(&**jet))->emEnergyFraction();
    std::cout << "emFrac: " << emFrac << std::endl;

  }

  // B hadron variables reconstruction
  // -------------------------------------------------------------------------------
  // WARNING:
  //  P0L seems to have a tiny negative tail, 
  //  probably owned to a not always right computation of cosTheta
  //  and so due to L3d rather than P3sumotrktag estimation
  //  we are expecting, indeed, a very small angle 
  //  between B hadron flight direction and tag jet axis,
  //  while a not perfect calculation of the primary and secondary vertex position 
  //  could explain this anomaly
  //
  //  However, _ptneu and _ptbhad are saved in SqzNtuple 
  // -------------------------------------------------------------------------------

  double P0L    = 0.;
  double ptBhad = 0.;

  // - secondary vertex information from "secondaryVertexBTagInfos"
  //    WARNING: available in the edm event only for IC5 jet collection
  // - primary vertex information from "offlinePrimaryVertices"
  // also available info about jet vertex
  TLorentzVector P4secVtx;

  const TVector3 secVtx = 0.;  
  const TVector3 priVtx = 0.;
  //  TVector3 L3d = operator - (secVtx,priVtx);
  TVector3 L3d = secVtx - priVtx;

  TVector3 P3sumTagTrk = P4secVtx.Vect();
  double sinThetaBHadron = TMath::Sin(L3d.Theta());  // angle between B hadron flight direction and detector z axis
  double cosTheta = 0.;                               // angle between B hadron flight direction and tag jet axis
  if ( P3sumTagTrk.Mag() != 0 && L3d.Mag() != 0 )
    cosTheta = P3sumTagTrk.Dot(L3d)/(P3sumTagTrk.Mag()*L3d.Mag()); 

  double PchL = P3sumTagTrk.Mag()*fabs(cosTheta);
  double PT   = P3sumTagTrk.Mag()*sqrt(1-pow(cosTheta,2));  // HP: B hadron rest frame/ P0T==PchT
  double Mch  = P4secVtx.M();
  double M0sq = MB2_ - 2*MB_*sqrt(Mch*Mch+PT*PT) + Mch*Mch;
  if ( M0sq < 0 ) M0sq = 0;
  if (Mch*Mch+PT*PT != 0) P0L = (MB2_ - (Mch*Mch + PT*PT) - (M0sq + PT*PT))*PchL/(2*(Mch*Mch+PT*PT));

  double PBrec = PchL + P0L;
  ptBhad   = PBrec*sinThetaBHadron;
  
// 	cout << "----------------------------------------------------------------------------------------------------" << endl;
// 	cout << "   SECONDARY VERTEX :: X: " << tmp_secondary_vertex.X() << "  Y: " << tmp_secondary_vertex.Y() << "  Z: " << tmp_secondary_vertex.Z() << endl;
// 	cout << "    PRIMARY VERTEX  :: X: " << _primary_vertex.X() << "  Y: " << _primary_vertex.Y() << "  Z: " << _primary_vertex.Z() << endl;
// 	cout << "----------------------------------------------------------------------------------------------------" << endl;
// 	if (P0L<0)
// 	  cout << "WARNING" << endl;
//      else
//        cout << "OK" << endl;
// 	cout << "----------------------------------------------------------------------------------------------------" << endl;
// 	cout << "ijet: " << ijet << endl;
// 	cout << "_TagJet[ije]: " << _TagJet[ijet] << endl;
// 	cout << "Pmod: " << P3sumTagTrk.Mag() << endl;
// 	cout << "Mch: " << Mch << endl;
// 	cout << "----------------------------------------------------------------------------------------------------" << endl;

//      ptBhad[ijet] = ptBhad;
//      ptBneu[ijet] = P0L;
//      L3dTag[ijet] = L3d;
}



// ------------ method called once each job just before starting event loop  ------------
void 
VBFHZZllbbBhadronReconstruction::beginJob(const edm::EventSetup&)
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
VBFHZZllbbBhadronReconstruction::endJob() {
}

//DEFINE_FWK_MODULE(VBFHZZllbbBhadronReconstruction);

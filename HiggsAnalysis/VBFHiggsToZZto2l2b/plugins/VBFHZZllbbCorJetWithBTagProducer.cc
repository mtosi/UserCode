// -*- C++ -*-
//
// Package:    VBFHZZllbbCorJetWithBTagProducer
// Class:      VBFHZZllbbCorJetWithBTagProducer
// 
/**\class VBFHZZllbbCorJetWithBTagProducer VBFHZZllbbCorJetWithBTagProducer.cc HiggsAnalysis/VBFHZZllbbCorJetWithBTagProducer/plugins/VBFHZZllbbCorJetWithBTagProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Mia TOSI
//         Created:  Wed Mar 18 16:33:51 CET 2009
// $Id: VBFHZZllbbCorJetWithBTagProducer.cc,v 1.3 2009/07/06 13:15:48 tosi Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbCorJetWithBTagProducer.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

#include "DataFormats/JetReco/interface/JetFloatAssociation.h"

#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/CorJetWithBTagDiscr.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/VBFHZZllbbUtils.h"

using namespace std;
using namespace edm;
using namespace reco;
// using namespace vbfhzz2l2b;

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//

VBFHZZllbbCorJetWithBTagProducer::VBFHZZllbbCorJetWithBTagProducer(const edm::ParameterSet& iConfig) :
  bTagConfigLabel_ ( iConfig.getParameter< std::vector<edm::ParameterSet> >( "bTagConfig" ) ) ,
  jetCorrectionService_ ( iConfig.getParameter<std::string>( "jetCorrectionService" ) )
{
  //if do put with a label
  produces<vbfhzz2l2b::CorJetWithBTagDiscrCollection>( "corJetWithBTagDiscr" );

  //now do what ever other initialization is needed
  std::cout << "[VBFHZZllbbCorJetWithBTagProducer::VBFHZZllbbCorJetWithBTagProducer] DONE" << std::endl;
}


VBFHZZllbbCorJetWithBTagProducer::~VBFHZZllbbCorJetWithBTagProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  std::cout << "VBFHZZllbbCorJetWithBTagProducer::~VBFHZZllbbCorJetWithBTagProducer] DONE" << std::endl;
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void VBFHZZllbbCorJetWithBTagProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //  std::cout << "[VBFHZZllbbCorJetWithBTagProducer::produce]" << std::endl;
  //  edm::Handle<reco::SecondaryVertexTagInfoCollection> secondaryVtxTagInfosHandle;
  //  iEvent.getByLabel( "secondaryVertexTagInfos", secondaryVtxTagInfosHandle );
  //  std::cout << "secondaryVtxTagInfosHandle" << std::endl;
  //  std::cout << "secondaryVtxTagInfosHandle->size(): " << secondaryVtxTagInfosHandle->size() << std::endl;

  std::auto_ptr<vbfhzz2l2b::CorJetWithBTagDiscrCollection> jetCollection;
  
  const JetCorrector* corrector = JetCorrector::getJetCorrector (jetCorrectionService_,iSetup);

  std::vector<edm::Handle<reco::JetTagCollection> > bTagHandle_vec;
  std::vector<reco::JetTagCollection>  bTagCollection_vec;
  for ( unsigned int bTagLabelIndex = 0; bTagLabelIndex != bJetTagInputTags_.size(); bTagLabelIndex++ ) {
    
    edm::Handle<reco::JetTagCollection> bTagHandle;
    iEvent.getByLabel(bJetTagInputTags_[bTagLabelIndex], bTagHandle);
    bTagHandle_vec.push_back(bTagHandle);
    
    const reco::JetTagCollection & bTagCollection = *(bTagHandle_vec[bTagLabelIndex].product());
    bTagCollection_vec.push_back(bTagCollection);
    //    std::cout << "Found " << bTagCollection_vec[bTagLabelIndex].size() 
    //	      << " B candidates in collection " << bJetTagInputTags_[bTagLabelIndex]
    //	      << std::endl;
  }
  

  std::vector<reco::JetBaseRef> jets = reco::JetFloatAssociation::allJets(*bTagHandle_vec[0]);
  //  std::cout << "[VBFHZZllbbCorJetWithBTagProducer::produce] found " << jets.size()
  //	    << " B candidates" << std::endl;
  if (jets.size() != 0) {
    edm::RefToBase<Jet> jj = jets[0];
    jetCollection.reset( new vbfhzz2l2b::CorJetWithBTagDiscrCollection(edm::RefToBaseProd<Jet>(jj)));
  }
  else jetCollection.reset( new vbfhzz2l2b::CorJetWithBTagDiscrCollection(edm::RefToBaseProd<Jet>()));

//  for ( unsigned int taggerIndex = 0; taggerIndex != bJetTagInputTags_.size(); taggerIndex++ ) {
//    for ( reco::SecondaryVertexTagInfoCollection::const_iterator secondaryVtxTagInfos_itr = secondaryVtxTagInfosHandle->begin();
//	  secondaryVtxTagInfos_itr != secondaryVtxTagInfosHandle->end(); ++secondaryVtxTagInfos_itr ) {
//      
//      RefToBase< Jet >        jetRef          = secondaryVtxTagInfos_itr -> jet();
//      std::cout << "jetRef et: " << jetRef->et() <<std::endl;
//      double discr = reco::JetFloatAssociation::getValue(*bTagHandle_vec[taggerIndex], jetRef);
//      std::cout << "[VBFHZZllbbCorJetWithBTagProducer::produce] taggerIndex: " << taggerIndex 
//		<< " --> discriminator: " << discr << std::endl;
//    }
//  }
//  std::cout << "[VBFHZZllbbCorJetWithBTagProducer::produce] DONE w/ SecondaryVertexTagInfoCollection" << std::endl;
//  std::cout << "******************" << std::endl;

  int jetIndex = 0;    
  for ( std::vector<reco::JetBaseRef>::const_iterator jet = jets.begin(); 
	jet != jets.end();
	++jet, jetIndex++ ) {
    //    std::cout << "[VBFHZZllbbCorJetWithBTagProducer::produce] jetIndex: " << jetIndex << std::endl;

    vbfhzz2l2b::DiscriminatorVector bTagDiscrs;
    for ( unsigned int taggerIndex = 0; taggerIndex != bJetTagInputTags_.size(); taggerIndex++ ) {
      double discr = reco::JetFloatAssociation::getValue(*bTagHandle_vec[taggerIndex], **jet);
      //      std::cout << "[VBFHZZllbbCorJetWithBTagProducer::produce] taggerIndex: " << taggerIndex 
      //      		<< " --> discriminator: " << discr << std::endl;
      bTagDiscrs.push_back(discr);
    }

    double uncorEt = (*jet)->et();
    double uncorPt = (*jet)->pt();
    double corScale = corrector->correction( **jet );
    double corEt = uncorEt*corScale;
    double corPt = uncorPt*corScale;
    //    std::cout << "[VBFHZZllbbCorJetWithBTagProducer::produce] corScale: " << corScale 
    //	      << " => uncorrected et: " << uncorEt
    //	      << " -----> corrected et: " << corEt 
    //	      << " => uncorrected pt: " << uncorPt
    //	      << " -----> corrected pt: " << corPt 
    //	      << std::endl;

    vbfhzz2l2b::CorJetBTagDiscrAssociation::CorBTagDiscrData corBTagDiscrData;
    corBTagDiscrData.corEt_  = corEt;
    corBTagDiscrData.corPt_  = corPt;
    corBTagDiscrData.discrVec_ = bTagDiscrs;
    corBTagDiscrData.highEffDiscr_ = bTagDiscrs[vbfhzz2l2b::bTaggerCode("HIGHEFF")];
    corBTagDiscrData.highPurDiscr_ = bTagDiscrs[vbfhzz2l2b::bTaggerCode("HIGHPUR")];
    corBTagDiscrData.compoSVDiscr_ = bTagDiscrs[vbfhzz2l2b::bTaggerCode("COMBSECVTX")];
    corBTagDiscrData.jetProbDiscr_ = bTagDiscrs[vbfhzz2l2b::bTaggerCode("JETPROB")];

    (*jetCollection)[*jet] = corBTagDiscrData;

    double emFrac = (dynamic_cast<const reco::CaloJet*>(&**jet))->emEnergyFraction();
    //    std::cout << "[VBFHZZllbbCorJetWithBTagProducer::produce] emFrac: " << emFrac << std::endl;
  }

  //  std::cout << "[VBFHZZllbbCorJetWithBTagProducer::produce] jetCollection->size(): " << jetCollection->size() << std::endl;

  // puts Corrected Jet Collection with BTags discriminator into event
  iEvent.put(jetCollection, "corJetWithBTagDiscr");

  std::cout << "[VBFHZZllbbCorJetWithBTagProducer::produce] DONE" << std::endl;

}

// ------------ method called once each job just before starting event loop  ------------
void VBFHZZllbbCorJetWithBTagProducer::beginJob(const edm::EventSetup&)
{
  for (unsigned int iModule = 0; iModule != bTagConfigLabel_.size(); ++iModule) {
    std::string dataFormatType = "JetTag";
    edm::InputTag bTagLabel = bTagConfigLabel_[iModule].getParameter<edm::InputTag>("label");
    bJetTagInputTags_.push_back( bTagLabel );
  }
  std::cout << "[VBFHZZllbbCorJetWithBTagProducer::beginJob] DONE" << std::endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void VBFHZZllbbCorJetWithBTagProducer::endJob() {
  std::cout << "[VBFHZZllbbCorJetWithBTagProducer::endJob] DONE" << std::endl;
}

DEFINE_FWK_MODULE(VBFHZZllbbCorJetWithBTagProducer);

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/UtilAlgos/interface/ObjectSelector.h"
#include "DataFormats/Common/interface/RefVector.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"


//define this as a plug-in
DEFINE_SEAL_MODULE () ;

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbBTagInfoAnalyzer.h"
DEFINE_ANOTHER_FWK_MODULE(VBFHZZllbbBTagInfoAnalyzer);

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbBhadronReconstruction.h"
DEFINE_ANOTHER_FWK_MODULE(VBFHZZllbbBhadronReconstruction);

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/newSimpleNtple.h"
DEFINE_ANOTHER_FWK_MODULE(newSimpleNtple); 

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/SimpleNtple.h"
DEFINE_ANOTHER_FWK_MODULE(SimpleNtple); 

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbDeltaRAnalyzer.h"
DEFINE_ANOTHER_FWK_MODULE(VBFHZZllbbDeltaRAnalyzer);

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbAnalyzer.h"
DEFINE_ANOTHER_FWK_MODULE(VBFHZZllbbAnalyzer);

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbJetMatching.h"
DEFINE_ANOTHER_FWK_MODULE(VBFHZZllbbJetMatching);

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbMCfilterValidation.h"
DEFINE_FWK_MODULE(VBFHZZllbbMCfilterValidation);

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbMCvalidation.h"
DEFINE_FWK_MODULE(VBFHZZllbbMCvalidation);

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbRAWCORBTAGtesting.h"
DEFINE_FWK_MODULE(VBFHZZllbbRAWCORBTAGtesting);


#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbMCzTObbFilter.h"
DEFINE_ANOTHER_FWK_MODULE (VBFHZZllbbMCzTObbFilter); 

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbMCbQuarkFilter.h"
DEFINE_ANOTHER_FWK_MODULE (VBFHZZllbbMCbQuarkFilter); 

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbMCprocessFilter.h"
DEFINE_ANOTHER_FWK_MODULE (VBFHZZllbbMCprocessFilter); 

//#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbValidationPlots.h"
//DEFINE_ANOTHER_FWK_MODULE (VBFHZZllbbValidationPlots); 

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbDisplay.h"
DEFINE_ANOTHER_FWK_MODULE (VBFHZZllbbDisplay) ; 

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbPreSelection.h"
DEFINE_ANOTHER_FWK_MODULE (VBFHZZllbbPreSelection) ;

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbCommonPreselection.h"
DEFINE_ANOTHER_FWK_MODULE (VBFHZZllbbCommonPreselection);

//#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbMuonSelector.h"
//typedef ObjectSelector<VBFHZZllbbMuonSelector> VBFHZZllbbMuonSelection;
//typedef ObjectSelector<VBFHZZllbbMuonSelector, edm::RefVector<reco::MuonCollection> > VBFHZZllbbMuonSelectionRef;
//DEFINE_ANOTHER_FWK_MODULE (VBFHZZllbbMuonSelection);
//DEFINE_ANOTHER_FWK_MODULE (VBFHZZllbbMuonSelectionRef);

//#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbElectronSelector.h"
//typedef ObjectSelector<VBFHZZllbbElectronSelector> VBFHZZllbbElectronSelection;
//typedef ObjectSelector<VBFHZZllbbElectronSelector, edm::RefVector<reco::PixelMatchGsfElectronCollection> > VBFHZZllbbElectronSelectionRef;
//DEFINE_ANOTHER_FWK_MODULE (VBFHZZllbbElectronSelection);
//DEFINE_ANOTHER_FWK_MODULE (VBFHZZllbbElectronSelectionRef);

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbCommonOfflineSelection.h"
DEFINE_ANOTHER_FWK_MODULE (VBFHZZllbbCommonOfflineSelection);



#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbMuonIsolationProducer.h"
DEFINE_ANOTHER_FWK_MODULE (VBFHZZllbbMuonIsolationProducer);

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbElectronIsolationProducer.h"
DEFINE_ANOTHER_FWK_MODULE (VBFHZZllbbElectronIsolationProducer);

//#include "HiggsAnalysis/VBFHiggsToZZto2l2b/plugins/VBFHZZllbbCorJetWithBTagProducer.cc"
//DEFINE_ANOTHER_FWK_MODULE(VBFHZZllbbCorJetWithBTagProducer);


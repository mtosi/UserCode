#ifndef classes_h
#define classes_h

#include "Rtypes.h" 

#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Common/interface/RefProd.h" 
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/RefHolder.h"
#include "DataFormats/Common/interface/Holder.h"

#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "DataFormats/JetReco/interface/BasicJet.h" 
#include "DataFormats/JetReco/interface/CaloJetCollection.h" 
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/GenericJetCollection.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/CorJetBTagDiscrAssociation.h"
#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/CorJetWithBTagDiscr.h"

//#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/SimpleNtpleObj.h"

using namespace vbfhzz2l2b;
using namespace std;
using namespace reco;
 
namespace {
  struct dictionary {

    edm::Wrapper<CorJetBTagDiscrAssociation::Container>  jea_c_w;
    CorJetBTagDiscrAssociation::Container       jea_c;
    CorJetBTagDiscrAssociation::Ref             jea_r;
    CorJetBTagDiscrAssociation::RefProd         jea_rp;
    CorJetBTagDiscrAssociation::RefVector       jea_rv;


//    SimpleNtpleObj::EVT evt1;
//    edm::Wrapper<SimpleNtpleObj::EVT> wevt1;
//    SimpleNtpleObj::JET jet1;
//    edm::Wrapper<SimpleNtpleObj::JET> wjet1;
//    SimpleNtpleObj::MUON muon1;
//    edm::Wrapper<SimpleNtpleObj::MUON> wmuon1;
//    SimpleNtpleObj::ELECTRON electron1;
//    edm::Wrapper<SimpleNtpleObj::ELECTRON> welectron1;
//    SimpleNtpleObj::ZHAD zhad1;
//    edm::Wrapper<SimpleNtpleObj::ZHAD> wzhad1;


  };
}

#endif // classes_h

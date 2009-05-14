#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "HiggsAnalysis/Skimming/interface/VBFHZZllbbSkim.h"
#include "HiggsAnalysis/Skimming/interface/VBFHZZllbbSkimEff.h"
#include "HiggsAnalysis/Skimming/interface/VBFHZZllbbPreFilter.h"

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(VBFHZZllbbSkim);
DEFINE_ANOTHER_FWK_MODULE(VBFHZZllbbSkimEff);
DEFINE_ANOTHER_FWK_MODULE(VBFHZZllbbPreFilter);



import FWCore.ParameterSet.Config as cms

from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbMCprocessFilter_cfi import *
vbfhzzllbbMCprocessFilter.vbfFlag = cms.bool(False)

#from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbZtoMuMuSelectorSequence_cff import *

from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbcorjetwithbtagproducer_cff import *
from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbanalyzer_cfi import *
#vbfHZZllbbAnalyzer.muonLabel = cms.InputTag('vbfhzzllbbMuonSelector')

vbfhzzllbbAnalysisSequence = cms.Sequence(
    vbfhzzllbbMCprocessFilter *
    vbfhzzllbbCorJetWithBTagSequence  *
#    vbfhzzllbbZtoMuMuSequence *
    vbfHZZllbbAnalyzer
)

########
# Path #
########
vbfhzzllbbAnalysisPath = cms.Path(vbfhzzllbbAnalysisSequence)

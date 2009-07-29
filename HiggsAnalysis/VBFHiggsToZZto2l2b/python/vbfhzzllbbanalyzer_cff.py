import FWCore.ParameterSet.Config as cms

from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbcorjetwithbtagproducer_cff import *
from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbanalyzer_cfi import *

vbfhzzllbbAnalysisSequence = cms.Sequence(
    vbfhzzllbbCorJetWithBTagSequence  *
    vbfHZZllbbAnalyzer
)

########
# Path #
########
vbfhzzllbbAnalysisPath = cms.Path(vbfhzzllbbAnalysisSequence)

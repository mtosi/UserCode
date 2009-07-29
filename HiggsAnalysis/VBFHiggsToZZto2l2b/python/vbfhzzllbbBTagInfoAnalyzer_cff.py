import FWCore.ParameterSet.Config as cms

from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbBTagInfoAnalyzer_cfi import *

from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbcorjetwithbtagproducer_cff import *

vbfhzzllbbBTagInfoAnalyzerSequence = cms.Sequence(
    vbfhzzllbbCorJetWithBTagSequence *
    vbfhzzllbbBTagInfoAnalyzer
)

#vbfhzzllbbBTagInfo = cms.Path(
#    vbfhzzllbbBTagInfoAnalyzerSequence
#)

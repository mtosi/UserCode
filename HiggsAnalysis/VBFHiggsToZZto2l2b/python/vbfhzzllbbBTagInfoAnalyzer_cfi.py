import FWCore.ParameterSet.Config as cms


vbfhzzllbbBTagInfoAnalyzer = cms.EDAnalyzer("VBFHZZllbbBTagInfoAnalyzer",
  impactParameterTagInfosLabel = cms.InputTag('impactParameterTagInfos'),
  secondaryVertexTagInfosLabel = cms.InputTag('secondaryVertexTagInfos')
)

from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbcorjetwithbtagproducer_cfi import *
vbfhzzllbbBTagInfoAnalyzerSequence = cms.Sequence(
    vbfhzzllbbCorJetWithBTagSequence *
    vbfhzzllbbBTagInfoAnalyzer
)

#vbfhzzllbbBTagInfo = cms.Path(
#    vbfhzzllbbBTagInfoAnalyzerSequence
#)

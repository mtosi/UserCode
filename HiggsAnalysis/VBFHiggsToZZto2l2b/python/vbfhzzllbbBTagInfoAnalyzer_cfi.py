import FWCore.ParameterSet.Config as cms

vbfhzzllbbBTagInfoAnalyzer = cms.EDAnalyzer("VBFHZZllbbBTagInfoAnalyzer",
  impactParameterTagInfosLabel = cms.InputTag('impactParameterTagInfos'),
  secondaryVertexTagInfosLabel = cms.InputTag('secondaryVertexTagInfos')
)

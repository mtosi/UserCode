import FWCore.ParameterSet.Config as cms

vbfhzzllbbCorJetWithBTagAnalyzer = cms.EDAnalyzer('VBFHZZllbbCorJetWithBTagAnalyzer',
                                             corJetWithBTagLabel = cms.string('vbfhzzllbbCorJetWithBTagProd')
)

vbfhzzllbbCorJetWithBTagAna = cms.Sequence(
    vbfhzzllbbCorJetWithBTagAnalyzer
)

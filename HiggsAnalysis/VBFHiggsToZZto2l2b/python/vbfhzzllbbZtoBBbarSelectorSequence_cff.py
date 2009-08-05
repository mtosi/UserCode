import FWCore.ParameterSet.Config as cms

goodBjetMCMatch = cms.EDFilter("CaloJetMatcher",
#    src = cms.InputTag('vbfhzzllbbCorJetWithBTagProd'),
    src = cms.InputTag('iterativeCone5CaloJets'),
    maxDPtRel = cms.double(1.0),
    mcPdgId = cms.vint32(5), ## b-quarks
    mcStatus = cms.vint32(1),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = cms.double(0.5),
    checkCharge = cms.bool(False),
    resolveAmbiguities = cms.bool(True),
    matched = cms.InputTag("genParticles")
)

# Z to bbbar
from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbZtoBBbar_cfi import *
vbfhzzllbbZtoBBbar.decay = cms.string('goodBjetMCMatch goodBjetMCMatch')

vbfHZZllbbZtoBBbarPath = cms.Path(
    goodBjetMCMatch    *
    vbfhzzllbbZtoBBbar
)

vbfhzzllbbZtoBBbarSequence = cms.Sequence(
    goodBjetMCMatch    *
    vbfhzzllbbZtoBBbar
)

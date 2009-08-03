import FWCore.ParameterSet.Config as cms

# muon selection
from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbMuonSequence_cff import *

goodMuonMCMatch = cms.EDFilter("MCMatcher",
    src = cms.InputTag("vbfhzzllbbMuonSelector"),
    maxDPtRel = cms.double(1.0),
    mcPdgId = cms.vint32(13), ## muons
    mcStatus = cms.vint32(1),
    resolveByMatchQuality = cms.bool(False),
    maxDeltaR = cms.double(0.15),
    checkCharge = cms.bool(True),
    resolveAmbiguities = cms.bool(False),
    matched = cms.InputTag("genParticles")
)

# 2 muon filter
from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbTwoMuonFilter_cfi import *
vbfhzzllbbTwoMuonFilter.src = cms.InputTag("vbfhzzllbbMuonSelector")

# Z to Muon@+ Muon@-
from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbZtoMuMu_cfi import *
vbfhzzllbbZtoMuMu.decay = cms.string('vbfhzzllbbMuonSelector@+ vbfhzzllbbMuonSelector@-')

vbfHZZllbbZtoMuMuPath = cms.Path(
    vbfhzzllbbMuonSequence  *
    vbfhzzllbbTwoMuonFilter *
    vbfhzzllbbZtoMuMu
)

vbfhzzllbbZtoMuMuSequence = cms.Sequence(
    vbfhzzllbbMuonSequence  *
    goodMuonMCMatch         *
    vbfhzzllbbTwoMuonFilter *
    vbfhzzllbbZtoMuMu
)

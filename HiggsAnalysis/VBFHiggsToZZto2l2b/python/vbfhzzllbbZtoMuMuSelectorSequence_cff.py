import FWCore.ParameterSet.Config as cms

# muon selection
from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfHZZllbbMuonSelector_cfi import *

MuonPtIsolation = cms.EDProducer("MuonPtIsolationProducer",
    src = cms.InputTag("globalMuons"),
    d0Max = cms.double(20.0),
    dRMin = cms.double(0.0),
    dRMax = cms.double(0.3),
    elements = cms.InputTag("generalTracks"),
    ptMin = cms.double(0.0),
    dzMax = cms.double(20.0)

)
MuonPtIsolation05 = cms.EDProducer("MuonPtIsolationProducer",
    dzMax = cms.double(20.0),
    src = cms.InputTag("globalMuons"),
    dRMin = cms.double(0.0),
    elements = cms.InputTag("generalTracks"),
    dRMax = cms.double(0.5)
)

goodMuonMCMatch = cms.EDFilter("MCMatcher",
    src = cms.InputTag("vbfHZZllbbMuonSelector"),
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
vbfhzzllbbTwoMuonFilter.src = cms.InputTag("vbfHZZllbbMuonSelector")

# Z to Muon@+ Muon@-
from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbZtoMuMu_cfi import *
vbfhzzllbbZtoMuMu.decay = cms.string('vbfHZZllbbMuonSelector@+ vbfHZZllbbMuonSelector@-')

vbfHZZllbbZtoMuMuPath = cms.Path(
    vbfHZZllbbMuonSelector  *
    vbfhzzllbbTwoMuonFilter *
    vbfhzzllbbZtoMuMu
)

vbfhzzllbbZtoMuMuSequence = cms.Sequence(
    vbfHZZllbbMuonSelector  *
    goodMuonMCMatch         *
    vbfhzzllbbTwoMuonFilter *
    vbfhzzllbbZtoMuMu
)

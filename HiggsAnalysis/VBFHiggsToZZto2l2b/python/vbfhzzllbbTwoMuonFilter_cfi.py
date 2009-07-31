import FWCore.ParameterSet.Config as cms

vbfhzzllbbTwoMuonFilter = cms.EDFilter("MuonCountFilter",
    src = cms.InputTag("muons"),
    minNumber = cms.uint32(2)
)

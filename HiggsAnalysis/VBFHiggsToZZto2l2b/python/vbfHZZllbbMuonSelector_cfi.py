import FWCore.ParameterSet.Config as cms

vbfHZZllbbMuonSelector = cms.EDProducer("VBFHZZllbbMuonSelector",
    sourceLabel = cms.InputTag("muons"),
    sourceMinPtBarrelCut = cms.double(5.0),
    sourceMinPtEndcapCut = cms.double(3.0),
    sourceMinPEndcapCut  = cms.double(9.0),
    sourceMaxEtaCut      = cms.double(2.4)
)

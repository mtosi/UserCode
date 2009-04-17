import FWCore.ParameterSet.Config as cms

overlapElectronResolver = cms.EDFilter("AmbResolver",
    filter = cms.bool(False),
    src = cms.InputTag("pixelMatchGsfElectrons")
)



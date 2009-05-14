import FWCore.ParameterSet.Config as cms

vbfhzzllbbMCbQuarkFilter = cms.EDFilter("VBFHZZllbbMCbQuarkFilter",
        genParticleLabel = cms.InputTag('genParticles'),
        signal = cms.int32(0)
)



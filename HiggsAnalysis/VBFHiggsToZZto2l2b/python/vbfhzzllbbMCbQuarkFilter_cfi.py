import FWCore.ParameterSet.Config as cms

vbfhzzllbbMCbQuarkFilter = cms.EDFilter("VBFHZZllbbMCbQuarkFilter",
        genParticleLabel = cms.InputTag('genParticles'),
        signal = cms.int32(1)  # 1:Signal,  0:Background
)



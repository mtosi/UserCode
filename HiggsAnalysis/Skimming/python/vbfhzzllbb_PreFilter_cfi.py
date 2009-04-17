import FWCore.ParameterSet.Config as cms

vbfhzzllbbPreFilter = cms.EDFilter('VBFHZZllbbPreFilter',
                                   genParticleLabel = cms.untracked.InputTag('genParticles'),
                                   VBFHZZllbbPreFilterDebug = cms.bool(False),
                                   VBFHZZllbbPreFilterLeptonFlavour = cms.int32(0)
                                   )

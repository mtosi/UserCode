import FWCore.ParameterSet.Config as cms

hTozzTo4leptons = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('zToMuMu zToBjetBjet'),
    cut = cms.string('0.0 < mass < 20000.0')
)


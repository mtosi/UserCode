import FWCore.ParameterSet.Config as cms

zToEE = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('VBFHZZllbbElectronSelector@+ VBFHZZllbbElectronSelector@-'),
    cut = cms.string('0.0 < mass < 20000.0')
)


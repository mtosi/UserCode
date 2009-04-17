import FWCore.ParameterSet.Config as cms

zToMuMu = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('vbfHZZllbbMuonSelector@+ vbfHZZllbbMuonSelector@-'),
    cut = cms.string('0.0 < mass < 20000.0')
)


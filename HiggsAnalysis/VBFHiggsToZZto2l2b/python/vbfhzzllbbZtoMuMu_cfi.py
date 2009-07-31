import FWCore.ParameterSet.Config as cms

vbfhzzllbbZtoMuMu = cms.EDProducer("CandViewShallowCloneCombiner",
    checkCharge = cms.bool(True),
    decay = cms.string('muons@+ muons@-'),
    cut = cms.string('20.0 < mass < 200.0')
)


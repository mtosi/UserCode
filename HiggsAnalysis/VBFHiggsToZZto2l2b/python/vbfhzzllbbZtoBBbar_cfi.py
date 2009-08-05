import FWCore.ParameterSet.Config as cms

vbfhzzllbbZtoBBbar = cms.EDProducer("CandViewShallowCloneCombiner",
    checkCharge = cms.bool(False),
    decay = cms.string('iterativeCone5CaloJets iterativeCone5CaloJets'),
    cut = cms.string('20.0 < mass < 200.0')
)


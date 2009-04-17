import FWCore.ParameterSet.Config as cms

vbfHZZllbbElectronSelector = cms.EDProducer("VBFHZZllbbElectronSelector",
    sourceLabel   = cms.InputTag("overlapElectronResolver"),
    sourceIDLabel = cms.InputTag("eidClassLoose"),
    sourceEtaMaxCut = cms.double(2.5),
    sourcePtMinCut  = cms.double(5.0)
)



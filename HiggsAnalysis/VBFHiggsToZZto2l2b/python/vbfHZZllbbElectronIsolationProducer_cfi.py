import FWCore.ParameterSet.Config as cms

VBFHZZllbbElectronIsolationProducer = cms.EDProducer("VBFHZZllbbElectronIsolationProducer",
    ElectronsLabel     = cms.InputTag("VBFHZZllbbElectronSelector"),
    ElectronsVetoLabel = cms.InputTag("overlapElectronResolver"),
    TracksLabel        = cms.InputTag("generalTracks"),
    MuonsLabel         = cms.InputTag("muons"),
    isolationCut         = cms.double(0.7),
    isolationConeCut     = cms.double(0.25),
    isolationConeVetoCut = cms.double(0.015)
)



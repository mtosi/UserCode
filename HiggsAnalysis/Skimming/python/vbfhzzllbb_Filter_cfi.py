import FWCore.ParameterSet.Config as cms

vbfhzzllbbFilter = cms.EDFilter("VBFHZZllbbSkim",
    VBFHZZllbbDebug = cms.bool(False),
    # Collection to be accessed
    muonLabel     = cms.InputTag('muons'),
    electronLabel = cms.InputTag('pixelMatchGsfElectrons'),
    # Minimum number of identified leptons above pt threshold
    tightLeptonMinimumNumber = cms.int32(1),
    softLeptonMinimumNumber  = cms.int32(2),
    # Pt threshold for leptons
    tightMinimumPt = cms.double(10.0),
    softMinimumPt  = cms.double(5.0)
)



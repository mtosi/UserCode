import FWCore.ParameterSet.Config as cms

vbfhzzllbbMCprocessFilter = cms.EDFilter("VBFHZZllbbMCprocessFilter",
  whichSim = cms.int32(0),
  signal   = cms.int32(0)
)



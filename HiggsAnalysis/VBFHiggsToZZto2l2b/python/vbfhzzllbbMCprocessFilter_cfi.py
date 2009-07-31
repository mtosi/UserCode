import FWCore.ParameterSet.Config as cms

vbfhzzllbbMCprocessFilter = cms.EDFilter("VBFHZZllbbMCprocessFilter",
  whichSim = cms.int32(1), # 0:FastSim, 1:FullSim
  signal   = cms.int32(1),
  vbfFlag  = cms.bool(True)                                       
)



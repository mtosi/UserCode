import FWCore.ParameterSet.Config as cms

from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfHZZllbbMuonSelector_cfi import *
vbfHZZllbbMuonSequence = cms.Sequence(vbfHZZllbbMuonSelector)


import FWCore.ParameterSet.Config as cms

from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbMuonSelector_cfi import *

vbfhzzllbbMuonSequence = cms.Sequence(vbfhzzllbbMuonSelector)

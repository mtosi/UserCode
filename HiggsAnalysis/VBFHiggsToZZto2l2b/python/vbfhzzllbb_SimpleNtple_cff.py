import FWCore.ParameterSet.Config as cms

from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbcorjetwithbtagproducer_cff import *
from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbb_SimpleNtple_cfi import *
from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbb_newSimpleNtple_cfi import *

vbfhzzllbbSimpleNtupleSequence = cms.Sequence(
    vbfhzzllbbCorJetWithBTagSequence  *
    vbfhzzllbbSimpleNtple
)

vbfhzzllbbNewSimpleNtupleSequence = cms.Sequence(
    vbfhzzllbbCorJetWithBTagSequence  *
    vbfhzzllbbNewSimpleNtple
)

########
# Path #
########
vbfhzzllbbSimpleNtuplePath = cms.Path(vbfhzzllbbSimpleNtupleSequence)
vbfhzzllbbNewSimpleNtuplePath = cms.Path(vbfhzzllbbNewSimpleNtupleSequence)

import FWCore.ParameterSet.Config as cms

from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbcorjetwithbtagproducer_cff import *
from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbb_SimpleNtple_cfi import *

vbfhzzllbbSimpleNtupleSequence = cms.Sequence(
    vbfhzzllbbCorJetWithBTagSequence  *
    vbfhzzllbbSimpleNtple
)

########
# Path #
########
vbfhzzllbbSimpleNtuplePath = cms.Path(vbfhzzllbbSimpleNtupleSequence)

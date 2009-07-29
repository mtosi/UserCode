import FWCore.ParameterSet.Config as cms

from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbcorjetwithbtagproducer_cff import *
from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbb_newSimpleNtple_cfi import *
#vbfhzzllbbNewSimpleNtple.corIC5CaloJetsWithBTagLabel = cms.string('iterativeCone5CaloJets')

vbfNewSimpleNtupleSequence = cms.Sequence(
    vbfhzzllbbCorJetWithBTagSequence  *
    vbfhzzllbbNewSimpleNtple
)

########
# Path #
########
pathTree = cms.Path(vbfNewSimpleNtupleSequence)

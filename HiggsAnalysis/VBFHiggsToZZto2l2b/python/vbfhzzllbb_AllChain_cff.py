#!/usr/bin/env cmsRun
import FWCore.ParameterSet.Config as cms

from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbb_Sequences_cff import *
from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbcorjetwithbtagproducer_cff import *
from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbb_SimpleNtple_cfi import *

VBFAllChainSequence = cms.Sequence(
    vbfhzzllbbSequence *
    vbfhzzllbbCorJetWithBTagSequence  *
    vbfhzzllbbSimpleNtple
)

########
# Path #
########
pathTree = cms.Path(VBFAllChainSequence)

#!/usr/bin/env cmsRun
import FWCore.ParameterSet.Config as cms

from HiggsAnalysis.Skimming.vbfhzzllbb_Sequences_cff import *
from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbcorjetwithbtagproducer_cfi import *
from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbb_SimpleNtple_cfi import *
#from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbb_newSimpleNtple_cfi import *

VBFAllChainSequence = cms.Sequence(
    vbfhzzllbbSequence *
    vbfhzzllbbCorJetWithBTagSequence  *
    vbfhzzllbbSimpleNtple 
#    vbfhzzllbbNewSimpleNtple
)

########
# Path #
########
pathTree = cms.Path(VBFAllChainSequence)

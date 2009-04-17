#!/usr/bin/env cmsRun
import FWCore.ParameterSet.Config as cms

from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbPreSelection_cff import *
from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbJetMatching_cfi import *

vbfhzzllbbJetMatchingSequence = cms.Sequence(
    vbfhzzllbbSequence *
    vbfhzzllbbJetMatching
)

########
# Path #
########
jetMatchinPath = cms.Path(vbfhzzllbbJetMatchingSequence)

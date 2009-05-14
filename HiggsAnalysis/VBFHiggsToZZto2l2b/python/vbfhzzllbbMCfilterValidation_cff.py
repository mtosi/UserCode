#!/usr/bin/env cmsRun
import FWCore.ParameterSet.Config as cms

#from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbPreSelection_cff import *
from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbMCfilterValidation_cfi import *

vbfhzzllbbMCfilterValidationSequence = cms.Sequence(
#    vbfhzzllbbPreSelectionSequence *
    vbfhzzllbbMCfilterValidation
)

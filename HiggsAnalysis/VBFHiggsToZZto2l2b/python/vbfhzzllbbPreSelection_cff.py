#!/usr/bin/env cmsRun
import FWCore.ParameterSet.Config as cms

from HiggsAnalysis.Skimming.vbfhzzllbb_HLTPaths_cfi import *
from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbcorjetwithbtagproducer_cfi import *
from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbMCprocessFilter_cfi import *
from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbPreSelection_cfi import *
#from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbb_SimpleNtple_cfi import *

vbfhzzllbbTrigReport = cms.EDFilter("HLTrigReport",
                                    HLTriggerResults = cms.InputTag('TriggerResults',
                                                                    '',
                                                                    'HLT'
                                                                    )
                                    )

vbfhzzllbbSequence = cms.Sequence(
    vbfhzzllbbTrigReport *
    vbfhzzllbbHLTFilter *
    vbfhzzllbbMCprocessFilter *
    vbfhzzllbbCorJetWithBTagSequence *
    vbfhzzllbbPreSelection 
#    vbfhzzllbbSimpleNtple
)

########
# Path #
########
path = cms.Path(vbfhzzllbbSequence)

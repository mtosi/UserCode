import FWCore.ParameterSet.Config as cms

from HiggsAnalysis.Skimming.vbfhzzllbb_HLTPaths_cfi import *
from HiggsAnalysis.Skimming.vbfhzzllbb_Filter_cfi import *
vbfhzzllbbTrigReport = cms.EDFilter("HLTrigReport",
                                    HLTriggerResults = cms.InputTag('TriggerResults',
                                                                    '',
                                                                    'HLT'
                                                                    )
                                    )
vbfhzzllbbSequence = cms.Sequence( vbfhzzllbbTrigReport
                                   * vbfhzzllbbHLTFilter
                                   * vbfhzzllbbFilter
                                   )


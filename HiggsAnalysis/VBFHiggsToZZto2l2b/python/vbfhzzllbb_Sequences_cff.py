import FWCore.ParameterSet.Config as cms

from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbb_HLTPaths_cfi import *
from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbb_Filter_cfi import *
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


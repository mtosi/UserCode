#!/usr/bin/env cmsRun
import FWCore.ParameterSet.Config as cms

process = cms.Process("VBFHZZllbbFILTEREFF")

process.load('Configuration/StandardSequences/Services_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
## to write every 10th event
process.MessageLogger.cerr.FwkReport.reportEvery = 10
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.cerr.threshold = 'ERROR'
process.MessageLogger.categories.append('VBFHZZllbbFilterEffSummary')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    default                   = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
    VBFHZZllbbFilterEffSummary = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True) # default values if False
)

## this defines the input files
from HiggsAnalysis.VBFHiggsToZZto2l2b.Data.H150_ZZ_qqllSummer08_IDEALV9v2_GENSIMRECO_FULL_Input_cfi import *
#from HiggsAnalysis.VBFHiggsToZZto2l2b.Data.H130_ZZ_mumuqqFastSim_Input_cfi import *

# this inputs the input files from the previous function
process.source = RecoInput()

## set the number of events
process.maxEvents = cms.untracked.PSet(
#    input = cms.untracked.int32(-1)
    input = cms.untracked.int32(20)
)

## debugging porpose
cms.Service('Tracer')
## if you get bad_alloc problems or you just observe increasing memory needs
## [for the real big problems use valgrind (http://valgrind.org) or igtools]
cms.Service('SimpleMemoryCheck', 
    ignoreTotal = cms.untracked.int32(1)
) 
## quickly identification on how fast your process and single modules are 
cms.Service('Timing')

#from HiggsAnalysis.Skimming.vbfhzzllbb_HLTPaths_cfi import *


import copy
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
process.vbfhzzllbbHLTFilter = copy.deepcopy(hltHighLevel)
# provide list of HLT paths (or patterns) you want
process.vbfhzzllbbHLTFilter.HLTPaths = ['HLT_IsoMu11',
                                        'HLT_Mu15_L1Mu7', 
                                        'HLT_DoubleMu3', 
                                        'HLT_IsoEle15_L1I', 
                                        'HLT_IsoEle18_L1R', 
                                        'HLT_DoubleIsoEle10_L1I', 
                                        'HLT_DoubleIsoEle12_L1R'
                                        ]
# how to deal with multiple triggers:
# - True (OR) accept if ANY is true,
# - False (AND) accept if ALL are
process.vbfhzzllbbHLTFilter.andOr = cms.bool(True)

process.vbfhzzllbbHLTReport = cms.EDFilter("HLTrigReport",
    HLTriggerResults = cms.InputTag('TriggerResults',
                                    '',
                                    'HLT'
                                    )
)

process.vbfhzzllbbFilterEff = cms.EDAnalyzer("VBFHZZllbbSkimEff",
                                             VBFHZZllbbDebug = cms.bool(False),
                                             genParticleLabel = cms.untracked.InputTag('genParticles'),
                                             muonLabel        = cms.InputTag( 'muons'                  ),
                                             electronLabel    = cms.InputTag( 'pixelMatchGsfElectrons' ),
                                             tightMinimumPt = cms.double(10.),
                                             softMinimumPt  = cms.double(5.),
                                             tightLeptonMinimumNumber = cms.int32(1),
                                             softLeptonMinimumNumber  = cms.int32(2)
                                             )

process.p = cms.Path( process.vbfhzzllbbHLTReport
                      * process.vbfhzzllbbHLTFilter
                      * process.vbfhzzllbbFilterEff
                      )

## define output event selection to be that which satisfies 'p'
process.vbfHZZllbbEventSelection = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p')
    )
)

from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfHZZllbbEventContent_cff import *

## output module configuration
process.out = cms.OutputModule("PoolOutputModule",
    VBFHZZ2l2bEventContent,
    process.vbfHZZllbbEventSelection,
    fileName = cms.untracked.string('VBFHZZllbbFilterEff.root'),
    verbose  = cms.untracked.bool(False)
)

# extend event content to include objects from EDNtuple
process.out.outputCommands.extend(["keep *_*_*_VBFHZZllbbFILTEREFF"])

# define output path
process.outpath = cms.EndPath(process.out)

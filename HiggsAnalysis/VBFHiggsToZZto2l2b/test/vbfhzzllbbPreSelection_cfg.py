#!/usr/bin/env cmsRun
import FWCore.ParameterSet.Config as cms

process = cms.Process("VBFHZZllbbPRESELECTION")

process.load('Configuration/StandardSequences/Services_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
## to write every 10th event
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.categories.append('VBFHZZllbbPreSelectionSummary')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    default                   = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
    VBFHZZllbbPreSelectionSummary = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
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

process.load("HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbPreSelection_cff")
process.vbfhzzllbbMCprocessFilter.whichSim = cms.int32(1) 
#process.vbfhzzllbbPreSelection.tightLeptonMinNumber = cms.int32(0)
#process.vbfhzzllbbPreSelection.softLeptonMinNumber = cms.int32(0)
#process.vbfhzzllbbPreSelection.tightJetMinNumber = cms.int32(0)
#process.vbfhzzllbbPreSelection.softJetMinNumber = cms.int32(0)

## define output event selection to be that which satisfies 'p'
process.vbfHZZllbbEventSelection = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('path')
    )
)

from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfHZZllbbEventContent_cff import *

## output module configuration
process.out = cms.OutputModule("PoolOutputModule",
    VBFHZZ2l2bEventContent,
    process.vbfHZZllbbEventSelection,
    fileName = cms.untracked.string('VBFHZZllbbPreSelection.root'),
    verbose  = cms.untracked.bool(False)
)

# extend event content to include objects from EDNtuple
process.out.outputCommands.extend(["keep *_*_*_VBFHZZllbbPRESELECTION"])

# define output path
process.outpath = cms.EndPath(process.out)

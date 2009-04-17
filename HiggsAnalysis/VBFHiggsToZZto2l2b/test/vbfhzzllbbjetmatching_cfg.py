import FWCore.ParameterSet.Config as cms

process = cms.Process("VBFHZZllbbJetMatching")

process.load("FWCore.MessageService.MessageLogger_cfi")
## to write every 10th event
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.categories.append('VBFHZZllbbJetMatchingSummary')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    default                      = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
    VBFHZZllbbJetMatchingSummary = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True) # default values if False
)

## this defines the input files
from HiggsAnalysis.VBFHiggsToZZto2l2b.Data.H150_ZZ_qqllSummer08_IDEALV9v2_GENSIMRECO_FULL_Input_cfi import *
#from HiggsAnalysis.VBFHiggsToZZto2l2b.Data.H130_ZZ_mumuqqFastSim_Input_cfi import *

# this inputs the input files from the previous function
process.source = RecoInput()

process.maxEvents = cms.untracked.PSet(
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

## talk to TFileService for output histograms
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('VBFHZZllbbJetMatchingHistos.root')
)

process.load("HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbJetMatching_cff")

## define output event selection to be that which satisfies 'p'
process.vbfHZZllbbEventSelection = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('jetMatchinPath')
    )
)

from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfHZZllbbEventContent_cff import *

## output module configuration
process.out = cms.OutputModule("PoolOutputModule",
    VBFHZZ2l2bEventContent,
    process.vbfHZZllbbEventSelection,
    fileName = cms.untracked.string('VBFHZZllbbJetMatchingSkim.root'),
    verbose  = cms.untracked.bool(False)
)

# extend event content to include objects from EDNtuple
process.out.outputCommands.extend(["keep *_*_*_VBFHZZllbbJetMatching"])

# define output path
process.outpath = cms.EndPath(process.out)


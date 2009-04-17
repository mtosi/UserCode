import FWCore.ParameterSet.Config as cms

process = cms.Process("VBFHZZllbbCorJetWithBTagANALYZER")

process.load("FWCore.MessageService.MessageLogger_cfi")
## to write every 10th event
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.categories.append('VBFHZZllbbCorJetWithBTagANALYZERSummary')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    default                   = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
    VBFHZZllbbCorJetWithBTagANALYZERSummary = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True) # default values if False
)

## this defines the input files
process.source = cms.Source(
    "PoolSource",
    debugFlag = cms.untracked.bool(True),
    debugVebosity = cms.untracked.uint32(1),
    fileNames = cms.untracked.vstring(
        'file:VBFHZZllbbCorJetWithBTagTESTSkim.root'
    )
)

## set the number of events
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

process.vbfhzzllbbCorJetWithBTagAnalyzer = cms.EDAnalyzer('VBFHZZllbbCorJetWithBTagAnalyzer',
    corJetWithBTagLabel = cms.untracked.string('vbfHZZllbbCorJetWithBTagProd')
)

process.p = cms.Path(process.vbfhzzllbbCorJetWithBTagAnalyzer)

## define output event selection to be that which satisfies 'p'
process.vbfHZZllbbEventSelection = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p')
    )
)

from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfHZZllbbEventContent_cff import *
process.out = cms.OutputModule("PoolOutputModule",
    VBFHZZ2l2bEventContent,
    process.vbfHZZllbbEventSelection,
    fileName = cms.untracked.string('VBFHZZllbbANALYZER.root'),
    verbose  = cms.untracked.bool(False)
)

# define output path
process.outpath = cms.EndPath(process.out)

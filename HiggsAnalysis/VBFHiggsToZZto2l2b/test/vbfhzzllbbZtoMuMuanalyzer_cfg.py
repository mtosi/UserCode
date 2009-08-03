import FWCore.ParameterSet.Config as cms

process = cms.Process("VBFHZZllbbZtoMuMuAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
## to write every 10th event
process.MessageLogger.cerr.FwkReport.reportEvery = 10
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.cerr.threshold = 'ERROR'
process.MessageLogger.categories.append('VBFHZZllbbAnalysisSummary')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    default                   = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
    VBFHZZllbbAnalysisSummary = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True) # default values if False
)

## this defines the input files
from HiggsAnalysis.VBFHiggsToZZto2l2b.Data.PYTHIA6_SM_H_ZZ_qqll_mH150_10TeV_RECO_IDEAL_legnaro_cfi import *

# this inputs the input files from the previous function
process.source = RecoInput()

## set the number of events
process.maxEvents = cms.untracked.PSet(
#    input = cms.untracked.int32(-1)
    input = cms.untracked.int32(1000)
)

# Muon selection
process.load("HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbMuonSelector_cfi")
#process.vbfhzzllbbMuonSelector.sourceMinPtBarrelCut = cms.double(5.0)
#process.vbfhzzllbbMuonSelector.sourceMinPtEndcapCut = cms.double(3.0)
#process.vbfhzzllbbMuonSelector.sourceMinPEndcapCut  = cms.double(9.0)
#process.vbfhzzllbbMuonSelector.sourceMaxEtaCut      = cms.double(2.4)

# zToMuMu
process.load("HiggsAnalysis.VBFHiggsToZZto2l2b.zToMuMu_cfi")
process.zToMuMu.decay = cms.string('vbfhzzllbbMuonSelector@+ vbfhzzllbbMuonSelector@-')

process.p = cms.Path(process.vbfhzzllbbMuonSelector*
                     process.zToMuMu)

## define output event selection to be that which satisfies 'p'
process.vbfHZZllbbEventSelection = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p')
    )
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

### talk to TFileService for output histograms
#process.TFileService = cms.Service("TFileService",
#    fileName = cms.string('VBFHZZllbbZtoMuMuAnalysisHistos.root')
#)
from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfHZZllbbEventContent_cff import *

## output module configuration
process.out = cms.OutputModule("PoolOutputModule",
    VBFHZZ2l2bEventContent,
    process.vbfHZZllbbEventSelection,
    fileName = cms.untracked.string('VBFHZZllbbZtoMuMuAnalysisSkim.root'),
    verbose  = cms.untracked.bool(False),
)
# extend event content to include objects from EDNtuple
process.out.outputCommands.extend(["keep *_*_*_VBFHZZllbbZtoMuMuAnalysis"      ])

# define output path
process.outpath = cms.EndPath(process.out)


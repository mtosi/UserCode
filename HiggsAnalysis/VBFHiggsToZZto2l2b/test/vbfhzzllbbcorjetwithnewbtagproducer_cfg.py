import FWCore.ParameterSet.Config as cms

process = cms.Process("VBFHZZllbbCorJetWithBTagTEST")

process.load("FWCore.MessageService.MessageLogger_cfi")
## to write every 10th event
process.MessageLogger.cerr.FwkReport.reportEvery = 10
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.cerr.threshold = 'ERROR'
process.MessageLogger.categories.append('VBFHZZllbbCorJetWithBTagTESTSummary')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    default                   = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
    VBFHZZllbbCorJetWithBTagTESTSummary = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
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

## talk to TFileService for output histograms
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('VBFHZZllbbCorJetWithBTagTESTHistosNEW.root')
)

# from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbNewIterCone5PFBTagSequence_cff import *
# from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbNewSisCone5CaloBTagSequence_cff import *
# from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbcorjetwithbtagproducer_cfi import *

process.load("HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbNewIterCone5PFBTagSequence_cff")
process.load("HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbNewSisCone5CaloBTagSequence_cff")
process.load("HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbNewSisCone5PFBTagSequence_cff")
process.load("HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbcorjetwithbtagproducer_cfi")

process.pSequence = cms.Sequence(
    process.sisC5CaloBtaggingSequence *
    process.iC5PFBtaggingSequence *
    process.sisC5PFBtaggingSequence *
    process.vbfhzzllbbCorJetWithBTagSequence
)

process.p = cms.Path(
    process.pSequence
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
    fileName = cms.untracked.string('VBFHZZllbbCorJetWithBTagTESTSkimNEW.root'),
    verbose  = cms.untracked.bool(False)
)

# extend event content to include objects from EDNtuple
process.out.outputCommands.extend(["keep *_*_*_VBFHZZllbbCorJetWithBTagTEST"])
process.out.outputCommands.extend(["keep corJetWithBTagsDiscr_*_*_*"])
process.out.outputCommands.extend(["keep *_corJetWithBTagsDiscr_*_*"])
process.out.outputCommands.extend(["keep *_*_corJetWithBTagsDiscr_*"])
process.out.outputCommands.extend(["keep *_*_*_corJetWithBTagsDiscr"])
process.out.outputCommands.extend(["keep L2L3CorJetIC5Calo_*_*_*"])
process.out.outputCommands.extend(["keep *_L2L3CorJetIC5Calo_*_*"])
process.out.outputCommands.extend(["keep *_*_L2L3CorJetIC5Calo_*"])
process.out.outputCommands.extend(["keep *_*_*_L2L3CorJetIC5Calo"])

# define output path
process.outpath = cms.EndPath(process.out)

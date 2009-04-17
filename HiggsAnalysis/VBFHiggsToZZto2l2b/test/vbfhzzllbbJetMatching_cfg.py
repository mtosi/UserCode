import FWCore.ParameterSet.Config as cms

process = cms.Process("VBFHZZllbbAnalysis")

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
from PhysicsTools.PatAlgos.H150_ZZ_qqllSummer08_IDEALV9v2_GENSIMRECO_FULL_Input_cfi import *

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
    fileName = cms.string('VBFHZZllbbDeltaRAnalysisHistos.root')
)

# this is to get a CandidateCollection from a JetCollection
process.caloJetCollectionClone = cms.EDProducer("CaloJetShallowCloneProducer",
    src = cms.InputTag("iterativeCone5CaloJets")
)

process.genJetCollectionClone = cms.EDProducer("GenJetShallowCloneProducer",
    src = cms.InputTag("iterativeCone5GenJets")
)

# Selection of jets
process.caloJetSele = cms.EDFilter("PtMinCandSelector",
    src = cms.InputTag("caloJetCollectionClone"),
    ptMin = cms.double(20.0)
)

process.genJetSele = cms.EDFilter("PtMinCandSelector",
    src = cms.InputTag("genJetCollectionClone"),
    ptMin = cms.double(20.0)
)

# The matching routine OneToOne
process.jetMatchOne = cms.EDFilter("CandOneToOneDeltaRMatcher",
    src = cms.InputTag("iterativeCone5GenJets"),
    algoMethod = cms.string('SwitchMode'),
    matched = cms.InputTag("iterativeCone5CaloJets")
)

# The matching routine OneToMany
process.jetMatchMany = cms.EDFilter("CandOneToManyDeltaRMatcher",
    printDebug = cms.untracked.bool(True),
    src = cms.InputTag("genJetSele"),
    matched = cms.InputTag("caloJetSele")
)

process.printJet = cms.EDFilter("jetMatch",
    src = cms.InputTag("genJetSele"),
    matchMapMany = cms.InputTag("jetMatchMany"),
    HistOutFile = cms.untracked.string('myPlots.root'),
    matchMapOne = cms.InputTag("jetMatchOne","src2mtc"),
    matched = cms.InputTag("caloJetSele")
)

process.printEventNumber = cms.OutputModule("AsciiOutputModule")

process.p = cms.Path(process.caloJetCollectionClone*
                     process.genJetCollectionClone*
                     process.caloJetSele*
                     process.genJetSele*
                     process.jetMatchOne*
                     process.jetMatchMany*
                     process.printJet)

process.outpath = cms.EndPath(process.printEventNumber)






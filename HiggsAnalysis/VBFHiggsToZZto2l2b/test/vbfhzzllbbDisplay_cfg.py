import FWCore.ParameterSet.Config as cms

process = cms.Process("VBFHZZllbbDisplay")

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
    fileName = cms.string('VBFHZZllbbDisplayHisto.root')
)

#from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbanalyzer_cfi import *

process.vbfHZZllbbDisplay = cms.EDAnalyzer('VBFHZZllbbDisplay',
    genParticleLabel = cms.untracked.InputTag( 'genParticles'                    ),
    jetLabel         = cms.untracked.InputTag( 'iterativeCone5CaloJets'          ),
    muonLabel        = cms.untracked.InputTag( 'muons'                           ),
    electronLabel    = cms.untracked.InputTag( 'pixelMatchGsfElectrons'          ),
    tagJetLabel      = cms.untracked.InputTag( 'iterativeCone5CaloJets'          ),
    metLabel         = cms.untracked.InputTag( 'met'                             )
)

process.p = cms.Path(process.vbfHZZllbbDisplay)


import FWCore.ParameterSet.Config as cms

process = cms.Process("VBFHZZllbbRAWCORBTAGtest")

process.load("FWCore.MessageService.MessageLogger_cfi")
## to write every 10th event
process.MessageLogger.cerr.FwkReport.reportEvery = 10
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.cerr.threshold = 'ERROR'
process.MessageLogger.categories.append('VBFHZZllbbRAW_COR_BTAGtestingSummary')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    default                   = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
    VBFHZZllbbRAW_COR_BTAGtestingSummary = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
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
    fileName = cms.string('VBFHZZllbbRAWCORBTAGtestingHistos.root')
)

process.load("HiggsAnalysis.VBFHiggsToZZto2l2b.vbfHZZllbbIC5CaloCorrections_cff")

process.vbfHZZllbbTest = cms.EDAnalyzer('VBFHZZllbbRAWCORBTAGtesting',
    jetLabel      = cms.untracked.InputTag( 'iterativeCone5CaloJets'          ),
    corrJetLabel  = cms.untracked.InputTag( 'L2L3CorJetIC5Calo'               ),
    muonLabel     = cms.untracked.InputTag( 'muons'                           ),
    electronLabel = cms.untracked.InputTag( 'pixelMatchGsfElectrons'          ),
    metLabel      = cms.untracked.InputTag( 'met'                             ),
    corrMetLabel  = cms.untracked.InputTag( 'corMetType1Icone5Muons'          ),
#    bTagLabel     = cms.untracked.InputTag( 'combinedSecondaryVertexBJetTags' ),
    bTagConfig = cms.VPSet(
        cms.PSet(
#            bTagGenericAnalysisBlock,
#            label = cms.InputTag("combinedSecondaryVertexBJetTags")
            label = cms.InputTag("combinedSecondaryVertexBJetTags")
        ), 
        cms.PSet(
#            bTagTrackCountingAnalysisBlock,
            label = cms.InputTag("trackCountingHighEffBJetTags")
        ), 
        cms.PSet(
#            bTagTrackCountingAnalysisBlock,
            label = cms.InputTag("trackCountingHighPurBJetTags")
        ), 
        cms.PSet(
#            bTagProbabilityAnalysisBlock,
            label = cms.InputTag("jetProbabilityBJetTags")
#         ), 
#         cms.PSet(
# #            bTagBProbabilityAnalysisBlock,
#             label = cms.InputTag("jetBProbabilityBJetTags")
#         ), 
#         cms.PSet(
# #            bTagSimpleSVAnalysisBlock,
#             label = cms.InputTag("simpleSecondaryVertexBJetTags")
#         ), 
#         cms.PSet(
# #            bTagGenericAnalysisBlock,
#             label = cms.InputTag("combinedSecondaryVertexMVABJetTags")
        )
    ),                                        
     jetCorrectionService = cms.string('L2L3JetCorrectorIC5Calo')                                        
 )
#    metLabel         = cms.untracked.InputTag( 'corMetType1Icone5'               ),

process.p = cms.Path(
    process.ic5CaloJetMETCorrections *
    process.vbfHZZllbbTest
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
    fileName = cms.untracked.string('VBFHZZllbbRAWCORBTAGtestingSkim.root'),
    verbose  = cms.untracked.bool(False)
)

# extend event content to include objects from EDNtuple
process.out.outputCommands.extend(["keep *_*_*_VBFHZZllbbTest"])

# define output path
process.outpath = cms.EndPath(process.out)


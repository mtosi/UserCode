import FWCore.ParameterSet.Config as cms

process = cms.Process("VBFHZZllbbBTAGtest")

process.load("FWCore.MessageService.MessageLogger_cfi")
## to write every 10th event
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.categories.append('VBFHZZllbbBTAGtestSummary')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    default                   = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
    VBFHZZllbbBTAGtestSummary = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
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
    fileName = cms.string('VBFHZZllbbBTAGtestingHistos.root')
)

process.load("HiggsAnalysis.VBFHiggsToZZto2l2b.vbfHZZllbbIC5CaloCorrections_cff")

# default configuration with frontier conditions
process.load("Configuration.StandardSequences.Reconstruction_cff")

# b-tagging general configuration
process.load("RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi")
process.load("RecoBTag.Configuration.RecoBTag_cff")

#    caloJetLabel     = cms.untracked.InputTag( 'iterativeCone5CaloJets'          ),
#    corrCaloJetLabel = cms.untracked.InputTag( 'L2L3CorJetIC5Calo'               ),
#    bTagLabel        = cms.untracked.InputTag( 'combinedSecondaryVertexBJetTags' ),
#    muonLabel        = cms.untracked.InputTag( 'muons'                           ),
#    electronLabel    = cms.untracked.InputTag( 'pixelMatchGsfElectrons'          ),
#    metLabel         = cms.untracked.InputTag( 'met'                             )
#    metLabel         = cms.untracked.InputTag( 'corMetType1Icone5'               ),

process.p = cms.Path(
    process.ic5CaloJetMETCorrections   *

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
    fileName = cms.untracked.string('VBFHZZllbbBTAGtestingSkim.root'),
    verbose  = cms.untracked.bool(False)
)

# extend event content to include objects from EDNtuple
process.out.outputCommands.extend(["keep *_*_*_VBFHZZllbbBTAGtest"])

# define output path
process.outpath = cms.EndPath(process.out)


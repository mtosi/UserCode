import FWCore.ParameterSet.Config as cms

process = cms.Process("VBFHZZllbbJetMatching")

process.load("FWCore.MessageService.MessageLogger_cfi")
## to write every 10th event
process.MessageLogger.cerr.FwkReport.reportEvery = 10
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.cerr.threshold = 'ERROR'
process.MessageLogger.categories.append('VBFHZZllbbJetMatchingSummary')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    default                      = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
    VBFHZZllbbJetMatchingSummary = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True) # default values if False
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)


## this defines the input files
from PhysicsTools.PatAlgos.H150_ZZ_qqllSummer08_IDEALV9v2_GENSIMRECO_FULL_Input_cfi import *
#from PhysicsTools.PatAlgos.PYTHIA6_SM_H_ZZ_2l_2jets_mH800_10TeV_GEN_FASTSIM_Input_cfi import *
# this inputs the input files from the previous function
process.source = RecoInput()

## talk to TFileService for output histograms
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('VBFHZZllbbJetMatchingHistos.root')
)

process.load("HiggsAnalysis.VBFHiggsToZZto2l2b.vbfHZZllbbIC5CaloCorrections_cff")

process.vbfHZZllbbJetMatching = cms.EDAnalyzer('VBFHZZllbbJetMatching',
    genParticleLabel = cms.untracked.InputTag( 'genParticles'                    ),
#    jetLabel         = cms.untracked.InputTag( 'L2L3CorJetIC5Calo'               ),
    jetLabel          = cms.untracked.InputTag( 'iterativeCone5CaloJets'          ),
    muonLabel        = cms.untracked.InputTag( 'muons'                           ),
    electronLabel    = cms.untracked.InputTag( 'pixelMatchGsfElectrons'          ),
    bTagLabel        = cms.untracked.InputTag( 'combinedSecondaryVertexBJetTags' ),
#    tauLabel         = cms.untracked.InputTag( '' ),
#    metLabel         = cms.untracked.InputTag( 'corMetType1Icone5Muons'          ),
#    metLabel         = cms.untracked.InputTag( 'corMetType1Icone5'               ),
    metLabel         = cms.untracked.InputTag( 'met'                             ),
    jetEtCut = cms.untracked.double(10)                           
)

process.p = cms.Path(
#    process.ic5CaloJetMETCorrections_withoutMuonCorr *
#    process.ic5CaloJetMETCorrections *
    process.vbfHZZllbbJetMatching
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
    fileName = cms.untracked.string('VBFHZZllbbJetMatchingSkim.root'),
    verbose  = cms.untracked.bool(False)
)

# extend event content to include objects from EDNtuple
process.out.outputCommands.extend(["keep *_*_*_VBFHZZllbbJetMatching"])

# define output path
process.outpath = cms.EndPath(process.out)


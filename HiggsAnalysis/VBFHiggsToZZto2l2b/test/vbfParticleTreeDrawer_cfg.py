import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

## ## this defines the input files
## from PhysicsTools.PatAlgos.H150_ZZ_qqllSummer08_IDEALV9v2_GENSIMRECO_FULL_Input_cfi import *
## #from PhysicsTools.PatAlgos.PYTHIA6_SM_H_ZZ_2l_2jets_mH800_10TeV_GEN_FASTSIM_Input_cfi import *
## # this inputs the input files from the previous function
## process.source = RecoInput()

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/data/tosi/VQQ-madgraph/VQQ-madgraphFile/VQQ-madgraphFall08_IDEAL_GEN-SIM-RECO_1.root',
        'file:/data/tosi/VQQ-madgraph/VQQ-madgraphFile/VQQ-madgraphFall08_IDEAL_GEN-SIM-RECO_2.root',
        'file:/data/tosi/VQQ-madgraph/VQQ-madgraphFile/VQQ-madgraphFall08_IDEAL_GEN-SIM-RECO_3.root',
        'file:/data/tosi/VQQ-madgraph/VQQ-madgraphFile/VQQ-madgraphFall08_IDEAL_GEN-SIM-RECO_4.root',
        'file:/data/tosi/VQQ-madgraph/VQQ-madgraphFile/VQQ-madgraphFall08_IDEAL_GEN-SIM-RECO_5.root',
        'file:/data/tosi/VQQ-madgraph/VQQ-madgraphFile/VQQ-madgraphFall08_IDEAL_GEN-SIM-RECO_6.root',
        'file:/data/tosi/VQQ-madgraph/VQQ-madgraphFile/VQQ-madgraphFall08_IDEAL_GEN-SIM-RECO_7.root',
        'file:/data/tosi/VQQ-madgraph/VQQ-madgraphFile/VQQ-madgraphFall08_IDEAL_GEN-SIM-RECO_8.root',
        'file:/data/tosi/VQQ-madgraph/VQQ-madgraphFile/VQQ-madgraphFall08_IDEAL_GEN-SIM-RECO_9.root',
        'file:/data/tosi/VQQ-madgraph/VQQ-madgraphFile/VQQ-madgraphFall08_IDEAL_GEN-SIM-RECO_10.root'
    )
)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.printList = cms.EDFilter("ParticleListDrawer",
    src              = cms.InputTag("genParticles"),
    maxEventsToPrint = cms.untracked.int32(1)
)

process.printTree = cms.EDFilter("ParticleTreeDrawer",
    src           = cms.InputTag("genParticles"),
    printP4       = cms.untracked.bool(False),
    printPtEtaPhi = cms.untracked.bool(True),
    printVertex   = cms.untracked.bool(False),
    printStatus   = cms.untracked.bool(False),
    printIndex    = cms.untracked.bool(True),
    status        = cms.untracked.vint32(1,3),
)

process.printDecay = cms.EDAnalyzer("ParticleDecayDrawer",
    src           = cms.InputTag("genParticles"),
    printP4       = cms.untracked.bool(False),
    printPtEtaPhi = cms.untracked.bool(False),
    printVertex   = cms.untracked.bool(True),
)

process.processAnalyzer = cms.EDAnalyzer("VBFHZZllbbMCprocessAnalyzer",
   SIGNAL = cms.untracked.bool(True)
)

process.p = cms.Path(
#    process.printList *
    process.printTree
#    * process.printDecay
    * process.processAnalyzer
)


import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

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

process.demo = cms.EDAnalyzer('VBFHZZllbbMCvalidation',
                              OutputName  = cms.untracked.string('VQQ-madgraphFall08_IDEAL_GEN-SIM-RECOOutput.root'),
                              MCParticles = cms.untracked.string('genParticles'),
                              SIGNAL      = cms.untracked.bool(False)
                              )

process.p = cms.Path(process.demo)

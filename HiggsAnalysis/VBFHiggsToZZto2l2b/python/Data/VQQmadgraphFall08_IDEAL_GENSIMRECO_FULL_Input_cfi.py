## import configurations
import FWCore.ParameterSet.Config as cms

# from 

def RecoInput() : 
 return cms.Source("PoolSource",
                   debugVerbosity = cms.untracked.uint32(200),
                   debugFlag = cms.untracked.bool(True),
                   fileNames = cms.untracked.vstring(
                        'file:/data/tosi/VQQ-madgraph/VQQ-madgraphFile/VQQ-madgraphFall08_IDEAL_GEN-SIM-RECO_FULL1.root',
                        'file:/data/tosi/VQQ-madgraph/VQQ-madgraphFile/VQQ-madgraphFall08_IDEAL_GEN-SIM-RECO_FULL2.root',
                        'file:/data/tosi/VQQ-madgraph/VQQ-madgraphFile/VQQ-madgraphFall08_IDEAL_GEN-SIM-RECO_FULL3.root',
                        'file:/data/tosi/VQQ-madgraph/VQQ-madgraphFile/VQQ-madgraphFall08_IDEAL_GEN-SIM-RECO_FULL4.root',
                        'file:/data/tosi/VQQ-madgraph/VQQ-madgraphFile/VQQ-madgraphFall08_IDEAL_GEN-SIM-RECO_FULL5.root',
                        'file:/data/tosi/VQQ-madgraph/VQQ-madgraphFile/VQQ-madgraphFall08_IDEAL_GEN-SIM-RECO_FULL6.root',
                        'file:/data/tosi/VQQ-madgraph/VQQ-madgraphFile/VQQ-madgraphFall08_IDEAL_GEN-SIM-RECO_FULL7.root',
                        'file:/data/tosi/VQQ-madgraph/VQQ-madgraphFile/VQQ-madgraphFall08_IDEAL_GEN-SIM-RECO_FULL8.root',
                        'file:/data/tosi/VQQ-madgraph/VQQ-madgraphFile/VQQ-madgraphFall08_IDEAL_GEN-SIM-RECO_FULL9.root',
                        'file:/data/tosi/VQQ-madgraph/VQQ-madgraphFile/VQQ-madgraphFall08_IDEAL_GEN-SIM-RECO_FULL10.root'
                    )
                   )

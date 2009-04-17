import FWCore.ParameterSet.Config as cms

vbfhzzllbbJetMatching = cms.EDFilter("VBFHZZllbbJetMatching",
    # reconstructed objects
     electronLabel        = cms.InputTag('pixelMatchGsfElectrons'),
     muonLabel            = cms.InputTag('muons'),                 
     metLabel             = cms.InputTag('met'),
#     metLabel             = cms.InputTag('corMetType1Icone5Muons'),
     jetLabel             = cms.InputTag('iterativeCone5CaloJets'),
     corJetsWithBTagLabel = cms.string('vbfhzzllbbCorJetWithBTagProd'),
     mcParticleLabel      = cms.InputTag('genParticles'),
     genJetLabel          = cms.InputTag('iterativeCone5GenJets'),
     genMetLabel          = cms.InputTag('genMetNoNuBSM')
)



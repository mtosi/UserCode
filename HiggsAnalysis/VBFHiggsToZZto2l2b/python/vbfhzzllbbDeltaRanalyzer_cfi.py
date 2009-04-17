import FWCore.ParameterSet.Config as cms

vbfHZZllbbDeltaRAnalyzer = cms.EDAnalyzer('VBFHZZllbbDeltaRAnalyzer',
    genParticleLabel = cms.untracked.InputTag( 'genParticles'                    ),
    caloJetLabel     = cms.untracked.InputTag( 'iterativeCone5CaloJets'          ),
    muonLabel        = cms.untracked.InputTag( 'muons'                           ),
    electronLabel    = cms.untracked.InputTag( 'pixelMatchGsfElectrons'          ),
    bTagLabel        = cms.untracked.InputTag( 'combinedSecondaryVertexBJetTags' ),
    metLabel         = cms.untracked.InputTag( 'caloMET'                         ), # needs to be checked!!!
    ptMax = cms.untracked.double(10)                           
)

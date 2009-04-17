import FWCore.ParameterSet.Config as cms

vbfhzzllbbPreSelection = cms.EDFilter("VBFHZZllbbPreSelection",
    # reconstructed objects
    electronLabel        = cms.InputTag ( 'pixelMatchGsfElectrons'      ),
    muonLabel            = cms.InputTag ( 'muons'                       ),
    metLabel             = cms.InputTag ( 'met'                         ),
#    metLabel             = cms.InputTag ( 'corMetType1Icone5Muons'      ),
    jetLabel             = cms.InputTag ( 'iterativeCone5CaloJets'      ),
    corJetsWithBTagLabel = cms.string   ( 'vbfhzzllbbCorJetWithBTagProd'),
    mcParticleLabel      = cms.InputTag ( 'genParticles'                ),
    genJetLabel          = cms.InputTag ( 'iterativeCone5GenJets'       ),
    genMetLabel          = cms.InputTag ( 'genMetNoNuBSM'               ),
    # minimum number of identified leptons above pt threshold
    tightLeptonMinNumber = cms.int32(1),
    softLeptonMinNumber  = cms.int32(2),
    # minimum number of identified jets above pt threshold
    tightJetMinNumber = cms.int32(3),
    softJetMinNumber  = cms.int32(4),
    # pt threshold for leptons
    tightLeptonMinPt = cms.double(10.0),
    softLeptonMinPt  = cms.double( 5.0),
    # pt threshold for jets
    tightJetMinPt = cms.double(20.0),
    softJetMinPt  = cms.double(15.0)
)



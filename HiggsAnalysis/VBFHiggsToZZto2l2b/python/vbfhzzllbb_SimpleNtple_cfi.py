import FWCore.ParameterSet.Config as cms


vbfhzzllbbSimpleNtple = cms.EDAnalyzer("SimpleNtple",
        whichSim = cms.int32(1),                               
	trackLabel    = cms.InputTag('generalTracks'),
	muonLabel     = cms.InputTag('muons'),
	electronLabel = cms.InputTag('pixelMatchGsfElectrons'),
	metLabel      = cms.InputTag('met'),
# 	metLabel      = cms.InputTag('corMetType1Icone5Muons'),
	tagJetLabel   = cms.InputTag('iterativeCone5CaloJets'),
        corIC5CaloJetsWithBTagLabel = cms.string('vbfhzzllbbCorJetWithBTagProd'),
#        corIC5CaloJetsWithBTagLabel = cms.string('iterativeCone5CaloJets'),
        corIC5PFJetsWithBTagLabel   = cms.string(''),  
        corIC5PFJetsWithBTagFlag = cms.bool(False),
	genParticleLabel = cms.InputTag('genParticles'),
	genJetLabel      = cms.InputTag('iterativeCone5GenJets'),
	genMetLabel      = cms.InputTag('genMet')
)


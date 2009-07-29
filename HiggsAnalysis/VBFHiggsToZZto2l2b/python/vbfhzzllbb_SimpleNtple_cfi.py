import FWCore.ParameterSet.Config as cms


vbfhzzllbbSimpleNtple = cms.EDAnalyzer("SimpleNtple",
        whichSim = cms.int32(1), # 0:FastSim, 1:FullSim                              
        vertexLabel   = cms.InputTag('offlinePrimaryVertexWithBS'),
	trackLabel    = cms.InputTag('generalTracks'),
	muonLabel     = cms.InputTag('muons'),
	electronLabel = cms.InputTag('pixelMatchGsfElectrons'),
#	metLabel      = cms.InputTag('met'),
 	metLabel      = cms.InputTag('corMetType1Icone5Muons'),
	tagJetLabel   = cms.InputTag('iterativeCone5CaloJets'),
        corIC5CaloJetsWithBTagLabel = cms.string('vbfhzzllbbCorJetWithBTagProd'),
#        corIC5CaloJetsWithBTagLabel = cms.string('iterativeCone5CaloJets'),
        corIC5PFJetsWithBTagLabel   = cms.string(''),  
        corIC5PFJetsWithBTagFlag = cms.bool(False),
	genParticleLabel  = cms.InputTag('genParticles'),
	genJetLabel       = cms.InputTag('iterativeCone5GenJets'),
	genMetLabel       = cms.InputTag('genMet'),
        eleTrkIsoAlgoFlag = cms.bool(True),
           coneRadius      = cms.double(0.3),
           vetoRadius      = cms.double(0.1),
           otherVetoRadius = cms.double(0.0002),
           ptMin           = cms.double(100.),
           lipMax          = cms.double(100.),
           useTkQuality    = cms.untracked.bool(True)

)


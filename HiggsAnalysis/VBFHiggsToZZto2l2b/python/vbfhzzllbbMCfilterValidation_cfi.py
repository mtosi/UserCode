import FWCore.ParameterSet.Config as cms

vbfhzzllbbMCfilterValidation = cms.EDFilter("VBFHZZllbbMCfilterValidation",
        whichSim = cms.int32(1),
        signal   = cms.int32(1),
        genJetLabel      = cms.InputTag('iterativeCone5GenJets'),
        genParticleLabel = cms.InputTag('genParticles'),
        jetNumberCut       = cms.int32(2), 
        firstJetPtCut      = cms.double(20.),
        secondJetPtCut     = cms.double(15.),
        jetpairInvMassCut  = cms.double(300.),
        jetPairDeltaEtaCut = cms.double(1.),
        leptonEtaCut   = cms.double(2.7),
        leptonPtCut    = cms.double(9.)
)



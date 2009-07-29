import FWCore.ParameterSet.Config as cms

vbfhzzllbbCorJetWithBTagProd = cms.EDProducer('VBFHZZllbbCorJetWithBTagProducer',
    bTagConfig = cms.VPSet(
        cms.PSet(
            label = cms.InputTag("combinedSecondaryVertexBJetTags")
        ), 
        cms.PSet(
            label = cms.InputTag("trackCountingHighEffBJetTags")
        ), 
        cms.PSet(
            label = cms.InputTag("trackCountingHighPurBJetTags")
        ), 
        cms.PSet(
            label = cms.InputTag("jetProbabilityBJetTags")
        )
    ),                                        
    jetCorrectionService = cms.string('L2L3JetCorrectorIC5Calo')                                        
)



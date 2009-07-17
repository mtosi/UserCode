import FWCore.ParameterSet.Config as cms

# Full Event content 
VBFHZZ2l2bEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring(
#        'keep *'
        'drop *',

        'keep edmTriggerResults_*_*_*',
        'keep edmHepMCProduct_*_*_*',
        'keep *_genEventProcID_*_*',
        
        # MC collections
        'keep edmGenInfoProduct_*_*_*',
        'keep *_genMetNoNuBSM_*_*',
        'keep *_genParticles_*_*',
        'keep *_source_*_*',

        # muon collections
        'keep *_muons_*_*',
        'keep *_globalMuons_*_*',  

        # electron collections
        'keep *_pixelMatchGsfElectrons_*_*',

        # jet collections
        'keep *_iterativeCone5GenJets_*_*',
        'keep *_kt4GenJets_*_*',
        'keep *_kt6GenJets_*_*',
        'keep *_sisCone5GenJets_*_*',
        'keep *_sisCone7GenJets_*_*',        
        'keep *_iterativeCone5CaloJets_*_*',
        'keep *_kt4CaloJets_*_*',
        'keep *_kt6CaloJets_*_*',
        'keep *_sisCone5CaloJets_*_*',
        'keep *_sisCone7CaloJets_*_*',
        
        # particle-flow collections
#        'keep *recoPF*_*_*_*',
        'keep *_particleFlow_*_*',
        'keep *_iterativeCone5PFJets_*_*',
        'keep *_kt4PFJets_*_*',
        'keep *_kt6PFJets_*_*',
        'keep *_sisCone5PFJets_*_*',
        'keep *_sisCone7PFJets_*_*',
#        'keep *_pfMet_*_*',

        # met collections
        'keep *_genMet_*_*',
        'keep *_met_*_*',

        # b-tagging JetTags collections
        'keep *_simpleSecondaryVertexBJetTags_*_*',
        'keep *_combinedSecondaryVertexBJetTags_*_*',
        'keep *_combinedSecondaryVertexMVABJetTags_*_*',
        'keep *_impactParameterMVABJetTags_*_*',
        'keep *_jetProbabilityBJetTags_*_*',
        'keep *_jetBProbabilityBJetTags_*_*',
        'keep *_trackCountingHighEffBJetTags_*_*',
        'keep *_trackCountingHighPurBJetTags_*_*',
        # b-tagging TagInfos collections
        'keep *_impactParameterTagInfos_*_*',
        'keep *_secondaryVertexTagInfos_*_*',

        # vertices collection
        'keep *_offlinePrimaryVertices_*_*',
        'keep *_offlinePrimaryVerticesWithBS_*_*',
        'keep *_pixelVertices_*_*',


        'keep *_offlineBeamSpot_*_*',

        'keep recoGsfTracks_*_*_*',
        'keep recoGsfTrackExtras_*_*_*',
        'keep *_pixelMatchGsfFit_*_*',                       
        'keep recoTracks_*_*_*',
        'keep recoTrackExtras_*_*_*',             

        'keep recoBasicClusters_*_*_*',
        'keep recoSuperClusters_*_*_*',
        'keep *_towerMaker_*_*',
        
        'keep *_hbhereco_*_*',
        'keep *_hfreco_*_*',
        'keep *_horeco_*_*', 
        'keep *_ecalRecHit_*_*',
        'keep EcalRecHitsSorted_*_*_*',  

        'keep *_generalTracks_*_*',  

        'keep *_eleIsoFromDepsEcalFromHits_*_*',
        'keep *_eleIsoFromDepsHcalFromHits_*_*',
        'keep *_eleIsoFromDepsTk_*_*',
        
        'keep *_eleIsoDepositEcalFromHits_*_*',
        'keep *_eleIsoDepositHcalFromHits_*_*',
        'keep *_eleIsoDepositTk_*_*',
       
        'keep *_eidRobustLoose_*_*',
        'keep *_eidRobustTight_*_*',
        'keep *_eidLoose_*_*',
        'keep *_eidTight_*_*',
       
        'keep recoIsoDepositedmValueMap_*_*_*'
              
        )
    )

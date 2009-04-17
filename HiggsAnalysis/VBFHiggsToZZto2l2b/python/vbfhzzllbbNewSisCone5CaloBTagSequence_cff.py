####################################################################################
## In order to rerun all the algorithms on a different jet (or track) collection, ##
## one needs to duplicate the original b-tag sequence.                            ##
## Only the EDProducers needs to be cloned.                                       ##
####################################################################################

import FWCore.ParameterSet.Config as cms

# default configuration with frontier conditions
from Configuration.StandardSequences.MagneticField_cff import *
from Configuration.StandardSequences.Geometry_cff import *
from Configuration.StandardSequences.Reconstruction_cff import *
from Configuration.StandardSequences.FrontierConditions_GlobalTag_cff import *
GlobalTag.globaltag = 'IDEAL_V9::All'

# b-tagging general configuration
from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import *

# create a new jets and tracks association
sisC5CaloJetTracksAssociatorAtVertex = ic5JetTracksAssociatorAtVertex.clone()
sisC5CaloJetTracksAssociatorAtVertex.jets = "sisCone5CaloJets"
sisC5CaloJetTracksAssociatorAtVertex.tracks = "generalTracks"


# The one needs to clone the b-tag producers and instruct them to use this new collection
from RecoBTag.Configuration.RecoBTag_cff import *

# impact parameter b-tag
sisC5CaloImpactParameterTagInfos = impactParameterTagInfos.clone()
sisC5CaloImpactParameterTagInfos.jetTracks = "sisC5CaloJetTracksAssociatorAtVertex"
sisC5CaloTrackCountingHighEffBJetTags = trackCountingHighEffBJetTags.clone()
sisC5CaloTrackCountingHighEffBJetTags.tagInfos = cms.VInputTag( cms.InputTag("sisC5CaloImpactParameterTagInfos") )
sisC5CaloTrackCountingHighPurBJetTags = trackCountingHighPurBJetTags.clone()
sisC5CaloTrackCountingHighPurBJetTags.tagInfos = cms.VInputTag( cms.InputTag("sisC5CaloImpactParameterTagInfos") )
sisC5CaloJetProbabilityBJetTags = jetProbabilityBJetTags.clone()
sisC5CaloJetProbabilityBJetTags.tagInfos = cms.VInputTag( cms.InputTag("sisC5CaloImpactParameterTagInfos") )
sisC5CaloJetBProbabilityBJetTags = jetBProbabilityBJetTags.clone()
sisC5CaloJetBProbabilityBJetTags.tagInfos = cms.VInputTag( cms.InputTag("sisC5CaloImpactParameterTagInfos") )

# secondary vertex b-tag
sisC5CaloSecondaryVertexTagInfos = secondaryVertexTagInfos.clone()
sisC5CaloSecondaryVertexTagInfos.trackIPTagInfos = "sisC5CaloImpactParameterTagInfos"
sisC5CaloSimpleSecondaryVertexBJetTags = simpleSecondaryVertexBJetTags.clone()
sisC5CaloSimpleSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("sisC5CaloSecondaryVertexTagInfos") )
sisC5CaloCombinedSecondaryVertexBJetTags = combinedSecondaryVertexBJetTags.clone()
sisC5CaloCombinedSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("sisC5CaloImpactParameterTagInfos"), cms.InputTag("sisC5CaloSecondaryVertexTagInfos") )
sisC5CaloCombinedSecondaryVertexMVABJetTags = combinedSecondaryVertexMVABJetTags.clone()
sisC5CaloCombinedSecondaryVertexMVABJetTags.tagInfos = cms.VInputTag( cms.InputTag("sisC5CaloImpactParameterTagInfos"), cms.InputTag("sisC5CaloSecondaryVertexTagInfos") )

# soft electron b-tag
sisC5CaloSoftElectronTagInfos = softElectronTagInfos.clone()
sisC5CaloSoftElectronTagInfos.jets = "sisCone5CaloJets"
sisC5CaloSoftElectronBJetTags = softElectronBJetTags.clone()
sisC5CaloSoftElectronBJetTags.tagInfos = cms.VInputTag( cms.InputTag("sisC5CaloSoftElectronTagInfos") )

# soft muon b-tag
sisC5CaloSoftMuonTagInfos = softMuonTagInfos.clone()
sisC5CaloSoftMuonTagInfos.jets = "sisCone5CaloJets"
sisC5CaloSoftMuonBJetTags = softMuonBJetTags.clone()
sisC5CaloSoftMuonBJetTags.tagInfos = cms.VInputTag( cms.InputTag("sisC5CaloSoftMuonTagInfos") )
sisC5CaloSoftMuonByIP3dBJetTags = softMuonByIP3dBJetTags.clone()
sisC5CaloSoftMuonByIP3dBJetTags.tagInfos = cms.VInputTag( cms.InputTag("sisC5CaloSoftMuonTagInfos") )
sisC5CaloSoftMuonByPtBJetTags = softMuonByPtBJetTags.clone()
sisC5CaloSoftMuonByPtBJetTags.tagInfos = cms.VInputTag( cms.InputTag("sisC5CaloSoftMuonTagInfos") )


# prepare a path running the new modules
sisC5CaloJetTracksAssociator = cms.Sequence(
    sisC5CaloJetTracksAssociatorAtVertex
)

sisC5CaloJetBtaggingIP = cms.Sequence(
    sisC5CaloImpactParameterTagInfos * (
        sisC5CaloTrackCountingHighEffBJetTags +
        sisC5CaloTrackCountingHighPurBJetTags +
        sisC5CaloJetProbabilityBJetTags +
        sisC5CaloJetBProbabilityBJetTags
    )
)

sisC5CaloJetBtaggingSV = cms.Sequence(
    sisC5CaloImpactParameterTagInfos *
    sisC5CaloSecondaryVertexTagInfos * (
        sisC5CaloSimpleSecondaryVertexBJetTags +
        sisC5CaloCombinedSecondaryVertexBJetTags +
        sisC5CaloCombinedSecondaryVertexMVABJetTags
    )
)

sisC5CaloJetBtaggingEle = cms.Sequence(
    btagSoftElectrons *
    sisC5CaloSoftElectronTagInfos *
    sisC5CaloSoftElectronBJetTags
)

sisC5CaloJetBtaggingMu = cms.Sequence(
    sisC5CaloSoftMuonTagInfos * (
        sisC5CaloSoftMuonBJetTags +
        sisC5CaloSoftMuonByIP3dBJetTags +
        sisC5CaloSoftMuonByPtBJetTags
    )
)

sisC5CaloJetBtagging = cms.Sequence(
    sisC5CaloJetBtaggingIP +
    sisC5CaloJetBtaggingSV +
    sisC5CaloJetBtaggingEle +
    sisC5CaloJetBtaggingMu
)

sisC5CaloBtaggingSequence = cms.Sequence(
    sisC5CaloJetTracksAssociator *
    sisC5CaloJetBtagging
)

sisC5CaloBtaggingPath = cms.Path(
    sisC5CaloBtaggingSequence
)

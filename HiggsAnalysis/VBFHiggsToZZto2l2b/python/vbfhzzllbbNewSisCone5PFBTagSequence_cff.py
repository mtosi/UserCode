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
sisC5PFJetTracksAssociatorAtVertex = ic5JetTracksAssociatorAtVertex.clone()
sisC5PFJetTracksAssociatorAtVertex.jets = "sisCone5PFJets"
sisC5PFJetTracksAssociatorAtVertex.tracks = "generalTracks"


# The one needs to clone the b-tag producers and instruct them to use this new collection
from RecoBTag.Configuration.RecoBTag_cff import *

# impact parameter b-tag
sisC5PFImpactParameterTagInfos = impactParameterTagInfos.clone()
sisC5PFImpactParameterTagInfos.jetTracks = "sisC5PFJetTracksAssociatorAtVertex"
sisC5PFTrackCountingHighEffBJetTags = trackCountingHighEffBJetTags.clone()
sisC5PFTrackCountingHighEffBJetTags.tagInfos = cms.VInputTag( cms.InputTag("sisC5PFImpactParameterTagInfos") )
sisC5PFTrackCountingHighPurBJetTags = trackCountingHighPurBJetTags.clone()
sisC5PFTrackCountingHighPurBJetTags.tagInfos = cms.VInputTag( cms.InputTag("sisC5PFImpactParameterTagInfos") )
sisC5PFJetProbabilityBJetTags = jetProbabilityBJetTags.clone()
sisC5PFJetProbabilityBJetTags.tagInfos = cms.VInputTag( cms.InputTag("sisC5PFImpactParameterTagInfos") )
sisC5PFJetBProbabilityBJetTags = jetBProbabilityBJetTags.clone()
sisC5PFJetBProbabilityBJetTags.tagInfos = cms.VInputTag( cms.InputTag("sisC5PFImpactParameterTagInfos") )

# secondary vertex b-tag
sisC5PFSecondaryVertexTagInfos = secondaryVertexTagInfos.clone()
sisC5PFSecondaryVertexTagInfos.trackIPTagInfos = "sisC5PFImpactParameterTagInfos"
sisC5PFSimpleSecondaryVertexBJetTags = simpleSecondaryVertexBJetTags.clone()
sisC5PFSimpleSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("sisC5PFSecondaryVertexTagInfos") )
sisC5PFCombinedSecondaryVertexBJetTags = combinedSecondaryVertexBJetTags.clone()
sisC5PFCombinedSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("sisC5PFImpactParameterTagInfos"), cms.InputTag("sisC5PFSecondaryVertexTagInfos") )
sisC5PFCombinedSecondaryVertexMVABJetTags = combinedSecondaryVertexMVABJetTags.clone()
sisC5PFCombinedSecondaryVertexMVABJetTags.tagInfos = cms.VInputTag( cms.InputTag("sisC5PFImpactParameterTagInfos"), cms.InputTag("sisC5PFSecondaryVertexTagInfos") )

# soft electron b-tag
sisC5PFSoftElectronTagInfos = softElectronTagInfos.clone()
sisC5PFSoftElectronTagInfos.jets = "sisCone5PFJets"
sisC5PFSoftElectronBJetTags = softElectronBJetTags.clone()
sisC5PFSoftElectronBJetTags.tagInfos = cms.VInputTag( cms.InputTag("sisC5PFSoftElectronTagInfos") )

# soft muon b-tag
sisC5PFSoftMuonTagInfos = softMuonTagInfos.clone()
sisC5PFSoftMuonTagInfos.jets = "sisCone5PFJets"
sisC5PFSoftMuonBJetTags = softMuonBJetTags.clone()
sisC5PFSoftMuonBJetTags.tagInfos = cms.VInputTag( cms.InputTag("sisC5PFSoftMuonTagInfos") )
sisC5PFSoftMuonByIP3dBJetTags = softMuonByIP3dBJetTags.clone()
sisC5PFSoftMuonByIP3dBJetTags.tagInfos = cms.VInputTag( cms.InputTag("sisC5PFSoftMuonTagInfos") )
sisC5PFSoftMuonByPtBJetTags = softMuonByPtBJetTags.clone()
sisC5PFSoftMuonByPtBJetTags.tagInfos = cms.VInputTag( cms.InputTag("sisC5PFSoftMuonTagInfos") )


# prepare a path running the new modules
sisC5PFJetTracksAssociator = cms.Sequence(
    sisC5PFJetTracksAssociatorAtVertex
)

sisC5PFJetBtaggingIP = cms.Sequence(
    sisC5PFImpactParameterTagInfos * (
        sisC5PFTrackCountingHighEffBJetTags +
        sisC5PFTrackCountingHighPurBJetTags +
        sisC5PFJetProbabilityBJetTags +
        sisC5PFJetBProbabilityBJetTags
    )
)

sisC5PFJetBtaggingSV = cms.Sequence(
    sisC5PFImpactParameterTagInfos *
    sisC5PFSecondaryVertexTagInfos * (
        sisC5PFSimpleSecondaryVertexBJetTags +
        sisC5PFCombinedSecondaryVertexBJetTags +
        sisC5PFCombinedSecondaryVertexMVABJetTags
    )
)

sisC5PFJetBtaggingEle = cms.Sequence(
    btagSoftElectrons *
    sisC5PFSoftElectronTagInfos *
    sisC5PFSoftElectronBJetTags
)

sisC5PFJetBtaggingMu = cms.Sequence(
    sisC5PFSoftMuonTagInfos * (
        sisC5PFSoftMuonBJetTags +
        sisC5PFSoftMuonByIP3dBJetTags +
        sisC5PFSoftMuonByPtBJetTags
    )
)

sisC5PFJetBtagging = cms.Sequence(
    sisC5PFJetBtaggingIP +
    sisC5PFJetBtaggingSV +
    sisC5PFJetBtaggingEle +
    sisC5PFJetBtaggingMu
)

sisC5PFBtaggingSequence = cms.Sequence(
    sisC5PFJetTracksAssociator *
    sisC5PFJetBtagging
)

sisC5PFBtaggingPath = cms.Path(
    sisC5PFBtaggingSequence
)

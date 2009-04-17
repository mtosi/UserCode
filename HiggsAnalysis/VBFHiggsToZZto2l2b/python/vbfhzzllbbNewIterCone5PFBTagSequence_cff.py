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
iC5PFJetTracksAssociatorAtVertex = ic5JetTracksAssociatorAtVertex.clone()
iC5PFJetTracksAssociatorAtVertex.jets = "iterativeCone5PFJets"
iC5PFJetTracksAssociatorAtVertex.tracks = "generalTracks"


# The one needs to clone the b-tag producers and instruct them to use this new collection
from RecoBTag.Configuration.RecoBTag_cff import *

# impact parameter b-tag
iC5PFImpactParameterTagInfos = impactParameterTagInfos.clone()
iC5PFImpactParameterTagInfos.jetTracks = "iC5PFJetTracksAssociatorAtVertex"
iC5PFTrackCountingHighEffBJetTags = trackCountingHighEffBJetTags.clone()
iC5PFTrackCountingHighEffBJetTags.tagInfos = cms.VInputTag( cms.InputTag("iC5PFImpactParameterTagInfos") )
iC5PFTrackCountingHighPurBJetTags = trackCountingHighPurBJetTags.clone()
iC5PFTrackCountingHighPurBJetTags.tagInfos = cms.VInputTag( cms.InputTag("iC5PFImpactParameterTagInfos") )
iC5PFJetProbabilityBJetTags = jetProbabilityBJetTags.clone()
iC5PFJetProbabilityBJetTags.tagInfos = cms.VInputTag( cms.InputTag("iC5PFImpactParameterTagInfos") )
iC5PFJetBProbabilityBJetTags = jetBProbabilityBJetTags.clone()
iC5PFJetBProbabilityBJetTags.tagInfos = cms.VInputTag( cms.InputTag("iC5PFImpactParameterTagInfos") )

# secondary vertex b-tag
iC5PFSecondaryVertexTagInfos = secondaryVertexTagInfos.clone()
iC5PFSecondaryVertexTagInfos.trackIPTagInfos = "iC5PFImpactParameterTagInfos"
iC5PFSimpleSecondaryVertexBJetTags = simpleSecondaryVertexBJetTags.clone()
iC5PFSimpleSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("iC5PFSecondaryVertexTagInfos") )
iC5PFCombinedSecondaryVertexBJetTags = combinedSecondaryVertexBJetTags.clone()
iC5PFCombinedSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("iC5PFImpactParameterTagInfos"), cms.InputTag("iC5PFSecondaryVertexTagInfos") )
iC5PFCombinedSecondaryVertexMVABJetTags = combinedSecondaryVertexMVABJetTags.clone()
iC5PFCombinedSecondaryVertexMVABJetTags.tagInfos = cms.VInputTag( cms.InputTag("iC5PFImpactParameterTagInfos"), cms.InputTag("iC5PFSecondaryVertexTagInfos") )

# soft electron b-tag
iC5PFSoftElectronTagInfos = softElectronTagInfos.clone()
iC5PFSoftElectronTagInfos.jets = "iterativeCone5PFJets"
iC5PFSoftElectronBJetTags = softElectronBJetTags.clone()
iC5PFSoftElectronBJetTags.tagInfos = cms.VInputTag( cms.InputTag("iC5PFSoftElectronTagInfos") )

# soft muon b-tag
iC5PFSoftMuonTagInfos = softMuonTagInfos.clone()
iC5PFSoftMuonTagInfos.jets = "iterativeCone5PFJets"
iC5PFSoftMuonBJetTags = softMuonBJetTags.clone()
iC5PFSoftMuonBJetTags.tagInfos = cms.VInputTag( cms.InputTag("iC5PFSoftMuonTagInfos") )
iC5PFSoftMuonByIP3dBJetTags = softMuonByIP3dBJetTags.clone()
iC5PFSoftMuonByIP3dBJetTags.tagInfos = cms.VInputTag( cms.InputTag("iC5PFSoftMuonTagInfos") )
iC5PFSoftMuonByPtBJetTags = softMuonByPtBJetTags.clone()
iC5PFSoftMuonByPtBJetTags.tagInfos = cms.VInputTag( cms.InputTag("iC5PFSoftMuonTagInfos") )


# prepare a path running the new modules
iC5PFJetTracksAssociator = cms.Sequence(
    iC5PFJetTracksAssociatorAtVertex
)

iC5PFJetBtaggingIP = cms.Sequence(
    iC5PFImpactParameterTagInfos * (
        iC5PFTrackCountingHighEffBJetTags +
        iC5PFTrackCountingHighPurBJetTags +
        iC5PFJetProbabilityBJetTags +
        iC5PFJetBProbabilityBJetTags
    )
)

iC5PFJetBtaggingSV = cms.Sequence(
    iC5PFImpactParameterTagInfos *
    iC5PFSecondaryVertexTagInfos * (
        iC5PFSimpleSecondaryVertexBJetTags +
        iC5PFCombinedSecondaryVertexBJetTags +
        iC5PFCombinedSecondaryVertexMVABJetTags
    )
)

iC5PFJetBtaggingEle = cms.Sequence(
    btagSoftElectrons *
    iC5PFSoftElectronTagInfos *
    iC5PFSoftElectronBJetTags
)

iC5PFJetBtaggingMu = cms.Sequence(
    iC5PFSoftMuonTagInfos * (
        iC5PFSoftMuonBJetTags +
        iC5PFSoftMuonByIP3dBJetTags +
        iC5PFSoftMuonByPtBJetTags
    )
)

iC5PFJetBtagging = cms.Sequence(
    iC5PFJetBtaggingIP +
    iC5PFJetBtaggingSV +
    iC5PFJetBtaggingEle +
    iC5PFJetBtaggingMu
)

iC5PFBtaggingSequence = cms.Sequence(
    iC5PFJetTracksAssociator *
    iC5PFJetBtagging
)

iC5PFBtaggingPath = cms.Path(
    iC5PFBtaggingSequence
)

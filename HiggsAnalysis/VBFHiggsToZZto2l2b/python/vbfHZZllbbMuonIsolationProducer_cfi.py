import FWCore.ParameterSet.Config as cms

VBFHZZllbbMuonIsolationProducer = cms.EDProducer("VBFHZZllbbMuonIsolationProducer",
    HCALIsoDepositLabel    = cms.InputTag("muIsoDepositCalByAssociatorTowersNew","hcal"),
    HOCALIsoDepositLabel   = cms.InputTag("muIsoDepositCalByAssociatorTowersNew","ho"  ),
    ECALIsoDepositLabel    = cms.InputTag("muIsoDepositCalByAssociatorTowersNew","ecal"),
    TrackerIsoDepositLabel = cms.InputTag("muIsoDepositTkNew"),
    tracksLabel    = cms.InputTag("generalTracks"          ),
    muonsLabel     = cms.InputTag("VBFHZZllbbMuonSelector" ),
    electronsLabel = cms.InputTag("overlapElectronResolver"),
    isolationCut      = cms.double(60.0),
    isolationCone     = cms.double(0.3),
    isolationConeVeto = cms.double(0.015)
)



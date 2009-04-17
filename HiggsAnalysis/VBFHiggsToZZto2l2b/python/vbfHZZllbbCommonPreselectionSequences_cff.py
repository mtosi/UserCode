import FWCore.ParameterSet.Config as cms

vbfHZZmumubbCommonPreselection = cms.EDProducer("VBFHZZllbbCommonPreselection",
    HLabel               = cms.InputTag("vbfHZZmumubb"),
    ZllLabel             = cms.InputTag("zToMuMu"),
    decaychannel         = cms.string('2mu2b'),
    leptonLooseIsolLabel = cms.InputTag("vbfHZZllbbMuonIsolationProducer"),
    leptonLabel          = cms.InputTag("vbfHZZllbbMuonSelector"),
    cuts = cms.PSet(
        nLeptonCut             = cms.int32(2),
        nLooseLeptonCut        = cms.int32(2),
        llMassCut              = cms.double(12.0),
        fourBodyMassCut        = cms.double(100.0),
        numberOfLeptonCombsCut = cms.int32(1),
        numberOf2l2bCombsCut   = cms.int32(1)
    )
)

vbfHZZeebbCommonPreselection = cms.EDProducer("VBFHZZllbbCommonPreselection",
    HLabel               = cms.InputTag("vbfHZZeebb"),
    ZllLabel             = cms.InputTag("zToEE"),
    decaychannel         = cms.string('2e2b'),
    leptonLooseIsolLabel = cms.InputTag("vbfHZZllbbElectronIsolationProducer"),
    leptonLabel          = cms.InputTag("vbfHZZllbbElectronSelector"),
    cuts = cms.PSet(
        nLeptonCut             = cms.int32(2),
        nLooseLeptonCut        = cms.int32(2),
        llMassCut              = cms.double(12.0),
        fourBodyMassCut        = cms.double(100.0),
        numberOfLeptonCombsCut = cms.int32(1),
        numberOf2l2bCombsCut   = cms.int32(1)
    )
)

vbfHZZeebbCommonPreselectionFilter = cms.EDFilter("VBFHZZllbbCommonPreselectionFilter",
    decaychannel = cms.string('2e2b'),
    preSelInst   = cms.string('vbfHZZeebbCommonPreselection'),
    preSelLabels = cms.vstring('PreSelAtLeast2Electron', 
                               'PreSelAtLeast1ZEE', 
                               'PreSelAtLeast1H', 
                               'PreSelLoose2IsolElectron' 
                               ),
    preSelectFileName = cms.string('preSelec2e2b.out'),
    rootFileName      = cms.string('preSelec2e2b.root')
)

vbfHZZmumubbCommonPreselectionFilter = cms.EDFilter("VBFHZZllbbCommonPreselectionFilter",
    decaychannel = cms.string('2mu2b'),
    preSelInst   = cms.string('vbfHZZmmumubbCommonPreselection'),
    preSelLabels = cms.vstring('PreSelAtLeast2Muon', 
                               'PreSelAtLeast1ZMuMu', 
                               'PreSelAtLeast1H', 
                               'PreSelLoose2IsolMuon'
                               ),
    preSelectFileName = cms.string('preSelec2mu2b.out'),
    rootFileName = cms.string('preSelec2mu2b.root')
)

vbfHZZmumubbCommonPreselectionSequence = cms.Sequence(vbfHZZmumubbCommonPreselection
                                                      +vbfHZZmumubbCommonPreselectionFilter)

vbfHZZeebbCommonPreselectionSequence = cms.Sequence(vbfHZZeebbCommonPreselection
                                                    +vbfHZZeebbCommonPreselectionFilter)


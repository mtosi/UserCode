import FWCore.ParameterSet.Config as cms

vbfHZZmumubbCommonOfflineSelection = cms.EDProducer("VBFHZZllbbCommonOfflineSelection",
    decaychannel = cms.string('2mu2b'),
    bestCandidatesLeptons = cms.InputTag("vbfHZZllbbBestCandidateProducer:vbfHZZllbbBestCandidateLeptons"),
    useBestCandidate = cms.bool(True),                                                     
    # vertexing
    leptonLabelVertexLabel  = cms.InputTag("vbfHZZllbbMuonIsolationProducer"      ),
    leptonMapLabelVertLabel = cms.InputTag("vbfHZZllbbIpToVtxProducer:VertexMuMap"),
    vertexVarCut            = cms.vdouble(12.0, 8.0),
    # tight isolation
    leptonLabel    = cms.InputTag("vbfHZZllbbMuonIsolationProducer"      ),
    leptonMapLabel = cms.InputTag("vbfHZZllbbMuonIsolationProducerOffsel"),
    leptonIsoVarLabel = cms.InputTag("MuonIsolationX"),
    leptonIsoVarCut   = cms.vdouble(30.0)                                                       
)

vbfHZZeebbCommonOfflineSelection = cms.EDProducer("VBFHZZllbbCommonOfflineSelection",
    decaychannel = cms.string('2e2b'),
    bestCandidatesLeptons = cms.InputTag("vbfHZZllbbBestCandidateProducer:vbfHZZllbbBestCandidateLeptons"),
    useBestCandidate      = cms.bool(True),                                                     
    # vertexing
    leptonLabelVertLabel    = cms.InputTag("vbfHZZllbbHadIsolationProducer"),
    leptonMapLabelVertLabel = cms.InputTag("vbfHZZllbbIpToVtxProducer:VertexEleMap"),
    vertexVarCut            = cms.vdouble(12.0, 8.0),
    # tight isolation
    leptonLabel    = cms.InputTag("vbfHZZllbbElectronIsolationProducer"),
    leptonMapLabel = cms.InputTag("vbfHZZllbbHadIsolationProducer"),
    leptonIsoVarLabel = cms.InputTag("HadIsolationXele"),
    leptonIsoVarCut   = cms.vdouble(0.35)                                                  
)

vbfHZZmumubbCommonOfflineSelectionFilter = cms.EDFilter("VBFHZZllbbCommonOfflineSelectionFilter",
    decaychannel = cms.string('2mu2b'),
    offSelLabels = cms.vstring('OffSelTightCombIsolMuon', 
                               'OffSelVertexComb'),
    offSelectFileName = cms.string('offSelect2mu2b.out'),
    offSelInst   = cms.string('vbfHZZmumubbCommonOfflineSelection'),
    rootFileName = cms.string('offSelect2mu2b.root')
)

vbfHZZeebbCommonOfflineSelectionFilter = cms.EDFilter("VBFHZZllbbCommonOfflineSelectionFilter",
    decaychannel = cms.string('2e2b'),
    offSelLabels = cms.vstring('OffSelTightCombIsolElectron', 
                               'OffSelVertexComb'),
    offSelectFileName = cms.string('offSelect2e2b.out'),
    offSelInst = cms.string('vbfHZZeebbCommonOfflineSelection'),
    rootFileName = cms.string('offSelect2e2b.root')
)

vbfHZZmumubbCommonOffSelSequence = cms.Sequence(vbfHZZmumubbCommonOfflineSelection
                                                +vbfHZZmumubbCommonOfflineSelectionFilter)

vbfHZZeebbCommonOffSelSequence = cms.Sequence(vbfHZZeebbCommonOfflineSelection
                                              +vbfHZZeebbCommonOfflineSelectionFilter)

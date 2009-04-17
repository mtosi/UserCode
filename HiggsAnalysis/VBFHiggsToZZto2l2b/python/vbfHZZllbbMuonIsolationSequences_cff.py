import FWCore.ParameterSet.Config as cms

from RecoMuon.MuonIsolationProducers.muIsolation_cff import *

from RecoMuon.MuonIsolationProducers.muIsoDepositTk_cfi import *
muIsoDepositTkNew=RecoMuon.MuonIsolationProducers.muIsoDepositTk_cfi.muIsoDepositTk.clone()
muIsoDepositTkNew.IOPSet.inputMuonCollection = cms.InputTag("VBFHZZllbbMuonSelector")

from RecoMuon.MuonIsolationProducers.muIsoDepositCalByAssociatorTowers_cfi  import *
muIsoDepositCalByAssociatorTowersNew=RecoMuon.MuonIsolationProducers.muIsoDepositCalByAssociatorTowers_cfi.muIsoDepositCalByAssociatorTowers.clone()
muIsoDepositCalByAssociatorTowersNew.IOPSet.inputMuonCollection = cms.InputTag("VBFHZZllbbMuonSelector")

from RecoMuon.MuonIsolationProducers.muIsoDepositJets_cfi  import *
muIsoDepositJetsNew=RecoMuon.MuonIsolationProducers.muIsoDepositJets_cfi.muIsoDepositJets.clone()
muIsoDepositJetsNew.IOPSet.inputMuonCollection = cms.InputTag("VBFHZZllbbMuonSelector")

muIsoDeposits_muonsNew = cms.Sequence(muIsoDepositTkNew+muIsoDepositCalByAssociatorTowersNew+muIsoDepositJetsNew)
muIsolation_muonsNew = cms.Sequence(muIsoDeposits_muonsNew)
muIsolationNew = cms.Sequence(muIsolation_muonsNew)


from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfHZZllbbMuonIsolationProducer_cfi import *

VBFHZZllbbMuonIsolationSequence = cms.Sequence( muIsolationNew 
						+ VBFHZZllbbMuonIsolationProducer )

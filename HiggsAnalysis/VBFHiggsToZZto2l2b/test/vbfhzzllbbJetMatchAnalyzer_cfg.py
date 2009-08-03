import FWCore.ParameterSet.Config as cms

process = cms.Process("testJET")

process.load('Configuration/StandardSequences/Services_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.default.limit = 10
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.categories.append('testJETSummary')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    default        = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
    testJETSummary = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True) # default values if False
)


## debugging porpose
cms.Service('Tracer')
## if you get bad_alloc problems or you just observe increasing memory needs
## [for the real big problems use valgrind (http://valgrind.org) or igtools]
cms.Service('SimpleMemoryCheck', 
    ignoreTotal = cms.untracked.int32(1)
) 
## quickly identification on how fast your process and single modules are 
cms.Service('Timing')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

## this defines the input files
from HiggsAnalysis.VBFHiggsToZZto2l2b.Data.H150_ZZ_qqllSummer08_IDEALV9v2_GENSIMRECO_Input_cfi import * 
#from HiggsAnalysis.VBFHiggsToZZto2l2b.Data.PYTHIA6_SM_H_ZZ_qqll_mH150_10TeV_RECO_IDEAL_legnaro_cfi import *

# this inputs the input files from the previous function
process.source = RecoInput()

process.load("HiggsAnalysis.VBFHiggsToZZto2l2b.vbfHZZllbbIC5CaloCorrections_cff")

process.load("HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbMCprocessFilter_cfi")

process.load("PhysicsTools.RecoAlgos.TrackWithVertexRefSelector_cfi")
process.trackWithVertexRefSelector.src = cms.InputTag("generalTracks")

process.muonFilter = cms.EDFilter("MuonCountFilter",
    src = cms.InputTag("muons"),
#    ptMin = cms.double(10.0),
    minNumber = cms.uint32(2)
)

process.load("HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbMuonSelector_cfi")

process.load("HiggsAnalysis.VBFHiggsToZZto2l2b.zToMuMu_cfi")
#process.zToMuMu.decay = cms.string('muons@+ muons@-')
process.zToMuMu.cut = cms.string('mass > 20.0')

process.threeJetsFilter = cms.EDFilter("EtMinCaloJetCountFilter",
    src = cms.InputTag("iterativeCone5CaloJets"),
    etMin = cms.double(15.0),
    minNumber = cms.uint32(3)
)

process.twoForwardJetsFilter = cms.EDFilter("EtaEtMinCaloJetCountFilter",
    src = cms.InputTag("iterativeCone5CaloJets"),
    etMin = cms.double(10.0),
    etaMin = cms.double(-2.),
    etaMax = cms.double(2.),
    minNumber = cms.uint32(2)
)                                            

process.twoBJetsFilter = cms.EDFilter("JetTagCountFilter",
    minDiscriminator = cms.double(2.0),
    src = cms.InputTag("trackCountingHighEffBJetTags"),
    maxJetEta = cms.double(2.5),
    minNumber = cms.uint32(1),
    minJetEt = cms.double(15.0)
)

process.caloJetCollectionClone = cms.EDProducer("CaloJetShallowCloneProducer",
#    src = cms.InputTag("iterativeCone5CaloJets")
    src = cms.InputTag("L2L3CorJetIC5Calo")                                                
)

process.genJetCollectionClone = cms.EDProducer("GenJetShallowCloneProducer",
    src = cms.InputTag("iterativeCone5GenJets")
)

process.caloJetSele = cms.EDFilter("PtMinCandSelector",
    src = cms.InputTag("caloJetCollectionClone"),
    ptMin = cms.double(15.0)
)

process.genJetSele = cms.EDFilter("PtMinCandSelector",
    src = cms.InputTag("genJetCollectionClone"),
    ptMin = cms.double(15.0)
)

process.jetMatchOne = cms.EDFilter("CandOneToOneDeltaRMatcher",
    src = cms.InputTag("iterativeCone5GenJets"),
    algoMethod = cms.string('SwitchMode'),
#    matched = cms.InputTag("iterativeCone5CaloJets")
    matched = cms.InputTag("L2L3CorJetIC5Calo")
)

process.jetMatchMany = cms.EDFilter("CandOneToManyDeltaRMatcher",
    printDebug = cms.untracked.bool(True),
    src = cms.InputTag("genJetSele"),
    matched = cms.InputTag("caloJetSele")
)

process.printJet = cms.EDFilter("VBFHZZllbbJetMatchAnalyzer",
    source       = cms.InputTag("genJetSele"),
    matchMapMany = cms.InputTag("jetMatchMany"),
    matchMapOne  = cms.InputTag("jetMatchOne","src2mtc"),
    matched      = cms.InputTag("caloJetSele"),
    histomakerFlag = cms.bool(True),                            
    dRcut = cms.double(0.5)
)

process.effMatchedJet03 = cms.EDFilter("VBFHZZllbbJetMatchAnalyzer",
    source       = cms.InputTag("genJetSele"),
    matchMapMany = cms.InputTag("jetMatchMany"),
    matchMapOne  = cms.InputTag("jetMatchOne","src2mtc"),
    matched      = cms.InputTag("caloJetSele"),
    histomakerFlag = cms.bool(True),                            
    dRcut = cms.double(0.3)
)
process.effMatchedJet04 = cms.EDFilter("VBFHZZllbbJetMatchAnalyzer",
    source       = cms.InputTag("genJetSele"),
    matchMapMany = cms.InputTag("jetMatchMany"),
    matchMapOne  = cms.InputTag("jetMatchOne","src2mtc"),
    matched      = cms.InputTag("caloJetSele"),
    histomakerFlag = cms.bool(True),                            
    dRcut = cms.double(0.4)
)
process.effMatchedJet06 = cms.EDFilter("VBFHZZllbbJetMatchAnalyzer",
    source       = cms.InputTag("genJetSele"),
    matchMapMany = cms.InputTag("jetMatchMany"),
    matchMapOne  = cms.InputTag("jetMatchOne","src2mtc"),
    matched      = cms.InputTag("caloJetSele"),
    histomakerFlag = cms.bool(True),                            
    dRcut = cms.double(0.6)
)

## talk to TFileService for output histograms
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('myPlots.root')
)

## define output event selection to be that which satisfies 'p'
process.vbfHZZllbbEventSelection = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p')
    )
)

process.printEventNumber = cms.OutputModule("AsciiOutputModule")

from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfHZZllbbEventContent_cff import *
process.out = cms.OutputModule("PoolOutputModule",
    VBFHZZ2l2bEventContent,
    process.vbfHZZllbbEventSelection,
    fileName = cms.untracked.string('VBFHZZllbb.root'),
    verbose  = cms.untracked.bool(False)
)                               

# extend event content to include objects from EDNtuple
process.out.outputCommands.extend(["keep *_*_*_testJET"])
process.out.outputCommands.extend(["keep *_TrackWithVertexRefSelector_*_"])


process.p = cms.Path(
    process.vbfhzzllbbMCprocessFilter  *
    process.ic5CaloJetMETCorrections  *
    process.vbfhzzllbbMuonSelector     *
    process.muonFilter                 *
    process.zToMuMu                    *
#    process.dimuons                    *
#    process.zToMuMuGolden              *
#    ~process.twoForwardJetsFilter       *
    process.threeJetsFilter            *
    process.twoBJetsFilter             *
#    process.trackWithVertexRefSelector *
    process.caloJetCollectionClone    *
    process.genJetCollectionClone     *
    process.caloJetSele               *
    process.genJetSele                *
    process.jetMatchOne               *
    process.jetMatchMany              *
    process.printJet                  *
    process.effMatchedJet03           *
    process.effMatchedJet04           *
    process.effMatchedJet06           
)

process.outpath = cms.EndPath(
    process.out
#    +process.printEventNumber
)




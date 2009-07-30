from FWCore.ParameterSet.Config import *

process = Process("CandSelectorTest")

process.include("FWCore/MessageLogger/data/MessageLogger.cfi")

process.maxEvents = untracked.PSet( input = untracked.int32(10) )

## this defines the input files
from HiggsAnalysis.VBFHiggsToZZto2l2b.Data.H150_ZZ_qqllSummer08_IDEALV9v2_GENSIMRECO_Input_cfi import * 
#from HiggsAnalysis.VBFHiggsToZZto2l2b.Data.PYTHIA6_SM_H_ZZ_qqll_mH150_10TeV_RECO_IDEAL_legnaro_cfi import *

# this inputs the input files from the previous function
process.source = RecoInput()

process.goodMuonMCMatch = cms.EDFilter("MCMatcher",
  src = cms.InputTag("muons"),
  maxDPtRel = cms.double(1.0),
  mcPdgId = cms.vint32(13), ## muons

  mcStatus = cms.vint32(1),
  resolveByMatchQuality = cms.bool(False),
  maxDeltaR = cms.double(0.15),
  checkCharge = cms.bool(True),
  resolveAmbiguities = cms.bool(False),
  matched = cms.InputTag("genParticles")
)

# select the 10 particles with the larget Pt
process.largestPtCands = EDProducer("LargestPtCandSelector",
  src = InputTag('genParticles'),
  maxNumber = uint32( 10 )
)

## select only electrons, and save a vector of references 
#process.electronRefs = EDProducer("PdgIdCandRefVectorSelector",
#  src = InputTag("genParticles"),
#  pdgId = vint32( 11 )
#)

# select only electrons, and save clones
process.electrons = EDProducer("PdgIdCandSelector",
  src = InputTag("genParticles"),
  pdgId = vint32( 11 )
)

## select only muons, and save a vector of references 
#process.muonRefs = EDProducer("PdgIdCandRefVectorSelector",
#  src = InputTag("genParticles"),
#  pdgId = vint32( 13 )
#)

# select only muons, and save clones
process.muons = EDProducer("PdgIdCandSelector",
  src = InputTag("genParticles"),
  pdgId = vint32( 13 )
)

# select only b-quark, and save clones
process.bquarks = EDProducer("PdgIdCandSelector",
  src = InputTag("genParticles"),
  pdgId = vint32( 5 )
)

# select only quark, and save clones
process.quarks = EDProducer("PdgIdCandSelector",
  src = InputTag("genParticles"),
  pdgId = vint32( 1,2,3,4,5,6 )
)

# select only gluons, and save clones
process.quarks = EDProducer("PdgIdCandSelector",
  src = InputTag("genParticle"),
  pdgId = vint32( 21 )
)

# select only Z, and save clones
process.quarks = EDProducer("PdgIdCandSelector",
  src = InputTag("genParticles"),
  pdgId = vint32( 23 )
)

## select only electrons within eta and Pt cuts 
#process.bestElectrons = EDProducer("EtaPtMinCandViewSelector",
#  src = InputTag("electronRefs"),
#  ptMin = double( 20 ),
#  etaMin = double( -2.5 ),
#  etaMax = double( 2.5 )
#)

# make Z->e+e-
process.zTOeeCands = EDProducer("CandViewShallowCloneCombiner",
  decay = string("electrons@+ electrons@-"),
  cut = string("20 < mass < 200")
)

# make Z->mu+mu-
process.zTOmumuCands = EDProducer("CandViewShallowCloneCombiner",
  decay = string("muons@+ muons@-"),
  cut = string("20 < mass < 200")
)

# make Z->bbbar
process.zTObbCands = EDProducer("CandViewShallowCloneCombiner",
  decay = string("bquark@+ bquarks@-"),
  cut = string("20 < mass < 200")
)

# make exotic decay to three electron
process.HCands = EDProducer("CandViewShallowCloneCombiner",
  decay = string("zTOmumuCands@+ zTObbCands@-"),
  cut = string("20 < mass < 400")
)

# merge muons and electrons into leptons
process.leptons = EDProducer("CandMerger",
  src = VInputTag("electrons", "muons")
)

process.out = OutputModule("PoolOutputModule",
  fileName = untracked.string("cands.root")
)

process.printEventNumber = OutputModule("AsciiOutputModule")

process.select = Path(
   process.goodMuonMCMatch

#  process.goodMuons *
#  process.largestPtCands *
#  process.electrons *
#  process.muons *
#  process.leptons *
#  process.zTOeeCands *
#  process.zTOmumuCands *
#  process.zTObbCands *
#  process.HCands
)
 	
process.ep = EndPath(
  process.printEventNumber *
  process.out
)

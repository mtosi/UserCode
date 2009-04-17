import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

## this defines the input files
from PhysicsTools.PatAlgos.H150_ZZ_qqllSummer08_IDEALV9v2_GENSIMRECO_FULL_Input_cfi import *
#from PhysicsTools.PatAlgos.PYTHIA6_SM_H_ZZ_2l_2jets_mH800_10TeV_GEN_FASTSIM_Input_cfi import *
# this inputs the input files from the previous function
process.source = RecoInput()

process.demo = cms.EDAnalyzer('VBFHZZllbbMCvalidation',
                              OutputName  = cms.untracked.string('VBFHZZllbbValidationOutput.root'),
                              MCParticles = cms.untracked.string('genParticles'),
                              SIGNAL      = cms.untracked.bool(True)
                              )

process.p = cms.Path(process.demo)

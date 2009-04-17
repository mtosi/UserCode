# The following comments couldn't be translated into the new config version:

#! /bin/env cmsRun

import FWCore.ParameterSet.Config as cms

process = cms.Process("validation")
##keep the logging output to a nice level
process.load("FWCore.MessageLogger.MessageLogger_cfi")
## to write every 10th event
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    default = cms.untracked.PSet( limit = cms.untracked.int32(0)  )
)
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

## please note: this should not be needed; just a temporary hack
process.load("RecoBTag.Configuration.RecoBTag_cff")

## Analysis and validation plots
process.load("RecoBTag.Analysis.bTagAnalysis_cfi")
process.bTagCommonBlock.jetMCSrc = cms.InputTag("jetFlavour")
process.bTagAnalysis.rootfile = cms.string('testJetTagAnalysis.root')

process.maxEvents = cms.untracked.PSet(
#    input = cms.untracked.int32(-1)
    input = cms.untracked.int32(2)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)

process.plots = cms.Path(process.bTagAnalysis)
## this defines the input files
#from PhysicsTools.PatAlgos.PYTHIA6_SM_H_ZZ_2l_2jets_mH800_10TeV_GEN_FASTSIM_Input_cfi import *
from PhysicsTools.PatAlgos.H150_ZZ_qqllSummer08_IDEALV9v2_GENSIMRECO_Input_cfi import *
# this inputs the input files from the previous function
process.source = RecoInput()

process.bTagAnalysis.producePs = False

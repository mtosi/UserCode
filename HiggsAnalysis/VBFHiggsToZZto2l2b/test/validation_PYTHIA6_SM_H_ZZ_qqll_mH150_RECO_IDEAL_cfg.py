# The following comments couldn't be translated into the new config version:

#! /bin/env cmsRun

import FWCore.ParameterSet.Config as cms

process = cms.Process("validation")
process.load("DQMServices.Components.DQMEnvironment_cfi")

#keep the logging output to a nice level
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load("DQMServices.Core.DQM_cfg")

process.load("RecoBTag.Configuration.RecoBTag_cff")

process.load("Validation.RecoB.bTagAnalysis_cfi")


process.load("PhysicsTools.JetMCAlgos.CaloJetsMCFlavour_cfi")  
process.bTagValidation.jetMCSrc = 'IC5byValAlgo'
process.bTagValidation.allHistograms = True 
process.bTagValidation.fastMC = True
process.bTagValidation.ptRecJetMin = 15.0
process.bTagCommonBlock.ptRecJetMin = 15.0


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)

process.plots = cms.Path(process.myPartons* process.iterativeCone5Flavour * process.bTagValidation*process.dqmSaver)
process.dqmEnv.subSystemFolder = 'BTAG'
process.dqmSaver.producer = 'DQM'
process.dqmSaver.workflow = '/POG/BTAG/BJET'
process.dqmSaver.convention = 'RelVal'

## this defines the input files
from HiggsAnalysis.VBFHiggsToZZto2l2b.Data.PYTHIA6_SM_H_ZZ_qqll_mH150_10TeV_RECO_IDEAL_legnaro_cfi import *
# this inputs the input files from the previous function
process.source = RecoInput()


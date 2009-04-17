#! /usr/bin/env cmsRun
import FWCore.ParameterSet.Config as cms

process = cms.Process("validation")

# keep the logging output to a nice level
process.load("FWCore.MessageLogger.MessageLogger_cfi")


# please note: this should not be needed; just a temporary hack
#process.load("RecoBTag.Configuration.RecoBTag_cff")

# default configuration with frontier conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.GlobalTag.globaltag = cms.string('IDEAL_V9::All')
process.load("PhysicsTools.JetMCAlgos.CaloJetsMCFlavour_cfi")

## this defines the input files
from PhysicsTools.PatAlgos.H150_ZZ_qqllSummer08_IDEALV9v2_GENSIMRECO_FULL_Input_cfi import *

# this inputs the input files from the previous function
process.source = RecoInput()


## set the number of events
process.maxEvents = cms.untracked.PSet(
#    input = cms.untracked.int32(-1)
    input = cms.untracked.int32(20)
)


# validation plots
process.load("RecoBTag.Analysis.bTagAnalysis_cfi")
#process.bTagAnalysis.fastMC = True
process.bTagAnalysis.fastMC = False
process.bTagAnalysis.jetMCSrc = 'IC5byValAlgo'
#process.bTagAnalysis.jetMCSrc = 'IC5byValPhys'
#process.bTagAnalysis.jetMCSrc = 'IC5byRef'
#process.bTagAnalysis.allHistograms = False
process.bTagAnalysis.allHistograms = True
process.bTagAnalysis.rootfile = cms.string('IC5byValAlgojetTagAnalysis.root')
#process.bTagAnalysis.rootfile = cms.string('IC5byValPhysjetTagAnalysis.root')
#process.bTagAnalysis.rootfile = cms.string('IC5byRefjetTagAnalysis.root')

# run btag, then validation
process.btag  = cms.Path(process.btagging)
#process.plots = cms.Path(process.myPartons + process.iterativeCone5Flavour + process.bTagAnalysis)
process.plots = cms.Path(process.caloJetMCFlavour + process.bTagAnalysis)

process.schedule = cms.Schedule(process.btag,
                                process.plots
)

# define output path
#process.outpath = cms.EndPath(process.bTagAnalysis)

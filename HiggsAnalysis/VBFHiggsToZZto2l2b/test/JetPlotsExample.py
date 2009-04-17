# PYTHON configuration file for class: JetPlotsExample
# Description:  Example of simple EDAnalyzer for jets.
# Author: K. Kousouris
# Date:  25 - August - 2008
import FWCore.ParameterSet.Config as cms

process = cms.Process("Ana")
process.load("FWCore.MessageService.MessageLogger_cfi")
#############   Set the number of events #############
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
#############   Define the source file ###############
## this defines the input files
from PhysicsTools.PatAlgos.H150_ZZ_qqllSummer08_IDEALV9v2_GENSIMRECO_FULL_Input_cfi import *

# this inputs the input files from the previous function
process.source = RecoInput()

#############   Calo Jets  ###########################
process.calo = cms.EDAnalyzer("CaloJetPlotsExample",
    JetAlgorithm  = cms.string('iterativeCone5CaloJets'),
    HistoFileName = cms.string('CaloJetPlotsExample.root'),
    NJets         = cms.int32(2)
)
#############   Gen Jets   ###########################
process.gen = cms.EDAnalyzer("GenJetPlotsExample",
    JetAlgorithm  = cms.string('iterativeCone5GenJets'),
    HistoFileName = cms.string('GenJetPlotsExample.root'),
    NJets         = cms.int32(2)
)
#############   PF Jets    ###########################
process.pf = cms.EDAnalyzer("PFJetPlotsExample",
    JetAlgorithm  = cms.string('iterativeCone5PFJets'),
    HistoFileName = cms.string('PFJetPlotsExample.root'),
    NJets         = cms.int32(2)
)
#############   Path       ###########################
process.p = cms.Path(process.calo*process.gen*process.pf)
#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 10


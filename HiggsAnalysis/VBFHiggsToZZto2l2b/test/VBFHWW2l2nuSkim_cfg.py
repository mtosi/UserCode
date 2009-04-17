#!/usr/bin/env cmsRun
import FWCore.ParameterSet.Config as cms

process = cms.Process("VBFHWW2l2nuSkim")



process.load('Configuration/StandardSequences/Services_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'

##########
# Source #
##########


process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source(
    "PoolSource",
    debugFlag = cms.untracked.bool(True),
    debugVebosity = cms.untracked.uint32(1),
    fileNames = cms.untracked.vstring(
        'file:/tmp/amassiro/VBFHWW2l2nuTest_8.root'
#         'rfio:/castor/cern.ch/cms/store/path/filename.root'
        )
    )

   
process.muonFilter = cms.EDFilter("PtMinMuonCountFilter",
  src = cms.InputTag("muons"),
  minNumber = cms.uint32(1),
  ptMin = cms.double(5.0)
)

process.electronFilter = cms.EDFilter("PtMinPixelMatchGsfElectronCountFilter",
  src = cms.InputTag("pixelMatchGsfElectrons"),
  minNumber = cms.uint32(1),
  ptMin = cms.double(5.0) 
)



########
# Path #
########

process.countMu= cms.Path (process.muonFilter)
process.countEle = cms.Path (process.electronFilter)


##########
# Output #
##########

from HiggsAnalysis.VBFHiggsToWWto2l2nu.VBFHWW2l2nuEventContent_cff import *

process.VBFHWW2l2nuOutputModule = cms.OutputModule(
    "PoolOutputModule",
    VBFHWW2l2nuEventContent,
    dataset = cms.untracked.PSet(dataTier = cms.untracked.string('USER')),
    fileName = cms.untracked.string('VBFHWW2l2nuTest.root'),
#     fileName = cms.untracked.string('/tmp/amassiro/VBFHWW2l2nuTest_12Feb09_nuovoFiltro.root'),
   
    SelectEvents = cms.untracked.PSet(
                SelectEvents = cms.vstring('countMu','countEle')
    )
)


process.o = cms.EndPath ( process.VBFHWW2l2nuOutputModule )

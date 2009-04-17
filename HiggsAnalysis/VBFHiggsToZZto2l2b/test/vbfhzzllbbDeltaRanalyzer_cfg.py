import FWCore.ParameterSet.Config as cms

process = cms.Process("VBFHZZllbbAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
## to write every 10th event
process.MessageLogger.cerr.FwkReport.reportEvery = 10
#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.cerr.threshold = 'ERROR'
process.MessageLogger.categories.append('VBFHZZllbbAnalysisSummary')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    default                   = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
    VBFHZZllbbAnalysisSummary = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True) # default values if False
)

## this defines the input files
from PhysicsTools.PatAlgos.H150_ZZ_qqllSummer08_IDEALV9v2_GENSIMRECO_FULL_Input_cfi import *

# this inputs the input files from the previous function
process.source = RecoInput()

## set the number of events
process.maxEvents = cms.untracked.PSet(
#    input = cms.untracked.int32(-1)
    input = cms.untracked.int32(20)
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

## talk to TFileService for output histograms
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('VBFHZZllbbDeltaRAnalysisHistos.root')
)

process.load("HiggsAnalysis.VBFHiggsToZZto2l2b.vbfHZZllbbIC5CaloCorrections_cff")

#from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbDeltaRanalyzer_cfi import *
process.vbfHZZllbbDeltaRAnalyzer = cms.EDAnalyzer('VBFHZZllbbDeltaRAnalyzer',
    genParticleLabel = cms.untracked.InputTag( 'genParticles'                    ),
    caloJetLabel     = cms.untracked.InputTag( 'L2L3CorJetIC5Calo'               ),
#    caloJetLabel     = cms.untracked.InputTag( 'iterativeCone5CaloJets'          ),
    muonLabel        = cms.untracked.InputTag( 'muons'                           ),
    electronLabel    = cms.untracked.InputTag( 'pixelMatchGsfElectrons'          ),
    bTagLabel        = cms.untracked.InputTag( 'combinedSecondaryVertexBJetTags' ),
    metLabel         = cms.untracked.InputTag( 'corMetType1Icone5'               ),
#    metLabel         = cms.untracked.InputTag( 'met'                             ),
    ptMax = cms.untracked.double(10)                           
)

process.corrections = cms.Sequence( process.L2L3CorJetIC5Calo *
                                    process.corMetType1Icone5
)


process.p = cms.Path(
    process.corrections *
    process.vbfHZZllbbDeltaRAnalyzer
)

## define output event selection to be that which satisfies 'p'
process.vbfHZZllbbEventSelection = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p')
    )
)

## output module configuration
process.out = cms.OutputModule("PoolOutputModule",
    process.vbfHZZllbbEventSelection,
    fileName = cms.untracked.string('VBFHZZllbbDeltaRAnalysisSkim.root'),
    verbose  = cms.untracked.bool(False),
#    outputCommands = cms.untracked.vstring('keep *')
    outputCommands = cms.untracked.vstring('drop *')
)

# extend event content to include objects from EDNtuple
# MC info
process.out.outputCommands.extend(["keep recoGenParticles_genParticles_*_*"       ])
process.out.outputCommands.extend(["keep *_source_*_*"                            ])
# muons                                                                           
process.out.outputCommands.extend(["keep *_muons_*_*"                             ])
# electrons!!!!!!!!
process.out.outputCommands.extend(["keep *_pixelMatchGsfElectrons_*_*"            ])
# jets                                                                            
process.out.outputCommands.extend(["keep *_iterativeCone5GenJets_*_*"             ])
process.out.outputCommands.extend(["keep *_iterativeCone5CaloJets_*_*"            ])
process.out.outputCommands.extend(["keep *_kt4GenJets_*_*"                        ])
process.out.outputCommands.extend(["keep *_kt4CaloJets_*_*"                       ])
process.out.outputCommands.extend(["keep *_kt6GenJets_*_*"                        ])
process.out.outputCommands.extend(["keep *_kt6CaloJets_*_*"                       ])
## # PF                                                                              
#process.out.outputCommands.extend(["keep *recoPF*_*_*_*"                      ])
## process.out.outputCommands.extend(["keep *_particleFlow_*_*"                      ])
process.out.outputCommands.extend(["keep *_iterativeCone5PFJets_*_*"              ]) 
process.out.outputCommands.extend(["keep *_kt4PFJets_*_*"                         ])
process.out.outputCommands.extend(["keep *_kt6PFJets_*_*"                         ])
# met                                                                             
process.out.outputCommands.extend(["keep *_genMet_*_*"                            ])
process.out.outputCommands.extend(["keep *_met_*_*"                               ])
process.out.outputCommands.extend(["keep *_corMetType1Icone5_*_*"                 ])
#process.out.outputCommands.extend(["keep *_htMet*_*_*"                            ])
# b-tagging
process.out.outputCommands.extend(["keep *_simpleSecondaryVertexBJetTags_*_*"     ])
process.out.outputCommands.extend(["keep *_combinedSecondaryVertexBJetTags_*_*"   ])
process.out.outputCommands.extend(["keep *_combinedSecondaryVertexMVABJetTags_*_*"])
process.out.outputCommands.extend(["keep *_impactParameterMVABJetTags_*_*"        ])
process.out.outputCommands.extend(["keep *_jetProbabilityBJetTags_*_*"            ])
process.out.outputCommands.extend(["keep *_jetBProbabilityBJetTags_*_*"           ])
process.out.outputCommands.extend(["keep *_trackCountingHighEffBJetTags_*_*"      ])
process.out.outputCommands.extend(["keep *_trackCountingHighPurBJetTags_*_*"      ])
process.out.outputCommands.extend(["keep *_impactParameterTagInfos_*_*"           ])
process.out.outputCommands.extend(["keep *_secondaryVertexTagInfos_*_*"           ]) 
# primary vertices
process.out.outputCommands.extend(["keep *_offlinePrimaryVertices_*_*"            ])
process.out.outputCommands.extend(["keep *_offlinePrimaryVerticesWithBS_*_*"      ])
process.out.outputCommands.extend(["keep *_pixelVertices_*_*"      ])


process.out.outputCommands.extend(["keep *_*_*_VBFHZZllbbAnalysis"      ])

# define output path
process.outpath = cms.EndPath(process.out)


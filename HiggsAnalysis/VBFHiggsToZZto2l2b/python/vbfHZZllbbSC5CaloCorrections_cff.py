import FWCore.ParameterSet.Config as cms

# ## load geometry
from Configuration.StandardSequences.Geometry_cff import *
from Configuration.StandardSequences.FrontierConditions_GlobalTag_cff import *
from Geometry.CommonDetUnit.globalTrackingGeometry_cfi import *
#process.GlobalTag.globaltag = cms.string('IDEAL_V11::All')
GlobalTag.globaltag = cms.string('IDEAL_V9::All')
from Configuration.StandardSequences.MagneticField_cff import *

from TrackingTools.TrackAssociator.default_cfi import *
from TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff import *

# default JetMET calibration on IC, KT and MC Jets

# define eta and pt correction services
from JetMETCorrections.Configuration.L2L3Corrections_Summer08Redigi_cff import *
# define jet flavour correction services
#process.load("JetMETCorrections.Configuration.L5FlavorCorrections_cff")
# define jet parton correction services
from JetMETCorrections.Configuration.L7PartonCorrections_cff import *

# pick theL2+L3 corrections
es_prefer_L2L3JetCorrectorIC5 = cms.ESPrefer("JetCorrectionServiceChain",
                                             "L2L3JetCorrectorSC5Calo")
# MET corrections from JES
from JetMETCorrections.Type1MET.MetType1Corrections_cff import *
# change corrector to L2+L3
corMetType1Scone5 = corMetType1Icone5.clone()
corMetType1Scone5.inputUncorJetsLabel = cms.string('sisCone5CaloJets')
corMetType1Scone5.corrector = cms.string('L2L3JetCorrectorSC5Calo')

# MET corrections from muons
from JetMETCorrections.Type1MET.MetMuonCorrections_cff import corMetGlobalMuons,goodMuonsforMETCorrection
# muon MET correction maker 
globalMuonsForMET     = goodMuonsforMETCorrection.clone(
   cut = cms.string('isGlobalMuon = 1')
)
goodGlobalMuonsForMET = goodMuonsforMETCorrection.clone(
    src = cms.InputTag("globalMuonsForMET")
)
corMetType1Scone5Muons = corMetGlobalMuons.clone(
    inputUncorMetLabel = cms.InputTag('corMetType1Scone5'),
    muonsInputTag      = cms.InputTag('goodGlobalMuonsForMET')
)

# It would be better to get this config to JetMETCorrections/Type1MET/data/ at some point
corMetType1Scone5Muons.TrackAssociatorParameters.useEcal    = False ## RecoHits
corMetType1Scone5Muons.TrackAssociatorParameters.useHcal    = False ## RecoHits
corMetType1Scone5Muons.TrackAssociatorParameters.useHO      = False ## RecoHits
corMetType1Scone5Muons.TrackAssociatorParameters.useCalo    = True  ## CaloTowers
corMetType1Scone5Muons.TrackAssociatorParameters.useMuon    = False ## RecoHits
corMetType1Scone5Muons.TrackAssociatorParameters.truthMatch = False

# default sequence for JetMET corrections
sc5CaloJetMETCorrections = cms.Sequence(L2L3CorJetSC5Calo      +
                                        globalMuonsForMET      *
                                        goodGlobalMuonsForMET  *
                                        corMetType1Scone5      *
                                        corMetType1Scone5Muons
)

# default sequence for JetMET corrections
sc5CaloJetMETCorrections_withoutMuonCorr = cms.Sequence(L2L3CorJetSC5Calo +
                                                        corMetType1Scone5   

)

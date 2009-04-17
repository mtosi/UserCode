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
from JetMETCorrections.Configuration.L2L3Corrections_Summer08_cff import *
# define jet flavour correction services
#process.load("JetMETCorrections.Configuration.L5FlavorCorrections_cff")
# define jet parton correction services
from JetMETCorrections.Configuration.L7PartonCorrections_cff import *

# pick theL2+L3 corrections
es_prefer_L2L3JetCorrectorIC5 = cms.ESPrefer("JetCorrectionServiceChain",
                                             "L2L3JetCorrectorIC5PF")

# default sequence for JetMET corrections
ic5PFJetMETCorrections = cms.Sequence(L2L3CorJetIC5PF)

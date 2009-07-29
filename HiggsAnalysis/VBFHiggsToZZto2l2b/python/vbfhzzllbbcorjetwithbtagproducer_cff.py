import FWCore.ParameterSet.Config as cms

from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfHZZllbbIC5CaloCorrections_cff import *
from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfhzzllbbcorjetwithbtagproducer_cfi import *

vbfhzzllbbCorJetWithBTagSequence = cms.Sequence(
    ic5CaloJetMETCorrections     *
    vbfhzzllbbCorJetWithBTagProd
)

vbfhzzllbbCorJetWithBTag = cms.Path(
    vbfhzzllbbCorJetWithBTagSequence
)

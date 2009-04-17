import FWCore.ParameterSet.Config as cms

from HiggsAnalysis.VBFHiggsToZZto2l2b.overlapElectronResolver_cfi import *
from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfHZZllbbElectronIdSequences_cff import *
from HiggsAnalysis.VBFHiggsToZZto2l2b.vbfHZZllbbElectronSelector_cfi import *

vbfHZZllbbElectronSequence = cms.Sequence(overlapElectronResolver
                                          + vbfHZZllbbElectronIdSequence
                                          + vbfHZZllbbElectronSelector)


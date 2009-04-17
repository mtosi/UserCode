import FWCore.ParameterSet.Config as cms


from RecoEgamma.ElectronIdentification.electronIdCutBasedClassesExt_cfi import *
import RecoEgamma.ElectronIdentification.electronIdCutBasedClassesExt_cfi
eidClassLoose = RecoEgamma.ElectronIdentification.electronIdCutBasedClassesExt_cfi.eidCutBasedClassesExt.clone()
eidClassLoose.src = "overlapElectronResolver"
import RecoEgamma.ElectronIdentification.electronIdCutBasedClassesExt_cfi
eidClassMedium = RecoEgamma.ElectronIdentification.electronIdCutBasedClassesExt_cfi.eidCutBasedClassesExt.clone()
eidClassMedium.electronQuality = 'medium'
eidClassMedium.src = "overlapElectronResolver"
import RecoEgamma.ElectronIdentification.electronIdCutBasedClassesExt_cfi
eidClassTight = RecoEgamma.ElectronIdentification.electronIdCutBasedClassesExt_cfi.eidCutBasedClassesExt.clone()
eidClassTight.electronQuality = 'tight'
eidClassTight.src = "overlapElectronResolver"

vbfHZZllbbElectronIdSequence = cms.Sequence( eidClassLoose
                                             + eidClassMedium
                                             + eidClassTight )

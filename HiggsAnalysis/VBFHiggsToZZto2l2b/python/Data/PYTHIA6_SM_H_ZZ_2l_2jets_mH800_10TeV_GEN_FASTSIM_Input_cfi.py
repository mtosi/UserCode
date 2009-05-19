## import configurations
import FWCore.ParameterSet.Config as cms

# from 

def RecoInput() : 
 return cms.Source("PoolSource",
                   debugVerbosity = cms.untracked.uint32(0),
                   debugFlag = cms.untracked.bool(True),
                   fileNames = cms.untracked.vstring(
                        'file:/data/tosi/PYTHIA6_SM_H_ZZ_2l_2jets_10TeV_FASTSIM/PYTHIA6_SM_H_ZZ_2l_2jets_mH800_10TeV_cff_py_GEN_FASTSIM_1.root'
                    )
                   )

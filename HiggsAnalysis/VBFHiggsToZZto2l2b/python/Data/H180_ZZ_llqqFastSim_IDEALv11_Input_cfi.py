## import configurations
import FWCore.ParameterSet.Config as cms

# from 

def RecoInput() : 
 return cms.Source("PoolSource",
                   debugVerbosity = cms.untracked.uint32(0),
                   debugFlag = cms.untracked.bool(True),
                   fileNames = cms.untracked.vstring(
                        'castor:/castor/cern.ch/user/d/dorigo/HZZ/FastSim180/PYTHIA6_SM_H_ZZ_2l_2jets_mH180_10TeV_cff_py_GEN_FASTSIM_1.root',
                        'castor:/castor/cern.ch/user/d/dorigo/HZZ/FastSim180/PYTHIA6_SM_H_ZZ_2l_2jets_mH180_10TeV_cff_py_GEN_FASTSIM_2.root',
                        'castor:/castor/cern.ch/user/d/dorigo/HZZ/FastSim180/PYTHIA6_SM_H_ZZ_2l_2jets_mH180_10TeV_cff_py_GEN_FASTSIM_3.root',
                        'castor:/castor/cern.ch/user/d/dorigo/HZZ/FastSim180/PYTHIA6_SM_H_ZZ_2l_2jets_mH180_10TeV_cff_py_GEN_FASTSIM_4.root',
                        'castor:/castor/cern.ch/user/d/dorigo/HZZ/FastSim180/PYTHIA6_SM_H_ZZ_2l_2jets_mH180_10TeV_cff_py_GEN_FASTSIM_5.root',
                        'castor:/castor/cern.ch/user/d/dorigo/HZZ/FastSim180/PYTHIA6_SM_H_ZZ_2l_2jets_mH180_10TeV_cff_py_GEN_FASTSIM_6.root',
                        'castor:/castor/cern.ch/user/d/dorigo/HZZ/FastSim180/PYTHIA6_SM_H_ZZ_2l_2jets_mH180_10TeV_cff_py_GEN_FASTSIM_7.root',
                        'castor:/castor/cern.ch/user/d/dorigo/HZZ/FastSim180/PYTHIA6_SM_H_ZZ_2l_2jets_mH180_10TeV_cff_py_GEN_FASTSIM_8.root',
                        'castor:/castor/cern.ch/user/d/dorigo/HZZ/FastSim180/PYTHIA6_SM_H_ZZ_2l_2jets_mH180_10TeV_cff_py_GEN_FASTSIM_9.root',
                        'castor:/castor/cern.ch/user/d/dorigo/HZZ/FastSim180/PYTHIA6_SM_H_ZZ_2l_2jets_mH180_10TeV_cff_py_GEN_FASTSIM_10.root',
                        'castor:/castor/cern.ch/user/d/dorigo/HZZ/FastSim180/PYTHIA6_SM_H_ZZ_2l_2jets_mH180_10TeV_cff_py_GEN_FASTSIM_10K_1.root',
                        'castor:/castor/cern.ch/user/d/dorigo/HZZ/FastSim180/PYTHIA6_SM_H_ZZ_2l_2jets_mH180_10TeV_cff_py_GEN_FASTSIM_10K_2.root',
                        'castor:/castor/cern.ch/user/d/dorigo/HZZ/FastSim180/PYTHIA6_SM_H_ZZ_2l_2jets_mH180_10TeV_cff_py_GEN_FASTSIM_10K_3.root',
                        'castor:/castor/cern.ch/user/d/dorigo/HZZ/FastSim180/PYTHIA6_SM_H_ZZ_2l_2jets_mH180_10TeV_cff_py_GEN_FASTSIM_10K_4.root',
                        'castor:/castor/cern.ch/user/d/dorigo/HZZ/FastSim180/PYTHIA6_SM_H_ZZ_2l_2jets_mH180_10TeV_cff_py_GEN_FASTSIM_10K_5.root',
                        'castor:/castor/cern.ch/user/d/dorigo/HZZ/FastSim180/PYTHIA6_SM_H_ZZ_2l_2jets_mH180_10TeV_cff_py_GEN_FASTSIM_10K_7.root',
                        'castor:/castor/cern.ch/user/d/dorigo/HZZ/FastSim180/PYTHIA6_SM_H_ZZ_2l_2jets_mH180_10TeV_cff_py_GEN_FASTSIM_10K_8.root',
                        'castor:/castor/cern.ch/user/d/dorigo/HZZ/FastSim180/PYTHIA6_SM_H_ZZ_2l_2jets_mH180_10TeV_cff_py_GEN_FASTSIM_10K_9.root',
                        'castor:/castor/cern.ch/user/d/dorigo/HZZ/FastSim180/PYTHIA6_SM_H_ZZ_2l_2jets_mH180_10TeV_cff_py_GEN_FASTSIM_10K_10.root'
                        )
                   )                   

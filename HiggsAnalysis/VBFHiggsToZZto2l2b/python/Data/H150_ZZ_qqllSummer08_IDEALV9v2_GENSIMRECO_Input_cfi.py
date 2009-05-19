## import configurations
import FWCore.ParameterSet.Config as cms

# from 

def RecoInput() : 
 return cms.Source("PoolSource",
                   debugVerbosity = cms.untracked.uint32(0),
                   debugFlag = cms.untracked.bool(True),
                   fileNames = cms.untracked.vstring(
                        'file:/data/tosi/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2Files/H150_ZZ_qqllSummer08_IDEALV9v2_GEN-SIM-RECO_1.root',
                        'file:/data/tosi/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2Files/H150_ZZ_qqllSummer08_IDEALV9v2_GEN-SIM-RECO_2.root',
                        'file:/data/tosi/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2Files/H150_ZZ_qqllSummer08_IDEALV9v2_GEN-SIM-RECO_3.root',
                        'file:/data/tosi/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2Files/H150_ZZ_qqllSummer08_IDEALV9v2_GEN-SIM-RECO_4.root',
                        'file:/data/tosi/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2Files/H150_ZZ_qqllSummer08_IDEALV9v2_GEN-SIM-RECO_5.root',
                        'file:/data/tosi/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2Files/H150_ZZ_qqllSummer08_IDEALV9v2_GEN-SIM-RECO_6.root',
                        'file:/data/tosi/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2Files/H150_ZZ_qqllSummer08_IDEALV9v2_GEN-SIM-RECO_7.root',
                        'file:/data/tosi/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2Files/H150_ZZ_qqllSummer08_IDEALV9v2_GEN-SIM-RECO_8.root',
                        'file:/data/tosi/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2Files/H150_ZZ_qqllSummer08_IDEALV9v2_GEN-SIM-RECO_9.root',
                        'file:/data/tosi/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2Files/H150_ZZ_qqllSummer08_IDEALV9v2_GEN-SIM-RECO_10.root'
                        )
                   )

## import configurations
import FWCore.ParameterSet.Config as cms

# from 

def RecoInput() : 
 return cms.Source("PoolSource",
                   debugVerbosity = cms.untracked.uint32(0),
                   debugFlag = cms.untracked.bool(True),
                   fileNames = cms.untracked.vstring(
                        'file:/data/tosi/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2Files/H150_ZZ_qqllSummer08_IDEALV9v2_GEN-SIM-RECO_FULL1.root',
                        'file:/data/tosi/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2Files/H150_ZZ_qqllSummer08_IDEALV9v2_GEN-SIM-RECO_FULL2.root',
                        'file:/data/tosi/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2Files/H150_ZZ_qqllSummer08_IDEALV9v2_GEN-SIM-RECO_FULL3.root',
                        'file:/data/tosi/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2Files/H150_ZZ_qqllSummer08_IDEALV9v2_GEN-SIM-RECO_FULL4.root',
                        'file:/data/tosi/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2Files/H150_ZZ_qqllSummer08_IDEALV9v2_GEN-SIM-RECO_FULL5.root',
                        'file:/data/tosi/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2Files/H150_ZZ_qqllSummer08_IDEALV9v2_GEN-SIM-RECO_FULL6.root',
                        'file:/data/tosi/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2Files/H150_ZZ_qqllSummer08_IDEALV9v2_GEN-SIM-RECO_FULL7.root',
                        'file:/data/tosi/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2Files/H150_ZZ_qqllSummer08_IDEALV9v2_GEN-SIM-RECO_FULL8.root',
                        'file:/data/tosi/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2Files/H150_ZZ_qqllSummer08_IDEALV9v2_GEN-SIM-RECO_FULL9.root',
                        'file:/data/tosi/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2/PYTHIA6_SM_H_ZZ_qqll_Summer08IDEALV9v2Files/H150_ZZ_qqllSummer08_IDEALV9v2_GEN-SIM-RECO_FULL10.root'
                        )
                   )

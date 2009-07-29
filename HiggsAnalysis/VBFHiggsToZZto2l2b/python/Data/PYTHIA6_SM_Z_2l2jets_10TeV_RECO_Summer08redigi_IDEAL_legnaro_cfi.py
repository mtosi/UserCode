## import configurations
import FWCore.ParameterSet.Config as cms

# from 

def RecoInput() : 
 return cms.Source("PoolSource",
                   debugVerbosity = cms.untracked.uint32(0),
                   debugFlag = cms.untracked.bool(True),
                   fileNames = cms.untracked.vstring(
                '/store/mc/Summer08/Z_2l2jets/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0017/C672F4A7-70F3-DD11-B661-00093D120FD9.root',
                '/store/mc/Summer08/Z_2l2jets/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0016/04C4C143-05F3-DD11-A720-001A4BA97406.root',
                '/store/mc/Summer08/Z_2l2jets/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0014/FC40AA20-DAF2-DD11-8FC0-001EC9D7FA3C.root',
                '/store/mc/Summer08/Z_2l2jets/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0014/ECA70B66-DCF2-DD11-ADAB-001EC9D26F7D.root',
                '/store/mc/Summer08/Z_2l2jets/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0014/B6D225E1-D9F2-DD11-AFD3-001EC9D7FF5B.root',
                '/store/mc/Summer08/Z_2l2jets/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0014/B0D1FC09-DCF2-DD11-B207-001EC9D83165.root',
                '/store/mc/Summer08/Z_2l2jets/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0014/966FFA84-DAF2-DD11-9B85-001EC9D8BDEF.root',
                '/store/mc/Summer08/Z_2l2jets/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0014/8A49960A-DCF2-DD11-80BD-001EC9D80AB5.root',
                '/store/mc/Summer08/Z_2l2jets/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0014/864AADA7-DAF2-DD11-9F33-001EC9D8D085.root',
                '/store/mc/Summer08/Z_2l2jets/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0014/6E99D0D5-D9F2-DD11-AD87-0030487C6A2C.root',
                '/store/mc/Summer08/Z_2l2jets/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0014/44C38767-DCF2-DD11-80E8-001EC9D80A95.root',
                '/store/mc/Summer08/Z_2l2jets/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0014/3E5B987B-DAF2-DD11-98FA-00093D00932C.root',
                '/store/mc/Summer08/Z_2l2jets/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0014/324F4049-DCF2-DD11-B9BB-00093D120700.root',
                '/store/mc/Summer08/Z_2l2jets/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0014/2A2C6419-DCF2-DD11-93FD-0019BB3FD4AA.root',
                '/store/mc/Summer08/Z_2l2jets/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0014/1E00CC7F-DCF2-DD11-8D2C-001E4F3F165E.root',
                '/store/mc/Summer08/Z_2l2jets/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0014/163FBBB0-DAF2-DD11-95C5-001A4BA81F0A.root',
                '/store/mc/Summer08/Z_2l2jets/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0013/9E99686D-A2F2-DD11-AD7D-001EC9D8BDC3.root',
                '/store/mc/Summer08/Z_2l2jets/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0013/5645BA49-A2F2-DD11-BB51-001EC9D7FF57.root',
                '/store/mc/Summer08/Z_2l2jets/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0013/340E570B-9DF2-DD11-963F-00093D11B2AE.root',
                '/store/mc/Summer08/Z_2l2jets/GEN-SIM-RECO/IDEAL_V11_redigi_v1/0013/1E072690-A2F2-DD11-97CF-0030487CDA4C.root'

                        )
                   )

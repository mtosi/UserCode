#!/bin/bash

cd $CMSSW_BASE/src/Configuration/GenProduction/python/

HIGGSMASSLIST="130 150 160 170 180 190 200 220 240 260 280 300 320 340 360 370 380 400 420 440 460 480 500 520 540 560 580 600 700 720 740 760 800"

echo '********' GENERATING PYTHIA6_SM_H_ZZ_2l_2jets_mHxxx_10TeV_cff.py '********'
echo w/ xxx'='$HIGGSMASSLIST
for k in ${HIGGSMASSLIST}; do 
    sed -e 's/<HIGGSMASS>/'${k}'/g' PYTHIA6_SM_H_ZZ_2l_2jets_mHTEMPLATE_10TeV_cff.py > PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_10TeV_cff.py
done
echo '**************************************************************************'

scramv1 b --reset
echo '**************************************************************************'

CONDITION="FrontierConditions_GlobalTag,IDEAL_31X::All"
EVENTCONTENT="RECO" #"RAWSIM" "AODSIM"

echo condition $CONDITION
echo eventcontent $EVENTCONTENT

echo '********' GENERATING PYTHIA6_SM_H_ZZ_2l_2jets_mHxxx_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_IDEAL.py '********'
for k in ${HIGGSMASSLIST}; do 

        echo `ls PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_10TeV_cff.py`

### FullSim ###
	cmsDriver.py Configuration/GenProduction/python/PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_10TeV_cff.py \
	-s GEN,SIM,DIGI,L1,DIGI2RAW,HLT \
	--eventcontent $EVENTCONTENT \
	--datatier GEN-SIM-RAW \
	--conditions $CONDITION \
	-n 10 \
	--no_exec

	sed -i -e 's/_py_/_/' PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_IDEAL.py
done
echo '****************************************************************************************************************'

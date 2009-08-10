#!/bin/bash

if [ "$1" == "" ]; then
    echo 'no Higgs mass value is been specified'
    echo '   => got the default list of values' 
    HIGGSMASSLIST="130 150 160 170 180 190 200 220 240 260 280 300 320 340 360 370 380 400 420 440 460 480 500 520 540 560 580 600 700 720 740 760 800"
else
    HIGGSMASSLIST=$1
fi
echo 'Higgs mass value: ' $HIGGSMASSLIST

CONDITION="FrontierConditions_GlobalTag,MC_31X_V1::All"
EVENTCONTENT="RAWSIM"
DATATIER="GEN-SIM-RAW"
echo condition $CONDITION
echo eventcontent $EVENTCONTENT
echo datatier $DATATIER

cd $CMSSW_BASE/src/Configuration/GenProduction/python/

scramv1 b --reset
echo '**************************************************************************'

for k in 130; do 
#for k in ${HIGGSMASSLIST}; do 
    DIROUT="H${k}ZZllqq/"

    echo '********' GENERATING PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_7TeV_SIM_DIGI_L1_DIGI2RAW_HLT_MC31X.py '********'
    echo `ls ${DIROUT}PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_7TeV_GEN.root`

    cmsDriver.py SIM_DIGI_L1_DIGI2RAW_HLT_${k}_7TeV.py \
	-s SIM,DIGI,L1,DIGI2RAW,HLT \
	--eventcontent $EVENTCONTENT \
	--datatier $DATATIER \
	--conditions $CONDITION \
	--processName HLT \
	--filein file:${DIROUT}PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_7TeV_GEN_MC31X.root \
	--fileout PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_7TeV_SIM_RAW_MC31X.root \
	--python_filename ${DIROUT}PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_7TeV_SIM_RAW_MC31X_cfg.py \
	--mc \
	--dirout $DIROUT 
	-n 10 \
	--no_exec

    sed -i -e 's/SIM_DIGI_L1_DIGI2RAW_HLT/PYTHIA6_SM_H_ZZ_2l_2jets_mH/g' ${DIROUT}PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_7TeV_SIM_RAW_MC31X_cfg.py

done
echo '****************************************************************************************************************'


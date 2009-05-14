#!/bin/bash

if [ "$1" == "" ]; then
    echo 'no Higgs mass value is been specified'
    echo '   => got the default list of values' 
    HIGGSMASSLIST="130 150 160 170 180 190 200 220 240 260 280 300 320 340 360 370 380 400 420 440 460 480 500 520 540 560 580 600 700 720 740 760 800"
else
    HIGGSMASSLIST=$1
fi
echo 'Higgs mass value: ' $HIGGSMASSLIST

CONDITION="FrontierConditions_GlobalTag,IDEAL_31X::All"
EVENTCONTENT="RECOSIM"
DATATIER="GEN-SIM-RECO"
echo condition $CONDITION
echo eventcontent $EVENTCONTENT
echo datatier $DATATIER

cd $CMSSW_BASE/src/Configuration/GenProduction/python/

scramv1 b --reset
echo '**************************************************************************'

for k in 130; do 
#for k in ${HIGGSMASSLIST}; do 
    DIROUT="H${k}ZZllqq/"

    echo '********' GENERATING PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_10TeV_RECO_IDEAL.py '********'
    echo `ls ${DIROUT}PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_10TeV_SIM_RAW_IDEAL.root`

    cmsDriver.py Reconstruction_${k}_10TeV.py \
	-s RAW2DIGI,RECO \
	--eventcontent $EVENTCONTENT \
	--datatier $DATATIER \
	--conditions $CONDITION \
	--filein file:${DIROUT}PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_10TeV_SIM_RAW_IDEAL.root \
	--fileout file:PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_10TeV_RECO_IDEAL.root \
	--python_filename ${DIROUT}PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_10TeV_RECO_IDEAL_cfg.py \
	--mc \
	--dirout $DIROUT \
	-n 10 \
	--no_exec

    sed -i -e 's/Reconstruction_/PYTHIA6_SM_H_ZZ_2l_2jets_mH/g' ${DIROUT}PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_10TeV_RECO_IDEAL_cfg.py
done
echo '****************************************************************************************************************'


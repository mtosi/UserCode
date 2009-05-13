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
EVENTCONTENT="RAWSIM"
DATATIER="GEN-SIM-RAW"
echo condition $CONDITION
echo eventcontent $EVENTCONTENT
echo datatier $DATATIER

cd $CMSSW_BASE/src/Configuration/GenProduction/python/


for k in ${HIGGSMASSLIST}; do 
    echo '********' GENERATING PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_10TeV_cff.py '********'
    sed -e 's/<HIGGSMASS>/'${k}'/g' PYTHIA6_SM_H_ZZ_2l_2jets_mHTEMPLATE_10TeV_cff.py > PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_10TeV_cff.py
done
echo '**************************************************************************'

scramv1 b --reset
echo '**************************************************************************'

for k in 130; do 
#for k in ${HIGGSMASSLIST}; do 
    echo '********' GENERATING PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_10TeV_cff_py_GEN_SIM_RAW_IDEAL.py '********'
    echo `ls PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_10TeV_cff.py`

### FullSim ###
	cmsDriver.py Configuration/GenProduction/python/PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_10TeV_cff.py \
	-s GEN:ProductionFilterSequence,SIM,DIGI,L1,DIGI2RAW,HLT \
	--eventcontent $EVENTCONTENT \
	--datatier $DATATIER \
	--conditions $CONDITION \
	-n 10 \
	--no_exec

	mv PYTHIA6_SM_H_ZZ_2l_2jets_mH130_10TeV_cff_py_GEN_SIM_DIGI_L1_DIGI2RAW_HLT_IDEAL.py PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_10TeV_cff_py_GEN_SIM_RAW_IDEAL_cfg.py
	sed -i -e 's/_py_/_/' PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_10TeV_cff_py_GEN_SIM_RAW_IDEAL_cfg.py
done
echo '****************************************************************************************************************'


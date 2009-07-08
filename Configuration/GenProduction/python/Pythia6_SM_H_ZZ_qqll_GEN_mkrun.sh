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
DATATIER="GEN"
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

#for k in 130; do 
for k in ${HIGGSMASSLIST}; do 
    echo '********' GENERATING PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_10TeV_GEN_MC31X.py '********'
    echo `ls PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_10TeV_cff.py`

    DIROUT="H${k}ZZllqq/"

    mkdir $DIROUT

    cmsDriver.py Configuration/GenProduction/python/PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_10TeV_cff.py \
	-s GEN:ProductionFilterSequence \
	--eventcontent $EVENTCONTENT \
	--datatier $DATATIER \
	--conditions $CONDITION \
        --fileout PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_10TeV_GEN_MC31X.root \
	--python_filename ${DIROPUT}PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_10TeV_GEN_MC31X_cfg.py \
	--mc \
	--dirout $DIROUT \
	-n 10 \
	--no_exec 

#    cvs add PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_10TeV_GEN_MC31X_cfg.py
    cp PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_10TeV_GEN_MC31X_cfg.py ${DIROUT}
    cp PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_10TeV_cff.py* ${DIROUT}
done
echo '****************************************************************************************************************'

#cvs commit

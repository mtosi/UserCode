#!/bin/bash

if [ "$1" == "" ]; then
    echo 'no Higgs mass value is been specified'
    echo '   => got the default list of values' 
    HIGGSMASSLIST="130 150 160 170 180 190 200 220 240 260 280 300 320 340 360 370 380 400 420 440 460 480 500 520 540 560 580 600 700 720 740 760 800"
else
    HIGGSMASSLIST=$1
fi
echo 'Higgs mass value: ' $HIGGSMASSLIST

OUTPUT=output

for k in ${HIGGSMASSLIST}; do 
    DIR="H${k}ZZllqq/"
    FILE=${DIR}${OUTPUT}

    cmsRun ${DIR}PYTHIA6_SM_H_ZZ_2l_2jets_mH${k}_10TeV_GEN_IDEAL_cfg.py > ${FILE}

    xSec=`cat ${FILE} | grep "All included subprocesses" | cut -d 'I' -f 4`
    exp=`cat ${FILE} | grep "All included subprocesses" | cut -d 'I' -f 4 | cut -d '-' -f 2`
    num=`cat ${FILE} | grep "All included subprocesses" | cut -d 'I' -f 4 | cut -d 'E' -f 1`
    echo Higgs mass ${k} xSec: $xSec mb '(xE-9 pb)' 
done

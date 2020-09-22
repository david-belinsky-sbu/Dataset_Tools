#!/bin/bash

source config.sh

echo 'extracting data'
python extract_patches.py ${DATAMANIFEST} ${SVSLOC} ${OUTPUTDIR} ${REMOTEUSER} ${REMOTEKEY} ${ANNOTDATA} ${TUMORTYPE} ${HIGH} ${MED} ${LOW} ${SEED} ${VALIDDIR}
#echo 'splitting data'
#python Folder_Arrangement_Like_CoNSeP.py ${OUTPUTDIR}

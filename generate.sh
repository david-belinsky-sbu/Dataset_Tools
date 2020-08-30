#!/bin/bash

source config.sh

echo 'extracting data'
python extract_patches.py ${DATAMANIFEST} ${SVSLOC} ${OUTPUTDIR}
#echo 'splitting data'
#python Folder_Arrangement_Like_CoNSeP.py ${OUTPUTDIR}

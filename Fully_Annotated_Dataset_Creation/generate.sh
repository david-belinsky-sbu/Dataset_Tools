#!/bin/bash

source config.sh

echo 'generating data'
python feature_extraction_from_manual_annotation.py ${DATADUMP} ${DATAMANIFEST} ${SVSLOC} ${SEGDATA} ${OUTPUTDIR}
echo 'splitting data'
python Folder_Arrangement_Like_CoNSeP.py ${OUTPUTDIR}

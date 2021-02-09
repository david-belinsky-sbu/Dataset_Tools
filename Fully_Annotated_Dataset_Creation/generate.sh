#!/bin/bash

source config.sh
set -o noglob
echo 'generating data'
python feature_extraction_from_manual_annotation.py ${DATADUMP} ${DATAMANIFEST} ${SVSLOC} ${SVS} ${SEGDATA} ${OUTPUTDIR} ${REMOTE} ${REMOTEUSER} ${REMOTEKEY}
echo 'splitting data'
python Folder_Arrangement_Like_CoNSeP.py ${OUTPUTDIR}
set +o noglob

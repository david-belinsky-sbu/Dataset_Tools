#!/bin/bash

MANIFEST_FILE='/nfs/data01/shared/mahmudul/UPenn:ColonCancer-2020-10-7-14-24-33/manifest.csv' 
WSI_FOLDER_PATH='/nfs/data02/shared/tcga_analysis/upenn-colon/'
OUTPUT_DIR='/nfs/data01/shared/mahmudul/UPennColonTrainPatchV2/'
RESOLUTION=40.0
PATCH_SIZE=400
JSON_FOLDER='/nfs/data01/shared/mahmudul/UPenn:ColonCancer-2020-10-7-14-24-33/'

python patch_extractor.py ${MANIFEST_FILE} ${WSI_FOLDER_PATH} ${OUTPUT_DIR} ${RESOLUTION} ${PATCH_SIZE} ${JSON_FOLDER}

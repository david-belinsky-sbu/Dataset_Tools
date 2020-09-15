#!/bin/bash

source config.sh

for files in ${IN_FOLDER}/*.*; do
    SVS=`echo ${files} | awk -F'/' '{print $NF}'`

    python save_svs_to_tiles.py ${SVS} ${IN_FOLDER} ${OUT_FOLDER} ${PATCH_SIZE}
    
    if [ $? -ne 0 ]; then
        echo "failed extracting patches for " ${SVS}
        rm -rf ${OUT_FOLDER}/${SVS}
    else
        touch ${OUT_FOLDER}/${SVS}/extraction_done.txt
    fi

done

exit 0;


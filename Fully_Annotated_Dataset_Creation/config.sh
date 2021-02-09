#!/bin/bash

DATADUMP=../../seer:rutgers:lung-2021-1-16-6-38-34/
DATAMANIFEST=manifest.csv

SEGMENT=1

#SEERLUNGCONFIG
SVS=/data/images/
SVSLOC=/data/quip_distro/images/
SEGDATA='/data11/bwang/run*/batch0*/seg_tiles/'

#TCGABRCACONFIG
#SVS=/data/images/tcga_data/brca/
#SVSLOC=/data03/tcga_data/tumor/brca/
#SEGDATA=../../segdata/

OUTPUTDIR=../../seertest3/

REMOTE=1
REMOTEUSER=belinsky
REMOTEKEY=~/quipkey

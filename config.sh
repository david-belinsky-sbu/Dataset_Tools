#!/bin/bash

WORKDIR=../super-patches-evaluation-merged/
DATAFILE=super-patches-label_m.csv
DATAMANIFEST=$WORKDIR$DATAFILE
SVSLOC="/data{}/tcga_data/tumor/"
REMOTEUSER=dbelinsky
REMOTEKEY=~/215key
#SEGDATA=../segdata/

OUTPUTFOLDER=super-patches/
OUTPUTDIR=$WORKDIR$OUTPUTFOLDER

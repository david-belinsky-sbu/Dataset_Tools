#!/bin/bash

WORKDIR=../../tumor_superpatches2/
DATAFILE=tumor_manifest.csv
DATAMANIFEST=$WORKDIR$DATAFILE
SVSLOC="/data{}/tcga_data/tumor/"
REMOTEUSER=dbelinsky
REMOTEKEY=~/215key

#ANNOTDATA should be a folder containing the binary tumor mask pngs at 2x magnification.
ANNOTDATA=../../tumordata/
TUMORTYPE='brca'

HIGH=10
MED=10
LOW=10
SEED=1

OUTPUTFOLDER=super-patches/
VALIDFOLDER=annot-location/
OUTPUTDIR=$WORKDIR$OUTPUTFOLDER
VALIDDIR=$WORKDIR$VALIDFOLDER

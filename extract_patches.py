import numpy as np
import openslide
import sys
import os
from PIL import Image
import datetime
import time
import cv2
from shutil import copyfile as cp
import multiprocessing as mp
import random
import pandas as pd
import json
import glob
import scipy.io as sio
import matplotlib.pyplot as plt
import fnmatch


svsloc = sys.argv[2]
print(svsloc)
outputDir = sys.argv[3]

datafile = sys.argv[1]
manifest = pd.read_csv(datafile)

def find(pattern, path):
    #print(path)
    #print(pattern)
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                return os.path.join(path, name)
    raise NameError('SVS not found')

patch_size_40X = 1600
level = 0
missing_slides = []
missing_tumors = []

if not os.path.exists(outputDir): os.mkdir(outputDir)
    
for index, row in manifest.iterrows():    
    slide_name = row["Slide"]
    patchloc = row["Patch"]
    patchlocs = patchloc.split('_')
    tumor = row["Tumor"]
    output_file = os.path.join(outputDir, slide_name + '_' + patchloc + '.png' )
    tumorloc = os.path.join(svsloc, tumor)
    for i in range(12):
        tumordir = tumorloc.format(str(i).zfill(2))
        if os.path.exists(tumordir):
            #tumordir = os.path.join(tumordir, tumor)
            break
        



    try:
        if not os.path.exists(tumordir):
            if tumor not in missing_tumors:
                missing_tumors.append(tumor)
            raise ValueError("Tumor location nonexistant")

        svsfile = find(slide_name + '*.svs', tumordir)
        oslide = openslide.OpenSlide(svsfile);
        if openslide.PROPERTY_NAME_MPP_X in oslide.properties:
            mag = 10.0 / float(oslide.properties[openslide.PROPERTY_NAME_MPP_X]);
        elif "XResolution" in oslide.properties:
            mag = 10.0 / float(oslide.properties["XResolution"]);
        elif 'tiff.XResolution' in oslide.properties:   # for Multiplex IHC WSIs, .tiff images
            Xres = float(oslide.properties["tiff.XResolution"])
            if Xres < 10:
                mag = 10.0 / Xres;
            else:
                mag = 10.0 / (10000/Xres)       # SEER PRAD
        else:
            print('[WARNING] mpp value not found. Assuming it is 40X with mpp=0.254!', slide_name);
            mag = 10.0 / float(0.254);
    
        scale = 40.0 / mag;  # scale patch size from 'mag' to 40x
        pw = int(patch_size_40X / scale)
    
        width = oslide.dimensions[0];
        height = oslide.dimensions[1];

        if int(patchlocs[0]) + pw > width:
            pw_x = width - int(patchlocs[0])
        else:
            pw_x = pw
        if int(patchlocs[1]) + pw > height:
            pw_y = height - int(patchlocs[1])
        else:
            pw_y = pw

        patch = oslide.read_region((int(patchlocs[0]), int(patchlocs[1])), 0, (pw_x, pw_y));
        patch.resize((int(patch_size_40X * pw_x / pw), int(patch_size_40X * pw_y / pw)), Image.ANTIALIAS)
        patch.save(output_file)
        oslide.close()
        print('extracting {}'.format(output_file));
    except NameError:
        missing_slides.append(slide_name)
        print('Find Failed: {}'.format(slide_name))
    except:
        print('{}: exception caught'.format(slide_name));


with open('errors.txt', 'w') as f:
    f.write('Missing Tumors: ')
    f.write(str(missing_tumors))
    f.write('\n')
    f.write('Missing Slides: ')
    f.write(str(missing_slides))
    f.write('\n')

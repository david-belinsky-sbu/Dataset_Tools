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
import subprocess


svsloc = sys.argv[2]
print(svsloc)
outputDir = sys.argv[3]
remotefile = False

datafile = sys.argv[1]
remoteuser = sys.argv[4]
remotekey = sys.argv[5]
annotdata = sys.argv[6]
tumor = sys.argv[7]
high = int(sys.argv[8])
med = int(sys.argv[9])
low = int(sys.argv[10])
seed = int(sys.argv[11])
random.seed(seed)
#manifest = {'Slide': [], 'Patch': [], 'Tumor': [], 'Annotation': [], 'Annotation Discrete': [], 'Native Scale': []}
manifest = pd.DataFrame(columns = ['Slide', 'Patch', 'Tumor', 'Annotation', 'Annotation Discrete', 'Native Scale'])
#manifest = pd.read_csv(datafile)
manifest["Native_Scale"] = 0

def downRemoteFile(path, remotefile):
    basename = path.split('/')
    basename = basename[-1]
    print(path)
    print("Downloading File: {}".format(path), end=' ')
    p = subprocess.Popen(["scp", '-i', remotekey, '{}@129.49.254.215:{}'.format(remoteuser, path), basename])
    sts = os.waitpid(p.pid, 0)
    print("Done")
    remotefile = True
    return basename, remotefile

def find(pattern, path, remotefile):
    #print(path)
    #print(pattern)
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                return os.path.join(path, name), remotefile
    return downRemoteFile('/home/tcga_all/{}/{}'.format(path.split('/')[-1], pattern), remotefile)

patch_size_40X = 1600
level = 0
missing_slides = []
missing_tumors = []
tumorloc = os.path.join(svsloc, tumor)
for i in range(12):
    tumordir = tumorloc.format(str(i).zfill(2))
    if os.path.exists(tumordir):
        break


if not os.path.exists(outputDir): os.makedirs(outputDir)

file_list = glob.glob(annotdata + '*.png')
    
for fname in file_list:    
    slide_name = os.path.basename(fname)
    slide_name = slide_name[:-4] + '.svs'
    #output_file = os.path.join(outputDir, slide_name + '_' + patchloc + '.png' )
    #tumorloc = os.path.join(svsloc, tumor)

    try:
        if not os.path.exists(tumordir):
            if tumor not in missing_tumors:
                missing_tumors.append(tumor)
            svsfile, remotefile = downRemoteFile('/home/tcga_all/{}/{}'.format(tumor, slide_name), remotefile)
        else:
            svsfile, remotefile = find(slide_name, tumordir, remotefile)
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

        highcount = 0
        medcount = 0
        lowcount = 0
        xmax = width // pw
        ymax = height // pw

        annot = cv2.imread(fname)
        pwannot = pw / 20
        print('annotation shape: {}'.format(str(annot.shape)))
        print('original shape: {}, {}'.format(width, height))
        print('annotpwidth: {}, patchwidth: {}'.format(pwannot, pw))

        while highcount < high or medcount < med or lowcount < low:
            #Get random Patches
            x = random.randint(0, xmax - 1)
            y = random.randint(0, ymax - 1)
            

        #if int(patchlocs[0]) + pw > width:
        #    pw_x = width - int(patchlocs[0])
        #else:
        #    pw_x = pw
        #if int(patchlocs[1]) + pw > height:
        #    pw_y = height - int(patchlocs[1])
        #else:
        #    pw_y = pw

        patch = oslide.read_region((int(patchlocs[0]), int(patchlocs[1])), 0, (pw_x, pw_y));
        if pw_x != pw and pw_y != pw:
            patch.resize((int(patch_size_40X * pw_x / pw), int(patch_size_40X * pw_y / pw)), Image.ANTIALIAS)
        else:
            patch.resize((int(patch_size_40X), int(patch_size_40X)), Image.ANTIALIAS)
        patch.save(output_file)
        manifest.at[index, 'Native_Scale'] = mag
        oslide.close()
        print('extracting {}'.format(output_file));
        if remotefile:
            print("Removing {}".format(svsfile))
            os.remove(svsfile)
            remotefile = False
    #except NameError:
    #    missing_slides.append(slide_name)
    #    print('Find Failed: {}'.format(slide_name))
    except Exception as e:
        print('{}: exception caught\nMessage: {}'.format(slide_name, str(e)));


with open('errors.txt', 'w') as f:
    f.write('Missing Tumors: ')
    f.write(str(missing_tumors))
    f.write('\n')
    f.write('Missing Slides: ')
    f.write(str(missing_slides))
    f.write('\n')

manifest.to_csv(datafile)

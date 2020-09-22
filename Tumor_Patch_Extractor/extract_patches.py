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
validDir = sys.argv[12]
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
manifest.set_index(['Slide', 'Patch'], drop=False)
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

def whiteness(patch):
    wh = np.mean(np.std(patch, axis=(0,1)))
    return wh

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
if not os.path.exists(validDir): os.makedirs(validDir)

file_list = glob.glob(annotdata + '*.png')
    
for fname in file_list:
    #seed += 1    
    slide_name = os.path.basename(fname)
    slide_name = slide_name[:-4] + '.svs'
    print(slide_name)
    #output_file = os.path.join(outputDir, slide_name + '_' + patchloc + '.png' )
    #tumorloc = os.path.join(svsloc, tumor)

    try:
    #if True:
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
        pw = int(round((patch_size_40X / scale)))
        width = oslide.dimensions[0];
        height = oslide.dimensions[1];
        highcount = 0
        medcount = 0
        lowcount = 0
        xmax = width // pw
        ymax = height // pw
        slidename = slide_name.split('.')[-3]
        annot = cv2.imread(fname, 0)
        annotrgb = np.stack((annot, annot, annot), axis=-1)
        annotscale = width / annot.shape[1]
        pwannot = int(round(pw / annotscale))
        #print('annotation shape: {}'.format(str(annot.shape)))
        #print('original shape: {}, {}'.format(width, height))
        #print('annotpwidth: {}, patchwidth: {}'.format(pwannot, pw))
        #random.seed(seed)
        dupecount = 0
        randomcount = 0
        while highcount < high or medcount < med or lowcount < low:
            #Get random Patches
            if dupecount > 50:
                break
            if randomcount > 10000:
                break
            randomcount += 1
            x = random.randint(0, xmax - 1)
            y = random.randint(0, ymax - 1)
            annoty = (y * pwannot) + 1
            annotx = (x * pwannot) + 1
            annotpatch = annot[annoty: annoty + pwannot, annotx: annotx + pwannot]
            totaltumor = np.sum(annotpatch != 0)
            annotport = (totaltumor * 1.0) / (pwannot * pwannot)
            try:
            #if True:
                if annotport > 0.75:
                    if highcount < high:
                        patch = oslide.read_region(((x*pw)+1, (y*pw)+1), 0, (pw, pw));
                        patcharr = np.array(patch)
                        if (patcharr[:,:,3].max() == 0):
                            #print('blank: Patch: {}_{}'.format((x*pw)+1, (y*pw)+1))
                            continue
                        elif whiteness(patcharr[:,:,:3]) < 12:
                            #print('white: Patch: {}_{}'.format((x*pw)+1, (y*pw)+1))
                            continue
                        patch.resize((int(patch_size_40X), int(patch_size_40X)), Image.ANTIALIAS)
                        manifest = manifest.append(pd.DataFrame({'Slide': [slidename], 'Patch': ['{}_{}'.format((x*pw)+1, (y*pw)+1)], 'Tumor': [tumor], 'Annotation': [annotport], 'Annotation Discrete': [1], 
                                                            'Native Scale': [mag]}).set_index(["Slide", "Patch"], drop=False), verify_integrity=True)
                        patch.save('{}{}_{}_{}.png'.format(outputDir, slidename, (x*pw)+1, (y*pw)+1))
                        validpatch = np.copy(annotrgb)
                        validpatch[annoty: annoty + pwannot, annotx: annotx + pwannot] = np.array([255, 0, 0])
                        cv2.imwrite('{}{}_{}_{}.png'.format(validDir, slidename, (x*pw)+1, (y*pw)+1), validpatch)
                        highcount += 1
                        print('highcount: {}, Patch: {}_{}'.format(highcount, (x*pw)+1, (y*pw)+1))
                elif annotport < 0.25:
                    if lowcount < low:
                        patch = oslide.read_region(((x*pw)+1, (y*pw)+1), 0, (pw, pw));
                        patcharr = np.array(patch)
                        if (patcharr[:,:,3].max() == 0):
                            #print('blank: Patch: {}_{}'.format((x*pw)+1, (y*pw)+1))
                            continue
                        elif whiteness(patcharr[:,:,:3]) < 12:
                            #print('white: Patch: {}_{}'.format((x*pw)+1, (y*pw)+1))
                            continue
                        patch.resize((int(patch_size_40X), int(patch_size_40X)), Image.ANTIALIAS)
                        manifest = manifest.append(pd.DataFrame({'Slide': [slidename], 'Patch': ['{}_{}'.format((x*pw)+1, (y*pw)+1)], 'Tumor': [tumor], 'Annotation': [annotport], 'Annotation Discrete': [0],
                                                            'Native Scale': [mag]}).set_index(["Slide", "Patch"], drop=False), verify_integrity=True)
                        patch.save('{}{}_{}_{}.png'.format(outputDir, slidename, (x*pw)+1, (y*pw)+1))
                        validpatch = np.copy(annotrgb)
                        validpatch[annoty: annoty + pwannot, annotx: annotx + pwannot] = np.array([255, 0, 0])
                        cv2.imwrite('{}{}_{}_{}.png'.format(validDir, slidename, (x*pw)+1, (y*pw)+1), validpatch)
                        lowcount += 1
                        print('lowcount: {}, Patch: {}_{}'.format(lowcount, (x*pw)+1, (y*pw)+1))
                else:
                    if medcount < med:
                        patch = oslide.read_region(((x*pw)+1, (y*pw)+1), 0, (pw, pw));
                        patcharr = np.array(patch)
                        if (patcharr[:,:,3].max() == 0):
                            #print('blank: Patch: {}_{}'.format((x*pw)+1, (y*pw)+1))
                            continue
                        elif whiteness(patcharr[:,:,:3]) < 12:
                            #print('white: Patch: {}_{}'.format((x*pw)+1, (y*pw)+1))
                            continue
                        patch.resize((int(patch_size_40X), int(patch_size_40X)), Image.ANTIALIAS)
                        manifest = manifest.append(pd.DataFrame({'Slide': [slidename], 'Patch': ['{}_{}'.format((x*pw)+1, (y*pw)+1)], 'Tumor': [tumor], 'Annotation': [annotport], 'Annotation Discrete': [0.5],
                                                            'Native Scale': [mag]}).set_index(["Slide", "Patch"], drop=False), verify_integrity=True)
                        patch.save('{}{}_{}_{}.png'.format(outputDir, slidename, (x*pw)+1, (y*pw)+1))
                        validpatch = np.copy(annotrgb)
                        validpatch[annoty: annoty + pwannot, annotx: annotx + pwannot] = np.array([255, 0, 0])
                        cv2.imwrite('{}{}_{}_{}.png'.format(validDir, slidename, (x*pw)+1, (y*pw)+1), validpatch)
                        medcount += 1
                        print('medcount: {}, Patch: {}_{}'.format(medcount, (x*pw)+1, (y*pw)+1))
            except Exception as e:
                print('Duplicate patch {}'.format(str(e)))
                dupecount += 1
                continue
        oslide.close()
        if remotefile:
            print("Removing {}".format(svsfile))
            os.remove(svsfile)
            remotefile = False
    except Exception as e:
        print('{}: exception caught\nMessage: {}'.format(slide_name, str(e)));
        remotefile = False


manifest.to_csv(datafile, index=False)

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
import skimage.draw as skd
from misc.viz_utils import visualize_instances, class_colors
from misc.utils import get_inst_centroid, bounding_box
from misc.patch_extractor import PatchExtractor
from metrics.stats_utils import remap_label
import glob
from shapely.geometry import Polygon
import scipy.io as sio
import matplotlib.pyplot as plt


brca = '/data/images/tcga_data/brca/'
brcaloc = sys.argv[3]

outputDir = sys.argv[5]

extract_type = 'mirror'
step_size = [80, 80]
win_size = [540, 540]
xtractor = PatchExtractor(win_size, step_size)
dumppath = sys.argv[1]
polygonpath = sys.argv[4] + 'brca_polygon/{}/*.csv'
manifest = pd.read_csv(os.path.join(dumppath, sys.argv[2]))
manifest["imagepath"] = manifest["imagepath"].str.replace(brca, brcaloc, regex=False)

def annotDict(annot):
    if annot == "DOTS-Tumor":
        return 2
    elif annot == "DOTS-Lymphs:plasmacells":
        return 1
    else:
        return 3

def inbound(bound, coord): # This function checks whether the dots from JSON lies into the boundary box. Here bound = box coordinates and coord =  dots coordinate
    return coord[0] > bound[0][0] and coord[1] > bound[0][1] and coord[0] <= bound[1][0] and coord[1] <= bound[1][1]

def filematch(filelist, corner): # This function finds the segmented polygon csv file where the box lies. Filelist =  contains the segmented polygons csv files and corner = top left corner 4k box which contains the 1k annotated box   
    for f in filelist:
        if corner in f:
            return f

        
def applypolygons(inst, region, mask, offset): #The function draws segmented polygons on the annotated box. Here, inst = 1k np.zero array for getting instance_map, region = contains the dots coordinates and types, mask = file_name of the segmented polygons, 
    masks = pd.read_csv(mask)
    bound = Polygon([(region["bound"][0][1] + 1, region["bound"][0][0] + 1), 
                     (region["bound"][0][1] + 1, region["bound"][1][0]), 
                     (region["bound"][1][1], region["bound"][1][0]), 
                     (region["bound"][1][1], region["bound"][0][0] + 1)])
    lastidx = 0
    for idx, instance in masks.iterrows():
        points = instance["Polygon"][1:-1]
        points = points.split(':')
        points = list(map(float, points))
        points = list(map(int, points))
        points = [(points[i+1], points[i]) for i in range(0, len(points), 2)]
        poly = Polygon(points)
        try:
            if poly.intersects(bound) and not poly.touches(bound):
                if poly.within(bound):
                    xx, yy = skd.polygon(np.array([int(x[0]) - (region["bound"][0][1] + 1) for x in points]), np.array([int(x[1]) - (region["bound"][0][0] + 1) for x in points]))
                    inst[xx, yy] = idx + offset
                    lastidx = idx + offset
                else:
                    if not poly.is_valid:
                        poly = poly.buffer(0)
                    poly = poly.intersection(bound)
                    if poly.geom_type != 'Polygon':
                        coll = list(poly)
                        for i in coll:
                            if i.geom_type == 'Polygon':
                                poly = i
                                break
                    if poly.geom_type != 'Polygon':
                        print(poly)
                    else:
                        points = list(poly.exterior.coords)
                        xx, yy = skd.polygon(np.array([int(x[0]) - (region["bound"][0][1] + 1) for x in points]), np.array([int(x[1]) - (region["bound"][0][0] + 1) for x in points]))
                        inst[xx, yy] = idx + offset
                        lastidx = idx + offset
        except:
            if not poly.is_valid:
                poly = poly.buffer(0)
            if poly.intersects(bound) and not poly.touches(bound):
                if poly.within(bound):
                    xx, yy = skd.polygon(np.array([int(x[0]) - (region["bound"][0][1] + 1) for x in points]), np.array([int(x[1]) - (region["bound"][0][0] + 1) for x in points]))
                    inst[xx, yy] = idx + offset
                    lastidx = idx + offset
                else:
                    if not poly.is_valid:
                        poly = poly.buffer(0)
                    poly = poly.intersection(bound)
                    if poly.geom_type != 'Polygon':
                        coll = list(poly)
                        for i in coll:
                            if i.geom_type == 'Polygon':
                                poly = i
                                break
                    if poly.geom_type != 'Polygon':
                        print(poly)
                    else:
                        points = list(poly.exterior.coords)
                        xx, yy = skd.polygon(np.array([int(x[0]) - (region["bound"][0][1] + 1) for x in points]), np.array([int(x[1]) - (region["bound"][0][0] + 1) for x in points]))
                        inst[xx, yy] = idx + offset
                        lastidx = idx + offset

    return inst, lastidx

pw = 1000
level = 0

if not os.path.exists(outputDir): os.mkdir(outputDir)
    
for index, row in manifest.iterrows():    
    slide_name = row["imagepath"].split('/')[-1]
    output_folder = os.path.join(outputDir, slide_name)
    
    try:
        oslide = openslide.OpenSlide(row["imagepath"]);
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

        scale = 20.0 / mag;  # scale patch size from 'mag' to 40x
        patch_size_40X = int(pw * scale)

        width = oslide.dimensions[0];
        height = oslide.dimensions[1];
        oslide.close()
    except:
        print('{}: exception caught'.format(slide_name));
        continue;

    polygonloc = polygonpath.format(slide_name)          
    file_list = glob.glob(polygonloc) # This filelist contains the coordinates of the segmentation polygon
    if len(file_list)<=0: continue
    file_list.sort() 
    
    if not os.path.exists(output_folder) : os.mkdir(output_folder)
    
    fdone = '{}/extraction_done.txt'.format(output_folder);
    if os.path.isfile(fdone):
        print('fdone {} exist, skipping'.format(fdone));
        continue;

    print('extracting {}'.format(output_folder));
    print('height/width: {}/{}'.format(height, width))

    
    with open(os.path.join(dumppath, row["path"])) as f: #Loading JSON from Manifest file. JSONs contain all annotations
        data = json.load(f)

    corrs = []
    regions = []
    badsegs = []
    nosegs = []
    splitsegs = []
    points = []

    for annot in data:
        if annot["properties"]["annotations"]["notes"] == "500p-Tumor": # This block loads the box coordinates (1000x1000)
            patch = annot["geometries"]["features"][0]["geometry"]["coordinates"][0]
            patch = [[int(coor[0] * width), int(coor[1] * height)] for coor in patch]
            maxx = max(patch, key=lambda x: x[0])[0] + 500
            maxy = max(patch, key=lambda x: x[1])[1] + 500
            minx = min(patch, key=lambda x: x[0])[0]
            miny = min(patch, key=lambda x: x[1])[1]
            patch = { "bound": [[minx, miny], [maxx, maxy]], "points": [], "badsegs": [], "nosegs": [], "splitsegs": []} # Mins are top left corner and Maxs are right bottom corner of the 1k box
            regions.append(patch)
        elif annot["properties"]["annotations"]["notes"] == "DOTS-Bad Segmentation":
            badseg = annot["geometries"]["features"][0]["geometry"]["coordinates"]
            badseg[0] = int(badseg[0] * width)
            badseg[1] = int(badseg[1] * height)
            badsegs.append(badseg)
        elif annot["properties"]["annotations"]["notes"] == "DOTS-split cell":
            splitseg = annot["geometries"]["features"][0]["geometry"]["coordinates"]
            splitseg[0] = int(splitseg[0] * width)
            splitseg[1] = int(splitseg[1] * height)
            splitsegs.append(splitseg)
        elif annot["properties"]["annotations"]["notes"] == "DOTS-non-seg object":
            noseg = annot["geometries"]["features"][0]["geometry"]["coordinates"]
            noseg[0] = int(noseg[0] * width)
            noseg[1] = int(noseg[1] * height)
            nosegs.append(noseg)
        elif 'Prostate' in annot["properties"]["annotations"]["notes"]:
            continue
        elif 'DOTS' not in annot["properties"]["annotations"]["notes"]:
            continue
        else: # This block loads the nuclei type annotation with coordinates
            point = { "type": annot["properties"]["annotations"]["notes"],
                      "coor": annot["geometries"]["features"][0]["geometry"]["coordinates"]
                    }
            try:
                point["coor"] = [int(point["coor"][0] * width), int(point["coor"][1] * height)]
            except:
                print(point["coor"])
                print(point["type"])
            points.append(point)
    
    for point in points:
        for region in regions:
            if inbound(region["bound"], point["coor"]):
                region["points"].append(point)
    for badseg in badsegs:
        for region in regions:
            if inbound(region["bound"], badseg):
                region["badsegs"].append(badseg)
    for noseg in nosegs:
        for region in regions:
            if inbound(region["bound"], noseg):
                region["nosegs"].append(noseg)
    for splitseg in splitsegs:
        for region in regions:
            if inbound(region["bound"], splitseg):
                region["splitsegs"].append(splitseg)
    
    matches = 0
    totals = 0
    b = 0
    segremoved = 0
    nucleicount = [0, 0, 0]
    inst_sizes = []
    inst_types = []
    inst_lengths = []
    for region in regions:
        instance_map = np.zeros((region["bound"][1][1] - region["bound"][0][1], region["bound"][1][0] - region["bound"][0][0]), np.int32)
        previd = 0
        corner = "/{}_{}_".format(((region["bound"][0][0]-1) // 4000)*4000 + 1, ((region["bound"][0][1]-1) // 4000)*4000 + 1)
        mask = filematch(file_list, corner) # Get the 4k box with segmented polygon
        instance_map, previd = applypolygons(instance_map, region, mask, previd)
        # Following ifs are for to handle 1k box that spans over multiple 4k box
        if (region["bound"][0][1]-1) // 4000 < (region["bound"][1][1]-1) // 4000:
            corner = "/{}_{}_".format(((region["bound"][0][0]-1) // 4000)*4000 + 1, ((region["bound"][1][1]-1) // 4000)*4000 + 1)
            mask = filematch(file_list, corner) 
            instance_map, previd = applypolygons(instance_map, region, mask, previd)
        if (region["bound"][0][0]-1) // 4000 < (region["bound"][1][0]-1) // 4000:
            corner = "/{}_{}_".format(((region["bound"][1][0]-1) // 4000)*4000 + 1, ((region["bound"][0][1]-1) // 4000)*4000 + 1)
            mask = filematch(file_list, corner)
            instance_map, previd = applypolygons(instance_map, region, mask, previd)
        if (region["bound"][0][0]-1) // 4000 < (region["bound"][1][0]-1) // 4000 and (region["bound"][0][1]-1) // 4000 < (region["bound"][1][1]-1) // 4000:
            corner = "/{}_{}_".format(((region["bound"][1][0]-1) // 4000)*4000 + 1, ((region["bound"][1][1]-1) // 4000)*4000 + 1)
            mask = filematch(file_list, corner)
            instance_map, previd = applypolygons(instance_map, region, mask, previd)

        splitmap = np.zeros(instance_map.shape, np.int32)
        for splitseg in region["splitsegs"]:
            inst_id = instance_map[splitseg[1] - (region["bound"][0][1] + 1)][splitseg[0] - (region["bound"][0][0] + 1)]
            if inst_id > 0:
                splitmap = np.logical_or(splitmap, instance_map == inst_id)
        splitinst = instance_map * splitmap
        downmap = np.pad(splitinst, ((3,0),(0,0)), mode='constant')[:-3, :]
        downmiss = splitinst != downmap
        downremap = np.logical_and(splitinst[downmiss] != 0, downmap[downmiss] != 0)
        orig = splitinst[downmiss][downremap]
        new = downmap[downmiss][downremap]
        while len(orig) != 0:
            for i in range(len(orig)):
                #print('merge {} into {}'.format(orig[i], new[i]))
                instance_map[instance_map == orig[i]] = new[i]
                orig[orig == orig[i]] = new[i]
            splitinst = instance_map * splitmap
            downmap = np.pad(splitinst, ((3,0),(0,0)), mode='constant')[:-3, :]
            downmiss = splitinst != downmap
            downremap = np.logical_and(splitinst[downmiss] != 0, downmap[downmiss] != 0)
            orig = splitinst[downmiss][downremap]
            new = downmap[downmiss][downremap]
        leftmap = np.pad(splitinst, ((0,0),(3,0)), mode='constant')[:, :-3]
        leftmiss = splitinst != leftmap
        leftremap = np.logical_and(splitinst[leftmiss] != 0, leftmap[leftmiss] != 0)
        orig = splitinst[leftmiss][leftremap]
        new = leftmap[leftmiss][leftremap]
        while len(orig) != 0:
            for i in range(len(orig)):
                #print('merge {} into {}'.format(orig[i], new[i]))
                instance_map[instance_map == orig[i]] = new[i]
                orig[orig == orig[i]] = new[i]
            splitinst = instance_map * splitmap
            leftmap = np.pad(splitinst, ((0,0),(3,0)), mode='constant')[:, :-3]
            leftmiss = splitinst != leftmap
            leftremap = np.logical_and(splitinst[leftmiss] != 0, leftmap[leftmiss] != 0)
            orig = splitinst[leftmiss][leftremap]
            new = leftmap[leftmiss][leftremap]
        rightmap = np.pad(splitinst, ((0,0),(0,2)), mode='constant')[:, 2:] 
        leftmap = np.pad(splitinst, ((0,0),(2,0)), mode='constant')[:, :-2]
        lrbridge = np.logical_and(instance_map == 0, np.logical_and(np.logical_and(rightmap != 0, leftmap != 0), rightmap == leftmap))
        instance_map[lrbridge] = rightmap[lrbridge]
        splitmap[lrbridge] = 1
        splitinst = instance_map * splitmap
        upmap = np.pad(splitinst, ((0,2),(0,0)), mode='constant')[2:, :]
        downmap = np.pad(splitinst, ((2,0),(0,0)), mode='constant')[:-2, :]
        udbridge = np.logical_and(instance_map == 0, np.logical_and(np.logical_and(upmap != 0, downmap != 0), upmap == downmap))
        instance_map[udbridge] = upmap[udbridge]
 
        instance_map = remap_label(instance_map, by_size=True) # This resets instance number from 0
        type_map = np.zeros(instance_map.shape, np.int32)
        instance_list = list(np.unique(instance_map))[1:] # Background 0 excluded from instance list
        inst_type = np.full(len(instance_list), 0, dtype=np.int32)
        inst_type2 = np.full(len(instance_list), 0, dtype=np.int32)
        inst_size = np.full(len(instance_list), 0, dtype=np.int32)
        inst_length = np.full(len(instance_list), 0, dtype=np.float32)
        inst_centroids = get_inst_centroid(instance_map)
        for point in region["points"]:
            annot = annotDict(point["type"])
            inst_id = instance_map[(point["coor"][1] - (region["bound"][0][1] + 1), point["coor"][0] - (region["bound"][0][0] + 1))]
            if inst_id > 0:
                #print("{}: {}".format(inst_id, annot))
                if inst_type[inst_id - 1] != annot:
                    nucleicount[annot - 1] += 1
                inst_type[inst_id - 1] = annot
                inst_type2[inst_id - 1] = annot
                type_map[instance_map == inst_id] = annot
                matches += 1
            else:
                inst_type = np.append(inst_type, [annotDict(point["type"])], 0)
                inst_centroids = np.append(inst_centroids, [[point["coor"][1] - (region["bound"][0][1] + 1), point["coor"][0] - (region["bound"][0][0] + 1)]], 0)
        totals += len(instance_list)
        b +=len(region["points"])
        
        
        save_dir = os.path.join(output_folder, 'Images')
        if not os.path.exists(save_dir) : os.mkdir(save_dir)
        oslide = openslide.OpenSlide(row["imagepath"])
        patch = oslide.read_region((region["bound"][0][0] + 1, region["bound"][0][1] + 1), 0, (region["bound"][1][0] - region["bound"][0][0], region["bound"][1][1] - region["bound"][0][1]));
        oslide.close()
        patch = np.array(patch)
        patch = cv2.cvtColor(np.float32(patch), cv2.COLOR_RGBA2BGR)
        patchjoin = np.copy(patch)
        patchcorr = np.copy(patch)
        #patchpoints = np.copy(patch)
        for badseg in region["badsegs"]:
            inst_id = instance_map[badseg[1] - (region["bound"][0][1] + 1)][badseg[0] - (region["bound"][0][0] + 1)]
            if inst_id > 0:
                bmatch = False
                inst_type2[inst_id - 1] = 0
                for point in region["points"]:
                    if instance_map[(point["coor"][1] - (region["bound"][0][1] + 1), point["coor"][0] - (region["bound"][0][0] + 1))] == inst_id:
                        if not bmatch:
                            bmatch = True
                            inst_type[inst_id - 1] = annotDict(point["type"])
                            inst_centroids[inst_id - 1] = (point["coor"][1] - (region["bound"][0][1] + 1), point["coor"][0] - (region["bound"][0][0] + 1))
                        else:
                            inst_type = np.append(inst_type, [annotDict(point["type"])], 0)
                            inst_centroids = np.append(inst_centroids, [[point["coor"][1] - (region["bound"][0][1] + 1), point["coor"][0] - (region["bound"][0][0] + 1)]], 0)
                patchjoin[instance_map == inst_id] = 255
                patchcorr[instance_map == inst_id] = 255
                type_map[instance_map == inst_id] = 0
                instance_map[instance_map == inst_id] = 0
            patchcorr = cv2.circle(patchcorr, (badseg[0] - (region["bound"][0][0] + 1), badseg[1] - (region["bound"][0][1] + 1)), 5, class_colors[6])
        for point in region["points"]:
            annot = annotDict(point["type"])
            patchcorr = cv2.circle(patchcorr, (point["coor"][0] - (region["bound"][0][0] + 1), point["coor"][1] - (region["bound"][0][1] + 1)), 5, class_colors[annot])
        for splitseg in region["splitsegs"]:
            patchcorr = cv2.circle(patchcorr, (splitseg[0] - (region["bound"][0][0] + 1), splitseg[1] - (region["bound"][0][1] + 1)), 5, class_colors[5])
        for noseg in region["nosegs"]:
            patchcorr = cv2.circle(patchcorr, (noseg[0] - (region["bound"][0][0] + 1), noseg[1] - (region["bound"][0][1] + 1)), 5, class_colors[7])
        px, py = (int(patch.shape[0] * scale), int(patch.shape[1] * scale))
        patch = cv2.resize(patch, (px, py), interpolation=cv2.INTER_NEAREST);
        patchcorr = cv2.resize(patchcorr, (px, py), interpolation=cv2.INTER_NEAREST);
        patchjoin = cv2.resize(patchjoin, (px, py), interpolation=cv2.INTER_NEAREST);
        instance_map = np.float32(instance_map)
        type_map = np.float32(type_map)
        #cv2.imwrite('{}/{}_inst.png'.format(save_dir, 'sample'), instance_map)
        #cv2.imwrite('{}/{}_type.png'.format(save_dir, 'sample'), type_map * 60)
        instance_map = cv2.resize(instance_map, (px, py), interpolation=cv2.INTER_NEAREST)
        type_map = cv2.resize(type_map, (px, py), interpolation=cv2.INTER_NEAREST)
        instance_map = np.int32(instance_map)
        type_map = np.int32(type_map)
        instances, instcounts = np.unique(instance_map, return_counts=True)
        for i in range(1, len(instances)):
            inst_size[instances[i]-1] = instcounts[i]
            inst_map = np.array(instance_map == instances[i], np.uint8)
            y1, y2, x1, x2  = bounding_box(inst_map)
            y1 = y1 - 2 if y1 - 2 >= 0 else y1
            x1 = x1 - 2 if x1 - 2 >= 0 else x1
            x2 = x2 + 2 if x2 + 2 <= inst_map.shape[1] - 1 else x2
            y2 = y2 + 2 if y2 + 2 <= inst_map.shape[0] - 1 else y2
            inst_map = inst_map[y1:y2, x1:x2]
            contours = cv2.findContours(inst_map, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
            #inst_size[instances[i]-1] = cv2.contourArea(contours[0])
            inst_length[instances[i]-1] = cv2.arcLength(contours[0][0], True)
        indexer = np.logical_and(inst_type2 != 0, inst_size != 0)
        inst_sizes = np.append(inst_sizes, inst_size[indexer])
        inst_lengths = np.append(inst_lengths, inst_length[indexer])
        inst_types = np.append(inst_types, inst_type2[indexer])
        
        basename = '{}_{}_{}_{}_{}'.format(region["bound"][0][0] + 1, region["bound"][0][1] + 1, region["bound"][1][0] - region["bound"][0][0], region["bound"][1][1] - region["bound"][0][1], scale)
        cv2.imwrite('%s/%s.png' % (save_dir, basename), patch)
        save_dir = os.path.join(output_folder, 'Labels')
        if not os.path.exists(save_dir) : os.mkdir(save_dir)
        inst_centroids *= scale
        sio.savemat('%s/%s.mat' % (save_dir, basename),
                        {'inst_map'  :     instance_map,
                            'type_map'  :     type_map,
                            'inst_type' :     inst_type[:, None],
                            'inst_centroid' : inst_centroids,
                        })
 
        save_dir = os.path.join(output_folder, 'Combined')
        if not os.path.exists(save_dir) : os.mkdir(save_dir)
        colormap = [class_colors[x] for x in inst_type]
        overlaid_output = visualize_instances(instance_map, patchcorr, colormap) 
        cv2.imwrite("{}/{}.png".format(save_dir, basename), overlaid_output)
        #type_map[np.logical_and(instance_map > 0, type_map == 0)] = -1
        unannot = np.logical_and(instance_map != 0, type_map == 0)
        patchjoin[unannot] = 255
        instance_map[unannot] = 0
        ann = np.dstack([instance_map, type_map])
        patchjoin = np.array(patchjoin, np.int32)
        print(ann.shape)
        print(patchjoin.shape)
        patchjoin = np.concatenate([patchjoin, ann], axis=-1)
        sub_patches = xtractor.extract(patchjoin, extract_type) # This extracts mirrored subpatch of 540x540 from 1k patch
        for idx, patch in enumerate(sub_patches):
            np.save("{}/{}_{}.npy".format(save_dir, basename, idx), patch)
    
    f = open(fdone, 'w')
    t = len(points)
    f.write("Total Number of Dots annotated are: {}\n".format(t))
    f.write("Number of Dots that fall outside of boxes are: {}\n".format(t-b))
    f.write("Percentage of the dots falling into the boxes: {}%\n\n".format((b/t)*100))
    f.write("Number of automated segmented polygon instance: {}\n".format(totals))
    f.write("Number of dots annotated falling inside the boxes: {}\n".format(b))
    f.write("Number of dots that matched polygons: {}\n".format(matches))
    f.write("Percentage of manual dots vs automated polygon: {}%\n\n".format((matches/totals)*100))
    f.write("Lymphocyte Count: {}\n".format(nucleicount[0]))
    f.write("Tumor Count: {}\n".format(nucleicount[1]))
    f.write("Misc. Count: {}\n\n".format(nucleicount[2]))
    f.write("Cell Measurements (Avg Area, Std. Dev Area, Avg Perimter, Std. Dev Perimeter)\n")
    f.write("Lymphocyte: ({},{},{},{})\n".format(np.mean(inst_sizes[inst_types==1]), np.std(inst_sizes[inst_types==1]),np.mean(inst_lengths[inst_types==1]),np.std(inst_lengths[inst_types==1])))
    f.write("Tumor: ({},{},{},{})\n".format(np.mean(inst_sizes[inst_types==2]), np.std(inst_sizes[inst_types==2]),np.mean(inst_lengths[inst_types==2]),np.std(inst_lengths[inst_types==2])))
    f.write("Misc.: ({},{},{},{})\n".format(np.mean(inst_sizes[inst_types==3]), np.std(inst_sizes[inst_types==3]),np.mean(inst_lengths[inst_types==3]),np.std(inst_lengths[inst_types==3])))
    #f.write("Stroma: ({},{},{},{})\n".format(np.mean(inst_sizes[inst_types==4]), np.std(inst_sizes[inst_types==4]),np.mean(inst_lengths[inst_types==4]),np.std(inst_lengths[inst_types==4])))
    #f.write("Misc: ({},{},{},{})\n".format(np.mean(inst_sizes[inst_types==5]), np.std(inst_sizes[inst_types==5]),np.mean(inst_lengths[inst_types==5]),np.std(inst_lengths[inst_types==5])))
    f.write("Overall: ({},{},{},{})\n".format(np.mean(inst_sizes), np.std(inst_sizes),np.mean(inst_lengths),np.std(inst_lengths)))
    f.close()

import sys
import pandas as pd
import os
import openslide
import json
from PIL import Image
from functools import reduce
import numpy as np
import cv2
import os.path
from os import path
import time
import multiprocessing as mp



manifest_file = sys.argv[1] #This csv contains which json file corresponds to which slides.
previous_imagepath = '/data/images/upenn-2020/colon/' 
wsi_folder = sys.argv[2]
outputDir = sys.argv[3]
if not os.path.exists(outputDir) : os.mkdir(outputDir)
resolution = float(sys.argv[4])
patch_size = int(sys.argv[5])
json_folder = sys.argv[6] #JSONS files contain manual annotation
manifest = pd.read_csv(manifest_file)
manifest["imagepath"] = manifest["imagepath"].str.replace(previous_imagepath, wsi_folder, regex=False)

# Following dictonary corresponds to different types of cell mapped to numerical values
labels_to_number_dict = {'Tumor':1,'Dysplastic epithelium':2,'NON-Tumor':3,'Lymphs:plasma cells':4,'Normal epithelium':5,'Stroma':6,'Neutrophils':7,'Vessels':8,'Nerves':9,'Necrosis':10,'Curly Collagen':11,'Straight Collagen':12}
unique_labels = [0 for i in range(13)]

for index, row in manifest.iterrows():
    start = time.time() 
    slide_name = row["imagepath"].split('/')[-1]
    output_folder = os.path.join(outputDir, slide_name) # This outputfolder is the folder named with the slide name where patches of that WSI are stored
    try:
        oslide = openslide.OpenSlide(row["imagepath"]);
        if openslide.PROPERTY_NAME_MPP_X in oslide.properties:
            mag = 10.0 / float(oslide.properties[openslide.PROPERTY_NAME_MPP_X]);
        elif "XResolution" in oslide.properties:
            mag = 10.0 / float(oslide.properties["XResolution"]);
        elif 'tiff.XResolution' in oslide.properties:
            Xres = float(oslide.properties["tiff.XResolution"])
            if Xres < 10:
                mag = 10.0 / Xres;
            else:
                mag = 10.0 / (10000/Xres)
        else:
            print('[WARNING] mpp value not found. Assuming it is 40X with mpp=0.254!', slide_name);
            mag = 10.0 / float(0.254);

#        scale = mag / resolution; # This is used for magnifying the slide according to user choice
        scale = 1 # This is used for native resolution
        pw = int(patch_size * scale)
        width = oslide.dimensions[0];
        height = oslide.dimensions[1];
    except:
        print('{}: exception caught'.format(slide_name));
        continue;

    if not os.path.exists(output_folder) : os.mkdir(output_folder)

    fdone = '{}/extraction_done.txt'.format(output_folder);
    if os.path.isfile(fdone):
        print('fdone {} exist, skipping'.format(fdone));
        continue;

    print('extracting {}'.format(output_folder));
    print('height/width: {}/{}'.format(height, width))

    coors = []
    final_coors = []
    final_coors_with_label = []
    for x in range(1, width, pw):
        for y in range(1, height, pw):
            if x + pw > width:
                pw_x = width - x
            else:
                pw_x = pw
            if y + pw > height:
                pw_y = height - y
            else:
                pw_y = pw

            if int(patch_size * pw_x / pw) > 50 and int(patch_size * pw_y / pw) > 50:
                coors.append((x,y,pw_x,pw_y))

    with open(os.path.join(json_folder, row["path"])) as f:
        data = json.load(f)

    for annot in data:
        unique_labels[labels_to_number_dict[annot["properties"]["annotations"]["notes"][15:]]] += 1 
        polygon_coord = annot["geometries"]["features"][0]["geometry"]["coordinates"][0]
        min_x = int(min(polygon_coord, key=lambda x:x[0])[0]*width)
        min_y = int(min(polygon_coord, key=lambda x:x[1])[1]*height)
        max_x = int(max(polygon_coord, key=lambda x:x[0])[0]*width)
        max_y = int(max(polygon_coord, key=lambda x:x[1])[1]*height)
        for coor  in coors:
            x1,y1,x2,y2 = coor
            x2 = x1+x2
            y2 = y1+y2
            if min_x>=x1 and min_y>=y1 and min_x<=x2 and min_y<=y2:
                min_x = x1
                min_y = y1
            if max_x>=x1 and max_y>=y1 and max_x<=x2 and max_y<=y2:
                max_x = x2
                max_y = y2

        polygon_coord = [[int(i[0]*width)-min_x,int(i[1]*height)-min_y] for i in polygon_coord]
        polygon_coord = np.array([polygon_coord],dtype=np.int32)
        big_patch = np.zeros((max_y-min_y,max_x-min_x),dtype=np.uint8)
        cv2.fillPoly(big_patch,polygon_coord,labels_to_number_dict[annot["properties"]["annotations"]["notes"][15:]])
        for coor in coors:
            x1,y1,x2,y2 = coor
            pw_x = x2
            pw_y = y2
            x2 = x1+x2
            y2 = y1+y2
            if x1>=min_x and y1>=min_y and x2<=max_x and y2<=max_y and (x1,y1,x2,y2) not in final_coors:
                m_x1 = x1 - min_x
                m_y1 = y1 - min_y
                m_x2 = x2 - min_x
                m_y2 = y2 - min_y
                m_min_x = min(m_x1,m_x2)
                m_min_y = min(m_y1,m_y2)
                m_max_x = max(m_x1,m_x2)
                m_max_y = max(m_y1,m_y2)
                mask = big_patch[m_min_y:m_max_y,m_min_x:m_max_x] 
                label = np.unique(mask,return_counts=True)
                # For Labeling we have labeled the patch with  whichever nuclei has the maximum frequency, if there is no specific type of nuclei we labeled it as a background.  
                if len(label[0])>1 and label[0][np.argmax(label[1])]==0:
                    label_numbers = label[0]
                    frequency = label[1]
                    frequency, label_numbers = zip(*sorted(zip(frequency,label_numbers)))
                    label = label_numbers[-2]
                else:
                    label = label[0][0] 
                final_coors.append((x1,y1,x2,y2))
                final_coors_with_label.append((x1,y1,x2,y2,label))
 
    def extract_patch(corr):
        x, y, x2, y2,label = corr
        pw_x = x2-x
        pw_y = y2-y
        fname = '{}/{}_{}_{}_{}_{}_{}_{}.png'.format(output_folder, x, y, pw, patch_size, pw_x, pw_y, label)
        patch = oslide.read_region((x, y), 0, (pw_x, pw_y));
        patch_arr = np.array(patch);
        patch = patch.resize((int(patch_size * pw_x / pw), int(patch_size * pw_y / pw)), Image.ANTIALIAS);
        patch.save(fname);
        fname = '{}_{}_{}_{}_{}_{}_{}.png'.format(x, y, pw, patch_size, pw_x, pw_y, label)
        return fname,label

    f = open(fdone, 'w')
    pool = mp.Pool(processes=8)
    for result in pool.imap(extract_patch, final_coors_with_label):
        f.write('%s %d\n'%result)

    oslide.close()
    end = time.time()
    print("{} took {}".format(slide_name,(end-start)/60)) 
    f.close()
print(unique_labels)

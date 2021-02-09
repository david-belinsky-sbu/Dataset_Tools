import glob
import os
import shutil
import random
import sys

random.seed(14)
base = sys.argv[1]
files = glob.glob(sys.argv[1] + "*.tif")
random.shuffle(files)

split_ratio = .8
split = int(len(files) * split_ratio)

train_files = files[:split]
valid_files = files[split:]


train_folder = os.path.join(base,'train/')
valid_folder = os.path.join(base,'valid/')
test_folder = os.path.join(base,'test/')



if not os.path.exists(train_folder) : os.mkdir(train_folder)
if not os.path.exists(valid_folder) : os.mkdir(valid_folder)
if not os.path.exists(test_folder) : os.mkdir(test_folder)
    
train_folder = os.path.join(train_folder,'540x540_80x80/')
valid_folder = os.path.join(valid_folder,'540x540_80x80/')
#test_folder = os.path.join(test_folder,'Images/')

if not os.path.exists(train_folder) : os.mkdir(train_folder)
if not os.path.exists(valid_folder) : os.mkdir(valid_folder)
if not os.path.exists(os.path.join(test_folder,'Images/')) : os.mkdir(os.path.join(test_folder,'Images/'))
if not os.path.exists(os.path.join(test_folder,'Labels/')) : os.mkdir(os.path.join(test_folder,'Labels/'))
if not os.path.exists(os.path.join(test_folder,'Combined/')) : os.mkdir(os.path.join(test_folder,'Combined/'))


def move_numpy(filelist, destFolder):
    for t in filelist:
        i = 0
        combined = glob.glob(os.path.join(t,'Combined/*.npy'))
        for src in combined:
            dest = destFolder + "{}_{}.npy".format(t[-19:-13],i)
            new_path = shutil.copy(src, dest)
            print(new_path)
            i+=1



def move_image(filelist, destFolder):
    for t in filelist:
        i = 0
        combined = glob.glob(os.path.join(t,'Images/*.png'))
        for src in combined:
            dest = destFolder+ 'Images/' + "{}.png".format(t[-19:-13])
            new_path = shutil.copy(src, dest)
            print(new_path)
            i+=1
        combined = glob.glob(os.path.join(t, 'Labels/*.mat'))
        for src in combined:
            dest = destFolder + 'Labels/' + "{}.mat".format(t[-19:-13])
            new_path = shutil.copy(src, dest)
            print(new_path)
            i+=1
        combined = glob.glob(os.path.join(t, 'Combined/*.png'))
        for src in combined:
            dest = destFolder + 'Combined/' + "{}.png".format(t[-19:-13])
            new_path = shutil.copy(src, dest)
            print(new_path)
            i+=1


move_numpy(train_files, train_folder)
move_numpy(valid_files, valid_folder)
move_image(valid_files, test_folder)





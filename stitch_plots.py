import uproot
import pandas as pd
import numpy as np
import argparse
import itertools
import os, sys
from icecream import ic
import matplotlib
matplotlib.use('Agg') 

import matplotlib.pyplot as plt
from copy import copy
from utils.utils import dot
from utils.utils import mag
from utils.utils import mag2
from utils.utils import cosTheta
from utils.utils import angle
from utils.utils import cross
from utils.utils import vecAdd
from utils.utils import pi0Energy
from utils.utils import pi0InvMass
from utils.utils import getPhi
from utils.utils import getTheta
from utils.utils import getEnergy
from utils.utils import readFile
from utils import make_histos
from utils import histo_plotting
from utils import filestruct
from convert_root_to_pickle import convert_GEN_NORAD_root_to_pkl
from convert_root_to_pickle import convert_GEN_RAD_root_to_pkl
from convert_root_to_pickle import convert_RECON_NORAD_root_to_pkl
from convert_root_to_pickle import convert_RECON_RAD_root_to_pkl
import pickle_analysis
from root2pickleEpggRec import root2pickle
pd.set_option('mode.chained_assignment', None)

import random 
import sys
import os, subprocess
import argparse
import shutil
import time
from datetime import datetime 
import json

import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt 
import random 
import sys
import os, subprocess
from pdf2image import convert_from_path
import math
from icecream import ic
import shutil
from PIL import Image, ImageDraw, ImageFont


M = 0.938272081 # target mass
me = 0.5109989461 * 0.001 # electron mass
ebeam = 10.604 # beam energy
pbeam = np.sqrt(ebeam * ebeam - me * me) # beam electron momentum
beam = [0, 0, pbeam] # beam vector
target = [0, 0, 0] # target vector
alpha = 1/137 #Fund const
mp = 0.938 #Mass proton
prefix = alpha/(8*np.pi)
E = 10.6


def img_from_pdf(img_dir):
    image_files = []
    lists = os.listdir(img_dir)
    sort_list = sorted(lists)
    #for f in sort_list:
    #    print(f)
    #sys.exit()
    for img_file in sort_list:
        #print("On file " + img_file)
        image1 = Image.open(img_dir+img_file)
        image_files.append(image1)

    return image_files

def append_images(images, xb_counter, direction='horizontal', 
                  bg_color=(255,255,255), aligment='center'):
    
    # Appends images in horizontal/vertical direction.

    # Args:
    #     images: List of PIL images
    #     direction: direction of concatenation, 'horizontal' or 'vertical'
    #     bg_color: Background color (default: white)
    #     aligment: alignment mode if images need padding;
    #        'left', 'right', 'top', 'bottom', or 'center'

    # Returns:
    #     Concatenated image as a new PIL image object.
    
    widths, heights = zip(*(i.size for i in images))

    if direction=='horizontal':
        new_width = sum(widths)
        new_height = max(heights)
    else:
        new_width = max(widths)
        new_height = sum(heights)

    new_im = Image.new('RGB', (new_width, new_height), color=bg_color)

    if direction=='vertical':
        new_im = Image.new('RGB', (int(new_width+0), int(new_height+images[0].size[1]/2)), color=bg_color)


    offset = 0
    for im_counter,im in enumerate(reversed(images)):
        ic(im_counter)
        if direction=='horizontal':
            y = 0
            if aligment == 'center':
                y = int((new_height - im.size[1])/2)
            elif aligment == 'bottom':
                y = new_height - im.size[1]
            new_im.paste(im, (offset, y))
            offset += im.size[0]
        else:
            x = 0
            if aligment == 'center':
                x = int((new_width - im.size[0])/2)
            elif aligment == 'right':
                x = new_width - im.size[0]
            new_im.paste(im, (x, offset))
            offset += im.size[1]

    return new_im


def chunks(l, n):
	spits = (l[i:i+n] for i in range(0, len(l), n))
	return spits




fs = filestruct.fs()


allz = True
TEST = True
#All
if allz:
    DoWipe = True
    DoGen = True
    DoRecon = True
    DoInspect = True
    DoBin = True
    DoCombine = True
    DoMetaCombine = True
#Just ana
else:
    DoWipe = False
    DoGen = False
    DoRecon = False
    DoInspect = False
    DoBin = False
    DoCombine = False
    DoMetaCombine = True



base_dir = "/mnt/d/GLOBUS/CLAS12/simulations/production/F2018_In_Norad"
# To add data to this base dir, make a subdir /runs/<run_number>/roots and put the files there
# Note that dataset 101 is just for testing
if base_dir[-1] != '/':
    base_dir += '/'
if TEST:
    base_plot_dir = base_dir + "plots_test/"
else:
    base_plot_dir = base_dir + "plots/"

if TEST:
    plot_out_dir = base_dir + "combined_plots_test/"
else:
    plot_out_dir = base_dir + "combined_plots/"

if TEST:
    binned_pkl = "rec_and_gen_binned_events_meta_test2_binning.pkl"
else:
    binned_pkl = "rec_and_gen_binned_events_meta.pkl"


df = pd.read_pickle(base_dir + binned_pkl)
df = df[["Q2","xB","t1","phi","rec_sum","gen_sum"]]


t_bins = df['t1'].unique()
q_bins = df['Q2'].unique()
x_bins = df['xB'].unique()
p_bins = df['phi'].unique()

TEST2 = False
if TEST2:
    t_bins = [t_bins[0]]

print(t_bins)
for t_bin in t_bins:
    #figures = os.listdir(base_plot_dir+"/"+str(t_bin))

    images = img_from_pdf(base_plot_dir+"/"+str(t_bin)+"/")


    num_ver_slices = len(q_bins)
    num_hori_slices = len(x_bins)

    layers = []
    for i in range(0,num_hori_slices):
        layer = list(images[i*num_ver_slices:i*num_ver_slices+num_ver_slices])
        #print(layer)
        layers.append(layer)

    layers = reversed(layers)
    horimg = []

    for xb_counter,layer in enumerate(layers):
        print("len of layers is {}".format(len(layer)))
        print("counter is {}".format(xb_counter))
        print("On vertical layer {}".format(xb_counter))
        #print(layer)
        imglay = append_images(layer, -1, direction='vertical')
        #imglay.save("testing1.jpg")
        horimg.append(imglay)

    #print(horimg) 
    #ssreversed(horimg)   
    print("Joining images horizontally")
    final = append_images(horimg, 0,  direction='horizontal')
    final_name = "joined_pictures_t_{}.jpg".format(t_bin)
    final.save(plot_out_dir+final_name,optimize=True, quality=100)
    print("saved {}".format(plot_out_dir+final_name))





sys.exit()



#for i in range(0,int(len(images)/num_ver_slices)):

#print(layers)
#sys.exit()
#print(layers[0])








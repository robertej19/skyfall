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


fs = filestruct.fs()


allz = True
TEST = True
TEST2 = False
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



base_dir = "/mnt/d/GLOBUS/CLAS12/simulations/production/F2018_In_Norad/runs"
# To add data to this base dir, make a subdir /runs/<run_number>/roots and put the files there
# Note that dataset 101 is just for testing
if base_dir[-1] != '/':
    base_dir += '/'
if TEST:
    base_plot_dir = base_dir + "plots_test/"
else:
    base_plot_dir = base_dir + "plots/"

if TEST:
    binned_pkl = "rec_and_gen_binned_events_meta.pkl"
else:
    binned_pkl = "rec_and_gen_binned_events_meta.pkl"


df = pd.read_pickle(base_dir + binned_pkl)

df = df[["qmin","xmin","tmin","pmin","rec_sum","gen_sum"]]


df["acc"] = df["rec_sum"]/df["gen_sum"]
df['acc_inv'] = 1/df['acc']
df["acc_inv_err"] = 100*np.sqrt(1/df["rec_sum"] + 1/df["gen_sum"])
df.replace([np.inf, -np.inf], np.nan, inplace=True)
df.replace(np.nan, 0, inplace=True)


t_bins = df['tmin'].unique()
q_bins = df['qmin'].unique()
x_bins = df['xmin'].unique()
p_bins = df['pmin'].unique()

p_labels = []


for p_bin in p_bins:
    p_labels.append(str(p_bin+9))
p_pos = np.arange(len(p_labels))

if DoWipe:
    shutil.rmtree(base_plot_dir) if os.path.isdir(base_plot_dir) else None
for dir_ending in t_bins:
        if not os.path.isdir(base_plot_dir+str(dir_ending)):
            os.makedirs(base_plot_dir+str(dir_ending))

plt.rcParams["font.size"] = "20"

print(df)
for t_bin in t_bins:
    for q_bin in q_bins:
        for x_bin in x_bins:
            query = "tmin=={} and qmin=={} and xmin=={}".format(t_bin, q_bin, x_bin)
            print(query)
            df_over_phi = df.query(query)
            ic(df_over_phi)

            #    ax.errorbar(x_c6,tel_c6,fmt='k',yerr=tel_err_c6,marker='x',linestyle="None", ms=10,label="CLAS6 - t+l")

            # x = df_over_phi['pmin']
            # y = df_over_phi['acc_inv']
            # yerr = df_over_phi['acc_inv_err']

            # fig, ax = plt.subplots(figsize =(16, 10)) 
            # ax.bar(x, y, yerr=yerr, align='center', alpha=0.5, ecolor='black', capsize=10, width=16)
            # ax.set_ylabel('1/Acceptance Correction')
            # ax.set_ylim(0,100)
            # #ax.set_xticks(p_pos[::2])
            # #ax.set_xticklabels(p_labels[::2])
            # ax.set_title('1/Acceptance Correction, $Q^2$ bin: {} GeV$^2$, $x_B$ bin: {}, t bin: {} GeV$^2$'.format(q_bin, x_bin, t_bin))
            # #ax.yaxis.grid(True)

            # if float(q_bin)<10:
            #     q_label = "0"+str(q_bin)

            # print(q_label)

            # plt_title = base_plot_dir+str(t_bin)+"/x_{}_q_{}_acc_inv.png".format(str(x_bin),q_label)
            # plt.savefig(plt_title)

            x = df_over_phi['pmin']
            y = df_over_phi['rec_sum']
            yg = df_over_phi['gen_sum']

            fig, ax = plt.subplots(figsize =(16, 10)) 
            ax.bar(x, y, align='center', alpha=1, color="black", capsize=10, width=16)
            ax.bar(x, yg, align='center', alpha=0.5, color="red", capsize=10, width=16)

            ax.set_ylabel('Simualtions: Gen and Rec Counts')
            #ax.set_ylim(0,100)
            #ax.set_xticks(p_pos[::2])
            #ax.set_xticklabels(p_labels[::2])
            ax.set_title('Simualtions: Gen and Rec Counts, $Q^2$ bin: {} GeV$^2$, $x_B$ bin: {}, t bin: {} GeV$^2$'.format(q_bin, x_bin, t_bin))
            #ax.yaxis.grid(True)

            if float(q_bin)<10:
                q_label = "0"+str(q_bin)

            print(q_label)

            plt_title = base_plot_dir+str(t_bin)+"/x_{}_q_{}_acc_inv.png".format(str(x_bin),q_label)
            plt.savefig(plt_title)


            plt.close()


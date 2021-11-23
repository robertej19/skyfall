import uproot
import pandas as pd
import numpy as np
import argparse
import itertools
import os, sys
from icecream import ic
import matplotlib
#matplotlib.use('Agg') 

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

df = pd.read_pickle("/mnt/d/GLOBUS/CLAS12/simulations/production/F2018_In_Norad/runs/2003_20211112_1732/event_pickles/merged_Fall_2018_Inbending_recon_reconstructed_events_after_cuts.pkl")
#df = pd.read_pickle("/mnt/d/GLOBUS/CLAS12/simulations/production/F2018_In_Norad/runs/100_20211103_1524/event_pickles/100_20211103_1524_merged_Fall_2018_Inbending_recon_reconstructed_events_after_cuts.pkl")


ic(df.head())
print(df.shape)
#for i in df.columns:
#    print(i)

df = df.query("Q2 > 1.5 and Q2 < 2.0 and xB < 0.3 and xB>0.2 and t>0.2 and t<0.3")

print(df.shape)

x_data = df["phi1"]
plot_title = "F 2018 Inbending, epgg, all exclusivity cuts"

#plot_title = "F 2018 Inbending, epgg, no exclusivity cuts"

vars = ["XB (GeV)"]
make_histos.plot_1dhist(x_data,vars,ranges="none",second_x="none",logger=False,first_label="F18IN",second_label="norad",
            saveplot=False,pics_dir="none",plot_title=plot_title,first_color="blue",sci_on=False)
            
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




# xbs = np.arange(start=.2, stop=.3, step=.001)
# q2s = np.arange(start=2, stop=3.5, step=.01)

# xv, yv = np.meshgrid(xbs, q2s)

xv = .2
yv = 3.5
def get_ratios(x,q2):
    epsis = []
    gammas = []
    energies = [5.776, 10.604]
    for e in energies:
        s = e*e+0.938*0.938
        y = q2/(s*x)
        num = 1-y-q2/(4*e*e)
        denom = 1- y + y*y/2 + q2/(4*e*e)
        print(y,q2,e,num,denom)
        epsi = num/denom
        gamma = 1/(e*e)*(1/(1-epsi))
        epsis.append(epsi)
        gammas.append(gamma)
    gamr = gammas[0]/gammas[1]
    epsr = epsis[0]/epsis[1]

    return gamr, epsr

#for x in xbs:
#    for q2 in q2s:

g,e = get_ratios(xv, yv)

print(g)
# g[g<0] =0
# print(g)
# fig, ax = plt.subplots(figsize =(14, 10)) 
# h = plt.pcolormesh(xbs, q2s, g)
# plt.xlim(.2,.3)
# plt.ylim(2,3.5)
# plt.colorbar(h)
# #plt.axis('scaled')
# plt.show()
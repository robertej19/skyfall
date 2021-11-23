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




xbs = np.arange(start=.1, stop=.6, step=.001)
q2s = np.arange(start=1, stop=6, step=.01)

xv, yv = np.meshgrid(xbs, q2s)


def get_ratios(x,q2):
    epsis = []
    gammas = []
    a8p = 1/137*(1/(8*3.14159))
    #print(a8p)
    energies = [5.75, 10.604]
    for e in energies:
        s = 2*0.938*e+0.938*0.938
        y = q2/(s*x)
        num = 1-y-q2/(4*e*e)
        denom = 1- y + y*y/2 + q2/(4*e*e)
        #print(y,q2,e,num,denom)
        epsi = num/denom
        gamma = 1/(e*e)*(1/(1-epsi))*(1-x)/(x*x*x)*a8p*q2/(0.938*.938)
        epsis.append(epsi)
        #print(epsi)
        gammas.append(gamma)
        print(gamma)
    gamr = gammas[0]/gammas[1]
    epsr = epsis[0]/epsis[1]

    return 1/gamr, epsr

# g, e = get_ratios(.2,3.5)
# print(1/g)
#for x in [0.1,0.15,0.2,0.3,0.4,0.5,0.6]:
#     for q2 in [1,2,3,4,5,6]:
#         g,e = get_ratios(x,q2)
#     

#for x in xbs:
#    for q2 in q2s:

g,e = get_ratios(xv,yv)

g[g<0] =0
print(g)
plt.rcParams["font.size"] = "20"

fig, ax = plt.subplots(figsize =(14, 10)) 
h = plt.pcolormesh(xbs, q2s, g)
plt.title('Ratio of Gamma (Clas12/Clas6)')
plt.xlabel('Xb')
plt.ylabel('Q2 (GeV$^2$)')
#h = plt.contour(xv, yv, g, cmap='RdGy');
#plt.xlim(.2,.3)
#plt.ylim(2,3.5)
plt.colorbar(h)
#plt.axis('scaled')
plt.show()
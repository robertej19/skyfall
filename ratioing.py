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
import matplotlib as mpl

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

a = False
b = False
c = False
d = False
e = False
f = False
g = True


if g:
  
    df_minibin = pd.read_pickle("ratioed.pkl")

    #q2bins,xBbins, tbins, phibins = [fs.q2bins_test, fs.xBbins_test, fs.tbins_test, fs.phibins_test]
    q2bins,xBbins, tbins, phibins = [fs.q2bins, fs.xBbins, fs.tbins, fs.phibins]

    #q2bins = [2.0,2.5]
    #xBbins = [0.3,0.25]
    for qmin in q2bins[:-1]:
        for xmin in xBbins[:-1]:
            print(" ON q-{} x-{}".format(qmin, xmin))
            # qmin = 1.5
            # xmin = 0.25
            df = df_minibin.query("qmin==@qmin and xmin==@xmin")
            #for tmin,tmax in zip(tbins[0:-1],tbins[1:]):

            #pave_arr = []

            #for pmin,pmax in zip(phibins[0:-1],phibins[1:]):
            #    pave_arr.append((pmin+pmax)/2)


            for tmin in tbins[0:-1]:
                for pmin,pmax in zip(phibins[0:-1],phibins[1:]):
                    #pave = (pmin+pmax)/2
                    df_t = df.query("tmin==@tmin and pmin==@pmin")
                    #print(df_t)
                    if df_t.shape[0] == 0:
                        #print("APPENDING ZEROS")
                        #dict0 = {'qmin':[qmin],'xmin':[xmin],'tmin':[tmin],'pmin':[pmin],'ratio':['nan'],'ratio_uncert':['nan']}
                        dict0 = {'qmin':[qmin],'xmin':[xmin],'tmin':[tmin],'pmin':[pmin],'ratio':[0],'ratio_uncert':[0]}
                        df2 = pd.DataFrame(dict0)
            #            df = pd.concat([df,df2],ignore_index=True)
                        df = df.append(df2)#,ignore_index=True)

            t = []
            p = []
            r = []


            # for tmin,tmax in zip(tbins[0:-1],tbins[1:]):
            #     tave = (tmin+tmax)/2
            #     #print(tave)
            #     for pave in pave_arr:
            #         df_t = df.query("tmin==@tmin and pave==@pave")
            #         t.append(tave)
            #         p.append(pave)
            #         rval = df_t.ratio.values[0]
            #         print(tave,pave,rval)
            #         r.append(rval)

            for tind,tmin in enumerate(tbins):
                for pind,pmin in enumerate(phibins):
                    print(tmin,pmin)
                    if (tind<len(tbins)-1) and (pind<len(phibins)-1):
                        df_t = df.query("tmin==@tmin and pmin==@pmin")
                        rval = df_t.ratio.values[0]
                    else:
                        rval = 0
                    t.append(tmin)
                    p.append(pmin)
                    r.append(rval)

            # for tmin in tbins[0:-1]:
            #     for pmin in phibins[0:-1]:
            #         #print(tmin,pmin)
            #         df_t = df.query("tmin==@tmin and pmin==@pmin")
            #         print(df_t)
            #         #t.append(tmin)
            #         #p.append(pmin)
            #         rval = df_t.ratio.values[0]
            #         #print(tave,pave,rval)
            #         r.append(rval)
            

            print(t)
            print(p)
            print(len(t))
            print(len(p))
            print(len(r))



            x = np.reshape(p, (len(tbins), len(phibins)))
            y = np.reshape(t, (len(tbins), len(phibins)))
            z = np.reshape(r, (len(tbins), len(phibins)))
            z = np.ma.masked_where(z==0, z)
            cmap = mpl.cm.get_cmap("OrRd").copy()

            cmap.set_bad(color='black')

            print(x)
            print(y)
            print(z)
            fig, ax = plt.subplots(figsize =(36, 17)) 

            plt.rcParams["font.family"] = "Times New Roman"
            plt.rcParams["font.size"] = "20"
            plt.pcolormesh(x,y,z,cmap=cmap)#norm=mpl.colors.LogNorm())


            plt.title("Ratio of CLAS12 to CLAS6 Reduced Cross Sections, Q2 = {}, xB = {}".format(qmin,xmin))
            ax.set_xlabel('Lepton-Hadron Angle')
            ax.set_ylabel('-t (GeV$^2)$')

            plt.colorbar()

# data = np.random.random((4, 4))

# fig, ax = plt.subplots()
# # Using matshow here just because it sets the ticks up nicely. imshow is faster.
# ax.matshow(data, cmap='seismic')
# plt.show()
            
            print(
                "RPINTING Z"
            )
            print(z.shape)
            #z = z[:-1,:-1]
            #print(z.shape)

            plt.rcParams["font.family"] = "Times New Roman"
            plt.rcParams["font.size"] = "10"

            for (i, j), zz in np.ndenumerate(z[2:-5,:-1]):
                print(i,j)
                ii = x[i+2][j]+9
                jj = y[i+2][j]*1.2
                ax.text(ii, jj, '{:0.1f}'.format(zz), ha='center', va='center',
                        bbox=dict(facecolor='white', edgecolor='0.3'))


            plt.ylim([0,2])
            #plt.show()
            #sys.exit()
            plt.savefig("newratios/ratio_q2_{}_xB_{}.png".format(qmin,xmin))
            plt.close()

        #     sys.exit()

        #     plt.contourf(x,y,z, 20, plt.cm.nipy_spectral)
        #     plt.colorbar()
        #     # plt.show()


        #     # ax = fig.add_subplot(111, projection='3d')

        #     # ax.plot_surface(x, y, z)
            
        #     #ax.set_zlabel('CLAS12:CLAS6')
                    

            

        #     # print(df)
        #     # Y = df.tmin
        #     # X = df.pave
        #     # Z = df.ratio
        #     # plt.contourf(X, Y, Z, 20, cmap='RdGy')
        #     # plt.colorbar()
        #     # plt.show()








            





        # #     pd.concat(
        # #     objs,
        # #     axis=0,
        # #     join="outer",
        # #     ignore_index=False,
        # #     keys=None,
        # #     levels=None,
        # #     names=None,
        # #     verify_integrity=False,
        # #     copy=True,
        # # # )
        # #     print(dfout)
        # # sys.exit()


        # # #outname = recon_file.split(".")[0]
        # # #output_loc_event_pkl_after_cuts = dirname+run+"/binned_pickles/"+outname+"_reconstructed_events_after_cuts.pkl"
        # # df = pd.read_pickle(output_loc_event_pkl_after_cuts)
        # # #df = df.query("Q2 > 2 and Q2 < 2.5 and xB < 0.38 and xB>0.3 and t>0.2 and t<0.3")

        # # # print(df.shape)

        # # # x_data = df["phi1"]
        # # # plot_title = "F 2018 Inbending, epgg, all exclusivity cuts"

        # # # #plot_title = "F 2018 Inbending, epgg, no exclusivity cuts"

        # # # vars = ["XB (GeV)"]
        # # # make_histos.plot_1dhist(x_data,vars,ranges="none",second_x="none",logger=False,first_label="F18IN",second_label="norad",
        # # #             saveplot=False,pics_dir="none",plot_title=plot_title,first_color="blue",sci_on=False)

        # # # sys.exit()

        # # df_gen = pd.read_pickle(output_loc_event_pkl_all_gen_events)
        # # #df = pd.read_pickle(save_base_dir+"100_20211103_1524_merged_Fall_2018_Inbending_gen_all_generated_events_all_generated_events.pkl")
        # # for col in df.columns:
        # #     print(col)

        # # df['t1'] = df['t']
        # # orginial_sum = df.shape[0]



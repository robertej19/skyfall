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

a = False
b = False
c = False
d = False
e = False
f = True
g = False

if a:

    df = pd.read_pickle("merged_Fall_2018_Inbending_recon_reconstructed_events_after_cuts.pkl")

    # for col in df.columns:
    #     print(col)
    # print(df_0.shape)

    # df_rec_small = df_0.query('Q2>2.5 and Q2<3 and xB<0.3 and xB>0.25 and t>0.2 and t<0.3')
    # print(df_rec_small.shape)

    # df_rec_small.phi1.plot.hist()
    # plt.show()

    df['t1'] = df['t']
    df_recon = df[["Q2", "W", "xB", "t1", "phi1"]]


    dfs = [df_recon,]

    for index,df0 in enumerate(dfs):
        print("Binning df: {}".format(df))
        prefix = "Gen" if index == 1 else ""

        
        q2bins,xBbins, tbins, phibins = fs.q2bins, fs.xBbins, fs.tbins, fs.phibins
        if True:
                q2bins,xBbins, tbins, phibins = fs.q2bins_test, fs.xBbins_test, fs.tbins_test, fs.phibins_test

        num_counts = []

        qrange = [q2bins[0], q2bins[-1]]
        xBrange = [xBbins[0], xBbins[-1]]
        trange = [tbins[0], tbins[-1]]

        total_num = df0.query('Q2>{} and Q2<{} and xB>{} and xB<{} and t1>{} and t1<{}'.format(*qrange, *xBrange, *trange)).shape[0]
        for qmin,qmax in zip(q2bins[0:-1],q2bins[1:]):
            print("Q2 bin: {} to {}".format(qmin,qmax))                    
            for xmin,xmax in zip(xBbins[0:-1],xBbins[1:]):
                print("xB bin: {} to {}".format(xmin,xmax))
                for tmin,tmax in zip(tbins[0:-1],tbins[1:]):
                    for pmin,pmax in zip(phibins[0:-1],phibins[1:]):
                        if prefix=="Gen":
                            query = "GenQ2 > {} and GenQ2 < {} and GenxB > {} and GenxB < {} and Gent1 > {} and Gent1 < {} and Genphi1 > {} and Genphi1 < {}".format(qmin,qmax,xmin,xmax,tmin,tmax,pmin,pmax)
                        else:  
                            query = "Q2 > {} and Q2 <= {} and xB > {} and xB <= {} and t1 > {} and t1 <= {} and phi1>{} and phi1<={}".format(qmin,qmax,xmin,xmax,tmin,tmax,pmin,pmax)
                        
                        df_bin = df0.query(query)
                        num_counts.append([qmin,xmin,tmin,pmin,len(df_bin.index)])

        df_minibin = pd.DataFrame(num_counts, columns = ['qmin','xmin','tmin','pmin',prefix+'counts'])
        print("Total number of binned events: {}".format(df_minibin[prefix+'counts'].sum()))
        print("Total number of events: {}".format(total_num))

        #df_minibin.to_pickle(dirname+run+"/binned_pickles/"+rec_out_name+prefix+"_reconstructed_events_binned_NEW2.pkl")
        df_minibin.to_pickle("binned_events.pkl")

if b:

    df = pd.read_pickle("merged_Fall_2018_Inbending_recon_reconstructed_events_after_cuts.pkl")
    df_rec_small = df.query('Q2>2.5 and Q2<3 and xB<0.3 and xB>0.25 and t>0.2 and t<0.3')

    pb = [0,18,36,54,72,90,108,126,144,162,180,198,216,234,252,270,288,306,324,342,360]

    for pmin,pmax in zip(pb[0:-1],pb[1:]):
        df_new = df_rec_small.query('phi1>{} and phi1<{}'.format(pmin,pmax))
        print("p bin: {} to {} : {}".format(pmin,pmax,df_new.shape[0]))

        

    x_data = df_rec_small.phi1.values
    vars = ["Phi","Counts"]
    ranges = [0,360,20]
    make_histos.plot_1dhist(x_data,vars,ranges=ranges,second_x="none",logger=False,first_label="rad",second_label="norad",
            saveplot=False,pics_dir="none",plot_title="none",first_color="blue",sci_on=False)

    binned = pd.read_pickle("binned_events.pkl")

    dfb = binned.query('qmin==2.5 and xmin==0.25 and tmin==0.2')
    print(dfb)
    x = dfb.pmin
    y = dfb.counts
    fig, ax = plt.subplots(figsize =(14, 10)) 
    ax.bar(x, y, align='center', alpha=0.5, ecolor='black', capsize=10, width=17)

    plt.show()

if c:
    df = pd.read_pickle("merged_Fall_2018_Inbending_gen_all_generated_events.pkl")

    # for col in df.columns:
    #     print(col)
    print(df.shape)

    df_rec_small = df.query('GenQ2>2.5 and GenQ2<3 and GenxB<0.3 and GenxB>0.25 and Gent1>0.2 and Gent1<0.3')
    print(df_rec_small.shape)

    df_rec_small.Genphi1.plot.hist()
    plt.show()

    df_gen = df[["GenQ2", "GenW", "GenxB", "Gent1", "Genphi1"]]

    dfs = [df_gen,]

    for index,df0 in enumerate(dfs):
        print("Binning df: {}".format(df))
        prefix = "Gen" if index == 0 else ""

        
        q2bins,xBbins, tbins, phibins = fs.q2bins, fs.xBbins, fs.tbins, fs.phibins
        if True:
                q2bins,xBbins, tbins, phibins = fs.q2bins_test, fs.xBbins_test, fs.tbins_test, fs.phibins_test

        num_counts = []

        qrange = [q2bins[0], q2bins[-1]]
        xBrange = [xBbins[0], xBbins[-1]]
        trange = [tbins[0], tbins[-1]]

        total_num = df0.query('GenQ2>{} and GenQ2<{} and GenxB>{} and GenxB<{} and Gent1>{} and Gent1<{}'.format(*qrange, *xBrange, *trange)).shape[0]
        for qmin,qmax in zip(q2bins[0:-1],q2bins[1:]):
            print("Q2 bin: {} to {}".format(qmin,qmax))                    
            for xmin,xmax in zip(xBbins[0:-1],xBbins[1:]):
                print("xB bin: {} to {}".format(xmin,xmax))
                for tmin,tmax in zip(tbins[0:-1],tbins[1:]):
                    for pmin,pmax in zip(phibins[0:-1],phibins[1:]):
                        if prefix=="Gen":
                            query = "GenQ2 > {} and GenQ2 < {} and GenxB > {} and GenxB < {} and Gent1 > {} and Gent1 < {} and Genphi1 > {} and Genphi1 < {}".format(qmin,qmax,xmin,xmax,tmin,tmax,pmin,pmax)
                        else:  
                            query = "Q2 > {} and Q2 <= {} and xB > {} and xB <= {} and t1 > {} and t1 <= {} and phi1>{} and phi1<={}".format(qmin,qmax,xmin,xmax,tmin,tmax,pmin,pmax)
                        
                        df_bin = df0.query(query)
                        num_counts.append([qmin,xmin,tmin,pmin,len(df_bin.index)])

        df_minibin = pd.DataFrame(num_counts, columns = ['qmin','xmin','tmin','pmin',prefix+'counts'])
        print("Total number of binned events: {}".format(df_minibin[prefix+'counts'].sum()))
        print("Total number of events: {}".format(total_num))

        #df_minibin.to_pickle(dirname+run+"/binned_pickles/"+rec_out_name+prefix+"_reconstructed_events_binned_NEW2.pkl")
        df_minibin.to_pickle("binned_events_gen.pkl")


if d:

    df = pd.read_pickle("merged_Fall_2018_Inbending_gen_all_generated_events.pkl")
    df_rec_small = df.query('GenQ2>2.5 and GenQ2<3 and GenxB<0.3 and GenxB>0.25 and Gent1>0.2 and Gent1<0.3')


    pb = [0,18,36,54,72,90,108,126,144,162,180,198,216,234,252,270,288,306,324,342,360]

    for pmin,pmax in zip(pb[0:-1],pb[1:]):
        df_new = df_rec_small.query('Genphi1>{} and Genphi1<{}'.format(pmin,pmax))
        print("p bin: {} to {} : {}".format(pmin,pmax,df_new.shape[0]))

        

    x_data = df_rec_small.Genphi1.values
    vars = ["Phi","Counts"]
    ranges = [0,360,20]
    make_histos.plot_1dhist(x_data,vars,ranges=ranges,second_x="none",logger=False,first_label="rad",second_label="norad",
            saveplot=False,pics_dir="none",plot_title="none",first_color="blue",sci_on=False)

    binned = pd.read_pickle("binned_events_gen.pkl")

    dfb = binned.query('qmin==2.5 and xmin==0.25 and tmin==0.2')
    print(dfb)
    x = dfb.pmin
    y = dfb.Gencounts
    fig, ax = plt.subplots(figsize =(14, 10)) 
    ax.bar(x, y, align='center', alpha=0.5, ecolor='black', capsize=10, width=17)

    plt.show()


if e:
    df = pd.read_pickle("/mnt/d/GLOBUS/CLAS12/simulations/production/F2018_In_Norad/runs/1904_20211112_1057/rec_and_gen_binned_events.pkl")

    print(df)

if f:
    df = pd.read_pickle("F18_All_DVPi0_Events.pkl")
    df["t1"] = df["t"]

    ic(df)
    for col in df.columns:
        print(col)

    df_recon = df[["Q2", "W", "xB", "t1", "phi1"]]


    dfs = [df_recon,]

    for index,df0 in enumerate(dfs):
        print("Binning df: {}".format(df))
        prefix = "Gen" if index == 1 else ""

        
        q2bins,xBbins, tbins, phibins = fs.q2bins, fs.xBbins, fs.tbins, fs.phibins
        #if True:
        #        q2bins,xBbins, tbins, phibins = fs.q2bins_test, fs.xBbins_test, fs.tbins_test, fs.phibins_test

        num_counts = []

        qrange = [q2bins[0], q2bins[-1]]
        xBrange = [xBbins[0], xBbins[-1]]
        trange = [tbins[0], tbins[-1]]

        total_num = df0.query('Q2>{} and Q2<{} and xB>{} and xB<{} and t1>{} and t1<{}'.format(*qrange, *xBrange, *trange)).shape[0]



        for qmin,qmax in zip(q2bins[0:-1],q2bins[1:]):
                    if prefix=="Gen":
                        query = "GenQ2 > {} and GenQ2 < {}".format(qmin,qmax)
                    else:  
                        query = "Q2 > {} and Q2 <= {}".format(qmin,qmax)            
                    df_q = df0.query(query)
                    print("Q2 bin: {} to {}".format(qmin,qmax))                    
                    for xmin,xmax in zip(xBbins[0:-1],xBbins[1:]):
                       # print("xB bin: {} to {}".format(xmin,xmax))
                        if prefix=="Gen":
                            query = "GenQ2 > {} and GenQ2 < {} and GenxB > {} and GenxB".format(qmin,qmax,xmin,xmax)
                        else:  
                            query = "Q2 > {} and Q2 <= {} and xB > {} and xB <= {}".format(qmin,qmax,xmin,xmax)
                        df_qx = df_q.query(query)
                        for tmin,tmax in zip(tbins[0:-1],tbins[1:]):
                            if prefix=="Gen":
                                query = "GenQ2 > {} and GenQ2 < {} and GenxB > {} and GenxB < {} and Gent1 > {} and Gent1 < {}".format(qmin,qmax,xmin,xmax,tmin,tmax)
                            else:  
                                query = "Q2 > {} and Q2 <= {} and xB > {} and xB <= {} and t1 > {} and t1 <= {}".format(qmin,qmax,xmin,xmax,tmin,tmax)
                            df_qxt = df_qx.query(query)
                            for pmin,pmax in zip(phibins[0:-1],phibins[1:]):
                                if prefix=="Gen":
                                    query = "GenQ2 > {} and GenQ2 < {} and GenxB > {} and GenxB < {} and Gent1 > {} and Gent1 < {} and Genphi1 > {} and Genphi1 < {}".format(qmin,qmax,xmin,xmax,tmin,tmax,pmin,pmax)
                                else:  
                                    query = "Q2 > {} and Q2 <= {} and xB > {} and xB <= {} and t1 > {} and t1 <= {} and phi1>{} and phi1<={}".format(qmin,qmax,xmin,xmax,tmin,tmax,pmin,pmax)
                                
                                df_bin = df_qxt.query(query)
                                num_counts.append([qmin,xmin,tmin,pmin,len(df_bin.index)])


        df_minibin = pd.DataFrame(num_counts, columns = ['qmin','xmin','tmin','pmin',prefix+'counts'])
        print("Total number of binned events: {}".format(df_minibin[prefix+'counts'].sum()))
        print("Total number of events: {}".format(total_num))

        #df_minibin.to_pickle(dirname+run+"/binned_pickles/"+rec_out_name+prefix+"_reconstructed_events_binned_NEW2.pkl")
        df_minibin.to_pickle("F18_In_binned_events.pkl")



def get_gamma(x,q2):
    a8p = 1/137*(1/(8*3.14159))
    #print(a8p)
    energies = [10.604]
    for e in energies:
        y = q2/(2*x*e*mp)
        num = 1-y-q2/(4*e*e)
        denom = 1- y + y*y/2 + q2/(4*e*e)
        #print(y,q2,e,num,denom)
        epsi = num/denom
        gamma = 1/(e*e)*(1/(1-epsi))*(1-x)/(x*x*x)*a8p*q2/(0.938*.938)

    return [gamma, epsi]

if g:
    # df6 = pd.read_csv("xs_clas6.csv")

    # df6 = df6.query("q<4 and q>3.5 and x>0.4 and x<0.5 and t>0.6 and t<1")
    # df6 = df6[["p","dsdtdp","stat","sys"]]
    # print(df6)
    #q2bins,xBbins, tbins, phibins = [fs.q2bins_test, fs.xBbins_test, fs.tbins_test, fs.phibins_test]
    q2bins,xBbins, tbins, phibins = [fs.q2bins, fs.xBbins, fs.tbins, fs.phibins]

    num_counts = []

    df60 = pd.read_pickle('xs_clas6_binned.pkl')
    sample0 = pd.read_pickle("/mnt/d/GLOBUS/CLAS12/simulations/production/F2018_In_Norad/rec_and_gen_binned_events_meta.pkl")
    df00 = pd.read_pickle("F18_In_binned_events.pkl")

    for qmin,qmax in zip(q2bins[0:-1],q2bins[1:]):
        for xmin,xmax in zip(xBbins[0:-1],xBbins[1:]):
            for tmin,tmax in zip(tbins[0:-1],tbins[1:]):
                qave = round((qmin+qmax)/2,2)
                xave = round((xmin+xmax)/2,2)
                tave = round((tmin+tmax)/2,2)
                
                #print("On bins: q-{} x-{} t-{}".format(qave,xave,tave))

                

                bin_query = "qmin=={} and xmin=={} and tmin=={}".format(qmin,xmin,tmin)
                print("HERE IS THE QUERY")
                print(bin_query)
                df6 = df60.query(bin_query)

                if df6.shape[0]==0:
                    print("CLAS6 is EMPTY")

                else:

                    df = df00.query(bin_query)

                    sample = sample0.query(bin_query)

                    sample = sample[["qmin","xmin","tmin","pmin","rec_sum","gen_sum"]]

                    try:


                        dfout = pd.merge(df,sample,how='inner', on=['qmin','xmin','tmin','pmin'])

                        dfout.loc[:,"acc_corr"] = dfout["rec_sum"]/dfout["gen_sum"]
                        dfout.loc[:,"N_corr"] = dfout["counts"]/dfout["acc_corr"]

                        dfout.loc[:,"Lumi"] = 5.5E40
                        dfout.loc[:,"qmax"] = dfout["qmin"]+0.5
                        dfout.loc[:,"xmax"] = dfout["xmin"]+0.05
                        dfout.loc[:,"tmax"] = dfout["tmin"]+0.2
                        dfout.loc[:,"pmax"] = dfout["pmin"]+18
                        dfout.loc[:,"binvol"] = (dfout["qmax"]-dfout["qmin"])*(dfout["xmax"]-dfout["xmin"])*(dfout["tmax"]-dfout["tmin"])*(dfout["pmax"]-dfout["pmin"])*3.14159/180
                        dfout.loc[:,"xsec"] = dfout["N_corr"]/dfout["Lumi"]/dfout["binvol"]
                        dfout.loc[:,"gamma"] = get_gamma((dfout["xmin"]+dfout["xmax"])/2,(dfout["qmin"]+dfout["qmax"])/2)[0]
                        dfout.loc[:,"epsi"] = get_gamma((dfout["xmin"]+dfout["xmax"])/2,(dfout["qmin"]+dfout["qmax"])/2)[1]
                        dfout.loc[:,"red_xsec"] = dfout["xsec"]/dfout["gamma"]
                        dfout.loc[:,"red_xsecnb"] = dfout["red_xsec"]*1E33
                        dfout.loc[:,"red_xsecnb_err"] = dfout["red_xsecnb"]*np.sqrt(1/dfout["counts"]+1/dfout["rec_sum"] + 1/dfout["gen_sum"])
                        dfout.loc[:,"frac_err"] = dfout["red_xsecnb_err"]/dfout["red_xsecnb"]
                        dfout.loc[:,"p"] = dfout["pmin"]+9
                        dfout = dfout[["p","acc_corr","red_xsecnb","red_xsecnb_err"]]


                        dffinal = pd.merge(dfout,df6,how='inner', on='p')
                        dffinal = dffinal[dffinal['acc_corr'] > 0.005]


                        fig, ax = plt.subplots(figsize =(16, 10)) 

                        x = dffinal["p"]
                        x2min = dffinal["p"]-9

                        y1 = dffinal["red_xsecnb"]
                        y_err1 =   dffinal["red_xsecnb_err"]
                        y2 = dffinal["dsdtdp"]
                        y_err2 =   dffinal["stat"]

                        ratio = y1/y2
                        ratio_uncert = np.sqrt((y_err1/y1)**2+(y_err2/y2)**2)*ratio

                        #print("RATIO IS")
                        #print(ratio)

                        for p,r,r_uncert in zip(x2min,ratio,ratio_uncert):
                            num_counts.append([qmin,xmin,tmin,p,r,r_uncert])


                        plt.rcParams["font.size"] = "20"


                        ax.bar(x, y2,  yerr=y_err2, align='center', alpha=0.5, color="black", ecolor="black", capsize=10, width=16, label="CLAS6")
                        ax.bar(x, y1, yerr=y_err1, align='center', alpha=0.5, color="red", ecolor='red', capsize=10, width=16, label="CLAS12")

                        

                        #ax.set_ylim(0,100)
                        #ax.set_xticks(p_pos[::2])
                        #ax.set_xticklabels(p_labels[::2])
                        #ax.set_title('Simualtions: Gen and Rec Counts, $Q^2$ bin: {} GeV$^2$, $x_B$ bin: {}, t bin: {} GeV$^2$'.format(q_bin, x_bin, t_bin))
                        #ax.yaxis.grid(True)

                        #if float(q_bin)<10:
                        #    q_label = "0"+str(q_bin)

                        #print(q_label)


                        plt_title = "Reduced Cross Section, Q2= {}, xB = {}, -t = {}".format(qave, xave, tave)
                        ax.set_ylabel('Reduced Cross Section (nb/GeV$^2$)')
                        plt.rcParams["font.size"] = "20"
                        ax.yaxis.label.set_size(20) 
                        ax.set_title(plt_title)

                        ax.set_xlabel('Lepton-Hadron Plane Angle')
                        plt.rcParams["font.size"] = "20"
                        ax.xaxis.label.set_size(20) 

                        ax.legend(loc="best")
                        #plt_title = base_plot_dir+str(t_bin)+"/x_{}_q_{}_acc_inv.png".format(str(x_bin),q_label)
                        #plt.show()
                        plt.savefig("ratio_plots/"+plt_title)

                    except:
                        print("PLOTTING FAILED")


            # ax.bar(x, y, align='center', alpha=1, color="black", capsize=10, width=16)
            # ax.bar(x, yg, align='center', alpha=0.5, color="red", capsize=10, width=16)

            # ax.set_ylabel('Simualtions: Gen and Rec Counts')
            # #ax.set_ylim(0,100)
            # #ax.set_xticks(p_pos[::2])
            # #ax.set_xticklabels(p_labels[::2])
            # ax.set_title('Simualtions: Gen and Rec Counts, $Q^2$ bin: {} GeV$^2$, $x_B$ bin: {}, t bin: {} GeV$^2$'.format(q_bin, x_bin, t_bin))
            # #ax.yaxis.grid(True)

            # if float(q_bin)<10:
            #     q_label = "0"+str(q_bin)

            # print(q_label)

            # plt_title = base_plot_dir+str(t_bin)+"/x_{}_q_{}_acc_inv.png".format(str(x_bin),q_label)
            # plt.savefig(plt_title)

    print("making DF")
    df_minibin = pd.DataFrame(num_counts, columns = ['qmin','xmin','tmin','pmin','ratio','ratio_uncert'])
    print(df_minibin)
    df_minibin.to_pickle("ratioed.pkl")




    





#     pd.concat(
#     objs,
#     axis=0,
#     join="outer",
#     ignore_index=False,
#     keys=None,
#     levels=None,
#     names=None,
#     verify_integrity=False,
#     copy=True,
# # )
#     print(dfout)
# sys.exit()


# #outname = recon_file.split(".")[0]
# #output_loc_event_pkl_after_cuts = dirname+run+"/binned_pickles/"+outname+"_reconstructed_events_after_cuts.pkl"
# df = pd.read_pickle(output_loc_event_pkl_after_cuts)
# #df = df.query("Q2 > 2 and Q2 < 2.5 and xB < 0.38 and xB>0.3 and t>0.2 and t<0.3")

# # print(df.shape)

# # x_data = df["phi1"]
# # plot_title = "F 2018 Inbending, epgg, all exclusivity cuts"

# # #plot_title = "F 2018 Inbending, epgg, no exclusivity cuts"

# # vars = ["XB (GeV)"]
# # make_histos.plot_1dhist(x_data,vars,ranges="none",second_x="none",logger=False,first_label="F18IN",second_label="norad",
# #             saveplot=False,pics_dir="none",plot_title=plot_title,first_color="blue",sci_on=False)

# # sys.exit()

# df_gen = pd.read_pickle(output_loc_event_pkl_all_gen_events)
# #df = pd.read_pickle(save_base_dir+"100_20211103_1524_merged_Fall_2018_Inbending_gen_all_generated_events_all_generated_events.pkl")
# for col in df.columns:
#     print(col)

# df['t1'] = df['t']
# orginial_sum = df.shape[0]



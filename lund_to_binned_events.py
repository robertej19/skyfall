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

#This project

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def vec_angle(v1, v2):
    """stolen from https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python/13849249#13849249"""
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))*180/np.pi

def vec_subtract(vec1,vec2):
    res = tuple(map(lambda i, j: i - j, vec1, vec2)) 
    return res

def vec_add(vec1,vec2):
    res = tuple(map(lambda i, j: i + j, vec1, vec2)) 
    return res

def calc_inv_mass_squared(four_vector):
    fv = four_vector
    inv_mass2 = fv[0]**2-fv[1]**2-fv[2]**2-fv[3]**2
    return inv_mass2

def calculate_kinematics(event_df):
    ele = event_df.query("particleID == 11")
    pro = event_df.query("particleID == 2212")
    ic(ele)
    ic(pro)
    photons = event_df.query("particleID == 22")
    photon1 = photons.head(n=1)#This will only work for two photons!
    photon2 = photons.tail(n=1)#This will only work for two photons!
    ic(photons)
    ic(photon1)
    

    e_mass = 0.000511
    pro_mass = 0.938
    Ebeam_4mom = (10.6,0,0,10.6)
    e_4mom = (ele["E_GeV"].values[0],ele["mom_x"].values[0],ele["mom_y"].values[0],ele["mom_z"].values[0])
    #e_energy = np.sqrt(e_mass**2+np.sum(np.square(e_mom)))
    pro_4mom = (pro["E_GeV"].values[0],pro["mom_x"].values[0],pro["mom_y"].values[0],pro["mom_z"].values[0])
    target_4mom = (pro_mass,0,0,0)
    pho1_4mom = (photon1["E_GeV"].values[0],photon1["mom_x"].values[0],photon1["mom_y"].values[0],photon1["mom_z"].values[0])
    pho2_4mom = (photon2["E_GeV"].values[0],photon2["mom_x"].values[0],photon2["mom_y"].values[0],photon2["mom_z"].values[0])



    ic(e_4mom)
    ic(Ebeam_4mom)
    ic(pro_4mom)
    Eprime = e_4mom[0]
    ic(Eprime)

    virtual_gamma = vec_subtract(Ebeam_4mom,e_4mom)
    ic(virtual_gamma)
    
    #Calculate kinematic quantities of interest
    Q2 = -1*calc_inv_mass_squared(virtual_gamma)
    nu = virtual_gamma[0]
    xB = Q2/(2*pro_mass*nu)
    W2 = calc_inv_mass_squared(vec_subtract(vec_add(Ebeam_4mom,(pro_mass,0,0,0)),e_4mom))
    #print("W2 is: {}".format(W2))

    pion_4mom = vec_add(pho1_4mom,pho2_4mom)
    ic(pion_4mom)

    #Calculate t
    #t = (P - P')^2
    #t[i] = 2*0.938*(p4_proton[i].E() - 0.938);
    t = -1*calc_inv_mass_squared(vec_subtract(target_4mom,pro_4mom)) 
    #Could also calculate this using (target+e_beam-e'-pion)

    #Calculate phi (trento angle)

    e3 = e_4mom[1:]
    ic(e3)
    v_lepton = np.cross(Ebeam_4mom[1:],e_4mom[1:])
    v_hadron = np.cross(pro_4mom[1:],virtual_gamma[1:])
    v_hadron2 = np.cross(pro_4mom[1:],pion_4mom[1:])

    phi = vec_angle(v_lepton,v_hadron)

    if (np.dot(v_lepton,pro_4mom[1:])>0):
        phi = 360 - phi
    
    
    ic.disable()


    return Q2, xB, W2,nu,t,phi,Eprime, e_4mom, pro_4mom, pho1_4mom, pho2_4mom



lund_header_labels =["num_particles",
"target_mass",
"target_atom_num",
"target_pol",
"beam_pol",
"beam_type",
"beam_energy",
"interaction_nuc_id",
"process_id",
"event_weight"]

lund_particle_labels = ["sub_index",
    "lifetime",
    "type_1active",
    "particleID",
    "ind_parent",
    "ind_daughter",
    "mom_x",
    "mom_y",
    "mom_z",
    "E_GeV",
    "Mass_GeV",
    "Vx",
    "Vy",
    "Vz"]


def convert_lund_file_to_df(filename):
    print("Converting file {}".format(filename))
    events = []
    event_ind = -1
    with open(filename,"r") as f:
        for line in f:
            #ic(line)
            line_str = str(line)
            #ic(line_str[1])
            if line_str[1] != ' ': #LUND format has the number of particles in the second character of a header
                event_ind += 1
                #print("header")
                values = [event_ind,]
                events.append([])
                
                cols = line.split()  
                for ind, val in enumerate(cols):
                    values.append(float(val))
                #print(values)
                events[event_ind].append(values)
            
                #print(events)
                ###Write to header
            else:
                #ic(events)
                #ic(event_ind)
                values = []
                #print("particle content")
                cols = line.split()
                for ind, val in enumerate(cols):
                    values.append(float(val))
                events[event_ind].append(values)
    #ic.enable()
    events_repacked = []
    for event in events:
        for particle_ind in range(1,len(event)):
            #ic(event[particle_ind])   
            events_repacked.append(event[0]+event[particle_ind])    

    #ic(events_repacked)

    df_labels = ["event_num"]+lund_header_labels+lund_particle_labels
    df = pd.DataFrame(events_repacked, columns=df_labels)
    
    return df



    
def process_lund_into_events(df,run_num):
    events_list = []
    num_events = df["event_num"].max()
    ic(num_events)
    for ind in range(0,num_events+1):
        if ind % 100 ==0:
            ic.enable()
            ic(ind)
            ic.disable()
        event_dataframe = df.query("event_num == {}".format(ind))
        
        event_num = ind
        lumi = 0
        heli = 0
        Ebeam = 10.6
        
        q2,xb,W2,nu,t,phi,Eprime, e_4mom, pro_4mom, pho1_4mom, pho2_4mom = calculate_kinematics(event_dataframe)


        events_list.append([run_num,event_num,lumi,heli,
            Ebeam,Eprime,q2,xb,W2,nu,t,phi,
            e_4mom[0],e_4mom[1],e_4mom[2],e_4mom[3], 
            pro_4mom[0], pro_4mom[1], pro_4mom[2], pro_4mom[3], 
            pho1_4mom[0], pho1_4mom[1], pho1_4mom[2], pho1_4mom[3], 
            pho2_4mom[0],pho2_4mom[1],pho2_4mom[2],pho2_4mom[3]])
    
    return events_list


def get_events_from_lunds(data_dir,out_dir):
    data_list = os.listdir(data_dir)

    out_labels = ["run","event","luminosity","helicity","Ebeam","Eprime","Q2","xB","W2","nu","t1","phi1",
                    "GenEe","GenEpx","GenEpy","GenEpz",
                    "GenPe","GenPpx","GenPpy","GenPpz",
                    "GenGe","GenGpx","GenGpy","GenGpz",
                    "GenG2e","GenG2px","GenG2py","GenG2pz"]
    for count,lund_pickle in enumerate(data_list):
        print("on file {} of {}".format(count,len(data_list)))

        #run_num = str(lund_pickle).split(".dat")[0].split("_")[-1]
        run_num = 1
        filename = data_dir+lund_pickle
        print(filename)
        df = pd.read_pickle(filename)


        events_list = process_lund_into_events(df,run_num)
        df_out = pd.DataFrame(events_list, columns=out_labels)

        print(df_out.columns)

        # make_histos.plot_2dhist(x_data,y_data,var_names,ranges,colorbar=True,
        #     saveplot=False,pics_dir="none",plot_title="none",
        #     filename="ExamplePlot",units=["",""])

        df_out.to_pickle(out_dir+lund_pickle.split(".")[0]+"_events.pkl")



lund_dir = "EnergyDep/lunds/"
pkl_dir = "EnergyDep/raw_pkls/"
event_dir = "EnergyDep/evented_pkls/"
binned_dir = "EnergyDep/binned_pkls/"

high_E = "CLAS12_test"
low_E = "CLAS6_test"

files = [high_E,low_E]


convert_lund_to_pkl = True
convert_pkl_to_events = True
convert_events_to_bins = True
get_ratio = True

fs = filestruct.fs()


if convert_lund_to_pkl:
    for lund_file in files:
        df = convert_lund_file_to_df(lund_dir+lund_file+".lund")
        print(df)
        df.to_pickle(pkl_dir+lund_file+".pkl")

if convert_pkl_to_events:
    #for event_file in files:
    #    df = pd.read_pickle(event_dir+lund_file+"_evented.pkl")
    #data_dir = "angle_studies/20_pandas/"
    #out_dir = "angle_studies/20_evented/"

    get_events_from_lunds(pkl_dir,event_dir)

if convert_events_to_bins:

    for filename in files:
        df0 = pd.read_pickle(event_dir+filename+"_events.pkl")

        q2bins,xBbins, tbins, phibins = fs.q2bins, fs.xBbins, fs.tbins, fs.phibins
        #if True:
        #        q2bins,xBbins, tbins, phibins = fs.q2bins_test, fs.xBbins_test, fs.tbins_test, fs.phibins_test

        num_counts = []

        qrange = [q2bins[0], q2bins[-1]]
        xBrange = [xBbins[0], xBbins[-1]]
        trange = [tbins[0], tbins[-1]]

        total_num = df0.query('Q2>{} and Q2<{} and xB>{} and xB<{} and t1>{} and t1<{}'.format(*qrange, *xBrange, *trange)).shape[0]

        prefix = ""

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
        df_minibin.to_pickle(binned_dir+filename+"_binned.pkl")
        print(df_minibin)

if get_ratio:
    df_low = pd.read_pickle(binned_dir+low_E+"_binned.pkl")
    df_low = df_low.rename(columns={"counts":"counts_low"})
    #                    dfout.loc[:,"qmax"] = dfout["qmin"]+0.5

    #df_low.loc[:,"counts_low"] = df_low["counts"]
    df_high = pd.read_pickle(binned_dir+high_E+"_binned.pkl")
    df_high = df_high.rename(columns={"counts":"counts_high"})

    
    df_out = pd.merge(df_high,df_low,how='inner', on=['qmin','xmin','tmin','pmin'])
    df_out.loc[:,"Energy_ratio"] = df_out["counts_high"]/df_out["counts_low"]

    print(df_out)
    df_out.to_pickle("EnergyDependenceRatio.pkl")







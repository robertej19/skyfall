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

#This project


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
            ic(line_str[1])
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
                ic(events)
                ic(event_ind)
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
            ic(event[particle_ind])   
            events_repacked.append(event[0]+event[particle_ind])    

    ic(events_repacked)

    df_labels = ["event_num"]+lund_header_labels+lund_particle_labels
    df = pd.DataFrame(events_repacked, columns=df_labels)
    
    return df

def convert_lund_dir_to_dfs(data_dir,out_dir):
    data_list = os.listdir(data_dir)

    for lund_file in data_list:
        ic.disable()
        df = convert_lund_file_to_df(data_dir+lund_file)
        print("DF IS ----------")
        print(df)  
        df.to_pickle(out_dir+lund_file+".pkl")
    
    print("\nProcessed {} files from {}\n".format(len(data_list),data_dir))
    print("Saved pkl files to {}\n".format(out_dir))

if __name__ == "__main__":
    #data_dir = "/mnt/d/GLOBUS/CLAS12/simulations/production/Fall_2018_Inbending/Test/lunds/"
    #out_dir = "/mnt/d/GLOBUS/CLAS12/simulations/production/Fall_2018_Inbending/Test/panda_lunds/"
  
    #data_dir = "ex_rad_lund/"
    #data_dir = "ex_norad_lund/"
    #out_dir = "ex_rad_pd/"
    #out_dir = "ex_norad_pd/"
    data_dir = "angle_studies/2_lund/"
    out_dir = "angle_studies/2_evented/"

    convert_lund_dir_to_dfs(data_dir,out_dir)


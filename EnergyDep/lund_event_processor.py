import matplotlib
matplotlib.use('Agg')
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
#import utils.make_histos
import make_histos

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


def plotter(filename):
    df_out = pd.read_pickle(filename)

    print(df_out.columns)

    x_data = df_out['xb']
    y_data = df_out['q2']
    var_names = ['xb','q2']
    ranges = [[0,1,100],[0,12,120]]

    make_histos.plot_2dhist(x_data,y_data,var_names,ranges,colorbar=True,
            saveplot=False,pics_dir="none",plot_title="none",
            filename="ExamplePlot",units=["",""])


def get_events_from_lunds(data_dir,out_dir):
    data_list = os.listdir(data_dir)

    out_labels = ["run","event","luminosity","helicity","Ebeam","Eprime","q2","xb","W2","nu","t","phi",
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

        x_data = df_out['xb']
        y_data = df_out['q2']
        var_names = ['xb','q2']
        ranges = [[0,1,10],[0,12,12]]

        # make_histos.plot_2dhist(x_data,y_data,var_names,ranges,colorbar=True,
        #     saveplot=False,pics_dir="none",plot_title="none",
        #     filename="ExamplePlot",units=["",""])

        df_out.to_pickle(out_dir+lund_pickle+"_events.pkl")


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

    return df



if __name__ == "__main__":


    # data_dir = "lunds_pickle/"
    # out_dir = "./"


    #data_dir = "/mnt/d/GLOBUS/CLAS12/simulations/production/Fall_2018_Inbending/Test/lunds/"
    #data_dir = "/mnt/d/GLOBUS/CLAS12/simulations/production/Fall_2018_Inbending/Test/filts/"
    #out_dir = "/mnt/d/GLOBUS/CLAS12/simulations/production/Fall_2018_Inbending/Test/filts/"
    #df = convert_lund_dir_to_dfs(data_dir,out_dir)
    
    # data_dir = "ex_rad_pd/EXAMPLERAD4000.lund.pkl"
    # data_dir2 = "ex_norad_pd/EXAMPLENORAD5000.lund.pkl"

    # #mass_statement = "Mass_GeV < 0.00052 and Mass_GeV > 0.000499"
    # mass_statement = "Mass_GeV < 0.95 and Mass_GeV > 0.9"
    # df = pd.read_pickle(data_dir)
    # ic(df.head(12))
    # df = df.query(mass_statement)

    # #ic(df.columns)
    # #ic(df.head(4))

    # vars = ["Proton Energy Distribution"]
    # x_data = df['E_GeV']


    # df2 = pd.read_pickle(data_dir2)
    # ic(df2.head(12))
    # df2 = df2.query(mass_statement)
    # df2 = df2.sample(df.shape[0])
    # ic(df2.columns)
    # ic(df2.head(4))
    # x_data_2 = df2['E_GeV']


    # make_histos.plot_1dhist(x_data,vars,ranges="none",second_x=x_data_2,
    #             saveplot=True,pics_dir="./",plot_title=vars[0],first_color="blue",sci_on=False)
    

    # sys.exit()
    #data_dir = "ex_norad_pd/"
    #out_dir = "evented_norad/"
    data_dir = "angle_studies/20_pandas/"
    out_dir = "angle_studies/20_evented/"

    get_events_from_lunds(data_dir,out_dir)
    #filename = "TEST3_real_output_aao_norad.lund_events.pkl"
    #plotter(filename)




        
        
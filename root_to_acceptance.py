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

def get_files(dirname):
    """
    Get all files in base_dir
    """

    print("looking for files in {}".format(dirname))
    jobs_list = []
    print(os.listdir(dirname))
    for f in os.listdir(dirname):
        f_path = dirname + f
        if os.path.isfile(f_path):
            jobs_list.append(f)
    jobs_list = sorted(jobs_list)
    print("Found {} files in jobs directory".format(len(jobs_list)))

    for fname in jobs_list:
        if "recon" in fname:
            recon_file = fname
        if "gen" in fname:
            gen_file = fname

    print("Generator file: {}".format(gen_file))
    print("Reconstructed file: {}".format(recon_file))
            
    return gen_file, recon_file

fs = filestruct.fs()


allz = True
QuickTesting = False
#All
if allz:
    DoWipe = True
    DoGen = True
    DoRecon = True
    DoInspect = True
#Just ana
else:
    DoWipe = False
    DoGen = False
    DoRecon = True
    DoInspect = False
    DoBin = True

base_dir = "/mnt/d/GLOBUS/CLAS12/simulations/production/F2018_In_Norad"
# To add data to this base dir, make a subdir /runs/<run_number>/roots and put the files there
# Note that dataset 101 is just for testing
if base_dir[-1] != '/':
    base_dir += '/'

dirname = base_dir+"runs/"
#print(os.listdir(dirname))
runs_list = os.listdir(dirname)


dir_ending_list = ["/event_pickles/", "/plots/", "/binned_pickles/"]
for run in runs_list:
    for dir_ending in dir_ending_list:
        if DoWipe:
            shutil.rmtree(dirname+run+dir_ending) if os.path.isdir(dirname+run+dir_ending) else None
        if not os.path.isdir(dirname+run+dir_ending):
            os.makedirs(dirname+run+dir_ending)

generator_type = "norad"

converter_gen = convert_GEN_NORAD_root_to_pkl
converter_recon = convert_RECON_NORAD_root_to_pkl

if generator_type == "rad":
    converter_gen = convert_GEN_RAD_root_to_pkl
    converter_recon = convert_RECON_RAD_root_to_pkl

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-f","--fname", help="a single root file to convert into pickles", default="/Users/sangbaek/Dropbox (MIT)/data/project/merged_9628_files.root")
    parser.add_argument("-o","--out", help="a single pickle file name as an output", default="goodbyeRoot.pkl")
    parser.add_argument("-s","--entry_stop", help="entry_stop to stop reading the root file", default = None)
    parser.add_argument("-p","--polarity", help="polarity", default = "inbending")
    
    args = parser.parse_args()

    if QuickTesting:
        runs_list = [runs_list[0]]
        print("Using first run directory only as test: {}".format(runs_list))
    for run in runs_list:
        root_file_list = os.listdir(dirname+run+"/roots/")
        gen_file, recon_file = get_files(dirname+run+"/roots/")

        if DoGen:
            outname = gen_file.split(".")[0]+"_all_generated_events"
            output_loc_event_pkl = dirname+run+"/event_pickles/"+outname+"_all_generated_events.pkl"
            output_loc_plots = dirname+run+"/plots/"+outname+"/"

            tree = converter_gen.readFile(dirname+run+"/roots/" + gen_file)
            df_gen_all = converter_gen.readEPGG(tree)
            #ic("saving file to: {}".format(output_loc_event_pkl))
            df_gen_all.to_pickle(output_loc_event_pkl)
        
            histo_plotting.make_all_histos(df_gen_all,datatype="Gen",
                hists_2d=True,hists_1d=True,hists_overlap=False,
                saveplots=True,output_dir=output_loc_plots)

            print(df_gen_all.columns)
            print(df_gen_all.head(5))
            print("Number of events: {}".format(df_gen_all.shape[0]))

        if DoRecon:
            outname = recon_file.split(".")[0]
            args.fname = dirname+run+"/roots/" + recon_file
            output_loc_event_pkl_after_cuts = dirname+run+"/event_pickles/"+outname+"_reconstructed_events_after_cuts.pkl"
            output_loc_plots_after_cuts = dirname+run+"/plots/"+outname+"_after_cuts/"
            output_loc_plots_truth_after_cuts = dirname+run+"/plots/"+outname+"_truth_after_cuts/"

            output_loc_event_pkl_before_cuts = dirname+run+"/event_pickles/"+outname+"_reconstructed_events_before_cuts.pkl"
            output_loc_plots_before_cuts = dirname+run+"/plots/"+outname+"_before_cuts/"
            output_loc_plots_truth_before_cuts = dirname+run+"/plots/"+outname+"_truth_before_cuts/"



            converter = root2pickle(args.fname, entry_stop = args.entry_stop, pol = args.polarity)
            df_after_cuts = converter.df_after_cuts
            df_after_cuts.to_pickle(output_loc_event_pkl_after_cuts)

            df_before_cuts = converter.df_before_cuts
            df_before_cuts.to_pickle(output_loc_event_pkl_before_cuts)

            print("NOW PLOTTING RECON AFTER CUTS")
            histo_plotting.make_all_histos(df_after_cuts,datatype="Recon",
                                hists_2d=True,hists_1d=True,hists_overlap=False,
                                saveplots=True,output_dir=output_loc_plots_after_cuts)

            print("NOW PLOTTING RECON TRUTH AFTER CUTS")

            histo_plotting.make_all_histos(df_after_cuts,datatype="Truth",
                                hists_2d=True,hists_1d=False,hists_overlap=False,
                                saveplots=True,output_dir=output_loc_plots_truth_after_cuts)

            print("NOW PLOTTING RECON BEFORE CUTS")

            histo_plotting.make_all_histos(df_before_cuts,datatype="Recon",
                                hists_2d=True,hists_1d=True,hists_overlap=False,
                                saveplots=True,output_dir=output_loc_plots_before_cuts)

            print("NOW PLOTTING RECON TRUTH BEFORE CUTS")

            histo_plotting.make_all_histos(df_before_cuts,datatype="Truth",
                                hists_2d=True,hists_1d=False,hists_overlap=False,
                                saveplots=True,output_dir=output_loc_plots_truth_before_cuts)

            #print(df_after_cuts.columns)
            #print(df_after_cuts.head(5))
            print("Number of events: {}".format(df_after_cuts.shape[0]))

    if DoBin:
        outname = recon_file.split(".")[0]
        #df = pd.read_pickle(save_base_dir+"100_20211103_1524_merged_Fall_2018_Inbending_gen_all_generated_events_all_generated_events.pkl")

        df = pd.read_pickle(save_base_dir+outname+"_reconstructed_events.pkl")

        four_squared = df[["Q2", "W", "xB", "t", "phi1"]]#.head(10)
        ic(four_squared)



        bins = [fs.q2bins, fs.xBbins, fs.tbins, fs.phibins]
        if args.test:
                bins = [fs.q2bins_test, fs.xBbins_test, fs.tbins_test, fs.phibins_test]

        qlabels, xBlabels, tlabels, philabels = [],[] ,[],[]

        labels = [qlabels, xBlabels, tlabels, philabels]

        for label_item,bin_item in zip(labels,bins):
            for count in range(1,len(qbins)):
                label_item.append(str(bin_item[count-1])+"-"+str(bin_item[count]))
        
        four_squared['qbin'] = pd.cut(four_squared['Q2'], qbins, labels=qlabels)
        four_squared['tbin'] = pd.cut(four_squared['t'], tbins, labels=tlabels)
        four_squared['xBbin'] = pd.cut(four_squared['xB'], xBbins, labels=xBlabels)
        four_squared['phibin'] = pd.cut(four_squared['phi1'], phibins, labels=philabels)

        rude_sum = 0

        num_counts = []

        for qval in qlabels:
            df_min = four_squared.query("qbin==@qval")
            if len(df_min.index) == 0:
                    num_counts.append([0]*len(xBlabels)*len(tlabels)*len(philabels))
                    print([0]*len(xBlabels)*len(tlabels)*len(philabels))
                    print("made a triple")
            else:
                for xval in xBlabels:
                    df_min2 = df_min.query("xBbin==@xval")
                    if len(df_min2.index) == 0:
                        num_counts.append([0]*len(tlabels)*len(philabels))
                        print([0]*len(tlabels)*len(philabels))
                        print("made a triple")

                    else:
                        for tval in tlabels:
                            df_min3 = df_min2.query("tbin==@tval")
                            if len(df_min3.index) == 0:
                                for i in philabels:
                                    num_counts.append(0)
                            else:
                                for phival in philabels:
                                    df_min4 = df_min3.query("phibin==@phival")
                                    print(len(df_min4.index))
                                    rude_sum += len(df_min4.index)
                                    num_counts.append(len(df_min4.index))

        print(num_counts)

        a1 = pd.DataFrame(tlabels,columns=["t"])
        b1 = pd.DataFrame(philabels,columns=["phi"])
        f1 = pd.merge(a1,b1,how='cross')

        a = pd.DataFrame(qlabels,columns=["Q"])
        b = pd.DataFrame(xBlabels,columns=["x"])
        f = pd.merge(a,b,how='cross')

        ff = pd.merge(f,f1,how='cross')
        ic(a)
        ic(b)
        #f = pd.DataFrame(qlabels, xBlabels, tlabels, philabels,num_counts, columns=['Month','Day','Year','Hour','date'])

        ff['counts'] = num_counts
        ic(ff)
        ic(ff.sum())

        ff.to_pickle(save_base_dir+outname+"_reconstructed_events_binned.pkl")

sys.exit()



if DoWipe:
    if not os.path.exists(save_base_dir):
        os.makedirs(save_base_dir)
    else:
        shutil.rmtree(save_base_dir)
        os.makedirs(save_base_dir)

sys.exit()

base_dir = "/mnt/d/GLOBUS/CLAS12/simulations/production/new_100/raw_root/"
#base_dir = "/mnt/d/GLOBUS/CLAS12/simulations/production/Fall_2018_Inbending/Raw_Root_Files/Norad/testana/"










#############!!!!!!!!!!!!!!!!!!
    

if True:
    """-----------------------------------------"""
    
    if base_dir[-1] != '/':
        base_dir += '/'

    save_base_dir =  base_dir + "analyzed/"



    gen_file, recon_file = get_files(base_dir)
    


    # Process gen_all


    # Process recon
    # outname = recon_file.split(".")[0]

    # tree = converter_gen.readFile(base_dir + recon_file)

    # df_gen, df_rec = converter_recon.readEPGG(tree)
    # ic("saving file to: {}".format(save_base_dir+outname+"_reconstructed_events.pkl"))
    # df_rec.to_pickle(save_base_dir+outname+"_reconstructed_events.pkl")
    # df_gen.to_pickle(save_base_dir+outname+"_detected_gen_events.pkl")

    # df_rec2 = pickle_analysis.makeDVpi0vars(df_rec)
    # ic(df_rec2.columns)
    # print("RYIN RECON PLOT")                    
    # histo_plotting.make_all_histos(df_rec2,datatype="Recon",
    #                     hists_2d=True,hists_1d=True,hists_overlap=False,
    #                     saveplots=True,output_dir=save_base_dir+outname+"/")

    # print(df_rec2.columns)
    # print(df_rec2.head(5))
    # print("Number of events: {}".format(df_rec2.shape[0]))



    if DoInspect:
        pkl_name = base_dir + "analyzed/100_20211103_1524_merged_Fall_2018_Inbending_recon_reconstructed_events.pkl"
        #pkl_name = base_dir + "analyzed/5000_20210731_2317_norad_recon_reconstructed_events.pkl"

        df = pd.read_pickle(pkl_name)

        mins = df.min()
        ic(df.GenEtheta.min())#for minz,name in mins:
        ic(df.Etheta.min())#for minz,name in mins:

        

        pkl_name = base_dir + "analyzed/100_20211103_1524_merged_Fall_2018_Inbending_gen_all_generated_events_all_generated_events.pkl"
        #pkl_name = base_dir + "analyzed/5000_20210731_2317_norad_gen_all_generated_events_all_generated_events.pkl"

        df = pd.read_pickle(pkl_name)

        mins = df.min()
        ic(df.GenEtheta.min())#for minz,name in mins:
        df0=df
        #df0 = df.query("GenEtheta<10 and GenW>2")        #    ic(minz,name)

        ic(df0.columns)

        x_data = df0["GenW"]
        plot_title = "F 2018 Inbending, epgg, all exclusivity cuts"

        #plot_title = "F 2018 Inbending, epgg, no exclusivity cuts"

        vars = ["XB (GeV)"]
        make_histos.plot_1dhist(x_data,vars,ranges="none",second_x="none",logger=False,first_label="F18IN",second_label="norad",
                    saveplot=False,pics_dir="none",plot_title=plot_title,first_color="blue",sci_on=False)
            

        y_data = df0["GenQ2"]
        var_names = ["Etheta","Q2"]
        ranges = [[0,10,100],[0,2,100]]
        make_histos.plot_2dhist(x_data,y_data,var_names,ranges,colorbar=True,
                    saveplot=False,pics_dir="none",plot_title=plot_title,logger=False,first_label="rad",
                    filename="ExamplePlot",units=["GeV","GeV^2"])


        
        #ic(four_squared)

        #for index, row in four_squared.iterrows():
        #    for q_ind, q_val in enumerat(qbins):
        #        if q_val>row["Q2"]:
                    
                
        #    print(row['Q2'], row['xB'])




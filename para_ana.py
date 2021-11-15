import uproot
import pandas as pd
import numpy as np
import argparse
import itertools
import os, sys
from icecream import ic
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


base_dir = "/mnt/d/GLOBUS/CLAS12/simulations/production/new_100/raw_root/"
#base_dir = "/mnt/d/GLOBUS/CLAS12/simulations/production/Fall_2018_Inbending/Raw_Root_Files/Norad/testana/"

allz = False
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
    DoRecon = False
    DoInspect = False
    DoBin = True


generator_type = "norad"

converter_gen = convert_GEN_NORAD_root_to_pkl
converter_recon = convert_RECON_NORAD_root_to_pkl

if generator_type == "rad":
    converter_gen = convert_GEN_RAD_root_to_pkl
    converter_recon = convert_RECON_RAD_root_to_pkl


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


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-f","--fname", help="a single root file to convert into pickles", default="/Users/sangbaek/Dropbox (MIT)/data/project/merged_9628_files.root")
    parser.add_argument("-o","--out", help="a single pickle file name as an output", default="goodbyeRoot.pkl")
    parser.add_argument("-s","--entry_stop", help="entry_stop to stop reading the root file", default = None)
    parser.add_argument("-p","--polarity", help="polarity", default = "inbending")
    
    args = parser.parse_args()

    
    """-----------------------------------------"""


    
    if base_dir[-1] != '/':
        base_dir += '/'

    save_base_dir =  base_dir + "analyzed/"

    if DoWipe:
        if not os.path.exists(save_base_dir):
            os.makedirs(save_base_dir)
        else:
            shutil.rmtree(save_base_dir)
            os.makedirs(save_base_dir)


    gen_file, recon_file = get_files(base_dir)
    


    # Process gen_all
    if DoGen:
        outname = gen_file.split(".")[0]+"_all_generated_events"

        tree = converter_gen.readFile(base_dir + gen_file)
        df_gen_all = converter_gen.readEPGG(tree)
        ic("saving file to: {}".format(save_base_dir+outname+".pkl"))
        df_gen_all.to_pickle(save_base_dir+outname+"_all_generated_events.pkl")
    
        histo_plotting.make_all_histos(df_gen_all,datatype="Gen",
            hists_2d=True,hists_1d=True,hists_overlap=False,
            saveplots=True,output_dir=save_base_dir+outname+"/")

        print(df_gen_all.columns)
        print(df_gen_all.head(5))
        print("Number of events: {}".format(df_gen_all.shape[0]))


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

    if DoRecon:
        outname = recon_file.split(".")[0]
        args.fname = base_dir + recon_file

        converter = root2pickle(args.fname, entry_stop = args.entry_stop, pol = args.polarity)
        df = converter.df
        df.to_pickle(save_base_dir+outname+"_reconstructed_events.pkl")

        histo_plotting.make_all_histos(df,datatype="Recon",
                            hists_2d=True,hists_1d=True,hists_overlap=False,
                            saveplots=True,output_dir=save_base_dir+outname+"/")

        print(df.columns)
        print(df.head(5))
        print("Number of events: {}".format(df.shape[0]))


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

    if DoBin:
        outname = recon_file.split(".")[0]
        #df = pd.read_pickle(save_base_dir+"100_20211103_1524_merged_Fall_2018_Inbending_gen_all_generated_events_all_generated_events.pkl")

        df = pd.read_pickle(save_base_dir+outname+"_reconstructed_events.pkl")

        four_squared = df[["Q2", "W", "xB", "t", "phi1"]]#.head(10)
        ic(four_squared)

        qbins = [0, 4,123]
        qlabels = ['0-1', '4-5']
        tbins = [0,.2,63]
        tlabels = ['0-1', '1-12']
        xBbins = [0, .6, 1.8]
        xBlabels = ['0-.6', '.6-1']
        phibins = [0, 120, 360]
        philabels = ['0-120', '120-360']

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

        #c = list(itertools.product(qlabels, xBlabels, tlabels, philabels))
        #print(c)


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
        #ic(four_squared)

        #for index, row in four_squared.iterrows():
        #    for q_ind, q_val in enumerat(qbins):
        #        if q_val>row["Q2"]:
                    
                
        #    print(row['Q2'], row['xB'])


















sys.exit()
if True:

    #For rad and norad
    # data_paths = ["Rad/Gen/","Rad/Recon/","Norad/Gen/","Norad/Recon/",]
    # converters = [convert_GEN_RAD_root_to_pkl,convert_RECON_RAD_root_to_pkl,
    #                convert_GEN_NORAD_root_to_pkl,convert_RECON_NORAD_root_to_pkl]
    # data_paths = ["Rad/Gen/","Rad/Recon/","Norad/Gen/","Norad/Recon/"]
    #For just rad
    #data_paths = ["Rad/Gen/","Rad/Recon/"]
    #For just norad
    # data_paths = ["Norad/Gen/","Norad/Recon/"]
    # converters = [convert_GEN_NORAD_root_to_pkl,convert_RECON_NORAD_root_to_pkl]

    data_paths = ["Gen/","Recon/"]


    for index,dpath in enumerate(data_paths):
        data_step = 0


        #print("On path {}".format(dpath))
        #dirname = '/mnt/d/GLOBUS/CLAS12/simulations/production/Fall_2018_Inbending/Raw_Root_Files/'+dpath    
        # dirname = '/mnt/d/GLOBUS/CLAS12/simulations/production/Fall_2018_Inbending/Before_Cuts/'+dpath
        # dirname = '/mnt/d/GLOBUS/CLAS12/simulations/production/Fall_2018_Inbending/After_Cuts/'+dpath

        #Inbending rad
        # dirname = '/mnt/d/GLOBUS/CLAS12/simulations/W_LessThan_1900MeV_rad/4_Final_Output_Files/Fall_2018_Inbending/Raw_Root_Files/'+dpath
        # dirname = '/mnt/d/GLOBUS/CLAS12/simulations/W_LessThan_1900MeV_rad/4_Final_Output_Files/Fall_2018_Inbending/Before_Cuts/'+dpath
        # dirname = '/mnt/d/GLOBUS/CLAS12/simulations/W_LessThan_1900MeV_rad/4_Final_Output_Files/Fall_2018_Inbending/After_Cuts/'+dpath

        #Inbending Norad
        # dirname = '/mnt/d/GLOBUS/CLAS12/simulations/W_LessThan_1900MeV_norad/4_Final_Output_Files/Fall_2018_Inbending/Raw_Root_Files/'+dpath
        # dirname = '/mnt/d/GLOBUS/CLAS12/simulations/W_LessThan_1900MeV_norad/4_Final_Output_Files/Fall_2018_Inbending/Before_Cuts/'+dpath
        # dirname = '/mnt/d/GLOBUS/CLAS12/simulations/W_LessThan_1900MeV_norad/4_Final_Output_Files/Fall_2018_Inbending/After_Cuts/'+dpath

        # #Outbending Norad
        # dirname = '/mnt/d/GLOBUS/CLAS12/simulations/W_LessThan_1900MeV_norad/4_Final_Output_Files/Fall_2018_Outbending_100/Raw_Root_Files/'+dpath
        # dirname = '/mnt/d/GLOBUS/CLAS12/simulations/W_LessThan_1900MeV_norad/4_Final_Output_Files/Fall_2018_Outbending_100/Before_Cuts/'+dpath
        # dirname = '/mnt/d/GLOBUS/CLAS12/simulations/W_LessThan_1900MeV_norad/4_Final_Output_Files/Fall_2018_Outbending_100/After_Cuts/'+dpath

        #Outbending rad
        #dirname = '/mnt/d/GLOBUS/CLAS12/simulations/W_LessThan_1900MeV_rad/4_Final_Output_Files/Fall_2018_Outbending_100/Raw_Root_Files/'+dpath
        #dirname = '/mnt/d/GLOBUS/CLAS12/simulations/W_LessThan_1900MeV_rad/4_Final_Output_Files/Fall_2018_Outbending_100/Before_Cuts/'+dpath
        #dirname = '/mnt/d/GLOBUS/CLAS12/simulations/W_LessThan_1900MeV_rad/4_Final_Output_Files/Fall_2018_Outbending_100/After_Cuts/'+dpath
        dirname = "/mnt/d/GLOBUS/CLAS12/simulations/production/new_100/"+dpath

        # First we convert from Root to Pkl

        if data_step == 0:

            jobs_list = []
            for f in os.listdir(dirname):
                f_path = dirname + f
                if os.path.isfile(f_path):
                    jobs_list.append(f)
            jobs_list = sorted(jobs_list)
            print("Found {} files in jobs directory".format(len(jobs_list)))

            for fname in jobs_list:
                print("Processing file {}".format(fname))
                tree = converters[index].readFile(dirname+fname)

                outname = fname.split(".")[0]
                if "Recon" in dpath:
                    df_gen, df_rec = converters[index].readEPGG(tree)
                    ic("saving file to: {}".format(outname))
                    df_rec.to_pickle(dirname+outname+"_reconstructed_events.pkl")
                    df_gen.to_pickle(dirname+"../Gen/"+outname+"_recon_generated_events.pkl")

                    df_rec2 = pickle_analysis.makeDVpi0vars(df_rec)
                    ic(df_rec2.columns)
                    print("RYIN RECON PLOT")                    
                    histo_plotting.make_all_histos(df_rec2,datatype="Recon",
                        hists_2d=True,hists_1d=True,hists_overlap=False,
                        saveplots=True,output_dir=dirname+outname+fname.split(".")[0]+"/")

                elif "Gen" in dpath:
                    df_gen = converters[index].readEPGG(tree)
                    ic("saving file to: {}".format(dirname+outname+"_all_generated_events.pkl"))
                    df_gen.to_pickle(dirname+outname+"_all_generated_events.pkl")

                    ic(df_gen.columns)
                    print("RYIN GEN PLOT")   
                    histo_plotting.make_all_histos(df_gen,datatype="Gen",
                        hists_2d=True,hists_1d=True,hists_overlap=False,
                        saveplots=True,output_dir=dirname+outname+fname.split(".")[0]+"/")

                    print(df_gen.shape)
                    print(df_gen.columns)
                    print(df_gen.head(5))
            data_step += 1
        
        # Now we implement DVEP cuts
        #print("Now implementing DVEP cuts")
        if data_step == 1:
            jobs_list = []
            for f in os.listdir(dirname):
                f_path = dirname + f
                #print(f_path)
                if os.path.isfile(f_path):
                    jobs_list.append(f)
                    print(f)
            jobs_list = sorted(jobs_list)

            ic(jobs_list)

            size_gen_chunks = 4000000 #Change this, depending on how big the gen dataset is


            for fname in jobs_list:
                dfs = []

                ic(fname)
                datatype = 'Gen' if 'gen' in fname else 'Recon'
                df = pd.read_pickle(dirname+fname)
                n = size_gen_chunks  #chunk row size
                ic(df.shape)
                list_df = []
                for i in range(0,df.shape[0],n):
                    ic(i)
                    #ic(i+n)
                    list_df.append(df[i:i+n])

                for index, df_chunk in enumerate(list_df):
                    #print("On DF chunk {}".format(index))
                    ic(df_chunk.shape)

                    outname = fname.split(".")[0]+"_{}.pkl".format(index)
                    print(outname)
                    if datatype == "Gen":
                        df_gen = pickle_analysis.makeGenDVpi0vars(df_chunk)
                        ic(df_gen.shape)

                        #df_gen.to_pickle(dirname+outname)
                        dfs.append(df_gen)
                    elif datatype == "Recon":
                        df_recon_pi0vars = pickle_analysis.makeDVpi0vars(df_chunk)
                        df_recon = pickle_analysis.cutDVpi(df_recon_pi0vars)
                        df_recon.to_pickle(dirname+fname)
                        ic(df_recon)
                        dfs.append(df_recon)


                df = dfs[0]

                histo_plotting.make_all_histos(df,datatype=datatype,
                    hists_2d=True,hists_1d=True,hists_overlap=False,
                    saveplots=True,output_dir=dirname+fname.split(".")[0]+"/")
            # histo_plotting.make_all_histos(df,datatype=datatype,
            #     hists_2d=True,hists_1d=True,hists_overlap=False,
            #     saveplots=True,output_dir=dirname+fname.split(".")[0]+"/hist_1D/")
            #print("Found {} files in jobs directory".format(len(jobs_list)))



    #File structure:
    # /mnt/d/GLOBUS/CLAS12/simulations/producion
    # ├── Fall_2018_Inbending
    # │   └── Raw_Root_Files
    # │   │   └── Rad
    # │   │   │   └── Gen
    # │   │   │   │   └── merged_gen_5000.root
    # │   │   │   │   └── merged_gen_6000.root
    # │   │   │   └── Recon
    # │   │   │   │   └── merged_recon_5000.root
    # │   │   │   │   └── merged_recon_6000.root
    # │   │   └── Norad
    # │   │   │   │   └── merged_gen_4000.root
    # │   │   │   │   └── merged_gen_5600.root
    # │   │   │   └── Recon
    # │   │   │   │   └── merged_recon_4000.root
    # │   │   │   │   └── merged_recon_5600.root
    # │   └── converted files before cuts
                # file
                # corresponding histogram pdf
    # │   └── after cuts
                # file
                # corresponding histogram pdf
    # │   └── binned data
                # file
                # corresponding histogram pdf


    # The following is needed since an executable does not have __file__ defined, but when working in interpreted mode,
    # __file__ is needed to specify the relative file path of other packages. In principle strict relative 
    # path usage should be sufficient, but it is easier to debug / more robust if absolute.
    try:
        __file__
    except NameError:
        full_file_path = sys.executable #This sets the path for compiled python
    else:
        full_file_path = os.path.abspath(__file__) #This sets the path for interpreted python


#     main_source_dir = "/".join(full_file_path.split("/")[:-3])

#     now = datetime.now()
#     dt_string = now.strftime("%Y%m%d_%H%M")
#     subdirs = ["0_JSub_Factory","1_Generated_Events",
#             "2_GEMC_DSTs","3_Filtered_Converted_Root_Files","4_Final_Output_Files"]


#         # Jsub file creator for norad and rad generator
#     location_of_jsub_factory_aao_generator = main_source_dir + "/aao_gen/gen_wrapper/batch_farm_executables/src/jsub_aao_generator.py"
#         # Jsub submitting tool
#     location_of_jsubmitter = main_source_dir+"/jlab_farm_tools/src/jsubmitter.py"
#         # aao_(no)rad generator wrapper aka aao_gen
#     location_of_aao_gen = main_source_dir+"/aao_gen/gen_wrapper/batch_farm_executables/src/aao_gen.py"
#         # actual generator location: aao_norad
#     location_of_aao_norad = main_source_dir+"/aao_gen/aao_norad/build/aao_norad.exe"
#         # actual generator location: aao_rad
#     location_of_aao_rad = main_source_dir+"/aao_gen/aao_rad/build/aao_rad.exe"


#         # input file maker: aao_norad and rad
#     location_of_input_file_maker = main_source_dir+"/aao_gen/gen_wrapper/batch_farm_executables/src/aao_input_file_maker.py"
#         # filter path for aao_rad and norad
#     location_of_lund_event_filter = main_source_dir+"/aao_gen/gen_wrapper/batch_farm_executables/src/lund_filter.py"


#         # dst copier path
#     location_of_dst_copier = main_source_dir+"/jlab_farm_tools/src/dst_copier_from_gemc_output.py"
#         # filter convert jsub machine
#     location_of_fc_jsub_machine = main_source_dir + "/jlab_farm_tools/src/jsub_filter_convert_machine.py"
#         #filter exe path
#     location_of_filter_exe = main_source_dir + "/filter/fiducial-filtering/filterEvents/"
#         #root combiner path
#     location_of_root_merger = main_source_dir + "/jlab_farm_tools/src/root_combiner.py"

#     location_of_converter_gen_exe = main_source_dir +  "/convertingHipo/minimal/convertGen"
#     location_of_converter_recon_exe = main_source_dir + "/convertingHipo/minimal/convertRec"

#     parser = argparse.ArgumentParser(description="""Need to write the description \n
#                     This script: \n
#                     1.) \n
#                     2.) \n
#                     3.) """,formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#     #Directory Structure
#     parser.add_argument("--base_dir",help="Location for directory containing all files",default='/volatile/clas12/robertej')
        
#         #Executables locations
#             #Jsub file creator for norad generator
#     parser.add_argument("--source_aao_jsub_generator",help="Location for norad generator jsub creator",default=location_of_jsub_factory_aao_generator)
#             #Jsub submitting tool
#     parser.add_argument("--jsubmitter",help="Location for jsubmission utility",default=location_of_jsubmitter)
#             # aao_(no)rad generator wrapper aka aao_gen
#     parser.add_argument("--pi0_gen_exe_path",help="Path to lund filter executable",default=location_of_aao_gen)

#             # actual generator location: aao_norad
#     parser.add_argument("--aao_norad_exe_location",help="Path to generator executable",default=location_of_aao_norad)
#             # actual generator location: aao_rad
#     parser.add_argument("--aao_rad_exe_location",help="Path to generator executable",default=location_of_aao_rad)



#     parser.add_argument("--generator_input_exe_path",help="Path to input file maker executable for aao_rad",default=location_of_input_file_maker)
#             # filter path for aao_norad
#     parser.add_argument("--lund_filter_exe_path",help="Path to lund filter executable",default=location_of_lund_event_filter)



#             # dst copier path
#     parser.add_argument("--dst_copier_path",help="Path to dst copier",default=location_of_dst_copier)
#             # Filter & Convert machine jsub path
#     parser.add_argument("--filt_conv_jsub_path",help="Location for filt-convert jsub creator",default=location_of_fc_jsub_machine)
#                 # Filter executable path
#     parser.add_argument("--filter_exe_path",help="Location for filter executable path",default=location_of_filter_exe)
#                 # gen converter executable path
#     parser.add_argument("--converter_gen_exe_path",help="Location for converter for gen executable path",default=location_of_converter_gen_exe)
#                     #recon converter executable path
#     parser.add_argument("--converter_recon_exe_path",help="Location for converter for recon executable path",default=location_of_converter_recon_exe)
#         #root merger path
#     parser.add_argument("--root_merger_path",help="Location for root merger script ",default=location_of_root_merger)
    
    

#             #This arguement can be ignored and should be deleted
#     parser.add_argument("--outdir",help="Location of intermediate return files between generation and filtering, can be ignored for batch farm",default="output/")


#     parser.add_argument("--generator_type",help="rad | norad, lets you build input for either aao_rad or aao_norad generators",default="norad")
#     parser.add_argument("--input_filename_rad",help="filename for aao_rad",default="aao_rad_input.inp")
#     parser.add_argument("--input_filename_norad",help="filename for aao_norad",default="aao_norad_input.inp")


#     with open('../../aao_gen/gen_wrapper/batch_farm_executables/src/default_generator_args.json') as fjson:
#         d = json.load(fjson)

#     norad = Dict2Class(d["aao_norad"][0])
#     rad = Dict2Class(d["aao_rad"][0])

   #File structure:
    # simulations_20210341014_2302
    # ├── 0_Jsub_factory
        # │   ├── Generation
            # │   │   └── jsub_gen_#.txt
        # │   ├── Filtering_Converting
            # │   ├── config_1
                # │   │   ├── filt_conv_recon
                    # │   │   └── jsub_fc_config_#_recon_#.txt
                # │   │   ├── filt_conv_gen
                    # │   │   └── jsub_fc_config_#_recon_#.txt
            # │   ├── config_2
                # │   │   ├── filt_conv_recon
                    # │   │   └── jsub_fc_config_#_recon_#.txt
                # │   │   ├── filt_conv_gen
                    # │   │   └── jsub_fc_config_#_recon_#.txt
            # │   ├── config_3
                # │   │   ├── filt_conv_recon
                    # │   │   └── jsub_fc_config_#_recon_#.txt
                # │   │   ├── filt_conv_gen
                    # │   │   └── jsub_fc_config_#_recon_#.txt
    # ├── 1_Lund files
        # └── xbmin_qtmin_#_norad.lund
    # ├── 2_GEMC DSTs
        # │   ├── config_1
            # │   │   └── dst_#.hipo
        # │   ├── config_2
            # │   │   └── dst_#.hipo
        # │   ├── config_3
        # │   │   └── dst_#.hipo
    # ├── 3_Filtered & Converted root files
        # │   ├── config_1
            # │   │   ├── filt_conv_recon root files
                # │   │   │   └── dst_fc_config_#_fc_recon_#.root
            # │   │   ├── filt_conv_recon root files
                # │   │   │   └── dst_fc_config_#__gen_#.root
        # │   ├── config_2
            # │   │   ├── filt_conv_recon root files
                # │   │   │   └── dst_fc_config_#_fc_recon_#.root
            # │   │   ├── filt_conv_recon root files
                # │   │   │   └── dst_fc_config_#__gen_#.root
        # │   ├── config_3
            # │   │   ├── filt_conv_recon root files
                # │   │   │   └── dst_fc_config_#_fc_recon_#.root
            # │   │   ├── filt_conv_recon root files
                # │   │   │   └── dst_fc_config_#__gen_#.root
    # ├── 4_Final files
        # │   ├── config_1
                # │   │   └── merged_recon.root
                # │   │   └── merged_gen.root
        # │   ├── config_2
            # │   │   └── merged_recon.root
            # │   │   └── merged_gen.root
        # │   ├── config_3
                # │   │   └── merged_recon.root
                # │   │   └── merged_gen.root
    
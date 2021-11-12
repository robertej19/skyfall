import uproot
import pandas as pd
import numpy as np
import argparse
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
pd.set_option('mode.chained_assignment', None)

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

def makeGenDVpi0vars(df_epgg):

    # useful objects
    ele = [df_epgg['GenEpx'], df_epgg['GenEpy'], df_epgg['GenEpz']]
    df_epgg.loc[:, 'GenEp'] = mag(ele)
    df_epgg.loc[:, 'GenEe'] = getEnergy(ele, me)
    df_epgg.loc[:, 'GenEtheta'] = getTheta(ele)
    df_epgg.loc[:, 'GenEphi'] = getPhi(ele)

    pro = [df_epgg['GenPpx'], df_epgg['GenPpy'], df_epgg['GenPpz']]
    df_epgg.loc[:, 'GenPp'] = mag(pro)
    df_epgg.loc[:, 'GenPe'] = getEnergy(pro, M)
    df_epgg.loc[:, 'GenPtheta'] = getTheta(pro)
    df_epgg.loc[:, 'GenPphi'] = getPhi(pro)

    gam = [df_epgg['GenGpx'], df_epgg['GenGpy'], df_epgg['GenGpz']]
    df_epgg.loc[:, 'GenGp'] = mag(gam)
    df_epgg.loc[:, 'GenGe'] = getEnergy(gam, 0)
    df_epgg.loc[:, 'GenGtheta'] = getTheta(gam)
    df_epgg.loc[:, 'GenGphi'] = getPhi(gam)

    gam2 = [df_epgg['GenGpx2'], df_epgg['GenGpy2'], df_epgg['GenGpz2']]
    df_epgg.loc[:, 'GenGp2'] = mag(gam2)
    df_epgg.loc[:,'GenGe2'] = getEnergy(gam2, 0)
    df_epgg.loc[:, 'GenGtheta2'] = getTheta(gam2)
    df_epgg.loc[:, 'GenGphi2'] = getPhi(gam2)

    pi0 = vecAdd(gam, gam2)
    VGS = [-df_epgg['GenEpx'], -df_epgg['GenEpy'], pbeam - df_epgg['GenEpz']]
    v3l = cross(beam, ele)
    v3h = cross(pro, VGS)
    v3g = cross(VGS, gam)
    VmissPi0 = [-df_epgg["GenEpx"] - df_epgg["GenPpx"], -df_epgg["GenEpy"] -
                df_epgg["GenPpy"], pbeam - df_epgg["GenEpz"] - df_epgg["GenPpz"]]
    VmissP = [-df_epgg["GenEpx"] - df_epgg["GenGpx"] - df_epgg["GenGpx2"], -df_epgg["GenEpy"] -
                df_epgg["GenGpy"] - df_epgg["GenGpy2"], pbeam - df_epgg["GenEpz"] - df_epgg["GenGpz"] - df_epgg["GenGpz2"]]
    Vmiss = [-df_epgg["GenEpx"] - df_epgg["GenPpx"] - df_epgg["GenGpx"] - df_epgg["GenGpx2"],
                -df_epgg["GenEpy"] - df_epgg["GenPpy"] - df_epgg["GenGpy"] - df_epgg["GenGpy2"],
                pbeam - df_epgg["GenEpz"] - df_epgg["GenPpz"] - df_epgg["GenGpz"] - df_epgg["GenGpz2"]]

    df_epgg.loc[:, 'GenMpx'], df_epgg.loc[:, 'GenMpy'], df_epgg.loc[:, 'GenMpz'] = Vmiss

    #print(2)
    # binning kinematics
    df_epgg.loc[:,'GenQ2'] = -((ebeam - df_epgg['GenEe'])**2 - mag2(VGS))
    df_epgg.loc[:,'Gennu'] = (ebeam - df_epgg['GenEe'])
    df_epgg.loc[:,'GenxB'] = df_epgg['GenQ2'] / 2.0 / M / df_epgg['Gennu']
    df_epgg.loc[:,'Gent'] = 2 * M * (df_epgg['GenPe'] - M)
    df_epgg.loc[:,'GenW'] = np.sqrt(np.maximum(0, (ebeam + M - df_epgg['GenEe'])**2 - mag2(VGS)))
    df_epgg.loc[:,'GenMPt'] = np.sqrt((df_epgg["GenEpx"] + df_epgg["GenPpx"] + df_epgg["GenGpx"] + df_epgg["GenGpx2"])**2 +
                                (df_epgg["GenEpy"] + df_epgg["GenPpy"] + df_epgg["GenGpy"] + df_epgg["GenGpy2"])**2)

    #print(3)
    # trento angles
    df_epgg.loc[:,'Genphi1'] = angle(v3l, v3h)
    #print(3.1)

    df_epgg.loc[:,'Genphi1'] = np.where(dot(v3l, pro) > 0, 360.0 -
                                df_epgg['Genphi1'], df_epgg['Genphi1'])
    #print(3.2)
    
    df_epgg.loc[:,'Genphi2'] = angle(v3l, v3g)
    #print(3.3)

    df_epgg.loc[:,'Genphi2'] = np.where(dot(VGS, cross(v3l, v3g)) <
                                0, 360.0 - df_epgg['Genphi2'], df_epgg['Genphi2'])

    #print(3.4)

    return df_epgg

def makeDVpi0vars(df_epgg):

    # useful objects
    ele = [df_epgg['Epx'], df_epgg['Epy'], df_epgg['Epz']]
    df_epgg.loc[:, 'Ep'] = mag(ele)
    df_epgg.loc[:, 'Ee'] = getEnergy(ele, me)
    df_epgg.loc[:, 'Etheta'] = getTheta(ele)
    df_epgg.loc[:, 'Ephi'] = getPhi(ele)

    pro = [df_epgg['Ppx'], df_epgg['Ppy'], df_epgg['Ppz']]
    df_epgg.loc[:, 'Pp'] = mag(pro)
    df_epgg.loc[:, 'Pe'] = getEnergy(pro, M)
    df_epgg.loc[:, 'Ptheta'] = getTheta(pro)
    df_epgg.loc[:, 'Pphi'] = getPhi(pro)

    gam = [df_epgg['Gpx'], df_epgg['Gpy'], df_epgg['Gpz']]
    df_epgg.loc[:, 'Gp'] = mag(gam)
    df_epgg.loc[:, 'Ge'] = getEnergy(gam, 0)
    df_epgg.loc[:, 'Gtheta'] = getTheta(gam)
    df_epgg.loc[:, 'Gphi'] = getPhi(gam)

    gam2 = [df_epgg['Gpx2'], df_epgg['Gpy2'], df_epgg['Gpz2']]
    df_epgg.loc[:, 'Gp2'] = mag(gam2)
    df_epgg.loc[:,'Ge2'] = getEnergy(gam2, 0)
    df_epgg.loc[:, 'Gtheta2'] = getTheta(gam2)
    df_epgg.loc[:, 'Gphi2'] = getPhi(gam2)

    pi0 = vecAdd(gam, gam2)
    VGS = [-df_epgg['Epx'], -df_epgg['Epy'], pbeam - df_epgg['Epz']]
    v3l = cross(beam, ele)
    v3h = cross(pro, VGS)
    v3g = cross(VGS, gam)
    VmissPi0 = [-df_epgg["Epx"] - df_epgg["Ppx"], -df_epgg["Epy"] -
                df_epgg["Ppy"], pbeam - df_epgg["Epz"] - df_epgg["Ppz"]]
    VmissP = [-df_epgg["Epx"] - df_epgg["Gpx"] - df_epgg["Gpx2"], -df_epgg["Epy"] -
                df_epgg["Gpy"] - df_epgg["Gpy2"], pbeam - df_epgg["Epz"] - df_epgg["Gpz"] - df_epgg["Gpz2"]]
    Vmiss = [-df_epgg["Epx"] - df_epgg["Ppx"] - df_epgg["Gpx"] - df_epgg["Gpx2"],
                -df_epgg["Epy"] - df_epgg["Ppy"] - df_epgg["Gpy"] - df_epgg["Gpy2"],
                pbeam - df_epgg["Epz"] - df_epgg["Ppz"] - df_epgg["Gpz"] - df_epgg["Gpz2"]]

    df_epgg.loc[:, 'Mpx'], df_epgg.loc[:, 'Mpy'], df_epgg.loc[:, 'Mpz'] = Vmiss

    # binning kinematics
    df_epgg.loc[:,'Q2'] = -((ebeam - df_epgg['Ee'])**2 - mag2(VGS))
    df_epgg.loc[:,'nu'] = (ebeam - df_epgg['Ee'])
    df_epgg.loc[:,'xB'] = df_epgg['Q2'] / 2.0 / M / df_epgg['nu']
    df_epgg.loc[:,'t'] = 2 * M * (df_epgg['Pe'] - M)
    df_epgg.loc[:,'W'] = np.sqrt(np.maximum(0, (ebeam + M - df_epgg['Ee'])**2 - mag2(VGS)))
    df_epgg.loc[:,'MPt'] = np.sqrt((df_epgg["Epx"] + df_epgg["Ppx"] + df_epgg["Gpx"] + df_epgg["Gpx2"])**2 +
                                (df_epgg["Epy"] + df_epgg["Ppy"] + df_epgg["Gpy"] + df_epgg["Gpy2"])**2)

    # trento angles
    df_epgg['phi1'] = angle(v3l, v3h)
    df_epgg['phi1'] = np.where(dot(v3l, pro) > 0, 360.0 -
                                df_epgg['phi1'], df_epgg['phi1'])
    df_epgg['phi2'] = angle(v3l, v3g)
    df_epgg['phi2'] = np.where(dot(VGS, cross(v3l, v3g)) <
                                0, 360.0 - df_epgg['phi2'], df_epgg['phi2'])


    # exclusivity variables
    df_epgg.loc[:,'MM2_ep'] = (-M - ebeam + df_epgg["Ee"] +
                            df_epgg["Pe"])**2 - mag2(VmissPi0)
    df_epgg.loc[:,'MM2_egg'] = (-M - ebeam + df_epgg["Ee"] +
                            df_epgg["Ge"] + df_epgg["Ge2"])**2 - mag2(VmissP)
    df_epgg.loc[:,'MM2_epgg'] = (-M - ebeam + df_epgg["Ee"] + df_epgg["Pe"] +
                            df_epgg["Ge"] + df_epgg["Ge2"])**2 - mag2(Vmiss)
    df_epgg.loc[:,'ME_epgg'] = (M + ebeam - df_epgg["Ee"] - df_epgg["Pe"] - df_epgg["Ge"] - df_epgg["Ge2"])
    df_epgg.loc[:,'Mpi0'] = pi0InvMass(gam, gam2)
    df_epgg.loc[:,'reconPiAngleDiff'] = angle(VmissPi0, pi0)
    df_epgg.loc[:,"Pie"] = df_epgg['Ge'] + df_epgg['Ge2']
    
    df_math_epgg = df_epgg

    return df_math_epgg

def cutDVpi(df_epgg):
    #make dvpi0 pairs

    df_epgg.loc[:, "closeness"] = np.abs(df_epgg.loc[:, "Mpi0"] - .1349766)

    cut_xBupper = df_epgg.loc[:, "xB"] < 1  # xB
    cut_xBlower = df_epgg.loc[:, "xB"] > 0  # xB
    cut_Q2 = df_epgg.loc[:, "Q2"] > 1  # Q2
    cut_W = df_epgg.loc[:, "W"] > 2  # W

    # Exclusivity cuts
    cut_mmep = df_epgg.loc[:, "MM2_ep"] < 0.7  # mmep
    cut_meepgg = df_epgg.loc[:, "ME_epgg"] < 0.7  # meepgg
    cut_mpt = df_epgg.loc[:, "MPt"] < 0.2  # mpt
    cut_recon = df_epgg.loc[:, "reconPiAngleDiff"] < 2  # recon gam angle
    cut_pi0upper = df_epgg.loc[:, "Mpi0"] < 0.2
    cut_pi0lower = df_epgg.loc[:, "Mpi0"] > 0.07
    cut_sector = (df_epgg.loc[:, "Esector"]!=df_epgg.loc[:, "Gsector"]) & (df_epgg.loc[:, "Esector"]!=df_epgg.loc[:, "Gsector2"])

    df_dvpi0 = df_epgg.loc[cut_xBupper & cut_xBlower & cut_Q2 & cut_W & cut_mmep & cut_meepgg &
                        cut_mpt & cut_recon & cut_pi0upper & cut_pi0lower & cut_sector, :]

    #For an event, there can be two gg's passed conditions above.
    #Take only one gg's that makes pi0 invariant mass
    #This case is very rare.
    #For now, duplicated proton is not considered.
    df_dvpi0.sort_values(by='closeness', ascending=False)
    df_dvpi0.sort_values(by='event')        
    df_dvpi0 = df_dvpi0.loc[~df_dvpi0.event.duplicated(), :]

    #df_x = df_dvpi0.loc[:, ["event", "Epx", "Epy", "Epz", "Ep", "Ephi", "Etheta", "Ppx", "Ppy", "Ppz", "Pp", "Pphi", "Ptheta", "Gpx", "Gpy", "Gpz", "Gp", "Gtheta", "Gphi", "Gpx2", "Gpy2", "Gpz2", "Gp2", "Gtheta2", "Gphi2"]]
    #self.df_x = df_x #done with saving x

    return df_dvpi0


def get_counts(df_base,tmin=0,tmax=1,xbmin=0,xbmax=1,q2min=0,q2max=12,datatype="Recon"):
    cut_q = "xB>{} & xB<{} & Q2>{} & Q2<{} & t>{} & t<{}".format(xbmin,xbmax,q2min,q2max,tmin,tmax)
    x_var_name = "phi1"
    if datatype == "Gen":
        cut_q = "GenxB>{} & GenxB<{} & GenQ2>{} & GenQ2<{} & Gent>{} & Gent<{}".format(xbmin,xbmax,q2min,q2max,tmin,tmax)
        x_var_name ="Genphi1"

    df = df_base.query(cut_q)
    ic(df)

    x_data = df[x_var_name]
    
    var_names = ["$\phi$"]
    ranges = [0,360,20]
    output_dir = "pics/"
    title = "$\phi$, Sim, {}<t<{} GeV$^2$,{}<$x_B$<{}, {}<$Q^2$<{}".format(tmin,tmax,xbmin,xbmax,q2min,q2max)
    if not (datatype=="Gen"):
        make_histos.plot_1dhist(x_data,var_names,ranges,
                    saveplot=True,pics_dir=output_dir,plot_title=title.replace("/",""),first_color="darkslateblue")

    count, division = np.histogram(x_data, bins = fs.phibins)
    tmin_arr = tmin*np.ones(len(count))

    mean_g, mean_epsi = 0,0
    if datatype =="Real":
        mean_g = df['gamma'].mean()*np.ones(len(count))
        mean_epsi = df['epsi'].mean()*np.ones(len(count))
    
    return count, tmin_arr, division, mean_g, mean_epsi



if __name__ == "__main__":

    #To process suite:
    ## If processing data for the first time:
    # python pickle_analysis.py -c
    # python pickle_analysis.py -g -c
    # python pickle_analysis.py -g
    # python pickle_analysis.py -r -c

    ## If running after the first time:
    # python pickle_analysis.py
    # python pickle_analysis.py -g
    # python pickle_analysis.py -r

    size_gen_chunks = 4000000 #Change this, depending on how big the gen dataset is


    parser = argparse.ArgumentParser(description="Get args",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-f","--fname", help="a single root file to convert into pickles", default="infile.root")
    parser.add_argument("-o","--out", help="a single pickle file name as an output", default="outfile.pkl")
    parser.add_argument("-s","--entry_stop", help="entry_stop to stop reading the root file", default = None)
    parser.add_argument("-c","--cut", help="use this flag to cut out non-DVPiP events", default=False, action='store_true')
    parser.add_argument("-g","--gen", help="enable to use gen events instead of recon", default=False, action='store_true')
    parser.add_argument("-r","--real", help="enable to use real events instead of recon", default=False, action='store_true')
    parser.add_argument("-p","--plot", help="enable to use real events instead of recon", default=False, action='store_true')
    
    args = parser.parse_args()

    dfs = []
    gen_path = "data/after_cuts/gen/"

    datatype = "Recon"
    if args.gen:
        datatype = "Gen"
    if args.real:
        datatype = "Real"


    if args.gen:
        #fname0 = "df_7999_genONLY"
        #fname0 = "df_gentest3"
        #fname0 = "test34gen"
        #fname0 = "df_in_18"
        #fname0 = "df_test_10M_gen"
        #fname0 = "df_gen_pisamp"
        #fname0 = "merged_Fall_2018_Inbending_gen_10radtest_genONLY"
        #fname0 = "merged_Fall_2018_Inbending_recon_10radtest_recon_gen"
        #fname0 = "df_gen_TEST2"
        #fname0 = "df_in_18_recon_gen"
        if args.cut:
            #df_gen = pd.read_pickle("data/before_cuts/df_gen.pkl")
            #df = pd.read_pickle("data/before_cuts/{}.pkl".format(fname0))
            df = pd.read_pickle("/mnt/d/GLOBUS/CLAS12/simulations/production/Fall_2018_Inbending/Test/ACTUAL_Gen_test_norad_gen_genONLY.pkl")
            n = size_gen_chunks  #chunk row size
            ic(df.shape)
            list_df = []
            for i in range(0,df.shape[0],n):
                ic(i)
                ic(i+n)
                list_df.append(df[i:i+n])

            for index, df_chunk in enumerate(list_df):
                print("On DF chunk {}".format(index))
                df_gen = makeGenDVpi0vars(df_chunk)
                #df_gen.to_pickle("data/after_cuts/gen/{}_with_cuts_{}.pkl".format(fname0,index))
                dfs.append(df_gen)
        else:
            gen_files = os.listdir(gen_path)
            for gf in gen_files:
                dfs.append(pd.read_pickle(gen_path+gf))
    elif args.real:
        if args.cut:
            df_after_cuts = pd.read_pickle("data/after_cuts/F18_All_DVPi0_Events.pkl")
            df_after_cuts.loc[:, "y"] = (E-df_after_cuts.loc[:, "nu"])/E    
            df_after_cuts.loc[:, "q24E2"] = df_after_cuts.loc[:, "Q2"]/(4*E*E)
            df_after_cuts.loc[:, "epsi"] = (1-df_after_cuts.loc[:, "y"]-df_after_cuts.loc[:, "q24E2"])/(1-df_after_cuts.loc[:, "y"]+(df_after_cuts.loc[:, "q24E2"]*df_after_cuts.loc[:, "q24E2"])/2+df_after_cuts.loc[:, "q24E2"])
            df_after_cuts.loc[:, "gamma"] = prefix*df_after_cuts.loc[:, "Q2"]/(mp*mp*E*E)*(1-df_after_cuts.loc[:,"xB"])/(df_after_cuts.loc[:,"xB"]**3)*(1/(1- df_after_cuts.loc[:, "epsi"]))/(2*np.pi)
            df_after_cuts.to_pickle("data/after_cuts/F18_All_DVPi0_Events_with_extra_vars.pkl")
        else:
            df_after_cuts = pd.read_pickle("data/after_cuts/F18_All_DVPi0_Events_with_extra_vars.pkl")
        dfs.append(df_after_cuts)         
    else:
        #fname0 = "data/before_cuts/df_recon_7999.pkl"
        #fname1 = "data/after_cuts/df_recon_7999_with_cuts.pkl"
        #fname0 = "df_radtest7_recon.pkl"
        #fname1 = "df_radtest7_recon_with_cuts.pkl"
        #fname0 = "data/before_cuts/df_in_18_recon_gen.pkl"
        #fname0 = "df_recon_TEST2.pkl"
        #fname0 = "df_test_outbending_recon.pkl"
        fname0 = "merged_Fall_2018_Inbending_recon_10radtest_recon_recon.pkl"

        if args.cut:
            df_recon_0 = pd.read_pickle("data/before_cuts/{}".format(fname0))
            #Calculate pi0 parameters    
            df_recon_pi0vars = makeDVpi0vars(df_recon_0)
            #Apply exclusivity cuts    
            df_recon = cutDVpi(df_recon_pi0vars)
            df_recon.to_pickle("data/after_cuts/{}".format(fname0))
        else:
            #df_recon = pd.read_pickle("data/after_cuts/df_recon_with_cuts.pkl")
            print('getting df recon')
            df_recon = pd.read_pickle(fname1)
        dfs.append(df_recon)

    ic(dfs)
    df = dfs[0]


    print(df.columns.values)
    if args.plot:
        df = dfs[0]

        # df = df.query("GenW>2")
        # df = df.query("GenQ2>1")
        ic(df)


        histo_plotting.make_all_histos(df,datatype=datatype,hists_2d=False,hists_1d=True,hists_overlap=False,saveplots=False)
        sys.exit()

    for df_index, df in enumerate(dfs):
        print("Processing df {}".format(df_index))
        ic(df)
        binned_dfs = []
        xbmin = 0.3
        xbmax = 0.38
        q2min = 2
        q2max = 2.5
        for t_set in fs.t_ranges:
            tmin = t_set[0]
            tmax = t_set[1]
            count, tmin_arr, division, mean_g, mean_epsi = get_counts(df_base=df,tmin=tmin,tmax=tmax,xbmin=xbmin,xbmax=xbmax,q2min=q2min,q2max=q2max,datatype=datatype)
            #binned = pd.DataFrame(data=tmin_arr,index=division[:-1],columns=['tmin'])
            binned = pd.DataFrame({'phi_min':division[:-1],'tmin':tmin_arr})
            if datatype =="Real":
                binned['real_counts'] = count
                binned['gamma'] = mean_g
                binned['epsi'] = mean_epsi
                binned['xb_min'] = xbmin
                binned['xb_max'] = xbmax
                binned['q2_min'] = q2min
                binned['q2_max'] = q2max
                binned['t_max'] = tmax
            elif datatype =="Gen":
                binned['gen_counts_{}'.format(df_index)] = count
            else:
                binned['recon_counts'] = count
            ic(binned)
            binned_dfs.append(binned)

        real_out = pd.concat(binned_dfs,axis='index',ignore_index=True,join="outer")
        ic(real_out)
        real_out.to_pickle("data/binned/{}_phi_binned_{}.pkl".format(datatype,df_index))
        

    email = True
    if email:
        import os
        from datetime import datetime
        from pytools import circle_emailer

        now = datetime.now()
        script_end_time = now.strftime("%H:%M:%S")
        s_name = os.path.basename(__file__)
        subject = "Completion of {}".format(s_name)
        body = "Your script {} finished running at {}".format(s_name,script_end_time)
        circle_emailer.send_email(subject,body)



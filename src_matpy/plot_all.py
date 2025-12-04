# Plot all data
# Version: Mar-04-2025
# Import modules
# Note the experimental fraction of 0.25 is actually 0.16
# So the files are named with 0.25 at the end. But the legends are
# adjusted acoordingly
# Similarly fraction of 0.11 is actually 0.09
import os
import math
import glob
import numpy as np
import analyze_rdfs as ardf
import pandas as pd
import aux_plot as paux
from env_plot import *
from matplotlib.ticker import FuncFormatter, LogLocator  # add these
#---------------Directory paths-------------------------------------
main_dir = os.getcwd() # current dir

all_dir  = '../../new_results' # analysis dir
if not os.path.isdir(all_dir):
    print("FATAL ERROR: ", all_dir, " not found")
    exit("Check scratch directory path")
    
anaout_dir = '../../analyzed_results'
if not os.path.isdir(anaout_dir):
    os.mkdir(anaout_dir)
    
figout_dir = '../../figures'
if not os.path.isdir(figout_dir):
    os.mkdir(figout_dir)

#---------Input Data-----------------------------------------------
syspref  = 'qanion_'
respref  = 'allresults_'
senspref = 'sens_' 
q_anion  = '0.2' # should be the exact string value
sys_arr  = ['S2'] # Give one at a time
fran_arr = ['0.05','0.09','0.20']
f_unpoly = [0, 0.6, 0.6]
deltat   = [0.006, 0.006, 0.006]
case_arr = [15]
molwt    = 40
nchains  = 100
S2_vol_tot  = [9926.68,9884.76,10032.43]
S3_vol_tot  = [9856.87,9807.30,9809.80]               

#---------File prefixes------------------------------------------
rdf_pref        = 'rdf_config_'
rg_pref         = 'rgavg_config_'
neigh_pref      = 'catanneigh_config_'
restime_pref    = 'autocorrcion_config'
diff_pref       = ['iondiff_config_', 'countiondiff_config_']
clust_pref      = 'clust_config_'
bfrdf_pref      = 'freeboundrdf_config_'
neigh_sens_pref = 'catanneigh_'
clust_sens_pref = 'clustavg_'

#---------Flags---------------------------------------------------
rdf_flag     = 0; rdf_keys    = {'LiAn':4,'LiPoly':3,'LiMon':7}
rg_flag      = 0; rg_keys     = {'rgt': 2} 
neigh_flag   = 0; neigh_keys  = {'Li-An':3,'An-Li':5}; maxneigh = 6
clust_flag   = 0; clust_keys  = {'Li-An':2}
bfrdf_flag   = 0
restime_flag = 0; restime_keys = {'Li-An':2}
simdiff_flag = 1 # Simulation diffusivities and transference numbers
simdiffavg_flag = 1 # Avg diffusivities with errbar
n_sens_flag  = 0 # sensitivity analysis on neighbors
c_sens_flag  = 0 # sensitivity analysis on clusters
s_clus_flag  = 0 # sum of clusters - reviewer's comment
n_surf_flag  = 0 # surface plot of neighbors - reviewer's comment
msd_flag     = 0 # MSD from simulations and MSD/t vs t
tg_flag      = 0 # Tg calculation
cross_flag   = 0 # Plot cross-MSD
expnmr_flag  = 0 # Experimental NMR data plot in semi-log scale
expdiff_flag = 0 # Experimental diffusivities and transference numbers
exppoly_flag = 0 # Experimental polymerization/unpolymerization data
expcond_flag = 0 # Experimental conductivities
exptnum_Ndep = 0 # Experimental dependence on transference numbers


#---------Experimental inputs---------------------------------------
fanfrac_values = [0.047, 0.091, 0.167]
fanfrac_labels = ['0.047', '0.091', '0.167']
expinp_dir = '../../ExperimentalResults'
if not os.path.isdir(expinp_dir) and expcond_flag == 1:
    print("FATAL ERROR: ", expinp_dir, " not found")
    exit("Check scratch directory path")

#---------Calculations----------------------------------------------
nanions  = [nchains*int(molwt*float(x)) for x in fran_arr]
npols    = [nchains*molwt - nanions[i]  for i in range(len(nanions))]
totbeads = [int(2.0*(npols[i] + nanions[i])) for i in range(len(npols))]
npol_ans = [int(0.5*float(x)) for x in npols]

#---------Initialize global arrays--------------------------------
if restime_flag:
    meantau_arr = np.empty(len(sys_arr),len(fran_arr))            

#-----Sanity checks------------------------------------------------
if len(fran_arr) != len(deltat):
    raise RuntimeError('Mismatch in lengths between delta_t and fran_arr')

#-------Main analysis----------------------------------------------           
if(restime_flag):

    print("Plotting average relaxation time")
    fig,ax = plt.subplots() 
    paux.set_axes(ax,plt,r'$f_{an}$',r'Mean Relaxation Time$($\tau$)')

    # Bar settings
    x = np.arange(len(sys_arr))  # the label locations
    bar_width = 0.2

    for i, label in enumerate(fran_arr):
        ax.bar(x + i * bar_width, [row[i] for row in meantau_arr],\
               bar_width,label='$f_{an}$ = ' + label,color=clr_arr[i])
        
    # Labels, legend, and formatting
    ax.set_xticks(x + bar_width)
    ax.set_xticklabels(categories)
    ax.legend()
    fig.savefig(f'{figout_dir}/meanrelax_ion.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/meanrelax_ion.eps',format = 'eps')

                
if(simdiff_flag):

    column_names = [
        "System", "f_an", "Case#", "Traj_#", "Timestep", "deltaT", "D_Li", "D_STFSI", 
        "t_Li_z=1", "t_Li_t=f_an*N", "Li_Power", "TFSI_Power"
    ]

    # Load the data into a DataFrame
    df = pd.read_csv('../../analyzed_results/computed_diffusivity_' + \
                     str(case_arr[0]) + '.dat', sep="\t", skipinitialspace=True)
    df.columns = column_names

    # Extract relevant columns and ensure correct datatypes
    df["f_an"] = df["f_an"].astype(float)
    df["D_Li"] = df["D_Li"].astype(float)
    df["D_STFSI"] = df["D_STFSI"].astype(float)
    df["t_Liz1"] = df["t_Li_z=1"].astype(float)
    df["t_LizfN"] = df["t_Li_t=f_an*N"].astype(float)


    # Rename columns for plotting
    plot_labels = {"S1": "RCP", "S2": "RCP+M", "S3": "HP+M"}
    
    # Pivot the table for grouped bar chart
    lipivot = df.pivot(index="System", columns="f_an", values="D_Li")
    lipivot = lipivot.sort_index()
    lipivot = lipivot.rename(index=plot_labels)

    stfsipivot = df.pivot(index="System", columns="f_an", values="D_STFSI")
    stfsipivot = stfsipivot.sort_index()
    stfsipivot = stfsipivot.rename(index=plot_labels)
    
    tz1lipivot = df.pivot(index="System", columns="f_an", values="t_Liz1")
    tz1lipivot = tz1lipivot.sort_index()
    tz1lipivot = tz1lipivot.rename(index=plot_labels)

    tzfNlipivot = df.pivot(index="System", columns="f_an", values="t_LizfN")
    tzfNlipivot = tzfNlipivot.sort_index()
    tzfNlipivot = tzfNlipivot.rename(index=plot_labels)

   
    # Bar settings
    x = np.arange(len(lipivot.index))  # the label locations
    bar_width = 0.2

    #------------Plotting cation and anion stacked-------------------
    print("Plotting stacked cation and anion data...")
    fig, ax = plt.subplots()
    paux.set_axes(ax,plt,r'Model',r'Diffusion Coefficient (${\sigma^2}/{\tau}$)')
    for i, label in enumerate(lipivot.columns):
        bars = ax.bar(x + i * bar_width, lipivot[label],width = bar_width,\
                      label='$D_{Li^{+}}$, $f_{an}$ = ' + str(label),color=clr_arr[i])


    for i, label in enumerate(stfsipivot.columns):
        bars = ax.bar(x + i * bar_width, stfsipivot[label],width = bar_width,\
                      label='$D_{STFSI^{-}}$, $f_{an}$ = ' + str(label),color=clr_arr[i],\
                      hatch='//', edgecolor='black', fill=False)

    ax.set_xticks(x + bar_width)
    ax.set_xticklabels(stfsipivot.index)
    ax.set_ylim([0, 0.06])
#    ax.legend(loc='upper center', ncol=3, frameon=False)
#    ax.tick_params(axis='x',pad=40)
    ax.legend()
    fig.savefig(f'{figout_dir}/diff_all_stacked_case' + str(case_arr[0]) + '.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/diff_all_stacked_case' + str(case_arr[0]) + '.eps',format = 'eps')

    
    #--------------Plotting ion-diffusivity------------------
    print("Plotting ion-diffusivity")    
    fig, ax = plt.subplots()
    paux.set_axes(ax,plt,r'Model',r'$D_{Li^{+}}$ (${\sigma^2}/{\tau}$)')

    for i, label in enumerate(lipivot.columns):
        bars = ax.bar(x + i * bar_width, lipivot[label],width = bar_width,\
                      label='$f_{an}$ = ' + str(label),color=clr_arr[i])
        
        # Annotate each bar with the corresponding f_an value
        # for j, bar in enumerate(bars):
        #     f_an_val = lipivot.columns[i]
        #     height = bar.get_height()
        #    ax.text(bar.get_x() + bar.get_width()/2, height,
        #            f"{f_an_val:.2f}", ha='center', va='bottom', fontsize=9)
        
    # Labels, legend, and formatting
    ax.set_xticks(x + bar_width)
    ax.set_xticklabels(lipivot.index)
    ax.tick_params(axis='x',pad=20)
    ax.legend()
    fig.savefig(f'{figout_dir}/diff_ion_case' + str(case_arr[0]) + '.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/diff_ion_case' + str(case_arr[0]) + '.eps',format = 'eps')

    #--------Plotting counter-ion-diffusivity------------------
    print("Plotting counter-ion diffusivity")    
    fig, ax = plt.subplots()
    paux.set_axes(ax,plt,r'Model',r'$D_{STFSI^{-}}$ (${\sigma^2}/{\tau}$)')
    for i, label in enumerate(stfsipivot.columns):
        bars = ax.bar(x + i * bar_width, stfsipivot[label],width = bar_width,\
                      label='$f_{an}$ = ' + str(label),color=clr_arr[i])

        # Annotate each bar with the corresponding f_an value
        # for j, bar in enumerate(bars):
        #     f_an_val = stfsipivot.columns[i]
        #     height = bar.get_height()
        #     ax.text(bar.get_x() + bar.get_width()/2, height,
        #             f"{f_an_val:.2f}", ha='center', va='bottom', fontsize=9)

        
    # Labels, legend, and formatting
    ax.set_xticks(x + bar_width)
    ax.set_xticklabels(stfsipivot.index)
    ax.set_yscale('log')
    ax.legend()
    fig.savefig(f'{figout_dir}/diff_countion_case' + str(case_arr[0]) + '.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/diff_countion_case' + str(case_arr[0]) + '.eps',format = 'eps')

    #--------Plotting Transference Number (z = 1) --------------------
    print("Plotting transference number (z=1)")    
    fig, ax = plt.subplots()
    paux.set_axes(ax,plt,r'Model',r'$t_{+}$')
    for i, label in enumerate(tz1lipivot.columns):
        bars = ax.bar(x + i * bar_width, tz1lipivot[label],width = bar_width,\
                      label='$f_{an}$ = ' + str(label),color=clr_arr[i])

        # Annotate each bar with the corresponding f_an value
        #for j, bar in enumerate(bars):
        #    f_an_val = tz1lipivot.columns[i]
        #    height = bar.get_height()
        #    ax.text(bar.get_x() + bar.get_width()/2, height,
        #            f"{f_an_val:.2f}", ha='center', va='bottom', fontsize=9)
        
    # Labels, legend, and formatting
    ax.set_ylim([0.6, 1.05])
    ax.set_xticks(x + bar_width)
    ax.set_xticklabels(tz1lipivot.index)
    ax.legend()
    fig.savefig(f'{figout_dir}/tplus_case' + str(case_arr[0]) + '.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/tplus_case' + str(case_arr[0]) + '.eps',format = 'eps')

    #--------Plotting Transference Number (z = fN) --------------------
    print("Plotting transference number (z = fN)")    
    fig, ax = plt.subplots()
    paux.set_axes(ax,plt,r'Model',r'$t_{+}$')

    for i, label in enumerate(tzfNlipivot.columns):
        ax.bar(x + i * bar_width, tzfNlipivot[label],width = bar_width,\
               label='$f_{an}$ = ' + str(label),color=clr_arr[i])

        
    # Labels, legend, and formatting
    ax.set_ylim([0.3, 1.0])
    ax.set_xticks(x + bar_width)
    ax.set_xticklabels(tzfNlipivot.index)
    ax.legend()
    fig.savefig(f'{figout_dir}/tplus_scaledfN_case' + str(case_arr[0]) + '.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/tplus_scaledfN_case' + str(case_arr[0]) + '.eps',format = 'eps')

if(simdiffavg_flag):

    column_names = [
        "System", "f_an", "Case#", "Traj_#", "Timestep", "deltaT", "D_Li", "D_STFSI", 
        "t_Li_z=1", "t_Li_t=f_an*N", "Li_Power", "TFSI_Power","ebar_D_Li","ebar_D_STFSI","ebar_t_fN"

    ]

    # Load the data into a DataFrame
    df = pd.read_csv('../../analyzed_results/computed_diffusivity_avg.txt', \
                     sep="\t", skipinitialspace=True)
    df.columns = column_names

    # Extract relevant columns and ensure correct datatypes
    df["f_an"] = df["f_an"].astype(float)
    df["D_Li"] = df["D_Li"].astype(float)
    df["D_STFSI"] = df["D_STFSI"].astype(float)
    df["t_Liz1"] = df["t_Li_z=1"].astype(float)
    df["t_LizfN"] = df["t_Li_t=f_an*N"].astype(float)
    df["ebar_D_Li"] = df["ebar_D_Li"].astype(float)
    df["ebar_D_STFSI"] = df["ebar_D_STFSI"].astype(float)
    df["ebar_t_fN"] = df["ebar_t_fN"].astype(float)
        
    # Rename columns for plotting
    plot_labels = {"S1": "RCP", "S2": "RCP+M", "S3": "HP+M"}
    
    # Pivot the table for grouped bar chart
    lipivot = df.pivot(index="System", columns="f_an", values="D_Li")
    lipivot = lipivot.sort_index()
    lipivot = lipivot.rename(index=plot_labels)

    stfsipivot = df.pivot(index="System", columns="f_an", values="D_STFSI")
    stfsipivot = stfsipivot.sort_index()
    stfsipivot = stfsipivot.rename(index=plot_labels)
    
    tz1lipivot = df.pivot(index="System", columns="f_an", values="t_Liz1")
    tz1lipivot = tz1lipivot.sort_index()
    tz1lipivot = tz1lipivot.rename(index=plot_labels)

    tzfNlipivot = df.pivot(index="System", columns="f_an", values="t_LizfN")
    tzfNlipivot = tzfNlipivot.sort_index()
    tzfNlipivot = tzfNlipivot.rename(index=plot_labels)

    ebDlipivot = df.pivot(index="System", columns="f_an", values="ebar_D_Li")
    ebDlipivot = ebDlipivot.sort_index()
    ebDlipivot = ebDlipivot.rename(index=plot_labels)

    ebDstfsipivot = df.pivot(index="System", columns="f_an", values="ebar_D_STFSI")
    ebDstfsipivot = ebDstfsipivot.sort_index()
    ebDstfsipivot = ebDstfsipivot.rename(index=plot_labels)

    ebtzfNpivot = df.pivot(index="System", columns="f_an", values="ebar_t_fN")
    ebtzfNpivot = ebtzfNpivot.sort_index()
    ebtzfNpivot = ebtzfNpivot.rename(index=plot_labels)
    


   
    # Bar settings
    x = np.arange(len(lipivot.index))  # the label locations
    bar_width = 0.2

    #------------Plotting cation and anion stacked-------------------
    print("Plotting stacked cation and anion data...")
    fig, ax = plt.subplots()
    paux.set_axes(ax,plt,r'Model',r'Diffusion Coefficient (${\sigma^2}/{\tau}$)')
    for i, label in enumerate(lipivot.columns):
        bars = ax.bar(x + i * bar_width, lipivot[label], yerr = ebDlipivot[label], width = bar_width,\
                      label='$D_{Li^{+}}$, $f_{an}$ = ' + str(label),color=clr_arr[i],capsize=5)


    for i, label in enumerate(stfsipivot.columns):
        bars = ax.bar(x + i * bar_width, stfsipivot[label],yerr = ebDstfsipivot[label],width = bar_width,\
                      label='$D_{STFSI^{-}}$, $f_{an}$ = ' + str(label),color=clr_arr[i],\
                      hatch='//', edgecolor='black', fill=False,capsize=5)

    ax.set_xticks(x + bar_width)
    ax.set_xticklabels(stfsipivot.index)
    ax.set_ylim([0, 0.06])
#    ax.legend(loc='upper center', ncol=3, frameon=False)
#    ax.tick_params(axis='x',pad=40)
    ax.legend()
    fig.savefig(f'{figout_dir}/diff_all_stacked_avg.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/diff_all_stacked_avg.eps',format = 'eps')

    
    #--------------Plotting ion-diffusivity------------------
    print("Plotting ion-diffusivity")    
    fig, ax = plt.subplots()
    paux.set_axes(ax,plt,r'Model',r'$D_{Li^{+}}$ (${\sigma^2}/{\tau}$)')

    for i, label in enumerate(lipivot.columns):
        bars = ax.bar(x + i * bar_width, lipivot[label], width = bar_width,\
                      label='$f_{an}$ = ' + str(label),color=clr_arr[i])
        
        # Annotate each bar with the corresponding f_an value
        # for j, bar in enumerate(bars):
        #     f_an_val = lipivot.columns[i]
        #     height = bar.get_height()
        #    ax.text(bar.get_x() + bar.get_width()/2, height,
        #            f"{f_an_val:.2f}", ha='center', va='bottom', fontsize=9)
        
    # Labels, legend, and formatting
    ax.set_xticks(x + bar_width)
    ax.set_xticklabels(lipivot.index)
    ax.tick_params(axis='x',pad=20)
    ax.legend()
    fig.savefig(f'{figout_dir}/diff_ion_avg.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/diff_ion_avg.eps',format = 'eps')

    #--------Plotting counter-ion-diffusivity------------------
    print("Plotting counter-ion diffusivity")    
    fig, ax = plt.subplots()
    paux.set_axes(ax,plt,r'Model',r'$D_{STFSI^{-}}$ (${\sigma^2}/{\tau}$)')
    for i, label in enumerate(stfsipivot.columns):
        bars = ax.bar(x + i * bar_width, stfsipivot[label],width = bar_width,\
                      label='$f_{an}$ = ' + str(label),color=clr_arr[i])

        # Annotate each bar with the corresponding f_an value
        # for j, bar in enumerate(bars):
        #     f_an_val = stfsipivot.columns[i]
        #     height = bar.get_height()
        #     ax.text(bar.get_x() + bar.get_width()/2, height,
        #             f"{f_an_val:.2f}", ha='center', va='bottom', fontsize=9)

        
    # Labels, legend, and formatting
    ax.set_xticks(x + bar_width)
    ax.set_xticklabels(stfsipivot.index)
    ax.set_yscale('log')
    ax.legend()
    fig.savefig(f'{figout_dir}/diff_countion_avg.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/diff_countion_avg.eps',format = 'eps')

    #--------Plotting Transference Number (z = 1) --------------------
    print("Plotting transference number (z=1)")    
    fig, ax = plt.subplots()
    paux.set_axes(ax,plt,r'Model',r'$t_{+}$')
    for i, label in enumerate(tz1lipivot.columns):
        bars = ax.bar(x + i * bar_width, tz1lipivot[label],width = bar_width,\
                      label='$f_{an}$ = ' + str(label),color=clr_arr[i])

        # Annotate each bar with the corresponding f_an value
        #for j, bar in enumerate(bars):
        #    f_an_val = tz1lipivot.columns[i]
        #    height = bar.get_height()
        #    ax.text(bar.get_x() + bar.get_width()/2, height,
        #            f"{f_an_val:.2f}", ha='center', va='bottom', fontsize=9)
        
    # Labels, legend, and formatting
    ax.set_ylim([0.6, 1.05])
    ax.set_xticks(x + bar_width)
    ax.set_xticklabels(tz1lipivot.index)
    ax.legend()
    fig.savefig(f'{figout_dir}/tplus_avg.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/tplus_avg.eps',format = 'eps')

    #--------Plotting Transference Number (z = fN) --------------------
    print("Plotting transference number (z = fN)")    
    fig, ax = plt.subplots()
    paux.set_axes(ax,plt,r'Model',r'$t_{+}$')

    for i, label in enumerate(tzfNlipivot.columns):
        ax.bar(x + i * bar_width, tzfNlipivot[label], yerr = ebtzfNpivot[label],width = bar_width,\
               label='$f_{an}$ = ' + str(label),color=clr_arr[i],capsize=5)

        
    # Labels, legend, and formatting
    ax.set_ylim([0.3, 1.0])
    ax.set_xticks(x + bar_width)
    ax.set_xticklabels(tzfNlipivot.index)
    ax.legend()
    fig.savefig(f'{figout_dir}/tplus_scaledfN_avg.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/tplus_scaledfN_avg.eps',format = 'eps')

    
    
if(exppoly_flag):

    print ("Plotting experimental (un)polymerization data")
    
    # Bar settings
    fanfrac_bar = np.arange(len(fanfrac_labels))  # the label locations
    bar_width = 0.4

    #----------Unpolymerized fraction-----------------------------------
    print("Unpolymerized polymer fraction data")

    #Categories and values
    unpoly_values = [0.707, 0.609, 0.884]
        
    fig, ax = plt.subplots()
    paux.set_axes(ax,plt,r'$f_{an}$',r'Unpolymerized Fraction')

    ax.bar(fanfrac_bar, unpoly_values, bar_width,color='skyblue',edgecolor='black')
        
    # Labels, legend, and formatting
    ax.set_ylim([0.5, 1])
    ax.set_xticks(fanfrac_bar)
    ax.set_xticklabels(fanfrac_labels)
    fig.savefig(f'{figout_dir}/unpoly_expdata.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/unpoly_expdata.eps',format = 'eps')

if(expdiff_flag):

    print ("Plotting experimental diffusivity and transference number data")
   
    # Bar settings
    fanfrac_bar = np.arange(len(fanfrac_labels))  # the label locations
    bar_width = 0.4
    
    #----------Ion diffusivity-----------------------------------
    print("Ion diffusivity data")
    #Categories and values
    dliexp_values = [4.51E-11, 1.02E-11, 3.29E-12]
        
    fig, ax = plt.subplots()
    paux.set_axes(ax,plt,r'$f_{an}$',r'$D_{Li^{+}}$ (m$^2$/s)')
    ax.bar(fanfrac_bar, dliexp_values, bar_width,color='skyblue',edgecolor='black')
        
    # Labels, legend, and formatting
    ax.set_xticks(fanfrac_bar)
    ax.set_xticklabels(fanfrac_labels)
    ax.set_yscale('log')
    fig.savefig(f'{figout_dir}/dli_expdata.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/dli_expdata.eps',format = 'eps')

    #----------Counterion diffusivity-----------------------------------
    print("Counterion unaveraged diffusivity data")
    #Categories and values
    dstfsiexp_values = [1.55E-11,3.60E-13,3.44E-13]
    
    fig, ax = plt.subplots()
    paux.set_axes(ax,plt,r'$f_{an}$',r'$D_{STFSI^{-}}$ (m$^2$/s)')
    ax.bar(fanfrac_bar, dstfsiexp_values, bar_width,color='skyblue',edgecolor='black')
        
    # Labels, legend, and formatting
    ax.set_xticks(fanfrac_bar)
    ax.set_xticklabels(fanfrac_labels)
    ax.set_yscale('log')
    fig.savefig(f'{figout_dir}/dstfsi_expdata.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/dstfsi_expdata.eps',format = 'eps')

    #----------Avged Counterion diffusivity-----------------------------------
    print("Counterion AVERAGED diffusivity data")
    #Categories and values
    
    f1 = np.array([0.933, 0.917, 0.854])
    f2 = 1-f1
    
    d1 = np.array([1.55E-11,3.60E-13,3.44E-13])
    d2 = np.array([1.04E-12,1.64E-11,1.05E-11])
    
    dstfsiavgexp_values = f1*d1 + f2*d2
    
    fig, ax = plt.subplots()
    paux.set_axes(ax,plt,r'$f_{an}$',r'$D_{STFSI^{-}}$ (m$^2$/s)')
    ax.bar(np.array(fanfrac_bar), dstfsiavgexp_values, bar_width,color='skyblue',edgecolor='black')
        
    # Labels, legend, and formatting
    ax.set_xticks(fanfrac_bar)
    ax.set_xticklabels(fanfrac_labels)
    ax.set_yscale('log')
    fig.savefig(f'{figout_dir}/avgdstfsi_expdata.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/avgdstfsi_expdata.eps',format = 'eps')

    #-----Combined diffusivities (lithium and averaged anion diffusivities)--------------

    fig, ax = plt.subplots()
    paux.set_axes(ax,plt,r'$f_{an}$',r'Diffusion Coefficients (m$^2$/s)')

    # Bars
    bwidth = 0.35
    rects1 = ax.bar(fanfrac_bar - bwidth/2, dliexp_values, bwidth, label="D$_{Li^{+}}$")
    rects2 = ax.bar(fanfrac_bar + bwidth/2, dstfsiavgexp_values, bwidth, label="D$_{STFSI^{-}}$")
    
    # Labels, legend, and formatting
    ax.set_xticks(fanfrac_bar)
    ax.set_xticklabels(fanfrac_labels)
    ax.set_yscale('log')
    ax.legend()
    ax.set_yticks([5e-12, 1e-11, 5e-11])
    ax.yaxis.set_major_formatter(FuncFormatter(paux.sci_mathtext))
    fig.savefig(f'{figout_dir}/combined_diff_expdata.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/combined_diff_expdata.eps',format = 'eps')

    
    #----------Transference number (z=1)-----------------------------------
    print("Transference number (z=1) data")
    #Categories and values
    tz1exp_values = [0.756, 0.858, 0.644]
      
    fig, ax = plt.subplots()
    paux.set_axes(ax,plt,r'$f_{an}$',r'$t_{+}$')

    ax.bar(fanfrac_bar, tz1exp_values, bar_width,color='skyblue',edgecolor='black')
        
    # Labels, legend, and formatting
    ax.set_ylim([0.6, 1])
    ax.set_xticks(fanfrac_bar)
    ax.set_xticklabels(fanfrac_labels)
    fig.savefig(f'{figout_dir}/tliz1exp_data.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/tliz1exp_data.eps',format = 'eps')

    #----------Avg Transference number (z=fN with N=100)-----------------------------------
    print("Transference number (z=fN) data")
    #Categories and values
    Ndp_exp = 100
    fNdp = np.array(fanfrac_values)*Ndp_exp
    tzfNexpavg_values = np.array(dliexp_values)/(np.array(dliexp_values) + fNdp*dstfsiavgexp_values)
    fig, ax = plt.subplots()
    paux.set_axes(ax,plt,r'$f_{an}$',r'$t_{+}$')

    ax.bar(fanfrac_bar, tzfNexpavg_values, bar_width,color='skyblue',edgecolor='black')
        
    # Labels, legend, and formatting
#    ax.set_ylim([0.6, 1])
    ax.set_xticks(fanfrac_bar)
    ax.set_xticklabels(fanfrac_labels)
    fig.savefig(f'{figout_dir}/tlizfNexpavg_data.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/tlizfNexpavg_data.eps',format = 'eps')

if(exptnum_Ndep):

    print ("Plotting experimental transference number as a function of N")

    # Nvalues
    Ndeg_poly = np.array([20,40,75,100,500,1000])

    # Bar settings
    N_barpos  = np.arange(len(Ndeg_poly))  # the label locations
    n_series  = len(fanfrac_values)
    bar_width = 0.8/n_series
    
    #Diffusivity values
    dliexp_values = np.array([4.51E-11, 1.02E-11, 3.29E-12])


    f1 = np.array([0.933, 0.917, 0.854])
    f2 = 1-f1    
    d1 = np.array([1.55E-11,3.60E-13,3.44E-13])
    d2 = np.array([1.04E-12,1.64E-11,1.05E-11])
    dstfsiavgexp_values = f1*d1 + f2*d2

    print(dliexp_values)
    print(dstfsiavgexp_values)
    
    # Transference number as a function of N
    fNalldp = Ndeg_poly.T.reshape(-1,1)*np.array(fanfrac_values).reshape(1,-1) # Transpose and convert to 2D for broadcasting
    tzfNexp_values_set = dliexp_values/(dliexp_values + fNalldp*dstfsiavgexp_values)
    
    print(tzfNexp_values_set)
    
    fig, ax = plt.subplots()
    for j in range(len(fanfrac_values)):
        ax.bar(N_barpos + (j - (n_series-1)/2)*bar_width, tzfNexp_values_set[:,j], bar_width,
               label='$f_{an}$ = ' + fanfrac_labels[j])
    
    # Labels, legend, and formatting
    paux.set_axes(ax,plt,r'$N_{poly}$',r'$t_{+}$')
    ax.set_xticks(N_barpos)
    ax.set_xticklabels([str(val) for val in Ndeg_poly.T])
    ax.legend()
    fig.savefig(f'{figout_dir}/exptnum_Ndep.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/exptnum_Ndep.eps',format = 'eps')
    

if(expcond_flag):

    print ("Plotting experimental conductivity data")
    
    # Bar settings
    fanfrac_bar = np.arange(len(fanfrac_labels))  # the label locations
    bar_width = 0.4
        
    #----------Ion conductivity-----------------------------------
    # Plot data
    fig, ax = plt.subplots()
    bar_width = 0.2
    
    all_data = np.zeros((maxneigh,len(fran_arr)+1))
    expcond_data = expinp_dir + '/expcond_data_jy.xlsx'

    heating = pd.read_excel(expcond_data, sheet_name="Alldata_heating")
    cooling = pd.read_excel(expcond_data, sheet_name="Alldata_cooling")

    # clean column names (strip spaces/newlines)
    heating.columns = [str(c).strip() for c in heating.columns]
    cooling.columns = [str(c).strip() for c in cooling.columns]

    # X-axis: 1000/T from Temp_degK
    for df in (heating, cooling):
        df["Temp_degK"] = pd.to_numeric(df.get("Temp_degK", np.nan), errors="coerce")
        df["x"] = 1000.0 / df["Temp_degK"]
        # Coerce series to numeric if present
        for c in ["5to1_heating", "10to1_heating", "20to1_heating",
                  "5to1_cooling", "10to1_cooling", "20to1_cooling"]:
            if c in df.columns:
                df[c] = pd.to_numeric(df[c], errors="coerce")


    # Map columns to legend labels
    heating_map = {
        "5to1_heating":  "f$_{an}$ = 0.167 (heating)",
        "10to1_heating": "f$_{an}$ = 0.091 (heating)",
        "20to1_heating": "f$_{an}$ = 0.047 (heating)",
    }
    cooling_map = {
        "5to1_cooling":  "f$_{an}$ = 0.167 (cooling)",
        "10to1_cooling": "f$_{an}$ = 0.091 (cooling)",
        "20to1_cooling": "f$_{an}$ = 0.047 (cooling)",
    }


    # Heating (solid, circles)
    fr_id = 0
    for col, label in heating_map.items():
        if col in heating.columns:
            plt.plot(heating["x"].to_numpy(), heating[col].to_numpy(), marker="o", linestyle="--", \
                     color = clr_arr[fr_id], markerfacecolor = clr_arr[fr_id],markersize=10,label=label)
            fr_id +=1
            
    # Cooling (dashed, squares)
    fr_id = 0
    for col, label in cooling_map.items():
        if col in cooling.columns:
            plt.plot(cooling["x"].to_numpy(), cooling[col].to_numpy(), marker="s", linestyle="--", \
                     color = clr_arr[fr_id],markersize=10,label=label)
            fr_id +=1

    paux.set_axes(ax,plt,r'1000/T (K$^{-1}$)',r'$\sigma$ (S/cm)')
    ax.set_yscale('log')
    ax.legend()
    fig.savefig(f'{figout_dir}/expcond_data.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/expcond_data.eps',format = 'eps')


if(s_clus_flag):

    print("Plotting cumulative cluster sums...")
    df = pd.read_excel('../../analyzed_results/avg_cluster_S2_f_0p2.xlsx',sheet_name='Cluster_details')

    fig, ax = plt.subplots()
    paux.set_axes(ax,plt,r'Cluster size, $s$',r'Average number of ions, ($sH(s)$)')
    plt.scatter(df["s"],df["sH_s"],color = dark_g)
    ax.set_xlim([1350, 1600])
    fig.savefig(f'{figout_dir}/sHs_data.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/sHs_data.eps',format = 'eps')

    fig, ax = plt.subplots()
    paux.set_axes(ax,plt,r'Cluster size, $s$',r'Cumulative sum, ($\Sigma_{i=s}^{s_{max}} i H(i)$)')
    plt.scatter(df["s"],df["Rev_Cum_Sum"],color = dark_g)
    ax.set_xlim([1350, 1600])
    ax.set_ylim([0, 1600])
    fig.savefig(f'{figout_dir}/rev_cum_sum_clust.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/rev_cum_sum_clust.eps',format = 'eps')


if(expnmr_flag):

    sheet_pref = ['0p05','0p11','0p25'];
    leg_vals   = ['0.047','0.091','0.167'];
    
    print("Plotting Li-NMR data...")

    for sid,sname in enumerate(sheet_pref):

        print("Plotting " + sname)
        
        df = pd.read_excel('../../analyzed_results/Li_NMRdata.xlsx',sheet_name='Li_'+sname)
        
        fig, ax = plt.subplots()
        paux.set_axes(ax,plt,r'Gradient (s/m$^{2}$)',r'Area')
        ax.set_title(r'$f_{an}$ = ' +  leg_vals[sid])
        plt.plot(df["Gradient"].to_numpy(),df["Experiment"].to_numpy(),\
                 color='red',marker='o',markerfacecolor='red',ls='None',\
                 markersize=8,label = 'Experimental data')
        plt.plot(df["Gradient"].to_numpy().reshape(-1,1),df["Fit"].to_numpy(),\
                 color='black',linestyle='-',linewidth=2,label = 'Fit')

        ax.set_yscale('log')
        ax.legend()
        fig.savefig(f'{figout_dir}/Li_{sname}.png',dpi = fig.dpi)
        fig.savefig(f'{figout_dir}/Li_{sname}.eps',format = 'eps')

        df = pd.read_excel('../../analyzed_results/F_NMRdata.xlsx',sheet_name='F_'+sname)

        fig, ax = plt.subplots()
        paux.set_axes(ax,plt,r'Gradient (s/m$^{2}$)',r'Area')
        ax.set_title(r'$f_{an}$ = ' +  leg_vals[sid])
        plt.plot(df["Gradient"].to_numpy(),df["Experiment"].to_numpy(),\
                 color='red',marker='o',markerfacecolor='red',ls='None',\
                 markersize=8,label = 'Experimental data')
        plt.plot(df["Gradient"].to_numpy().reshape(-1,1),df["Fit"].to_numpy(),\
                 color='black',linestyle='-',linewidth=2,label = 'Fit')

        ax.set_yscale('log')
        ax.legend()
        fig.savefig(f'{figout_dir}/F_{sname}.png',dpi = fig.dpi)
        fig.savefig(f'{figout_dir}/F_{sname}.eps',format = 'eps')


if(n_surf_flag):

    neigh_vals = ['0p05','0p11','0p25'];
    leg_vals   = ['0.05','0.09','0.16'];
    step_val   = 2500
    time_step  = [0.006, 0.006, 0.006]
    
    print("Plotting Li cluster surface plot...")

    for sid,sname in enumerate(neigh_vals):

        print("Plotting " + sname)
        surf_dir = '../../new_results/qanion_0.2/surf_plots'
        fname = surf_dir + '/clustime_' + sname + '.txt'
        if not os.path.exists(fname):
            print("ERR: " + fname + " does not exist")
            continue
        
        times,x,Y = paux.compute_neigharr(fname,step_val,time_step[sid])
        Tg, Xg = np.meshgrid(times,x)
        
        fig, ax = plt.subplots()
        #plt.figure(figsize=(9, 6))
        pcm = ax.pcolormesh(Tg,Xg, Y.T/100, shading="auto") # Convert to fraction
        
        plt.ylabel(r'Anion neighbors, $s$')
        plt.xlabel(r'Time ($\tau$)')
        cbar = plt.colorbar(pcm)
        cbar.set_label(r'$f_{\mathrm{neigh}}(s)$')
        plt.ylim(-0.5,6.5)
        plt.tight_layout()

        fig.savefig(f'{figout_dir}/neighsurf_{sname}.png',dpi = fig.dpi)
        fig.savefig(f'{figout_dir}/neighsurf_{sname}.eps',format = 'eps')

if(n_surf_flag):

    neigh_vals = ['0p05','0p11','0p25'];
    leg_vals   = ['0.05','0.09','0.16'];
    step_val   = 2500
    time_step  = [0.006, 0.006, 0.006]
    
    print("Plotting Li neighbor surface plot...")

    for sid,sname in enumerate(neigh_vals):

        print("Plotting " + sname)
        surf_dir = '../../new_results/qanion_0.2/surf_plots'
        fname = surf_dir + '/catanneightime_' + sname + '.txt'
        if not os.path.exists(fname):
            print("ERR: " + fname + " does not exist")
            continue
        
        times,x,Y = paux.compute_neigharr(fname,step_val,time_step[sid])
        Tg, Xg = np.meshgrid(times,x)
        
        fig, ax = plt.subplots()
        #plt.figure(figsize=(9, 6))
        pcm = ax.pcolormesh(Tg,Xg, Y.T/100, shading="auto") # Convert to fraction
        
        plt.ylabel(r'Anion neighbors, $s$')
        plt.xlabel(r'Time ($\tau$)')
        cbar = plt.colorbar(pcm)
        cbar.set_label(r'$f_{\mathrm{neigh}}(s)$')
        plt.ylim(-0.5,6.5)
        plt.tight_layout()

        fig.savefig(f'{figout_dir}/neighsurf_{sname}.png',dpi = fig.dpi)
        fig.savefig(f'{figout_dir}/neighsurf_{sname}.eps',format = 'eps')

if(tg_flag):

    print("Plotting Tg data...")

    tg_dir = '../../analyzed_results'
    fname = tg_dir + '/tgVals_S3_0p5_15.xlsx'
    
    if not os.path.exists(fname):
        raise RuntimeError("FATAL ERROR: " + fname + " not found")
        
    df = pd.read_excel(fname, sheet_name="AvgData")

    t_x   = df['Avg_T'].to_numpy()
    rho_y = df['Avg_rho'].to_numpy()

    
    idx = np.argsort(t_x) # Sort by temp
    T = t_x[idx]
    rho = rho_y[idx]


    # Compute best fit slopes and ssr
    i_star, p1, p2, ssr = paux.best_two_line_fit(T, rho)
    
    a1, b1 = p1  # glassy side
    a2, b2 = p2  # liquid side

    Tg = (b2 - b1) / (a1 - a2)
    rho_g = a1 * Tg + b1

    print(f"Tg: {Tg:.4f}, rho(Tg) {rho_g:.4f}")


    #Plot data
    fig, ax = plt.subplots()
    T1_line = np.linspace(T[0], T[i_star-1] + 0.1, 50)
    T2_line = np.linspace(T[i_star] - 0.1, T[-1], 50)

    rho1_line = np.polyval(p1, T1_line)
    rho2_line = np.polyval(p2, T2_line)

    # Data
    plt.scatter(T, rho, label="Data", zorder=3, marker='d',\
                s = 30, color = orange)

    # Fits
    plt.plot(T1_line, rho1_line, label="Glassy fit",\
             linewidth=2,linestyle='--',color = 'red')
    plt.plot(T2_line, rho2_line, label="Liquid fit",\
             linewidth=2,linestyle='--',color = dark_g)

    # Tg
    plt.axvline(Tg, linestyle="--", linewidth=1.5,
                label=None,color=violet)
    plt.scatter([Tg], [rho_g], marker="o", s=80, zorder=4, \
                label=fr"$T_g \approx {Tg:.3f}$",color = violet)
    

    plt.ylabel(r'Reduced number density, $\rho^{*}$')
    plt.xlabel(r'Reduced temperature, $T^{*}$')
    plt.tight_layout()
    ax.legend()
    
    fig.savefig(f'{figout_dir}/tg_sim.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/tg_sim.eps',format = 'eps')

if(cross_flag):

    print("Plotting Cross data...")
    delt = 0.005
    msd_dir = '../../new_results/qanion_0.2/cross_diff/'
    fname = glob.glob(msd_dir + 'iondiff*')   
    if fname == []:
        raise RuntimeError("FATAL ERROR: " + fname + " not found")
    ion_df = pd.read_csv(fname[0], sep=" ",skiprows = 1,skipinitialspace=True)
    
    fname = glob.glob(msd_dir + 'countiondiff*')
    if fname == []:
        raise RuntimeError("FATAL ERROR: " + fname + " not found")
    countion_df = pd.read_csv(fname[0], sep=" ",skiprows = 1,skipinitialspace=True)    

    fname = glob.glob(msd_dir + 'crossdiff*')
    if fname == []:
        raise RuntimeError("FATAL ERROR: " + fname + " not found")
    cross_df = pd.read_csv(fname[0], sep=" ",skiprows = 1,skipinitialspace=True)

    #Plot data
    fig, ax = plt.subplots()
    t1 = ion_df.iloc[0,0]
    plt.plot(delt*(ion_df.iloc[:,0]-t1).to_numpy(), ion_df.iloc[:,1:4].sum(axis=1).to_numpy(),\
             label=r"$r^{2}_{++}$",linewidth=2,linestyle='--',color = 'red')
    
    plt.plot(delt*(countion_df.iloc[:,0]-t1).to_numpy(), countion_df.iloc[:,1:4].sum(axis=1).to_numpy(),\
             label=r"$r^2_{--}$",linewidth=2,linestyle='--',color = violet)

    plt.plot(delt*(cross_df.iloc[:,0]-t1).to_numpy(), cross_df.iloc[:,1:4].sum(axis=1).to_numpy(),\
             label=r"$r^{2}_{+-}$",linewidth=2,linestyle='--',color = dark_g)


    ax.set_xlim([0,12000])
    plt.ylabel(r'Mean Squared Displacements ($\sigma^{2}$)')
    plt.xlabel(r'Time ($\tau$)')
    plt.tight_layout()
    ax.legend()
    
    fig.savefig(f'{figout_dir}/crossmsd_all.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/crossmsd_all.eps',format = 'eps')

    #Plot data
    fig, ax = plt.subplots()
    plt.plot(delt*(cross_df.iloc[:,0]-t1).to_numpy(), cross_df.iloc[:,1:4].sum(axis=1).to_numpy(),\
             label=r"$r^{2}_{+-}$",linewidth=2,linestyle='--',color = dark_g)

    ax.set_xlim([0,10000])
    plt.ylabel(r'Mean Squared Displacements ($\sigma^{2}$)')
    plt.xlabel(r'Time ($\tau$)')
    plt.tight_layout()
    ax.legend()
    
    fig.savefig(f'{figout_dir}/crossmsd_zoom.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/crossmsd_zoom.eps',format = 'eps')

        
# All simulation analysis    
for sysid,sysname in enumerate(sys_arr):
    
    sys_dir = all_dir + '/' + syspref + q_anion
    if not os.path.isdir(sys_dir):
        raise RuntimeError("FATAL ERROR: " + sys_dir + " not found")

    for casenum in case_arr:
        casestr = str(casenum)

        if casestr == 'S2':
            frho_an  = [nanions[i]/S2_vol_tot[i] for i in range(len(totbeads))]
            frho_pol = [npol_ans[i]/S2_vol_tot[i] for i in range(len(totbeads))]
        else:
            frho_an  = [nanions[i]/S3_vol_tot[i] for i in range(len(totbeads))]
            frho_pol = [npol_ans[i]/S3_vol_tot[i] for i in range(len(totbeads))]

        if(rdf_flag):

            
            for cid, (rdfkey,colval) in enumerate(rdf_keys.items()):

                # Plot data
                fig, ax = plt.subplots()
                paux.set_axes(ax,plt,r'$r$ ($\sigma$)',r'g(r), n(r)')
                maxgr = 2
                for fr_id, frac_an in enumerate(fran_arr):

                    if not 'Poly' in rdfkey:
                        rhof = frho_an[fr_id]
                    else:
                        rhof = frho_pol[fr_id]*(1-f_unpoly[sysid])

                    result_dir = sys_dir + '/' + respref + sysname + '_' + frac_an + \
                        '_' + str(casenum)

                    if not os.path.isdir(result_dir):
                        raise RuntimeError("FATAL ERROR: " + result_dir + " not found")

                    print (f'RDF analysis for {rdfkey} in {result_dir}')

                    rarr,gofr,nofr = ardf.ana_rdf(result_dir,rdf_pref,colval-1)
                    plt.plot(rarr, gofr, color = clr_arr[fr_id], linestyle = '--', \
                             label = '$f_{an}$ = ' + str(frac_an))
                    plt.plot(rarr, rhof*nofr, color = clr_arr[fr_id], linestyle='-')
                    maxgrt = max(gofr)
                    if maxgrt > maxgr: maxgr = maxgrt

                
                ax.set_xlim([0, 3.5])
                ax.set_ylim([0, maxgr+1])
                ax.legend()
                fig.savefig(f'{figout_dir}/rdf_{sysname}_{rdfkey}_{casestr}.png',dpi = fig.dpi)
                fig.savefig(f'{figout_dir}/rdf_{sysname}_{rdfkey}_{casestr}.eps',format = 'eps')

        if(neigh_flag):
            
            for nid, (neighkey,colval) in enumerate(neigh_keys.items()):

                # Plot data
                fig, ax = plt.subplots()
                paux.set_axes(ax,plt,r'Anion neighbor size, $s$',r'$f_{\mathrm{neigh}}(s)$')
                bar_width = 0.2
                
                all_data = np.zeros((maxneigh,len(fran_arr)+1))
                for fr_id, frac_an in enumerate(fran_arr):
                    
                    result_dir = sys_dir + '/' + respref + sysname + '_' + frac_an + \
                        '_' + str(casenum)

                    if not os.path.isdir(result_dir):
                        raise RuntimeError("FATAL ERROR: " + result_dir + " not found")
                    
                    print (f'Neighbor analysis for {neighkey} in {result_dir}')

                    xarr,yarr = paux.return_neigh_arrays(result_dir,\
                                                         neigh_pref,[0,colval-1],\
                                                         0,maxneigh)
                    if fr_id == 0:
                        all_data[:,fr_id] = xarr - 1
                    all_data[:,fr_id+1]   = yarr

                for i in range(len(fran_arr)):
                    ax.bar(all_data[:,0]+i*bar_width,all_data[:,i+1]/100,\
                           bar_width,label='$f_{an}$ = ' + fran_arr[i],\
                           color=clr_arr[i])


                ax.set_xlim([-0.2, maxneigh+1])
                xtick_centers = [r + bar_width / 2 for r in range(maxneigh+1)]
                plt.xticks(xtick_centers, np.arange(0,maxneigh+1,1))
                ax.legend()
                fig.savefig(f'{figout_dir}/neigh_{sysname}_{neighkey}_{casestr}.png',dpi = fig.dpi)
                fig.savefig(f'{figout_dir}/neigh_{sysname}_{neighkey}_{casestr}.eps',format = 'eps')

        if(clust_flag):
            
            for nid, (clustkey,colval) in enumerate(clust_keys.items()):

                bar_width = 0.2
                maxclust = 1600
                
                # Plot data
                fig, ax = plt.subplots()
                paux.set_axes(ax,plt,r'Cluster size, $s$',r'$H$($s$)')

                fig2, ax2 = plt.subplots()
                paux.set_axes(ax,plt,r'Cluster size, $s$',r'$H$($s$)')
                
                for fr_id, frac_an in enumerate(fran_arr):

                    result_dir = sys_dir + '/' + respref + sysname + '_' + frac_an + \
                        '_' + str(casenum)

                    if not os.path.isdir(result_dir):
                        raise RuntimeError("FATAL ERROR: " + result_dir + " not found")

                    print (f'Cluster analysis for {clustkey} in {result_dir}')

                    xarr,yarr = paux.return_clust_arrays(result_dir,\
                                                         clust_pref,[0,colval-1],\
                                                         0,0)

                    ax.bar(xarr+fr_id*bar_width,yarr,bar_width,\
                           label='$f_{an}$ = ' + fran_arr[fr_id],color=clr_arr[fr_id])
                    ax2.bar(xarr+fr_id*bar_width,yarr,bar_width,\
                            label='$f_{an}$ = ' + fran_arr[fr_id],color=clr_arr[fr_id])
                    
                
                xtick_centers = [r + bar_width / 2 for r in range(1,10,1)] 
                ax.set_xlim([0.5, 10])
                plt.xticks(xtick_centers, np.arange(1,10,1))
                ax.legend()
                fig.savefig(f'{figout_dir}/clust_{sysname}_{clustkey}_{casestr}.png',dpi = fig.dpi)
                fig.savefig(f'{figout_dir}/clust_{sysname}_{clustkey}_{casestr}.eps',format = 'eps')

                xtick_centers = [r + bar_width / 2 for r in range(1300,maxclust+1,100)]
                ax2.set_xlim([1300, maxclust+5])
                ax2.set_ylim([0, 0.02])
                plt.xticks(xtick_centers, np.arange(1300,maxclust+1,100))
                ax2.legend()
                fig2.savefig(f'{figout_dir}/clustinset_{sysname}_{clustkey}_{casestr}.png',dpi = fig.dpi)
                fig2.savefig(f'{figout_dir}/clustinset_{sysname}_{clustkey}_{casestr}.eps',format = 'eps')
                
        if(restime_flag):
            
            # Plot data
            fig, ax = plt.subplots()
            paux.set_axes(ax,plt,r'Time ($\tau$)',r'$TACF(t)$')

            for fr_id, frac_an in enumerate(fran_arr):
                
                result_dir = sys_dir + '/' + respref + sysname + '_' + frac_an + \
                    '_' + str(casenum)
                
                if not os.path.isdir(result_dir):
                    raise RuntimeError("FATAL ERROR: " + result_dir + " not found")
                
                print (f'Autocorrelation analysis for {reskey} in {result_dir}')
                    
                tarr,cfun,tau,beta,meantau = atcf.ana_tacf(result_dir,tacf_pref,colval-1)
                plt.plot(tarr, cfun, color = clr_arr[fr_id], linestyle = '--', \
                         label = '$f_{an}$ = ' + str(frac_an))

            meantau_arr[sysid,fr_id] = meantau
                    
            ax.set_xlim([0.001, 20])
            ax.set_ylim([0, 1])
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.legend()
            fig.savefig(f'{figout_dir}/tacf_{sysname}_{rdfkey}.png',dpi = fig.dpi)
            fig.savefig(f'{figout_dir}/tacf_{sysname}_{rdfkey}.eps',format = 'eps')
            fig.close()

        if(n_sens_flag):

            rcut_vals = ['15','16','17','18','19'] # corresponds to TEN-TIMES cutoff distance
            rcut_leg  = ['1.5','1.6','1.7','1.8','1.9'] # correct r-cut legends
            
            for nid, (neighkey,colval) in enumerate(neigh_keys.items()):

                # Plot data
                fig, axes = plt.subplots(2,3,figsize=(12,6))
                fig.supxlabel(r'Anion neighbor size, $s$', fontsize=16)
                fig.supylabel(r'$f_{\mathrm{neigh}}(s)$', fontsize=16)
                bar_width = 0.2

                axes = axes.ravel()
                
                for rind,rcut in enumerate(rcut_vals):

                    print (f'Neighbor analysis for rcut = {rcut}')
                    
                    all_data = np.zeros((maxneigh,len(fran_arr)+1))
                    sensfile_pref = neigh_sens_pref + rcut
                    ax = axes[rind]
                    ax.set_title(r'$r_{c}$ = ' + rcut_leg[rind] + ' $\sigma$')
                    
                    for fr_id, frac_an in enumerate(fran_arr):
                    
                        result_dir = sys_dir + '/' + senspref + sysname + '_' + frac_an + \
                            '_' + str(casenum)
                        
                        if not os.path.isdir(result_dir):
                            raise RuntimeError("FATAL ERROR: " + result_dir + " not found")
                        
                        print (f'{neighkey} in {result_dir}')

                        
                        xarr,yarr = paux.return_neigh_arrays(result_dir,\
                                                         sensfile_pref,[0,colval-1],\
                                                         0,maxneigh)
                        if fr_id == 0:
                            all_data[:,fr_id] = xarr - 1
                        all_data[:,fr_id+1]   = yarr

                    
                    for i in range(len(fran_arr)):
                        ax.plot(all_data[:,0],all_data[:,i+1]/100,\
                                label='$f_{an}$ = ' + fran_arr[i],\
                                color=clr_arr[rind],marker=mrk_arr[i],\
                                markersize=8, markerfacecolor=clr_arr[rind],\
                                linestyle='--',linewidth=2)
                    ax.set_xlim([-0.2, maxneigh+1])                            
                    ax.legend()

                axes[5].set_visible(False)
                fig.savefig(f'{figout_dir}/sens_neigh_{sysname}_{neighkey}_{casestr}.png',dpi = fig.dpi)
                fig.savefig(f'{figout_dir}/sens_neigh_{sysname}_{neighkey}_{casestr}.eps',format = 'eps')

        if(c_sens_flag):

            rcut_vals = ['15','16','17','18','19'] # corresponds to TEN-TIMES cutoff distance
            rcut_leg  = ['1.5','1.6','1.7','1.8','1.9'] # correct r-cut legends
            
            for nid, (clustkey,colval) in enumerate(clust_keys.items()):

                # Plot data
                fig, axes = plt.subplots(2,3,figsize=(12,6))
                fig.supxlabel(r'Cluster size, $s$', fontsize=16)
                fig.supylabel(r'$H$($s$)', fontsize=16)
                bar_width = 0.2

                axes = axes.ravel()
                
                for rind,rcut in enumerate(rcut_vals):

                    print (f'Cluster analysis for rcut = {rcut}')
                    
                    sensfile_pref = clust_sens_pref + rcut
                    ax = axes[rind]
                    ax.set_title(r'$r_{c}$ = ' + rcut_leg[rind] + ' $\sigma$')
                    maxclust = 1600
                    
                    for fr_id, frac_an in enumerate(fran_arr):
                    
                        result_dir = sys_dir + '/' + senspref + sysname + '_' + frac_an + \
                            '_' + str(casenum)
                        
                        if not os.path.isdir(result_dir):
                            raise RuntimeError("FATAL ERROR: " + result_dir + " not found")
                        
                        print (f'{clustkey} in {result_dir}')

                        
                        xarr,yarr = paux.return_clust_arrays(result_dir,\
                                                             sensfile_pref,[0,colval-1],\
                                                             0,0)
                                        
                        ax.plot(xarr,yarr,label='$f_{an}$ = ' + fran_arr[fr_id],\
                                color=clr_arr[rind],marker=mrk_arr[fr_id],\
                                markersize=8, markerfacecolor=clr_arr[rind],\
                                linestyle='--',linewidth=2)
                    ax.set_xlim([0.5, 10])                            
                    ax.legend()

                axes[5].set_visible(False)
                fig.savefig(f'{figout_dir}/sens_clust_{sysname}_{clustkey}_{casestr}.png',dpi = fig.dpi)
                fig.savefig(f'{figout_dir}/sens_clust_{sysname}_{clustkey}_{casestr}.eps',format = 'eps')
                

        if(rg_flag):

            for nid, (rgkey,colval) in enumerate(rg_keys.items()):

                
                # Plot data
                fig, ax = plt.subplots()
                paux.set_axes(ax,plt,r'Time ($\tau$)',r'R$_{g}$ ($\sigma$)')
                deltvals = np.array([0.005,0.006,0.005])
                models = ['S1','S2','S3']
                leg_vals = ['RCP','RCP+M','HP+M']
                result_dir = '../../new_results/qanion_0.2/rg_data'
                
                if not os.path.isdir(result_dir):
                    raise RuntimeError("FATAL ERROR: " + result_dir + " not found")

                
                for mid,modname in enumerate(models):


                    print (f'Rg analysis for {modname} in {result_dir}')
                    pattern = result_dir + '/rgavg_' + modname + '*'
                    all_files = glob.glob(pattern)
                    
                    if all_files:
                        latest_file = max(all_files, key=os.path.getmtime)
                        print(f"The most recent file matching '{rg_pref}' is: {latest_file}")
                    else:
                        print(f"No files found matching the pattern: {pattern}")

                    df = pd.read_csv(latest_file, sep=" ", skipinitialspace=True)
                    tvals = deltvals[mid]*(df['tstep'] - df['tstep'].iloc[0]).to_numpy()
                    rgall = df['rgt'].to_numpy()
                    
                    plt.plot(tvals, rgall, color = clr_arr[mid], linestyle = '--', \
                             label = 'Model: ' + leg_vals[mid])

                ax.set_xlim([0,20000])
                ax.legend()
                fig.savefig(f'{figout_dir}/rg_all.png',dpi = fig.dpi)
                fig.savefig(f'{figout_dir}/rg_all.eps',format = 'eps')

# Plot all data
# Version: Mar-04-2025
# Import modules
import os
import math
import numpy as np
import analyze_rdfs as ardf
#import analyze_tacf as atcf
#from compute_props import *
import pandas as pd
import aux_plot as paux
from env_plot import *

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
q_anion  = '0.2' # should be the exact string value
sys_arr  = ['S3']
fran_arr = ['0.05','0.09','0.20']
f_unpoly = [0, 0.6, 0.6]
deltat   = [0.006, 0.006, 0.006] 
case_arr = [1]
molwt    = 40
nchains  = 100
rho_tot  = 0.8

#---------File prefixes------------------------------------------
rdf_pref     = 'rdf_config_'
rg_pref      = 'rgavg_config_'
neigh_pref   = 'catanneigh_config_'
restime_pref = 'autocorrcion_config'
diff_pref    = ['iondiff_config_', 'countiondiff_config_']
clust_pref   = 'clust_config_'
bfrdf_pref   = 'freeboundrdf_config_'


#---------Flags---------------------------------------------------
rdf_flag     = 0; rdf_keys    = {'LiAn':3,'LiPoly':2}
rg_flag      = 0
neigh_flag   = 0; neigh_keys  = {'Li-An':3,'An-Li':5}; maxneigh = 6
clust_flag   = 0; clust_keys  = {'Li-An':2}
bfrdf_flag   = 0
restime_flag = 0; restime_keys = {'Li-An':2}
simdiff_flag = 1 # Simulation diffusivities and transference numbers
expdiff_flag = 0 # Experimental diffusivities and transference numbers
expcond_flag = 0 # Experimental conductivities

#---------Experimental inputs---------------------------------------
fanfrac_values = [0.05, 0.11, 0.25]
fanfrac_labels = ['0.05', '0.11', '0.25']
expinp_dir = '../../ExperimentalResults'
if not os.path.isdir(expinp_dir) and expcond_flag == 1:
    print("FATAL ERROR: ", expinp_dir, " not found")
    exit("Check scratch directory path")

#---------Calculations----------------------------------------------
nanions  = [nchains*int(molwt*float(x)) for x in fran_arr]
npols    = [nchains*molwt - nanions[i]  for i in range(len(nanions))]
totbeads = [int(2.0*(npols[i] + nanions[i])) for i in range(len(npols))]
npol_ans = [int(0.5*float(x)) for x in npols]
frho_an  = [rho_tot*(nanions[i]/totbeads[i]) for i in range(len(totbeads))]
frho_pol = [rho_tot*(npol_ans[i]/totbeads[i]) for i in range(len(totbeads))]  

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
        "t_Li_z=1", "t_Li_t=f_an*N"
    ]

    # Load the data into a DataFrame
    df = pd.read_csv('../../analyzed_results/computed_diffusivity_1.dat', sep="\t", skipinitialspace=True)
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
    fig.savefig(f'{figout_dir}/diff_ion.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/diff_ion.eps',format = 'eps')

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
    fig.savefig(f'{figout_dir}/diff_countion.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/diff_countion.eps',format = 'eps')

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
    fig.savefig(f'{figout_dir}/tplus.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/tplus.eps',format = 'eps')

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
    fig.savefig(f'{figout_dir}/tplus_scaledfN.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/tplus_scaledfN.eps',format = 'eps')

if(expdiff_flag):

    print ("Plotting experimental polymerization data")
    
    # Bar settings
    fanfrac_bar = np.arange(len(fanfrac_labels))  # the label locations
    bar_width = 0.4

    #----------Unpolymerized fraction-----------------------------------
    print("Unpolymerized polymer fraction data")

    #Categories and values
    unpoly_values = [70.7, 60.9, 88.4]
        
    fig, ax = plt.subplots()
    paux.set_axes(ax,plt,r'$f_{an}$',r'Unpolymerized Fraction')

    ax.bar(fanfrac_bar, unpoly_values, bar_width,color='skyblue',edgecolor='black')
        
    # Labels, legend, and formatting
    ax.set_ylim([0.6, 1])
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
    print("Counterion counterion diffusivity data")
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
    print("Counterion AVERAGED counterion diffusivity data")
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

    # Combined diffusivities (lithium and averaged anion diffusivities)

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
    fig.savefig(f'{figout_dir}/combined_diff_expdata.png',dpi = fig.dpi)
    fig.savefig(f'{figout_dir}/combined_diff_expdata.eps',format = 'eps')

    
    #----------Transference number (z=1)-----------------------------------
    print("Transference number (z=1) data")
    #Categories and values
    tz1exp_values = [0.744, 0.966, 0.905]
      
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


if(expcond_flag):

    print ("Plotting experimental conductivity data")
    
    # Bar settings
    fanfrac_bar = np.arange(len(fanfrac_labels))  # the label locations
    bar_width = 0.4
        
    #----------Ion conductivity-----------------------------------
    # Plot data
    fig, ax = plt.subplots()
    paux.set_axes(ax,plt,r'$n_{neigh}$',r'$f$($n_{neigh}$)')
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
        "5to1_heating":  "f$_{an}$ = 0.25 (heating)",
        "10to1_heating": "f$_{an}$ = 0.11 (heating)",
        "20to1_heating": "f$_{an}$ = 0.05 (heating)",
    }
    cooling_map = {
        "5to1_cooling":  "f$_{an}$ = 0.25 (cooling)",
        "10to1_cooling": "f$_{an}$ = 0.11 (cooling)",
        "20to1_cooling": "f$_{an}$ = 0.05 (cooling)",
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


    
    
for sysid,sysname in enumerate(sys_arr):
    
    sys_dir = all_dir + '/' + syspref + q_anion
    if not os.path.isdir(sys_dir):
        raise RuntimeError("FATAL ERROR: " + sys_dir + " not found")

    for casenum in case_arr:
               
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
                fig.savefig(f'{figout_dir}/rdf_{sysname}_{rdfkey}.png',dpi = fig.dpi)
                fig.savefig(f'{figout_dir}/rdf_{sysname}_{rdfkey}.eps',format = 'eps')

        if(neigh_flag):
            
            for nid, (neighkey,colval) in enumerate(neigh_keys.items()):

                # Plot data
                fig, ax = plt.subplots()
                paux.set_axes(ax,plt,r'$n_{neigh}$',r'$f$($n_{neigh}$)')
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
                fig.savefig(f'{figout_dir}/neigh_{sysname}_{neighkey}.png',dpi = fig.dpi)
                fig.savefig(f'{figout_dir}/neigh_{sysname}_{neighkey}.eps',format = 'eps')

        if(clust_flag):
            
            for nid, (clustkey,colval) in enumerate(clust_keys.items()):
                # Plot data
                fig, ax = plt.subplots()
                paux.set_axes(ax,plt,r'$C_{size}$',r'$f$($C_{size}$)')
                bar_width = 0.2
                maxclust = 1600 
                
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
                    
                
                ax.set_xlim([1, 10])
                #ax.set_ylim([0, 0.02])
                xtick_centers = [r + bar_width / 2 for r in range(1,10,1)] #maxclust+1,100)]
                plt.xticks(xtick_centers, np.arange(1,10,1))#np.arange(1300,maxclust+1,100))
                ax.legend()
                fig.savefig(f'{figout_dir}/clust_{sysname}_{clustkey}.png',dpi = fig.dpi)
                fig.savefig(f'{figout_dir}/clust_{sysname}_{clustkey}.eps',format = 'eps')

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

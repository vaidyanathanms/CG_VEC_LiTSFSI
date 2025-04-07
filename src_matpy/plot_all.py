# Plot all data
# Version: Mar-04-2025
# Import modules
import os
import math
import numpy as np
import analyze_rdfs as ardf
#from compute_props import *
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
syspref  = '_qanion_'
respref  = 'allresults_'
q_anion  = '0.2' # should be the exact string value
sys_arr  = ['S1','S2','S3']
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
restime_pref = ['autocorrcion_config','autocorrpol_config_']
diff_pref    = ['iondiff_config_', 'countiondiff_config_']
clust_pref   = 'clust_config_'
bfrdf_pref   = 'freeboundrdf_config_'

#---------Flags---------------------------------------------------
rdf_flag     = 0; rdf_keys    = {'LiAn':2,'LiPoly':3}
rg_flag      = 0
neigh_flag   = 1; neigh_keys  = {'Li-An':3,'An-Li':5}; maxneigh = 6
restime_flag = 0
diff_flag    = 0
clust_flag   = 0
bfrdf_flag   = 0

#---------Calculations----------------------------------------------
nanions  = [nchains*int(molwt*float(x)) for x in fran_arr]
npols    = [nchains*molwt - nanions[i]  for i in range(len(nanions))]
totbeads = [int(2.0*(npols[i] + nanions[i])) for i in range(len(npols))]
npol_ans = [int(0.5*float(x)) for x in npols]
frho_an  = [rho_tot*(nanions[i]/totbeads[i]) for i in range(len(totbeads))]
frho_pol = [rho_tot*(npol_ans[i]/totbeads[i]) for i in range(len(totbeads))]  

#-----Sanity checks------------------------------------------------
if len(fran_arr) != len(deltat):
    raise RuntimeError('Mismatch in lengths between delta_t and fran_arr')

#-------Main analysis----------------------------------------------
for sysid,sysname in enumerate(sys_arr):
    
    sys_dir = all_dir + '/' + sysname + syspref + q_anion
    if not os.path.isdir(sys_dir):
        raise RuntimeError("FATAL ERROR: " + sys_dir + " not found")

    for casenum in case_arr:
               
        if(rdf_flag):
                
            for cid, (rdfkey,colval) in enumerate(rdf_keys.items()):

                # Plot data
                fig, ax = plt.subplots()
                paux.set_axes(ax,plt,r'$r$ ($\sigma$)',r'g(r), n(r)')

                for fr_id, frac_an in enumerate(fran_arr):

                    if not 'Poly' in rdfkey:
                        rhof = frho_an[fr_id]
                    else:
                        rhof = frho_pol[fr_id]*(1-f_unpoly[sysid])

                    result_dir = sys_dir + '/' + respref + sysname + '_' + frac_an + \
                        '_' + str(casenum)

                    if not os.path.isdir(result_dir):
                        raise RuntimeError("FATAL ERROR: " + result_dir + " not found")

                    print ("RDF analysis for ", result_dir)

                    rarr,gofr,nofr = ardf.ana_rdf(result_dir,rdf_pref,colval-1)
                    plt.plot(rarr, gofr, color = clr_arr[fr_id], linestyle = '--', \
                             label = '$f_{an}$ = ' + str(frac_an))
                    plt.plot(rarr, rhof*nofr, color = clr_arr[fr_id], linestyle='-')

                ax.set_xlim([0, 3])
                ax.legend()
                fig.savefig(f'{figout_dir}/rdf_{sysname}_{rdfkey}.png',dpi = fig.dpi)
                fig.savefig(f'{figout_dir}/rdf_{sysname}_{rdfkey}.eps',format = 'eps')

        if(neigh_flag):
            
            for nid, (neighkey,colval) in enumerate(neigh_keys.items()):

                # Plot data
                fig, ax = plt.subplots()
                paux.set_axes(ax,plt,r'$n_{neigh}$',r'$p$($n_{neigh}$)')
                width = 0.2
                
                all_neigh_data = np.zeros((maxneigh,len(fran_arr)+1))
                for fr_id, frac_an in enumerate(fran_arr):
                    
                    result_dir = sys_dir + '/' + respref + sysname + '_' + frac_an + \
                        '_' + str(casenum)

                    if not os.path.isdir(result_dir):
                        raise RuntimeError("FATAL ERROR: " + result_dir + " not found")
                    
                    print ("Neighbor analysis for ", result_dir)

                    xarr,yarr = paux.return_neigh_arrays(result_dir,\
                                                         neigh_pref,[0,colval-1],\
                                                         0,maxneigh)
                    if nid == 0:
                        all_neigh_data[:,0] = xarr - 1
                    all_neigh_data[:,fr_id+1] = yarr

                ax.set_xlim([-0.2, maxneigh+1])
                for i in range(len(fran_arr)):
                    ax.bar(all_neigh_data[:,0]+i*width,all_neigh_data[:,i+1]/100,\
                           width,label='$f_{an}$ = ' + fran_arr[i],\
                           color=clr_arr[i])
                ax.legend()
                plt.xticks(np.arange(0,maxneigh+1,1))
                fig.savefig(f'{figout_dir}/neigh_{sysname}_{neighkey}.png',dpi = fig.dpi)
                fig.savefig(f'{figout_dir}/neigh_{sysname}_{neighkey}.eps',format = 'eps')

        if(clust_flag):
            
            for nid, (neighkey,colval) in enumerate(neigh_keys.items()):

                # Plot data
                fig, ax = plt.subplots()
                paux.set_axes(ax,plt,r'$n_{neigh}$',r'$p$($n_{neigh}$)')
                width = 0.2
                
                all_neigh_data = np.zeros((maxneigh,len(fran_arr)+1))
                for fr_id, frac_an in enumerate(fran_arr):
                    
                    result_dir = sys_dir + '/' + respref + sysname + '_' + frac_an + \
                        '_' + str(casenum)

                    if not os.path.isdir(result_dir):
                        raise RuntimeError("FATAL ERROR: " + result_dir + " not found")
                    
                    print ("Neighbor analysis for ", result_dir)

                    xarr,yarr = paux.return_neigh_arrays(result_dir,\
                                                         neigh_pref,[0,colval-1],\
                                                         0,maxneigh)
                    if nid == 0:
                        all_neigh_data[:,0] = xarr - 1
                    all_neigh_data[:,fr_id+1] = yarr

                ax.set_xlim([-0.2, maxneigh+1])
                for i in range(len(fran_arr)):
                    ax.bar(all_neigh_data[:,0]+i*width,all_neigh_data[:,i+1]/100,\
                           width,label='$f_{an}$ = ' + fran_arr[i],\
                           color=clr_arr[i])
                ax.legend()
                plt.xticks(np.arange(0,maxneigh+1,1))
                fig.savefig(f'{figout_dir}/neigh_{sysname}_{neighkey}.png',dpi = fig.dpi)
                fig.savefig(f'{figout_dir}/neigh_{sysname}_{neighkey}.eps',format = 'eps')

            
        


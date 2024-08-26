import numpy
import os
import shutil
import subprocess
import sys
import glob
from subprocess import call


#---------mypython functions------------------------------

from my_python_functions import cpy_main_files
from my_python_functions import compile_anafiles
from my_python_functions import find_datafyle
from my_python_functions import find_latest_trajfyle
from my_python_functions import edit_generate_anainp_files
from my_python_functions import run_analysis

#---------input flags------------------------------------------
analysis_type = 2 #1-polydisp_pe,2-newparams,3-mwchange,#4-oldsize

#---------input details----------------------------------------
free_chains  = [150]
free_avg_mw  = 30
graft_chains = 32
graft_avg_mw = 35 
tail_mons    = 5
nsalt        = 510
f_charge     = 0.5
archarr      = [3]
ncases_pdi   = [3]
pdi_free     = 1.0
pdi_graft    = 1.0
cutoff_dist  = 1.50 #use two decimal places

#--------file_lists--------------------------------------------

ana_files = ['pe_analyze.f90','pe_params.f90','anainp_var.txt']
job_files = ['jobana_var.sh']
traj_pref = 'config_*'

#---------directory info---------------------------------------
maindir = os.getcwd()
scratchdir = '/scratch.global/vaidya/'
jobdir = '/home/dorfmank/vsethura/allfiles/files_polydisperse/src_lmp'
   
#Main analysis
for ifree in range(len(free_chains)):
    
    print( "Free Number of Chains: ", free_chains[ifree])

    if analysis_type == 1:
        workdir1 = scratchdir + 'polydisp_pe'
    elif analysis_type == 2:
        workdir1 = scratchdir + 'newparams_polydisp_pe'
    elif analysis_type == 3:
        workdir1 = scratchdir + 'mwchange_polydisp_pe'
    elif analysis_type == 4:
        workdir1 = scratchdir + 'oldsize_polydisp_pe'
    else:
        print ("ERROR: Analysis_type unknown")
        continue
    
    if not os.path.isdir(workdir1):
        print("ERROR: ", workdir1, " not found")
        continue

    workdir2 = workdir1 + '/n_' + str(free_chains[ifree])
	
    if not os.path.isdir(workdir2):
        print("ERROR: ", workdir2, " not found")
        continue
        
    for iarch in range(len(archarr)):
             
        if archarr[iarch] == 1:
            print( "Archval: Block_Block")
            dirstr = 'bl_bl'
            fylstr = 'block_block'
        elif archarr[iarch] == 2:
            print( "Archval: Block_Alter")
            dirstr = 'bl_al'
            fylstr = 'block_alter'
        elif archarr[iarch] == 3:
            print( "Archval: Alter_Block")
            dirstr = 'al_bl'
            fylstr = 'alter_block'
        elif archarr[iarch] == 4:
            print( "Archval: Alter_Alter")
            dirstr = 'al_al'
            fylstr = 'alter_alter'
        else:
            print( "Unknown Architecture")
            
        workdir_arch = workdir2 + '/' + dirstr

        if not os.path.isdir(workdir_arch):
            print("ERROR: ", workdir_arch, " not found")
            continue

        workdir_freepdi = workdir_arch + '/pdi_free_' + str(pdi_free)
        if not os.path.isdir(workdir_freepdi):
            print("ERROR: ", workdir_freepdi, " not found")
            continue

        workdir_graftpdi = workdir_freepdi + '/pdi_graft_' + str(pdi_graft)
        if not os.path.isdir(workdir_graftpdi):
            print("ERROR: ", workdir_graftpdi, " not found")
            continue

        for casenum in range(len(ncases_pdi)):

            workdir_subpdi = workdir_graftpdi + '/Case_' + str(ncases_pdi[casenum])
            if not os.path.isdir(workdir_subpdi):
                print("ERROR: ", workdir_subpdi, " not found")
                continue

            os.chdir(workdir_subpdi)
            destdir = os.getcwd()

            print( "Starting analysis for case ", ncases_pdi[casenum], "in ",\
                       free_chains[ifree],dirstr)

            #---Copying files------
            print( "Current Dir ", destdir)
            print( "Copying Files")
                
            for fyllist in range(len(ana_files)):
                cpy_main_files(maindir,destdir,ana_files[fyllist])

            for fyllist in range(len(job_files)):
                cpy_main_files(jobdir,destdir,job_files[fyllist])

            dataname =  find_datafyle(free_chains[ifree],fylstr,\
                                      ncases_pdi[casenum],destdir)
            if dataname == 'ERROR':
                print("ERROR: No restart files found")
                continue

            
            #latest_traj = find_latest_trajfyle(traj_pref,destdir)
            #if latest_traj == 'ERROR':
            #    print("ERROR: No trajectory files found")
            #    continue
            
            compile_anafiles()            
            ntotch = free_chains[ifree] + graft_chains
            traj_arr = glob.glob(traj_pref)
            if traj_arr == []:
                print("ERROR: No trajectory files found")
                continue
            
            for fyllist in range(len(traj_arr)):
                print("Analyzing ", traj_arr[fyllist])
                edit_generate_anainp_files(dataname,traj_arr[fyllist],ntotch,\
                                           free_chains[ifree],graft_chains,\
                                           cutoff_dist,fyllist+1)
                outana = 'jobana_' + str(fyllist+1) + '.sh'
                run_analysis(free_chains[ifree],pdi_free,ncases_pdi[casenum],\
                             fylstr,'jobana_var.sh',outana,fyllist+1,destdir)
                



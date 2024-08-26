# To recreate datafiles for computing PDI and initial MW distributions


import numpy
import os
import shutil
import subprocess
import sys
import glob

#---------------Import Functions-----------------------------------

from subprocess import call
from my_python_functions import my_cpy_generic

#---------input flags------------------------------------------
analysis_type  = 2 #1-polydisp_pe,2-newparams,3-mwchange,#4-oldsize

#---------input details----------------------------------------
free_chains  = [150] #[16,32,64,128,150]
free_avg_mw  = 30
graft_chains = 32
graft_avg_mw = 35 
tail_mons    = 5
nsalt        = 510
f_charge     = 0.5
archarr      = [2,3]
ncases_pdi   = [1,2,3,4]
pdi_free     = 1.0
pdi_graft    = 1.0
cutoff_dist  = 1.50 #use two decimal places

#--------file_lists--------------------------------------------
#Give prefix for files to be copied followed by *
restfyle_name = 'archival*'

#---------directory info---------------------------------------
maindir = os.getcwd()
scratchdir = '/scratch.global/vaidya/'
lmp_dir = '/home/dorfmank/vsethura/mylammps/src'
lmp_fyle = 'lmp_mesabi'

if not os.path.isdir(scratchdir):
    print(scratchdir, "does not exist")
    sys.exit()

#--------------Main Analysis------------------------------------

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

        #------------------Make global analysis files output directory-----
        
        data_superdir = workdir1 + '/data_all_dir'
        if not os.path.isdir(data_superdir):
            os.mkdir(data_superdir)
        
        data_dirname = 'n_' + str(free_chains[ifree])
        data_main_dir = data_superdir + '/' + data_dirname
        if not os.path.isdir(data_main_dir):
            os.mkdir(data_main_dir)

        data_arch_dir = data_main_dir + '/' + dirstr
        if not os.path.isdir(data_arch_dir):
            print("making", data_arch_dir)
            os.mkdir(data_arch_dir)

        data_pdi_dir = data_arch_dir + '/pdifree' + \
                       str(pdi_free) + '_pdigraft_' + str(pdi_graft)
        if not os.path.isdir(data_pdi_dir):
            print("making", data_pdi_dir)
            os.mkdir(data_pdi_dir)

        #---------------------------------------------------

        for casenum in range(len(ncases_pdi)):

            workdir_subpdi = workdir_graftpdi + '/Case_' + str(ncases_pdi[casenum])
            if not os.path.isdir(workdir_subpdi):
                print("ERROR: ", workdir_subpdi, " not found")
                continue

            os.chdir(workdir_subpdi)
            destdir = os.getcwd()

            print( "Starting analysis for case ", ncases_pdi[casenum], "in ",\
                       free_chains[ifree],dirstr)


            #----Make case directory in pdi_directory----------

            data_dir = data_pdi_dir + '/' + \
                       'Case_' + str(ncases_pdi[casenum])
            if not os.path.isdir(data_dir):
                os.mkdir(data_dir)

            #------All copying/manipulations--------------------------

            print("Copying output files from", workdir_subpdi)

            #search for restart files in the directory
            os.chdir(workdir_subpdi)
            destdir = os.getcwd()
            list_of_files = glob.glob(restfyle_name)
            
            if list_of_files == []: # if list is empty here too
                print("Did not find files of type: ", restfyle_name)
                continue

            # All restart files should have same atom type details
            fylename = max(list_of_files, key=os.path.getctime)
            print( "Copying file: ", fylename)

            dataname = 'PEinitdata.txt'
            if not os.path.exists(lmp_fyle):
                my_cpy_generic(lmp_dir,destdir,lmp_fyle,lmp_fyle)

            subprocess.call(["mpirun","-np","48","./lmp_mesabi","-r"\
                             ,fylename,dataname])
            my_cpy_generic(destdir,data_dir,dataname,dataname)
            

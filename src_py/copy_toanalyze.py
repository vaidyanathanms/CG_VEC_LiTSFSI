# Version: May-22-2020
# To copy the results and restart files for further analysis

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
rest_job_flag  = 0 #copy restart/job files
fyl_flag       = 1 #copy output files

#---------input details----------------------------------------
free_chains  = [0]
free_avg_mw  = 30
graft_chains = 32
graft_avg_mw = 35 
tail_mons    = 5
nsalt        = 510
f_charge     = 0.5
archarr      = [4]
ncases_pdi   = [1]
pdi_free     = 1.5
pdi_graft    = 1.0
cutoff_dist  = 1.00 #use two decimal places


#--------file_lists--------------------------------------------
#Give prefix for files to be copied followed by *
#fyl_list = ['rgall*','rgavg*']
fyl_list     = ['adsfrac*','tether*','chainadsval*','log*',\
                'dens*','chdens*','PErdf*','chgrpdens*','grpdens*'\
                ,'polydens*','anainp*','rgall*','rgavg*']
restart_list = ['restart*','archival*','job*']

#---------directory info---------------------------------------
maindir = os.getcwd()
scratchdir = '/scratch.global/vaidya/'
jobdir = '/home/dorfmank/vsethura/allfiles/files_polydisperse/src_lmp'

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
        if fyl_flag == 1:

            out_dir = 'outresults_dir_n_' + str(free_chains[ifree])         
            anafyl_main_dir = workdir2 + '/' + out_dir
            if not os.path.isdir(anafyl_main_dir):
                print("making", anafyl_main_dir)
                os.mkdir(anafyl_main_dir)

            anafyl_arch_dir = anafyl_main_dir + '/' + dirstr
            if not os.path.isdir(anafyl_arch_dir):
                print("making", anafyl_arch_dir)
                os.mkdir(anafyl_arch_dir)

            anafyl_pdi_dir = anafyl_arch_dir + '/pdifree_' + \
                             str(pdi_free) + '_pdigraft_' + str(pdi_graft)

            if not os.path.isdir(anafyl_pdi_dir):
                print("making", anafyl_pdi_dir)
                os.mkdir(anafyl_pdi_dir)

        if rest_job_flag == 1:

            restart_dirname = 'restart_all_dir_n_' + str(free_chains[ifree])
            restart_main_dir = workdir2 + '/' + restart_dirname
            if not os.path.isdir(restart_main_dir):
                os.mkdir(restart_main_dir)

            restart_arch_dir = restart_main_dir + '/' + dirstr
            if not os.path.isdir(restart_arch_dir):
                print("making", restart_arch_dir)
                os.mkdir(restart_arch_dir)

            restart_pdi_dir = restart_arch_dir + '/pdifree' + \
                              str(pdi_free) + '_pdigraft_' + str(pdi_graft)
            if not os.path.isdir(restart_pdi_dir):
                print("making", restart_pdi_dir)
                os.mkdir(restart_pdi_dir)

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
            if fyl_flag == 1:
                anafyl_dir = anafyl_pdi_dir + '/' + \
                             'Case_' + str(ncases_pdi[casenum])
                if not os.path.isdir(anafyl_dir):
                    os.mkdir(anafyl_dir)
            if rest_job_flag == 1:
                restart_dir = restart_pdi_dir + '/' + \
                              'Case_' + str(ncases_pdi[casenum])
                if not os.path.isdir(restart_dir):
                    os.mkdir(restart_dir)

            #------------------------------------------------------------


            #------All copying/manipulations--------------------------
            result_dir = 'results_'+ str(free_chains[ifree]) + '_' \
                         + str(pdi_free) + '_' + str(ncases_pdi[casenum])
            workdir_results = workdir_subpdi + '/' + result_dir
            if not os.path.isdir(workdir_results):
                print(workdir_results, " does not exist")
                continue

            print("Copying output files from", workdir_results)

            if fyl_flag == 1:

                for fylcnt in range(len(fyl_list)):

                    #search in results_directory
                    os.chdir(workdir_results)
                    destdir = os.getcwd()
                    list_of_files = glob.glob(fyl_list[fylcnt])

                    if list_of_files == []: #if list is empty in results dir

                        #search in previous directory
                        os.chdir(workdir_subpdi)
                        destdir = os.getcwd()
                        list_of_files = glob.glob(fyl_list[fylcnt])

                        if list_of_files == []: #if list is empty here too
                            print("Did not find files of type", \
                                  fyl_list[fylcnt])
                            continue

                    print( "Copying files of type: ", fyl_list[fylcnt])

                    for filenum in range(len(list_of_files)):

                        fylname = list_of_files[filenum]
                        anafylname = fylname
                        my_cpy_generic(destdir,anafyl_dir,fylname,anafylname)

            if rest_job_flag == 1:

                for fylcnt in range(len(restart_list)):

                    #search only in main directory
                    os.chdir(workdir_subpdi)
                    destdir = os.getcwd()
                    list_of_files = glob.glob(restart_list[fylcnt])

                    if list_of_files == []: #
                        print("Did not find files of type", \
                              restart_list[fylcnt])
                        continue

                    print( "Copying files of type: ", restart_list[fylcnt])
                    for filenum in range(len(list_of_files)):

                        fylname = list_of_files[filenum]
                        anafylname = fylname
                        my_cpy_generic(destdir,restart_dir,fylname,anafylname)

                

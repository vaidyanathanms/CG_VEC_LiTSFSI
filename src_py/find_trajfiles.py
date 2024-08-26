# Version: July-07-2020
# To find all the trajectory files within directories

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
free_chains  = [64]
free_avg_mw  = 30
graft_chains = 32
graft_avg_mw = 35 
tail_mons    = 5
nsalt        = 510
f_charge     = 0.5
archarr      = [1,4]
ncases_pdi   = [1,2,3,4]
pdi_free     = 1.5
pdi_graft    = 1.0

#--------file_lists--------------------------------------------
#Give prefix for files to be copied followed by *
fyl_list     = ['config*']

#---------directory info---------------------------------------
maindir = os.getcwd()
scratchdir = '/scratch.global/vaidya/'
jobdir = '/home/dorfmank/vsethura/allfiles/files_polydisperse/src_lmp'

if not os.path.isdir(scratchdir):
    print(scratchdir, "does not exist")
    sys.exit()

#--------create output file------------------------------------
trajfyl = maindir + '/all_trajfiles_pdi_' + str(pdi_free) + '.txt'
fyle_id = open(trajfyl,'w')
fyle_id.write('%s\t %s\t %s\t %s\t %s\n' %('nfree','arch',\
                                           'pdi_free','case_num',\
                                           'filename'))

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



        #------------------Make global analysis file------------------------

        for casenum in range(len(ncases_pdi)):

            print("Finding trajecory files in ", ncases_pdi[casenum])
            workdir_subpdi = workdir_graftpdi + '/Case_' + str(ncases_pdi[casenum])
            if not os.path.isdir(workdir_subpdi):
                print("ERROR: ", workdir_subpdi, " not found")
                continue

            os.chdir(workdir_subpdi)
            destdir = os.getcwd()

            for fylcnt in range(len(fyl_list)):

                #search in previous directory
                os.chdir(workdir_subpdi)
                destdir = os.getcwd()
                list_of_files = glob.glob(fyl_list[fylcnt])
                
                if list_of_files == []: #if list is empty here too
                    print("Did not find files of type", \
                          fyl_list[fylcnt], "in casenum", \
                          ncases_pdi[casenum])
                    continue

                for filenum in range(len(list_of_files)):

                    traj_fylename = list_of_files[filenum]
                    fyle_id.write('%d\t%s\t%g\t%d\t%s\n' %(free_chains[ifree],\
                                  dirstr,pdi_free,ncases_pdi[casenum],\
                                  traj_fylename))
                    


                
fyle_id.close()

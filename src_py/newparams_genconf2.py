# Edited to avoid restart flags
# New Version: April-14-2020
# The path to the source file for LAMMPS (search keyword "srcfyl") needs to be updated
import numpy
import os
import shutil
import subprocess
import sys
import glob
import re
from subprocess import call


#---------mypython functions------------------------------

from my_python_functions import cpy_main_files
from my_python_functions import init_pdi_write
from my_python_functions import compile_and_run_pdi
from my_python_functions import check_pdi_files
from my_python_functions import create_paramfyl_for_datafyl
from my_python_functions import compile_and_run_inpgenfyles
from my_python_functions import edit_generate_input_lmp_files
from my_python_functions import run_lammps
from my_python_functions import clean_backup_initfiles
from my_python_functions import find_recent_file

#---------input flags------------------------------------------
#0-initial run  1- production
num_hrs   = 48 # Total number of hours for run
num_nodes = 2  # Number of nodes
num_procs = 24 # Number of procs per node
maxtime   = 40000000 # Maximum timesteps
equiltime = 15000000 # Equilibrium timesteps

#---------input details----------------------------------------
free_chains  = [0]
free_avg_mw  = 30
graft_chains = 32
graft_avg_mw = 35 
tail_mons    = 5
nsalt        = 510
f_charge     = 0.5
archarr      = [4]
ncases_pdi   = [1,2]
pdi_free     = 1.5
pdi_graft    = 1.0

box_data     = [35.0, 35.0, 120.0] #x-box, y-box, z-box

#--------file_lists--------------------------------------------

f90_files = ['ran_numbers.f90','lammps_inp.f90','lmp_params.f90'\
             ,'SZDist2.f90','init_pdi.txt','polyinp_fordiffbox_var.dat']
lmp_files = ['in.longrun','in.init_var','in.run1','jobmain_var.sh']
lmp_long  = ['in.longrun','jobmain_long_var.sh']
lmp_equil = ['in.run2','in.longrun','jobmain_equil_var.sh']

pdi_files = ['FreeChains.dat','GraftChains.dat'] #Order: freepdi,graftpdi
par_files = ['polyinp_fordiffbox_var.dat','polyinp.txt'] #Order: varfile,oufyle
archive_file_style = 'archival*'
#---------directory info---------------------------------------
maindir = os.getcwd()
scratchdir = '/scratch.global/vaidya/'
lmpdir = '/home/dorfmank/vsethura/allfiles/files_polydisperse/src_lmp'

   
#Main analysis
for ifree in range(len(free_chains)):
    
    print( "Free Number of Chains: ", free_chains[ifree])
    workdir1 = scratchdir + 'newparams_polydisp_pe'
    
    if not os.path.isdir(workdir1):
        os.mkdir(workdir1)

    workdir2 = workdir1 + '/n_' + str(free_chains[ifree])
	
    if not os.path.isdir(workdir2):
        os.mkdir(workdir2)
        
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
            os.mkdir(workdir_arch)

        workdir_freepdi = workdir_arch + '/pdi_free_' + str(pdi_free)
        if not os.path.isdir(workdir_freepdi):
            os.mkdir(workdir_freepdi)

        workdir_graftpdi = workdir_freepdi + '/pdi_graft_' + str(pdi_graft)
        if not os.path.isdir(workdir_graftpdi):
            os.mkdir(workdir_graftpdi)

        for caselen in range(len(ncases_pdi)):

            archive_flag = 1 # archived files are present

            workdir_subpdi = workdir_graftpdi + '/Case_' + str(ncases_pdi[caselen])

            if not os.path.isdir(workdir_subpdi):
                archive_flag = 0 #if dir is not there, file is not there
                print("Creating new simulation")
                os.mkdir(workdir_subpdi)


            if archive_flag == 1:

                os.chdir(workdir_subpdi)
                print("Current path: ", workdir_subpdi)
                list_of_files = glob.glob(archive_file_style)
                if list_of_files == []:
                    print("No archive files are found: restarting simulations")
                    archive_flag = 0 #no files are found. reset flag
                else:
                    print("Archive files found: RESTARTING from archives")


            if archive_flag == 0:

                os.chdir(workdir_subpdi)
                destdir = os.getcwd()

                print( "Starting case", ncases_pdi[caselen], "for ",\
                       free_chains[ifree],dirstr)

                #---Copying files------
                print( "Current Dir ", destdir)
                print( "Copying Files")
                
                for fyllist in range(len(f90_files)):
                    cpy_main_files(maindir,destdir,f90_files[fyllist])

                for fyllist in range(len(lmp_files)):
                    cpy_main_files(lmpdir,destdir,lmp_files[fyllist])

                srcfyl = '/home/dorfmank/vsethura/mylammps/src/lmp_mesabi'
                desfyl = destdir + '/lmp_mesabi'
                shutil.copy2(srcfyl, desfyl)

                #----Run PDI code------
                print("Running PDI calculation file")
                tot_chains = free_chains[ifree] + graft_chains
                init_pdi_write(pdi_free,free_avg_mw,free_chains[ifree],\
                               pdi_graft,graft_avg_mw,graft_chains,destdir)
                compile_and_run_pdi(destdir)

                pdiflag = 1
                check_pdi_files(destdir,pdi_files,pdiflag)

                if pdiflag == -1:
                    print("Check PDI files")
                    continue

                    
                #----Generate input files-----
                print( "Copy Successful - Generating Input Files")
                
                lmp_data_fyle = create_paramfyl_for_datafyl(destdir,\
                                                            par_files,free_chains[ifree]\
                                                            ,fylstr,ncases_pdi[caselen],free_avg_mw,\
                                                            graft_chains,graft_avg_mw,nsalt,\
                                                            f_charge,tail_mons,archarr[iarch],\
                                                            pdi_files,box_data)

                compile_and_run_inpgenfyles(par_files[1],destdir)            

                #---Run LAMMPS files-------------

                edit_generate_input_lmp_files('in.init_var',lmp_data_fyle)
                run_lammps(free_chains[ifree],pdi_free,ncases_pdi[caselen],fylstr,\
                           'jobmain_var.sh','jobmain.sh',num_hrs,num_nodes,num_procs)
                
                #----Copy/Backup initial files---
                clean_backup_initfiles(f90_files,pdi_files,par_files,destdir)
                os.chdir(maindir)

            else:

                print("%s\t %g\t %d\t %s\t %d\t" 
                      %("Restarting simulation for pdi/nfree/arch/casenum",
                        pdi_free,free_chains[ifree],fylstr,ncases_pdi[caselen]))
                if not os.path.isdir(workdir_subpdi):
                    print( workdir_subpdi, "not found")
                    continue

                os.chdir(workdir_subpdi)
                destdir = os.getcwd()

                os.chdir(workdir_subpdi)
                print("Current path: ", workdir_subpdi)
                latest_restfyl = find_recent_file(destdir,archive_file_style)
                delimited_vals = re.split("\W+|_",latest_restfyl)
                timeval = delimited_vals[len(delimited_vals)-2]
                print("timeval", timeval)

                if int(timeval) > maxtime:
                    print('Max time reached: ', timeval, maxtime)
                    continue

                elif int(timeval) > equiltime:

                    print('Running production cycle..')
                    for fyllist in range(len(lmp_long)):
                        cpy_main_files(lmpdir,destdir,lmp_long[fyllist])

                    fylename = destdir + '/lmp_mesabi'
                    if not fylename: 
                        srcfyl = '/home/dorfmank/vsethura/mylammps/src/lmp_mesabi'
                        desfyl = destdir + '/lmp_mesabi'
                        shutil.copy2(srcfyl, desfyl)

                    run_lammps(free_chains[ifree],pdi_free,ncases_pdi[caselen],fylstr,\
                               'jobmain_long_var.sh','jobmain_long.sh',num_hrs,\
                               num_nodes,num_procs)


                    os.chdir(maindir)


                else:

                    print('Running long equilibrium cycle..')
                    for fyllist in range(len(lmp_equil)):
                        cpy_main_files(lmpdir,destdir,lmp_equil[fyllist])

                    fylename = destdir + '/lmp_mesabi'
                    if not fylename: 
                        srcfyl = '/home/dorfmank/vsethura/mylammps/src/lmp_mesabi'
                        desfyl = destdir + '/lmp_mesabi'
                        shutil.copy2(srcfyl, desfyl)

                    run_lammps(free_chains[ifree],pdi_free,ncases_pdi[caselen],fylstr,\
                               'jobmain_equil_var.sh','jobmain_equil.sh',num_hrs,\
                               num_nodes,num_procs)


                    os.chdir(maindir)

                    

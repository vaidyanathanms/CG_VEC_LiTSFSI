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

#---------input details----------------------------------------
analyze_only_latest = 1 # Only analyze the latest trajectory file
frac_anions  = [1/10]#,1/15,1/10]#,1/20,1/10,1/6,1/5,1/3] # fraction of anions
tot_mons     = 6000 # total number of MONOMERS in the poly CHAIN
chain_mw     = [40]#60,40]#,60,90] # of monomer range per chain
num_chains   = [int(tot_mons/x) for x in chain_mw] # of polymerized ch
unpoly_farr  = [0.6] # fraction of unpolymerized mons
nrepeats     = 1 # number of replica
nframes      = 200 # total frames to be analyzed
skipfr       = 0 # skip frames
freqfr       = 5 # freq of anaylsis


#---------job details------------------------------------------
tottime   = 3 # in hours
nnodes    = 1 # number of nodes
ncores    = 36 # number of cores
hpc_sys   = 'cades'  # Opt: kestrel, cades

#--------file_lists--------------------------------------------
ana_files = ['analyze.f90','ana_params.f90','anainp_var.txt']
job_files = ['jobana_var.sh']
traj_pref = 'config_*'
data_pref = 'VECdata_'; data_ext = '.dat'

#---------directory info---------------------------------------
maindir = os.getcwd() #src_py dir
if hpc_sys == 'kestrel':
    home_path = '/home/vaidyams'
    scr_path  = '/scratch/vaidyams'
elif hpc_sys == 'cades':
    home_path = '/home/vm5'
    scr_path  = '/lustre/or-scratch/cades-birthright/vm5'
else:
    raise RuntimeError('Unknown HPC system ' + hpc_sys)

src_f90    = home_path + '/all_codes/files_CG-SIC/src_f90' #src_f90 dir
src_lmp    = home_path + '/all_codes/files_CG-SIC/src_lmp' #src_lmp dir
scratchdir = scr_path  + '/cg_sic' #output headdir
scr_head   = 'sic_mixedvec_listsfi' # head dir scratch'

#--------lammps executable-------------------------------------
if hpc_sys == 'kestrel':
    lmpexe_dir = 'None' # Using kestrel compiled lammps
    lmp_exe    = 'lmp' # lmp executable file
    f90_comp   = 'ifx' # FORTRAN compiler
elif hpc_sys == 'cades':
    lmpexe_dir = home_path + '/mylammps/src' #lmp executable path
    lmp_exe    = 'lmp_mpi' # lmp executable file
    f90_comp   = 'ifort' # FORTRAN compiler
else:
    raise RuntimeError('Unknown HPC system ' + hpc_sys)

if not os.path.isdir(scratchdir):
    raise RuntimeError(scratchdir + " not found!")

#---------main analysis---------------------------------------
for mw_ch in range(len(chain_mw)):
    
    print( "Analyzing MW/number of Chains: ", chain_mw[mw_ch],\
           num_chains[mw_ch])
    workdir1 = scratchdir + '/' + scr_head
    
    if not os.path.isdir(workdir1):
        print("ERROR: " + workdir1 + " not found!"); continue

    workdir2 = workdir1 + '/MW_' + str(chain_mw[mw_ch])
	
    if not os.path.isdir(workdir2):
        print("ERROR: " + workdir2 + " not found!"); continue
        
    for fr_an in range(len(frac_anions)):
             
        workdir_super = workdir2 + '/fanion_' + str(round(frac_anions[fr_an],2))

        if not os.path.isdir(workdir_super):
            print("ERROR: " + workdir_super + " not found!"); continue

        if len(unpoly_farr) == 1:
            unpoly_frac = unpoly_farr[0]
        else:
            unpoly_frac = unpoly_farr[fr_an]

        for casenum in range(1,nrepeats+1):

            workdir_main = workdir_super + '/Case_' + str(casenum)
            if not os.path.isdir(workdir_main):
                print("ERROR: " + workdir_main + " not found!"); continue

            os.chdir(workdir_main)
            destdir = os.getcwd()

            print( "Analyzing case ", casenum, "for f_anion/MW_chains/nchains: ",\
                   frac_anions[fr_an],chain_mw[mw_ch],num_chains[mw_ch])

            #---Copying files------
            print( "Current Dir ", destdir)
            print( "Copying Files")
                
            for fyllist in range(len(ana_files)):
                cpy_main_files(src_f90,destdir,ana_files[fyllist])

            for fyllist in range(len(job_files)):
                cpy_main_files(src_lmp,destdir,job_files[fyllist])

            #----Generate analysis input files
            print( "Copy Successful - Generating Input Files")
            tot_chains = num_chains[mw_ch]
            
            dataname = find_datafyle(data_pref,chain_mw[mw_ch],frac_anions[fr_an],\
                                     lmpexe_dir,lmp_exe,data_ext)
            if dataname == 'ERROR':
                print("ERROR: No restart files found"); continue

            #----Retrieve trajectory files
            if analyze_only_latest:
                traj_arr = glob.glob(traj_pref,key=os.path.getctime)
            else:
                traj_arr = glob.glob(traj_pref)
            if traj_arr == []:
                print("ERROR: No trajectory files found")
                continue

            compile_anafiles()            

            #----Iterate through trajectory files and submit-----            
            for fyllist in range(len(traj_arr)):
                print("Analyzing ", traj_arr[fyllist])
                anainp = edit_generate_anainp_files(dataname,traj_arr[fyllist],\
                                                    num_chains[mw_ch],nframes,skipfr,\
                                                    freqfr,fyllist+1)
                jobana = 'jobana_' + str(fyllist+1) + '.sh'
                jobstr = 'ana_' + str(frac_anions[fr_an]) + '_' + \
                         chain_mw[mw_ch]
                run_analysis(anainp, jobstr, num_chains[mw_ch],fyllist+1,\
                             'jobana_var.sh',jobana,tottime,nnodes,ncores)                

#-------Automated code to run SIC systems VEC-LiSTFSI----
#-------Date: Aug-26-2024--------------------------------
#-------Author: Vaidyanathan Sethuraman------------------
#-------Requirements: my_python_functions.py-------------
import numpy as np
import os
import shutil
import subprocess
import sys
import glob
from subprocess import call

#---------mypython functions------------------------------

from my_python_functions import cpy_main_files
from my_python_functions import init_pdi_write
from my_python_functions import create_paramfyl_for_datafyl
from my_python_functions import compile_and_run_inpgenfyles
from my_python_functions import create_angle_topo
from my_python_functions import edit_generate_input_lmp_files
from my_python_functions import gen_pair_coeff_file
from my_python_functions import gen_bond_coeff_file
from my_python_functions import gen_angl_coeff_file
from my_python_functions import run_lammps
from my_python_functions import clean_backup_initfiles

#---------input flags------------------------------------------
#0-initial run  1- production
restart   = 0  # For restarting from given configurations
num_hrs   = 47 # Total number of hours for run
num_nodes = 6  # Number of nodes
num_cores = 36 # Number of cores per node
part_type = 'hbw' #hbw or shared (partition type)
hpc_sys   = 'kestrel'  # Opt: kestrel, cades
sys_type  = 'S1' # S1; S2; S3; S4

#---------input details - Topology--------------------------------
frac_anions  = [1/5]#,1/20]#,1/15,1/20]#, 1/10, 1/5] # fraction of anions
tot_mons     = 4000 # total number of MONOMERS in the poly CHAIN
chain_mw     = [40] # of monomer range per chain
num_chains   = [int(tot_mons/x) for x in chain_mw] # of polymerized ch
unpoly_farr  = [0.0] # fraction of unpolymerized mons
density      = 0.8 # system density
cg_per_mon   = 2 # number of blobs per polymer monomer
same_mass    = 0 # 0-> different mass, 1 -> same mass for monomers
blob_charge  = 3 # charge per blob
case_arr     = [10] # case number

if sys_type == 'S3': #1-VEC-anions blend; 0-VEC-anion copolymer
    is_anion_sep = 1 
else:
    is_anion_sep = 0
    
if sys_type == 'S1' and max(unpoly_farr) > 0:
    raise RuntimeError('Cannot have unpolymerized frac for S1')

if sys_type != 'S1' and max(unpoly_farr) == 0:
    raise RuntimeError('Need unpolymerized frac for S2-S4')

#---------input details - Pair Coeff--------------------------------
gen_pair_lst = 1 # Generate pair list
gen_bond_lst = 1 # Generate bond list
gen_angl_lst = 1 # Generate angle list

ntypes       = 4 # Number of types
lj_cutoff    = float(2**(1.0/6.0)) # LJ cut-off in sigma
coul_cutoff  = 10.0 # Coulumb cut-off in sigma
name_list    = ['VEC without C=O','C=O of VEC','STFSI'] # Name of groups
eps_list     = [1, 1, 1] # epsilon values
sig_list     = [1,1,1] # sigma values
    

if sys_type.lower() != 's1': # unpolymerized VEC
    name_list.extend(['unpoly_VEC: VEC','unpoly_VEC: C=O'])
    eps_list.extend([1,1])
    sig_list.extend([1,1])

# Cation pair coefficient values
name_list.append('Li')
eps_list.append(1)
sig_list.append(1)

#-----------input details: pair cut-off-----------------------------    
ljcut_list   = [i*lj_cutoff for i in sig_list] # lj-cut values [2.7,3.2,6,2.5]
coulcut_rat  = [1.0,1.0,1.0,1.0] #coulumbic cut off ratio
coulcut_list = [i*coul_cutoff for  i in coulcut_rat] # coul-cut values

#---------input details - Bond Coeff--------------------------------
if sys_type.lower() != 's3':

    bname_list   = ['VEC-VEC (no C=O) [1-1]',' VEC - C=O [1-2]', \
                    'VEC - STFSI [1-3]' , 'STFSI - STFSI [3-3]'] #Name of groups
    kspr_list    = [30, 50, 30, 30] # spring constants
    bcon_list    = [[1,1],[1,2],[1,3],[3,3]] # connectivity list
else:
    bname_list   = ['VEC-VEC (no C=O) [1-1]',' VEC - C=O [1-2]', \
                    'STFSI - STFSI [3-3]'] #Name of groups
    kspr_list    = [30, 50, 30] # spring constants
    bcon_list    = [[1,1],[1,2],[3,3]] # connectivity list

if sys_type.lower() != 's1': # unpolymerized VEC
    ntypes = 6
    bname_list.append('unpoly_VEC: VEC - unpoly_VEC: C=O [4-5]')
    kspr_list.append(50) 
    bcon_list.append([4,5])

#---------input details - Angle Coeff--------------------------------
if sys_type.lower() == 's3':
    atype_list   = ['*', '1', '2']
    aname_list   = ['Default angle list',' 1 -> VEC1-VEC1-VEC1 [1-1-1]', \
                    '2 -> VEC1-VEC1-VEC2'] #Name of groups
    kang_list    = [3, 3, 5] # spring constants
    thet_list    = [110,110,70] # connectivity list
else:
    atype_list   = ['*', '1', '2', '6']
    aname_list   = ['Default angle list',' 1 -> VEC1-VEC1-VEC1 [1-1-1]', \
                    '2 -> VEC1-VEC1-VEC2', ' 6 -> VEC2-VEC1-STFSI [2-1-3]'] #Name of groups
    kang_list    = [3, 3, 5, 3] # spring constants
    thet_list    = [110,110,70,110] # connectivity list
    
#--------file_lists--------------------------------------------
f90_files = ['ran_numbers.f90','lmp_params_var.f90','lammps_inp.f90'] 
lmp_files = ['in.init_var','in.nve','in.nvt','in.nvt_main','in.npt','jobmain_var.sh']
tcl_files = ['guessangle_var.tcl'] 
lmp_long  = ['in.nvt','in.npt','in.nvt_main','in.rdf','jobmain_long_var.sh']

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
    
src_f90    = home_path + '/all_codes/files_CG-SIC/src_f90' #f90 dir
src_lmp    = home_path + '/all_codes/files_CG-SIC/src_lmp' #lmp dir
src_tcl    = home_path + '/all_codes/files_CG-SIC/src_tcl' #tcl dir
scratchdir = scr_path  + '/cg_sic' #output headdir
if same_mass:
    scr_head   = sys_type + '_sic_samemass_' + str(blob_charge) # head dir scratch'
else:
    scr_head   = sys_type + '_sic_mixedvec_listfsi_' + str(blob_charge) # head dir scratch'

#--------lammps executable-------------------------------------
if hpc_sys == 'kestrel':
    lmpexe_dir = 'None' # Using kestrel compiled lammps
    lmp_exe    = 'lmp' # lmp executable file
    f90_comp   = 'ifx' # FORTRAN compiler
elif hpc_sys == 'cades':
    lmpexe_dir = home_path + '/mylammps/src' #lmp executable path
    lmp_exe    = 'lmp_mpi' # lmp executable file
    f90_comp   = 'ifort'
else:
    raise RuntimeError('Unknown HPC system ' + hpc_sys)

num_cores = int(2*num_cores) if part_type == 'hbw' else num_cores

#--------Create scratchdir------------------------------------
if not os.path.isdir(scratchdir):
    os.mkdir(scratchdir)

#---------main analysis---------------------------------------
for mw_ch in range(len(chain_mw)):
    
    print( "MW/number of Chains: ", chain_mw[mw_ch],num_chains[mw_ch])
    workdir1 = scratchdir + '/' + scr_head
    
    if not os.path.isdir(workdir1):
        os.mkdir(workdir1)

    workdir2 = workdir1 + '/MW_' + str(chain_mw[mw_ch])
    
    if not os.path.isdir(workdir2):
        os.mkdir(workdir2)
        
    for fr_an in range(len(frac_anions)):
        
        workdir_super = workdir2 + '/fanion_' + str(round(frac_anions[fr_an],2))

        if not os.path.isdir(workdir_super):
            os.mkdir(workdir_super)

        if len(unpoly_farr) == 1:
            unpoly_frac = unpoly_farr[0]
        else:
            unpoly_frac = unpoly_farr[fr_an]

        for casenum in case_arr:

            workdir_main = workdir_super + '/Case_' + str(casenum)
            if not os.path.isdir(workdir_main):
                os.mkdir(workdir_main)

            if restart == 0:

                os.chdir(workdir_main)
                destdir = os.getcwd()

                print( "Starting case", casenum, "for f_anion/MW_chains: ",\
                       frac_anions[fr_an],chain_mw[mw_ch])

                #---Copying files------
                print( "Current Dir ", destdir)
                print( "Copying Files")
                
                for fyllist in range(len(f90_files)):
                    cpy_main_files(src_f90,destdir,f90_files[fyllist])

                for fyllist in range(len(lmp_files)):
                    cpy_main_files(src_lmp,destdir,lmp_files[fyllist])

                for fyllist in range(len(tcl_files)):
                    cpy_main_files(src_tcl,destdir,tcl_files[fyllist])

                if lmpexe_dir.lower() != 'None'.lower():
                    srcfyl = lmpexe_dir + '/' + lmp_exe
                    desfyl = destdir + '/' + lmp_exe
                    shutil.copy2(srcfyl, desfyl)

                #----Generate input files-----
                print( "Copy Successful - Generating Input Files")
                tot_chains = num_chains[mw_ch]
                lmp_par,lmp_data_fyle = create_paramfyl_for_datafyl(destdir,'lmp_params_var.f90',num_chains[mw_ch],\
                                                                    chain_mw[mw_ch],casenum,\
                                                                    round(frac_anions[fr_an],2),density,\
                                                                    cg_per_mon,blob_charge,unpoly_frac,\
                                                                    is_anion_sep,same_mass)

                compile_and_run_inpgenfyles(lmp_par,destdir,f90_comp=f90_comp)            


                #----Generate angle details-----
                print("Generating angle topology")
                ang_lmp_data_fyle = "topo_" + lmp_data_fyle
                create_angle_topo(destdir,tcl_files[0],lmp_data_fyle,ang_lmp_data_fyle)

                #----Generate pair/bond coefficient list----
                if gen_pair_lst:
                    gen_pair_coeff_file(destdir,ntypes,name_list,eps_list,sig_list,\
                                        ljcut_list,coulcut_list)

                if gen_bond_lst:
                    gen_bond_coeff_file(destdir,bname_list,kspr_list,sig_list,bcon_list)

                if gen_angl_lst:
                    gen_angl_coeff_file(destdir,atype_list,aname_list,kang_list,thet_list)
                    
                #---Run LAMMPS files-------------
                edit_generate_input_lmp_files('in.init_var',ang_lmp_data_fyle)
                run_lammps(sys_type,chain_mw[mw_ch],round(frac_anions[fr_an],2),casenum,\
                           'jobmain_var.sh','jobmain.sh',num_hrs,num_nodes,num_cores,part_type)
                
                #----Copy/Backup initial files---
                clean_backup_initfiles(f90_files,tcl_files,lmp_data_fyle,ang_lmp_data_fyle,destdir)
                os.chdir(maindir)


            else:

                print( "Starting case", casenum, "for f_anion/MW_chains: ",\
                       frac_anions[fr_an],chain_mw[mw_ch])

                if not os.path.isdir(workdir_main):
                    print( workdir_main, "not found")
                    continue

                os.chdir(workdir_main)
                destdir = os.getcwd()

                print('Current working directory: ' + destdir)

                #----Check for previous restart files------------
                archfiles = destdir + '/archival*'
                list_of_files = glob.glob(archfiles)

                if not list_of_files:
                    print("No archival files found in ", destdir)
                    continue

                #----Generate pair coefficient list--------------
                if gen_pair_lst:
                    gen_pair_coeff_file(destdir,ntypes,name_list,eps_list,sig_list,\
                                        ljcut_list,coulcut_list)
                if gen_bond_lst:
                    gen_bond_coeff_file(destdir,bname_list,kspr_list,sig_list,bcon_list)

                if gen_angl_lst:
                    gen_angl_coeff_file(destdir,atype_list,aname_list,kang_list,thet_list)

                #----Copy other required files-------------------
                for fyllist in range(len(lmp_long)):
                    cpy_main_files(src_lmp,destdir,lmp_long[fyllist])

                fylename = destdir + '/' + lmp_exe
                if not fylename:
                    if not lmpexe_dir.lower() != 'None'.lower():
                        srcfyl = lmpexe_dir + '/' + lmp_exe
                        desfyl = destdir + '/' + lmp_exe
                        shutil.copy2(srcfyl, desfyl)


                run_lammps(sys_type,chain_mw[mw_ch],round(frac_anions[fr_an],2),casenum,\
                           'jobmain_long_var.sh','jobmain_long.sh',num_hrs,num_nodes,num_cores)

                os.chdir(maindir)
	        

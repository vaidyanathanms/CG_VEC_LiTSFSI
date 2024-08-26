# Generic Function Definitions
# Version_1: V_May_31_2024

import numpy
import os
import shutil
import subprocess
import sys
import glob
import re

def cpy_main_files(dum_maindir,dum_destdir,fylname):

    srcfyl = dum_maindir + '/' + fylname

    if not os.path.exists(srcfyl):
        print('ERROR', srcfyl, 'not found')
        return

    desfyl = dum_destdir + '/' + fylname
    shutil.copy2(srcfyl, desfyl)

def init_pdi_write(freepdi, freemw,freechains,\
                   graftpdi,graftmw,graftch,destdir):

    pdi_fyl = destdir + '/init_pdi.txt'
    finit   = open(pdi_fyl,'w')
    finit.write('%s\t %s\n' %('free_data', '#pdi mw nchains'))
    finit.write('%g\t %g\t %g\n' %(freepdi,freemw,freechains))
    finit.write('%s\t %s\n' %('graft_data', '#pdi mw nchains'))
    finit.write('%g\t %g\t %g\n' %(graftpdi,graftmw,graftch))
    finit.close()


def compile_and_run_pdi(destdir):
    
    os.chdir(destdir)
    if not os.path.exists('init_pdi.txt'):
        print('init_pdi.txt not found')
        return

    if not os.path.exists('SZDist2.f90'):
        print('SZDist2.f90 not found')
        return
        
    # Generate PDI data

    print( "Running FORTRAN script for generating PDI files")
    subprocess.call(["ifort","-r8","-check","-traceback",
                     "SZDist2.f90","-o","pdiinp.o"])

    subprocess.call(["./pdiinp.o"])


def check_pdi_files(destdir,pdi_files,flagcheck):

    os.chdir(destdir)
    for fyl in pdi_files:
        if not os.path.exists(fyl):
            raise RuntimeError(fyl, 'not found')            


def create_paramfyl_for_datafyl(destdir,inpfyle,nchains,mw_chain\
                                ,casenum,fr_an,dens=0.85,cg_per_mon=2,\
                                blob_charge = 0.25):


    lmpinp = "lmpinp_" + str(mw_chain)+"_" + str("{:.2f}".format(fr_an)) + "_" +str(casenum)+".f90"
    logout = "log_" + str(mw_chain)+"_" + str("{:.2f}".format(fr_an)) + "_" +str(casenum)+".dat"
    datafyle = "VECdata_" + str(mw_chain)+"_" + str("{:.2f}".format(fr_an)) + ".dat"

    fr  = open(inpfyle,'r')
    fw  = open(destdir + '/' + lmpinp,'w')

    fid = fr.read().replace("py_nchains",str(nchains)).\
          replace("py_mwchain",str(mw_chain)).\
          replace("py_casenum",str(casenum)).\
          replace("py_fracanions",str(fr_an)).\
          replace("py_chargblob",str(blob_charge)).\
          replace("py_cgpermon",str(cg_per_mon)).\
          replace("py_density",str(dens))

    fw.write(fid)
    fw.close()
    fr.close()

    return lmpinp,datafyle

def compile_and_run_inpgenfyles(lmpinp,destdir):

    os.chdir(destdir)
    if not os.path.exists('ran_numbers.f90'):
        raise RuntimeError('ran_numbers.f90 not found in ' + destdir)
    if not os.path.exists(lmpinp):
        raise RuntimeError(lmpinp + ' not found in ' + destdir)
    if not os.path.exists('lammps_inp.f90'):
        raise RuntimeError('lammps_inp.f90 not found in ' + destdir)

    subprocess.call(['ifort','-r8','-qopenmp','-check','-traceback',\
                     'ran_numbers.f90',lmpinp,'lammps_inp.f90',\
                     '-o','inpgen.o'])
    subprocess.call('./inpgen.o')


def create_angle_topo(destdir,tcl_infyle,inpfyle,outfyle):
    if not os.path.exists(tcl_infyle):
        raise RuntimeError('ERROR: ' + tcl_infyle + ' not found')

    fr  = open(tcl_infyle,'r')
    fw  = open('gen_angle.tcl','w')
    fid = fr.read().replace("py_inpname",inpfyle).\
          replace("py_outname",outfyle)
    fw.write(fid)
    fw.close()
    fr.close()

    subprocess.call(['vmd','-dispdev','text','-e','gen_angle.tcl'])


def edit_generate_input_lmp_files(lmp_infyle,lmp_datafyle):

    if not os.path.exists(lmp_infyle):
        raise RuntimeError('ERROR: ' + lmp_infyle + ' not found')

    if not os.path.exists(lmp_datafyle):
        raise RuntimeError('ERROR: ' + lmp_datafyle + ' not found')

    fr  = open(lmp_infyle,'r')
    fw  = open('in.init','w')
    fid = fr.read().replace("py_dataname",lmp_datafyle)
    fw.write(fid)
    fw.close()
    fr.close()
    

def gen_pair_coeff_file(destdir,ntypes,name_list,eps_list,sig_list,ljcut_list,coulcut_list):

    with open(destdir + '/pair_info_file.dat','w') as fpair:
        fpair.write('# LJ-Coul Pair coefficient file\n')
        fpair.write('# Add after read_data datafile in LAMMPS input file\n')
        fpair.write('pair_style      lj/cut/coul/long 1.122462 10.0\n')
        for attype in range(ntypes):
            fpair.write('%s\t %d %d %g %g %g #%s\n' %('pair_coeff', attype+1, attype+1,\
                                                      eps_list[attype],sig_list[attype],\
                                                      round(ljcut_list[attype],10),\
                                                      name_list[attype]))

def gen_bond_coeff_file(destdir,bname_list,kspr_list,sig_list,bcon_list):

    with open(destdir + '/bond_info_file.dat','w') as fbond:
        fbond.write('# Bond coefficient file\n')
        fbond.write('# Add after bond_style command in LAMMPS input file\n')
        for bttype in range(len(bcon_list)):
            eqbm_dst = 0.5*(sig_list[bcon_list[bttype][0]-1] + sig_list[bcon_list[bttype][1]-1])
            fbond.write('%s\t %d %g %g #%s\n' %('bond_coeff', bttype+1,kspr_list[bttype],\
                                                eqbm_dst,bname_list[bttype]))


def run_lammps(mw_chain,fr_anion,casenum,inpjob,outjob,tot_hrs,tot_nodes,tot_cores):

    if not os.path.exists(inpjob):
        raise RuntimeError('ERROR: ' + inpjob + ' not found')
    
    jobstr = "job_" + str(mw_chain) + "_" + str(fr_anion) + "_" \
             + str(casenum)
    fr  = open(inpjob,'r')
    fw  = open(outjob,'w')
    fid = fr.read().replace("py_jobname",jobstr).\
          replace("py_tottime",str(tot_hrs)).\
          replace("py_nnodes",str(tot_nodes)).\
          replace("py_ncores",str(tot_cores)).\
          replace("py_nptot",str(tot_cores*tot_nodes))
    fw.write(fid)
    fw.close()
    fr.close()

    subprocess.call(["sbatch", outjob])
    
    
def clean_backup_initfiles(f90_files,tcl_files,lmpinp1,lmpinp2,destdir):
    
    initdir = destdir + '/init_files'
    if not os.path.isdir(initdir):
        os.mkdir(initdir)

    for fyl in f90_files:
        if os.path.exists(fyl):
            cpy_main_files(destdir,initdir,fyl)
            os.remove(fyl)

    for fyl in tcl_files:
        if os.path.exists(fyl):
            cpy_main_files(destdir,initdir,fyl)
            os.remove(fyl)

    if os.path.exists(lmpinp1):
        cpy_main_files(destdir,initdir,lmpinp1)

    if os.path.exists(lmpinp2):
        cpy_main_files(destdir,initdir,lmpinp2)

    files = glob.glob('init_*')
    for fyl in files:
        if not os.path.isdir(fyl):
            cpy_main_files(destdir,initdir,fyl)
            os.remove(fyl)

    files = glob.glob('*var*')
    for fyl in files:
        if os.path.exists(fyl):
            cpy_main_files(destdir,initdir,fyl)
            os.remove(fyl)

    files = glob.glob('*.o')
    for fyl in files:
        if os.path.exists(fyl):
            cpy_main_files(destdir,initdir,fyl)
            os.remove(fyl)
    
    files = glob.glob(destdir +'/*.mod')
    for fyl in files:
        os.remove(fyl)


def compile_anafiles():

    if not os.path.exists('pe_params.f90'):
        print('ERROR: pe_params.f90 not found')
        return

    if not os.path.exists('pe_analyze.f90'):
        print('ERROR: pe_analyze.f90 not found')
        return


    subprocess.call(['ifort','-r8','-qopenmp','-mkl',\
                     '-check','-traceback','pe_params.f90',\
                     'pe_analyze.f90','-o','anainp.o'])

def find_datafyle(nch_free,archstr,casenum,destdir):

    os.chdir(destdir)
    curr_dir = os.getcwd()
    datafyle = "PEdata_"+str(nch_free)+"_" +archstr+ \
               "_"+str(casenum)+".dat"

    if not os.path.exists(datafyle):
        print ("Data file not found ..")
        print ("Making datafile from restart files")
        restart_fyles = glob.glob('archival_*')
        
        if restart_fyles == []:
            return 'ERROR'

        if not os.path.exists('lmp_mesabi'):
        
            src_lmp = '/home/dorfmank/vsethura/mylammps/src/lmp_mesabi'
            destfyle = curr_dir + '/lmp_mesabi'
            shutil.copy2(src_lmp,destfyle)

        subprocess.call(['mpirun','-np','48','./lmp_mesabi','-r',restart_fyles[0],datafyle])
        
    return datafyle


def find_latest_trajfyle(pref,destdir):
    
    os.chdir(destdir)
    traj_arr = glob.glob(pref)
    if traj_arr == []:
        return 'ERROR'
    latest_fyle = max(traj_arr,key = os.path.getctime)
    return latest_fyle


def edit_generate_anainp_files(inpdata,inptraj,nch_tot,nch_free,\
                               nch_graft,cutoff,listnum):
    
    if not os.path.exists('anainp_var.txt'):
        print('ERROR: pe_params not found')
        return

    fr  = open('anainp_var.txt','r')
    outana = 'anainp_' + str(listnum) + '.txt'
    fw  = open(outana,'w')

    datafyle = re.split('/',inpdata)
    dataname = datafyle[len(datafyle)-1]

    trajfyle = re.split('/',inptraj)
    trajname = trajfyle[len(trajfyle)-1]
    

    fid = fr.read().replace("py_datafyl",dataname).\
          replace("py_trajfyl",trajname).\
          replace("py_ntotchains",str(nch_tot)).\
          replace("py_nfrchains", str(nch_free)).\
          replace("py_ngrchains",str(nch_graft)).\
          replace("py_cutoff",str(cutoff))
    fw.write(fid)
    fw.close()
    fr.close()

def run_analysis(nch_free,pdifree,casenum,dirstr,inpjob,outjob,listnum,destdir):

    if not os.path.exists(inpjob):
        print('ERROR: ', inpjob,'not found')
        return
    
    jobstr = "ana_" + str(nch_free) + "_" + str(pdifree) + "_" \
             + str(casenum) + "_" + dirstr
    fr  = open(inpjob,'r')
    fw  = open(outjob,'w')
    fid = fr.read().replace("py_jobname",jobstr).\
          replace("py_freech",str(nch_free)).\
          replace("py_pdifree",str(pdifree)).\
          replace("py_caselen",str(casenum)).\
          replace("pyfylval",str(listnum)).\
          replace("pyoutdir",str(destdir))
    fw.write(fid)
    fw.close()
    fr.close()

    subprocess.call(["qsub", outjob])


def find_recent_file(destdir,keyword): #A replica of find_recent_traj_file

    fylnames = destdir + '/' + keyword
    list_of_files = glob.glob(fylnames)
    if list_of_files != []:
        fyl_str = max(list_of_files, key=os.path.getctime)
        fyl_arr = fyl_str.split("/")
        print( "File Name: ", fyl_arr[len(fyl_arr)-1])
        return fyl_arr[len(fyl_arr)-1]
    else:
        return "nil"
    
def my_cpy_generic(srcdir,destdir,inpfylname,destfylname):
    src_fyl  = srcdir  + '/' + inpfylname
    dest_fyl = destdir + '/' + destfylname
    shutil.copy2(src_fyl, dest_fyl)

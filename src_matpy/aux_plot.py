# Auxiliary data for generic plot
#------------------------------------------------------------------
# Make labels and graphics
# Version: Mar-04-2025
#------------------------------------------------------------------

# Import modules
import os
import sys
import numpy as np
import re
import shutil
import glob
import math
import fileinput
import subprocess
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib.ticker as tck
from matplotlib.ticker import LogFormatterSciNotation

#------------------------------------------------------------------

# General copy script
def gencpy(dum_maindir,dum_destdir,fylname):

    srcfyl = dum_maindir + '/' + fylname

    if not os.path.exists(srcfyl):
        print('ERROR: ', srcfyl, 'not found')
        return

    desfyl = dum_destdir + '/' + fylname
    shutil.copy2(srcfyl,desfyl)
#------------------------------------------------------------------

def ret_temp_dir(scr_dir,inp_type,biomass,dispval,casenum,subdir,\
                 solv_type='None'):
    # Check all directories
    head_dir = scr_dir + '/' + inp_type
    if not os.path.isdir(head_dir):
        raise RuntimeError(head_dir + " does not exist!")
       
    poly_dir = head_dir + '/' + biomass
    if not os.path.isdir(poly_dir):
        raise RuntimeError(poly_dir + " does not exist!")
        
    if inp_type == 'melts':
        fig_dir  = poly_dir + '/results_figs'
        if not os.path.isdir(fig_dir):
            os.mkdir(fig_dir)

        if dispval == 'expts':
            poly_dir = poly_dir + '/expts'
        else:
            poly_dir = poly_dir + '/pdi_' + str(dispval)
        if not os.path.isdir(poly_dir):
            raise RuntimeError(poly_dir + " does not exist!")

    rundir = poly_dir + '/run_' + str(casenum)
    if not os.path.isdir(rundir):
        raise RuntimeError(rundir + " does not exist!")

    if inp_type == 'solvents':
        workdir1 = rundir + '/' + solv_type
    elif inp_type == 'cosolvents':
        h_workdir = rundir + '/' + solv_type
        workdir1 = h_workdir + '/' + solv_type + '_water'
    else:
        workdir1 = rundir #even for inp_type = melts

    if not os.path.isdir(workdir1):
        raise RuntimeError(workdir1, " does not exist!")

    anadir = workdir1 + '/T_' + str(subdir)
    if not os.path.isdir(anadir):
        raise RuntimeError(anadir, " does not exist!")

    return workdir1, anadir, fig_dir
#------------------------------------------------------------------

def ret_ana_dir(scr_dir,inp_type,biomass,dispval,casenum,subdir,\
                 solv_type='None'):
    # Check all directories
    head_dir = scr_dir + '/' + inp_type
    if not os.path.isdir(head_dir):
        raise RuntimeError(head_dir + " does not exist!")

    fig_dir  = head_dir + '/results_figs'
    if not os.path.isdir(fig_dir):
        os.mkdir(fig_dir)
        
    poly_dir = head_dir + '/' + biomass
    if not os.path.isdir(poly_dir):
        raise RuntimeError(poly_dir + " does not exist!")
        
    if inp_type == 'melts':
        poly_dir = poly_dir + '/pdi_' + str(dispval)
        if not os.path.isdir(poly_dir):
            raise RuntimeError(poly_dir + " does not exist!")

    rundir = poly_dir + '/run_' + str(casenum)
    if not os.path.isdir(rundir):
        raise RuntimeError(rundir + " does not exist!")

    if inp_type == 'solvents':
        workdir1 = rundir + '/' + solv_type
    elif inp_type == 'cosolvents':
        h_workdir = rundir + '/' + solv_type
        workdir1 = h_workdir + '/' + solv_type + '_water'
    else:
        workdir1 = rundir #even for inp_type = melts

    if not os.path.isdir(workdir1):
        raise RuntimeError(workdir1, " does not exist!")

    anadir = workdir1 + '/' + subdir
    if not os.path.isdir(workdir1):
        raise RuntimeError(anadir, " does not exist!")

    return workdir1, anadir, fig_dir
#------------------------------------------------------------------        
# Change width of seaborn barplot
def change_width(ax, new_value):
    for patch in ax.patches:
        current_width = patch.get_width()
        diff = current_width - new_value

        # Change the bar width
        patch.set_width(new_value)

        # Recenter the bar
        patch.set_x(patch.get_x() + diff * .5)
#------------------------------------------------------------------
# Return Chain ID from file name
def ret_chid(inpfyle):
    outstr = re.split(r'/|.xvg|_',inpfyle)
    chid = outstr[len(outstr)-2]
    if not chid.isdigit():
        print('ERR: not a numeric value', chid)
    else:
        return int(chid)
#------------------------------------------------------------------
# Return the monomers in a chain
def ret_mons(inpfyle):
    if not os.path.exists(inpfyle):
        raise RuntimeError('File not found ' + inpfyle)
    mon_arr = []
    with open(inpfyle) as finp:
        finp.readline()
        for line in finp:
            line = line.rstrip('/n')
            outstr = line.split()
            mon_arr.append(int(outstr[len(outstr)-2]))
    return mon_arr
#------------------------------------------------------------------
# All Rg plots
def plot_allrg(df,fig_dir,pdi_val):
    #Plot Rg^2 - T
    figa, axa = plt.subplots()
    axa.scatter(x=df['Temp'],y=df['<Rg^2>'])
    set_axes(axa,plt,r'Temperature ($K$)',r'$R_{g}^{2}$ ($nm^{2}$)') 
    figa.savefig(fig_dir + '/'+'rg2_' + str(pdi_val) + '.png',\
                 dpi=figa.dpi)
    plt.close(figa)
    
    #Plot Rg^4 - T
    figa, axa = plt.subplots()
    axa.scatter(x=df['Temp'],y=df['<Rg^4>'])
    set_axes(axa,plt,r'Temperature ($K$)',r'$R_{g}^{4}$ ($nm^{4}$)') 
    figa.savefig(fig_dir + '/'+'rg4_' + str(pdi_val) + '.png',\
                 dpi=figa.dpi)
    plt.close(figa)
    
    #Plot Alpha - T
    figa, axa = plt.subplots()
    axa.scatter(x=df['Temp'],y=df['Rg4/Rg2^2'])
    set_axes(axa,plt,r'Temperature ($K$)',\
             r'$\langle R_{g}^{4} \rangle/ \langle R_{g}^{2} \rangle^{2}$') 
    figa.savefig(fig_dir + '/'+'alpha' + str(pdi_val) + '.png',\
                 dpi=figa.dpi)
    plt.close(figa)
#------------------------------------------------------------------
# Plot Tg
def plot_tg(df,fig_dir,pdi_val):
    figa, axa = plt.subplots()
    set_axes(axa,plt,r'Temperature ($K$)',r'Specific Volume ($cm^3/g$)')
    axa.scatter(x=df['Temp'],y=df['SV_NPT'])
    figa.savefig(fig_dir + '/'+'sv_' + str(pdi_val) + '.png',\
                 dpi=figa.dpi)
    return figa, axa
#------------------------------------------------------------------
# All shape plots
def plot_allshape(df,fig_dir,pdi_val):
    #Plot lam_x - T
    figa, axa = plt.subplots()
    axa.scatter(x=df['Temp'],y=df['<lam_1>'])
    set_axes(axa,plt,r'Temperature ($K$)',r'$\lambda_{x}$ ($nm^{2}$)') 
    figa.savefig(fig_dir + '/'+'lamx_' + str(pdi_val) + '.png',\
                 dpi=figa.dpi)
    plt.close(figa)

    #Plot lam_y - T
    figa, axa = plt.subplots()
    axa.scatter(x=df['Temp'],y=df['<lam_2>'])
    set_axes(axa,plt,r'Temperature ($K$)',r'$\lambda_{y}$ ($nm^{2}$)') 
    figa.savefig(fig_dir + '/'+'lamy_' + str(pdi_val) + '.png',\
                 dpi=figa.dpi)
    plt.close(figa)

    #Plot lam_z - T
    figa, axa = plt.subplots()
    axa.scatter(x=df['Temp'],y=df['<lam_3>'])
    set_axes(axa,plt,r'Temperature ($K$)',r'$\lambda_{z}$ ($nm^{2}$)') 
    figa.savefig(fig_dir + '/'+'lamz_' + str(pdi_val) + '.png',\
                 dpi=figa.dpi)
    plt.close(figa)

    #Plot kappa - T
    figa, axa = plt.subplots()
    axa.scatter(x=df['Temp'],y=df['<\kappa>'])
    set_axes(axa,plt,r'Temperature ($K$)',r'$\kappa$ ($nm^{4}$)') 
    figa.savefig(fig_dir + '/'+'kappa_' + str(pdi_val) + '.png',\
                 dpi=figa.dpi)
    plt.close(figa)
#------------------------------------------------------------------
# Plot b-nu data
def plot_rgscaling(df,fig_dir,pdi_val):
    #Plot b - T
    figa, axa = plt.subplots()
    axa.scatter(x=df['Temp'],y=df['b'])
    set_axes(axa,plt,r'Temperature ($K$)',r'$C_{l}$ ($nm$)') 
    figa.savefig(fig_dir + '/'+'lamx_' + str(pdi_val) + '.png',\
                 dpi=figa.dpi)
    plt.close(figa)

    #Plot nu-T
    figa, axa = plt.subplots()
    axa.scatter(x=df['Temp'],y=df['nu'])
    set_axes(axa,plt,r'Temperature ($K$)',r'$\nu$') 
    figa.savefig(fig_dir + '/'+'lamx_' + str(pdi_val) + '.png',\
                 dpi=figa.dpi)
    plt.close(figa)
#------------------------------------------------------------------
# All shape plots
def plot_allhb(df,fig_dir,pdi_val):
    #Plot intra HB - T
    figa, axa = plt.subplots()
    axa.scatter(x=df['Temp'],y=df['Intra_HB'])
    set_axes(axa,plt,r'Temperature ($K$)',r'#Intra HB') 
    figa.savefig(fig_dir + '/'+'hbintra_' + str(pdi_val) + '.png',\
                 dpi=figa.dpi)
    plt.close(figa)

    #Plot inter HB - T
    figa, axa = plt.subplots()
    axa.scatter(x=df['Temp'],y=df['Inter_HB'])
    set_axes(axa,plt,r'Temperature ($K$)',r'Inter_HB') 
    figa.savefig(fig_dir + '/'+'hbinter_' + str(pdi_val) + '.png',\
                 dpi=figa.dpi)
    plt.close(figa)

    #Plot tot HB - T
    figa, axa = plt.subplots()
    axa.scatter(x=df['Temp'],y=df['Tot_HB'])
    set_axes(axa,plt,r'Temperature ($K$)',r'#HB') 
    figa.savefig(fig_dir + '/'+'hbtot_' + str(pdi_val) + '.png',\
                 dpi=figa.dpi)
    plt.close(figa)
#------------------------------------------------------------------    
# Plot axis details
def set_axes(axhdl,plt,xlabel,ylabel):
    plt.style.use('seaborn-colorblind')
    plt.tight_layout()
    axhdl.set_xlabel(xlabel)
    axhdl.set_ylabel(ylabel)
    change_width(axhdl,0.2)
#------------------------------------------------------------------
# find axes limits
def axlims(yminref,ymin,ymaxref,ymax):
    if ymax > ymaxref: ymaxref = ymax
    if ymin < yminref: yminref = ymin
    return yminref, ymaxref
#------------------------------------------------------------------
# Submit density job script
def submit_dens(currdir,workdir,f_name,tmin,tmax,tdt):
    os.chdir(workdir)
    rev_fname = f_name.replace('_pyinp','')
    fr  = open(workdir + '/' + f_name,'r')
    fw  = open(workdir + '/' + rev_fname,'w')
    fid = fr.read().replace("py_Tinit",str(tmin)).\
          replace("py_Tfin",str(tmax)).\
          replace("py_dT",str(tdt)).\
          replace("py_jobname","denscomp")
    fw.write(fid)
    fw.close()
    subprocess.call(["sbatch",rev_fname])
    os.chdir(currdir)
#------------------------------------------------------------------     
def set_axes(axhdl,plt,xlabel,ylabel):
    plt.style.use('seaborn-colorblind')
    plt.tight_layout()
    axhdl.set_xlabel(xlabel)
    axhdl.set_ylabel(ylabel)
    change_width(axhdl,0.2)
#------------------------------------------------------------------     
def return_neigh_arrays(simdir,strpref,colarr,nheaders=0,maxneigh=5):
    if len(colarr) != 2:
        raise RuntimeError("Expecting 2 columns for plotting neighbors")

    neigh_file = find_latest_file(glob.glob(simdir + '/' + strpref + '*'))
    if neigh_file == -1:
        return -1

    sarr,neigharr = extract_data(neigh_file,colarr,nheaders,maxneigh)
    return sarr,neigharr
#------------------------------------------------------------------         
def return_clust_arrays(simdir,strpref,colarr,nheaders=0,maxclust=20):
    if len(colarr) != 2:
        raise RuntimeError("Expecting 2 columns for plotting clusters")

    neigh_file = find_latest_file(glob.glob(simdir + '/' + strpref + '*'))
    if neigh_file == -1:
        return -1

    sarr,neigharr = extract_data(neigh_file,colarr,nheaders,maxclust)
    return sarr,neigharr
#------------------------------------------------------------------             
def extract_data(filename,colarr,skipl=1,refnr=0):
    data = np.loadtxt(filename, skiprows=skipl)  # Skips n header lines
    nrows = data.shape[0]
    if refnr == 0 or nrows < refnr:
        col1 = data[:, colarr[0]]  # Column 1
        col2 = data[:, colarr[1]]  # Column n
    else:
        col1 = data[:refnr, colarr[0]]  # Column 1
        col2 = data[:refnr, colarr[1]]  # Column n       
    return col1, col2
#------------------------------------------------------------------         
def find_latest_file(flist):
    if not len(flist):
        print('ERROR: No files found in ' + simdir)
        return -1
    outfile = max(flist, key=os.path.getmtime)
    return outfile
#------------------------------------------------------------------         
def sci_mathtext(x, pos):
    # format as $2\times10^{-12}$, $10^{-11}$, etc.
    if x <= 0:
        return ""
    exp = int(np.floor(np.log10(x)))
    coeff = x / 10**exp
    if abs(coeff - 1.0) < 1e-12:
        return rf"$10^{{{exp}}}$"
    else:
        # %g avoids trailing .0; switch to:.3g if you want fixed precision
        return rf"${coeff:g}\times10^{{{exp}}}$"
#------------------------------------------------------------------         
# if __name__
if __name__ == '__main__':
    main()
#------------------------------------------------------------------

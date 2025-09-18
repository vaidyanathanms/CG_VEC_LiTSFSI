# Computing different properties
# V2: Mar-04-2025
# V1: Oct-12-2021
# Current capabilities: Tg, Rgscaling, Bin averages, Nearest MSD
#-------------------------------------------------------------------

#Import modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.linear_model import LinearRegression
from sklearn import metrics
from scipy.optimize import curve_fit
from scipy.stats import sem
import os
import sys
import re
import shutil
import glob
import math
import subprocess
#------------------------------------------------------------------

# Compute Tg
def compute_tg(df,ax):

    tvals = np.array(df['Temp']); den_vals = np.array(df['SV_NPT'])
    wvals = np.array(df['SV_err']); wvals = 1/wvals
    swvals = np.sum(wvals); wvals /= swvals

    lden  = len(tvals); 
    if lden < 7:
        print("ERR:Requires at least 7 temperature points to fit")
        return 0

    
    clf = LinearRegression(fit_intercept = True) 
    toterr = 9.0*10**20; reflm = 0; reflc = 0; refhm = 0; refhc = 0
    refle = 0; refhe = 0

    for cnt in range(2,lden-1):

        low_fit = clf.fit(tvals[0:cnt].reshape(-1,1),den_vals[0:cnt],\
                          sample_weight = wvals[0:cnt])
        lc = low_fit.intercept_; lm = low_fit.coef_
        predlo = clf.predict(tvals[0:cnt].reshape(-1,1))
        errlo = np.sum(np.multiply(wvals[0:cnt],np.square((den_vals[0:cnt]-predlo))))
        
        hi_fit = clf.fit(tvals[cnt:lden].reshape(-1,1),den_vals[cnt:lden],\
                         sample_weight = wvals[cnt:lden])
        predhi = clf.predict(tvals[cnt:lden].reshape(-1,1))
        hc = low_fit.intercept_; hm = low_fit.coef_
        errhi = np.sum(np.multiply(wvals[cnt:lden],np.square((den_vals[cnt:lden]-predhi))))

        if errlo + errhi < toterr:
            refcnt = cnt; toterr = errlo + errhi
            reflm = lm; reflc = lc; refhm = hm; refhc = hc
            refle = errlo; refhe = errhi
            
    # Solve linear equations y-m1x = c1 & y-m2x = c2
    lhs = np.array([[1,-reflm[0]], [1,-refhm[0]]])
    rhs = np.array([reflc,refhc])
    tg  = np.linalg.solve(lhs,rhs)
    return tg[1],reflc, reflm, refle, refhc, refhm, refhe
#------------------------------------------------------------------ 

# Compute Rgscaling
def compute_rgscaling(df):
    N_all = np.array(df['N']); Rg_all = np.array(df['Rgmean'])
    lenarr = min(math.floor(0.5*len(N_all)),10)
    if lenarr == 0:
        print("ERR: No data found")
        return [float('NaN'),float('NaN')], float('NaN')
    xarr = N_all[:lenarr]; yarr = Rg_all[:lenarr]
    pars, cov = curve_fit(func_powerlaw,xarr,yarr,p0=np.asarray([10,0.5]))
    stdevs = np.sqrt(np.diag(cov))
    return pars, cov
#------------------------------------------------------------------ 

# Define power law
def func_powerlaw(x,a,b):
    return a*x**b
#------------------------------------------------------------------ 

# Compute theoretical SZ distribution
def comp_pdi_sz(pdi_val,sig,mwn):
    pinv = (1.0)/(pdi_val-1.0) #defn = 1/(PDI-1)
    term1 = pinv**pinv; term2 = sig**(pinv-1)
    term3 = np.exp(-pinv*sig); term4 = math.gamma(pinv)
    psz_theory = (term1*term2*term3)/term4
    return psz_theory, np.trapz(psz_theory,sig*mwn)
#------------------------------------------------------------------

# Return averages over bins
def compute_bin_averages(xdata,ydata):

    xoutarr = np.empty(0); ysumarr = np.empty(0);
    cntarr  = np.empty(0)

    if len(xdata) == [] or len(ydata) == []:
        print('Error: No Data Present\n')
        return xoutarr,ysumarr,cntarr,ysumarr
        
    for i in range(len(xdata)):
        xval = xdata[i]; yval = ydata[i]
        if math.isnan(xval) or math.isnan(yval):
            continue
        elif xval in xoutarr: 
            ysumarr[list(xoutarr).index(xval)] += yval
            cntarr[list(xoutarr).index(xval)]  += 1
        else:
            xoutarr = np.append(xoutarr,xval)
            ysumarr = np.append(ysumarr,yval)
            cntarr  = np.append(cntarr,1)

    return xoutarr, ysumarr/cntarr, cntarr, ysumarr
#------------------------------------------------------------------ 

# Return averages and SEM over bins
def comp_bin_ave_sem(xdata,ydata):

    xoutarr = []; yallarr = [] 

    if len(xdata) == [] or len(ydata) == []:
        print('Error: No Data Present\n')
        return [],[],[]
        
    for i in range(len(xdata)):
        xval = xdata[i]; yval = ydata[i]
        if math.isnan(xval) or math.isnan(yval):
            continue
        elif xval in xoutarr:
            index = list(xoutarr).index(xval)
            yallarr[index].append(yval) #as list
        else:
            xoutarr.append(xval) #just an element
            yallarr.append([yval]) #as list

    yavearr = np.empty(0); ysemarr = np.empty(0); cntarr = np.empty(0)
    for i in range(len(xoutarr)):
        yavearr = np.append(yavearr,np.mean(np.array(yallarr[i])))
        cntarr  = np.append(cntarr,np.size(np.array(yallarr[i])))
        if np.size(np.array(yallarr[i])) == 1: #if only one point is
            #there to average, SEM will be nan
            ysemarr = np.append(ysemarr,0)
        else:
            ysemarr = np.append(ysemarr,sem(np.array(yallarr[i])))
    return xoutarr, yavearr, ysemarr, cntarr
#------------------------------------------------------------------ 

# Return closest MSD at given tref
def msd_tref(data,reft):
    if data[len(data[:,0])-2,0] < reft:
        print('Simulation time less than reference time')
        print(data[len(data[:,0])-2,0],reft)
        return -1
    else:
        for tind in range(len(data[:,0])):
            if data[tind,0] < reft and data[tind+1,0] >= reft:
                sl = (data[tind+1,1]-data[tind,1])/\
                    (data[tind+1,0]-data[tind,0])
                return data[tind,1] + sl*(reft-data[tind,0])
#------------------------------------------------------------------ 



# if __name__
if __name__ == '__main__':
    main()
#------------------------------------------------------------------

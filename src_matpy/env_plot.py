# Environment variables for plotting

import matplotlib.pyplot as plt

# Color/line data
orange = '#FFA500'; dark_g = '#006400'; d_brown = '#8B4513'
lite_g = '#7FC97F'; violet = '#BEAED4'; lite_o = 'FDC086'
clr_arr = [violet,orange,dark_g,'m',lite_o,'b',d_brown]
mrk_arr = ['o','d','s','v','H','>']
lne_arr = ['-','--']
#------------------------------------------------------------------

# Default plot data
plt.rcParams['figure.dpi'] = 1200
plt.rc('legend',fontsize=12) # fontsize of the legends
plt.rcParams.update({'font.size': 16}) # other fonts
plt.rcParams.update({'figure.autolayout': True})
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.it'] = 'Arial:italic'
plt.rcParams['mathtext.rm'] = 'Arial'
#plt.rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
plt.rcParams['font.family'] = 'Arial' # font type
plt.rcParams['lines.markersize'] = 4 # marker size
plt.rcParams['lines.linewidth'] = 2 # line width
#plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
#plt.rcParams.update({'font.family': 'Times New Roman'})
#------------------------------------------------------------------

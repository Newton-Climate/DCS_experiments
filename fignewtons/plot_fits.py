import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib
import seaborn as sns
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc
rc('text', usetex=False)
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
from setup_retrieval_plots import *
from process_results import *
import logging
logging.getLogger('matplotlib.font_manager').disabled = True

print('Running em_plot_fits.py...')


font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 8}

matplotlib.rc('font', **font)

cp = sns.color_palette(palette="deep")

c = {'pic' : 'black', ### CHANGED TO PIC FROM PICARRO
    'hit08' : cp[0],
    'hit16' : cp[1], ### THIS MIGHT BE A MISTAKE?
    'hit16' : cp[2],
    'hit20' : cp[3],
    'tccon' : cp[4],
    'oco' : cp[5]}

markersize = 1

# on Fluo
datadir = '/net/fluo/data1/data/NIST/DCS_A/'

## read the data
#retrievals 
hit08 = read_data(datadir+'timeseries/data/hit08_results.jld2')
hit16 = read_data(datadir+'timeseries/data/hit16_results.jld2')
hit20 = read_data(datadir+'timeseries/data/hit20_results.jld2')
tccon = read_data(datadir+'timeseries/data/tccon_results.jld2')
OCO = read_OCO(datadir+'timeseries/data/OCO_results.jld2')


ch4_grid = np.arange(start=6050.0, stop=6120, step=(6120 - 6050)/len(hit08['ch4_measurement']))
co2_grid = np.arange(start=6180.0, stop=6260, step=(6260 - 6180)/len(hit08['co2_measurement']))


######### TOP SUBPLOT: CH4 Spectra #########
key = 'ch4_model'
obs = hit08['ch4_measurement']

fig1, (f1_ax1, f1_ax2)  = plt.subplots(nrows=2,ncols=1, sharex = False, figsize = (6,7), dpi = 300)

# CH4 spectra
f1_ax1.plot(ch4_grid, obs, '-', label='measurement', c='black', markersize = 2, linewidth = 0.5)
f1_ax1.set_xlabel(r'cm$^{-1}$')
f1_ax1.set_ylabel('Intensity')
f1_ax1.set_title(r'Fit over CH$_4$ window')
xtickslocs = f1_ax1.get_xticks()

######### BOTTOM SUBPLOT #########
gs = gridspec.GridSpec(2, 1)
f1_ax2 = plt.subplot(gs[1])
axMain = f1_ax2
plt.sca(axMain) #sca = set current axes

#all the axes, axMain is the bottom
divider = make_axes_locatable(axMain)
axShallow2 = divider.append_axes("top", size="100%", pad=0.3, sharex=axMain) #middle, needs to have cf = in front
axShallow3 = divider.append_axes("top", size="100%", pad=0.3, sharex=axMain) #top
###axShallow4 = divider.append_axes("top", size="100%", pad=0.1, sharex=axMain) #top

# residuals (model - measurement)
linewidth = 0.5#formerly 0.2
axMain.set_xlabel(r'Wavenumbers (cm$^{-1}$)')
axShallow2.set_ylabel('Model - Measurement')
axShallow2.plot(ch4_grid, obs - hit08[key], label='Hitran-16', c=c['hit16'], alpha = 1, markersize = 2, linewidth =linewidth)
axShallow3.plot(ch4_grid, obs - hit16[key], label='Hitran-08', c=c['hit08'], alpha = 1, markersize = 2, linewidth =linewidth)
axMain.plot(ch4_grid, obs - hit20[key], label='Hitran-20', c=c['hit20'], alpha = 1, markersize = 2, linewidth =linewidth)
f1_ax2.set_xlabel(r'cm$^{-1}$')

# Set labels and ylims and get rid of spines
axShallow2.tick_params(labelbottom=False, axis = 'x', colors = 'white')
axShallow3.tick_params(labelbottom=False, axis = 'x', colors = 'white')

offset_label = 0.0025 #6077
axMain.text(x = 6050, y = (obs-hit08[key]).mean() + offset_label, s = 'Hitran-2020')#, fontsize=10)
axShallow2.text(x = 6050, y = (obs-hit08[key]).mean() + offset_label, s = 'Hitran-2016')#, fontsize=10)
axShallow3.text(x = 6050, y = (obs-hit08[key]).mean() + offset_label, s = 'Hitran-2008')#, fontsize=10)

axShallow2.spines['top'].set_visible(False)
axShallow3.spines['bottom'].set_visible(False)
axMain.spines['top'].set_visible(False)
axShallow2.spines['bottom'].set_visible(False)

ymin = -0.005
ymax = 0.005
axMain.set_ylim([ymin, ymax])
axShallow2.set_ylim([ymin, ymax])
axShallow3.set_ylim([ymin, ymax])

plt.tight_layout()
plt.savefig(datadir+'figs/ch4_fit.pdf')
plt.close()
print('\tSAVED: figs/ch4_fit.pdf')





############################### Co2 ###############################
###############################

#fig1, (f1_ax0, f1_ax1, f1_ax2, f1_ax3, f1_ax4) = plt.subplots(nrows=5, sharex=True,figsize=[6, 9], dpi=300)            
#ymin = -0.005
#ymax = 0.005

key = 'co2_model'
obs = hit16['co2_measurement']

fig1, (f1_ax1, f1_ax2)  = plt.subplots(nrows=2,ncols=1, sharex = False, figsize = (6,7), dpi = 300)

######### TOP SUBPLOT: CO2 Spectra #########

# CH4 spectra
f1_ax1.plot(co2_grid, obs, '-', label='measurement', c='black', markersize = 2, linewidth = 0.5)
f1_ax1.set_xlabel(r'cm$^{-1}$')
f1_ax1.set_ylabel('Intensity')
f1_ax1.set_title(r'Fit over CH$_4$ window')
xtickslocs = f1_ax1.get_xticks()

######### BOTTOM SUBPLOT #########
gs = gridspec.GridSpec(2, 1)
f1_ax2 = plt.subplot(gs[1])
axMain = f1_ax2
plt.sca(axMain) #sca = set current axes

#all the axes, axMain is the bottom
divider = make_axes_locatable(axMain)
axShallow2 = divider.append_axes("top", size="100%", pad=0.3, sharex=axMain) #middle, needs to have cf = in front
axShallow3 = divider.append_axes("top", size="100%", pad=0.3, sharex=axMain) #top
axShallow4 = divider.append_axes("top", size="100%", pad=0.3, sharex=axMain) #top

# residuals (model - measurement)
linewidth = 0.5#formerly 0.2
axMain.set_xlabel(r'Wavenumbers (cm$^{-1}$)')
axShallow3.set_ylabel('Model - Measurement')
axShallow2.plot(co2_grid, obs - hit08[key], label='Hitran-16', c=c['hit16'], alpha = 1, markersize = 2, linewidth =linewidth)
axShallow3.plot(co2_grid, obs - hit16[key], label='Hitran-08', c=c['hit08'], alpha = 1, markersize = 2, linewidth =linewidth)
axShallow4.plot(co2_grid, obs - tccon[key], label='TCCON', c=c['tccon'], alpha = 1, markersize = 2, linewidth =linewidth)

axMain.plot(co2_grid, obs - hit20[key], label='Hitran-20', c=c['hit20'], alpha = 1, markersize = 2, linewidth =linewidth)
f1_ax2.set_xlabel(r'cm$^{-1}$')

# Set labels and ylims and get rid of spines
axShallow2.tick_params(labelbottom=False, axis = 'x', colors = 'white')
axShallow3.tick_params(labelbottom=False, axis = 'x', colors = 'white')
axShallow4.tick_params(labelbottom=False, axis = 'x', colors = 'white')

offset_label = -0.0035 #6077
axMain.text(x = 6180, y = (obs-hit08[key]).mean() + offset_label, s = 'Hitran-2020', fontsize=10)
axShallow2.text(x = 6180, y = offset_label, s = 'Hitran-2016', fontsize=10)
axShallow3.text(x = 6180, y = offset_label, s = 'Hitran-2008', fontsize=10)
axShallow4.text(x = 6180, y = offset_label, s = 'TCCON', fontsize=10)


axShallow2.spines['top'].set_visible(False)
axShallow3.spines['bottom'].set_visible(False)
axShallow4.spines['bottom'].set_visible(False)
axMain.spines['top'].set_visible(False)
axShallow2.spines['bottom'].set_visible(False)
axShallow3.spines['top'].set_visible(False)



ymin = -0.005
ymax = 0.005
axMain.set_ylim([ymin, ymax])
axShallow2.set_ylim([ymin, ymax])
axShallow3.set_ylim([ymin, ymax])
axShallow4.set_ylim([ymin, ymax])

plt.tight_layout()
plt.savefig(datadir+'figs/co2_fit.pdf')
plt.close()
print('\tSAVED: figs/co2_fit.pdf')



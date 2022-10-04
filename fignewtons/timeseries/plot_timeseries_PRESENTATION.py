from setup_retrieval_plots import *
from process_results import *
             

###################
plt.rcParams.update({'font.size': 14})
plt.rcParams['font.family'] = "serif"
#plt.rcParams['font.serif'] = "Times New Roman"
plt.rcParams['axes.grid'] = False
#plt.rcParams['axes.axisbelow'] = True
plt.rcParams['axes.labelpad'] = 6
figure_export_settings = {'dpi': 300, 'bbox_inches': 'tight'}
#####################

print('Running em_plot_timeseries.py for PRESENTATIONS...')

# How many days to plot??
PRINT_ALL = False # print all 4 subplots
change_xlim = True
x_max = 7


# on Fluo
datadir = '/net/fluo/data1/data/NIST/DCS_A/'

## read the data
#retrievals 
hit08 = read_data(datadir+'timeseries/data/hit08_results.jld2')
hit16 = read_data(datadir+'timeseries/data/hit16_results.jld2')
hit20 = read_data(datadir+'timeseries/data/hit20_results.jld2')
tccon = read_data(datadir+'timeseries/data/tccon_results.jld2')
OCO = read_OCO(datadir+'timeseries/data/OCO_results.jld2')

# observations
obs = read_obs(datadir+'timeseries/data/obs.jld2') # piccaro and tower pT data 
t = obs['t_pT'] 


print('\n\nOBS FILE WAS NOT UPDATED AS OF THE 12TH!\n')

# replace t, because time-stamps are inconsistant 
t = np.arange(start=0, stop=t[-1], step=(t[-1] - t[0])/len(t))

#################################### CH4 ########################
if PRINT_ALL == True:
    fig1, (f1_ax1, f1_ax2, f1_ax3, f1_ax4) = plt.subplots(nrows=4, sharex=True,figsize=[6, 9], dpi=300)            
else:
    fig1, (f1_ax1, f1_ax2) = plt.subplots(nrows=2, sharex=True,figsize=[6, 4], dpi=300)            

    
# CH4 concentrations
key = 'CH4_VMR'
f1_ax1.plot(t, hit08[key], c=c['hit08'], label='Hitran-08')
f1_ax1.plot(t, hit16[key], c=c['hit16'], label='Hitran-16')
f1_ax1.plot(t, hit20[key], c=c['hit20'], label='Hitran-20')
f1_ax1.plot(t, tccon[key], c=c['tccon'], label='TCCON')
f1_ax1.plot(obs['t_pic'], obs[key], 'o', markersize = 1, c=c['obs'], label='observation')
#f1_ax1.set_xlabel('days')
f1_ax1.set_ylabel('Concentration\n(ppb)')
#f1_ax1.legend(ncol = 2, frameon=False, loc='upper right')
if change_xlim: f1_ax1.set_xlim(0, x_max)
f1_ax1.set_ylim(1800, 2300)
f1_ax1.set_title('Retrieved $CH_4$')

# plot the column
key = 'CH4'
f1_ax2.plot(t, hit08[key], c=c['hit08'], label='Hitran-08')
f1_ax2.plot(t, hit16[key], c=c['hit16'], label='Hitran-16')
f1_ax2.plot(t, hit20[key], c=c['hit20'], label='Hitran-20')
f1_ax2.plot(t, tccon[key], c=c['tccon'], label='TCCON')
#f1_ax2.set_xlabel('days')
f1_ax2.set_ylim(0.7e19, 0.9e19)
f1_ax2.set_ylabel(r'$\frac{{\mathrm{molecules}}}{{\mathrm{cm}^2}}$')
if change_xlim: f1_ax2.set_xlim(0, x_max)

if PRINT_ALL == True:
    # plot retrieved pressure 
    key = 'p_ch4'
    f1_ax3.plot(t, hit08[key], c=c['hit08'], label='Hitran-08')
    f1_ax3.plot(t, hit16[key], c=c['hit16'], label='Hitran-16')
    f1_ax3.plot(t, hit20[key], c=c['hit20'], label='Hitran-20')
    f1_ax3.plot(t, tccon[key], c=c['tccon'], label='TCCON')
    f1_ax3.plot(t, obs['p'], 'o', markersize = 1, c=c['obs'], label='observation')
    #f1_ax3.set_xlabel('days')
    f1_ax3.set_ylabel('HPa')
    f1_ax3.set_ylim(780, 890)
    if change_xlim: f1_ax3.set_xlim(0, x_max)

    # plot retrieved temperature 
    key = 'T_ch4'
    f1_ax4.plot(t, hit08[key], c=c['hit08'], label='Hitran-08')
    f1_ax4.plot(t, hit16[key], c=c['hit16'], label='Hitran-16')
    f1_ax4.plot(t, hit20[key], c=c['hit20'], label='Hitran-20')
    f1_ax4.plot(t, tccon[key], c=c['tccon'], label='TCCON')
    f1_ax4.plot(t, obs['T'], 'o', markersize = 1, c=c['obs'], label='observation')

    f1_ax4.set_xlabel('days')
    f1_ax4.set_ylabel('Kelvin')
    f1_ax4.set_ylim(270,305)
    if change_xlim: f1_ax4.set_xlim(0, x_max)
else:
    f1_ax2.set_xlabel('days')
    
fig1.tight_layout()
fig1.savefig(datadir+'figs/presentation/retrieved_DCS_CH4_pres.pdf')
print('\tSAVED: figs/presentation/retrieved_DCS_CH4_pres.pdf')


#################################### C02 ########################
if PRINT_ALL == True:
    fig2, (f2_ax1, f2_ax2, f2_ax3, f2_ax4) = plt.subplots(nrows=4, sharex=True,figsize=[6, 9], dpi=300)            
else:
    fig2, (f2_ax1, f2_ax2) = plt.subplots(nrows=2, sharex=True,figsize=[6, 4], dpi=300)            

    
key = 'CO2_VMR'
f2_ax1.plot(t, hit16[key], c=c['hit16'], label='Hitran-16')
f2_ax1.plot(t, hit20[key], c=c['hit20'], label='Hitran-20')
f2_ax1.plot(t, tccon[key], c=c['tccon'], label='TCCON')
f2_ax1.plot(t, OCO[key], c=c['OCO'], label='OCO')
f2_ax1.plot(obs['t_pic'], obs[key], 'o', markersize = 1, c=c['obs'], label='observation')

#f2_ax1.set_xlabel('days')
f2_ax1.set_ylabel('Concentration\n(ppm)')
#f2_ax1.legend(ncol = 2, frameon=False, loc='upper left')
f2_ax1.set_ylim(385, 475)
if change_xlim: f2_ax1.set_xlim(0, x_max)
f2_ax1.set_title('Retrieved $CO_2$')


# plot the column
key = 'CO2'
f2_ax2.plot(t, hit16[key], c=c['hit16'], label='Hitran-16')
f2_ax2.plot(t, hit20[key], c=c['hit20'], label='Hitran-20')
f2_ax2.plot(t, tccon[key], c=c['tccon'], label='TCCON')
f2_ax2.plot(t, OCO[key], c=c['OCO'], label='OCO')
#f2_ax2.set_xlabel('days')
f2_ax2.set_ylim(1.5e21, 1.85e21)
if change_xlim: f2_ax2.set_xlim(0, x_max)

f2_ax2.set_ylabel(r'$\frac{{\mathrm{molecules}}}{{\mathrm{cm}^2}}$')

if PRINT_ALL == True:
    # plot retrieved pressure 
    key = 'p'
    f2_ax3.plot(t, hit16[key], c=c['hit16'], label='Hitran-16')
    f2_ax3.plot(t, hit20[key], c=c['hit20'], label='Hitran-20')
    f2_ax3.plot(t, tccon[key], c=c['tccon'], label='TCCON')
    f2_ax3.plot(t, OCO[key], c=c['OCO'], label='OCO')
    f2_ax3.plot(t, obs[key], 'o', markersize = 1, c=c['obs'], label='observation')


    #f2_ax3.set_xlabel('days')
    f2_ax3.set_ylabel('Pressure\n(HPa)')
    f2_ax3.set_ylim(805,855)
    if change_xlim: f2_ax3.set_xlim(0, x_max)

    # plot retrieved temperature 
    key = 'T'
    f2_ax4.plot(t, hit16[key], c=c['hit16'], label='Hitran-16')
    f2_ax4.plot(t, hit20[key], c=c['hit20'], label='Hitran-20')
    f2_ax4.plot(t, tccon[key], c=c['tccon'], label='TCCON')
    f2_ax4.plot(t, OCO[key], c=c['OCO'], label='OCO')
    f2_ax4.plot(t, obs[key], 'o', markersize = 1, c=c['obs'], label='observation')

    f2_ax4.set_xlabel('Days')
    f2_ax4.set_ylabel('Temperature\n(Kelvin)')
    f2_ax4.set_ylim(275,305)
    if change_xlim: f2_ax4.set_xlim(0, x_max)
else:
    f2_ax2.set_xlabel('Days')

fig2.tight_layout()
fig2.savefig(datadir+'figs/presentation/retrieved_DCS_CO2_pres.pdf')
print('\tSAVED: figs/presentation/retrieved_DCS_CO2_pres.pdf')

#################################### H20 ########################
if PRINT_ALL == True:
    fig3, (f3_ax1, f3_ax2, f3_ax3, f3_ax4) = plt.subplots(nrows=4, sharex=True,figsize=[6, 9], dpi=300)            
else:
    fig3, (f3_ax1, f3_ax2) = plt.subplots(nrows=2, sharex=True,figsize=[6, 4], dpi=300)            

key = 'H2O_VMR'
f3_ax1.plot(t, hit16[key], c=c['hit16'], label='Hitran-16')
f3_ax1.plot(t, hit20[key], c=c['hit20'], label='Hitran-20')
f3_ax1.plot(t, tccon[key], c=c['tccon'], label='TCCON')
f3_ax1.plot(t, OCO[key], c=c['OCO'], label='OCO')
f3_ax1.plot(obs['t_pic'], obs[key], 'o', markersize = 1, c=c['obs'], label='observation')

#f3_ax1.set_xlabel('days')
f3_ax1.set_ylabel('Percent')
#f3_ax1.legend(ncol = 1, frameon=False, loc='upper right')
if change_xlim: f3_ax1.set_xlim(0, x_max)
f3_ax1.set_ylim(0,3)
f3_ax1.set_title('Retrieved $H_2O$')

# plot the column
key = 'H2O'
f3_ax2.plot(t, hit16[key], c=c['hit16'], label='Hitran-16')
f3_ax2.plot(t, hit20[key], c=c['hit20'], label='Hitran-20')
f3_ax2.plot(t, tccon[key], c=c['tccon'], label='TCCON')
f3_ax2.plot(t, OCO[key], c=c['OCO'], label='OCO')
#f3_ax2.set_xlabel('days')
f3_ax2.set_ylim(1e22, 7e22)
f3_ax2.set_ylabel(r'$\frac{{\mathrm{molecules}}}{{\mathrm{cm}^2}}$')
if change_xlim: f3_ax2.set_xlim(0, x_max)

if PRINT_ALL == True:
    # plot retrieved pressure 
    key = 'p'
    f3_ax3.plot(t, hit16[key], c=c['hit16'], label='Hitran-16')
    f3_ax3.plot(t, hit20[key], c=c['hit20'], label='Hitran-20')
    f3_ax3.plot(t, tccon[key], c=c['tccon'], label='TCCON')
    f3_ax3.plot(t, OCO[key], c=c['OCO'], label='OCO')
    f3_ax3.plot(t, obs[key], 'o', markersize = 1,  c=c['obs'], label='obs')
    #f3_ax3.set_xlabel('days')
    f3_ax3.set_ylabel('Pressure\n(HPa)')
    f3_ax3.set_ylim(805,850)
    if change_xlim: f3_ax3.set_xlim(0, x_max)

    # plot retrieved temperature 
    key = 'T'
    f3_ax4.plot(t, hit16[key], c=c['hit16'], label='Hitran-16')
    f3_ax4.plot(t, hit20[key], c=c['hit20'], label='Hitran-20')
    f3_ax4.plot(t, tccon[key], c=c['tccon'], label='TCCON')
    f3_ax4.plot(t, OCO[key], c=c['OCO'], label='OCO')
    f3_ax4.plot(t, obs[key], 'o', markersize = 1, c=c['obs'], label='observation')
    f3_ax4.set_ylim(275,305)
    if change_xlim: f3_ax4.set_xlim(0, x_max)

    f3_ax4.set_xlabel('days')
    f3_ax4.set_ylabel('Temperature\n(K)')
else:
    f3_ax2.set_xlabel('days')
    
fig3.tight_layout()
fig3.savefig(datadir+'figs/presentation/retrieved_DCS_H2O_pres.pdf')
print('\tSAVED: figs/presentation/retrieved_DCS_H2O_pres.pdf')

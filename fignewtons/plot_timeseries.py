from setup_retrieval_plots import *
from process_results import *
             
print('Running em_plot_timeseries.py...')

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
obs = read_obs(datadir+'timeseries/data/obs.jld2') # piccaro tower pT data 
t = obs['t_pT'] 

print('\n\nOBS FILE WAS NOT UPDATED AS OF THE 12TH!\n')

# replace t, because time-stamps are inconsistant 
t = np.arange(start=0, stop=t[-1], step=(t[-1] - t[0])/len(t))

## Fig. 1: CH4 time-series and retrieval 
fig1, (f1_ax1, f1_ax2, f1_ax3, f1_ax4) = plt.subplots(nrows=4, sharex=True,figsize=[6, 9], dpi=300)            

# CH4 concentrations
key = 'CH4_VMR'
f1_ax1.plot(t, hit08[key], c=c['hit08'], label='Hitran-08')
f1_ax1.plot(t, hit16[key], c=c['hit16'], label='Hitran-16')
f1_ax1.plot(t, hit20[key], c=c['hit20'], label='Hitran-20')
f1_ax1.plot(t, tccon[key], c=c['tccon'], label='TCCON')
f1_ax1.plot(obs['t_pic'], obs[key], 'o', markersize = 1, c=c['obs'], label='observation')
#f1_ax1.set_xlabel('days')
f1_ax1.set_ylabel('ppb')
f1_ax1.legend(loc='upper right')

# plot the column
key = 'CH4'
f1_ax2.plot(t, hit08[key], c=c['hit08'], label='Hitran-08')
f1_ax2.plot(t, hit16[key], c=c['hit16'], label='Hitran-16')
f1_ax2.plot(t, hit20[key], c=c['hit20'], label='Hitran-20')
f1_ax2.plot(t, tccon[key], c=c['tccon'], label='TCCON')
#f1_ax2.set_xlabel('days')
#f1_ax2.set_ylim(0, 6e22)
f1_ax2.set_ylabel(r'$\frac{{\mathrm{molecules}}}{{\mathrm{cm}^2}}$')

# plot retrieved pressure 
key = 'p_ch4'
f1_ax3.plot(t, hit08[key], c=c['hit08'], label='Hitran-08')
f1_ax3.plot(t, hit16[key], c=c['hit16'], label='Hitran-16')
f1_ax3.plot(t, hit20[key], c=c['hit20'], label='Hitran-20')
f1_ax3.plot(t, tccon[key], c=c['tccon'], label='TCCON')
f1_ax3.plot(t, obs['p'], 'o', markersize = 1, c=c['obs'], label='observation')
#f1_ax3.set_xlabel('days')
f1_ax3.set_ylabel('HPa')
f1_ax3.set_ylim(800, 950)

# plot retrieved temperature 
key = 'T_ch4'
f1_ax4.plot(t, hit08[key], c=c['hit08'], label='Hitran-08')
f1_ax4.plot(t, hit16[key], c=c['hit16'], label='Hitran-16')
f1_ax4.plot(t, hit20[key], c=c['hit20'], label='Hitran-20')
f1_ax4.plot(t, tccon[key], c=c['tccon'], label='TCCON')
f1_ax4.plot(t, obs['T'], 'o', markersize = 1, c=c['obs'], label='observation')

f1_ax4.set_xlabel('days')
f1_ax4.set_ylabel('Kelvin')
f1_ax4.set_ylim(260, 310)
fig1.tight_layout()
fig1.savefig(datadir+'figs/retrieved_DCS_CH4.pdf')
print('\tSAVED: figs/retrieved_DCS_CH4.pdf')


## Fig. 1: CH4 time-series and retrieval 
fig2, (f2_ax1, f2_ax2, f2_ax3, f2_ax4) = plt.subplots(nrows=4, sharex=True,figsize=[6, 9], dpi=300)            

# CH4 concentrations
key = 'CO2_VMR'
f2_ax1.plot(t, hit16[key], c=c['hit16'], label='Hitran-16')
f2_ax1.plot(t, hit20[key], c=c['hit20'], label='Hitran-20')
f2_ax1.plot(t, tccon[key], c=c['tccon'], label='TCCON')
f2_ax1.plot(t, OCO[key], c=c['OCO'], label='OCO')
f2_ax1.plot(obs['t_pic'], obs[key], 'o', markersize = 1, c=c['obs'], label='obs')

#f2_ax1.set_xlabel('days')
f2_ax1.set_ylabel('Concentration (ppm)')
f2_ax1.legend(loc='upper right')
f2_ax1.set_ylim(360, 480)

# plot the column
key = 'CO2'
f2_ax2.plot(t, hit16[key], c=c['hit16'], label='Hitran-16')
f2_ax2.plot(t, hit20[key], c=c['hit20'], label='Hitran-20')
f2_ax2.plot(t, tccon[key], c=c['tccon'], label='TCCON')
f2_ax2.plot(t, OCO[key], c=c['OCO'], label='OCO')
#f2_ax2.set_xlabel('days')
#f2_ax2.set_ylim(1.5e21, 2e21)
f2_ax2.set_ylabel(r'$\frac{{\mathrm{molecules}}}{{\mathrm{cm}^2}}$')
# plot retrieved pressure 
key = 'p'
f2_ax3.plot(t, hit16[key], c=c['hit16'], label='Hitran-16')
f2_ax3.plot(t, hit20[key], c=c['hit20'], label='Hitran-20')
f2_ax3.plot(t, tccon[key], c=c['tccon'], label='TCCON')
f2_ax3.plot(t, OCO[key], c=c['OCO'], label='OCO')
f2_ax3.plot(t, obs[key], 'o', markersize = 1, c=c['obs'], label='observation')

#f2_ax3.set_xlabel('days')
f2_ax3.set_ylabel('Pressure (HPa)')
f2_ax3.set_ylim(800,900)

# plot retrieved temperature 
key = 'T'
f2_ax4.plot(t, hit16[key], c=c['hit16'], label='Hitran-16')
f2_ax4.plot(t, hit20[key], c=c['hit20'], label='Hitran-20')
f2_ax4.plot(t, tccon[key], c=c['tccon'], label='TCCON')
f2_ax4.plot(t, OCO[key], c=c['OCO'], label='OCO')
f2_ax4.plot(t, obs[key], 'o', markersize = 1, c=c['obs'], label='observation')

f2_ax4.set_xlabel('Days')
f2_ax4.set_ylabel('Temperature (Kelvin)')
f2_ax4.set_ylim(260,310)

fig2.tight_layout()
fig2.savefig(datadir+'figs/retrieved_DCS_CO2.pdf')
print('\tSAVED: figs/retrieved_DCS_CO2.pdf')

### Plot H2O
fig3, (f3_ax1, f3_ax2, f3_ax3, f3_ax4) = plt.subplots(nrows=4, sharex=True,figsize=[6, 9], dpi=300)            

# CH4 concentrations
key = 'H2O_VMR'
f3_ax1.plot(t, hit16[key], c=c['hit16'], label='Hitran-16')
f3_ax1.plot(t, hit20[key], c=c['hit20'], label='Hitran-20')
f3_ax1.plot(t, tccon[key], c=c['tccon'], label='TCCON')
f3_ax1.plot(t, OCO[key], c=c['OCO'], label='OCO')
f3_ax1.plot(obs['t_pic'], obs[key], 'o', markersize = 1, c=c['obs'], label='observation')

#f3_ax1.set_xlabel('days')
f3_ax1.set_ylabel('Percent')
f3_ax1.legend(loc='upper right')

# plot the column
key = 'H2O'
f3_ax2.plot(t, hit16[key], c=c['hit16'], label='Hitran-16')
f3_ax2.plot(t, hit20[key], c=c['hit20'], label='Hitran-20')
f3_ax2.plot(t, tccon[key], c=c['tccon'], label='TCCON')
f3_ax2.plot(t, OCO[key], c=c['OCO'], label='OCO')
#f3_ax2.set_xlabel('days')
f3_ax2.set_ylim(0, 7e22)
f3_ax2.set_ylabel(r'$\frac{{\mathrm{molecules}}}{{\mathrm{cm}^2}}$')


# plot retrieved pressure 
key = 'p'
f3_ax3.plot(t, hit16[key], c=c['hit16'], label='Hitran-16')
f3_ax3.plot(t, hit20[key], c=c['hit20'], label='Hitran-20')
f3_ax3.plot(t, tccon[key], c=c['tccon'], label='TCCON')
f3_ax3.plot(t, OCO[key], c=c['OCO'], label='OCO')
f3_ax3.plot(t, obs[key], 'o', markersize = 1,  c=c['obs'], label='obs')
#f3_ax3.set_xlabel('days')
f3_ax3.set_ylabel('Pressure (HPa)')
f3_ax3.set_ylim(800,900)
# plot retrieved temperature 
key = 'T'
f3_ax4.plot(t, hit16[key], c=c['hit16'], label='Hitran-16')
f3_ax4.plot(t, hit20[key], c=c['hit20'], label='Hitran-20')
f3_ax4.plot(t, tccon[key], c=c['tccon'], label='TCCON')
f3_ax4.plot(t, OCO[key], c=c['OCO'], label='OCO')
f3_ax4.plot(t, obs[key], 'o', markersize = 1, c=c['obs'], label='observation')
f3_ax4.set_ylim(225,325)

f3_ax4.set_xlabel('days')
f3_ax4.set_ylabel('Temperature (Kelvin)')

fig3.tight_layout()
#plt.subplots_adjust(hspace=0.05)
fig3.savefig(datadir+'figs/retrieved_DCS_H2O.pdf')
print('\tSAVED: figs/retrieved_DCS_H2O.pdf')

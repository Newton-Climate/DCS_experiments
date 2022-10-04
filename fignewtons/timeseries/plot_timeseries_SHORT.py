from setup_retrieval_plots import *
from process_results import *
             
print('Running em_plot_timeseries.py...')

# How many days to plot??
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
#t = np.arange(start=0, stop=t[-1], step=(t[-1] - t[0])/len(t))

#################################### CH4 ########################
fig1, (f1_ax1, f1_ax2, f1_ax3, f1_ax4) = plt.subplots(nrows=4, sharex=True,figsize=[6, 9], dpi=300)            

# CH4 concentrations
key = 'CH4_VMR'
f1_ax1.plot(t, hit08[key], 'o',markersize = 1, c=c['hit08'], label='Hitran-08')
f1_ax1.plot(t, hit16[key], 'o',markersize = 1,c=c['hit16'], label='Hitran-16')
f1_ax1.plot(t, hit20[key], 'o',markersize = 1,c=c['hit20'], label='Hitran-20')
f1_ax1.plot(t, tccon[key], 'o',markersize = 1,c=c['tccon'], label='TCCON')
f1_ax1.plot(obs['t_pic'], obs[key], 'o', markersize = 1, c=c['obs'], label='observation')
#f1_ax1.set_xlabel('days')
f1_ax1.set_ylabel('ppb')
f1_ax1.legend(ncol = 2, frameon=False, loc='upper right')
#if change_xlim: f1_ax1.set_xlim(0, x_max)
f1_ax1.set_ylim(1850, 2200)
f1_ax1.set_title('Retrieved CH$_4$')

# plot the column
key = 'CH4'
f1_ax2.plot(t, hit08[key], 'o',markersize = 2,c=c['hit08'], label='Hitran-08')
f1_ax2.plot(t, hit16[key], 'o',markersize = 2,c=c['hit16'], label='Hitran-16')
f1_ax2.plot(t, hit20[key], 'o',markersize = 2,c=c['hit20'], label='Hitran-20')
f1_ax2.plot(t, tccon[key], 'o',markersize = 2,c=c['tccon'], label='TCCON')
#f1_ax2.set_xlabel('days')
f1_ax2.set_ylim(0.7e19, 0.9e19)
f1_ax2.set_ylabel(r'$\frac{{\mathrm{molecules}}}{{\mathrm{cm}^2}}$')
if change_xlim: f1_ax2.set_xlim(0, x_max)

# plot retrieved pressure 
key = 'p_ch4'
f1_ax3.plot(t, hit08[key], 'o',markersize = 2,c=c['hit08'], label='Hitran-08')
f1_ax3.plot(t, hit16[key], 'o',markersize = 2,c=c['hit16'], label='Hitran-16')
f1_ax3.plot(t, hit20[key], 'o',markersize = 2,c=c['hit20'], label='Hitran-20')
f1_ax3.plot(t, tccon[key], 'o',markersize = 2,c=c['tccon'], label='TCCON')
f1_ax3.plot(t, obs['p'], 'x', markersize = 1, c=c['obs'], label='observation')
#f1_ax3.set_xlabel('days')
f1_ax3.set_ylabel('hPa')
f1_ax3.set_ylim(780, 890)
if change_xlim: f1_ax3.set_xlim(0, x_max)

# plot retrieved temperature 
key = 'T_ch4'
f1_ax4.plot(t, hit08[key], 'o',markersize = 2, c=c['hit08'], label='Hitran-08')
f1_ax4.plot(t, hit16[key], 'o',markersize = 2,c=c['hit16'], label='Hitran-16')
f1_ax4.plot(t, hit20[key], 'o',markersize = 2,c=c['hit20'], label='Hitran-20')
f1_ax4.plot(t, tccon[key], 'o',markersize = 2,c=c['tccon'], label='TCCON')
f1_ax4.plot(t, obs['T'], 'x', markersize = 2, c=c['obs'], label='observation')

f1_ax4.set_xlabel('days')
f1_ax4.set_ylabel('Kelvin')
f1_ax4.set_ylim(270,305)
if change_xlim: f1_ax4.set_xlim(0, x_max)

fig1.tight_layout()
fig1.savefig(datadir+'/figs/retrieved_DCS_CH4_short3.pdf')
print('\tSAVED: figs/retrieved_DCS_CH4_short.pdf')


#################################### C02 ########################
fig2, (f2_ax1, f2_ax2, f2_ax3, f2_ax4) = plt.subplots(nrows=4, sharex=True,figsize=[6, 9], dpi=300)            

key = 'CO2_VMR'
f2_ax1.plot(t, hit16[key],'o',markersize = 2, c=c['hit16'], label='Hitran-16')
f2_ax1.plot(t, hit20[key],'o',markersize = 2, c=c['hit20'], label='Hitran-20')
f2_ax1.plot(t, tccon[key],'o',markersize = 2, c=c['tccon'], label='TCCON')
f2_ax1.plot(t, OCO[key],'o',markersize = 2, c=c['OCO'], label='OCO')
f2_ax1.plot(obs['t_pic'], obs[key], 'o', markersize = 1, c=c['obs'], label='observation')

#f2_ax1.set_xlabel('days')
f2_ax1.set_ylabel('Concentration (ppm)')
f2_ax1.legend(ncol = 2, frameon=False, loc='upper left')
f2_ax1.set_ylim(385, 475)
if change_xlim: f2_ax1.set_xlim(0, x_max)
f2_ax1.set_title('Retrieved $CO_2$')


# plot the column
key = 'CO2'
f2_ax2.plot(t, hit16[key],'o',markersize = 2, c=c['hit16'], label='Hitran-16')
f2_ax2.plot(t, hit20[key],'o',markersize = 2, c=c['hit20'], label='Hitran-20')
f2_ax2.plot(t, tccon[key],'o',markersize = 2, c=c['tccon'], label='TCCON')
f2_ax2.plot(t, OCO[key],'o',markersize = 2, c=c['OCO'], label='OCO')
#f2_ax2.set_xlabel('days')
f2_ax2.set_ylim(1.5e21, 1.85e21)
if change_xlim: f2_ax2.set_xlim(0, x_max)

f2_ax2.set_ylabel(r'$\frac{{\mathrm{molecules}}}{{\mathrm{cm}^2}}$')
# plot retrieved pressure 
key = 'p'
f2_ax3.plot(t, hit16[key],'o',markersize = 2, c=c['hit16'], label='Hitran-16')
f2_ax3.plot(t, hit20[key],'o',markersize = 2, c=c['hit20'], label='Hitran-20')
f2_ax3.plot(t, tccon[key],'o',markersize = 2, c=c['tccon'], label='TCCON')
f2_ax3.plot(t, OCO[key],'o',markersize = 2, c=c['OCO'], label='OCO')
f2_ax3.plot(t, obs[key], 'o', markersize = 1, c=c['obs'], label='observation')

#f2_ax3.set_xlabel('days')
f2_ax3.set_ylabel('Pressure (hPa)')
f2_ax3.set_ylim(805,855)
if change_xlim: f2_ax3.set_xlim(0, x_max)

# plot retrieved temperature 
key = 'T'
f2_ax4.plot(t, hit16[key],'o',markersize = 2, c=c['hit16'], label='Hitran-16')
f2_ax4.plot(t, hit20[key],'o',markersize = 2, c=c['hit20'], label='Hitran-20')
f2_ax4.plot(t, tccon[key],'o',markersize = 2, c=c['tccon'], label='TCCON')
f2_ax4.plot(t, OCO[key],'o',markersize = 2, c=c['OCO'], label='OCO')
f2_ax4.plot(t, obs[key], 'o', markersize = 1, c=c['obs'], label='observation')

f2_ax4.set_xlabel('Days')
f2_ax4.set_ylabel('Temperature (Kelvin)')
f2_ax4.set_ylim(275,305)
if change_xlim: f2_ax4.set_xlim(0, x_max)

fig2.tight_layout()
fig2.savefig(datadir+'/figs/retrieved_DCS_CO2_short.pdf')
print('\tSAVED: figs/retrieved_DCS_CO2_short.pdf')

#################################### H20 ########################
fig3, (f3_ax1, f3_ax2, f3_ax3, f3_ax4) = plt.subplots(nrows=4, sharex=True,figsize=[6, 9], dpi=300)            

key = 'H2O_VMR'
f3_ax1.plot(t, hit16[key],'o',markersize = 2, c=c['hit16'], label='Hitran-16')
f3_ax1.plot(t, hit20[key],'o',markersize = 2, c=c['hit20'], label='Hitran-20')
f3_ax1.plot(t, tccon[key],'o',markersize = 2, c=c['tccon'], label='TCCON')
f3_ax1.plot(t, OCO[key],'o',markersize = 2, c=c['OCO'], label='OCO')
f3_ax1.plot(obs['t_pic'], obs[key], 'o', markersize = 1, c=c['obs'], label='observation')

#f3_ax1.set_xlabel('days')
f3_ax1.set_ylabel('Percent')
f3_ax1.legend(ncol = 1, frameon=False, loc='upper right')
if change_xlim: f3_ax1.set_xlim(0, x_max)
f3_ax1.set_ylim(0,3)
f3_ax1.set_title('Retrieved H$_2$O')

# plot the column
key = 'H2O'
f3_ax2.plot(t, hit16[key],'o',markersize = 2, c=c['hit16'], label='Hitran-16')
f3_ax2.plot(t, hit20[key],'o',markersize = 2, c=c['hit20'], label='Hitran-20')
f3_ax2.plot(t, tccon[key],'o',markersize = 2, c=c['tccon'], label='TCCON')
f3_ax2.plot(t, OCO[key],'o',markersize = 2, c=c['OCO'], label='OCO')
#f3_ax2.set_xlabel('days')
f3_ax2.set_ylim(1e22, 7e22)
f3_ax2.set_ylabel(r'$\frac{{\mathrm{molecules}}}{{\mathrm{cm}^2}}$')
if change_xlim: f3_ax2.set_xlim(0, x_max)


# plot retrieved pressure 
key = 'p'
f3_ax3.plot(t, hit16[key],'o',markersize = 2, c=c['hit16'], label='Hitran-16')
f3_ax3.plot(t, hit20[key],'o',markersize = 2, c=c['hit20'], label='Hitran-20')
f3_ax3.plot(t, tccon[key],'o',markersize = 2, c=c['tccon'], label='TCCON')
f3_ax3.plot(t, OCO[key],'o',markersize = 2, c=c['OCO'], label='OCO')
f3_ax3.plot(t, obs[key], 'o', markersize = 1,  c=c['obs'], label='obs')
#f3_ax3.set_xlabel('days')
f3_ax3.set_ylabel('Pressure (hPa)')
f3_ax3.set_ylim(805,850)
if change_xlim: f3_ax3.set_xlim(0, x_max)

# plot retrieved temperature 
key = 'T'
f3_ax4.plot(t, hit16[key],'o',markersize = 2, c=c['hit16'], label='Hitran-16')
f3_ax4.plot(t, hit20[key],'o',markersize = 2, c=c['hit20'], label='Hitran-20')
f3_ax4.plot(t, tccon[key],'o',markersize = 2, c=c['tccon'], label='TCCON')
f3_ax4.plot(t, OCO[key],'o',markersize = 2, c=c['OCO'], label='OCO')
f3_ax4.plot(t, obs[key], 'o', markersize = 1, c=c['obs'], label='observation')
f3_ax4.set_ylim(275,305)
if change_xlim: f3_ax4.set_xlim(0, x_max)

f3_ax4.set_xlabel('days')
f3_ax4.set_ylabel('Temperature (Kelvin)')

fig3.tight_layout()
#plt.subplots_adjust(hspace=0.05)
fig3.savefig(datadir+'figs/retrieved_DCS_H2O_short.pdf')
print('\tSAVED: figs/retrieved_DCS_H2O_short.pdf')

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import h5py
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as colors

print('Running plot_heatmap_ch4.py...')

# For resetting x and y ticks:
p = np.arange(start=400, step=20, stop=1001)
T = np.arange(start=200, step=3, stop=300)

t_num = 6
t_tick_labels = np.linspace(200, 300, num = t_num)
t_ticks = np.linspace(0, len(T), num = t_num)
t_tick_labels = [int(i) for i in t_tick_labels]

p_num = 5
p_tick_labels = np.linspace(400, 1000, num = p_num)
p_ticks = np.linspace(0, len(p)-1, num = p_num)
p_tick_labels = [int(i) for i in p_tick_labels]

def calc_errors(data_dict, p_true, T_true):
    p_retrieved = data_dict['p']
    T_retrieved = data_dict['T']
    (nT, nP) = p_retrieved.shape
    dp = np.empty((nT, nP))
    dT = dp.copy()

    # Negative values = True is bigger than retrieved
    # Positive values = True is smaller than retrieved
    for i in range(nP):
        dT[:,i] = T_retrieved[:,i] - T_true
    # end of loop

    for i in range(nT):
        dp[i,:] = p_retrieved[i,:] - p_true
    # End of for-loop

    data_dict['dp'] = dp
    data_dict['dT'] = dT
    data_dict['dch4'] = data_dict['ch4'] - 2000
    return data_dict
# end of function 

def load_data(datafile):
    with h5py.File(datafile, 'r') as f:
        read = lambda x: f[x][:,:]
        read_vec = lambda x: f[x][:]

        out = {
            'ch4_col' : read('ch4_col'),
            'ch4' : read('ch4'),
            'h2o_col' : read('h2o_col'),
            'h2o' : read('h2o'),
            'vcd' : read('vcd'),
            'p' : read('p'),
#            'p_true' : read_vec('p_true'),
            'T' : read('T')}

        out = calc_errors(out, p, T)

    return out

hit16 = load_data('heatmaps/ch4_hit16_results.jld2')
hit20 = load_data('heatmaps/ch4_hit20_results.jld2')
TCCON = load_data('heatmaps/ch4_TCCON_results.jld2')

var = ['dch4', 'dch4', 'dch4', 'dp', 'dp', 'dp','dT', 'dT', 'dT']
vmax = [250, 250, 250, 15, 15, 15, 0.5, 0.5, 0.5]
vmin = [-25, -25, -25, -75, -75, -75, -1, -1, -1]
array_to_plot = [hit16, hit20, TCCON] * 3
cmap = 'PiYG'


"""
print('\ndch4')
print(np.nanmin(hit16['dch4']))
print(np.nanmin(hit20['dch4']))
print(np.nanmin(TCCON['dch4']))

print(np.nanmax(hit16['dch4']))
print(np.nanmax(hit20['dch4']))
print(np.nanmax(TCCON['dch4']))
print('\ndT')
print(np.nanmin(hit16['dT']))
print(np.nanmin(hit20['dT']))
print(np.nanmin(TCCON['dT']))

print(np.nanmax(hit16['dT']))
print(np.nanmax(hit20['dT']))
print(np.nanmax(TCCON['dT']))
print('\ndp')
print(np.nanmin(hit16['dp']))
print(np.nanmin(hit20['dp']))
print(np.nanmin(TCCON['dp']))

print(np.nanmax(hit16['dp']))
print(np.nanmax(hit20['dp']))
print(np.nanmax(TCCON['dp']))
"""

fig, ax = plt.subplots(nrows=3, ncols=3, sharex = True, sharey = True, figsize = (12,10), dpi = 300)
ax = ax.flatten() 

for i in range(9):
    divnorm = colors.TwoSlopeNorm(vmin=vmin[i], vcenter=0, vmax=vmax[i])
    pcm = ax[i].pcolormesh(array_to_plot[i][var[i]], rasterized=True, norm=divnorm, cmap=cmap,)
    
    if i == 0:
        divnorm = colors.TwoSlopeNorm(vmin=vmin[i], vcenter=0, vmax=vmax[i])
        im1 = ax[0].pcolormesh(array_to_plot[0][var[0]], rasterized=True, norm=divnorm, cmap=cmap,)
        ax[i].set_ylabel('Temperature (K)')
        ax[i].text(-15, 17, '$CH_4$\n(ppb)', va='center', ha ='center', fontsize = 14)
    if i == 3:
        divnorm = colors.TwoSlopeNorm(vmin=vmin[i], vcenter=0, vmax=vmax[i])
        im2 = ax[3].pcolormesh(array_to_plot[3][var[3]], rasterized=True, norm=divnorm, cmap=cmap,)
        ax[i].set_ylabel('Temperature (K)')
        ax[i].text(-15, 17, 'Pressure\n(HPa)', va='center', ha="center", fontsize = 14)

    if i == 6:
        divnorm = colors.TwoSlopeNorm(vmin=vmin[i], vcenter=0, vmax=vmax[i])
        im3 = ax[6].pcolormesh(array_to_plot[6][var[6]], rasterized=True, norm=divnorm, cmap=cmap,)
        ax[i].set_ylabel('Temperature (K)')
        ax[i].text(-15, 17, 'Temperature\n(K)', va='center', ha="center", fontsize = 14)
    divider = make_axes_locatable(ax[i])
    cax = divider.append_axes("right", size = "7%", pad = 0.05)
    cax.remove()
    if i == 2:
            cax1 = divider.append_axes("right", size="7%", pad=.05)
            fig.colorbar(pcm, cax=cax1, orientation='vertical', label='Retrieval Error (ppb)', ticks = [-25, -20, -15, -10, -5, 0, 50, 100, 150, 200, 250])
  
    if i == 5:
            cax2 = divider.append_axes("right", size="7%", pad=.05)
            fig.colorbar(pcm, cax=cax2, orientation='vertical',  label='Retrieval Error (HPa)', ticks = [-75, -50, -25, 0, 5, 10, 15])
    if i == 8:
            cax3 = divider.append_axes("right", size="7%", pad=.05)
            fig.colorbar(pcm, cax=cax3, orientation='vertical',  label='Retrieval Error (K)', ticks = [-1, -0.5, 0, 0.25, 0.5])
    ax[0].set_title('Hitran-2016', fontsize = 14)
    ax[1].set_title('Hitran-2020', fontsize = 14)
    ax[2].set_title('TCCON', fontsize = 14)

    if i > 5: ax[i].set_xlabel('Pressure (HPa)')
    ax[i].grid(False)  
    ax[i].set_xticks(ticks = p_ticks, labels = p_tick_labels)
    ax[i].set_yticks(ticks = t_ticks, labels = t_tick_labels)

    

plt.subplots_adjust(wspace=0.005, hspace=0)
plt.tight_layout()
fig.savefig('figs/heatmap_ch4.pdf')
print('\tSAVED: figs/heatmap_ch4.pdf')




import h5py 
import matplotlib.pyplot as plt
import numpy as np
# nicer figures using ggg plot style.
#plt.style.use('ggplot')
from IPython.core.display import HTML

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=False)
from IPython.core.pylabtools import figsize

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import numpy as np
from matplotlib.figure import Figure
import matplotlib.pylab as plt

import scipy.io as spio

# import a color pallet with better colors 
import seaborn as sns
cp = sns.color_palette(palette="deep")

c = {
    'obs' : 'black',
    'hit08' : cp[0],
    'hit16' : cp[1],
    'hit20' : cp[2],
    'tccon' : cp[3],
    'OCO' : cp[4]}

import sys

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from cartopy import config
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime,timedelta
from datetime import date as date_datetime
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm
from os import listdir
from os.path import isfile, join
import matplotlib.dates as mdates
from matplotlib.colors import LinearSegmentedColormap
import mpl_scatter_density


# "Viridis-like" colormap with white background
white_viridis = LinearSegmentedColormap.from_list('white_viridis', [
    (0, '#ffffff'),
    (1e-20, '#440053'),
    (0.2, '#404388'),
    (0.4, '#2a788e'),
    (0.6, '#21a784'),
    (0.8, '#78d151'),
    (1, '#fde624'),
], N=256)

#===Load user libraries====
sys.path.append('./libraries/')
from datamanager import DataManager_GRIDDED
from matplotlib import cm


viridis = cm.get_cmap('viridis', 256)
RdYlGn=cm.get_cmap('Spectral', 256)
#===Create Color bar====
def make_Ramp( ramp_colors ): 
    from colour import Color
    from matplotlib.colors import LinearSegmentedColormap

    color_ramp = LinearSegmentedColormap.from_list( 'my_list', [ Color( c1 ).rgb for c1 in ramp_colors ] )
    return color_ramp


custom_ramp = make_Ramp( ["#ffffff","#aec8e0","#9dafb2","#e2b85f","#f0ae46","#ea9842","#d36f2b","#a04617","#6c423d"] ) 

emep_pm10=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_emep-pm10_2008.nc')

emep_pm10.hist2d(biascorr=[1e9,1],cmap=white_viridis,xlim=[0.01,1.05],ylim=[0.01,1.05],log=True,save=True,title='PM$_{10}$',xlabel='EMEP')

emep_pm10.graph_map(variable='yr',biascorr=1,vmin=0,vmax=50,n_levels_plot=10,cmap=custom_ramp,type_graph='mean',date=False,save=False,title='EMEP PM$_{10}$ [ug/m3]')
emep_pm10.graph_map(variable='ys',biascorr=1e9,vmin=0,vmax=50,n_levels_plot=10,cmap=custom_ramp,type_graph='mean',date=False,save=False,title='LE PM$_{10}$ [ug/m3]')


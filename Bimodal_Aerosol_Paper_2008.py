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
import met_brewer
from colour import Color

#===Load user libraries====
sys.path.append('./libraries/')
from datamanager import DataManager_GRIDDED
from datamanager import DataManager_LE
from datamanager import graph_map
from matplotlib import cm


viridis = cm.get_cmap('viridis', 256)
redgreen = cm.get_cmap('RdYlGn', 256)
redblue = cm.get_cmap('RdYlBu', 256)
RdYlGn=cm.get_cmap('Spectral', 256)
bwr = cm.get_cmap('bwr', 256)
redgreen = redgreen.reversed()
redblue = redblue.reversed()
demuth = met_brewer.met_brew('Demuth',n=256,brew_type='continuous')
demuth = LinearSegmentedColormap.from_list( 'my_list', [ Color( c1 ).rgb for c1 in demuth ] )
white_viridis = LinearSegmentedColormap.from_list('white_viridis', [
    (0, '#ffffff'),
    (1e-20, '#440053'),
    (0.2, '#404388'),
    (0.4, '#2a788e'),
    (0.6, '#21a784'),
    (0.8, '#78d151'),
    (1, '#fde624'),
], N=256)
def make_Ramp( ramp_colors ): 
    from colour import Color
    from matplotlib.colors import LinearSegmentedColormap

    color_ramp = LinearSegmentedColormap.from_list( 'my_list', [ Color( c1 ).rgb for c1 in ramp_colors ] )
    return color_ramp
custom_ramp = make_Ramp( ["#ffffff","#aec8e0","#9dafb2","#e2b85f","#f0ae46","#ea9842","#d36f2b","#a04617","#6c423d"] ) 

def getmeanofneighbors(matrix, i, j,r):
    region = matrix[max(0, i-r) : i+r+1,
                    max(0, j-r) : j+r+1]
    #print(region)
    #return region
    return (np.sum(region) - matrix[i, j])/(region.size-1) # Sum the region and subtract center

save = True

Bimodal_polder_565 = DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_polder-aod_bimodal_Bimodal_2008.nc')
Bimodal_aeronet_500 = DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_aeronet-aod-500_bimodal_2008.nc')
Bimodal_aeronet_440 = DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_aeronet-aod-440_bimodal_2008.nc')
Bimodal_aeronet_870 = DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_aeronet-aod-870_bimodal_2008.nc')
Bimodal_polder_ae = DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_polder-ae_bimodal_2008.nc')


Bimodal_polder_565.graph_map(variable='ys',vmin=0,vmax=0.6,n_levels_plot=10,grid=True,save=save,save_title='Bimodal_2008_LE_POLDER_AOD565',title='LE-AOD565nm 2008',cmap=demuth.reversed(),date=False,facecolor=(1,1,1),extend='max')
Bimodal_polder_565.graph_map(variable='yr',vmin=0,vmax=0.6,n_levels_plot=10,grid=True,save=save,save_title='Bimodal_2008_POLDER_AOD565',title='POLDER-AOD565nm 2008',cmap=demuth.reversed(),date=False,facecolor=(1,1,1),extend='max')

Bimodal_aeronet_500.graph_map(variable='ys',vmin=0,vmax=0.6,n_levels_plot=10,grid=True,save=save,save_title='Bimodal_2008_LE_Aeronet_AOD500',title='LE-AOD500nm 2008',cmap=demuth.reversed(),date=False,facecolor=(1,1,1),extend='max',stations=(True))
Bimodal_aeronet_500.graph_map(variable='yr',vmin=0,vmax=0.6,n_levels_plot=10,grid=True,save=save,save_title='Bimodal_2008_Aeronet_AOD500',title='Aeronet-AOD500nm 2008',cmap=demuth.reversed(),date=False,facecolor=(1,1,1),extend='max',stations=(True))
Bimodal_aeronet_500.graph_diff_map(vmin=-100,vmax=100,norm2=(True),grid=True,save=save,save_title='Bimodal_2008_Aeronet_AOD500_diff',title='Diff LE - Aeronet AOD500nm 2008',cmap=bwr,date=False,facecolor=(1,1,1),extend='max',stations=(True),porcentage_diff=(True))
Bimodal_aeronet_500.hist2d(cmap=white_viridis,log=False,xlabel='Aeronet-AOD500nm',ylabel='LE-AOD500nm',xlim=[0,2],ylim=[0,2],save=save,save_title='Scatter_Aeronet_aod_500')
Bimodal_aeronet_440.graph_diff_map(vmin=-100,vmax=100,norm2=(True),grid=True,save=save,save_title='Bimodal_2008_Aeronet_AOD440_diff',title='Diff LE - Aeronet AOD440nm 2008',cmap=bwr,date=False,facecolor=(1,1,1),extend='max',stations=(True),porcentage_diff=(True))
Bimodal_aeronet_440.hist2d(cmap=white_viridis,log=False,xlabel='Aeronet-AOD440nm',ylabel='LE-AOD440nm',xlim=[0,2],ylim=[0,2],save=save,save_title='Scatter_Aeronet_aod_440')
Bimodal_aeronet_870.graph_diff_map(vmin=-100,vmax=100,norm2=(True),grid=True,save=save,save_title='Bimodal_2008_Aeronet_AOD870_diff',title='Diff LE - Aeronet AOD870nm 2008',cmap=bwr,date=False,facecolor=(1,1,1),extend='max',stations=(True),porcentage_diff=(True))
Bimodal_aeronet_870.hist2d(cmap=white_viridis,log=False,xlabel='Aeronet-AOD870nm',ylabel='LE-AOD870nm',xlim=[0,2],ylim=[0,2],save=save,save_title='Scatter_Aeronet_aod_870')


Bimodal_polder_ae.graph_map(variable='ys',vmin=0,vmax=3,n_levels_plot=10,grid=True,save=save,save_title='Bimodal_2008_LE_POLDER_AE',title='LE-AE  443nm-865nm 2008',cmap=RdYlGn.reversed(),date=False,facecolor=(1,1,1),extend='max')
Bimodal_polder_ae.graph_map(variable='yr',vmin=0,vmax=3,n_levels_plot=10,grid=True,save=save,save_title='Bimodal_2008_POLDER_AE',title='POLDER-AE 443nm-865nm 2008',cmap=RdYlGn.reversed(),date=False,facecolor=(1,1,1),extend='max')




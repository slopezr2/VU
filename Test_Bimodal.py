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

Janot_radii_polder = DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_polder-aod_bimodal_Janot_radii_2008.nc')
Janot_radii_run = DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_Test_Janot_radii_aod.nc')
Bimodal_run = DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_Test_bimodal_aod.nc')





Janot_radii_polder.graph_map(variable='ys',vmin=0,vmax=0.6,n_levels_plot=10,grid=True,save=save,save_title='Janot_radii_polder_V1',title='LE-Janot_radii-AOD565nm 2008',cmap=demuth.reversed(),date=False,facecolor=(1,1,1),extend='max')

Janot_radii_polder.graph_map(variable='ys',vmin=0,vmax=0.6,n_levels_plot=10,grid=True,save=save,save_title='Janot_radii_polder_V1',title='LE-Janot_radii-AOD565nm 2008',cmap=demuth.reversed(),date=False,facecolor=(1,1,1),extend='max')
Janot_radii_run.graph_map(variable='angstrom_aeronet',vmin=0,vmax=2,n_levels_plot=10,grid=True,save=save,save_title='LE_Angstrom_Janot_radii_V1',title='LE-Janot_radii Angstrom 2008',cmap=RdYlGn.reversed(),date=False,facecolor=(1,1,1),extend='max')
Bimodal_run.graph_map(variable='angstrom_aeronet',vmin=0,vmax=2,n_levels_plot=10,grid=True,save=save,save_title='LE_Angstrom_Bimodal',title='LE-Bimodal Angstrom 2008',cmap=RdYlGn.reversed(),date=False,facecolor=(1,1,1),extend='max')




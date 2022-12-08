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

save = False

Natural_run = DataManager_LE('/Users/santiago/Documents/OSSE/OSSE_aod2_Natural_run_V3_200805.nc')
# Bimodal_run = DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_Test_bimodal_aod.nc')
# radii_01_run = DataManager_LE('/Users/santiago/Documents/LE_outputs/Different_radii/LE_m_Test_radii_01_aod.nc')


Natural_run.graph_map(variable='aod_563nm',vmin=0,vmax=0.7,n_levels_plot=10,grid=True,save=save,save_title='LE_Natural_Dust_Bounday_V1',
                          title='LE-Natural AOD 550nm 200805',cmap=demuth.reversed(),date=True
                          ,facecolor=(1,1,1),extend='max')

Dust_boundary_ens = {}
for i in range(1,11):
    ens = f'{i:02d}'
    print(ens)
    Dust_boundary_ens[ens] = DataManager_LE('/Users/santiago/Documents/LE_outputs/Smoother/LE_m_aod_Dust_Boundary_xi'+ens+'_V2.nc')
    Dust_boundary_ens[ens].graph_map(variable='aod_550nm',vmin=0,vmax=0.7,n_levels_plot=10,grid=True,save=save,save_title='LE_Dust_Bounday_Ens_'+ens+'_V2',
                          title='LE-Ens '+ens+' AOD 550nm 200805',cmap=demuth.reversed(),date=True
                          ,facecolor=(1,1,1),extend='max')















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

#==Save graphs==
save = True

#===Natural Run Version==
version = 'V2'
#==POLDER AOD and LE-AOD==
control_run = DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_polder-aod-563_2008.nc')
control_run_aeronet = DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_aeronet_500_2008.nc')
lat = control_run.nc.variables['latitude'][:]
lon = control_run.nc.variables['longitude'][:]
#==Map of relation between polder and LE==
control_run.graph_map(variable='yr',vmin=0,vmax=0.6,n_levels_plot=10,grid=True,save=save,save_title='OSSE_Polder_2008',title='POLDER-AOD565nm 2008',cmap=demuth.reversed(),date=False,facecolor=(1,1,1),extend='max')
control_run.graph_map(variable='ys',vmin=0,vmax=0.6,n_levels_plot=10,grid=True,save=save,save_title='OSSE_LE_CTRL_2008',title='LE-CTRL-AOD565nm 2008',cmap=demuth.reversed(),date=False,facecolor=(1,1,1),extend='max')

control_run.graph_diff_map(end_period=99,porcentage_diff=True,vmin=0,vmax=5,n_levels_plot=20,relation=True,norm2=True,vcenter=1,grid=True,save=save,save_title='OSSE_Relation_LE_Polder',extend='max')
control_run_aeronet.graph_diff_map(date=True,ini_period=1550,end_period=2040,porcentage_diff=True,vmin=0,vmax=5,n_levels_plot=20,relation=True,norm2=True,vcenter=1,grid=True,save=save,save_title='OSSE_Relation_LE_aeronet',extend='max',title='obs/LE ctrl aeronet-aod-500',stations=True,size_stations=25,facecolor=(1,1,1))


prueba = np.ones((lat.shape[0],lon.shape[0])) * 1.5
v1 = 10
v2 = 5
v3= 7
prueba[80:,:] = v1
prueba[:,0:20] = v1
prueba[8:40,10:35] = v2
prueba[30:65,20:47] = v2
prueba[0:8,15:55] = v2
prueba[40:50,45:60] = v2
prueba[55:80,65:] = v3
prueba[1,:] = 3000
for i in range(prueba.shape[0]):
    for j in range(prueba.shape[1]):
        prueba[i,j] = getmeanofneighbors(prueba, i, j,3)

graph_map(variable=prueba,lat=lat,lon=lon,vmax=1000,cmap=redblue,vmin=1,title='DC Correction Factors',title_on=True,stock_image=True,grid=True,n_levels_plot=10,save=save,save_title='OSSE_DC_Natural')


#===Experiments Natural run===
natural_ens = {}
if version=='V1':
    for i in range(1,11):
        ens = f'{i:02d}'
        print(ens)
        natural_ens[ens] = DataManager_GRIDDED('/Users/santiago/Documents/OSSE/OSSE_polder_Natural_run_2008_xi'+ens+'.nc')
        natural_ens[ens].graph_diff_map(end_period=-1,porcentage_diff=True,vmin=0,vmax=5,n_levels_plot=20,relation=True,norm2=True,vcenter=1,grid=True,date=True,save=save,save_title='OSSE_Relation_Natural_xi'+ens+'_Polder',extend='max',title='obs/LE xi'+ens)
elif version=='V2':
    for i in range(1,4):
        ens = f'{i:02d}'
        print(ens)
        natural_ens[ens] = DataManager_GRIDDED('/Users/santiago/Documents/OSSE/OSSE_polder_Natural_run_V2_2008_xi'+ens+'.nc')
        natural_ens[ens].graph_diff_map(end_period=-1,porcentage_diff=True,vmin=0,vmax=5,n_levels_plot=20,relation=True,norm2=True,vcenter=1,grid=True,date=True,save=save,save_title='OSSE_Relation_Natural_xi'+ens+'_Polder',extend='max',title='obs/LE xi'+ens)
        if i==3:
            natural_aeronet = DataManager_GRIDDED('/Users/santiago/Documents/OSSE/OSSE_aeronet-aod-500_Natural_run_V2_200805_xi'+ens+'.nc')
            natural_aeronet.graph_diff_map(end_period=-1,porcentage_diff=True,vmin=0,vmax=5,n_levels_plot=20,relation=True,norm2=True,vcenter=1,grid=True,date=True,save=save,save_title='OSSE_Relation_Natural_xi'+ens+'_Aeronet-aod-500',extend='max',title='obs/LE xi'+ens+' aeronet-aod-500',stations=True,size_stations=25,facecolor=(1,1,1))
           
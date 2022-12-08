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


mypath='/Users/santiago/Documents/CAMS_reanalysis/AOD_550nm_Daily/'
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
onlyfiles.sort()
i = 0
for file in onlyfiles:
    CAMS = Dataset(mypath+file)
    if i==0:
        aux = np.flip(np.mean(CAMS.variables['od550aer'][:,:,:],axis=0))
        lon = [f-180 for f in CAMS.variables['longitude'][:]]
        lat = np.flip(CAMS.variables['latitude'][:])
    else:
        aux += np.flip(np.mean(CAMS.variables['od550aer'][:,:,:],axis=0))
    i += 1
    
aod = aux/12
aod2 = np.zeros((181,360))
aod2[:,0:180] = aod[:,180:]
aod2[:,180:] = aod[:,0:180]

graph_map(variable=aod2,lon=lon,lat=lat,cmap=demuth.reversed(),n_levels_plot=10,subdomain=True,lat_lim=[100,165],lon_lim=[150,250],vmax=.6,grid=(True),title_on=(True),title='CAMS yearly average 2008',save=(True))

CAMS_daily = Dataset('/Users/santiago/Documents/CAMS_reanalysis/CAMS_AOD_2008.nc')
AOD_CAMS_daily = np.mean(CAMS_daily.variables['aod550'][13::8,:,:],axis=0)

lon_daily = CAMS_daily.variables['longitude'][:]
lat_daily = np.flip(CAMS_daily.variables['latitude'][:])
graph_map(variable=np.flip(AOD_CAMS_daily,axis=0),lon=lon_daily,lat=lat_daily,cmap=demuth.reversed(),vmax=.6,grid=(True),title_on=(True),title='CAMS-POLDER 2008',save=(True),  n_levels_plot=10)

import sys

from netCDF4 import Dataset
import numpy as np
import pandas as pd
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
import iris.plot as iplt


def moving_average_centered(a, n):
    return pd.Series(a).rolling(window=n, center=True).mean().to_numpy()

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
bwr=cm.get_cmap('bwr', 256)
#===Create Color bar====
def make_Ramp( ramp_colors ): 
    from colour import Color
    from matplotlib.colors import LinearSegmentedColormap

    color_ramp = LinearSegmentedColormap.from_list( 'my_list', [ Color( c1 ).rgb for c1 in ramp_colors ] )
    return color_ramp
custom_ramp = make_Ramp( ["#ffffff","#aec8e0","#9dafb2","#e2b85f","#f0ae46","#ea9842","#d36f2b","#a04617","#6c423d"] ) 


le_ext=Dataset('/Users/santiago/Documents/LE_outputs/LE_CALIPSO_V1_calipso_20080101_0200_state.nc')

calipso_ext=Dataset('/Users/santiago/Documents/CALIPSO/2008/CSO/CSO_calipso_2008-01-01T01-30-23.nc')

calipso=Dataset('/Users/santiago/Documents/CALIPSO/2008/CAL_LID_L2_05kmAPro-Standard-V4-20.2008-01-01T01-30-23ZN_Subset.hdf')
calipso.set_auto_mask(False)
elevation=moving_average_centered(calipso.variables['Surface_Elevation_Statistics'][:-1,2],20)

le_profile=le_ext.variables['ys'][:]
calipso_profile=calipso_ext.variables['AEC'][:-1,0:-1]
lon=calipso_ext.variables['longitude'][:-1]
lat=calipso_ext.variables['latitude'][:-1]

le_profile[np.isnan(calipso_profile)]=np.nan

nx=calipso_profile.shape[0]
nz=calipso_profile.shape[1]

cross=np.zeros([nx,nz+2])
cross[:,0]=0
cross[:,1]=elevation

xx=np.zeros([nx,nz+2])
lev_ht_1km=np.arange(0,nz)
for i in range(nz):
    cross[:,i+1]=elevation+lev_ht_1km[i]
    xx[:,i]=lon
vmin=0
vmax=1e-4
n_levels_plot=100
levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
norm = BoundaryNorm(levels, ncolors=256, clip=True)
fig, ax = plt.subplots()
ax.fill_between(xx[:,0],cross[:,0],cross[:,1],color='k')
#im=ax.contourf(xx[:,1:-1],cross[:,1:-1],column_1km,levels=levels,cmap='viridis',norm=norm)
im=ax.pcolor(xx[:,1:-1],cross[:,1:-1],calipso_profile,vmin=vmin,vmax=vmax,cmap='viridis',norm=norm)

ax.set_ylim(0,10)
ax.set_xlim(min(lon),max(lon))
lon_ticks=np.linspace(min(lon),max(lon), 10, dtype=int)
ax.set_xticks(lon_ticks)
ax.set_ylabel('Altitude Km')
ax.set_xlabel('Longitude ')
ax.set_title('AEC CALIOP')
ax.set_ylim(0,10)
# Set scond x-axis
ax2 = ax.twiny()

lat_ticks=np.linspace(min(lat), max(lat), 10, dtype=int)
ax2.set_xticks(lat_ticks)
ax2.set_xlim(min(lat),max(lat))

ax2.xaxis.set_ticks_position('bottom') # set the position of the second x-axis to bottom
ax2.xaxis.set_label_position('bottom') # set the position of the second x-axis to bottom
ax2.spines['bottom'].set_position(('outward', 36))
ax2.set_xlabel('Latitude')
#ax2.set_xlim(ax.get_xlim())

cax = fig.add_axes([ax.get_position().x1+0.02,
          ax.get_position().y0+0.09,
          0.08,
          ax.get_position().height-0.2])

plt.colorbar(im,cax=cax)        
plt.show()


vmin=0
vmax=1e-5
n_levels_plot=100
levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
norm = BoundaryNorm(levels, ncolors=256, clip=True)
fig, ax = plt.subplots()
ax.fill_between(xx[:,0],cross[:,0],cross[:,1],color='k')
#im=ax.contourf(xx[:,1:-1],cross[:,1:-1],column_1km,levels=levels,cmap='viridis',norm=norm)
im=ax.pcolor(xx[:,1:-1],cross[:,1:-1],le_profile,vmin=vmin,vmax=vmax,cmap='viridis',norm=norm)

ax.set_ylim(0,10)
ax.set_xlim(min(lon),max(lon))
lon_ticks=np.linspace(min(lon),max(lon), 10, dtype=int)
ax.set_xticks(lon_ticks)
ax.set_ylabel('Altitude Km')
ax.set_xlabel('Longitude ')
ax.set_title('AEC LE')
ax.set_ylim(0,10)
# Set scond x-axis
ax2 = ax.twiny()

lat_ticks=np.linspace(min(lat), max(lat), 10, dtype=int)
ax2.set_xticks(lat_ticks)
ax2.set_xlim(min(lat),max(lat))

ax2.xaxis.set_ticks_position('bottom') # set the position of the second x-axis to bottom
ax2.xaxis.set_label_position('bottom') # set the position of the second x-axis to bottom
ax2.spines['bottom'].set_position(('outward', 36))
ax2.set_xlabel('Latitude')
#ax2.set_xlim(ax.get_xlim())

cax = fig.add_axes([ax.get_position().x1+0.02,
          ax.get_position().y0+0.09,
          0.08,
          ax.get_position().height-0.2])

plt.colorbar(im,cax=cax)        
plt.show()

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
from datamanager import DataManager_LE
from matplotlib import cm


viridis = cm.get_cmap('viridis', 256)
blues = cm.get_cmap('Blues', 256)
#===Create Color bar====
def make_Ramp( ramp_colors ): 
    from colour import Color
    from matplotlib.colors import LinearSegmentedColormap

    color_ramp = LinearSegmentedColormap.from_list( 'my_list', [ Color( c1 ).rgb for c1 in ramp_colors ] )
    return color_ramp


custom_ramp = make_Ramp( ["#ffffff","#aec8e0","#9dafb2","#e2b85f","#f0ae46","#ea9842","#d36f2b","#a04617","#6c423d"] ) 

gridded_aod=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_POLDER_Complete/gridded/LE_2008.nc')
gridded_aaod=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_polder-aaod_2008.nc')

gridded_aaod.graph_map(variable='yr',biascorr=1,vmin=0,vmax=0.3,n_levels_plot=10,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='POLDER5 AAOD565')
gridded_aaod.graph_map(variable='ys',biascorr=1,vmin=0,vmax=0.03,n_levels_plot=10,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='LE AAOD565 ')

le_file=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_05_aaod/LE_ssa_200805.nc')
le_file.graph_map(variable='ssa_550nm',biascorr=1,vmin=0.7,vmax=1,n_levels_plot=10,cmap=viridis,type_graph='mean',date=False,save=True,title='LE output SSA 565')
le_file.graph_map(variable='aod_563nm',biascorr=1,vmin=0,vmax=0.5,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='LE output AOD 565')


ext_file=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_05_aaod/LE_ext_200805.nc')
ext_file.graph_map(variable='extinction_563nm',sum_level=True,biascorr=1,vmin=0,vmax=1e-3,n_levels_plot=20,cmap=blues,type_graph='mean',date=False,save=True,title='LE ext 565 nm')

lon=gridded_aod.nc.variables['longitude'][:]
lat=gridded_aod.nc.variables['latitude'][:]
polder_aod=gridded_aod.nc.variables['yr'][:]
polder_aod=np.nanmean( polder_aod[0:-1,:,:,0],axis=0)
polder_aod=np.squeeze(polder_aod)
polder_aaod=gridded_aaod.nc.variables['yr'][:]
polder_aaod=np.nanmean(polder_aaod[0:-1,:,:],axis=0)

polder_ssa=1-(polder_aaod/polder_aod)


le_aod=gridded_aod.nc.variables['ys'][:]
le_aod=np.nanmean( le_aod[0:-1,:,:],axis=0)
le_aod=np.squeeze(le_aod)

le_aaod=gridded_aaod.nc.variables['ys'][:]
le_aaod=np.nanmean( le_aaod[0:-1,:,:],axis=0)

le_ssa=1-(le_aaod/le_aod)
le_ssa.mask=polder_ssa.mask

vmin=0.7
vmax=1
n_levels_plot=10
cmap=viridis
biascorr=1
title='Polder SSA 565 nm'
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
ax.set_facecolor((1, 1, 1))

norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
plt.pcolormesh(lon, lat, polder_ssa, 
                   transform=ccrs.PlateCarree(),cmap=cmap, norm=norm,vmin=vmin, vmax=vmax,shading='auto')
european_borders=cfeature.NaturalEarthFeature(
         category='cultural',
         name='admin_0_countries',
         scale='50m',
         facecolor='none')


coastlines=cfeature.NaturalEarthFeature(
                 category='physical',
                 name='coastline',
                 scale='50m',
                 facecolor='none')
ax.add_feature(european_borders,edgecolor='black',linewidth=0.2)
ax.add_feature(coastlines,edgecolor='black',linewidth=1)
plt.xlim(lon[0],lon[-1])
plt.ylim(lat[0],lat[-1])
plt.title(title)
        

cax = fig.add_axes([ax.get_position().x1+0.02,
         ax.get_position().y0+0.09,
         0.08,
         ax.get_position().height-0.2])

plt.colorbar(cax=cax,extend='both')
plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")
plt.show()







fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
ax.set_facecolor((1, 1, 1))

norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
plt.pcolormesh(lon, lat,le_ssa, 
                   transform=ccrs.PlateCarree(),cmap=cmap, norm=norm,vmin=vmin, vmax=vmax,shading='auto')
european_borders=cfeature.NaturalEarthFeature(
         category='cultural',
         name='admin_0_countries',
         scale='50m',
         facecolor='none')


coastlines=cfeature.NaturalEarthFeature(
                 category='physical',
                 name='coastline',
                 scale='50m',
                 facecolor='none')
ax.add_feature(european_borders,edgecolor='black',linewidth=0.2)
ax.add_feature(coastlines,edgecolor='black',linewidth=1)
plt.xlim(lon[0],lon[-1])
plt.ylim(lat[0],lat[-1])
title='LE SSA 565 nm'
plt.title(title)
        

cax = fig.add_axes([ax.get_position().x1+0.02,
         ax.get_position().y0+0.09,
         0.08,
         ax.get_position().height-0.2])

plt.colorbar(cax=cax,extend='both')
plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")
plt.show()


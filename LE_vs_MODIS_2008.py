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
from matplotlib import cm


viridis = cm.get_cmap('viridis', 256)

#===Create Color bar====
def make_Ramp( ramp_colors ): 
    from colour import Color
    from matplotlib.colors import LinearSegmentedColormap

    color_ramp = LinearSegmentedColormap.from_list( 'my_list', [ Color( c1 ).rgb for c1 in ramp_colors ] )
    return color_ramp


custom_ramp = make_Ramp( ["#ffffff","#aec8e0","#9dafb2","#e2b85f","#f0ae46","#ea9842","#d36f2b","#a04617","#6c423d"] ) 


gridded_aqua=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_modis-aqua_2008.nc')
gridded_terra=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_modis-terra_2008.nc')
gridded_polder=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_polder-aod-563_2008.nc')

gridded_aqua.graph_map(variable='ys',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='LE-Aqua AOD550')
gridded_aqua.graph_map(variable='yr',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='MODIS-Aqua AOD550')

gridded_terra.graph_map(variable='ys',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='LE-Terra AOD550')
gridded_terra.graph_map(variable='yr',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='MODIS-Terra AOD550')


gridded_polder.graph_map(variable='ys',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='LE AOD565')
gridded_polder.graph_map(variable='yr',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='POLDER AOD565')




gridded_aqua.timeseries(variable=['ys','yr'],biascorr=[1,1],label=['LE ','MODIS-Aqua'],type_graph='mean',save=True,title='Aqua Domain Average time series',y_lim=[0,0.6],windows=4)
gridded_terra.timeseries(variable=['ys','yr'],biascorr=[1,1],label=['LE ','MODIS-Terra'],type_graph='mean',save=True,title='Terra Domain Average time series',y_lim=[0,0.6],windows=4)



gridded_polder.timeseries(variable=['ys','yr'],biascorr=[1,1],label=['LE ','POLDER'],type_graph='mean',save=True,title='POLDER Domain Average time series',y_lim=[0,2],windows=24)


gridded_terra.hist2d(biascorr=[1,1],cmap=white_viridis,xlim=[0.01,1.05],ylim=[0.01,1.05],log=True,save=True,title='Terra scatter')
gridded_aqua.hist2d(biascorr=[1,1],cmap=white_viridis,xlim=[0.01,1.05],ylim=[0.01,1.05],log=True,save=True,title='Agua scatter')
gridded_polder.hist2d(biascorr=[1,1],cmap=white_viridis,xlim=[0.01,1.05],ylim=[0.01,1.05],log=True,save=True,title='POLDER scatter')





lon=gridded_polder.nc.variables['longitude'][:]
lat=gridded_polder.nc.variables['latitude'][:]
polder_aod=gridded_polder.nc.variables['yr'][:]
polder_aod=np.nanmean( polder_aod[0:-1,:,:,0],axis=0)
modis_aod=gridded_aqua.nc.variables['yr'][:]
modis_aod=np.nanmean( modis_aod[0:-1,:,:],axis=0)
mask_2=np.logical_or(polder_aod.mask,modis_aod.mask)
polder_aod.mask=mask_2
modis_aod.mask=mask_2
vmin=0
vmax=0.7
n_levels_plot=20
cmap=custom_ramp
biascorr=1
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
ax.set_facecolor((0.6, 0.6, 0.6))

norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
plt.pcolormesh(lon, lat, polder_aod, 
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
title='Polder Colocated'
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
ax.set_facecolor((0.6, 0.6, 0.6))

norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
plt.pcolormesh(lon, lat, modis_aod, 
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
title='Modis Aqua Colocated'
plt.title(title)
        

cax = fig.add_axes([ax.get_position().x1+0.02,
         ax.get_position().y0+0.09,
         0.08,
         ax.get_position().height-0.2])

plt.colorbar(cax=cax,extend='both')
plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")

plt.show()

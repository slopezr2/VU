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
from matplotlib.colors import TwoSlopeNorm
from os import listdir
from os.path import isfile, join
import matplotlib.dates as mdates
from matplotlib.colors import LinearSegmentedColormap
import mpl_scatter_density
from matplotlib import cm
import pandas as pd
sys.path.append('./libraries/')
from datamanager import DataManager_GRIDDED
from datamanager import DataManager_LE
from matplotlib import cm


#====Graph map===

def graph_map(lon=0,lat=0,vmin=0,vmax=0.3,vcenter=0,n_levels_plot=20,cmap=1,variable=1,title='AOD',save=True,extend='max',norm2=False,grid=False):
    

    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
    ax.set_facecolor((0.6, 0.6, 0.6))
    if norm2==True:
        norm=TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    else:    
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    plt.pcolormesh(lon, lat, variable, 
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
    if grid==True:
        gl = ax.gridlines(draw_labels=True,color='black', alpha=0.5,linestyle='--')
        gl.top_labels = False
        gl.right_labels = False
    plt.xlim(lon[0],lon[-1])
    plt.ylim(lat[0],lat[-1])
    
    plt.title(title)
    cax = fig.add_axes([ax.get_position().x1+0.02,
              ax.get_position().y0+0.09,
              0.08,
              ax.get_position().height-0.2])
    cbar=plt.colorbar(cax=cax,extend=extend)
    cbar.formatter.set_powerlimits((-3, 3))
    cbar.update_ticks()
    cbar.ax.set_title('label', y=1.05, x=0.5, rotation=0)

    if save:        
        plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")
    
    plt.show()



#==Center histogram bins==
def bins_labels(bins,label, **kwargs):
    bin_w = (max(bins) - min(bins)) / (len(bins) - 1)
    plt.xticks(np.arange(min(bins)+bin_w/2, max(bins)+1, bin_w),labels=label, **kwargs)
    plt.xlim(bins[0], bins[-1])

#===Create Color bar====
def make_Ramp( ramp_colors ): 
    from colour import Color
    from matplotlib.colors import LinearSegmentedColormap

    color_ramp = LinearSegmentedColormap.from_list( 'my_list', [ Color( c1 ).rgb for c1 in ramp_colors ] )
    return color_ramp


custom_ramp = make_Ramp( ["#ffffff","#aec8e0","#9dafb2","#e2b85f","#f0ae46","#ea9842","#d36f2b","#a04617","#6c423d"] ) 
viridis = cm.get_cmap('viridis', 256)
RdYlGn=cm.get_cmap('Spectral', 256)
oranges=cm.get_cmap('Oranges', 256)
bwr=cm.get_cmap('bwr', 256)


aod_dry_post_datamanager=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_snellius_post_dry_2008.nc')
aod_dry_complete_datamanager=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_snellius_scheme_dry_2008.nc')
aod_dry_post_cartesius_datamanager=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_dry_aod_2008.nc')


aod_dry_post=aod_dry_post_datamanager.nc.variables['aod_563nm'][:]
aod_dry_complete=aod_dry_complete_datamanager.nc.variables['aod_563nm'][:]
aod_dry_post_cartesius=aod_dry_post_cartesius_datamanager.nc.variables['aod_563nm'][:]


lon=aod_dry_complete_datamanager.nc.variables['longitude'][:]
lat=aod_dry_complete_datamanager.nc.variables['latitude'][:]

graph_map(lon=lon,lat=lat,variable=((aod_dry_complete-aod_dry_post)/aod_dry_complete).mean(axis=0)*100,vmin=-0.5,vmax=0.5,n_levels_plot=20,cmap=bwr,title='LE complete-Post tool AOD 563nm 2008',save=True,extend='max',grid=True,norm2=True,vcenter=0)

graph_map(lon=lon,lat=lat,variable=((aod_dry_post-aod_dry_post_cartesius)/aod_dry_post).mean(axis=0)*100,vmin=-0.0005,vmax=0.0005,n_levels_plot=20,cmap=bwr,title='Snellius-Cartesius AOD 563nm 2008',save=True,extend='max',grid=True,norm2=True,vcenter=0)

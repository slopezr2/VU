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
bwr=cm.get_cmap('bwr', 256)
#===Create Color bar====
def make_Ramp( ramp_colors ): 
    from colour import Color
    from matplotlib.colors import LinearSegmentedColormap

    color_ramp = LinearSegmentedColormap.from_list( 'my_list', [ Color( c1 ).rgb for c1 in ramp_colors ] )
    return color_ramp

path_calipso='/Users/santiago/Documents/CALIPSO/2008/05/'
onlyfiles = [f for f in listdir(path_calipso) if isfile(join(path_calipso, f))]
onlyfiles.sort()
aeronet={}




vmin=0
vmax=0.5
n_levels_plot=10
#cmap=viridis
cmap=RdYlGn.reversed()
biascorr=1
title='Calipso Column Optical Depth Stratospheric'
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
ax.set_facecolor((1, 1, 1))

norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
for file in onlyfiles:
    if file[-3:]=='.nc':
        calipso=Dataset(path_calipso+file)
        calipso.set_auto_mask(False)
        custom_ramp = make_Ramp( ["#ffffff","#aec8e0","#9dafb2","#e2b85f","#f0ae46","#ea9842","#d36f2b","#a04617","#6c423d"] ) 
        
        column=calipso.variables['Column_Optical_Depth_Tropospheric_Aerosols_532'][:]
        column=np.ma.masked_array(column, mask=((column<=0)|(column>3)))
        
        
        lon=calipso.variables['Longitude'][:,1]
        lat=calipso.variables['Latitude'][:,1]
        
        plt.scatter(lon, lat,3 ,column, 's',
                           transform=ccrs.PlateCarree(),cmap=cmap, norm=norm,vmin=vmin, vmax=vmax)
        


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
#ax.stock_img()
plt.xlim(-15,35)
plt.ylim(35,70)
plt.title(title)
        

cax = fig.add_axes([ax.get_position().x1+0.02,
         ax.get_position().y0+0.09,
         0.08,
         ax.get_position().height-0.2])

plt.colorbar(cax=cax,extend='both')
#plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")
plt.show()






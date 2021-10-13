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


custom_ramp = make_Ramp( ["#ffffff","#aec8e0","#9dafb2","#e2b85f","#f0ae46","#ea9842","#d36f2b","#a04617","#6c423d"] ) 

gridded_440=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_aeronet_440_2008.nc')
gridded_500=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_aeronet_500_2008.nc')
gridded_870=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_aeronet_870_2008.nc')

gridded_440.hist2d(biascorr=[1,1],cmap=white_viridis,xlim=[0.01,1.05],ylim=[0.01,1.05],log=True,save=True,title='AERONET 440 scatter',xlabel='AERONET 440')
gridded_500.hist2d(biascorr=[1,1],cmap=white_viridis,xlim=[0.0,0.3],ylim=[0.0,0.3],log=False,save=False,title='AERONET 500 scatter',xlabel='AERONET 500')
# gridded_870.hist2d(biascorr=[1,1],cmap=white_viridis,xlim=[0.01,1.05],ylim=[0.01,1.05],log=True,save=True,title='AERONET 870 scatter',xlabel='AERONET 870')

gridded_440.graph_map(variable='yr',biascorr=1,vmin=0,vmax=0.5,n_levels_plot=10,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='AERONET AOD440')
# gridded_440.graph_map(variable='ys',biascorr=1,vmin=0,vmax=0.5,n_levels_plot=10,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='LE AOD440')

gridded_500.graph_map(variable='yr',biascorr=1,vmin=0,vmax=0.7,n_levels_plot=10,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='AERONET AOD500')
gridded_500.graph_diff_map(facecolor=(1,1,1),size_scatter=20,scatter=True,dif_porcentage=True,biascorr=1,vmin=-100,vmax=100,n_levels_plot=20,cmap=bwr,date=False,save=True,title='Percentage diff obs-mod AERONET AOD500')



# gridded_500.graph_map(variable='ys',biascorr=1,vmin=0,vmax=0.5,n_levels_plot=10,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='LE AOD500')

gridded_870.graph_map(variable='yr',biascorr=1,vmin=0,vmax=0.5,n_levels_plot=10,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='AERONET AOD870')
# gridded_870.graph_map(variable='ys',biascorr=1,vmin=0,vmax=0.5,n_levels_plot=10,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='LE AOD870')


lon=gridded_440.nc.variables['lon'][:]
lat=gridded_440.nc.variables['lat'][:]

aeronet_440_aod=gridded_440.nc.variables['yr'][:]
aeronet_870_aod=gridded_870.nc.variables['yr'][:]
aeronet_aexp_scatter=-np.log(aeronet_440_aod/aeronet_870_aod)/np.log(440/870)

aeronet_440_aod[aeronet_440_aod<0]=np.nan
aeronet_870_aod[aeronet_870_aod<0]=np.nan

aeronet_440_aod=np.nanmean( aeronet_440_aod,axis=0)
aeronet_870_aod=np.nanmean( aeronet_870_aod,axis=0)
aeronet_aexp=-np.log(aeronet_440_aod/aeronet_870_aod)/np.log(440/870)


le_440_aod=gridded_440.nc.variables['ys'][:]
le_870_aod=gridded_870.nc.variables['ys'][:]
le_aexp_scatter=-np.log(le_440_aod/le_870_aod)/np.log(440/870)

le_870_aod=np.nanmean( le_870_aod[0:-1,:,:],axis=0)
le_440_aod=np.nanmean( le_440_aod[0:-1,:,:],axis=0)
le_aexp=-np.log(le_440_aod/le_870_aod)/np.log(440/870)


gridded_500.graph_map(variable='yr',biascorr=(500/550)**-aeronet_aexp,vmin=0,vmax=0.7,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='AERONET AOD550 extrapolated')



vmin=0
vmax=3
n_levels_plot=10
#cmap=viridis
cmap=RdYlGn.reversed()
biascorr=1
title='Aeronet Agnstrom Exponent 440-870'
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
ax.set_facecolor((1, 1, 1))

norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
plt.pcolormesh(lon, lat, aeronet_aexp, 
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
plt.pcolormesh(lon, lat,le_aexp, 
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
title='LE Angstrom Exponent 440-870'
plt.title(title)
        

cax = fig.add_axes([ax.get_position().x1+0.02,
         ax.get_position().y0+0.09,
         0.08,
         ax.get_position().height-0.2])

plt.colorbar(cax=cax,extend='both')
plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")
plt.show()




le_aexp_scatter=le_aexp_scatter.reshape(le_aexp_scatter.shape[0]*le_aexp_scatter.shape[1]*le_aexp_scatter.shape[2])
aeronet_aexp_scatter=aeronet_aexp_scatter.reshape(aeronet_aexp_scatter.shape[0]*aeronet_aexp_scatter.shape[1]*aeronet_aexp_scatter.shape[2])

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection='scatter_density')
plt.grid(True, which="both", ls="-",alpha=0.3)
density = ax.scatter_density((aeronet_aexp_scatter), (le_aexp_scatter),cmap=white_viridis)
xx=np.arange(0,100)
plt.plot(xx, xx,'k-',linewidth=1)
plt.plot(xx, xx*0.5,'k-',linewidth=0.5)
plt.plot(xx, xx*2,'k-',linewidth=0.5)
fig.colorbar(density, label='Number of points per pixel')
ax.set_xlabel('Aeronet Aexp 440-870')
ax.set_ylabel('LE Aexp 440-870')
 
plt.xlim([0,3])
plt.ylim([0,3])
 
plt.savefig('./Figures/Aeronet_Aexp_scatter.png',format='png', dpi=1000,bbox_inches = "tight")

plt.show()
  
 


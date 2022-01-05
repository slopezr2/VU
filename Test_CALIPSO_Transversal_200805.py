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

#path_calipso='/Users/santiago/Documents/CALIPSO/2008/05/'
path_calipso='/Users/santiago/Documents/CALIPSO/2008/'
onlyfiles = [f for f in listdir(path_calipso) if isfile(join(path_calipso, f))]
onlyfiles.sort()

orbit=10

calipso=Dataset(path_calipso+onlyfiles[orbit])
calipso.set_auto_mask(False)


CAD=np.flip(calipso.variables['CAD_Score'][:,:,0],1)
CAD=CAD[:,:344]

qs=np.flip(calipso.variables['Extinction_Coefficient_Uncertainty_532'][:],1)
qs=qs[:,:344]
qs[qs<0]=10000
qc=np.flip(calipso.variables['Extinction_QC_Flag_532'][:],1)
qc=qc[:,:344,0]

column2=np.flip(calipso.variables['Extinction_Coefficient_532'][:],1) 
column2=column2[:,:344]#To take just the troposphere
#column2=np.ma.masked_array(column2, mask=((column2<=0)|(column2>1)| (CAD>-80) | (qs>column2*0.50)| ((qc!=1) & (qc!=0))))
column2=np.ma.masked_array(column2, mask=((column2<=0)|(column2>1)| (CAD>-80)| (CAD<-100) | (qs>column2*0.99)))
#column2=np.ma.masked_array(column2, mask=((column2<=0)|(column2>1.25)))

column=column2.data.copy()
column[column2.mask==True]=np.nan
pressure=np.flip(calipso.variables['Pressure'][:],1) 
pressure=pressure[:,:344]#To take just the troposphere
pressure=np.ma.masked_array(pressure, mask=((pressure<=0)))
inter_levels=17
column_1km=np.zeros([column.shape[0],column.shape[1]//inter_levels+1])
pressure_1km=np.zeros([pressure.shape[0],pressure.shape[1]//inter_levels+1])

for i in range(column.shape[0]):
    aux=pd.Series(column[i,:]).groupby(pd.Series(column[i,:]).index//inter_levels).mean().to_numpy()
    column_1km[i,:aux.shape[0]]= pd.Series(column[i,:]).groupby(pd.Series(column[i,:]).index//inter_levels).mean().to_numpy()
    pressure_1km[i,:aux.shape[0]]= pd.Series(pressure[i,:]).groupby(pd.Series(pressure[i,:]).index//inter_levels).mean().to_numpy()[:column_1km.shape[1]]
    

df=pd.DataFrame(pressure_1km)
#pressure_1km=df.fillna(df.mean(axis=0)).to_numpy()

lon=calipso.variables['Longitude'][:,1]
lat=calipso.variables['Latitude'][:,1]

lev_ht=np.arange(0,column.shape[1])#*60/1000
lev_ht_1km=np.arange(0,column_1km.shape[1])*60*inter_levels/1000 -0.5
elevation=calipso.variables['Surface_Elevation_Statistics'][:,2]

nx=column_1km.shape[0]
nz=column_1km.shape[1]

cross=np.zeros([nx,nz+2])
cross[:,0]=0
cross[:,1]=elevation

xx=np.zeros([nx,nz+2])
for i in range(nz):
    #cross[:,i+1]=elevation+lev_ht_1km[i]
    cross[:,i+1]=lev_ht_1km[i]
    xx[:,i]=lon
# plot

vmin=0
vmax=1e-1
n_levels_plot=256
levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
norm = BoundaryNorm(levels, ncolors=256, clip=True)
fig, ax = plt.subplots()

#im=ax.contourf(xx[:,1:-1],cross[:,1:-1],column_1km,levels=levels,cmap='viridis',norm=norm)
im=ax.pcolor(xx[:,1:-1],cross[:,1:-1],column_1km,vmin=vmin,vmax=vmax,cmap='viridis',norm=norm)
ax.fill_between(xx[:,0],cross[:,0],elevation,color='k')
#ax.set_ylim(0,10)
ax.set_xlim(min(lon),max(lon))
lon_ticks=np.linspace(min(lon),max(lon), 10, dtype=int)
ax.set_xticks(lon_ticks)
ax.set_ylabel('Altitude Km')
ax.set_xlabel('Longitude ')
ax.set_title('AEC CALIOP')
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
#plt.savefig('./Figures/Orbit_profile.png',format='png', dpi=1000,bbox_inches = "tight")      
plt.show()


# vmin=0
# vmax=1e-1
# n_levels_plot=256
# levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
# norm = BoundaryNorm(levels, ncolors=256, clip=True)
# fig, ax = plt.subplots()

# im=ax.pcolor(xx[:,1:-1],cross[:,1:-1],qs[:,:],cmap='viridis',norm=norm)
# ax.fill_between(xx[:,0],cross[:,0],elevation,color='k')
# #ax.set_ylim(0,10)
# ax.set_xlim(min(lon),max(lon))
# lon_ticks=np.linspace(min(lon),max(lon), 10, dtype=int)
# ax.set_xticks(lon_ticks)
# ax.set_ylabel('Altitude Km')
# ax.set_xlabel('Longitude ')
# ax.set_title('Pressure CALIOP')
# # Set scond x-axis
# ax2 = ax.twiny()

# lat_ticks=np.linspace(min(lat), max(lat), 10, dtype=int)
# ax2.set_xticks(lat_ticks)
# ax2.set_xlim(min(lat),max(lat))

# ax2.xaxis.set_ticks_position('bottom') # set the position of the second x-axis to bottom
# ax2.xaxis.set_label_position('bottom') # set the position of the second x-axis to bottom
# ax2.spines['bottom'].set_position(('outward', 36))
# ax2.set_xlabel('Latitude')
# #ax2.set_xlim(ax.get_xlim())

# cax = fig.add_axes([ax.get_position().x1+0.02,
#           ax.get_position().y0+0.09,
#           0.08,
#           ax.get_position().height-0.2])

# plt.colorbar(im,cax=cax)        
# plt.show()


vmin=0
vmax=0.5
n_levels_plot=10
#cmap=viridis
cmap=RdYlGn.reversed()
biascorr=1
title='Orbit'
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
ax.set_facecolor((1, 1, 1))
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
plt.scatter(lon, lat, 3 ,10*np.ones(lat.shape), 's',
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
ax.stock_img()
plt.xlim(-15,35)
plt.ylim(35,70)
plt.title(title)
        
# cax = fig.add_axes([ax.get_position().x1+0.02,
#          ax.get_position().y0+0.09,
#          0.08,
#          ax.get_position().height-0.2])

#plt.colorbar(cax=cax,extend='both')
#plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")
plt.show()



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
#===Create Color bar====
def make_Ramp( ramp_colors ): 
    from colour import Color
    from matplotlib.colors import LinearSegmentedColormap

    color_ramp = LinearSegmentedColormap.from_list( 'my_list', [ Color( c1 ).rgb for c1 in ramp_colors ] )
    return color_ramp


custom_ramp = make_Ramp( ["#ffffff","#aec8e0","#9dafb2","#e2b85f","#f0ae46","#ea9842","#d36f2b","#a04617","#6c423d"] ) 

emep_pm10=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_emep-pm10_2008.nc')

# emep_pm10.hist2d(biascorr=[1e9,1],cmap=white_viridis,xlim=[0.01,1.05],ylim=[0.01,1.05],log=True,save=True,title='PM$_{10}$',xlabel='EMEP')

# emep_pm10.graph_map(variable='yr',biascorr=1,vmin=0,vmax=50,n_levels_plot=10,cmap=custom_ramp,type_graph='mean',date=False,save=False,title='EMEP PM$_{10}$ [ug/m3]')
# emep_pm10.graph_map(variable='ys',biascorr=1e9,vmin=0,vmax=50,n_levels_plot=10,cmap=custom_ramp,type_graph='mean',date=False,save=False,title='LE PM$_{10}$ [ug/m3]')

lat=emep_pm10.nc.variables['lat'][:]
lon=emep_pm10.nc.variables['lon'][:]

i=0
conc_year=[]
for month in ['01','02','03','04','05','06','07','08','09','10','11','12']: #Abril a septiembre
    i=i+1
    print(month)
    file_conc='/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_conc_2008'+month+'.nc'
    conc_file_month=Dataset(file_conc)
    conc_month=conc_file_month.variables['tpm10'][:,1,:,:]
    conc_month=np.squeeze(np.nanmean(conc_month,axis=0))
    

    conc_year.append(conc_month)        
conc_year=np.array(conc_year)
conc_year=np.nanmean(conc_year,axis=0)*1e9 #En gramos

EMEP_LE={}
LE_pm10=[]
EMEP_pm10=[]
lon_stations=[]
lat_stations=[]
for station in EMEP.keys():
    lat_sta=EMEP[station]['latitude']
    lon_sta=EMEP[station]['longitude']
    for i in range(len(lat)-2):
        if lat_sta>=lat[i] and lat_sta<=lat[i+1]:
            lat_le_sta=i
    for i in range(len(lon)-2): 
        if lon_sta>=lon[i] and lon_sta<=lon[i+1]:
            lon_le_sta=i
    lon_stations.append(lon[lon_le_sta])
    lat_stations.append(lat[lat_le_sta])
    EMEP_LE[station]={}
    EMEP_LE[station]['LE']=conc_year[lat_le_sta,lon_le_sta]
    LE_pm10.append(conc_year[lat_le_sta,lon_le_sta])
    EMEP_LE[station]['EMEP']=EMEP[station]['values']['PM10'][EMEP[station]['values']['PM10']<230].mean()
    EMEP_pm10.append(EMEP[station]['values']['PM10'][(EMEP[station]['values']['PM10']<230)].mean())

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plt.grid(True, which="both", ls="-",alpha=0.3)
title='EMEP scatterplot'
plt.scatter((EMEP_pm10), (np.array(LE_pm10)))
xx=np.arange(0,100)
plt.xlim([0,40])
plt.ylim([0,40])
plt.plot(xx, xx,'k-',linewidth=1)
plt.plot(xx, xx*0.5,'k-',linewidth=0.5)
plt.plot(xx, xx*2,'k-',linewidth=0.5)
ax.set_xlabel(r'EMEP PM$_{10}$')
ax.set_ylabel(r'LE PM$_{10}$')
plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")
plt.show()



vmin=0
vmax=40
n_levels_plot=10
#cmap=viridis
cmap=RdYlGn.reversed()
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
ax.set_facecolor((1, 1, 1))

norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
plt.scatter(lon_stations, lat_stations,15,np.array(LE_pm10)
                   ,cmap=cmap, norm=norm,vmin=vmin, vmax=vmax)
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
title=r'LE PM$_{10}$ $\mu$g/m$^3$'
plt.title(title)
        

cax = fig.add_axes([ax.get_position().x1+0.02,
         ax.get_position().y0+0.09,
         0.08,
         ax.get_position().height-0.2])

plt.colorbar(cax=cax,extend='both')
plt.savefig('./Figures/Map_LE_pm10.png',format='png', dpi=1000,bbox_inches = "tight")
plt.show()


vmin=0
vmax=40
n_levels_plot=10
#cmap=viridis
cmap=RdYlGn.reversed()
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
ax.set_facecolor((1, 1, 1))

norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
plt.scatter(lon_stations, lat_stations,15,EMEP_pm10
                   ,cmap=cmap, norm=norm,vmin=vmin, vmax=vmax)
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
title=r'EMEP PM$_{10}$ $\mu$g/m$^3$'
plt.title(title)
        

cax = fig.add_axes([ax.get_position().x1+0.02,
         ax.get_position().y0+0.09,
         0.08,
         ax.get_position().height-0.2])

plt.colorbar(cax=cax,extend='both')
plt.savefig('./Figures/Map_EMEP_pm10.png',format='png', dpi=1000,bbox_inches = "tight")
plt.show()

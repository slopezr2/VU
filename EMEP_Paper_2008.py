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
import os
import matplotlib.dates as mdates
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import ListedColormap
import mpl_scatter_density
from matplotlib import cm
from matplotlib.colors import TwoSlopeNorm
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
bwr=cm.get_cmap('bwr', 256)
#===Load user libraries====
sys.path.append('./libraries/')
from datamanager import DataManager_GRIDDED
from matplotlib import cm



emep_pm10=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_emep-pm10_2008.nc')
lat=emep_pm10.nc.variables['lat'][:]
lon=emep_pm10.nc.variables['lon'][:]

i=0
conc_year_pm10=[]
conc_year_pm25=[]
for month in ['01','02','03','04','05','06','07','08','09','10','11','12']: #Abril a septiembre
    i=i+1
    print(month)
    file_conc='/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_conc_2008'+month+'.nc'
    conc_file_month=Dataset(file_conc)
    conc_month_pm10=conc_file_month.variables['tpm10'][:,1,:,:]
    conc_month_pm10=np.squeeze(np.nanmean(conc_month_pm10,axis=0))
    conc_month_pm25=conc_file_month.variables['tpm25'][:,1,:,:]
    conc_month_pm25=np.squeeze(np.nanmean(conc_month_pm25,axis=0))

    conc_year_pm10.append(conc_month_pm10) 
    conc_year_pm25.append(conc_month_pm25)        
conc_year_pm10=np.array(conc_year_pm10)
conc_year_pm10=np.nanmean(conc_year_pm10,axis=0)*1e9 #En gramos
conc_year_pm25=np.array(conc_year_pm25)
conc_year_pm25=np.nanmean(conc_year_pm25,axis=0)*1e9 #En gramos

EMEP_LE_pm10={}
LE_pm10=[]
EMEP_pm10=[]
EMEP_LE_pm25={}
LE_pm25=[]
EMEP_pm25=[]
lon_stations=[]
lat_stations=[]
lon_stations_pm25=[]
lat_stations_pm25=[]
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
    EMEP_LE_pm10[station]={}
    EMEP_LE_pm10[station]['LE']=conc_year_pm10[lat_le_sta,lon_le_sta]
    LE_pm10.append(conc_year_pm10[lat_le_sta,lon_le_sta])
    EMEP_LE_pm10[station]['EMEP']=EMEP[station]['values']['PM10'][EMEP[station]['values']['PM10']<230].mean()
    EMEP_pm10.append(EMEP[station]['values']['PM10'][(EMEP[station]['values']['PM10']<230)].mean())
    if 'PM2.5' in EMEP[station]['values'].columns:
        lon_stations_pm25.append(lon[lon_le_sta])
        lat_stations_pm25.append(lat[lat_le_sta])
        EMEP_LE_pm25[station]={}
        EMEP_LE_pm25[station]['LE']=conc_year_pm25[lat_le_sta,lon_le_sta]
        LE_pm25.append(conc_year_pm25[lat_le_sta,lon_le_sta])
        EMEP_LE_pm25[station]['EMEP']=EMEP[station]['values']['PM2.5'][EMEP[station]['values']['PM2.5']<230].mean()
        EMEP_pm25.append(EMEP[station]['values']['PM2.5'][(EMEP[station]['values']['PM2.5']<230)].mean())

EMEP_pm10 = np.array(EMEP_pm10)
LE_pm10 = np.array(LE_pm10)
EMEP_pm25 = np.array(EMEP_pm25)
LE_pm25 = np.array(LE_pm25)


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plt.grid(True, which="both", ls="-",alpha=0.3)
title=r'EMEP scatterplot PM$_{10}$'
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

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plt.grid(True, which="both", ls="-",alpha=0.3)
title=r'EMEP scatterplot PM$_{2.5}$'
plt.scatter((EMEP_pm25), (np.array(LE_pm25)))
xx=np.arange(0,100)
plt.xlim([0,20])
plt.ylim([0,20])
plt.plot(xx, xx,'k-',linewidth=1)
plt.plot(xx, xx*0.5,'k-',linewidth=0.5)
plt.plot(xx, xx*2,'k-',linewidth=0.5)
ax.set_xlabel(r'EMEP PM$_{2.5}$')
ax.set_ylabel(r'LE PM$_{2.5}$')
plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")
plt.show()


#===Diff EMEP-LE PM10===
vmin=-100
vmax=100
vcenter=0
norm=TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
size_stations = 20

cmap=bwr
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
ax.set_facecolor((1, 1, 1))


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
plt.scatter(lon_stations, lat_stations,size_stations,np.array((LE_pm10-EMEP_pm10)*100/EMEP_pm10)
                   ,cmap=cmap, norm=norm,edgecolors=(0.1,0.1,0.1),linewidths=0.5)
plt.xlim(lon[0],lon[-1])
plt.ylim(lat[0],lat[-1])
title=r'Relative difference LE-EMEP PM$_{10}$'
plt.title(title)
gl = ax.gridlines(draw_labels=True,color='black', alpha=0.5,linestyle='--',zorder=100000)
gl.top_labels = False
gl.right_labels = False        

cax = fig.add_axes([ax.get_position().x1+0.02,
         ax.get_position().y0+0.09,
         0.08,
         ax.get_position().height-0.2])

plt.colorbar(cax=cax,extend='both')
plt.savefig('./Figures/Map_LE_pm10.png',format='png', dpi=1000,bbox_inches = "tight")
plt.show()

#===Diff EMEP-LE PM25===
vmin=-100
vmax=100
vcenter=0
norm=TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
size_stations = 20

cmap=bwr
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
ax.set_facecolor((1, 1, 1))


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
plt.scatter(lon_stations_pm25, lat_stations_pm25,size_stations,np.array((LE_pm25-EMEP_pm25)*100/EMEP_pm25)
                   ,cmap=cmap, norm=norm,edgecolors=(0.1,0.1,0.1),linewidths=0.5)
plt.xlim(lon[0],lon[-1])
plt.ylim(lat[0],lat[-1])
title=r'Relative difference LE-EMEP PM$_{2.5}$'
plt.title(title)
gl = ax.gridlines(draw_labels=True,color='black', alpha=0.5,linestyle='--',zorder=100000)
gl.top_labels = False
gl.right_labels = False        

cax = fig.add_axes([ax.get_position().x1+0.02,
         ax.get_position().y0+0.09,
         0.08,
         ax.get_position().height-0.2])

plt.colorbar(cax=cax,extend='both')
plt.savefig('./Figures/Map_LE_pm25.png',format='png', dpi=2500,bbox_inches = "tight")
plt.show()

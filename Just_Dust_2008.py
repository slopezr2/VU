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
from matplotlib import cm

sys.path.append('./libraries/')
from datamanager import DataManager_GRIDDED
from datamanager import DataManager_LE
from matplotlib import cm

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


gridded_dust_aod=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_dust_aod_2008.nc')


gridded_dust_aod.graph_map(variable='aod_563nm',biascorr=1,vmin=0,vmax=0.1,n_levels_plot=20,level=0,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='AOD 563nm Just Dust')
aod_year=gridded_dust_aod.nc.variables['aod_563nm'][:,:,:] #Todo el ano
aod_year=np.squeeze(np.nanmean(aod_year,axis=0))
h_file=Dataset('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_h_2008.nc')
i=0

burden_year=[]
burden_year_ccc=[]
burden_year_cc=[]
burden_year_c=[]
burden_year_f=[]
burden_year_ff=[]
for month in ['01','02','03','04','05','06','07','08','09','10','11','12']: #Abril a septiembre
    i=i+1
    print(month)
    file_conc='/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_conc_2008'+month+'.nc'
    conc_file_month=Dataset(file_conc)
    conc_month=conc_file_month.variables['dust_ccc'][:]+conc_file_month.variables['dust_cc'][:]+conc_file_month.variables['dust_c'][:]+conc_file_month.variables['dust_f'][:]+conc_file_month.variables['dust_ff'][:]
    conc_month=np.squeeze(np.nanmean(conc_month,axis=0))
    
    conc_month_ccc=conc_file_month.variables['dust_ccc'][:]
    conc_month_ccc=np.squeeze(np.nanmean(conc_month_ccc,axis=0))
    conc_month_cc=conc_file_month.variables['dust_cc'][:]
    conc_month_cc=np.squeeze(np.nanmean(conc_month_cc,axis=0))
    conc_month_c=conc_file_month.variables['dust_c'][:]
    conc_month_c=np.squeeze(np.nanmean(conc_month_c,axis=0))
    conc_month_f=conc_file_month.variables['dust_f'][:]
    conc_month_f=np.squeeze(np.nanmean(conc_month_f,axis=0))
    conc_month_ff=conc_file_month.variables['dust_ff'][:]
    conc_month_ff=np.squeeze(np.nanmean(conc_month_ff,axis=0))
    
    if i==3: #Marzo, sumo solo 29 dias
        ini_dia=ini_dia+29
        end_dia=end_dia+31
    elif i==2:
        ini_dia=ini_dia+31
        end_dia=end_dia+29
    elif i==1:
        ini_dia=0
        end_dia=31
    elif i in [4,6,9,11]:
        ini_dia=ini_dia+31
        end_dia=end_dia+30
    elif i ==8:   
        ini_dia=ini_dia+31
        end_dia=end_dia+31
    else:
        ini_dia=ini_dia+30
        end_dia=end_dia+31
        
    if end_dia>357:
        end_dia=-1
    h_month=h_file.variables['h'][ini_dia*24:end_dia*24,:,:,:]
    h_month=np.squeeze(np.nanmean(h_month,axis=0))
    burden=0
    burden_ccc=0
    burden_cc=0
    burden_c=0
    burden_f=0
    burden_ff=0
    for j in range(12):
        if j>0:
            burden=burden+conc_month[j,:,:]*(h_month[j,:,:]-h_month[j-1,:,:])
            burden_ccc=burden_ccc+conc_month_ccc[j,:,:]*(h_month[j,:,:]-h_month[j-1,:,:])
            burden_cc=burden_cc+conc_month_cc[j,:,:]*(h_month[j,:,:]-h_month[j-1,:,:])
            burden_c=burden_c+conc_month_c[j,:,:]*(h_month[j,:,:]-h_month[j-1,:,:])
            burden_f=burden_f+conc_month_f[j,:,:]*(h_month[j,:,:]-h_month[j-1,:,:])
            burden_ff=burden_ff+conc_month_ff[j,:,:]*(h_month[j,:,:]-h_month[j-1,:,:])
        else:
            burden=burden+conc_month[j,:,:]*(h_month[j,:,:])
            burden_ccc=burden_ccc+conc_month_ccc[j,:,:]*(h_month[j,:,:])
            burden_cc=burden_cc+conc_month_cc[j,:,:]*(h_month[j,:,:])
            burden_c=burden_c+conc_month_c[j,:,:]*(h_month[j,:,:])
            burden_f=burden_f+conc_month_f[j,:,:]*(h_month[j,:,:])
            burden_ff=burden_ff+conc_month_ff[j,:,:]*(h_month[j,:,:])
    burden_year.append(burden)
    burden_year_ccc.append(burden_ccc)
    burden_year_cc.append(burden_cc)
    burden_year_c.append(burden_c)
    burden_year_f.append(burden_f)
    burden_year_ff.append(burden_ff)
        
burden_year=np.array(burden_year)
burden_year=np.nanmean(burden_year,axis=0)*1e3 #En gramos
burden_year_ccc=np.array(burden_year_ccc)
burden_year_ccc=np.nanmean(burden_year_ccc,axis=0)*1e3 #En gramos
burden_year_cc=np.array(burden_year_cc)
burden_year_cc=np.nanmean(burden_year_cc,axis=0)*1e3 #En gramos
burden_year_c=np.array(burden_year_c)
burden_year_c=np.nanmean(burden_year_c,axis=0)*1e3 #En gramos
burden_year_f=np.array(burden_year_f)
burden_year_f=np.nanmean(burden_year_f,axis=0)*1e3 #En gramos
burden_year_ff=np.array(burden_year_ff)
burden_year_ff=np.nanmean(burden_year_ff,axis=0)*1e3 #En gramos


MEC=aod_year/burden_year





lon=gridded_dust_aod.nc.variables['lon'][:]
lat=gridded_dust_aod.nc.variables['lat'][:]
vmin=0
vmax=0.15
n_levels_plot=20
cmap=custom_ramp
biascorr=1

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
ax.set_facecolor((0.6, 0.6, 0.6))

norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
plt.pcolormesh(lon, lat, burden_year, 
                    transform=ccrs.PlateCarree(),cmap=RdYlGn.reversed(), norm=norm,vmin=vmin, vmax=vmax,shading='auto')
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
title='Dust Burden g m-2'
plt.title(title)
cax = fig.add_axes([ax.get_position().x1+0.02,
          ax.get_position().y0+0.09,
          0.08,
          ax.get_position().height-0.2])

plt.colorbar(cax=cax,extend='both')        
plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")

plt.show()

vmin=0
vmax=0.05
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
ax.set_facecolor((0.6, 0.6, 0.6))

norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
plt.pcolormesh(lon, lat, burden_year_ccc, 
                    transform=ccrs.PlateCarree(),cmap=RdYlGn.reversed(), norm=norm,vmin=vmin, vmax=vmax,shading='auto')
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
title='Dust_ccc Burden g m-2'
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
plt.pcolormesh(lon, lat, burden_year_cc, 
                    transform=ccrs.PlateCarree(),cmap=RdYlGn.reversed(), norm=norm,vmin=vmin, vmax=vmax,shading='auto')
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
title='Dust_cc Burden g m-2'
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
plt.pcolormesh(lon, lat, burden_year_c, 
                    transform=ccrs.PlateCarree(),cmap=RdYlGn.reversed(), norm=norm,vmin=vmin, vmax=vmax,shading='auto')
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
title='Dust_c Burden g m-2'
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
plt.pcolormesh(lon, lat, burden_year_f, 
                    transform=ccrs.PlateCarree(),cmap=RdYlGn.reversed(), norm=norm,vmin=vmin, vmax=vmax,shading='auto')
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
title='Dust_f Burden g m-2'
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
plt.pcolormesh(lon, lat, burden_year_ff, 
                    transform=ccrs.PlateCarree(),cmap=RdYlGn.reversed(), norm=norm,vmin=vmin, vmax=vmax,shading='auto')
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
title='Dust_ff Burden g m-2'
plt.title(title)
cax = fig.add_axes([ax.get_position().x1+0.02,
          ax.get_position().y0+0.09,
          0.08,
          ax.get_position().height-0.2])

plt.colorbar(cax=cax,extend='both')        
plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")

plt.show()




vmin=0.2
vmax=0.4

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
ax.set_facecolor((0.6, 0.6, 0.6))

norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
plt.pcolormesh(lon, lat, MEC, 
                    transform=ccrs.PlateCarree(),cmap=viridis, norm=norm,vmin=vmin, vmax=vmax,shading='auto')
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
title='DUST MEC m2 g-1'
plt.title(title)
cax = fig.add_axes([ax.get_position().x1+0.02,
          ax.get_position().y0+0.09,
          0.08,
          ax.get_position().height-0.2])

plt.colorbar(cax=cax,extend='both')        
plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")

plt.show()


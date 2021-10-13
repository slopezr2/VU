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


aod_file=Dataset('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_dry_aod_2008.nc')
tau_file=Dataset('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_wet_tau_2008.nc')
aod_datamanager=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_dry_aod_2008.nc')
h_file=Dataset('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_h_2008.nc')
gridded_emis=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_emis_2008.nc')
gridded_drydepo=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_drydepo_2008.nc')
gridded_wetdepo=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_wetdepo_2008.nc')
gridded_dust_aod=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_dust_aod_2008.nc')


gridded_LE=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_polder-aod-563_2008.nc')
# gridded_LE.graph_map(variable='ys',biascorr=1,vmin=0,vmax=0.3,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',ini_period=86,end_period=266,date=False,save=True,title='LE AOD565 Apr-Sept')
# gridded_LE.graph_map(variable='yr',biascorr=1,vmin=0,vmax=0.3,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',ini_period=86,end_period=266,date=False,save=True,title='POLDER AOD565 Apr-Sept')
    
# gridded_LE.graph_map(variable='ys',biascorr=1,vmin=0,vmax=0.3,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='LE AOD565 Yearly')
# gridded_LE.graph_map(variable='yr',biascorr=1,vmin=0,vmax=0.3,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='POLDER AOD565 Yearly')

# gridded_emis.graph_map(variable='tpm10',biascorr=1,vmin=0,vmax=2e-10,n_levels_plot=20,level=0,cmap=oranges,type_graph='mean',date=False,save=True,title='Total Aerosol Emissions')
# gridded_drydepo.graph_map(variable='tpm10_dflux',biascorr=1,vmin=0,vmax=2e-10,n_levels_plot=20,level=0,cmap=oranges,type_graph='mean',date=False,save=True,title='Aerosol Dry Deposition')
# gridded_wetdepo.graph_map(variable='tpm10_wflux',biascorr=1,vmin=0,vmax=2e-10,n_levels_plot=20,level=0,cmap=oranges,type_graph='mean',date=False,save=True,title='Aerosol Wet Deposition')


emis_year=np.squeeze(np.nanmean(gridded_emis.nc.variables['tpm10'][:],axis=0))
depo_year=np.squeeze(np.nanmean(gridded_drydepo.nc.variables['tpm10_dflux'][:],axis=0))+np.squeeze(np.nanmean(gridded_wetdepo.nc.variables['tpm10_wflux'][:],axis=0))

#aod_year=aod_file.variables['aod_563nm'][91:273,:,:] #abri a septiembre
#aod_year=aod_file.variables['aod_563nm'][:,:,:] #Todo el ano
aod_year=tau_file.variables['tau_563nm'][:,0,:,:] #Todo el ano
aod_year=np.squeeze(np.nanmean(aod_year,axis=0))
#i=3
#ini_dia=86
#end_dia=117
i=0
burden_year=[]
#for month in ['04','05','06','07','08','09']: #Abril a septiembre
for month in ['01','02','03','04','05','06','07','08','09','10','11','12']: #Abril a septiembre
    i=i+1
    print(month)
    file_conc='/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_conc_2008'+month+'.nc'
    conc_file_month=Dataset(file_conc)
    conc_month=conc_file_month.variables['tpm10'][:]
    conc_month=np.squeeze(np.nanmean(conc_month,axis=0))
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
    # gridded_LE.graph_map(variable='ys',biascorr=1,vmin=0,vmax=0.3,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',ini_period=ini_dia,end_period=end_dia,date=False,save=True,title='LE AOD565 2008-'+month)
    # gridded_LE.graph_map(variable='yr',biascorr=1,vmin=0,vmax=0.3,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',ini_period=ini_dia,end_period=end_dia,date=False,save=True,title='POLDER AOD565 2008-'+month)
    # aod_datamanager.graph_map(variable='aod_563nm',biascorr=1,vmin=0,vmax=0.3,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',ini_period=ini_dia*24,end_period=end_dia*24,date=False,save=True,title='LE output AOD565 2008-'+month)    
    h_month=h_file.variables['h'][ini_dia*24:end_dia*24,:,:,:]
    h_month=np.squeeze(np.nanmean(h_month,axis=0))
    burden=0
    for j in range(1):
        if j>0:
            burden=burden+conc_month[j,:,:]*(h_month[j,:,:]-h_month[j-1,:,:])
        else:
            burden=burden+conc_month[j,:,:]*(h_month[j,:,:])
    burden_year.append(burden)        
burden_year=np.array(burden_year)
burden_time=burden_year.copy()
burden_year=np.nanmean(burden_year,axis=0)*1e3 #En gramos
MEC=aod_year/burden_year





lon=aod_file.variables['lon'][:]
lat=aod_file.variables['lat'][:]
vmin=0
vmax=0.3
n_levels_plot=20
cmap=custom_ramp
biascorr=1


fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
ax.set_facecolor((0.6, 0.6, 0.6))

norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
plt.pcolormesh(lon, lat, aod_year, 
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
title='Surface AOD 563 nm 2008'
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
title='Surface Burden g m-2 Apr-Sept'
plt.title(title)
cax = fig.add_axes([ax.get_position().x1+0.02,
         ax.get_position().y0+0.09,
         0.08,
         ax.get_position().height-0.2])

plt.colorbar(cax=cax,extend='both')        
plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")

plt.show()

vmin=0
vmax=12
n_levels_plot=20
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
title='Surface Ambient MEC m2 g-1 2008'
plt.title(title)
cax = fig.add_axes([ax.get_position().x1+0.02,
         ax.get_position().y0+0.09,
         0.08,
         ax.get_position().height-0.2])

plt.colorbar(cax=cax,extend='both')        
plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")

plt.show()


# for i in range(23):
#     aod_datamanager.graph_map(variable='aod_563nm',biascorr=1,vmin=0,vmax=0.3,n_levels_plot=20,time=i,cmap=custom_ramp,type_graph='mean',ini_period=i,step_period=24,end_period=31*24,date=False,save=False,title='LE output AOD565 2008-'+month+' mean at '+str(i)+':00')    
    
# aod_datamanager.graph_map(variable='aod_563nm',biascorr=1,vmin=0,vmax=0.3,n_levels_plot=20,time=i,cmap=custom_ramp,type_graph='mean',ini_period=0,step_period=1,end_period=31*24,date=False,save=False,title='LE output AOD565 2008-01')    
    


# burden_emis=(burden_year/(1e3*emis_year))/(3600*24)
# burden_depo=(burden_year/(1e3*depo_year))/(3600*24)


# vmin=0
# vmax=7
# n_levels_plot=20
# fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
# levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
# ax.set_facecolor((0.6, 0.6, 0.6))

# norm = BoundaryNorm(levels, ncolors=viridis.N, clip=True)
# plt.pcolormesh(lon, lat, burden_emis, 
#                    transform=ccrs.PlateCarree(),cmap=viridis, norm=norm,vmin=vmin, vmax=vmax,shading='auto')
# european_borders=cfeature.NaturalEarthFeature(
#          category='cultural',
#          name='admin_0_countries',
#          scale='50m',
#          facecolor='none')


# coastlines=cfeature.NaturalEarthFeature(
#                  category='physical',
#                  name='coastline',
#                  scale='50m',
#                  facecolor='none')
# ax.add_feature(european_borders,edgecolor='black',linewidth=0.2)
# ax.add_feature(coastlines,edgecolor='black',linewidth=1)
# plt.xlim(lon[0],lon[-1])
# plt.ylim(lat[0],lat[-1])
# title='Lifetime burden emis [days]'
# plt.title(title)
# cax = fig.add_axes([ax.get_position().x1+0.02,
#          ax.get_position().y0+0.09,
#          0.08,
#          ax.get_position().height-0.2])

# plt.colorbar(cax=cax,extend='both')        
# plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")

# plt.show()

# fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
# levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
# ax.set_facecolor((0.6, 0.6, 0.6))

# norm = BoundaryNorm(levels, ncolors=viridis.N, clip=True)
# plt.pcolormesh(lon, lat, burden_depo, 
#                    transform=ccrs.PlateCarree(),cmap=viridis, norm=norm,vmin=vmin, vmax=vmax,shading='auto')
# european_borders=cfeature.NaturalEarthFeature(
#          category='cultural',
#          name='admin_0_countries',
#          scale='50m',
#          facecolor='none')


# coastlines=cfeature.NaturalEarthFeature(
#                  category='physical',
#                  name='coastline',
#                  scale='50m',
#                  facecolor='none')
# ax.add_feature(european_borders,edgecolor='black',linewidth=0.2)
# ax.add_feature(coastlines,edgecolor='black',linewidth=1)
# plt.xlim(lon[0],lon[-1])
# plt.ylim(lat[0],lat[-1])
# title='Lifetime burden depo [days]'
# plt.title(title)
# cax = fig.add_axes([ax.get_position().x1+0.02,
#          ax.get_position().y0+0.09,
#          0.08,
#          ax.get_position().height-0.2])

# plt.colorbar(cax=cax,extend='both')        
# plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")

# plt.show()
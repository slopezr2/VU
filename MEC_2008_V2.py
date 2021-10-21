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
    
    plt.colorbar(cax=cax,extend=extend)
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

aod_dry_file=Dataset('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_dry_aod_2008.nc')
tau_dry_file=Dataset('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_dry_tau_2008.nc')
aod_dry_datamanager=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_dry_aod_2008.nc')

aod_wet_file=Dataset('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_aod_2008.nc')
tau_wet_file=Dataset('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_wet_tau_2008.nc')
aod_wet_datamanager=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_aod_2008.nc')



h_file=Dataset('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_h_2008.nc')
gridded_emis=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_emis_2008.nc')
gridded_drydepo=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_drydepo_2008.nc')
gridded_wetdepo=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_wetdepo_2008.nc')
gridded_dust_aod=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_dust_aod_2008.nc')


gridded_LE=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_polder-aod-563_2008.nc')

gridded_LE.graph_map(variable='ys',biascorr=1,vmin=0,vmax=0.3,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',ini_period=86,end_period=266,date=False,save=True,title='LE AOD565 Apr-Sept',orientation='horizontal')

# gridded_LE.graph_map(variable='yr',biascorr=1,vmin=0,vmax=0.3,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',ini_period=86,end_period=266,date=False,save=True,title='POLDER AOD565 Apr-Sept')
    
# gridded_LE.graph_map(variable='ys',biascorr=1,vmin=0,vmax=0.3,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='LE AOD565 Yearly')
# gridded_LE.graph_map(variable='yr',biascorr=1,vmin=0,vmax=0.3,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='POLDER AOD565 Yearly')

# gridded_emis.graph_map(variable='tpm10',biascorr=1,vmin=0,vmax=2e-10,n_levels_plot=20,level=0,cmap=oranges,type_graph='mean',date=False,save=True,title='Total Aerosol Emissions')
# gridded_drydepo.graph_map(variable='tpm10_dflux',biascorr=1,vmin=0,vmax=2e-10,n_levels_plot=20,level=0,cmap=oranges,type_graph='mean',date=False,save=True,title='Aerosol Dry Deposition')
# gridded_wetdepo.graph_map(variable='tpm10_wflux',biascorr=1,vmin=0,vmax=2e-10,n_levels_plot=20,level=0,cmap=oranges,type_graph='mean',date=False,save=True,title='Aerosol Wet Deposition')


emis_year=np.squeeze(np.nanmean(gridded_emis.nc.variables['tpm10'][:],axis=0))
depo_year=np.squeeze(np.nanmean(gridded_drydepo.nc.variables['tpm10_dflux'][:],axis=0))+np.squeeze(np.nanmean(gridded_wetdepo.nc.variables['tpm10_wflux'][:],axis=0))


aod_dry_year=aod_dry_file.variables['aod_563nm'][:,:,:] #Todo el ano
aod_wet_year=aod_wet_file.variables['aod_563nm'][:,:,:] #Todo el ano


tau_dry_year=tau_dry_file.variables['tau_563nm'][:,0,:,:] #Todo el ano
tau_wet_year=tau_wet_file.variables['tau_563nm'][:,0,:,:] #Todo el ano


i=0
burden_year=[]
burden_surface_year=[]
for month in ['01','02','03','04','05','06','07','08','09','10','11','12']: #Abril a septiembre
    i=i+1
    print(month)
    file_conc='/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_conc_2008'+month+'.nc'
    conc_file_month=Dataset(file_conc)
    conc_month=conc_file_month.variables['tpm10'][:]
    
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
    
    if end_dia==-1:
        h_month=h_file.variables['h'][ini_dia*24:,:,:,:]
    else:
        h_month=h_file.variables['h'][ini_dia*24:end_dia*24,:,:,:]
    
    burden=0
    
    for j in range(conc_month.shape[1]):
        if j>0:
            burden=burden+conc_month[:,j,:,:]*(h_month[:,j,:,:]-h_month[:,j-1,:,:])
        else:
            burden=burden+conc_month[:,j,:,:]*(h_month[:,j,:,:])
            burden_surface=conc_month[:,j,:,:]*(h_month[:,j,:,:])
    
    if i==1:
        burden_year=np.array(burden)
        burden_surface_year=np.array(burden_surface)
    else:
        burden_year=np.append(burden_year,burden,axis=0)
        burden_surface_year=np.append(burden_surface_year,burden_surface,axis=0)
burden_year=burden_year*1e3 #En gramos
burden_surface_year=burden_surface_year*1e3




MEC_dry=aod_dry_year/burden_year
MEC_wet=aod_wet_year/burden_year

MEC_dry_surface=tau_dry_year/burden_surface_year
MEC_wet_surface=tau_wet_year/burden_surface_year

lon=aod_dry_file.variables['lon'][:]
lat=aod_dry_file.variables['lat'][:]

#==MEC Dry===
graph_map(lon=lon,lat=lat,variable=MEC_dry.mean(axis=0),vmin=0,vmax=3,n_levels_plot=20,cmap=viridis,title='MEC dry',save=True,extend='max',grid=True)

#==MEC ambient===
graph_map(lon=lon,lat=lat,variable=MEC_wet.mean(axis=0),vmin=0,vmax=10,n_levels_plot=20,cmap=viridis,title='MEC ambient',save=True,extend='max',grid=True)

#==MEC Surface Dry===
graph_map(lon=lon,lat=lat,variable=MEC_dry_surface.mean(axis=0),vmin=0,vmax=3,n_levels_plot=20,cmap=viridis,title='MEC dry surface',save=True,extend='max',grid=True)

#==MEC Surface ambient===
graph_map(lon=lon,lat=lat,variable=MEC_wet_surface.mean(axis=0),vmin=0,vmax=10,n_levels_plot=20,cmap=viridis,title='MEC ambient surface',save=True,extend='max',grid=True)

#==MEC Surface ambient/dry===
graph_map(lon=lon,lat=lat,variable=MEC_wet_surface.mean(axis=0)/MEC_dry_surface.mean(axis=0),vmin=0,vmax=5,n_levels_plot=20,cmap=bwr,title='MEC relation ambient dry',save=True,extend='max',vcenter=1,norm2=True,grid=True)


#===Seasonal Surface Dry===
graph_map(lon=lon,lat=lat,variable=MEC_dry_surface[np.r_[335*24:366*24-1,1:59*24],:,:].mean(axis=0),vmin=0,vmax=3,n_levels_plot=20,cmap=viridis,title='MEC dry surface DJF',save=True,extend='max',vcenter=1,norm2=False,grid=True)
graph_map(lon=lon,lat=lat,variable=MEC_dry_surface[59*24:152*24,:,:].mean(axis=0),vmin=0,vmax=3,n_levels_plot=20,cmap=viridis,title='MEC dry surface MAM',save=True,extend='max',vcenter=1,norm2=False,grid=True)
graph_map(lon=lon,lat=lat,variable=MEC_dry_surface[152*24:244*24,:,:].mean(axis=0),vmin=0,vmax=3,n_levels_plot=20,cmap=viridis,title='MEC dry surface JJA',save=True,extend='max',vcenter=1,norm2=False,grid=True)
graph_map(lon=lon,lat=lat,variable=MEC_dry_surface[244*24:335*24,:,:].mean(axis=0),vmin=0,vmax=3,n_levels_plot=20,cmap=viridis,title='MEC dry surface SON',save=True,extend='max',vcenter=1,norm2=False,grid=True)

#===Seasonal Surface wet===
graph_map(lon=lon,lat=lat,variable=MEC_wet_surface[np.r_[335*24:366*24-1,1:59*24],:,:].mean(axis=0),vmin=0,vmax=10,n_levels_plot=20,cmap=viridis,title='MEC ambient surface DJF',save=True,extend='max',vcenter=1,norm2=False,grid=True)
graph_map(lon=lon,lat=lat,variable=MEC_wet_surface[59*24:152*24,:,:].mean(axis=0),vmin=0,vmax=10,n_levels_plot=20,cmap=viridis,title='MEC ambient surface MAM',save=True,extend='max',vcenter=1,norm2=False,grid=True)
graph_map(lon=lon,lat=lat,variable=MEC_wet_surface[152*24:244*24,:,:].mean(axis=0),vmin=0,vmax=10,n_levels_plot=20,cmap=viridis,title='MEC ambient surface JJA',save=True,extend='max',vcenter=1,norm2=False,grid=True)
graph_map(lon=lon,lat=lat,variable=MEC_wet_surface[244*24:335*24,:,:].mean(axis=0),vmin=0,vmax=10,n_levels_plot=20,cmap=viridis,title='MEC ambient surface SON',save=True,extend='max',vcenter=1,norm2=False,grid=True)

#====Seasonal Surface relationship===
graph_map(lon=lon,lat=lat,variable=MEC_wet_surface[np.r_[335*24:366*24-1,1:59*24],:,:].mean(axis=0)/MEC_dry_surface[np.r_[335*24:366*24-1,1:59*24],:,:].mean(axis=0),vmin=0,vmax=5,n_levels_plot=20,cmap=bwr,title='MEC relation ambient dry DJF',save=True,extend='max',vcenter=1,norm2=True,grid=True)
graph_map(lon=lon,lat=lat,variable=MEC_wet_surface[59*24:152*24,:,:].mean(axis=0)/MEC_dry_surface[59*24:152*24,:,:].mean(axis=0),vmin=0,vmax=5,n_levels_plot=20,cmap=bwr,title='MEC relation ambient dry MAM',save=True,extend='max',vcenter=1,norm2=True,grid=True)
graph_map(lon=lon,lat=lat,variable=MEC_wet_surface[152*24:244*24,:,:].mean(axis=0)/MEC_dry_surface[152*24:244*24,:,:].mean(axis=0),vmin=0,vmax=5,n_levels_plot=20,cmap=bwr,title='MEC relation ambient dry JJA',save=True,extend='max',vcenter=1,norm2=True,grid=True)
graph_map(lon=lon,lat=lat,variable=MEC_wet_surface[244*24:335*24,:,:].mean(axis=0)/MEC_dry_surface[244*24:335*24,:,:].mean(axis=0),vmin=0,vmax=5,n_levels_plot=20,cmap=bwr,title='MEC relation ambient dry SON',save=True,extend='max',vcenter=1,norm2=True,grid=True)

lat_cabauw=65
lon_cabauw=38
cabauw=pd.DataFrame()
cabauw['MEC_dry']=MEC_dry_surface[:,lat_cabauw,lon_cabauw]
cabauw['MEC_wet']=MEC_wet_surface[:,lat_cabauw,lon_cabauw]
cabauw.index=pd.date_range(start='1/1/2008',end='31/12/2008 23:00:00',freq='h')

#==Time series monthly average===
fig, ax = plt.subplots()
cabauw[['MEC_dry','MEC_wet']].resample('M').mean().plot(ax=ax)
ax.legend(['MEC surface dry','MEC surface ambient'])
plt.savefig('./Figures/MEC_temporal_series.png',format='png', dpi=1000,bbox_inches = "tight")

bins=np.arange(1,23)
bins_label=['[0,1]','(1,2]','(2,3]','(3,4]','(4,5]','(5,6]','(6,7]','(7,8]','(8,9]','(9,10]','(10,11]','(11,12]','(12,13]','(13,14]','(14,15]','(15,16]','(16,17]','(17,18]','(18,19]','(19,20]','(21,22]','(22,23]']
fig, ax = plt.subplots()
cabauw['MEC_dry'].resample('D').mean().hist(bins=bins,rwidth = 0.7 , ax=ax)
cabauw['MEC_wet'].resample('D').mean().hist(bins=bins,rwidth = 0.7 ,ax=ax,alpha=0.7)
bins_labels(bins,bins_label, fontsize=10)
ax.legend(['MEC surface dry','MEC surface ambient'])
plt.xticks(rotation=90)
ax.grid(False,axis='x')
plt.savefig('./Figures/MEC_histogram.png',format='png', dpi=1000,bbox_inches = "tight")

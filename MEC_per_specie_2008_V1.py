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

def graph_map(lon=0,lat=0,vmin=0,vmax=0.3,vcenter=0,n_levels_plot=20,cmap=1,variable=1,title='AOD',save=False,extend='max',norm2=False,grid=False):
    

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
    cbar.formatter.set_powerlimits((0, 1))
    cbar.update_ticks()

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


LE_ec=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_ec_aod_2008.nc')
LE_na=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_na_aod_2008.nc')
LE_nh4a=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_nh4a_aod_2008.nc')
LE_no3a=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_no3a_aod_2008.nc')
LE_pom=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_pom_aod_2008.nc')
LE_ppm=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_ppm_aod_2008.nc')
LE_so4a=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_so4a_aod_2008.nc')
LE_dust=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_dust_aod_2008.nc')
LE_aod=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_aod_2008.nc')




aod_ec_year=LE_ec.nc.variables['aod_563nm'][:,:,:] #Todo el ano
aod_na_year=LE_na.nc.variables['aod_563nm'][:,:,:] #Todo el ano
aod_nh4a_year=LE_nh4a.nc.variables['aod_563nm'][:,:,:] #Todo el ano
aod_no3a_year=LE_no3a.nc.variables['aod_563nm'][:,:,:] #Todo el ano
aod_pom_year=LE_pom.nc.variables['aod_563nm'][:,:,:] #Todo el ano
aod_ppm_year=LE_ppm.nc.variables['aod_563nm'][:,:,:] #Todo el ano
aod_so4a_year=LE_so4a.nc.variables['aod_563nm'][:,:,:] #Todo el ano
aod_dust_year=LE_dust.nc.variables['aod_563nm'][:,:,:] #Todo el ano
aod_total_year=LE_aod.nc.variables['aod_563nm'][:,:,:] #Todo el ano
aod_sum_year=aod_ec_year+aod_na_year+aod_nh4a_year+aod_no3a_year+aod_pom_year+aod_ppm_year+aod_so4a_year+aod_dust_year

h_file=Dataset('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_h_2008.nc')


i=0
burden_year=[]
burden_surface_year=[]
for month in ['01','02','03','04','05','06','07','08','09','10','11','12']: #Abril a septiembre
    i=i+1
    print(month)
    file_conc='/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_conc_2008'+month+'.nc'
    conc_file_month=Dataset(file_conc)
    conc_ec_month=(conc_file_month.variables['ec_f'][:]+conc_file_month.variables['ec_c'][:])*1e3
    conc_na_month=(conc_file_month.variables['na_ff'][:]+conc_file_month.variables['na_f'][:]+conc_file_month.variables['na_c'][:]+conc_file_month.variables['na_cc'][:]+conc_file_month.variables['na_ccc'][:])*1e3*3.26
    conc_nh4a_month=(conc_file_month.variables['nh4a_f'][:])*1e3
    conc_no3a_month=(conc_file_month.variables['no3a_f'][:]+conc_file_month.variables['no3a_c'][:])*1e3
    conc_pom_month=(conc_file_month.variables['pom_f'][:]+conc_file_month.variables['pom_c'][:])*1e3
    conc_ppm_month=(conc_file_month.variables['ppm_f'][:]+conc_file_month.variables['ppm_c'][:])*1e3
    conc_so4a_month=(conc_file_month.variables['so4a_f'][:]+conc_file_month.variables['so4a_c'][:])*1e3
    conc_dust_month=(conc_file_month.variables['dust_ff'][:]+conc_file_month.variables['dust_f'][:]+conc_file_month.variables['dust_c'][:]+conc_file_month.variables['dust_cc'][:]+conc_file_month.variables['dust_ccc'][:])*1e3
    
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
    
    if end_dia==-1:
        h_month=h_file.variables['h'][ini_dia*24:,:,:,:]
    else:
        h_month=h_file.variables['h'][ini_dia*24:end_dia*24,:,:,:]
    
    burden_ec=0
    burden_na=0
    burden_nh4a=0
    burden_no3a=0
    burden_pom=0
    burden_ppm=0
    burden_so4a=0
    burden_dust=0
    
    for j in range(conc_ec_month.shape[1]):
        if j>0:
            burden_ec=burden_ec+conc_ec_month[:,j,:,:]*(h_month[:,j,:,:]-h_month[:,j-1,:,:])
            burden_na=burden_na+conc_na_month[:,j,:,:]*(h_month[:,j,:,:]-h_month[:,j-1,:,:])
            burden_nh4a=burden_nh4a+conc_nh4a_month[:,j,:,:]*(h_month[:,j,:,:]-h_month[:,j-1,:,:])
            burden_no3a=burden_no3a+conc_no3a_month[:,j,:,:]*(h_month[:,j,:,:]-h_month[:,j-1,:,:])
            burden_pom=burden_pom+conc_pom_month[:,j,:,:]*(h_month[:,j,:,:]-h_month[:,j-1,:,:])
            burden_ppm=burden_ppm+conc_ppm_month[:,j,:,:]*(h_month[:,j,:,:]-h_month[:,j-1,:,:])
            burden_so4a=burden_so4a+conc_so4a_month[:,j,:,:]*(h_month[:,j,:,:]-h_month[:,j-1,:,:])
            burden_dust=burden_dust+conc_dust_month[:,j,:,:]*(h_month[:,j,:,:]-h_month[:,j-1,:,:])
        else:
            burden_ec=burden_ec+conc_ec_month[:,j,:,:]*(h_month[:,j,:,:])
            burden_na=burden_na+conc_na_month[:,j,:,:]*(h_month[:,j,:,:])
            burden_nh4a=burden_nh4a+conc_nh4a_month[:,j,:,:]*(h_month[:,j,:,:])
            burden_no3a=burden_no3a+conc_no3a_month[:,j,:,:]*(h_month[:,j,:,:])
            burden_pom=burden_pom+conc_pom_month[:,j,:,:]*(h_month[:,j,:,:])
            burden_ppm=burden_ppm+conc_ppm_month[:,j,:,:]*(h_month[:,j,:,:])
            burden_so4a=burden_so4a+conc_so4a_month[:,j,:,:]*(h_month[:,j,:,:])
            burden_dust=burden_dust+conc_dust_month[:,j,:,:]*(h_month[:,j,:,:])
            
    
    if i==1:
        burden_ec_year=np.array(burden_ec)
        burden_na_year=np.array(burden_na)
        burden_nh4a_year=np.array(burden_nh4a)
        burden_no3a_year=np.array(burden_no3a)
        burden_pom_year=np.array(burden_pom)
        burden_ppm_year=np.array(burden_ppm)
        burden_so4a_year=np.array(burden_so4a)
        burden_dust_year=np.array(burden_dust)
        
    else:
        burden_ec_year=np.append(burden_ec_year,burden_ec,axis=0)
        burden_na_year=np.append(burden_na_year,burden_na,axis=0)
        burden_nh4a_year=np.append(burden_nh4a_year,burden_nh4a,axis=0)
        burden_no3a_year=np.append(burden_no3a_year,burden_no3a,axis=0)
        burden_pom_year=np.append(burden_pom_year,burden_pom,axis=0)
        burden_ppm_year=np.append(burden_ppm_year,burden_ppm,axis=0)
        burden_so4a_year=np.append(burden_so4a_year,burden_so4a,axis=0)
        burden_dust_year=np.append(burden_dust_year,burden_dust,axis=0)
        


burden_total_year=burden_ec_year+burden_na_year+burden_nh4a_year+burden_no3a_year+burden_pom_year+burden_ppm_year+burden_so4a_year+burden_dust_year

MEC_ec=aod_ec_year/burden_ec_year
MEC_na=aod_na_year/burden_na_year
MEC_nh4a=aod_nh4a_year/burden_nh4a_year
MEC_no3a=aod_no3a_year/burden_no3a_year
MEC_pom=aod_pom_year/burden_pom_year
MEC_ppm=aod_ppm_year/burden_ppm_year
MEC_so4a=aod_so4a_year/burden_so4a_year
MEC_dust=aod_dust_year/burden_dust_year
MEC_total=aod_total_year/burden_total_year


MEC_ec[MEC_ec>20]=np.nan
MEC_na[MEC_na>20]=np.nan
MEC_nh4a[MEC_nh4a>20]=np.nan
MEC_no3a[MEC_no3a>20]=np.nan
MEC_pom[MEC_pom>20]=np.nan
MEC_ppm[MEC_ppm>20]=np.nan
MEC_dust[MEC_dust>20]=np.nan
MEC_total[MEC_total>20]=np.nan

lon=LE_ec.nc.variables['lon'][:]
lat=LE_ec.nc.variables['lat'][:]

#==aod per specie===
LE_ec.graph_map(variable='aod_563nm',vmin=0,vmax=0.1,n_levels_plot=20,cmap=custom_ramp,title='AOD 563nm EC 2008',save=True,extend='max',grid=True,date=False,type_graph='mean')
LE_na.graph_map(variable='aod_563nm',vmin=0,vmax=0.1,n_levels_plot=20,cmap=custom_ramp,title='AOD 563nm Na 2008',save=True,extend='max',grid=True,date=False,type_graph='mean')
LE_nh4a.graph_map(variable='aod_563nm',vmin=0,vmax=0.1,n_levels_plot=20,cmap=custom_ramp,title='AOD 563nm NH4 2008',save=True,extend='max',grid=True,date=False,type_graph='mean')
LE_no3a.graph_map(variable='aod_563nm',vmin=0,vmax=0.1,n_levels_plot=20,cmap=custom_ramp,title='AOD 563nm NO3 2008',save=True,extend='max',grid=True,date=False,type_graph='mean')
LE_pom.graph_map(variable='aod_563nm',vmin=0,vmax=0.1,n_levels_plot=20,cmap=custom_ramp,title='AOD 563nm pom 2008',save=True,extend='max',grid=True,date=False,type_graph='mean')
LE_ppm.graph_map(variable='aod_563nm',vmin=0,vmax=0.1,n_levels_plot=20,cmap=custom_ramp,title='AOD 563nm ppm 2008',save=True,extend='max',grid=True,date=False,type_graph='mean')
LE_so4a.graph_map(variable='aod_563nm',vmin=0,vmax=0.1,n_levels_plot=20,cmap=custom_ramp,title='AOD 563nm SO4 2008',save=True,extend='max',grid=True,date=False,type_graph='mean')
LE_dust.graph_map(variable='aod_563nm',vmin=0,vmax=0.1,n_levels_plot=20,cmap=custom_ramp,title='AOD 563nm Dust 2008',save=True,extend='max',grid=True,date=False,type_graph='mean')
LE_aod.graph_map(variable='aod_563nm',vmin=0,vmax=0.3,n_levels_plot=20,cmap=custom_ramp,title='AOD 563nm 2008',save=True,extend='max',grid=True,date=False,type_graph='mean')

graph_map(lon=lon,lat=lat,variable=aod_sum_year.mean(axis=0),vmin=0,vmax=0.3,n_levels_plot=20,cmap=custom_ramp,title='AOD 563nm sum 2008',save=True,extend='max',grid=True)

#==MEC ambient===
graph_map(lon=lon,lat=lat,variable=np.nanmean(MEC_ec,axis=0),vmin=0,vmax=10,n_levels_plot=20,cmap=viridis,title='MEC EC ambient',save=True,extend='max',grid=True)
graph_map(lon=lon,lat=lat,variable=np.nanmean(MEC_na,axis=0),vmin=0,vmax=10,n_levels_plot=20,cmap=viridis,title='MEC Na ambient',save=True,extend='max',grid=True)
graph_map(lon=lon,lat=lat,variable=np.nanmean(MEC_nh4a,axis=0),vmin=0,vmax=10,n_levels_plot=20,cmap=viridis,title='MEC NH4 ambient',save=True,extend='max',grid=True)
graph_map(lon=lon,lat=lat,variable=np.nanmean(MEC_no3a,axis=0),vmin=0,vmax=10,n_levels_plot=20,cmap=viridis,title='MEC NO3 ambient',save=True,extend='max',grid=True)
graph_map(lon=lon,lat=lat,variable=np.nanmean(MEC_pom,axis=0),vmin=0,vmax=10,n_levels_plot=20,cmap=viridis,title='MEC pom ambient',save=True,extend='max',grid=True)
graph_map(lon=lon,lat=lat,variable=np.nanmean(MEC_ppm,axis=0),vmin=0,vmax=10,n_levels_plot=20,cmap=viridis,title='MEC ppm ambient',save=True,extend='max',grid=True)
graph_map(lon=lon,lat=lat,variable=np.nanmean(MEC_so4a,axis=0),vmin=0,vmax=10,n_levels_plot=20,cmap=viridis,title='MEC SO4 ambient',save=True,extend='max',grid=True)
graph_map(lon=lon,lat=lat,variable=np.nanmean(MEC_dust,axis=0),vmin=0,vmax=10,n_levels_plot=20,cmap=viridis,title='MEC Dust ambient',save=True,extend='max',grid=True)

graph_map(lon=lon,lat=lat,variable=np.nanmean(MEC_total,axis=0),vmin=0,vmax=10,n_levels_plot=20,cmap=viridis,title='MEC All ambient',save=False,extend='max',grid=True)


lat_cabauw=65
lon_cabauw=38
MEC_series=pd.DataFrame()
MEC_series['MEC_ec_cabauw']=MEC_ec[:,lat_cabauw,lon_cabauw]
MEC_series['MEC_na_cabauw']=MEC_na[:,lat_cabauw,lon_cabauw]
MEC_series['MEC_nh4a_cabauw']=MEC_nh4a[:,lat_cabauw,lon_cabauw]
MEC_series['MEC_no3a_cabauw']=MEC_no3a[:,lat_cabauw,lon_cabauw]
MEC_series['MEC_pom_cabauw']=MEC_pom[:,lat_cabauw,lon_cabauw]
MEC_series['MEC_ppm_cabauw']=MEC_ppm[:,lat_cabauw,lon_cabauw]
MEC_series['MEC_so4a_cabauw']=MEC_so4a[:,lat_cabauw,lon_cabauw]
MEC_series['MEC_dust_cabauw']=MEC_dust[:,lat_cabauw,lon_cabauw]
MEC_series['MEC_total_cabauw']=MEC_total[:,lat_cabauw,lon_cabauw]

MEC_series['MEC_ec_domain']=MEC_ec.mean(axis=1).mean(axis=1)
MEC_series['MEC_na_domain']=MEC_na.mean(axis=1).mean(axis=1)
MEC_series['MEC_nh4a_domain']=MEC_nh4a.mean(axis=1).mean(axis=1)
MEC_series['MEC_no3a_domain']=MEC_no3a.mean(axis=1).mean(axis=1)
MEC_series['MEC_pom_domain']=MEC_pom.mean(axis=1).mean(axis=1)
MEC_series['MEC_ppm_domain']=MEC_ppm.mean(axis=1).mean(axis=1)
MEC_series['MEC_so4a_domain']=MEC_so4a.mean(axis=1).mean(axis=1)
MEC_series['MEC_dust_domain']=MEC_dust.mean(axis=1).mean(axis=1)
MEC_series['MEC_total_domain']=MEC_total.mean(axis=1).mean(axis=1)


MEC_series['aod_ec_domain']=aod_ec_year.mean(axis=1).mean(axis=1)
MEC_series['aod_na_domain']=aod_na_year.mean(axis=1).mean(axis=1)
MEC_series['aod_nh4a_domain']=aod_nh4a_year.mean(axis=1).mean(axis=1)
MEC_series['aod_no3a_domain']=aod_no3a_year.mean(axis=1).mean(axis=1)
MEC_series['aod_pom_domain']=aod_pom_year.mean(axis=1).mean(axis=1)
MEC_series['aod_ppm_domain']=aod_ppm_year.mean(axis=1).mean(axis=1)
MEC_series['aod_so4a_domain']=aod_so4a_year.mean(axis=1).mean(axis=1)
MEC_series['aod_dust_domain']=aod_dust_year.mean(axis=1).mean(axis=1)
MEC_series['aod_total_domain']=aod_total_year.mean(axis=1).mean(axis=1)
MEC_series['aod_sum_domain']=aod_sum_year.mean(axis=1).mean(axis=1)


MEC_series['burden_ec_domain']=aod_ec_year.mean(axis=1).mean(axis=1)
MEC_series['burden_na_domain']=burden_na_year.mean(axis=1).mean(axis=1)
MEC_series['burden_nh4a_domain']=burden_nh4a_year.mean(axis=1).mean(axis=1)
MEC_series['burden_no3a_domain']=burden_no3a_year.mean(axis=1).mean(axis=1)
MEC_series['burden_pom_domain']=burden_pom_year.mean(axis=1).mean(axis=1)
MEC_series['burden_ppm_domain']=burden_ppm_year.mean(axis=1).mean(axis=1)
MEC_series['burden_so4a_domain']=burden_so4a_year.mean(axis=1).mean(axis=1)
MEC_series['burden_dust_domain']=burden_dust_year.mean(axis=1).mean(axis=1)
MEC_series['burden_total_domain']=burden_total_year.mean(axis=1).mean(axis=1)


MEC_series.index=pd.date_range(start='1/1/2008',end='31/12/2008 23:00:00',freq='h')

labels=['EC','pom',r'NH$_{4}$',r'NO$_{3}$','Na','ppm',r'SO$_{4}$','Dust','All','Sum']
# #==Time series MEC monthly average===
fig, ax = plt.subplots()
MEC_series[['MEC_ec_domain','MEC_pom_domain','MEC_nh4a_domain','MEC_no3a_domain','MEC_na_domain','MEC_ppm_domain','MEC_so4a_domain','MEC_dust_domain','MEC_total_domain',]].resample('M').mean().plot(ax=ax,linewidth=2,marker='s')
ax.legend(labels,loc='upper center', bbox_to_anchor=(0.5, 1.15),ncol=5, fancybox=True, shadow=True)
plt.rcParams['axes.titley'] = 1.13    # y is in axes-relative coordinates.
plt.title(r'Average ambient MEC value m$^{2}$ g$^{-1}$')
plt.savefig('./Figures/MEC_species_temporal_series.png',format='png', dpi=1000,bbox_inches = "tight")



# #==Time series MEC monthly average===
fig, ax = plt.subplots()
MEC_series[['MEC_ec_cabauw','MEC_pom_cabauw','MEC_nh4a_cabauw','MEC_no3a_cabauw','MEC_na_cabauw','MEC_ppm_cabauw','MEC_so4a_cabauw','MEC_dust_cabauw','MEC_total_cabauw',]].resample('M').mean().plot(ax=ax,linewidth=2,marker='s')
ax.legend(labels,loc='upper center', bbox_to_anchor=(0.5, 1.15),ncol=5, fancybox=True, shadow=True)
plt.rcParams['axes.titley'] = 1.13    # y is in axes-relative coordinates.
plt.title(r'Cabauw ambient MEC value m$^{2}$ g$^{-1}$')
plt.savefig('./Figures/MEC_species_temporal_series_cabauw.png',format='png', dpi=1000,bbox_inches = "tight")


# #==Time series aod monthly average===
fig, ax = plt.subplots()
MEC_series[['aod_ec_domain','aod_pom_domain','aod_nh4a_domain','aod_no3a_domain','aod_na_domain','aod_ppm_domain','aod_so4a_domain','aod_dust_domain','aod_total_domain','aod_sum_domain']].resample('M').mean().plot(ax=ax,linewidth=2,marker='s')
ax.legend(labels,loc='upper center', bbox_to_anchor=(0.5, 1.15),ncol=5, fancybox=True, shadow=True)
plt.rcParams['axes.titley'] = 1.13    # y is in axes-relative coordinates.
plt.title(r'Average ambient AOD 563nm')
plt.savefig('./Figures/aod_species_temporal_series.png',format='png', dpi=1000,bbox_inches = "tight")


# #==Time series burden monthly average===
fig, ax = plt.subplots()
MEC_series[['burden_ec_domain','burden_pom_domain','burden_nh4a_domain','burden_no3a_domain','burden_na_domain','burden_ppm_domain','burden_so4a_domain','burden_dust_domain','burden_total_domain',]].resample('M').mean().plot(ax=ax,linewidth=2,marker='s')
ax.legend(labels,loc='upper center', bbox_to_anchor=(0.5, 1.15),ncol=5, fancybox=True, shadow=True)
plt.rcParams['axes.titley'] = 1.13    # y is in axes-relative coordinates.
plt.title(r'Average burden g m$^{-2}$')
plt.savefig('./Figures/burden_species_temporal_series.png',format='png', dpi=1000,bbox_inches = "tight")


# #==Pie aod chart average===
fig, ax = plt.subplots()

explode=(0,0.5,0,0,0,0.5,0,0)
MEC_series[['aod_ec_domain','aod_pom_domain','aod_nh4a_domain','aod_no3a_domain','aod_na_domain',
            'aod_ppm_domain','aod_so4a_domain','aod_dust_domain']].mean().plot.pie(ax=ax,normalize=True,
       labels=labels, wedgeprops={'linewidth': 2.0, 'edgecolor': 'white'}, textprops={'size': 13},
       autopct='%.1f%%',pctdistance=0.6,labeldistance=1.1,explode=explode,counterclock=False,startangle=0)

plt.rcParams['axes.titley'] = 1                                                                                                                                                                    
plt.title('Ambient AOD 563nm composition',fontsize=12)
plt.ylabel('')
plt.savefig('./Figures/aod_species_composition.png',format='png', dpi=1000,bbox_inches = "tight")

# #==Pie burden chart average===
fig, ax = plt.subplots()
labels=['EC','pom',r'NH$_{4}$',r'NO$_{3}$','Na','ppm',r'SO$_{4}$','Dust','All']
explode=(0,0.5,0,0,0.2,0.5,0,0)
h=MEC_series[['burden_ec_domain','burden_pom_domain','burden_nh4a_domain','burden_no3a_domain','burden_na_domain',
            'burden_ppm_domain','burden_so4a_domain','burden_dust_domain']].mean().plot.pie(ax=ax,normalize=True,
       labels=labels, wedgeprops={'linewidth': 2.0, 'edgecolor': 'white'}, textprops={'size': 13},
       autopct='%.1f%%',pctdistance=0.6,labeldistance=1.1,explode=explode,counterclock=False,startangle=0)

plt.rcParams['axes.titley'] = 1                                                                                                                                                              
plt.title('Burden composition',fontsize=12)
plt.ylabel('')
plt.savefig('./Figures/burden_species_composition.png',format='png', dpi=1000,bbox_inches = "tight")


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
from datamanager import graph_map
from matplotlib import cm


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
jet=cm.get_cmap('jet', 256)
nipy=cm.get_cmap('nipy_spectral', 256)
sizes=['ec_c','ec_f','na_ccc','na_cc','na_c','na_ff','na_f','dust_ccc','dust_cc','dust_c','dust_ff','dust_f','no3a_c','no3a_f','ppm_c','ppm_f','pom_c','pom_f','so4a_c','so4a_f','nh4a_f']
#sizes=['ec_c','ec_f','na_ccc','na_cc','na_c','na_ff','na_f','dust_ccc','dust_cc','dust_c','dust_ff','dust_f','no3a_c','no3a_f','ppm_f','pom_c','pom_f','so4a_f','nh4a_f']
#sizes=['ppm_f','pom_c','pom_f','so4a_f','nh4a']
#sizes=['nh4a_f']
#sizes=['so4a_c']


h_file=Dataset('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_h_2008.nc')
lon=h_file.variables['lon']
lat=h_file.variables['lat']
burden_year={}
MEC={}
aod_year={}
nested_dict={}
nested_dict['aod_domain']={}
nested_dict['burden_domain']={}
nested_dict['MEC_dry_domain']={}
nested_dict['MEC_dry_cabauw']={}

lat_cabauw=65
lon_cabauw=38
for size in sizes:
    if size[0:4]=='nh4a':
        LE=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_dry_nh4_2008.nc')
    else:
        LE=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_'+size+'_dry_2008.nc')
    
    aod_year[size]=LE.nc.variables['aod_563nm'][:]
    for month in ['01','02','03','04','05','06','07','08','09','10','11','12']: #Abril a septiembre
        i=int(month)
        print(month)
        file_conc='/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_conc_2008'+month+'.nc'
        conc_file_month=DataManager_LE(file_conc)
        conc_month=conc_file_month.nc.variables[size][:]*1e3
        if size[0:2]=='na':
            conc_month=conc_month*3.26
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
        burden=0    
        for j in range(conc_month.shape[1]):
            if j>0:
                burden=burden+conc_month[:,j,:,:]*(h_month[:,j,:,:]-h_month[:,j-1,:,:]) 
            else:
                burden=burden+conc_month[:,j,:,:]*(h_month[:,j,:,:])
                
        if i==1:
            burden_year[size]=np.array(burden)
        else:
            burden_year[size]=np.append(burden_year[size],burden,axis=0)
       
    MEC[size]=aod_year[size]/burden_year[size]
    MEC[size][MEC[size]>20]=np.nan
    
    graph_map(lon=lon,lat=lat,variable=np.nanmean(MEC[size],axis=0),vmin=0,vmax=10,n_levels_plot=20,cmap=viridis,title='MEC '+size+' dry',save=True,date=False,extend='max',grid=True,dpi=300)
    graph_map(lon=lon,lat=lat,variable=np.nanmean(burden_year[size],axis=0),vmin=0,vmax=0.01,n_levels_plot=20,cmap=RdYlGn,title='Burden '+size,save=True,date=False,extend='max',grid=True,dpi=300)
    graph_map(lon=lon,lat=lat,variable=np.nanmean(aod_year[size],axis=0),vmin=0,vmax=0.01,n_levels_plot=20,cmap=custom_ramp,title='AOD 563nm '+size+' dry',save=True,date=False,extend='max',grid=True,dpi=300)
    
    nested_dict['MEC_dry_cabauw'][size]=MEC[size][:,lat_cabauw,lon_cabauw]
    nested_dict['MEC_dry_domain'][size]=MEC[size].mean(axis=1).mean(axis=1)
    nested_dict['aod_domain'][size]=aod_year[size].mean(axis=1).mean(axis=1)
    nested_dict['burden_domain'][size]=burden_year[size].mean(axis=1).mean(axis=1)
    
    

reformed_dict = {}
for outerKey, innerDict in nested_dict.items():
    for innerKey, values in innerDict.items():
        reformed_dict[(outerKey,
                       innerKey)] = values

MEC_series=pd.DataFrame(reformed_dict)
MEC_series.index=pd.date_range(start='1/1/2008',end='31/12/2008 23:00:00',freq='h')



fig, ax = plt.subplots()
MEC_series['MEC_dry_domain'].resample('M').mean().plot(ax=ax,linewidth=3,marker='s')
ax.legend(sizes,loc='upper center', bbox_to_anchor=(0.5, 1.2),ncol=5, fancybox=True, shadow=True)
ax.set_ylabel('MEC dry [m$^{-2}$/g]')
plt.savefig('./Figures/MEC_domain_dry_sizes_series.png',format='png', dpi=300,bbox_inches = "tight")


fig, ax = plt.subplots()
MEC_series['aod_domain'].resample('M').mean().plot(ax=ax,linewidth=3,marker='s')
ax.legend(sizes,loc='upper center', bbox_to_anchor=(0.5, 1.2),ncol=5, fancybox=True, shadow=True)
ax.set_ylabel('AOD 563nm')
plt.savefig('./Figures/AOD_dry_sizes_series.png',format='png', dpi=300,bbox_inches = "tight")


fig, ax = plt.subplots()
MEC_series['MEC_dry_cabauw'].resample('M').mean().plot(ax=ax,linewidth=3,marker='s')
ax.legend(sizes,loc='upper center', bbox_to_anchor=(0.5, 1.2),ncol=5, fancybox=True, shadow=True)
ax.set_ylabel('AOD 563nm')
plt.savefig('./Figures/MEC_cabauw_dry_sizes_series.png',format='png', dpi=300,bbox_inches = "tight")

fig, ax = plt.subplots()
MEC_series['burden_domain'].resample('M').mean().plot(ax=ax,linewidth=3,marker='s')
ax.legend(sizes,loc='upper center', bbox_to_anchor=(0.5, 1.2),ncol=5, fancybox=True, shadow=True)
ax.set_ylabel('Burden $\mu$g/m$^{2}$')
plt.savefig('./Figures/burden_sizes_series.png',format='png', dpi=300,bbox_inches = "tight")



#for size in sizes:


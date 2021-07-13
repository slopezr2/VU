import sys
#===Load user libraries====
sys.path.append('./libraries/')
from datamanager import DataManager_GRIDDED
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
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
from matplotlib import cm


viridis = cm.get_cmap('viridis', 256)
orange=cm.get_cmap('Oranges', 256)

#===Create Color bar====
def make_Ramp( ramp_colors ): 
    from colour import Color
    from matplotlib.colors import LinearSegmentedColormap

    color_ramp = LinearSegmentedColormap.from_list( 'my_list', [ Color( c1 ).rgb for c1 in ramp_colors ] )
    return color_ramp


custom_ramp = make_Ramp( ["#ffffff","#aec8e0","#9dafb2","#e2b85f","#f0ae46","#ea9842","#d36f2b","#a04617","#6c423d"] ) 


no3a_yearly=[]
na_yearly=[]
so4a_yearly=[]
nh4a_yearly=[]
ppm_yearly=[]
pom_yearly=[]
ec_yearly=[]
dust_yearly=[]
total_yearly=[]

std_total=[]

no3a_monthly=[]
na_monthly=[]
so4a_monthly=[]
nh4a_monthly=[]
ppm_monthly=[]
pom_monthly=[]
ec_monthly=[]
dust_monthly=[]
total_monthly=[]
for month in ['01','02','03','04','05','06','07','08','09','10','11','12']:
    print(month)
    file='/Users/santiago/Documents/LE_outputs/2008_POLDER_Complete/gridded/LE_2008'+month+'.nc'
    file_conc='/Users/santiago/Documents/LE_outputs/2008_Cocentration/LE_conc_2008'+month+'.nc'
    
    le_month=Dataset(file_conc)
    # conc_month={}
    
    # conc_month['no3a']=le_month.variables['no3a_c'][:]+le_month.variables['no3a_f'][:]
    # conc_month['no3a']=np.sum(conc_month['no3a'],axis=1)
    # no3a_yearly=no3a_yearly+list((np.mean(np.mean(conc_month['no3a'],axis=1),axis=1)))
    
    # no3a_monthly.append(np.mean(np.mean(np.mean(conc_month['no3a'],axis=1),axis=1)))
    
    # conc_month['na']=le_month.variables['na_ccc'][:]+le_month.variables['na_cc'][:]+le_month.variables['na_c'][:]+le_month.variables['na_f'][:]+le_month.variables['na_ff'][:]
    # conc_month['na']=np.sum(conc_month['na'],axis=1)
    # na_yearly=na_yearly+list((np.mean(np.mean(conc_month['na'],axis=1),axis=1)))
    # na_monthly.append(np.mean(np.mean(np.mean(conc_month['na'],axis=1),axis=1)))
    
    # conc_month['so4a']=le_month.variables['so4a_c'][:]+le_month.variables['so4a_f'][:]
    # conc_month['so4a']=np.sum(conc_month['so4a'],axis=1)
    # so4a_yearly=so4a_yearly+list((np.mean(np.mean(conc_month['so4a'],axis=1),axis=1)))
    # so4a_monthly.append(np.mean(np.mean(np.mean(conc_month['so4a'],axis=1),axis=1)))
    
    # conc_month['nh4a']=le_month.variables['nh4a_f'][:]
    # conc_month['nh4a']=np.sum(conc_month['nh4a'],axis=1)
    # nh4a_yearly=nh4a_yearly+list((np.mean(np.mean(conc_month['nh4a'],axis=1),axis=1)))
    # nh4a_monthly.append(np.mean(np.mean(np.mean(conc_month['nh4a'],axis=1),axis=1)))
    
    # conc_month['ppm']=le_month.variables['ppm_c'][:]+le_month.variables['ppm_f'][:]
    # conc_month['ppm']=np.sum(conc_month['ppm'],axis=1)
    # ppm_yearly=ppm_yearly+list((np.mean(np.mean(conc_month['ppm'],axis=1),axis=1)))
    # ppm_monthly.append(np.mean(np.mean(np.mean(conc_month['ppm'],axis=1),axis=1)))
    
    # conc_month['pom']=le_month.variables['pom_c'][:]+le_month.variables['pom_f'][:]
    # conc_month['pom']=np.sum(conc_month['pom'],axis=1)
    # pom_yearly=pom_yearly+list((np.mean(np.mean(conc_month['pom'],axis=1),axis=1)))
    # pom_monthly.append(np.mean(np.mean(np.mean(conc_month['pom'],axis=1),axis=1)))
    
    # conc_month['ec']=le_month.variables['ec_c'][:]+le_month.variables['ec_f'][:]
    # conc_month['ec']=np.sum(conc_month['ec'],axis=1)
    # ec_yearly=ec_yearly+list((np.mean(np.mean(conc_month['ec'],axis=1),axis=1)))
    # ec_monthly.append(np.mean(np.mean(np.mean(conc_month['ec'],axis=1),axis=1)))
    
    # conc_month['dust']=le_month.variables['dust_ccc'][:]+le_month.variables['dust_cc'][:]+le_month.variables['dust_c'][:]+le_month.variables['dust_f'][:]+le_month.variables['dust_ff'][:]
    # conc_month['dust']=np.sum(conc_month['dust'],axis=1)
    # dust_yearly=dust_yearly+list((np.mean(np.mean(conc_month['dust'],axis=1),axis=1)))
    # dust_monthly.append(np.mean(np.mean(np.mean(conc_month['dust'],axis=1),axis=1)))
    
    # conc_month['total']=conc_month['no3a']+conc_month['na']+conc_month['so4a']+conc_month['nh4a']+conc_month['ppm']+conc_month['pom']+conc_month['ec']+conc_month['dust']
    # total_yearly=total_yearly+list((np.mean(np.mean(conc_month['total'],axis=1),axis=1)))
    # std_total=std_total+list((np.std(np.std(conc_month['total'],axis=1),axis=1)))
    # total_monthly.append(np.mean(np.mean(np.mean(conc_month['total'],axis=1),axis=1)))
    
    
    
    # vmin=0
    
    # n_levels_plot=20
    # cmap=orange
    # biascorr=1
    # lon=le_month.variables['lon'][:]
    # lat=le_month.variables['lat'][:]
    # for specie in ['total']:#['no3a','na','so4a','nh4a','ppm','pom','ec','dust','total']:
    #     title='Average Total Aerosol 2008-'+month
    #     if specie=='total':
    #         vmax=200
    #     else:
    #         vmax=20
    #     fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    #     levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
    #     ax.set_facecolor((0.6, 0.6, 0.6))
        
    #     norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    #     plt.pcolormesh(lon, lat, np.mean(conc_month[specie]*1e9,axis=0), 
    #                         transform=ccrs.PlateCarree(),cmap=cmap, norm=norm,vmin=vmin, vmax=vmax,shading='auto')
    #     european_borders=cfeature.NaturalEarthFeature(
    #               category='cultural',
    #               name='admin_0_countries',
    #               scale='50m',
    #               facecolor='none')
        
        
    #     coastlines=cfeature.NaturalEarthFeature(
    #                       category='physical',
    #                       name='coastline',
    #                       scale='50m',
    #                       facecolor='none')
    #     ax.add_feature(european_borders,edgecolor='black',linewidth=0.2)
    #     ax.add_feature(coastlines,edgecolor='black',linewidth=1)
    #     plt.xlim(lon[0],lon[-1])
    #     plt.ylim(lat[0],lat[-1])
    #     plt.title(title)
                
        
    #     cax = fig.add_axes([ax.get_position().x1+0.02,
    #               ax.get_position().y0+0.09,
    #               0.08,
    #               ax.get_position().height-0.2])
        
    #     plt.colorbar(cax=cax,extend='both')
    #     plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")
    #     plt.show()
    
    
    gridded_LE=DataManager_GRIDDED(file)
    gridded_LE.graph_map(variable='ys',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='LE AOD565 2008-'+month)
    gridded_LE.graph_map(variable='yr',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='POLDER AOD565 2008-'+month)
    # gridded_LE.graph_map(variable='ys',biascorr=2.5,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=True,save=False,title='LE AOD565 Biascorr=2.5')



# no3a_yearly=np.array(no3a_yearly)*1e9
# na_yearly=np.array(na_yearly)*1e9
# so4a_yearly=np.array(so4a_yearly)*1e9
# nh4a_yearly=np.array(nh4a_yearly)*1e9
# ppm_yearly=np.array(ppm_yearly)*1e9
# pom_yearly=np.array(pom_yearly)*1e9
# ec_yearly=np.array(ec_yearly)*1e9
# dust_yearly=np.array(dust_yearly)*1e9
# total_yearly=np.array(total_yearly)*1e9


# no3a_monthly=np.array(no3a_monthly)*1e9
# na_monthly=np.array(na_monthly)*1e9
# so4a_monthly=np.array(so4a_monthly)*1e9
# nh4a_monthly=np.array(nh4a_monthly)*1e9
# ppm_monthly=np.array(ppm_monthly)*1e9
# pom_monthly=np.array(pom_monthly)*1e9
# ec_monthly=np.array(ec_monthly)*1e9
# dust_monthly=np.array(dust_monthly)*1e9
# total_monthly=np.array(total_monthly)*1e9

# fig, ax = plt.subplots()
# dates=pd.date_range(start="2008-01-01 00:00:00",end="2008-12-31 23:00:00",freq='60T')
# plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b'))
# plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=1))
# plt.plot(dates,total_yearly)
# title='Total aerosol domain average'
# plt.title(title)
# plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")
# plt.show()


# dates_month=pd.date_range(start="2008-01-01 01:00:00",end="2008-12-31 01:00:00",freq='M')
# bin_width=20
# fig, ax = plt.subplots()
# plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b'))
# plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=1))

# plt.bar(dates_month,dust_monthly,bin_width,label=r'Dust')
# plt.bar(dates_month,na_monthly,bin_width,label=r'Na',bottom=dust_monthly)
# plt.bar(dates_month,so4a_monthly,bin_width,label=r'SO$_4$',bottom=dust_monthly+na_monthly)
# plt.bar(dates_month,nh4a_monthly,bin_width,label=r'NH$_4$',bottom=dust_monthly+na_monthly+so4a_monthly)
# plt.bar(dates_month,ppm_monthly,bin_width,label=r'ppm',bottom=dust_monthly+na_monthly+so4a_monthly+nh4a_monthly)
# plt.bar(dates_month,pom_monthly,bin_width,label=r'pom',bottom=dust_monthly+na_monthly+so4a_monthly+nh4a_monthly+ppm_monthly)
# plt.bar(dates_month,ec_monthly,bin_width,label=r'BC',bottom=dust_monthly+na_monthly+so4a_monthly+nh4a_monthly+ppm_monthly+pom_monthly)
# plt.bar(dates_month,no3a_monthly,bin_width,label=r'NO$_3$',bottom=dust_monthly+na_monthly+so4a_monthly+nh4a_monthly+ppm_monthly+pom_monthly+ec_monthly)
# title='Aerosol composition'
# plt.title(title)
# plt.legend()
# plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")
# plt.show()
# gridded_LE_year=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_POLDER_Complete/gridded/LE_2008.nc')
# gridded_LE_year.timeseries(variable=['ys','yr'],biascorr=[2.5,1],label=['LE Biascorr=2.5','POLDER'],type_graph='mean',save=True,title='Domain Average time series',y_lim=[0,0.6])


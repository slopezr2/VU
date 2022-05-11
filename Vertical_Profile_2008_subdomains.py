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
from datamanager import graph_map
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
RdYlBu_r=cm.get_cmap('RdYlBu_r', 256)


#===Subdominio===
subdominios=['All','Netherlands','Benelux+Germany','Scandinavia','Iberia','East Europe']
#subdominios=['Netherlands']


for subdominio in subdominios:
    
    if subdominio=='Netherlands':
        lat_ini=60
        lat_end=76
        lon_ini=37
        lon_end=45
    elif subdominio=='All':
        lat_ini=0
        lat_end=-1
        lon_ini=0
        lon_end=-1    
    elif subdominio=='Benelux+Germany':
        lat_ini=50
        lat_end=85
        lon_ini=37
        lon_end=60
    elif subdominio=='Scandinavia':
        lat_ini=85
        lat_end=140
        lon_ini=37
        lon_end=100
    elif subdominio=='Iberia':
        lat_ini=2
        lat_end=45
        lon_ini=10
        lon_end=45
    elif subdominio=='East Europe':
        lat_ini=20
        lat_end=80
        lon_ini=55
        lon_end=100
    
    print(subdominio)
    h_file=Dataset('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_h_2008.nc')
    h_year=h_file.variables['h'][0:30,:,:,:]
    h_year=np.squeeze(np.nanmean(h_year,axis=0))
    h_year=np.squeeze(np.nanmean(h_year,axis=1))
    h_year=np.squeeze(np.nanmean(h_year,axis=1))
    
    
    no3a_yearly=np.zeros(12)
    na_yearly=np.zeros(12)
    so4a_yearly=np.zeros(12)
    nh4a_yearly=np.zeros(12)
    ppm_yearly=np.zeros(12)
    pom_yearly=np.zeros(12)
    ec_yearly=np.zeros(12)
    dust_yearly=np.zeros(12)
    total_yearly=np.zeros(12)
    
    std_total=np.zeros(12)
    
    
    i=0
    conc_year=np.zeros([12,12])
    conc_std_year=np.zeros([12,12])
    conc_month_all={}
    for month in (['01','02','03','04','05','06','07','08','09','10','11','12']):
    
        print(month)
        file_conc='/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_conc_2008'+month+'.nc'
        conc_file_month=Dataset(file_conc)
        conc_month=conc_file_month.variables['tpm10'][:,:,lat_ini:lat_end,lon_ini:lon_end]
        conc_std_month=np.squeeze(np.std(conc_month,axis=0))
        conc_std_month=np.squeeze(np.nanmean(conc_std_month,axis=1))
        conc_std_month=np.squeeze(np.nanmean(conc_std_month,axis=1))
        conc_month=np.squeeze(np.nanmean(conc_month,axis=0))
        conc_month=np.squeeze(np.nanmean(conc_month,axis=1))
        conc_month=np.squeeze(np.nanmean(conc_month,axis=1))
        conc_year[:,i]=conc_month
        conc_std_year[:,i]=conc_std_month
        i=i+1
    
        le_month=Dataset(file_conc)
        conc_month={}
        
        conc_month['no3a']=le_month.variables['no3a_c'][:,:,lat_ini:lat_end,lon_ini:lon_end]+le_month.variables['no3a_f'][:,:,lat_ini:lat_end,lon_ini:lon_end]
        conc_month['no3a']=np.squeeze(np.nanmean(conc_month['no3a'],axis=0))
        conc_month['no3a']=np.squeeze(np.nanmean(conc_month['no3a'],axis=1))
        conc_month['no3a']=np.squeeze(np.nanmean(conc_month['no3a'],axis=1))
        no3a_yearly=no3a_yearly+conc_month['no3a']
        
    
        conc_month['na']=le_month.variables['na_ccc'][:,:,lat_ini:lat_end,lon_ini:lon_end]+le_month.variables['na_cc'][:,:,lat_ini:lat_end,lon_ini:lon_end]+le_month.variables['na_c'][:,:,lat_ini:lat_end,lon_ini:lon_end]+le_month.variables['na_f'][:,:,lat_ini:lat_end,lon_ini:lon_end]+le_month.variables['na_ff'][:,:,lat_ini:lat_end,lon_ini:lon_end]
        conc_month['na']=np.squeeze(np.nanmean(conc_month['na'],axis=0))*3.26
        conc_month['na']=np.squeeze(np.nanmean(conc_month['na'],axis=1))
        conc_month['na']=np.squeeze(np.nanmean(conc_month['na'],axis=1))
        na_yearly=na_yearly+conc_month['na']
        
        conc_month['so4a']=le_month.variables['so4a_c'][:,:,lat_ini:lat_end,lon_ini:lon_end]+le_month.variables['so4a_f'][:,:,lat_ini:lat_end,lon_ini:lon_end]
        conc_month['so4a']=np.squeeze(np.nanmean(conc_month['so4a'],axis=0))
        conc_month['so4a']=np.squeeze(np.nanmean(conc_month['so4a'],axis=1))
        conc_month['so4a']=np.squeeze(np.nanmean(conc_month['so4a'],axis=1))
        so4a_yearly=so4a_yearly+conc_month['so4a']
        
        conc_month['nh4a']=le_month.variables['nh4a_f'][:,:,lat_ini:lat_end,lon_ini:lon_end]
        conc_month['nh4a']=np.squeeze(np.nanmean(conc_month['nh4a'],axis=0))
        conc_month['nh4a']=np.squeeze(np.nanmean(conc_month['nh4a'],axis=1))
        conc_month['nh4a']=np.squeeze(np.nanmean(conc_month['nh4a'],axis=1))
        nh4a_yearly=nh4a_yearly+conc_month['nh4a']
                                                      
        conc_month['ppm']=le_month.variables['ppm_c'][:,:,lat_ini:lat_end,lon_ini:lon_end]+le_month.variables['ppm_f'][:,:,lat_ini:lat_end,lon_ini:lon_end]
        conc_month['ppm']=np.squeeze(np.nanmean(conc_month['ppm'],axis=0))
        conc_month['ppm']=np.squeeze(np.nanmean(conc_month['ppm'],axis=1))
        conc_month['ppm']=np.squeeze(np.nanmean(conc_month['ppm'],axis=1))
        ppm_yearly=ppm_yearly+conc_month['ppm']
        
        conc_month['pom']=le_month.variables['pom_c'][:,:,lat_ini:lat_end,lon_ini:lon_end]+le_month.variables['pom_f'][:,:,lat_ini:lat_end,lon_ini:lon_end]
        conc_month['pom']=np.squeeze(np.nanmean(conc_month['pom'],axis=0))
        conc_month['pom']=np.squeeze(np.nanmean(conc_month['pom'],axis=1))
        conc_month['pom']=np.squeeze(np.nanmean(conc_month['pom'],axis=1))
        pom_yearly=pom_yearly+conc_month['pom']
        
        conc_month['ec']=le_month.variables['ec_c'][:,:,lat_ini:lat_end,lon_ini:lon_end]+le_month.variables['ec_f'][:,:,lat_ini:lat_end,lon_ini:lon_end]
        conc_month['ec']=np.squeeze(np.nanmean(conc_month['ec'],axis=0))
        conc_month['ec']=np.squeeze(np.nanmean(conc_month['ec'],axis=1))
        conc_month['ec']=np.squeeze(np.nanmean(conc_month['ec'],axis=1))
        ec_yearly=ec_yearly+conc_month['ec']
        
        conc_month['dust']=le_month.variables['dust_ccc'][:,:,lat_ini:lat_end,lon_ini:lon_end]+le_month.variables['dust_cc'][:,:,lat_ini:lat_end,lon_ini:lon_end]+le_month.variables['dust_c'][:,:,lat_ini:lat_end,lon_ini:lon_end]+le_month.variables['dust_f'][:,:,lat_ini:lat_end,lon_ini:lon_end]+le_month.variables['dust_ff'][:,:,lat_ini:lat_end,lon_ini:lon_end]
        conc_month['dust']=np.squeeze(np.nanmean(conc_month['dust'],axis=0))
        conc_month['dust']=np.squeeze(np.nanmean(conc_month['dust'],axis=1))
        conc_month['dust']=np.squeeze(np.nanmean(conc_month['dust'],axis=1))
        dust_yearly=dust_yearly+conc_month['dust']
        
        conc_month['total']=conc_month['no3a']+conc_month['na']+conc_month['so4a']+conc_month['nh4a']+conc_month['ppm']+conc_month['pom']+conc_month['ec']+conc_month['dust']
        total_yearly=total_yearly+conc_month['total']
        std_total=std_total+(np.std(conc_month['total']))
        conc_month_all[month]=conc_month
        
    
    no3a_yearly=np.array(no3a_yearly)*1e9/12
    na_yearly=np.array(na_yearly)*1e9/12
    so4a_yearly=np.array(so4a_yearly)*1e9/12
    nh4a_yearly=np.array(nh4a_yearly)*1e9/12
    ppm_yearly=np.array(ppm_yearly)*1e9/12
    pom_yearly=np.array(pom_yearly)*1e9/12
    ec_yearly=np.array(ec_yearly)*1e9/12
    dust_yearly=np.array(dust_yearly)*1e9/12
    total_yearly=np.array(total_yearly)*1e9/12
    
    
    layers=np.arange(1,13)
    layers_ticks= [str(numeric_string)[0:-2] for numeric_string in np.round(h_year)]
    bin_width=20
    fig, ax = plt.subplots()
    
    plt.barh(layers*bin_width*1.5,dust_yearly,height=bin_width,label=r'Dust')
    plt.barh(layers*bin_width*1.5,na_yearly,height=bin_width,label=r'SS',left=dust_yearly)
    plt.barh(layers*bin_width*1.5,so4a_yearly,height=bin_width,label=r'SO$_4$',left=dust_yearly+na_yearly)
    plt.barh(layers*bin_width*1.5,nh4a_yearly,height=bin_width,label=r'NH$_4$',left=dust_yearly+na_yearly+so4a_yearly)
    plt.barh(layers*bin_width*1.5,ppm_yearly,height=bin_width,label=r'ppm',left=dust_yearly+na_yearly+so4a_yearly+nh4a_yearly)
    plt.barh(layers*bin_width*1.5,pom_yearly,height=bin_width,label=r'pom',left=dust_yearly+na_yearly+so4a_yearly+nh4a_yearly+ppm_yearly)
    plt.barh(layers*bin_width*1.5,ec_yearly,height=bin_width,label=r'BC',left=dust_yearly+na_yearly+so4a_yearly+nh4a_yearly+ppm_yearly+pom_yearly)
    plt.barh(layers*bin_width*1.5,no3a_yearly,height=bin_width,label=r'NO$_3$',left=dust_yearly+na_yearly+so4a_yearly+nh4a_yearly+ppm_yearly+pom_yearly+ec_yearly)
    ax.set_yticks(layers*bin_width*1.5)
    ax.set_yticklabels(layers_ticks, rotation=0)
    ax.set_xlabel(r'Mean concentration [$\mu$g/m$^3$]',fontsize=12)
    ax.set_ylabel(r'Mean layer altitude [m]',fontsize=12)
    title='Aerosol composition per layer for '+subdominio
    plt.title(title)
    plt.legend()
    plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")
    plt.show()
    # gridded_LE_year=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_POLDER_Complete/gridded/LE_2008.nc')
    # gridded_LE_year.timeseries(variable=['ys','yr'],biascorr=[2.5,1],label=['LE Biascorr=2.5','POLDER'],type_graph='mean',save=True,title='Domain Average time series',y_lim=[0,0.6])
    
    
    conc_year=np.squeeze(np.mean(conc_year,axis=1))
    conc_std_year=np.squeeze(np.mean(conc_std_year,axis=1))
    
    
    fig, ax = plt.subplots()
    ax.plot(conc_year*1e9,h_year/1000,'*-',linewidth=3,markersize=12)
    #ax.fill_between(h_year,conc_year*1e9-conc_std_year*1e9,conc_year*1e9+conc_std_year*1e9,alpha=0.3)
    ax.set_xlabel(r'Mean concentration [$\mu$g/m$^3$]',fontsize=12)
    ax.set_ylabel(r'Mean layer altitude [km]',fontsize=12)
    #plt.ylim([0, 27])
    plt.savefig('./Figures/Vertical_profile_Aerosols for '+subdominio+'.png',format='png', dpi=1000,bbox_inches = "tight")
    plt.show()
    
    lat=conc_file_month.variables['lat']
    lon=conc_file_month.variables['lon']


    nan_conc=np.empty(conc_file_month.variables['tpm10'][:].shape)
    nan_conc[:]=np.NaN
    graph_map(variable=np.mean(nan_conc[:,0,lat_ini:lat_end,lon_ini:lon_end],axis=0)*1e9,lat=lat[lat_ini:lat_end],lon=lon[lon_ini:lon_end],vmin=0,vmax=200,cmap=RdYlGn.reversed(),date=False,save=True,title=subdominio,dpi=300,grid=True,save_title='Map '+subdominio,stock_image=True)
 
    months=['January','February','March','April','May','June','July','August','September','October','November','December']
    for imonth,month in enumerate(['01','02','03','04','05','06','07','08','09','10','11','12']):
        layers=np.arange(1,13)
        layers_ticks= [str(numeric_string)[0:-2] for numeric_string in np.round(h_year)]
        bin_width=20
        dust_monthly=conc_month_all[month]['dust']*1e9
        na_monthly=conc_month_all[month]['na']*1e9
        so4a_monthly=conc_month_all[month]['so4a']*1e9
        nh4a_monthly=conc_month_all[month]['nh4a']*1e9
        ppm_monthly=conc_month_all[month]['ppm']*1e9
        pom_monthly=conc_month_all[month]['pom']*1e9
        ec_monthly=conc_month_all[month]['ec']*1e9
        no3a_monthly=conc_month_all[month]['no3a']*1e9
        
        
        fig, ax = plt.subplots()
        
        plt.barh(layers*bin_width*1.5,dust_monthly,height=bin_width,label=r'Dust')
        plt.barh(layers*bin_width*1.5,na_monthly,height=bin_width,label=r'SS',left=dust_monthly)
        plt.barh(layers*bin_width*1.5,so4a_monthly,height=bin_width,label=r'SO$_4$',left=dust_monthly+na_monthly)
        plt.barh(layers*bin_width*1.5,nh4a_monthly,height=bin_width,label=r'NH$_4$',left=dust_monthly+na_monthly+so4a_monthly)
        plt.barh(layers*bin_width*1.5,ppm_monthly,height=bin_width,label=r'ppm',left=dust_monthly+na_monthly+so4a_monthly+nh4a_monthly)
        plt.barh(layers*bin_width*1.5,pom_monthly,height=bin_width,label=r'pom',left=dust_monthly+na_monthly+so4a_monthly+nh4a_monthly+ppm_monthly)
        plt.barh(layers*bin_width*1.5,ec_monthly,height=bin_width,label=r'BC',left=dust_monthly+na_monthly+so4a_monthly+nh4a_monthly+ppm_monthly+pom_monthly)
        plt.barh(layers*bin_width*1.5,no3a_monthly,height=bin_width,label=r'NO$_3$',left=dust_monthly+na_monthly+so4a_monthly+nh4a_monthly+ppm_monthly+pom_monthly+ec_monthly)
        ax.set_yticks(layers*bin_width*1.5)
        ax.set_yticklabels(layers_ticks, rotation=0)
        ax.set_xlabel(r'Mean concentration [$\mu$g/m$^3$]',fontsize=12)
        ax.set_ylabel(r'Mean layer altitude [m]',fontsize=12)
        title='Aerosol composition per layer for '+subdominio+' '+months[imonth]
        plt.title(title)
        plt.legend()
        plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")
        plt.show()
  
         
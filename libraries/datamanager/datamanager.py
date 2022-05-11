from netCDF4 import Dataset
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from cartopy import config
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy as cart
from datetime import datetime,timedelta
from datetime import date as date_datetime
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm
from os import listdir
from os.path import isfile, join
import matplotlib.dates as mdates
import mpl_scatter_density
import scipy.stats
from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
from matplotlib.ticker import ScalarFormatter
from matplotlib.colors import TwoSlopeNorm
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
from matplotlib import cm


viridis = cm.get_cmap('viridis', 256)
orange=cm.get_cmap('Oranges', 256)
bwr=cm.get_cmap('bwr', 256)


def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

def graph_map(lon=0,lat=0,vmin=0,vmax=0.3,date=True,vcenter=0,date_plot='',level=0,n_levels_plot=20,cmap=viridis,variable=1,variable_title='',title='default',title_on=True,save=False,extend='max',norm2=False,grid=False,ocean=False, orientation='vertical',lat_lim=[0,-1],lon_lim=[0,-1],subdomain=False,ocean_color=[1,1,1],shape_file=False,shape_name='',units='',dpi=1000,save_title=None,stock_image=False,stations=False,size_stations=15,facecolor=(0.6,0.6,0.6)):
    
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
    levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
    ax.set_facecolor(facecolor)
    if norm2==True:
        norm=TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)
    else:    
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
    
    if stations==False:
        plt.pcolormesh(lon, lat, variable, 
                        transform=ccrs.PlateCarree(),cmap=cmap, norm=norm,vmin=vmin, vmax=vmax,shading='auto')
    else:
        
        lat_sca=[]
        lon_sca=[]
        aux_sca=[]
        mask_aux=variable.mask
        for i in range(variable.shape[0]):
            for j in range(variable.shape[1]):
                if (mask_aux[i,j]==False and not(np.isnan(variable[i,j]))):
                    lat_sca.append(lat[i])
                    lon_sca.append(lon[j])
                    aux_sca.append(variable[i,j])
        plt.scatter(lon_sca, lat_sca,size_stations,aux_sca
        ,cmap=cmap, norm=norm,edgecolors=(0.1,0.1,0.1),linewidths=0.5)
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
    if shape_file==True:
        
        shape_feature = ShapelyFeature(Reader(shape_name).geometries(),
                                ccrs.PlateCarree(), facecolor='none')
        ax.add_feature(shape_feature,edgecolor='black',linewidth=0.8)


    if ocean==True:
           ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k',facecolor=ocean_color)

    if stock_image==True:
        ax.stock_img()
    if grid==True:
        gl = ax.gridlines(draw_labels=True,color='black', alpha=0.5,linestyle='--',zorder=100000)
        gl.top_labels = False
        gl.right_labels = False
    
    if subdomain==True:
        plt.xlim(lon[lon_lim[0]],lon[lon_lim[1]])
        plt.ylim(lat[lat_lim[0]],lat[lat_lim[1]])

    else:
        plt.xlim(lon[0],lon[-1])
        plt.ylim(lat[0],lat[-1])       
    if title=='default':
        title=variable
        if variable.ndim==4:
            title=title+'_level_'+str(level+1)
 
    title_plot=title
           
    if date==True:
        title_plot=title_plot+'  '+date_plot
        
    if title_on==True:    
        plt.title(title_plot)
    
    if orientation=='vertical': 
        cax = fig.add_axes([ax.get_position().x1+0.02,
                 ax.get_position().y0+0.09,
                 0.08,
                 ax.get_position().height-0.2])
    elif orientation=='horizontal':
        if units=='':
            deltay=0.12
        else:
            deltay=0.18
        cax = fig.add_axes([ax.get_position().x0,
                 ax.get_position().y0-deltay,
                 ax.get_position().width+0.02,
                 0.05])
    
    cbar=plt.colorbar(cax=cax,extend=extend,orientation=orientation)
    cbar.formatter.set_powerlimits((-3, 4))
    cbar.update_ticks()

    cbar.ax.set_title(units, y=1.05, x=0.5, rotation=0)
    
    
    
    if save==True:
        if save_title==None:
            if date==True:
                plt.savefig('./Figures/'+title_plot+'.png',format='png', dpi=dpi,bbox_inches = "tight")
            else:
                plt.savefig('./Figures/'+title+'.png',format='png', dpi=dpi,bbox_inches = "tight")
        else:
                plt.savefig('./Figures/'+save_title+'.png',format='png', dpi=dpi,bbox_inches = "tight")
    plt.show()



class DataManager_LE:

   def __init__(self,file):
        self.nc = Dataset(file, "r", format="NETCDF4")
        
   def graph_map(self,biascorr=1,variable='aod_550nm',time=0,sum_level=False,level=1,vmin=0, vmax=1,vcenter=0.5,norm2=False, n_levels_plot=10,cmap=viridis, date=True,title='default',title_on=True,type_graph='mean',ini_period=0,end_period=-1,step_period=1,save=False,ocean=False, extend='both',grid=False,orientation='vertical',lat_lim=[0,-1],lon_lim=[0,-1],subdomain=False,ocean_color=[1,1,1],shape_file=False,shape_name='',units='',dpi=1000,save_title=None,stations=False,size_stations=15,facecolor=(0.6,0.6,0.6)):     
       
       if 'lon' in self.nc.variables.keys():
           lon=self.nc.variables['lon'][:]
       else:
           lon=self.nc.variables['longitude'][:]
       if 'lat' in self.nc.variables.keys():
           lat=self.nc.variables['lat'][:]
       else:
           lat=self.nc.variables['latitude'][:]
       
       aux=self.nc.variables[variable][:]
       timestamp_2012 = datetime.timestamp(datetime.strptime('2008-1-1 00:00:00', '%Y-%m-%d %H:%M:%S'))
       if type_graph=='instant':
           
           timestamp=timestamp_2012+self.nc.variables['time'][time]
           date_plot=datetime.fromtimestamp(timestamp)
           date_plot=date_plot.strftime("%d-%b-%Y (%H:%M:%S)")
           if aux.ndim==3:
               var_graph= biascorr*aux[time,:,:]
               
           elif aux.ndim==2:
               var_graph=biascorr*aux[:,:]
              
           elif aux.ndim==4:
               
               if sum_level==True:
                   var_graph=biascorr*np.sum(aux[time,level,:,:],1)
                   
               else:
                   var_graph=biascorr*aux[time,level,:,:]
                  
        
       elif type_graph=='mean':
           
           if aux.ndim==3:
               var_graph=np.nanmean( biascorr*aux[ini_period:end_period:step_period,:,:],axis=0)
              
           elif aux.ndim==4:
               if sum_level==True:
                   aux=np.sum(aux,axis=1)
                   var_graph= np.nanmean(biascorr*aux[ini_period:end_period:step_period,:,:],axis=0)
                   
               else:
                   var_graph=np.nanmean( biascorr*aux[ini_period:end_period:step_period,level,:,:],axis=0)
                   
           elif aux.ndim==2:
                   var_graph= np.nanmean( biascorr*aux[:,:],axis=0)
                   
           
           date_plot_ini=datetime.fromtimestamp(timestamp_2012+self.nc.variables['time'][ini_period])
           date_plot_ini=date_plot_ini.strftime("%d-%b-%Y")
           date_plot_end=datetime.fromtimestamp(timestamp_2012+self.nc.variables['time'][end_period])
           date_plot_end=date_plot_end.strftime("%d-%b-%Y")
           
           date_plot=date_plot_ini+' - '+date_plot_end
           
           
       graph_map(lon=lon,lat=lat,vmin=vmin,vmax=vmax,vcenter=vcenter,n_levels_plot=n_levels_plot,date=date,date_plot=date_plot,title_on=title_on,cmap=cmap,variable_title=variable,variable=var_graph,title=title,save=save,extend=extend,norm2=norm2,grid=grid,ocean=ocean, orientation=orientation,lat_lim=lat_lim,lon_lim=lon_lim,subdomain=subdomain,ocean_color=ocean_color,shape_file=shape_file,shape_name=shape_name,units=units,dpi=dpi,save_title=save_title,stations=stations,size_stations=size_stations,facecolor=facecolor)
           
       
       


class DataManager_GRIDDED:

   def __init__(self,file):
        self.nc = Dataset(file, "r", format="NETCDF4")
     
   def metrics(self,biascorr=[2.5,1],show=False):
       aod_sat=self.nc.variables['yr'][:]
       aod_le=self.nc.variables['ys'][:]
       
       if aod_sat.ndim==4:
           aod_sat=biascorr[0]*aod_sat.reshape(aod_sat.shape[0]*aod_sat.shape[1]*aod_sat.shape[2]*aod_sat.shape[3])
           aod_le=biascorr[1]*aod_le.reshape(aod_le.shape[0]*aod_le.shape[1]*aod_le.shape[2]*aod_le.shape[3])
       else:
           aod_sat=biascorr[0]*aod_sat.reshape(aod_sat.shape[0]*aod_sat.shape[1]*aod_sat.shape[2])
           aod_le=biascorr[1]*aod_le.reshape(aod_le.shape[0]*aod_le.shape[1]*aod_le.shape[2])
       aux=aod_sat
       aod_sat=aod_sat[aux>0]
       aod_le=aod_le[aux>0]
       
       pearson=pearsonr(aod_sat,aod_le)
       mfb=2*((np.nanmean(aod_sat)-np.nanmean(aod_le))/(np.nanmean(aod_sat)+np.nanmean(aod_le)))
       rmse=(mean_squared_error(aod_sat,aod_le))**(1/2)
       r2=r2_score(aod_sat,aod_le)
     
       if show==True:
            print('Pearsons correlation: '  + str(pearson[0]))
            print('Mean Fractional Bias: '+ str(mfb))
            print('Root Mean Square Error: '+ str(rmse))
            print('R^2: ' +str(r2))
       return pearson[0],mfb,rmse,r2    
   def hist2d(self,biascorr=[2.5,1],cmap='Oranges',xlabel='Satellite',ylabel='LE',xlim=[0,1],ylim=[0,1],log=False,save=True,title=''):
       aod_sat=self.nc.variables['yr'][:]
       aod_le=self.nc.variables['ys'][:]
       if aod_sat.ndim==4:
           aod_sat=biascorr[0]*aod_sat.reshape(aod_sat.shape[0]*aod_sat.shape[1]*aod_sat.shape[2]*aod_sat.shape[3])
           aod_le=biascorr[1]*aod_le.reshape(aod_le.shape[0]*aod_le.shape[1]*aod_le.shape[2]*aod_le.shape[3])
       else:
           aod_sat=biascorr[0]*aod_sat.reshape(aod_sat.shape[0]*aod_sat.shape[1]*aod_sat.shape[2])
           aod_le=biascorr[1]*aod_le.reshape(aod_le.shape[0]*aod_le.shape[1]*aod_le.shape[2])
       aux=aod_sat
       aod_sat=aod_sat[aux>0]
       aod_le=aod_le[aux>0]
       pearson,mfb,rmse,r2=self.metrics(biascorr=biascorr,show=False)
       N=len(aod_sat)
       fig = plt.figure()
       ax = fig.add_subplot(1, 1, 1, projection='scatter_density')
       plt.grid(True, which="both", ls="-",alpha=0.3)
       density = ax.scatter_density((aod_sat), (aod_le),cmap=cmap)
       xx=np.arange(0,100)
       plt.plot(xx, xx,'k-',linewidth=1)
       plt.plot(xx, xx*0.5,'k-',linewidth=0.5)
       plt.plot(xx, xx*2,'k-',linewidth=0.5)
       fig.colorbar(density, label='Number of points per pixel')
       ax.set_xlabel(xlabel)
       ax.set_ylabel(ylabel)
       if log==True:
           ax.set_yscale('log')
           ax.set_xscale('log')
           for axis in [ax.xaxis, ax.yaxis]:
            formatter = ScalarFormatter()
            formatter.set_scientific(False)
            axis.set_major_formatter(formatter)
            
       #ax.text(xlim[0]*1.2, ylim[1]*0.2, )
      # ax.annotate('MFB= %2.2f\nRMSE= %2.2f\nPcorr= %2.2f\nR^2= %2.2f\nN=%1i'%(mfb,rmse,pearson,r2,N), xy=(0.02, 0.97), xycoords='axes fraction', fontsize=12,horizontalalignment='left', verticalalignment='top')
       ax.annotate('MFB= %2.2f'%(mfb), xy=(0.02, 0.97), xycoords='axes fraction', fontsize=12,horizontalalignment='left', verticalalignment='top')
       
       plt.xlim(xlim)
       plt.ylim(ylim)
       


       
            
       if save==True:
           plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")
       plt.show()

       return aod_sat,aod_le
   def graph_map(self,biascorr=1,variable='aod_550nm',level=0,facecolor=(0.6,0.6,0.6),time=100,vmin=0, vcenter=0.5,vmax=1,norm2=False, n_levels_plot=10,cmap=viridis, date=True,title='default',title_on=True,type_graph='mean',ini_period=0,end_period=-1,step_period=1,save=False,ocean=False, extend='both',grid=False,orientation='vertical',lat_lim=[0,-1],lon_lim=[0,-1],subdomain=False,ocean_color=[1,1,1],shape_file=False,shape_name='',units='',dpi=1000,save_title=None,stations=False,size_stations=15):     
       
       if 'lon' in self.nc.variables.keys():
           lon=self.nc.variables['lon'][:]
       else:
           lon=self.nc.variables['longitude'][:]
       if 'lat' in self.nc.variables.keys():
           lat=self.nc.variables['lat'][:]
       else:
           lat=self.nc.variables['latitude'][:]
           
           
           
       aux=np.ma.masked_invalid(self.nc.variables[variable][:])
       if aux.ndim==4:
           aux=np.ma.masked_invalid(self.nc.variables[variable][:,:,:,level])
           aux_yr=np.ma.masked_invalid(self.nc.variables['yr'][:])
           aux[aux<0]=np.nan
           aux[aux_yr<0]=np.nan
           aux_mask=np.ma.masked_invalid(self.nc.variables['yr'][:,:,:,level])
       
       else:
           aux=np.ma.masked_invalid(self.nc.variables[variable][:,:,:])
           aux_yr=np.ma.masked_invalid(self.nc.variables['yr'][:])
           aux[aux<0]=np.nan
           aux[aux_yr<0]=np.nan
           aux_mask=np.ma.masked_invalid(self.nc.variables['yr'][:,:,:])
       aux.mask=aux_mask.mask
       

       timestamp_2012 = datetime.timestamp(datetime.strptime('2008-1-1 00:00:00', '%Y-%m-%d %H:%M:%S'))
       if type_graph=='instant':
           
           timestamp=timestamp_2012+self.nc.variables['time'][time]
           date_plot=datetime.fromtimestamp(timestamp)
           date_plot=date_plot.strftime("%d-%b-%Y (%H:%M:%S)")
           if aux.ndim==3:
               var_graph=biascorr*aux[time,:,:]
             
           elif aux.ndim==4:
               var_graph=biascorr*aux[time,:,:,level]
               
       elif type_graph=='mean':
           if aux.ndim==3:
               var_graph=np.nanmean( biascorr*aux[ini_period:end_period:step_period,:,:],axis=0)
               
           elif aux.ndim==4:
               var_graph= np.nanmean( biascorr*aux[ini_period:end_period:step_period,:,:,level],axis=0)
               
           date_plot_ini=datetime.fromtimestamp(timestamp_2012+self.nc.variables['time'][ini_period])
           date_plot_ini=date_plot_ini.strftime("%d-%b-%Y")
           date_plot_end=datetime.fromtimestamp(timestamp_2012+self.nc.variables['time'][end_period])
           date_plot_end=date_plot_end.strftime("%d-%b-%Y")
           
           date_plot=date_plot_ini+' - '+date_plot_end
       graph_map(lon=lon,lat=lat,vmin=vmin,vmax=vmax,vcenter=vcenter,n_levels_plot=n_levels_plot,date=date,date_plot=date_plot,title_on=title_on,cmap=cmap,variable_title=variable,variable=var_graph,title=title,save=save,extend=extend,norm2=norm2,grid=grid,ocean=ocean, orientation=orientation,lat_lim=lat_lim,lon_lim=lon_lim,subdomain=subdomain,ocean_color=ocean_color,shape_file=shape_file,shape_name=shape_name,units=units,dpi=dpi,save_title=save_title,stations=stations,size_stations=size_stations,facecolor=facecolor)
           
       
   
    
   def transversal(self,biascorr=1,extern=False,calipso_mean=1,le_mean=1,layer_meters=600,variable='aod_550nm',facecolor=(0.6,0.6,0.6),time=100,vmin=0, vmax=1, n_levels_plot=10,cmap=viridis, date=True,title='default',title_on=True,type_graph='mean',ini_period=0,end_period=-1,step_period=1,save=False,ocean=False, extend='both',direction='row',position=10,cut=True,normalize=False):
       
       if 'lon' in self.nc.variables.keys():
           lon=self.nc.variables['lon'][:]
       else:
           lon=self.nc.variables['longitude'][:]
       if 'lat' in self.nc.variables.keys():
           lat=self.nc.variables['lat'][:]
       else:
           lat=self.nc.variables['latitude'][:]
           
       aux=np.ma.masked_invalid(self.nc.variables[variable][:])  
       
       
       if normalize==True:
           calipso_profile=np.ma.masked_invalid(self.nc.variables['yr'][:])
           if extern==False:
               le_mean=aux.mean()
               calipso_mean=calipso_profile.mean()
           
       if aux.ndim<=3:
           print('La variable no tiene dimension de altitud')
           return 0
       
       if normalize==True:
           aux=aux*calipso_mean/le_mean
       if type_graph=='mean':
           aux=np.nanmean(aux[ini_period:end_period:step_period,:,:,:],axis=0)
           if normalize==True:
               calipso_profile=np.nanmean(calipso_profile[ini_period:end_period:step_period,:,:,:-1],axis=0)
       elif type_graph=='instant': 
           aux=aux[time,:,:,:]
           if normalize==True:
               calipso_profile=calipso_profile[time,:,:,:-1]
       
       
       if variable=='ys':
           if direction=='row':
               aux=aux[position,:,:]
               
           elif direction=='col':
               aux=aux[:,position,:]
               
       elif variable=='yr':
           if direction=='row':
               aux=aux[position,:,:-1]

           elif direction=='col':
               aux=aux[:,position,:-1]
               
       orography=Dataset('/Users/santiago/Documents/LE_outputs/LE_Orography_meteo_20080101.nc')
       if direction=='row':
               if normalize==True:
                   calipso_profile=calipso_profile[position,:,:]
               orography=orography.variables['oro'][position,:]/1000
       elif direction=='col':
               if normalize==True:
                   calipso_profile=calipso_profile[:,position,:]
               orography=orography.variables['oro'][:,position]/1000
       
       
       if variable=='ys':
           aux[~np.isfinite(calipso_profile)]=np.nan
       
       
       nx=aux.shape[0]
       nz=aux.shape[1]
        
       cross=np.zeros([nx,nz+2])
       cross[:,0]=0
       cross[:,1]=orography
        
       xx=np.zeros([nx,nz+2])
       lev_ht_1km=np.arange(0,nz)*layer_meters/1000
       for i in range(nz):
            cross[:,i+1]=lev_ht_1km[i]
            if direction=='row':
                xx[:,i]=lon
            elif direction=='col':
                xx[:,i]=lat
       levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
       norm = BoundaryNorm(levels, ncolors=256, clip=True)
       fig, ax = plt.subplots()
       ax.set_facecolor(facecolor)
       
       timestamp_2012 = datetime.timestamp(datetime.strptime('2008-1-1 00:00:00', '%Y-%m-%d %H:%M:%S'))
       im=ax.pcolor(xx[:,1:-1],cross[:,1:-1],aux[:,:],vmin=vmin,vmax=vmax,cmap=cmap,norm=norm)
       ax.fill_between(xx[:,0],cross[:,0],orography,color='k')
       if type_graph=='instant':  
           
           timestamp=timestamp_2012+self.nc.variables['time'][time]
           date_plot=datetime.fromtimestamp(timestamp)
           date_plot=date_plot.strftime("%d-%b-%Y (%H:%M:%S)")
       
       elif type_graph=='mean':
           date_plot_ini=datetime.fromtimestamp(timestamp_2012+self.nc.variables['time'][ini_period])
           date_plot_ini=date_plot_ini.strftime("%d-%b-%Y")
           date_plot_end=datetime.fromtimestamp(timestamp_2012+self.nc.variables['time'][end_period])
           date_plot_end=date_plot_end.strftime("%d-%b-%Y")
           date_plot=date_plot_ini+' - '+date_plot_end
       
        
       
       ax.set_ylim(0,10)
       if direction=='row':
           ax.set_xlim(min(lon),max(lon))
           lon_ticks=np.linspace(min(lon),max(lon), 10, dtype=int)
           ax.set_xticks(lon_ticks)
           ax.set_ylabel('Altitude Km')
           ax.set_xlabel('Longitude ')
           
           # Set scond x-axis
           ax2 = ax.twiny()
            
           lat_ticks=np.linspace(min(lat), max(lat), 10, dtype=int)
           ax2.set_xticks(lat_ticks)
           ax2.set_xlim(min(lat),max(lat))
            
           ax2.xaxis.set_ticks_position('bottom') # set the position of the second x-axis to bottom
           ax2.xaxis.set_label_position('bottom') # set the position of the second x-axis to bottom
           ax2.spines['bottom'].set_position(('outward', 36))
           ax2.set_xlabel('Latitude')
       elif direction=='col':
           ax.set_xlim(min(lat),max(lat))
           lat_ticks=np.linspace(min(lat),max(lat), 10, dtype=int)
           ax.set_xticks(lat_ticks)
           ax.set_ylabel('Altitude Km')
           ax.set_xlabel('Latitude ')
           
           # Set scond x-axis
           ax2 = ax.twiny()
            
           lon_ticks=np.linspace(min(lon), max(lon), 10, dtype=int)
           ax2.set_xticks(lon_ticks)
           ax2.set_xlim(min(lon),max(lon))
            
           ax2.xaxis.set_ticks_position('bottom') # set the position of the second x-axis to bottom
           ax2.xaxis.set_label_position('bottom') # set the position of the second x-axis to bottom
           ax2.spines['bottom'].set_position(('outward', 36))
           ax2.set_xlabel('Longitude')
        
       if title=='default':
           title=variable
           
 
       title_plot=title
           
       if date==True:
           title_plot=title_plot+'  '+date_plot
           
       if title_on==True:    
           plt.title(title_plot)
       
       cax = fig.add_axes([ax.get_position().x1+0.02,
                    ax.get_position().y0+0.09,
                    0.08,
                    ax.get_position().height-0.2])
       
       cbar=plt.colorbar(im,cax=cax,extend=extend)
       cbar.formatter.set_powerlimits((-3, 4))
       cbar.update_ticks()
       if save==True:
           if date==True:
               plt.savefig('./Figures/'+title_plot+'.png',format='png', dpi=1000,bbox_inches = "tight")
           else:
               plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")
                   
       plt.show()
       
       if cut==True:
           aux_cut=np.zeros([lat.shape[0],lon.shape[0]])
           aux_cut[aux_cut==0]=np.nan
           if direction=='row':
               aux_cut[position,:]=1e6
           elif direction=='col':
               aux_cut[:,position]=1e6
           fig2, ax2 = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
           levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
           ax2.set_facecolor((1, 1, 1))
           norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
           
           plt.pcolormesh(lon, lat, aux_cut, 
                                  transform=ccrs.PlateCarree(),cmap=cmap, norm=norm,shading='auto')
            
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
           ax2.add_feature(european_borders,edgecolor='black',linewidth=0.2)
           ax2.add_feature(coastlines,edgecolor='black',linewidth=1)
           ax2.stock_img()
           plt.xlim(-15,35)
           plt.ylim(35,70)
           plt.title(title)
                    
            # cax = fig.add_axes([ax.get_position().x1+0.02,
            #          ax.get_position().y0+0.09,
            #          0.08,
            #          ax.get_position().height-0.2])
            
            #plt.colorbar(cax=cax,extend='both')
           if save==True:
               if date==True:
                   plt.savefig('./Figures/Cut_'+title_plot+'.png',format='png', dpi=1000,bbox_inches = "tight")
               else:
                   plt.savefig('./Figures/Cut_'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")
           plt.show()
       
       
       return fig,ax
    
   def graph_diff_map(self,biascorr=1,level=0,facecolor=(0.6,0.6,0.6),time=100,vmin=0, vcenter=0.5,vmax=1,norm2=False, n_levels_plot=10,cmap=bwr, date=False,title='default',title_on=True,type_graph='mean',ini_period=0,end_period=-1,step_period=1,save=False,ocean=False, extend='both',grid=False,orientation='vertical',lat_lim=[0,-1],lon_lim=[0,-1],subdomain=False,ocean_color=[1,1,1],shape_file=False,shape_name='',units='',dpi=1000,save_title=None,stations=False,size_stations=15,porcentage_diff=False):     
       
       if 'lon' in self.nc.variables.keys():
           lon=self.nc.variables['lon'][:]
       else:
           lon=self.nc.variables['longitude'][:]
       if 'lat' in self.nc.variables.keys():
           lat=self.nc.variables['lat'][:]
       else:
           lat=self.nc.variables['latitude'][:]

       aux=np.ma.masked_invalid(self.nc.variables['ys'][:])
       if aux.ndim==4:
           aux=np.ma.masked_invalid(self.nc.variables['ys'][:,:,:,level])
           aux_yr=np.ma.masked_invalid(self.nc.variables['yr'][:])
           aux[aux<0]=np.nan
           aux[aux_yr<0]=np.nan
           aux = aux_yr-aux
           aux_mask=np.ma.masked_invalid(self.nc.variables['yr'][:,:,:,level])
       
       else:
           aux=np.ma.masked_invalid(self.nc.variables['ys'][:,:,:])
           aux_yr=np.ma.masked_invalid(self.nc.variables['yr'][:])
           aux[aux<0]=np.nan
           aux[aux_yr<0]=np.nan
           aux = aux_yr-aux
           aux_mask=np.ma.masked_invalid(self.nc.variables['yr'][:,:,:])
       aux.mask=aux_mask.mask
       
       if porcentage_diff==True:
           aux = (aux/aux_yr)*100
       timestamp_2012 = datetime.timestamp(datetime.strptime('2008-1-1 00:00:00', '%Y-%m-%d %H:%M:%S'))
       if type_graph=='instant':
           
           timestamp=timestamp_2012+self.nc.variables['time'][time]
           date_plot=datetime.fromtimestamp(timestamp)
           date_plot=date_plot.strftime("%d-%b-%Y (%H:%M:%S)")
           if aux.ndim==3:
               var_graph=biascorr*aux[time,:,:]
             
           elif aux.ndim==4:
               var_graph=biascorr*aux[time,:,:,level]
               
       elif type_graph=='mean':
           if aux.ndim==3:
               var_graph=np.nanmean( biascorr*aux[ini_period:end_period:step_period,:,:],axis=0)
               
           elif aux.ndim==4:
               var_graph= np.nanmean( biascorr*aux[ini_period:end_period:step_period,:,:,level],axis=0)
               
           date_plot_ini=datetime.fromtimestamp(timestamp_2012+self.nc.variables['time'][ini_period])
           date_plot_ini=date_plot_ini.strftime("%d-%b-%Y")
           date_plot_end=datetime.fromtimestamp(timestamp_2012+self.nc.variables['time'][end_period])
           date_plot_end=date_plot_end.strftime("%d-%b-%Y")
           
           date_plot=date_plot_ini+' - '+date_plot_end
       graph_map(lon=lon,lat=lat,vmin=vmin,vmax=vmax,vcenter=vcenter,n_levels_plot=n_levels_plot,date=date,date_plot=date_plot,title_on=title_on,cmap=cmap,variable_title='diff obs-LE',variable=var_graph,title='diff obs-LE',save=save,extend=extend,norm2=norm2,grid=grid,ocean=ocean, orientation=orientation,lat_lim=lat_lim,lon_lim=lon_lim,subdomain=subdomain,ocean_color=ocean_color,shape_file=shape_file,shape_name=shape_name,units=units,dpi=dpi,save_title=save_title,stations=stations,size_stations=size_stations,facecolor=facecolor)
           
       
   
    
   
   def timeseries(self,biascorr=[1,1],variable=['yr','ys'],label=['POLDER','LE'],title='default',title_on=True,type_graph='sum',ini_period=0,end_period=-1,step_period=1,save=False,y_lim=[0,1],windows=24):
       fig, ax = plt.subplots()
       if type_graph=='sum':
           
           for variable_aux in range(len(variable)):
               aux=self.nc.variables[variable[variable_aux]][ini_period:end_period:step_period,:,:]
               time=self.nc.variables['time'][ini_period:end_period:step_period].data
               timestamp_2012 = datetime.timestamp(datetime.strptime('2008-1-1 00:00:00', '%Y-%m-%d %H:%M:%S'))
               timestamp=np.round(time).astype(int)
               date_plot=[datetime.fromtimestamp(fs+timestamp_2012) for fs in timestamp]
               plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b'))
               plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=1))
               aux=pd.Series(np.squeeze(biascorr[variable_aux]*np.nansum(np.nansum(aux,axis=2),axis=1)))
               plt.plot(date_plot,aux.rolling(window=windows,min_periods=1).mean(),label=label[variable_aux])
               ax.set_ylim(y_lim)
               plt.gcf().autofmt_xdate()
           if title_on==True:
               plt.title(title)
           plt.legend()
           if save==True:
               plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")
           plt.show()
           
       elif type_graph=='mean':
           
           for variable_aux in range(len(variable)):
               aux=self.nc.variables[variable[variable_aux]][ini_period:end_period:step_period,:,:]
               time=self.nc.variables['time'][ini_period:end_period:step_period].data
               aux[aux<0] = np.NaN
               timestamp_2012 = datetime.timestamp(datetime.strptime('2008-1-1 00:00:00', '%Y-%m-%d %H:%M:%S'))
               timestamp=np.round(time).astype(int)
               date_plot=[datetime.fromtimestamp(fs+timestamp_2012) for fs in timestamp]
               aux=pd.Series(np.squeeze(biascorr[variable_aux]*np.nanmedian(np.nanmedian(aux,axis=2),axis=1)))
               plt.plot(date_plot,aux.rolling(window=windows,min_periods=1).mean(),label=label[variable_aux])
               plt.gcf().autofmt_xdate()
               ax.set_ylim(y_lim)
           if title_on==True:
               plt.title(title)
           plt.legend()
           
           if save==True:
               plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")
           plt.show()

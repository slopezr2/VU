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
class DataManager_LE:

   def __init__(self,file):
        self.nc = Dataset(file, "r", format="NETCDF4")
        
   def graph_map(self,biascorr=1,variable='aod_550nm',time=100,level=1,vmin=0, vmax=1, n_levels_plot=10,cmap='Oranges', date=True,title='default',title_on=True,type_graph='instant',ini_period=0,end_period=-1,step_period=1,save=False):     
       
       lon=self.nc.variables['lon'][:]
       lat=self.nc.variables['lat'][:]
       aux=self.nc.variables[variable][:]
       fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
       levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
    
    
        # pick the desired colormap, sensible levels, and define a normalization
        # instance which takes data values and translates those into levels.
        
       norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
       if vmax==None or vmin==None:
           levels = n_levels_plot
       else:
           levels = np.linspace(vmin, vmax, n_levels_plot)
       timestamp_2012 = datetime.timestamp(datetime.strptime('2008-1-1 00:00:00', '%Y-%m-%d %H:%M:%S'))
       if type_graph=='instant':
           
           timestamp=timestamp_2012+self.nc.variables['time'][time]
           date_plot=datetime.fromtimestamp(timestamp)
           date_plot=date_plot.strftime("%d-%b-%Y (%H:%M:%S)")
           if aux.ndim==3:
               plt.pcolormesh(lon, lat, biascorr*aux[time,:,:], 
                          transform=ccrs.PlateCarree(),cmap=cmap, norm=norm,vmin=vmin, vmax=vmax,shading='auto')
           elif aux.ndim==4:
               plt.pcolormesh(lon, lat,  biascorr*aux[time,level,:,:], 
                          transform=ccrs.PlateCarree(),cmap=cmap, norm=norm,vmin=vmin, vmax=vmax,shading='auto')
        
       elif type_graph=='mean':
           if aux.ndim==3:
               plt.pcolormesh(lon, lat, np.nanmean( biascorr*aux[ini_period:end_period:step_period,:,:],axis=0), 
                              transform=ccrs.PlateCarree(),cmap=cmap, norm=norm,vmin=vmin, vmax=vmax,shading='auto')
           elif aux.ndim==4:
               plt.pcolormesh(lon, lat, np.nanmean( biascorr*aux[ini_period:end_period:step_period,level,:,:],axis=0), 
                              transform=ccrs.PlateCarree(),cmap=cmap, norm=norm,vmin=vmin, vmax=vmax,shading='auto')
           date_plot_ini=datetime.fromtimestamp(timestamp_2012+self.nc.variables['time'][ini_period])
           date_plot_ini=date_plot_ini.strftime("%d-%b-%Y")
           date_plot_end=datetime.fromtimestamp(timestamp_2012+self.nc.variables['time'][end_period])
           date_plot_end=date_plot_end.strftime("%d-%b-%Y")
           
           date_plot=date_plot_ini+' - '+date_plot_end

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
       if title=='default':
           title=variable
           if aux.ndim==4:
               title=title+'_level_'+str(level+1)
 
       title_plot=title
           
       if date==True:
           title_plot=title_plot+'  '+date_plot
           
       if title_on==True:    
           plt.title(title_plot)
       
       cax = fig.add_axes([ax.get_position().x1+0.02,
                    ax.get_position().y0+0.09,
                    0.08,
                    ax.get_position().height-0.2])
       
       plt.colorbar(cax=cax,extend='both')
       
       if save==True:
           if date==True:
               plt.savefig('./Figures/'+title_plot+'.png',format='png', dpi=1000,bbox_inches = "tight")
           else:
               plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")
                   
       plt.show()

class DataManager_POLDER:   
    def __init__(self,path):
        self.path = path
    def graph_map(self,variable='AOD565',time=100,vmin=0, vmax=1, n_levels_plot=10,cmap='Oranges',title='default',title_on=True,type_graph='instant',ini_period=0,end_period=-1,step_period=1,save=False):     
        onlyfiles = [f for f in listdir(self.path) if isfile(join(self.path, f)) and f[0:12]=='GRASP_POLDER']

        onlyfiles.sort()

        #polder=Dataset('/Users/santiago/Documents/POLDER/GRASP_POLDER_L3_20120813.01degree.nc')
        polder=Dataset(self.path+onlyfiles[0])
        aux=polder.variables[variable][:]
        aod_polder=np.empty([aux.shape[0],aux.shape[1],len(onlyfiles)])
        i=0
        aod_polder[:,:,i]=aux.filled(fill_value=np.nan)
        
        for file_name in onlyfiles[1:]:
            if file_name[0:12]=='GRASP_POLDER':
                i=i+1
                polder=Dataset(self.path+file_name)
                aux=polder.variables[variable][:]
                aod_polder[:,:,i]=aux.filled(fill_value=np.nan)
            
                latitude=polder.variables['Latitude'][:]
                longitude=polder.variables['Longitude'][:]
                
            
        fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
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
        n_levels_plot=n_levels_plot
        vmin=vmin
        vmax=vmax
        levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
        
        
        # pick the desired colormap, sensible levels, and define a normalization
        # instance which takes data values and translates those into levels.
        #cmap = cmap
        #cmap.set_bad([0.7, 0.7, 0.7],1.)
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        # im=ax.contour(longitude[0,1652:2148],latitude[0:350,0],aod_polder[0:350,1652:2148,i],levels=levels,cmap=custom_ramp,extend='both')
        if not(title=='default') and title_on==True:
            plt.title(title)      
        im=plt.pcolor(longitude[0,1652:2148],latitude[0:350,0],np.nanmean(aod_polder[0:350,1652:2148,:],axis=2),cmap=cmap,norm=norm,vmin=vmin,vmax=vmax)
        cax = fig.add_axes([ax.get_position().x1+0.02,
                            ax.get_position().y0+0.09,
                            0.08,
                            ax.get_position().height-0.2])
         
        plt.colorbar(im,cax=cax,extend='both')
        if save==True:
            plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")
        plt.show()

class DataManager_GRIDDED:

   def __init__(self,file):
        self.nc = Dataset(file, "r", format="NETCDF4")
        
   def graph_map(self,biascorr=1,variable='aod_550nm',time=100,vmin=0, vmax=1, n_levels_plot=10,cmap='Oranges', date=True,title='default',title_on=True,type_graph='instant',ini_period=0,end_period=-1,step_period=1,save=False):     
       
       lon=self.nc.variables['longitude'][:]
       lat=self.nc.variables['latitude'][:]
       aux=self.nc.variables[variable][:]
       aux_mask=self.nc.variables['yr'][:]
       aux.mask=aux_mask.mask
       fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
       levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
       ax.set_facecolor((0.6, 0.6, 0.6))
    
        # pick the desired colormap, sensible levels, and define a normalization
        # instance which takes data values and translates those into levels.
        
       norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
       if vmax==None or vmin==None:
           levels = n_levels_plot
       else:
           levels = np.linspace(vmin, vmax, n_levels_plot)
       timestamp_2012 = datetime.timestamp(datetime.strptime('2008-1-1 00:00:00', '%Y-%m-%d %H:%M:%S'))
       if type_graph=='instant':
           
           timestamp=timestamp_2012+self.nc.variables['time'][time]
           date_plot=datetime.fromtimestamp(timestamp)
           date_plot=date_plot.strftime("%d-%b-%Y (%H:%M:%S)")
           if aux.ndim==3:
               plt.pcolormesh(lon, lat, biascorr*aux[time,:,:], 
                          transform=ccrs.PlateCarree(),cmap=cmap, norm=norm,vmin=vmin, vmax=vmax,shading='auto')
           elif aux.ndim==4:
               plt.pcolormesh(lon, lat,  biascorr*aux[time,:,:,0], 
                          transform=ccrs.PlateCarree(),cmap=cmap, norm=norm,vmin=vmin, vmax=vmax,shading='auto')
        
       elif type_graph=='mean':
           if aux.ndim==3:
               plt.pcolormesh(lon, lat, np.nanmean( biascorr*aux[ini_period:end_period:step_period,:,:],axis=0), 
                              transform=ccrs.PlateCarree(),cmap=cmap, norm=norm,vmin=vmin, vmax=vmax,shading='auto')
           elif aux.ndim==4:
               plt.pcolormesh(lon, lat, np.nanmean( biascorr*aux[ini_period:end_period:step_period,:,:,0],axis=0), 
                              transform=ccrs.PlateCarree(),cmap=cmap, norm=norm,vmin=vmin, vmax=vmax,shading='auto')
           date_plot_ini=datetime.fromtimestamp(timestamp_2012+self.nc.variables['time'][ini_period])
           date_plot_ini=date_plot_ini.strftime("%d-%b-%Y")
           date_plot_end=datetime.fromtimestamp(timestamp_2012+self.nc.variables['time'][end_period])
           date_plot_end=date_plot_end.strftime("%d-%b-%Y")
           
           date_plot=date_plot_ini+' - '+date_plot_end

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
       
       plt.colorbar(cax=cax,extend='both')
       
       if save==True:
           if date==True:
               plt.savefig('./Figures/'+title_plot+'.png',format='png', dpi=1000,bbox_inches = "tight")
           else:
               plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")
                   
       plt.show()
       return fig,ax
   def timeseries(self,biascorr=[1,1],variable=['yr','ys'],label=['POLDER','LE'],title='default',title_on=True,type_graph='sum',ini_period=0,end_period=-1,step_period=1,save=False):
       
       if type_graph=='sum':
           
           for variable_aux in range(len(variable)):
               aux=self.nc.variables[variable[variable_aux]][ini_period:end_period:step_period,:,:,0]
               time=self.nc.variables['time'][ini_period:end_period:step_period].data
               timestamp_2012 = datetime.timestamp(datetime.strptime('2008-1-1 00:00:00', '%Y-%m-%d %H:%M:%S'))
               timestamp=np.round(time).astype(int)
               date_plot=[datetime.fromtimestamp(fs+timestamp_2012) for fs in timestamp]
               plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%b'))
               plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=1))
               
               plt.plot(date_plot,biascorr[variable_aux]*np.nansum(np.nansum(aux,axis=2),axis=1),label=label[variable_aux])
               plt.gcf().autofmt_xdate()
           if title_on==True:
               plt.title(title)
           plt.legend()
           if save==True:
               plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")
           plt.show()
           
       elif type_graph=='mean':
           
           for variable_aux in range(len(variable)):
               aux=self.nc.variables[variable[variable_aux]][ini_period:end_period:step_period,:,:,0]
               time=self.nc.variables['time'][ini_period:end_period:step_period].data
               timestamp_2012 = datetime.timestamp(datetime.strptime('2008-1-1 00:00:00', '%Y-%m-%d %H:%M:%S'))
               timestamp=np.round(time).astype(int)
               date_plot=[datetime.fromtimestamp(fs+timestamp_2012) for fs in timestamp]
               
               plt.plot(date_plot,biascorr[variable_aux]*np.nanmean(np.nanmean(aux,axis=2),axis=1),label=label[variable_aux])
           if title_on==True:
               plt.title(title)
           plt.legend()
           if save==True:
               plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")
           plt.show()
class DataManager_CAMS:

   def __init__(self,file):
        self.nc = Dataset(file, "r", format="NETCDF4")
        
   def graph_map(self,biascorr=1,variable='aod_550nm',time=100,level=1,vmin=0, vmax=1, n_levels_plot=10,cmap='Oranges', date=True,title='default',title_on=True,type_graph='instant',ini_period=0,end_period=-1,step_period=1,save=False,lon_lim=[-25,44],lat_lim=[30,69]):     
       
       lon=self.nc.variables['longitude'][:]
       lat=self.nc.variables['latitude'][:]
       lat=np.flip(lat)
       aux=self.nc.variables[variable][:]
       aux=np.flip(aux,axis=1)
       fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
       levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
    
    
        # pick the desired colormap, sensible levels, and define a normalization
        # instance which takes data values and translates those into levels.
        
       norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
       base = date_datetime(1900, 1, 1)
       if vmax==None or vmin==None:
           levels = n_levels_plot
       else:
           levels = np.linspace(vmin, vmax, n_levels_plot)
       if type_graph=='instant':
           
           timestamp=self.nc.variables['time'][time].data
           
           date_plot=base+timedelta(hours=int(timestamp))
           
           date_plot=date_plot.strftime("%d-%b-%Y (%H:%M:%S)")
           if aux.ndim==3:
               plt.pcolormesh(lon, lat, biascorr*aux[time,:,:], 
                          transform=ccrs.PlateCarree(),cmap=cmap, norm=norm,vmin=vmin, vmax=vmax,shading='auto')
           elif aux.ndim==4:
               plt.pcolormesh(lon, lat,  biascorr*aux[time,level,:,:], 
                          transform=ccrs.PlateCarree(),cmap=cmap, norm=norm,vmin=vmin, vmax=vmax,shading='auto')
        
       elif type_graph=='mean':
           if aux.ndim==3:
               plt.pcolormesh(lon, lat, np.nanmean( biascorr*aux[ini_period:end_period:step_period,:,:],axis=0), 
                              transform=ccrs.PlateCarree(),cmap=cmap, norm=norm,vmin=vmin, vmax=vmax,shading='auto')
           elif aux.ndim==4:
               plt.pcolormesh(lon, lat, np.nanmean( biascorr*aux[ini_period:end_period:step_period,level,:,:],axis=0), 
                              transform=ccrs.PlateCarree(),cmap=cmap, norm=norm,vmin=vmin, vmax=vmax,shading='auto')
           date_plot_ini=base+timedelta(hours=int(self.nc.variables['time'][ini_period].data))
           date_plot_ini=date_plot_ini.strftime("%d-%b-%Y")
           date_plot_end=base+timedelta(hours=int(self.nc.variables['time'][end_period].data))
           date_plot_end=date_plot_end.strftime("%d-%b-%Y")
           
           date_plot=date_plot_ini+' - '+date_plot_end

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
       plt.xlim(lon_lim[0],lon_lim[-1])
       plt.ylim(lat_lim[0],lat_lim[-1])
       if title=='default':
           title=variable
           if aux.ndim==4:
               title=title+'_level_'+str(level+1)
 
       title_plot=title
           
       if date==True:
           title_plot=title_plot+'  '+date_plot
           
       if title_on==True:    
           plt.title(title_plot)
       
       cax = fig.add_axes([ax.get_position().x1+0.02,
                    ax.get_position().y0+0.09,
                    0.08,
                    ax.get_position().height-0.2])
       
       plt.colorbar(cax=cax,extend='both')
       
       if save==True:
           if date==True:
               plt.savefig('./Figures/'+title_plot+'.png',format='png', dpi=1000,bbox_inches = "tight")
           else:
               plt.savefig('./Figures/'+title+'.png',format='png', dpi=1000,bbox_inches = "tight")
                   
       plt.show()        
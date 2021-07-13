
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from cartopy import config
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm
from os import listdir,scandir
from os.path import isfile, join
import pandas as pd
import datetime


#Create listing table
template=pd.read_csv('/Users/santiago/scratch/CAMS-61/WP6130/CSO-data/S5p/listing-NO2-CAMS.csv',delimiter=';')
list_CSO={name: [] for name in list(template.columns)}

TERRA='MOD04_L2'
AQUA='MYD04_L2'
products=[TERRA,AQUA]
mypath='/Users/santiago/Documents/MODIS/2008'
CSO_path='/Users/santiago/Documents/MODIS/CSO'
subfolders = [ f.path for f in scandir(mypath) if f.is_dir() ]
subfolders.sort()
a=0
for i in subfolders:
    for product in products:
        month_path=i+'/'+product+'/2008/'
        dayfolders = [ f.path for f in scandir(month_path) if f.is_dir() ]
        dayfolders.sort()
        for day in dayfolders:
            onlyfiles = [f for f in listdir(day) if isfile(join(day, f))]

            onlyfiles.sort()
            for file_name in onlyfiles[1:]:
                if file_name[0:8]==product:
                    # Reading modis data
                    modis=Dataset(day+'/'+file_name)
                    aod_modis=modis.variables['Optical_Depth_Land_And_Ocean'][:]
                    latitude_modis=modis.variables['Latitude'][:]
                    longitude_modis=modis.variables['Longitude'][:]
                    latitude_modis=latitude_modis.reshape((latitude_modis.shape[0]*latitude_modis.shape[1]))
                    longitude_modis=longitude_modis.reshape((longitude_modis.shape[0]*longitude_modis.shape[1]))
                    aod_modis=aod_modis.reshape((aod_modis.shape[0]*aod_modis.shape[1]))
                    
                    
                    #Creating netcdf4
                    day_year=file_name[10:17]
                    day_year=datetime.datetime.strptime(day_year, '%Y%j')
                    day_listing=day_year.strftime('%Y-%m-%d')
                    day_year=day_year.strftime('%Y%m%d')
                    
                    hour_file='T'+file_name[18:20]+':'+file_name[20:22]+':00.00'
                    orbit=file_name[27:-4]
                    
                    cso_file='CSO_modis_'+product+'_'+day_year+file_name[17:-4]+'.nc'
                    root_grp = Dataset(CSO_path+'/'+cso_file, 'w', format='NETCDF4')
                    root_grp.description = 'modis observation following CSO estructure'
                    
                    
                    # dimensions
                    npixel = latitude_modis[np.logical_not(aod_modis.mask)].size 
                    root_grp.createDimension('pixel', npixel)
                    root_grp.createDimension('retr', 1)
                    root_grp.createDimension('corner', 4)
                    # variables
                    pixel = root_grp.createVariable('pixel', 'i4', ('pixel',),fill_value=-2147483647)
                    retr = root_grp.createVariable('retr', 'i4', ('retr',))
                    longitude = root_grp.createVariable('longitude', 'f4', ('pixel',),fill_value=3.4028235e+38)
                    longitude_bounds = root_grp.createVariable('longitude_bounds', 'f4', ('pixel','corner',),fill_value=3.4028235e+38)
                    latitude = root_grp.createVariable('latitude', 'f4', ('pixel',),fill_value=3.4028235e+38)
                    latitude_bounds = root_grp.createVariable('latitude_bounds', 'f4', ('pixel','corner',),fill_value=3.4028235e+38)
                    AOD550 = root_grp.createVariable('AOD550', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
                    
                    
                    #Attributes
                    pixel.long_name= 'pixel index'
                    pixel.units= '1'
                   
                    longitude.long_name= 'pixel center longitude'
                    longitude.units= 'degrees_east'
                    longitude.standard_name= 'longitude'
                    longitude.valid_min: -180.0
                    longitude.valid_max: 180.0
                    
                    longitude_bounds.long_name= 'longitude_bounds'
                    longitude_bounds.units= 'degrees_east'
                    
                    latitude.long_name= 'pixel center latitude'
                    latitude.units= 'degrees_north'
                    latitude.standard_name= 'latitude'
                    latitude.valid_min: -90.0
                    latitude.valid_max: 90.0
                    
                    latitude_bounds.long_name= 'latitude_bounds'
                    latitude_bounds.units= 'degrees_north'
                    
                    AOD550.long_name= 'aod 550 nm'
                    AOD550.units= '1'
                    AOD550.standard_name= 'aod_550'
                   
                    #values
                    pixel[:]=list(range(1,npixel+1))
                    retr[:]=1
                    longitude[:]=longitude_modis[np.logical_not(aod_modis.mask)]
                    longitude_bounds_aux=np.zeros((npixel,4))
                    pixel_resolution=0.1
                    
                    longitude_bounds_aux[:,0]=longitude_modis[np.logical_not(aod_modis.mask)] -pixel_resolution/2
                    longitude_bounds_aux[:,1]=longitude_modis[np.logical_not(aod_modis.mask)] +pixel_resolution/2
                    longitude_bounds_aux[:,2]=longitude_modis[np.logical_not(aod_modis.mask)] +pixel_resolution/2
                    longitude_bounds_aux[:,3]=longitude_modis[np.logical_not(aod_modis.mask)] -pixel_resolution/2
                    longitude_bounds[:]=longitude_bounds_aux
                    latitude[:]=latitude_modis[np.logical_not(aod_modis.mask)]
                    latitude_bounds_aux=np.zeros((npixel,4))
                    latitude_bounds_aux[:,0]=latitude_modis[np.logical_not(aod_modis.mask)] -pixel_resolution/2
                    latitude_bounds_aux[:,1]=latitude_modis[np.logical_not(aod_modis.mask)] -pixel_resolution/2
                    latitude_bounds_aux[:,2]=latitude_modis[np.logical_not(aod_modis.mask)] +pixel_resolution/2
                    latitude_bounds_aux[:,3]=latitude_modis[np.logical_not(aod_modis.mask)] +pixel_resolution/2
                    latitude_bounds[:]=latitude_bounds_aux
                    AOD550[:]=aod_modis[np.logical_not(aod_modis.mask)]
                    
                    
                            
                    
                    root_grp.close()
                    modis.close()
                
            
            
                    list_CSO['filename'].extend(['CSO/'+cso_file])
                    list_CSO['start_time'].extend([day_listing+hour_file])
                    list_CSO['end_time'].extend([day_listing+hour_file])
                    list_CSO['orbit'].extend([orbit])
                    a=a+1
                    print(a)
pd.DataFrame.from_dict(list_CSO).to_csv('/Users/santiago/Documents/MODIS/listing_modis.csv',sep=';',index=False)
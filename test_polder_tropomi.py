
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
from os import listdir
from os.path import isfile, join
import pandas as pd
#Create listing table
template=pd.read_csv('/Users/santiago/scratch/CAMS-61/WP6130/CSO-data/S5p/listing-NO2-CAMS.csv',delimiter=';')
list_CSO={name: [] for name in list(template.columns)}
for i in ['01','02','03','04','05','06','07','08','09','10','11','12']:
    mypath='/Users/santiago/Documents/POLDER/2008/'+i+'/'
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    onlyfiles.sort()
    cso_path=mypath+'/CSO/'
    
    
    for file_name in onlyfiles[1:]:
        if file_name[0:12]=='GRASP_POLDER':
            # Reading POLDER data
            polder=Dataset(mypath+file_name)
            aod_polder=polder.variables['AOD565'][0:350,1652:2148]
            latitude_polder=polder.variables['Latitude'][0:350,1652:2148]
            longitude_polder=polder.variables['Longitude'][0:350,1652:2148]
            latitude_polder=latitude_polder.reshape((latitude_polder.shape[0]*latitude_polder.shape[1]))
            longitude_polder=longitude_polder.reshape((longitude_polder.shape[0]*longitude_polder.shape[1]))
            aod_polder=aod_polder.reshape((aod_polder.shape[0]*aod_polder.shape[1]))
            
            # dimensions
            npixel = latitude_polder[np.logical_not(aod_polder.mask)].size
            if npixel<10:
                continue
            #Creating netcdf4
            cso_file='CSO_POLDER_'+file_name[16:24]+'.nc'
            root_grp = Dataset(cso_path+cso_file, 'w', format='NETCDF4')
            root_grp.description = 'POLDER observation following CSO estructure'
            
            
            
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
            AOD565 = root_grp.createVariable('AOD565', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
            
            
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
            
            AOD565.long_name= 'aod 565 nm'
            AOD565.units= '1'
            AOD565.standard_name= 'aod_565'
           
            #values
            pixel[:]=list(range(1,npixel+1))
            retr[:]=1
            longitude[:]=longitude_polder[np.logical_not(aod_polder.mask)]
            longitude_bounds_aux=np.zeros((npixel,4))
            longitude_bounds_aux[:,0]=longitude_polder[np.logical_not(aod_polder.mask)] -0.005
            longitude_bounds_aux[:,1]=longitude_polder[np.logical_not(aod_polder.mask)] +0.005
            longitude_bounds_aux[:,2]=longitude_polder[np.logical_not(aod_polder.mask)] +0.005
            longitude_bounds_aux[:,3]=longitude_polder[np.logical_not(aod_polder.mask)] -0.005
            longitude_bounds[:]=longitude_bounds_aux
            latitude[:]=latitude_polder[np.logical_not(aod_polder.mask)]
            latitude_bounds_aux=np.zeros((npixel,4))
            latitude_bounds_aux[:,0]=latitude_polder[np.logical_not(aod_polder.mask)] -0.005
            latitude_bounds_aux[:,1]=latitude_polder[np.logical_not(aod_polder.mask)] -0.005
            latitude_bounds_aux[:,2]=latitude_polder[np.logical_not(aod_polder.mask)] +0.005
            latitude_bounds_aux[:,3]=latitude_polder[np.logical_not(aod_polder.mask)] +0.005
            latitude_bounds[:]=latitude_bounds_aux
            AOD565[:]=aod_polder[np.logical_not(aod_polder.mask)]
            
            
                    
            
            root_grp.close()
            polder.close()
        
    
    csofiles = [f for f in listdir(cso_path) if isfile(join(cso_path, f))]
    csofiles.sort()
    for filename in csofiles:
        if filename[-3:]=='.nc':
            list_CSO['filename'].extend(['CSO/'+filename])
            list_CSO['start_time'].extend([filename[11:15]+'-'+filename[15:17]+'-'+filename[17:19]+'T13:30:00.00'])
            list_CSO['end_time'].extend([filename[11:15]+'-'+filename[15:17]+'-'+filename[17:19]+'T15:30:00.00'])
            list_CSO['orbit'].extend([filename[11:19]])
pd.DataFrame.from_dict(list_CSO).to_csv('/Users/santiago/Documents/POLDER/2008/listing_POLDER.csv',sep=';',index=False)
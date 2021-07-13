
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
valid_files=[]
tropomi=Dataset('/Users/santiago/scratch/CAMS-61/WP6130/CSO-data/S5p/RPRO/NO2/CAMS/2018/06/S5p_RPRO_NO2_03279.nc')
for i in ['01','02','03','04','05','06','07','08','09','10','11','12']:
    mypath='/Users/santiago/Documents/IASI/SO2/'+i+'/'
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    onlyfiles.sort()
    cso_path=mypath+'/CSO/'
    
    
    for file_name in onlyfiles[1:]:
        if file_name[0:11]=='IASI_METOPA':
            # Reading iasi data
            iasi=Dataset(mypath+file_name)
            so2_iasi=iasi.variables['SO2_interpolated'][:]
           
            
            latitude_iasi=iasi.variables['latitude'][:]
            longitude_iasi=iasi.variables['longitude'][:]
            
            
            # dimensions
            npixel = latitude_iasi[np.logical_and(longitude_iasi>=-15,np.logical_and(longitude_iasi<=35,np.logical_and(latitude_iasi>=35 ,latitude_iasi<=70)))].size
            if npixel<10:
                continue
            #Creating netcdf4
            
            cso_file='CSO_IASI_'+file_name[19:27]+'.nc'
            valid_files.append(cso_file)
            root_grp = Dataset(cso_path+cso_file, 'w', format='NETCDF4')
            root_grp.description = 'IASI observation following CSO estructure'
            
            
            
            root_grp.createDimension('pixel', npixel)
            root_grp.createDimension('layer', 34)
            root_grp.createDimension('layeri', 35)
            root_grp.createDimension('retr', 1)
            root_grp.createDimension('corner', 4)
            
            # variables
            pixel = root_grp.createVariable('pixel', 'i4', ('pixel',),fill_value=-2147483647)
            retr = root_grp.createVariable('retr', 'i4', ('retr',))
            longitude = root_grp.createVariable('longitude', 'f4', ('pixel',),fill_value=9.969209968386869e+36)
            longitude_bounds = root_grp.createVariable('longitude_bounds', 'f4', ('pixel','corner',),fill_value=9.969209968386869e+36)
            latitude = root_grp.createVariable('latitude', 'f4', ('pixel',),fill_value=9.969209968386869e+36)
            latitude_bounds = root_grp.createVariable('latitude_bounds', 'f4', ('pixel','corner',),fill_value=9.969209968386869e+36)
            SO2 = root_grp.createVariable('SO2', 'f4', ('pixel', 'retr',),fill_value=9.969209968386869e+36)
            pressure = root_grp.createVariable('pressure', 'f4', ('pixel', 'layeri',),fill_value=9.969209968386869e+36)
            Kernel = root_grp.createVariable('kernel', 'f4', ('pixel', 'layer','retr',),fill_value=9.969209968386869e+36)
            
            
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
            
            SO2.long_name= 'SO2 total column'
            SO2.units= 'mol m-2'
            SO2.standard_name= 'SO2'
            SO2.multiplication_factor_to_convert_to_molecules_percm2=6.02214179e+19
            
            pressure.long_name= 'half level pressures'
            pressure.units= 'Pa'
            pressure.standard_name= 'pressure'
            
            Kernel.long_name= 'Averaging kernel'
            Kernel.units= '1'
            Kernel.standard_name= 'kernel'
           
            #values
            resolution=0.12
            pixel[:]=list(range(1,npixel+1))
            retr[:]=1
            longitude[:]=longitude_iasi[np.logical_and(longitude_iasi>=-15,np.logical_and(longitude_iasi<=35,np.logical_and(latitude_iasi>=35 ,latitude_iasi<=70)))]
            longitude_bounds_aux=np.zeros((npixel,4))
            longitude_bounds_aux[:,0]=longitude_iasi[np.logical_and(longitude_iasi>=-15,np.logical_and(longitude_iasi<=35,np.logical_and(latitude_iasi>=35 ,latitude_iasi<=70)))] -resolution/2
            longitude_bounds_aux[:,1]=longitude_iasi[np.logical_and(longitude_iasi>=-15,np.logical_and(longitude_iasi<=35,np.logical_and(latitude_iasi>=35 ,latitude_iasi<=70)))] +resolution/2
            longitude_bounds_aux[:,2]=longitude_iasi[np.logical_and(longitude_iasi>=-15,np.logical_and(longitude_iasi<=35,np.logical_and(latitude_iasi>=35 ,latitude_iasi<=70)))] +resolution/2
            longitude_bounds_aux[:,3]=longitude_iasi[np.logical_and(longitude_iasi>=-15,np.logical_and(longitude_iasi<=35,np.logical_and(latitude_iasi>=35 ,latitude_iasi<=70)))] -resolution/2
            longitude_bounds[:]=longitude_bounds_aux
            latitude[:]=latitude_iasi[np.logical_and(longitude_iasi>=-15,np.logical_and(longitude_iasi<=35,np.logical_and(latitude_iasi>=35 ,latitude_iasi<=70)))]
            latitude_bounds_aux=np.zeros((npixel,4))
            latitude_bounds_aux[:,0]=latitude_iasi[np.logical_and(longitude_iasi>=-15,np.logical_and(longitude_iasi<=35,np.logical_and(latitude_iasi>=35 ,latitude_iasi<=70)))] -resolution/2
            latitude_bounds_aux[:,1]=latitude_iasi[np.logical_and(longitude_iasi>=-15,np.logical_and(longitude_iasi<=35,np.logical_and(latitude_iasi>=35 ,latitude_iasi<=70)))] -resolution/2
            latitude_bounds_aux[:,2]=latitude_iasi[np.logical_and(longitude_iasi>=-15,np.logical_and(longitude_iasi<=35,np.logical_and(latitude_iasi>=35 ,latitude_iasi<=70)))] +resolution/2
            latitude_bounds_aux[:,3]=latitude_iasi[np.logical_and(longitude_iasi>=-15,np.logical_and(longitude_iasi<=35,np.logical_and(latitude_iasi>=35 ,latitude_iasi<=70)))] +resolution/2
            latitude_bounds[:]=latitude_bounds_aux
            Kernel[:,:]=1
            SO2[:]=so2_iasi[np.logical_and(longitude_iasi>=-15,np.logical_and(longitude_iasi<=35,np.logical_and(latitude_iasi>=35 ,latitude_iasi<=70)))]
            
            
            pressure[:,:]=tropomi.variables['pressure'][1,:]
            
            
            
            
            root_grp.close()
            iasi.close()
        
    
    
for filename in valid_files:
        if filename[-3:]=='.nc':
            list_CSO['filename'].extend(['CSO/'+filename])
            list_CSO['start_time'].extend([filename[9:13]+'-'+filename[13:15]+'-'+filename[15:17]+'T09:30:00.00'])
            list_CSO['end_time'].extend([filename[9:13]+'-'+filename[13:15]+'-'+filename[15:17]+'T11:30:00.00'])
            list_CSO['orbit'].extend([filename[9:17]])
pd.DataFrame.from_dict(list_CSO).to_csv('/Users/santiago/Documents/IASI/SO2/listing_IASI.csv',sep=';',index=False)
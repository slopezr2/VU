
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



products=['SO2','NO2']
mypath='/Users/santiago/Documents/OMI/'
CSO_general_path='/Users/santiago/Documents/OMI/CSO/'
subfolders = [ f.path for f in scandir(mypath) if f.is_dir() ]
subfolders.sort()
a=0
for product in products:
    CSO_path=CSO_general_path+product
    OMI_Path=mypath+product
    template=pd.read_csv('/Users/santiago/scratch/CAMS-61/WP6130/CSO-data/S5p/listing-NO2-CAMS.csv',delimiter=';')
    list_CSO={name: [] for name in list(template.columns)}
    
    onlyfiles = [f for f in listdir(OMI_Path) if isfile(join(OMI_Path, f))]
    onlyfiles.sort()
    for file_name in onlyfiles:
        if file_name[0:8]=='OMI-Aura':
            # Reading omi data
            omi=Dataset(OMI_Path+'/'+file_name)
            if product=='SO2':
                specie_omi=omi.groups['HDFEOS'].groups['SWATHS'].groups['OMI Total Column Amount SO2'].groups['Data Fields'].variables['ColumnAmountSO2'][:]*4.4615e-4
                latitude_omi=omi.groups['HDFEOS'].groups['SWATHS'].groups['OMI Total Column Amount SO2'].groups['Geolocation Fields'].variables['Latitude'][:]
                latitude_bounds_omi = omi.groups['HDFEOS'].groups['SWATHS'].groups['OMI Total Column Amount SO2'].groups['Geolocation Fields'].variables['TiledCornerLatitude'][:,:]
                longitude_omi=omi.groups['HDFEOS'].groups['SWATHS'].groups['OMI Total Column Amount SO2'].groups['Geolocation Fields'].variables['Longitude'][:]
                longitude_bounds_omi = omi.groups['HDFEOS'].groups['SWATHS'].groups['OMI Total Column Amount SO2'].groups['Geolocation Fields'].variables['TiledCornerLongitude'][:,:]
                flaq_anomaly = omi.groups['HDFEOS'].groups['SWATHS'].groups['OMI Total Column Amount SO2'].groups['Data Fields'].variables['Flag_RowAnomaly'][:]
                kernel_omi = np.flip(omi.groups['HDFEOS'].groups['SWATHS'].groups['OMI Total Column Amount SO2'].groups['Data Fields'].variables['ScatteringWeight'][:-1])
                pressure_omi = np.flip(omi.groups['HDFEOS'].groups['SWATHS'].groups['OMI Total Column Amount SO2'].groups['Data Fields'].variables['LayerBottomPressure'][:])*100
            else:
                specie_omi=omi.groups['HDFEOS'].groups['SWATHS'].groups['ColumnAmountNO2'].groups['Data Fields'].variables['ColumnAmountSO2'][:]*4.4615e-4
                latitude_omi=omi.groups['HDFEOS'].groups['SWATHS'].groups['ColumnAmountNO2'].groups['Geolocation Fields'].variables['Latitude'][:]
                latitude_bounds_omi = omi.groups['HDFEOS'].groups['SWATHS'].groups['ColumnAmountNO2'].groups['Geolocation Fields'].variables['FoV75CornerLatitude'][:,:]
                longitude_omi=omi.groups['HDFEOS'].groups['SWATHS'].groups['ColumnAmountNO2'].groups['Geolocation Fields'].variables['Longitude'][:]
                longitude_bounds_omi = omi.groups['HDFEOS'].groups['SWATHS'].groups['ColumnAmountNO2'].groups['Geolocation Fields'].variables['FoV75CornerLongitude'][:,:]
                flaq_anomaly = omi.groups['HDFEOS'].groups['SWATHS'].groups['ColumnAmountNO2'].groups['Data Fields'].variables['Flag_RowAnomaly'][:]
                flaq_anomaly = flaq_anomaly & 1
                kernel_omi = omi.groups['HDFEOS'].groups['SWATHS'].groups['ColumnAmountNO2'].groups['Data Fields'].variables['ScatteringWeight'][:-1]
                pressure_omi = omi.groups['HDFEOS'].groups['SWATHS'].groups['ColumnAmountNO2'].groups['Data Fields'].variables['ScatteringWtPressure'][-1]*100
            
            specie_omi[specie_omi>2000] = np.nan
            specie_omi[specie_omi<0] = np.nan
            specie_omi[flaq_anomaly==1] = np.nan
            npixel = latitude_omi[np.logical_not(specie_omi.mask)].size
            
            if npixel<10:
                continue
            
            #Creating netcdf4
            day_year=file_name[18:22]+file_name[23:25]+file_name[25:27]
            
            day_listing=file_name[18:22]+'-'+file_name[23:25]+'-'+file_name[25:27]
            
            hour_file='T'+file_name[28:30]+':'+file_name[30:32]+':00.00'
            orbit=file_name[55:-8]
            
            cso_file='CSO_OMI_'+product+'_'+day_year+file_name[17:-4]+'.nc'
            root_grp = Dataset(CSO_path+'/'+cso_file, 'w', format='NETCDF4')
            root_grp.description = 'OMI '+product+' observation following CSO estructure'
            
            
            # dimensions
            n_layer = pressure_omi.shape[0]
            root_grp.createDimension('pixel', npixel)
            root_grp.createDimension('layer', n_layer-1)
            root_grp.createDimension('layeri', n_layer)
            root_grp.createDimension('retr', 1)
            root_grp.createDimension('corner', 4)
            # variables
            pixel = root_grp.createVariable('pixel', 'i4', ('pixel',),fill_value=-2147483647)
            retr = root_grp.createVariable('retr', 'i4', ('retr',))
            longitude = root_grp.createVariable('longitude', 'f4', ('pixel',),fill_value=3.4028235e+38)
            longitude_bounds = root_grp.createVariable('longitude_bounds', 'f4', ('pixel','corner',),fill_value=3.4028235e+38)
            latitude = root_grp.createVariable('latitude', 'f4', ('pixel',),fill_value=3.4028235e+38)
            latitude_bounds = root_grp.createVariable('latitude_bounds', 'f4', ('pixel','corner',),fill_value=3.4028235e+38)
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
            
            SO2.long_name= product+' total column'
            SO2.units= 'mol m-2'
            SO2.standard_name= product
            SO2.multiplication_factor_to_convert_to_molecules_percm2=6.02214179e+19
            
            pressure.long_name= 'half level pressures'
            pressure.units= 'Pa'
            pressure.standard_name= 'pressure'
            
            Kernel.long_name= 'Averaging kernel'
            Kernel.units= '1'
            Kernel.standard_name= 'kernel'
           
            #values
            pixel[:]=list(range(1,npixel+1))
            retr[:]=1
            valid_pixels = np.logical_not(specie_omi.mask)
            longitude[:]=longitude_omi[valid_pixels]
            longitude_bounds[:]=longitude_bounds_omi[valid_pixels,:]
            latitude[:]=latitude_omi[valid_pixels]
            latitude_bounds[:]=latitude_bounds_omi[valid_pixels,:]
            SO2[:]=specie_omi[valid_pixels]*4.4615e-4
            Kernel[:,:,0] = kernel_omi[valid_pixels,:]
            aux_pressure = np.ones((npixel,n_layer))
            aux_pressure[:,:] = pressure_omi
            pressure[:,:] = aux_pressure[valid_pixels,:]
            
                    
            
            root_grp.close()
            omi.close()
        
    
    
            list_CSO['filename'].extend([product+'/'+cso_file])
            list_CSO['start_time'].extend([day_listing+hour_file])
            list_CSO['end_time'].extend([day_listing+hour_file])
            list_CSO['orbit'].extend([orbit])
            a=a+1
            #print(a)
    pd.DataFrame.from_dict(list_CSO).to_csv('/Users/santiago/Documents/omi/CSO/listing_omi_'+product+'.csv',sep=';',index=False)
    
    
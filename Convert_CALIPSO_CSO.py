
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
mypath='/Users/santiago/Documents/calipso/2008/'
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
onlyfiles.sort()
cso_path=mypath+'/CSO/'

    
for file_name in onlyfiles[:]:
    if file_name[0:11]=='CAL_LID_L2_':
        # Reading calipso data
        calipso=Dataset(mypath+file_name)
        calipso.set_auto_mask(False)

        
        
        CAD=np.flip(calipso.variables['CAD_Score'][:,:,0],1)
        CAD=CAD[:,:344]

        qs=np.flip(calipso.variables['Extinction_Coefficient_Uncertainty_532'][:],1)
        qs=qs[:,:344]
        
        qc=np.flip(calipso.variables['Extinction_QC_Flag_532'][:],1)
        qc=qc[:,:344,0]
        
        
        AEC_calipso_aux=np.flip(calipso.variables['Extinction_Coefficient_532'][:],1)
        AEC_calipso_aux=AEC_calipso_aux[:,:344]#To take just the troposphere
        AEC_2=np.flip(calipso.variables['Extinction_Coefficient_532'][:],1)
        AEC_2=AEC_2[:,:344]
        #AEC_calipso_aux=np.ma.masked_array(AEC_calipso_aux, mask=((AEC_calipso_aux<=0)|(AEC_calipso_aux>1.25)))
        AEC_calipso_aux=np.ma.masked_array(AEC_calipso_aux, mask=((AEC_calipso_aux<=0)|(AEC_calipso_aux>1)| (CAD>-80)| (qs>AEC_calipso_aux*0.50) | ((qc!=1) & (qc!=0)) ))
        AEC_sa=AEC_calipso_aux
        AEC_calipso_aux=AEC_calipso_aux/1000
        
        latitude_calipso=calipso.variables['Latitude'][:]
        longitude_calipso=calipso.variables['Longitude'][:]
        
        pressure_calipso_aux=np.flip(calipso.variables['Pressure'][:],1)*100
        pressure_calipso_aux=pressure_calipso_aux[:,:344]#To take just the troposphere
        pressure_calipso_aux=np.ma.masked_array(pressure_calipso_aux, mask=((pressure_calipso_aux<=0)))
        
        df=pd.DataFrame(pressure_calipso_aux)
        pressure_calipso_aux=df.fillna(df.mean(axis=0)).to_numpy()
        
        inter_levels=17
        AEC_calipso=np.zeros([AEC_calipso_aux.shape[0],AEC_calipso_aux.shape[1]//inter_levels+1])
        pressure_calipso=np.zeros([pressure_calipso_aux.shape[0],pressure_calipso_aux.shape[1]//inter_levels+1])
        
        for i in range(AEC_calipso.shape[0]):
            AEC_calipso[i,:]= pd.Series(AEC_calipso_aux[i,:]).groupby(pd.Series(AEC_calipso_aux[i,:]).index//inter_levels).mean().to_numpy()
            pressure_calipso[i,:]= pd.Series(pressure_calipso_aux[i,:]).groupby(pd.Series(pressure_calipso_aux[i,:]).index//inter_levels).mean().to_numpy()
        
        fill_value=9.969209968386869e+36
        AEC_calipso[np.isnan(AEC_calipso)]=fill_value
        
        df=pd.DataFrame(pressure_calipso)
        pressure_calipso=df.fillna(df.mean(axis=0)).to_numpy()
        
        
        # dimensions
        npixel = latitude_calipso.shape[0]
        if npixel<10:
            continue
        #Creating netcdf4
        
        cso_file='CSO_calipso_'+file_name[35:54]+'.nc'
        valid_files.append(cso_file)
        print(cso_file)
        root_grp = Dataset(cso_path+cso_file, 'w', format='NETCDF4')
        root_grp.description = 'calipso observation following CSO estructure'
        
        
        
        root_grp.createDimension('pixel', npixel)
        root_grp.createDimension('layer', AEC_calipso.shape[1]-1)
        root_grp.createDimension('layeri', AEC_calipso.shape[1])
        root_grp.createDimension('retr', 1)
        root_grp.createDimension('corner', 4)
        
        # variables
        
        pixel = root_grp.createVariable('pixel', 'i4', ('pixel',),fill_value=-2147483647)
        retr = root_grp.createVariable('retr', 'i4', ('retr',))
        longitude = root_grp.createVariable('longitude', 'f4', ('pixel',),fill_value=fill_value)
        longitude_bounds = root_grp.createVariable('longitude_bounds', 'f4', ('pixel','corner',),fill_value=fill_value)
        latitude = root_grp.createVariable('latitude', 'f4', ('pixel',),fill_value=fill_value)
        latitude_bounds = root_grp.createVariable('latitude_bounds', 'f4', ('pixel','corner',),fill_value=fill_value)
        AEC = root_grp.createVariable('AEC', 'f4', ('pixel','layeri',),fill_value=fill_value)
        pressure = root_grp.createVariable('pressure', 'f4', ('pixel', 'layeri',),fill_value=fill_value)
        Kernel = root_grp.createVariable('kernel', 'f4', ('pixel', 'layer','retr',),fill_value=fill_value)
        
        
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
        
        AEC.long_name= 'AEC profile'
        AEC.units= 'm-1'
        AEC.standard_name= 'AEC'
        
        pressure.long_name= 'half level pressures'
        pressure.units= 'Pa'
        pressure.standard_name= 'pressure'
        
        Kernel.long_name= 'Averaging kernel'
        Kernel.units= '1'
        Kernel.standard_name= 'kernel'
       
        #values
    
        pixel[:]=list(range(1,npixel+1))
        retr[:]=1
        longitude[:]=longitude_calipso[:,1]
        longitude_bounds_aux=np.zeros((npixel,4))
        longitude_bounds_aux[:,0]=longitude_calipso[:,2]
        longitude_bounds_aux[:,1]=longitude_calipso[:,0] 
        longitude_bounds_aux[:,2]=longitude_calipso[:,0]
        longitude_bounds_aux[:,3]=longitude_calipso[:,2] 
        longitude_bounds[:]=longitude_bounds_aux
        latitude[:]=latitude_calipso[:,1]
        latitude_bounds_aux=np.zeros((npixel,4))
        latitude_bounds_aux[:,0]=latitude_calipso[:,2]
        latitude_bounds_aux[:,1]=latitude_calipso[:,2] 
        latitude_bounds_aux[:,2]=latitude_calipso[:,0] 
        latitude_bounds_aux[:,3]=latitude_calipso[:,0] 
        latitude_bounds[:]=latitude_bounds_aux
        Kernel[:,:]=1
        AEC[:,:]=AEC_calipso
        
        
        pressure[:,:]=pressure_calipso
        
        
        
        
        root_grp.close()
        calipso.close()
        
    
    
for filename in valid_files:
        if filename[-3:]=='.nc':
            list_CSO['filename'].extend(['CSO/'+filename])
            list_CSO['start_time'].extend([filename[12:16]+'-'+filename[17:19]+'-'+filename[20:22]+filename[22:25]+':'+filename[26:28]+':'+filename[29:31]+'.00'])
            list_CSO['end_time'].extend([filename[12:16]+'-'+filename[17:19]+'-'+filename[20:22]+filename[22:25]+':59:59.00'])
            list_CSO['orbit'].extend([filename[12:31]])
pd.DataFrame.from_dict(list_CSO).to_csv('/Users/santiago/Documents/CALIPSO/2008/listing_calipso.csv',sep=';',index=False)
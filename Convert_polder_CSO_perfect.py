import sys
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
#===Load user libraries====
sys.path.append('./libraries/')
from datamanager import DataManager_GRIDDED
from matplotlib import cm

#Create listing table
template=pd.read_csv('/Users/santiago/scratch/CAMS-61/WP6130/CSO-data/S5p/listing-NO2-CAMS.csv',delimiter=';')
list_CSO={name: [] for name in list(template.columns)}
mypath='/Users/santiago/Documents/POLDER/Perfect'
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
onlyfiles.sort()
cso_path=mypath+'/CSO/'
valid_files=[]
gridded_aqua=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_modis-aqua_2008.nc')
lat=gridded_aqua.nc.variables['lat'][:]
lon=gridded_aqua.nc.variables['lon'][:]
lat_perfect=[]
lon_perfect=[]
for i in lat:
    for j in lon:
        lat_perfect.append(i)
        lon_perfect.append(j)
lat_perfect=np.array(lat_perfect)
lon_perfect=np.array(lon_perfect)
aod_perfect=np.ones(lon_perfect.shape[0])

for dia in range(1,367):
    
    day, year = dia, 2008
    fecha=(pd.to_datetime(day-1, unit='D', origin=str(year)))
    print(fecha)
    fecha=fecha.strftime('%Y%m%d')
    fecha=fecha[0:8]
    print(fecha)
    npixel=aod_perfect.shape[0]
    cso_file='CSO_POLDER_Perfect_'+fecha+'.nc'
    valid_files.append(cso_file)
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
    AOD865 = root_grp.createVariable('AOD865', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
    AOD443 = root_grp.createVariable('AOD443', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
    AAOD565 = root_grp.createVariable('AAOD565', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
    AAOD865 = root_grp.createVariable('AAOD865', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
    AAOD443 = root_grp.createVariable('AAOD443', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
    SSA565 = root_grp.createVariable('SSA565', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
    SSA865 = root_grp.createVariable('SSA865', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
    SSA443 = root_grp.createVariable('SSA443', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
    AExp = root_grp.createVariable('AExp', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
    
    
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
    
    AOD443.long_name= 'aod 443 nm'
    AOD443.units= '1'
    AOD443.standard_name= 'aod_443'
    
    AOD865.long_name= 'aod 865 nm'
    AOD865.units= '1'
    AOD865.standard_name= 'aod_865'
    
    AAOD565.long_name= 'aod 565 nm'
    AAOD565.units= '1'
    AAOD565.standard_name= 'aod_565'
    
    AAOD443.long_name= 'aod 443 nm'
    AAOD443.units= '1'
    AAOD443.standard_name= 'aod_443'
    
    AAOD865.long_name= 'aod 865 nm'
    AAOD865.units= '1'
    AAOD865.standard_name= 'aod_865'
    
    SSA565.long_name= 'ssa 565 nm'
    SSA565.units= '1'
    SSA565.standard_name= 'ssa_565'
    
    SSA443.long_name= 'ssa 443 nm'
    SSA443.units= '1'
    SSA443.standard_name= 'ssa_443'
    
    SSA865.long_name= 'ssa 865 nm'
    SSA865.units= '1'
    SSA865.standard_name= 'ssa_865'
    
    
    AExp.long_name= 'Angstrom Exponent 443-865 nm'
    AExp.units= '1'
    AExp.standard_name= 'aexp'
   
    #values
    pixel[:]=list(range(1,npixel+1))
    retr[:]=1
    longitude[:]=lon_perfect
    longitude_bounds_aux=np.zeros((npixel,4))
    longitude_bounds_aux[:,0]=lon_perfect -0.005
    longitude_bounds_aux[:,1]=lon_perfect +0.005
    longitude_bounds_aux[:,2]=lon_perfect +0.005
    longitude_bounds_aux[:,3]=lon_perfect -0.005
    longitude_bounds[:]=longitude_bounds_aux
    latitude[:]=lat_perfect
    latitude_bounds_aux=np.zeros((npixel,4))
    latitude_bounds_aux[:,0]=lat_perfect -0.005
    latitude_bounds_aux[:,1]=lat_perfect -0.005
    latitude_bounds_aux[:,2]=lat_perfect +0.005
    latitude_bounds_aux[:,3]=lat_perfect +0.005
    latitude_bounds[:]=latitude_bounds_aux
    
    AOD565[:]=aod_perfect
    AOD443[:]=aod_perfect
    AOD865[:]=aod_perfect
    
    AAOD565[:]=aod_perfect
    AAOD443[:]=aod_perfect
    AAOD865[:]=aod_perfect
    
    SSA565[:]=aod_perfect
    SSA443[:]=aod_perfect
    SSA865[:]=aod_perfect
    
    AExp[:]=aod_perfect
    root_grp.close()

    
    
for filename in valid_files:
        if filename[-3:]=='.nc':
            list_CSO['filename'].extend(['CSO/'+filename])
            list_CSO['start_time'].extend([filename[19:23]+'-'+filename[23:25]+'-'+filename[25:27]+'T13:00:00.00'])
            list_CSO['end_time'].extend([filename[19:23]+'-'+filename[23:25]+'-'+filename[25:27]+'T13:30:00.00'])
            list_CSO['orbit'].extend([filename[19:27]])
pd.DataFrame.from_dict(list_CSO).to_csv('/Users/santiago/Documents/POLDER/Perfect/listing_POLDER_Perfect.csv',sep=';',index=False)
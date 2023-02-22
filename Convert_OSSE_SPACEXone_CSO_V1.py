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
import datetime
#===Load user libraries====
sys.path.append('./libraries/')
from datamanager import IF_NOT_OK
import pathlib

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
#Verificar que las variables existan en el archivo de CSO o el LE segun el caso.
#Ver si la variable dentro del for de CSO esta en el vector de nombres a cambiar, si si, consultar el indice para leer el nombre de la variable correspondiente en el vector de LE
#Hacer funcion para leer cada lat, lon y tomar el valor de LE.

def pixel_bounds(lat,lon,x_reso_meters,y_reso_meters):
    earth = 6378.137 #radius of the earth in kilometer
    m = (1 / ((2 * np.pi / 360) * earth)) / 1000 # 1 meter in degree
    bounds_lat=np.zeros(2)
    bounds_lon=np.zeros(2)
    bounds_lat[0] = lat - (y_reso_meters * m);
    bounds_lat[1] = lat + (y_reso_meters * m);

    bounds_lon[0] = lon - (x_reso_meters * m) / np.cos(lat * (np.pi / 180));
    bounds_lon[1] = lon + (x_reso_meters * m) / np.cos(lat * (np.pi / 180));
    return bounds_lat,bounds_lon
#=== Experiment Name ===
experiment_name = 'Test_SPACEOnex_V1'


#===LE AOD values===
#=A file for all the year=
LE_file = Dataset('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_aod_2008.nc')

#=== Sampling satellite ===
#=CSO file of the satellite used as sampling==
satellite_path =  '/Users/santiago/Documents/SPEXone_OSSEs/data'

#===Template CSO file ===
cso_template = Dataset('/Users/santiago/Documents/POLDER/Model/complete/CSO_POLDER_20080101.nc')

#===Variables satellite==
cso_variables = ['AOD565','AExp']

#===Variables LE==
LE_variables = ['aod_563nm','angstrom_polder']
Dim_variables = [3,3,4] #To improve the performance of the script
#==Observation error==
observation_error = [0.1,0.15,0.23,0.16,0.1]
#==========================================================================================================================================================================================================================
status = 0
error_message = ''

if len(cso_variables) != len(LE_variables):
    status = 1
    error_message = 'The number of variables on CSO file and LE file does not agree.'

IF_NOT_OK(status, error_message)

#===Timestamp 2008===
timestamp_2008 = datetime.datetime.timestamp(datetime.datetime.strptime('2008-1-1 00:00:00', '%Y-%m-%d %H:%M:%S'))

#Create listing table
template=pd.read_csv('/Users/santiago/scratch/CAMS-61/WP6130/CSO-data/S5p/listing-NO2-CAMS.csv',delimiter=';')
list_CSO={name: [] for name in list(template.columns)}


onlyfiles = [f for f in listdir(satellite_path) if isfile(join(satellite_path, f))]
onlyfiles.sort()
osse_path = '/Users/santiago/Documents/OSSE/'+experiment_name

pathlib.Path(osse_path).mkdir(parents=True,exist_ok=True)
pathlib.Path(osse_path+'/CSO').mkdir(parents=True,exist_ok=True)

valid_files = []
for file in onlyfiles:
    if file[0:3]!='Cut':
        continue
    src_file = join(satellite_path, file)
    new_name = 'CSO/CSO_OSSE_'+experiment_name+'_2008'+file[-7:] 
    trg_file = join(osse_path,new_name)
    idx_time = (datetime.datetime.strptime('2008'+file[-7:-3]+'13', "%Y%m%d%H").timestamp() - timestamp_2008)
    if idx_time<LE_file.variables['time'][0] or idx_time>LE_file.variables['time'][-1]:
        continue
    idx_time = idx_time/3600 - LE_file.variables['time'][0]/3600
    valid_files.append(new_name)
    src = Dataset(src_file)
    trg = Dataset(trg_file, mode='w')
    print('Creating the file: '+trg_file)
    # Create the dimensions of the file   
    dim_sat = src.variables['lat'][:].shape[0]
    for name, dim in cso_template.dimensions.items():
        if name=='pixel':
            trg.createDimension(name, dim_sat if not dim.isunlimited() else None) 
        else:
            trg.createDimension(name, len(dim) if not dim.isunlimited() else None)
    
    # Copy the global attributes
    trg.setncatts({a:cso_template.getncattr(a) for a in cso_template.ncattrs()})
    trg.description = 'SPEXOne observations following CSO estructure from OSSEs'
    # Create the variables in the file
    for name, var in cso_template.variables.items():
        trg.createVariable(name, var.dtype, var.dimensions) 
        # Copy the variable attributes
        trg.variables[name].setncatts({a:var.getncattr(a) for a in var.ncattrs()})
        # Put in zero all variables except non_zero_vars
        trg.variables[name][:] = 0
    
    aux_variables = {key: [] for key in cso_variables}
    aux_latitude_bounds = np.zeros((dim_sat,4))
    aux_longitude_bounds = np.zeros((dim_sat,4))
    for i in range(src.variables['lat'][:].shape[0]):
        if 'lat' in LE_file.variables.keys():
            idx_lat = find_nearest(LE_file.variables['lat'][:], src.variables['lat'][i])
        else:
            idx_lat = find_nearest(LE_file.variables['latitude'][:], src.variables['lat'][i])
        if 'lon' in LE_file.variables.keys():
            idx_lon = find_nearest(LE_file.variables['lon'][:], src.variables['lon'][i])
        else:
            idx_lon = find_nearest(LE_file.variables['longitude'][:], src.variables['lon'][i])
        i_variable = 0
        lat_bounds,lon_bounds = pixel_bounds(src.variables['lat'][i],src.variables['lon'][i],5400,4600)
        aux_latitude_bounds[i] = [lat_bounds[0],lat_bounds[0],lat_bounds[1],lat_bounds[1]]
        aux_longitude_bounds[i] = [lon_bounds[0],lon_bounds[1],lon_bounds[0],lon_bounds[1]]
        for cso_variable,LE_variable in zip(cso_variables,LE_variables):
            status = 0
            if not(LE_variable in LE_file.variables.keys()):
                status = 1
                error_message = 'The LE variable is not in the LE file.'
            IF_NOT_OK(status,error_message)
            if Dim_variables[i_variable]==3:
                LE_value = LE_file.variables[LE_variable][idx_time,idx_lat,idx_lon]
                LE_AOD = LE_file.variables['aod_563nm'][idx_time,idx_lat,idx_lon]
            else:
                LE_value = np.mean(LE_file.variables[LE_variable][idx_time,:,idx_lat,idx_lon])
                LE_AOD = np.mean(LE_file.variables['aod_563nm'][idx_time,:,idx_lat,idx_lon])
            
            if LE_AOD<=0.4:
                rand_error = observation_error[0]*np.random.randn(1)
            elif LE_AOD>0.4 and LE_AOD<=0.6:
                rand_error = observation_error[1]*np.random.randn(1)
            elif LE_AOD>0.6 and LE_AOD<=0.8:
                rand_error = observation_error[2]*np.random.randn(1)
            elif LE_AOD>0.8 and LE_AOD<=1.2:
                rand_error = observation_error[3]*np.random.randn(1)
            else:
                rand_error = observation_error[4]*np.random.randn(1)
                
            aux_variables[cso_variable].append(max(0,LE_value+rand_error))
            i_variable = i_variable + 1
    for cso_variable in cso_variables:       
        trg.variables[cso_variable][:] = aux_variables[cso_variable][:]
    trg.variables['latitude_bounds'][:] = aux_latitude_bounds
    trg.variables['longitude_bounds'][:] = aux_longitude_bounds
    trg.variables['latitude'][:] = src.variables['lat'][:]
    trg.variables['longitude'][:] = src.variables['lon'][:]
    # Save the file
    trg.close()
    src.close()
    
    
    
   
    
for filename in valid_files:
        if filename[-3:]=='.nc':
            list_CSO['filename'].extend([filename])
            list_CSO['start_time'].extend([filename[-11:-7]+'-'+filename[-7:-5]+'-'+filename[-5:-3]+'T13:00:00.00'])
            list_CSO['end_time'].extend([filename[-11:-7]+'-'+filename[-7:-5]+'-'+filename[-5:-3]+'T13:30:00.00'])
            list_CSO['orbit'].extend([filename[-11:-3]])
pd.DataFrame.from_dict(list_CSO).to_csv(osse_path+'/listing_'+experiment_name+'.csv',sep=';',index=False)
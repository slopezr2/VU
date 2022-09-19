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



#=== Experiment Name ===
experiment_name = 'Natural_Run_200805_POLDER'


#===LE AOD values===
#=A file for all the year=
LE_file = Dataset(po)

#=== Sampling satellite ===
#=CSO file of the satellite used as sampling==
satellite_path =  '/Users/santiago/Documents/POLDER/Model/complete'


#===Variables satellite==
cso_variables = ['AOD565','AExp','SSA565']

#===Variables LE==
LE_variables = ['aod_563nm','angstrom_polder','ssa_550nm']
Dim_variables = [3,3,4] #To improve the performance of the script
#==Observation error==
observation_error = [0.002,0.002,0.002]
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
    src_file = join(satellite_path, file)
    new_name = 'CSO/CSO_OSSE_'+experiment_name+'_'+file[-11:] 
    trg_file = join(osse_path,new_name)
    idx_time = (datetime.datetime.strptime(file[-11:-3]+'13', "%Y%m%d%H").timestamp() - timestamp_2008)
    if idx_time<LE_file.variables['time'][0] or idx_time>LE_file.variables['time'][-1]:
        continue
    idx_time = idx_time/3600 - LE_file.variables['time'][0]/3600
    valid_files.append(new_name)
    src = Dataset(src_file)
    trg = Dataset(trg_file, mode='w')
    print('Creating the file: '+trg_file)
    # Create the dimensions of the file   
    for name, dim in src.dimensions.items():
        trg.createDimension(name, len(dim) if not dim.isunlimited() else None) 
    # Copy the global attributes
    trg.setncatts({a:src.getncattr(a) for a in src.ncattrs()})
     
    # Create the variables in the file
    for name, var in src.variables.items():
        trg.createVariable(name, var.dtype, var.dimensions) 
        # Copy the variable attributes
        trg.variables[name].setncatts({a:var.getncattr(a) for a in var.ncattrs()})
        # Put in zero all variables except non_zero_vars
        trg.variables[name][:] = src.variables[name][:]
    
    aux_variables = {key: [] for key in cso_variables}
    for i in range(src.variables['latitude'][:].shape[0]):
        if 'lat' in LE_file.variables.keys():
            idx_lat = find_nearest(LE_file.variables['lat'][:], src.variables['latitude'][i])
        else:
            idx_lat = find_nearest(LE_file.variables['latitude'][:], src.variables['latitude'][i])
        if 'lon' in LE_file.variables.keys():
            idx_lon = find_nearest(LE_file.variables['lon'][:], src.variables['longitude'][i])
        else:
            idx_lon = find_nearest(LE_file.variables['longitude'][:], src.variables['longitude'][i])
        i_variable = 0
        for cso_variable,LE_variable in zip(cso_variables,LE_variables):
            status = 0
            error_message = ''
            if not(cso_variable in src.variables.keys()):
                status = 1
                error_message = 'The CSO variable is not in the CSO files.'
            
            if not(LE_variable in LE_file.variables.keys()):
                status = 1
                error_message = 'The LE variable is not in the LE file.'
            IF_NOT_OK(status,error_message)
            if Dim_variables[i_variable]==3:
                LE_value = LE_file.variables[LE_variable][idx_time,idx_lat,idx_lon]
            else:
                LE_value = np.mean(LE_file.variables[LE_variable][idx_time,:,idx_lat,idx_lon])
            
            rand_error = observation_error[i_variable]*np.random.randn(1)
            aux_variables[cso_variable].append(max(0,LE_value+rand_error))
            i_variable = i_variable + 1
    for cso_variable in cso_variables:       
        trg.variables[cso_variable][:] = aux_variables[cso_variable][:]
            
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
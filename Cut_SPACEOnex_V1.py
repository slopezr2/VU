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

#=== Sampling satellite ===
#=CSO file of the satellite used as sampling==
satellite_path =  '/Users/santiago/Documents/SPEXone_OSSEs/data'
#===LE AOD values===
#=A file for all the year=
LE_file = Dataset('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_aod_2008.nc')


#==========================================================================================================================================================================================================================
status = 0
error_message = ''



onlyfiles = [f for f in listdir(satellite_path) if isfile(join(satellite_path, f))]
onlyfiles.sort()

valid_files = []
for file in onlyfiles:
    if file[0]=='.' or file[0:3]=='Cut':
        continue
    src_file = join(satellite_path, file)
    new_name = 'Cut_'+file 
    trg_file = join(satellite_path,new_name)
    src = Dataset(src_file)
    trg = Dataset(trg_file, mode='w')
    print('Creating the file: '+trg_file)

    # Get dimension of pixel into de domain
    domain = (src.variables['lat'][:]>=LE_file.variables['lat'][:].min()) & (src.variables['lat'][:]<=LE_file.variables['lat'][:].max()) & (src.variables['lon'][:]>=LE_file.variables['lon'][:].min()) & (src.variables['lon'][:]<=LE_file.variables['lon'][:].max())
    lat_domain = src.variables['lat'][domain]
    lon_domain = src.variables['lon'][domain]
    var_domain = src.variables['cloud_mask'][domain]
    dim_domain=len(lat_domain)
    
    # Create the dimensions of the file   
    for name, dim in src.dimensions.items():
        trg.createDimension(name, dim_domain if not dim.isunlimited() else None) 
    # Copy the global attributes
    trg.setncatts({a:src.getncattr(a) for a in src.ncattrs()})
     
    # Create the variables in the file
    for name, var in src.variables.items():
        trg.createVariable(name, var.dtype, var.dimensions) 
        # Copy the variable attributes
        trg.variables[name].setncatts({a:var.getncattr(a) for a in var.ncattrs()})
        # Put in zero all variables except non_zero_vars
    
    trg.variables['lat'][:] = lat_domain
    trg.variables['lon'][:] = lon_domain
    trg.variables['cloud_mask'][:] = var_domain
    
   
    # Save the file
    trg.close()
    src.close()
    
    
    

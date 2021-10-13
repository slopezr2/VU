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


def zeros_variables(src_file, trg_file,non_zero_vars):
    """create_file_from_source('in.nc', 'out.nc') with setting all variable to zero except for the non_zero_vars """
    src = Dataset(src_file)
    trg = Dataset(trg_file, mode='w')
    # Variables that are always required to be a proper LE output
    neededvars = ['longitude','latitude','longitude_bnds','latitude_bnds','longitude_crnr','latitude_crnr','level','time','altitude', 'T', 'Q', 'hp', 'halt']
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
        if name in non_zero_vars or name in neededvars:
          trg.variables[name][:] = src.variables[name][:]
        else:
          trg.variables[name][:] = 0.*src.variables[name][:]
    # Save the file
    trg.close()
    src.close()

mypath='/Users/santiago/Documents/LOTOS-EUROS/conc_zero/'
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
onlyfiles.sort()

for filename in onlyfiles:
    if filename[-3:] == '.nc':
        zeros_variables(mypath+filename, mypath+filename[0:-3]+'_nao3.nc',['no3a_f', 'no3a_c'])

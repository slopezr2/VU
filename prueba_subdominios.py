import sys
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import cm
sys.path.append('./libraries/')
from datamanager import graph_map

RdBu = cm.get_cmap('RdBu', 256)
domain=np.loadtxt('./dc_bien.txt')
dc=Dataset('/Users/santiago/Downloads/LE_Noise_Thanos_V1_dc_20080502_xi01a.nc')

graph_map(variable=domain,lat=dc.variables['latitude'][:],lon=dc.variables['longitude'][:],cmap=RdBu.reversed(),date=False,title='Perturbation factor',vmax=2, dpi=300,save=True)

splitx=4
splity=4
process=splitx*splity
gnx,gny=domain.shape
nx=int(gnx/splitx)
ny=int(gny/splity)

subdomains={}
subdomains2={}
i=1
j=1
for iproces in range(process):
    subdomains[iproces]=np.empty((gnx,gny))
    subdomains2[iproces]=np.empty((nx,ny))
    subdomains[iproces][:]=np.nan
    subdomains2[iproces][:]=np.nan
    subdomains[iproces][(i-1)*nx:i*nx,(j-1)*ny:j*ny]=domain[(i-1)*nx:i*nx,(j-1)*ny:j*ny]
    subdomains2[iproces]=domain[(i-1)*nx:i*nx,(j-1)*ny:j*ny]
    i=i+1
    if i>splitx:
        i=1
        j=j+1
    
    graph_map(variable=subdomains[iproces],lat=dc.variables['latitude'][:],lon=dc.variables['longitude'][:],cmap=RdBu.reversed(),date=False,title='Perturbation factor',vmax=2, dpi=300,save=True)

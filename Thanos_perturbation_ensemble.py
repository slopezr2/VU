import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.stats as sts
from netCDF4 import Dataset
from matplotlib import cm
sys.path.append('./libraries/')
from datamanager import graph_map

RdBu = cm.get_cmap('RdBu', 256)

def getmeanofneighbors(matrix, i, j,r):
    region = matrix[max(0, i-r) : i+r+1,
                    max(0, j-r) : j+r+1]
    #print(region)
    #return region
    return (np.sum(region) - matrix[i, j])/(region.size-1) # Sum the region and subtract center



nx=100
ny=110
r=2
repeat=5
new_mean=1
new_std=0.65
ensembles=300


ensemble_dc = np.zeros((ensembles,nx,ny))
for ens in range(ensembles):
    domain=1+np.random.randn(nx,ny)*10
    domain_backup2=domain.copy()

    for k in range(repeat):
        domain_backup=domain.copy()
        for i in range(domain.shape[0]):
            for j in range(domain.shape[1]):
                domain[i,j]=getmeanofneighbors(domain_backup, i, j,r)
                
    
    domain=np.exp(domain)

    
    domain=new_std*((domain-domain.mean())/domain.std())+new_mean
    ensemble_dc[ens,:,:] = domain
    
   
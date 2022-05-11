import sys
#===Load user libraries====
sys.path.append('./libraries/')
from datamanager import DataManager_GRIDDED
from datamanager import DataManager_LE
from datamanager import graph_map
from netCDF4 import Dataset
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
from cartopy import config
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime,timedelta
from datetime import date as date_datetime
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm
from os import listdir
from os.path import isfile, join
import matplotlib.dates as mdates
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import cm
from matplotlib.colors import ListedColormap

#===Color definitions==
upper = cm.RdYlGn_r(np.arange(256))
lower = np.ones((1,4))
#lower = cm.viridis(np.arange(256))
value_union=1
for i in range(3):
  lower[-value_union:,i] = np.linspace(lower[-value_union,i], upper[0,i], value_union)


# combine parts of colormap
RdYlGn_w = np.vstack(( lower, upper ))

# convert to matplotlib colormap
RdYlGn_w = ListedColormap(RdYlGn_w, name='myColorMap', N=RdYlGn_w.shape[0])


upper = cm.Spectral_r(np.arange(256))
lower = np.ones((1,4))
#lower = cm.viridis(np.arange(256))
value_union=1
for i in range(3):
  lower[-value_union:,i] = np.linspace(lower[-value_union,i], upper[0,i], value_union)


# combine parts of colormap
Spectral_w = np.vstack(( lower, upper ))

# convert to matplotlib colormap
Spectral_w = ListedColormap(Spectral_w, name='myColorMap', N=Spectral_w.shape[0])

viridis = cm.get_cmap('viridis', 256)
orange=cm.get_cmap('Oranges', 256)
bwr=cm.get_cmap('bwr', 256)
seismic=cm.get_cmap('seismic', 256)
spectral=cm.get_cmap('Spectral', 256)
plasma=cm.get_cmap('RdYlGn', 256)


#====Parameters===
save=True


#====Concentrations===
tpm10_Xb=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm10_200805_S_T_all_0001_xb.nc')
tpm10_Xa_S_T_all_0001=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm10_200805_S_T_all_0001_xa.nc')


#DC_names=['S_T_all_0001','S_T_aer_0001']
DC_names=['S_T_all_0001']
i=0
DC={}
Emis={}
nnoises={}

for DC_name in DC_names:
    mypath='/Users/santiago/Documents/LE_outputs/LETKF/DC_'+DC_name
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f)) and f[0]!='.']
    onlyfiles.sort()
    DC_aux=Dataset(mypath+'/'+onlyfiles[0])
    nnoises[DC_name]=DC_aux.variables['dc'].shape[2]
    DC[DC_name]=DC_aux.variables['dc'][:,0,:,:,:]
    Emis[DC_name]=Dataset('/Users/santiago/Documents/LE_outputs/LETKF/LE_m_emis_'+DC_name+'_200805.nc')
    for file in onlyfiles[1:]:
       DC_aux=Dataset(mypath+'/'+file) 
       DC[DC_name]=np.concatenate((DC[DC_name],DC_aux.variables['dc'][:,0,:,:,:]))

noise_names={}
noise_names['S_T_all_0001']=[['SOx_emis'],['PPM_emis'],['BC_emis'],['POM_emis'],['dust_emis'],['SS_emis']]
noise_names['S_T_aer_0001']=[['aer_emis']]

emis_names={}
emis_names['S_T_all_0001']=[['so4a_f','so2'],['ppm_f','ppm_c'],['ec_f','ec_c'],['pom_f','pom_c'],['tdust'],['tss']]
emis_names['S_T_aer_0001']=[['so4a_f','so2'],['ppm_f','ppm_c'],['ec_f','ec_c'],['pom_f','pom_c'],['tdust'],['tss']]


#===Graph Noise====
lat=DC_aux.variables['latitude'][:]
lon=DC_aux.variables['longitude'][:]

for DC_name in DC_names:
    inoise=0
    for noise in noise_names[DC_name]:
        graph_map(variable=np.mean(DC[DC_name][:,inoise,:,:],axis=0),lat=lat,lon=lon,vmin=0,vmax=10,vcenter=0,cmap=Spectral_w,date=False,save=save,title=noise[0]+' DC',dpi=300,grid=True)
        emission=np.zeros(Emis[DC_name].variables[emis_names[DC_name][0][0]][:,0,:,:].shape)
        for var in emis_names[DC_name][inoise]:
            aux=Emis[DC_name].variables[var][:,0,:,:]
            emission=emission+aux
        new_emission=emission*DC[DC_name][:,inoise,:,:]
        vmax=emission.mean()*1e9*5
        graph_map(variable=np.mean(emission,axis=0)*1e9,lat=lat,lon=lon,vmin=0,vmax=vmax,vcenter=0,cmap=RdYlGn_w,date=False,save=save,title=noise[0]+r' Xb emissions $\mu$g m$^{-2}$ $s^{-1}$',dpi=300,grid=True,save_title=noise[0]+r' Xb emissions')
        graph_map(variable=np.mean(new_emission,axis=0)*1e9,lat=lat,lon=lon,vmin=0,vmax=vmax,vcenter=0,cmap=RdYlGn_w,date=False,save=save,title=noise[0]+r' Xa emissions $\mu$g m$^{-2}$ $s^{-1}$',dpi=300,grid=True,save_title=noise[0]+r' Xa emissions')
        
        inoise=inoise+1
        
        
#==== Graph concentrations maps===
#burden_tpm10_Xb=np.zeros(tpm10_Xa_S_T_all_0001.nc.variables['tpm10'][:,0,:,:].shape)
#burden_tpm10_Xa_S_T_all_0001=np.zeros(tpm10_Xa_S_T_all_0001.nc.variables['tpm10'][:,0,:,:].shape)
for layer in range(12):
    tpm10_Xa_S_T_all_0001.graph_map(variable='tpm10',level=layer,biascorr=1e9,vmin=0,vmax=100,n_levels_plot=20,cmap=orange,type_graph='mean',date=False,save=save,title='LE-Xa PM$_{10}$ layer '+str(layer) +' 2008-05-01_14 $\mu$g m$^{-3}$',dpi=300,grid=True,save_title='LE_Xa_PM10_layer_'+str(layer))
    tpm10_Xb.graph_map(variable='tpm10',level=layer,biascorr=1e9,vmin=0,vmax=100,n_levels_plot=20,cmap=orange,type_graph='mean',date=False,save=save,title='LE-Xb PM$_{10}$ layer '+str(layer) +' 2008-05-01_14 $\mu$g m$^{-3}$',dpi=300,grid=True,save_title='LE_Xb_PM10_layer_'+str(layer))
    graph_map(variable=np.mean(tpm10_Xa_S_T_all_0001.nc.variables['tpm10'][:,layer,:,:]/tpm10_Xb.nc.variables['tpm10'][:,layer,:,:],axis=0),lat=lat,lon=lon,vmin=0,vmax=10,vcenter=1,norm2=True,cmap=seismic,date=False,save=save,title='LE-Xa/LE-Xb PM$_{10}$ layer '+str(layer) +' 2008-05-01_14',dpi=300,grid=True,save_title='LE_Xa_LE_Xb_PM10'+str(layer))
        
    
    #burden_tpm10_Xb=burden_tpm10_Xb+tpm10_Xb.nc.variables['tpm10'][:,layer,:,:]
    #burden_tpm10_Xa_S_T_all_0001=burden_tpm10_Xa_S_T_all_0001+tpm10_Xa_S_T_all_0001.nc.variables['tpm10'][:,layer,:,:]

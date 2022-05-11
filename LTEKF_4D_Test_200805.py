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


viridis = cm.get_cmap('viridis', 256)
orange=cm.get_cmap('Oranges', 256)
bwr=cm.get_cmap('bwr', 256)
#===Create Color bar====
def make_Ramp( ramp_colors ): 
    from colour import Color
    from matplotlib.colors import LinearSegmentedColormap

    color_ramp = LinearSegmentedColormap.from_list( 'my_list', [ Color( c1 ).rgb for c1 in ramp_colors ] )
    return color_ramp


custom_ramp = make_Ramp( ["#ffffff","#aec8e0","#9dafb2","#e2b85f","#f0ae46","#ea9842","#d36f2b","#a04617","#6c423d"] ) 


#====Save===
save=True

#==Period for larger simulations==
ini_period=0 #02-May
end_period=12 #14-May
#===Without R factor 2008-05-01 --- 2008-07-31====
gridded_Xb=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_polder_200805_07_S_T_all_xb.nc')
gridded_Xa_S_F_all=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_polder_200805_07_S_F_all_xa.nc')
gridded_Xa_S_F_aer=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_polder_200805_07_S_F_aer_xa.nc')
gridded_Xa_S_T_all=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_polder_200805_07_S_T_all_xa.nc')
gridded_Xa_S_T_aer=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_polder_200805_07_S_T_aer_xa.nc')
gridded_Xa_Just_emiss=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_polder_200805_Just_emiss_V1_xa.nc')


#===With R factor, shorter period 2008-05-01 --- 2008-05-14=== 
gridded_Xb_R=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_polder_200805_S_T_all_001_xb.nc')
gridded_Xa_S_T_all_001=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_polder_200805_S_T_all_001_xa.nc')
gridded_Xa_S_T_all_0001=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_polder_200805_S_T_all_0001_xa.nc')
gridded_Xa_S_T_aer_001=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_polder_200805_S_T_aer_001_xa.nc')
gridded_Xa_S_T_aer_0001=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_polder_200805_S_T_aer_0001_xa.nc')

# LETKF_PM10_Xb=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm10_200805_07_S_F_all_xb.nc')
# LETKF_PM10_Xa_S_F_all=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm10_200805_07_S_F_all_xa.nc')
# LETKF_PM10_Xa_S_F_aer=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm10_200805_07_S_F_aer_xa.nc')
# LETKF_PM10_Xa_S_T_all=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm10_200805_07_S_T_all_xa.nc')
# LETKF_PM10_Xa_S_T_aer=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm10_200805_07_S_T_aer_xa.nc')

# LETKF_PM25_Xb=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm25_200805_07_S_F_all_xb.nc')
# LETKF_PM25_Xa_S_F_all=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm25_200805_07_S_F_all_xa.nc')
# LETKF_PM25_Xa_S_F_aer=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm25_200805_07_S_F_aer_xa.nc')
# LETKF_PM25_Xa_S_T_all=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm25_200805_07_S_T_all_xa.nc')
# LETKF_PM25_Xa_S_T_aer=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm25_200805_07_S_T_aer_xa.nc')


#LETKF_aod_Xa=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LE_m_aod563_200805_xa.nc')
#LETKF_aod_Xb=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LE_m_aod563_200805_xb.nc')

#dc_Xa=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LE_CSO_V5_dc_20080509_xa.nc')

#====AOD Polder maps===
gridded_Xb.graph_map(variable='yr',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=save,title='POLDER AOD565 2008-05_07',dpi=300,grid=True,ini_period=ini_period,end_period=end_period)
gridded_Xb.graph_map(variable='ys',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=save,title='LE-Xb AOD565 2008-05_07',dpi=300,grid=True,ini_period=ini_period,end_period=end_period)
#gridded_Xa_S_F_all.graph_map(variable='ys',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=save,title='LE-Xa S_F_all AOD565 2008-05_07',dpi=300,grid=True)
#gridded_Xa_S_F_aer.graph_map(variable='ys',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=save,title='LE-Xa S_F_aer AOD565 2008-05_07',dpi=300,grid=True)
#gridded_Xa_S_T_all.graph_map(variable='ys',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=save,title='LE-Xa S_T_all AOD565 2008-05_07',dpi=300,grid=True)
#gridded_Xa_S_T_aer.graph_map(variable='ys',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=save,title='LE-Xa S_T_aer AOD565 2008-05_07',dpi=300,grid=True)

#gridded_Xb_R.graph_map(variable='yr',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=save,title='POLDER AOD565 2008-05-01_14',dpi=300,grid=True)
#gridded_Xb_R.graph_map(variable='ys',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=save,title='LE-Xb AOD565 2008-05-01_14',dpi=300,grid=True)
#gridded_Xa_S_T_all_001.graph_map(variable='ys',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=save,title='LE-Xa S_T_all_001 AOD565 2008-05-01_14',dpi=300,grid=True)
gridded_Xa_S_T_all_0001.graph_map(variable='ys',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=save,title='LE-Xa Conc+Emiss AOD565 2008-05-01_14',dpi=300,grid=True)
#gridded_Xa_S_T_aer_001.graph_map(variable='ys',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=save,title='LE-Xa S_T_aer_001 AOD565 2008-05-01_14',dpi=300,grid=True)
#gridded_Xa_S_T_aer_0001.graph_map(variable='ys',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=save,title='LE-Xa S_T_aer_0001 AOD565 2008-05-01_14',dpi=300,grid=True)

gridded_Xa_Just_emiss.graph_map(variable='ys',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=save,title='LE-Xa Emiss AOD565 2008-05_07',dpi=300,grid=True,ini_period=ini_period,end_period=end_period)


#===Calculate MAPE and MAE====
mape_Xb=100*(-np.mean(gridded_Xb.nc.variables['yr'][ini_period:end_period,:,:,0],axis=0)+np.mean(gridded_Xb.nc.variables['ys'][ini_period:end_period,:,:,0],axis=0))/np.mean(gridded_Xb.nc.variables['yr'][ini_period:end_period,:,:,0],axis=0)
mape_Xa_S_F_all=100*(-np.mean(gridded_Xb.nc.variables['yr'][ini_period:end_period,:,:,0],axis=0)+np.mean(gridded_Xa_S_F_all.nc.variables['ys'][ini_period:end_period,:,:,0],axis=0))/np.mean(gridded_Xb.nc.variables['yr'][ini_period:end_period,:,:,0],axis=0)
mape_Xa_S_F_aer=100*(-np.mean(gridded_Xb.nc.variables['yr'][ini_period:end_period,:,:,0],axis=0)+np.mean(gridded_Xa_S_F_aer.nc.variables['ys'][ini_period:end_period,:,:,0],axis=0))/np.mean(gridded_Xb.nc.variables['yr'][ini_period:end_period,:,:,0],axis=0)
mape_Xa_S_T_all=100*(-np.mean(gridded_Xb.nc.variables['yr'][ini_period:end_period,:,:,0],axis=0)+np.mean(gridded_Xa_S_T_all.nc.variables['ys'][ini_period:end_period,:,:,0],axis=0))/np.mean(gridded_Xb.nc.variables['yr'][ini_period:end_period,:,:,0],axis=0)
mape_Xa_S_T_aer=100*(-np.mean(gridded_Xb.nc.variables['yr'][ini_period:end_period,:,:,0],axis=0)+np.mean(gridded_Xa_S_T_aer.nc.variables['ys'][ini_period:end_period,:,:,0],axis=0))/np.mean(gridded_Xb.nc.variables['yr'][ini_period:end_period,:,:,0],axis=0)

mae_Xb=(-np.mean(gridded_Xb.nc.variables['yr'][ini_period:end_period,:,:,0],axis=0)+np.mean(gridded_Xb.nc.variables['ys'][ini_period:end_period,:,:,0],axis=0))
mae_Xa_S_F_all=(-np.mean(gridded_Xb.nc.variables['yr'][ini_period:end_period,:,:,0],axis=0)+np.mean(gridded_Xa_S_F_all.nc.variables['ys'][ini_period:end_period,:,:,0],axis=0))
mae_Xa_S_F_aer=(-np.mean(gridded_Xb.nc.variables['yr'][ini_period:end_period,:,:,0],axis=0)+np.mean(gridded_Xa_S_F_aer.nc.variables['ys'][ini_period:end_period,:,:,0],axis=0))
mae_Xa_S_T_all=(-np.mean(gridded_Xb.nc.variables['yr'][ini_period:end_period,:,:,0],axis=0)+np.mean(gridded_Xa_S_T_all.nc.variables['ys'][ini_period:end_period,:,:,0],axis=0))
mae_Xa_S_T_aer=(-np.mean(gridded_Xb.nc.variables['yr'][ini_period:end_period,:,:,0],axis=0)+np.mean(gridded_Xa_S_T_aer.nc.variables['ys'][ini_period:end_period,:,:,0],axis=0))


mape_Xb_R=100*(-np.mean(gridded_Xb_R.nc.variables['yr'][:],axis=0)+np.mean(gridded_Xb_R.nc.variables['ys'][:],axis=0))/np.mean(gridded_Xb_R.nc.variables['yr'][:],axis=0)
mape_Xa_S_T_all_001=100*(-np.mean(gridded_Xa_S_T_all_001.nc.variables['yr'][:],axis=0)+np.mean(gridded_Xa_S_T_all_001.nc.variables['ys'][:],axis=0))/np.mean(gridded_Xa_S_T_all_001.nc.variables['yr'][:],axis=0)
mape_Xa_S_T_all_0001=100*(-np.mean(gridded_Xa_S_T_all_0001.nc.variables['yr'][:],axis=0)+np.mean(gridded_Xa_S_T_all_0001.nc.variables['ys'][:],axis=0))/np.mean(gridded_Xa_S_T_all_0001.nc.variables['yr'][:],axis=0)
mape_Xa_S_T_aer_001=100*(-np.mean(gridded_Xa_S_T_aer_001.nc.variables['yr'][:],axis=0)+np.mean(gridded_Xa_S_T_aer_001.nc.variables['ys'][:],axis=0))/np.mean(gridded_Xa_S_T_aer_001.nc.variables['yr'][:],axis=0)
mape_Xa_S_T_aer_0001=100*(-np.mean(gridded_Xa_S_T_aer_0001.nc.variables['yr'][:],axis=0)+np.mean(gridded_Xa_S_T_aer_0001.nc.variables['ys'][:],axis=0))/np.mean(gridded_Xa_S_T_aer_0001.nc.variables['yr'][:],axis=0)


mae_Xb_R=(-np.mean(gridded_Xb_R.nc.variables['yr'][:],axis=0)+np.mean(gridded_Xb_R.nc.variables['ys'][:],axis=0))
mae_Xa_S_T_all_001=(-np.mean(gridded_Xa_S_T_all_001.nc.variables['yr'][:],axis=0)+np.mean(gridded_Xa_S_T_all_001.nc.variables['ys'][:],axis=0))
mae_Xa_S_T_all_0001=(-np.mean(gridded_Xa_S_T_all_0001.nc.variables['yr'][:],axis=0)+np.mean(gridded_Xa_S_T_all_0001.nc.variables['ys'][:],axis=0))
mae_Xa_S_T_aer_001=(-np.mean(gridded_Xa_S_T_aer_001.nc.variables['yr'][:],axis=0)+np.mean(gridded_Xa_S_T_aer_001.nc.variables['ys'][:],axis=0))
mae_Xa_S_T_aer_0001=(-np.mean(gridded_Xa_S_T_aer_0001.nc.variables['yr'][:],axis=0)+np.mean(gridded_Xa_S_T_aer_0001.nc.variables['ys'][:],axis=0))



#===Graph MAPE and MAE maps====
lat=gridded_Xb_R.nc.variables['latitude'][:]
lon=gridded_Xb_R.nc.variables['longitude'][:]


graph_map(variable=mape_Xb[:,:],lat=lat,lon=lon,vmin=-100,vmax=100,norm2=True,vcenter=0,cmap=bwr,date=False,save=save,title='MAPE Xb mean='+str(mape_Xb.mean())[0:5],dpi=300,grid=True)
#graph_map(variable=mape_Xa_S_F_all[:,:],lat=lat,lon=lon,vmin=-100,vmax=100,norm2=True,vcenter=0,cmap=bwr,date=False,save=save,title='MAPE Xa_S_F_all mean='+str(mape_Xa_S_F_all.mean())[0:5],dpi=300,grid=True)
#graph_map(variable=mape_Xa_S_F_aer[:,:],lat=lat,lon=lon,vmin=-100,vmax=100,norm2=True,vcenter=0,cmap=bwr,date=False,save=save,title='MAPE Xa_S_F_aer mean='+str(mape_Xa_S_F_aer.mean())[0:5],dpi=300,grid=True)
#graph_map(variable=mape_Xa_S_T_all[:,:],lat=lat,lon=lon,vmin=-100,vmax=100,norm2=True,vcenter=0,cmap=bwr,date=False,save=save,title='MAPE Xa_S_T_all mean='+str(mape_Xa_S_T_all.mean())[0:5],dpi=300,grid=True)
#graph_map(variable=mape_Xa_S_T_aer[:,:],lat=lat,lon=lon,vmin=-100,vmax=100,norm2=True,vcenter=0,cmap=bwr,date=False,save=save,title='MAPE Xa_S_T_aer mean='+str(mape_Xa_S_T_aer.mean())[0:5],dpi=300,grid=True)

#graph_map(variable=mae_Xb[:,:],lat=lat,lon=lon,vmin=-2,vmax=2,norm2=True,vcenter=0,cmap=bwr,date=False,save=save,title='MAE Xb mean='+str(mae_Xb.mean())[0:5],dpi=300,grid=True)
#graph_map(variable=mae_Xa_S_F_all[:,:],lat=lat,lon=lon,vmin=-2,vmax=2,norm2=True,vcenter=0,cmap=bwr,date=False,save=save,title='MAE Xa_S_F_all mean='+str(mae_Xa_S_F_all.mean())[0:5],dpi=300,grid=True)
#graph_map(variable=mae_Xa_S_F_aer[:,:],lat=lat,lon=lon,vmin=-2,vmax=2,norm2=True,vcenter=0,cmap=bwr,date=False,save=save,title='MAE Xa_S_F_aer mean='+str(mae_Xa_S_F_aer.mean())[0:5],dpi=300,grid=True)
#graph_map(variable=mae_Xa_S_T_all[:,:],lat=lat,lon=lon,vmin=-2,vmax=2,norm2=True,vcenter=0,cmap=bwr,date=False,save=save,title='MAE Xa_S_T_all mean='+str(mae_Xa_S_T_all.mean())[0:5],dpi=300,grid=True)
#graph_map(variable=mae_Xa_S_T_aer[:,:],lat=lat,lon=lon,vmin=-2,vmax=2,norm2=True,vcenter=0,cmap=bwr,date=False,save=save,title='MAE Xa_S_T_aer mean='+str(mae_Xa_S_T_aer.mean())[0:5],dpi=300,grid=True)


#graph_map(variable=mape_Xb_R[:,:,0],lat=lat,lon=lon,vmin=-100,vmax=100,norm2=True,vcenter=0,cmap=bwr,date=False,save=save,title='MAPE Xb mean='+str(mape_Xb_R.mean())[0:5],dpi=300,grid=True)
#graph_map(variable=mape_Xa_S_T_all_001[:,:,0],lat=lat,lon=lon,vmin=-100,vmax=100,norm2=True,vcenter=0,cmap=bwr,date=False,save=save,title='MAPE Xa_S_T_all_001 mean='+str(mape_Xa_S_T_all_001.mean())[0:5],dpi=300,grid=True)
graph_map(variable=mape_Xa_S_T_all_0001[:,:,0],lat=lat,lon=lon,vmin=-100,vmax=100,norm2=True,vcenter=0,cmap=bwr,date=False,save=save,title='MAPE Xa Conc+Emiss mean='+str(mape_Xa_S_T_all_0001.mean())[0:5],dpi=300,grid=True)
#graph_map(variable=mape_Xa_S_T_aer_001[:,:,0],lat=lat,lon=lon,vmin=-100,vmax=100,norm2=True,vcenter=0,cmap=bwr,date=False,save=save,title='MAPE Xa_S_T_aer_001 mean='+str(mape_Xa_S_T_aer_001.mean())[0:5],dpi=300,grid=True)
#graph_map(variable=mape_Xa_S_T_aer_0001[:,:,0],lat=lat,lon=lon,vmin=-100,vmax=100,norm2=True,vcenter=0,cmap=bwr,date=False,save=save,title='MAPE Xa_S_T_aer_0001 mean='+str(mape_Xa_S_T_aer_0001.mean())[0:5],dpi=300,grid=True)


#graph_map(variable=mae_Xb_R[:,:,0],lat=lat,lon=lon,vmin=-1,vmax=1,norm2=True,vcenter=0,cmap=bwr,date=False,save=save,title='MAE Xb mean='+str(mae_Xb_R.mean())[0:5],dpi=300,grid=True)
#graph_map(variable=mae_Xa_S_T_all_001[:,:,0],lat=lat,lon=lon,vmin=-1,vmax=1,norm2=True,vcenter=0,cmap=bwr,date=False,save=save,title='MAE Xa_S_T_all_001 mean='+str(mae_Xa_S_T_all_001.mean())[0:5],dpi=300,grid=True)
graph_map(variable=mae_Xa_S_T_all_0001[:,:,0],lat=lat,lon=lon,vmin=-1,vmax=1,norm2=True,vcenter=0,cmap=bwr,date=False,save=save,title='MAE Xa Conc+Emiss mean='+str(mae_Xa_S_T_all_0001.mean())[0:5],dpi=300,grid=True)
#graph_map(variable=mae_Xa_S_T_aer_001[:,:,0],lat=lat,lon=lon,vmin=-1,vmax=1,norm2=True,vcenter=0,cmap=bwr,date=False,save=save,title='MAE Xa_S_T_aer_001 mean='+str(mae_Xa_S_T_aer_001.mean())[0:5],dpi=300,grid=True)
#graph_map(variable=mae_Xa_S_T_aer_0001[:,:,0],lat=lat,lon=lon,vmin=-1,vmax=1,norm2=True,vcenter=0,cmap=bwr,date=False,save=save,title='MAE Xa_S_T_aer_0001 mean='+str(mae_Xa_S_T_aer_0001.mean())[0:5],dpi=300,grid=True)

 
# LETKF_PM10_Xb.graph_map(variable='tpm10',level=1,biascorr=1e9,vmin=0,vmax=30,n_levels_plot=20,cmap=orange,type_graph='mean',date=False,save=save,title='LE-Xb PM10 2008-05_07',dpi=300,grid=True)
# LETKF_PM10_Xa_S_F_all.graph_map(variable='tpm10',level=1,biascorr=1e9,vmin=0,vmax=30,n_levels_plot=20,cmap=orange,type_graph='mean',date=False,save=save,title='LE-Xa S_F_all PM10 2008-05_07',dpi=300,grid=True)
# LETKF_PM10_Xa_S_F_aer.graph_map(variable='tpm10',level=1,biascorr=1e9,vmin=0,vmax=30,n_levels_plot=20,cmap=orange,type_graph='mean',date=False,save=save,title='LE-Xa S_F_aer PM10 2008-05_07',dpi=300,grid=True)
# LETKF_PM10_Xa_S_T_all.graph_map(variable='tpm10',level=1,biascorr=1e9,vmin=0,vmax=30,n_levels_plot=20,cmap=orange,type_graph='mean',date=False,save=save,title='LE-Xa S_T_all PM10 2008-05_07',dpi=300,grid=True)
# LETKF_PM10_Xa_S_T_aer.graph_map(variable='tpm10',level=1,biascorr=1e9,vmin=0,vmax=30,n_levels_plot=20,cmap=orange,type_graph='mean',date=False,save=save,title='LE-Xa S_T_aer PM10 2008-05_07',dpi=300,grid=True)


# LETKF_PM25_Xb.graph_map(variable='tpm25',level=1,biascorr=1e9,vmin=0,vmax=20,n_levels_plot=20,cmap=orange,type_graph='mean',date=False,save=save,title='LE-Xb PM25 2008-05_07',dpi=300,grid=True)
# LETKF_PM25_Xa_S_F_all.graph_map(variable='tpm25',level=1,biascorr=1e9,vmin=0,vmax=20,n_levels_plot=20,cmap=orange,type_graph='mean',date=False,save=save,title='LE-Xa S_F_all PM25 2008-05_07',dpi=300,grid=True)
# LETKF_PM25_Xa_S_F_aer.graph_map(variable='tpm25',level=1,biascorr=1e9,vmin=0,vmax=20,n_levels_plot=20,cmap=orange,type_graph='mean',date=False,save=save,title='LE-Xa S_F_aer PM25 2008-05_07',dpi=300,grid=True)     
# LETKF_PM25_Xa_S_T_all.graph_map(variable='tpm25',level=1,biascorr=1e9,vmin=0,vmax=20,n_levels_plot=20,cmap=orange,type_graph='mean',date=False,save=save,title='LE-Xa S_T_all PM25 2008-05_07',dpi=300,grid=True)
# LETKF_PM25_Xa_S_T_aer.graph_map(variable='tpm25',level=1,biascorr=1e9,vmin=0,vmax=20,n_levels_plot=20,cmap=orange,type_graph='mean',date=False,save=save,title='LE-Xa S_T_aer PM25 2008-05_07',dpi=300,grid=True)
                      

#LETKF_aod_Xa.graph_map(variable='aod_563nm',level=1,biascorr=1,vmin=0,vmax=0.5,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=save,title='LE-Xa AOD565 All 2008-05',dpi=300,grid=True)
#LETKF_aod_Xb.graph_map(variable='aod_563nm',level=1,biascorr=1,vmin=0,vmax=0.5,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=save,title='LE-Xb AOD565 All 2008-05',dpi=300,grid=True)

 
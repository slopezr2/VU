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

#===Create Color bar====
def make_Ramp( ramp_colors ): 
    from colour import Color
    from matplotlib.colors import LinearSegmentedColormap

    color_ramp = LinearSegmentedColormap.from_list( 'my_list', [ Color( c1 ).rgb for c1 in ramp_colors ] )
    return color_ramp


custom_ramp = make_Ramp( ["#ffffff","#aec8e0","#9dafb2","#e2b85f","#f0ae46","#ea9842","#d36f2b","#a04617","#6c423d"] ) 

gridded_Xb_R_1=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_polder_200805_R_1_xb.nc')
gridded_Xa_R_1=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_polder_200805_R_1_xa.nc')

gridded_Xb_R_001=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_polder_200805_R_001_xb.nc')
gridded_Xa_R_001=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_polder_200805_R_001_xa.nc')

gridded_Xb_R_0001=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_polder_200805_R_0001_xb.nc')
gridded_Xa_R_0001=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_polder_200805_R_0001_xa.nc')


LETKF_PM10_Xb_R_1=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm10_200805_R_1_xb.nc')
LETKF_PM10_Xa_R_1=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm10_200805_R_1_xa.nc')

LETKF_PM10_Xb_R_001=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm10_200805_R_001_xb.nc')
LETKF_PM10_Xa_R_001=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm10_200805_R_001_xa.nc')

LETKF_PM10_Xb_R_0001=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm10_200805_R_0001_xb.nc')
LETKF_PM10_Xa_R_0001=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm10_200805_R_0001_xa.nc')



LETKF_PM25_Xb_R_1=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm25_200805_R_1_xb.nc')
LETKF_PM25_Xa_R_1=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm25_200805_R_1_xa.nc')

LETKF_PM25_Xb_R_001=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm25_200805_R_001_xb.nc')
LETKF_PM25_Xa_R_001=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm25_200805_R_001_xa.nc')

LETKF_PM25_Xb_R_0001=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm25_200805_R_0001_xb.nc')
LETKF_PM25_Xa_R_0001=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm25_200805_R_0001_xa.nc')


#LETKF_aod_Xa=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LE_m_aod563_200805_xa.nc')
#LETKF_aod_Xb=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LE_m_aod563_200805_xb.nc')

#dc_Xa=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LE_CSO_V5_dc_20080509_xa.nc')

gridded_Xb_R_1.graph_map(variable='yr',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='POLDER AOD565 2008-05-01_03',dpi=300,grid=True)
gridded_Xb_R_1.graph_map(variable='ys',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='LE-Xb AOD565 2008-05-01_03',dpi=300,grid=True)
gridded_Xa_R_1.graph_map(variable='ys',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='LE-Xa R_1 AOD565 2008-05-01_03',dpi=300,grid=True)
gridded_Xa_R_001.graph_map(variable='ys',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='LE-Xa R_001 AOD565 2008-05-01_03',dpi=300,grid=True)
gridded_Xa_R_0001.graph_map(variable='ys',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='LE-Xa R_0001 AOD565 2008-05-01_03',dpi=300,grid=True)

 
LETKF_PM10_Xb_R_1.graph_map(variable='tpm10',level=1,biascorr=1e9,vmin=0,vmax=30,n_levels_plot=20,cmap=orange,type_graph='mean',date=False,save=True,title='LE-Xb PM10 2008-05-01_03',dpi=300,grid=True)
LETKF_PM10_Xa_R_1.graph_map(variable='tpm10',level=1,biascorr=1e9,vmin=0,vmax=30,n_levels_plot=20,cmap=orange,type_graph='mean',date=False,save=True,title='LE-Xa R_1 PM10 2008-05-01_03',dpi=300,grid=True)
LETKF_PM10_Xa_R_001.graph_map(variable='tpm10',level=1,biascorr=1e9,vmin=0,vmax=30,n_levels_plot=20,cmap=orange,type_graph='mean',date=False,save=True,title='LE-Xa R_001 PM10 2008-05-01_03',dpi=300,grid=True)
LETKF_PM10_Xa_R_0001.graph_map(variable='tpm10',level=1,biascorr=1e9,vmin=0,vmax=30,n_levels_plot=20,cmap=orange,type_graph='mean',date=False,save=True,title='LE-Xa R_0001 PM10 2008-05-01_03',dpi=300,grid=True)


LETKF_PM25_Xb_R_1.graph_map(variable='tpm25',level=1,biascorr=1e9,vmin=0,vmax=20,n_levels_plot=20,cmap=orange,type_graph='mean',date=False,save=True,title='LE-Xb PM25 2008-05-01_03',dpi=300,grid=True)
LETKF_PM25_Xa_R_1.graph_map(variable='tpm25',level=1,biascorr=1e9,vmin=0,vmax=20,n_levels_plot=20,cmap=orange,type_graph='mean',date=False,save=True,title='LE-Xa R_1 PM25 2008-05-01_03',dpi=300,grid=True)
LETKF_PM25_Xa_R_001.graph_map(variable='tpm25',level=1,biascorr=1e9,vmin=0,vmax=20,n_levels_plot=20,cmap=orange,type_graph='mean',date=False,save=True,title='LE-Xa R_001 PM25 2008-05-01_03',dpi=300,grid=True)     
LETKF_PM25_Xa_R_0001.graph_map(variable='tpm25',level=1,biascorr=1e9,vmin=0,vmax=20,n_levels_plot=20,cmap=orange,type_graph='mean',date=False,save=True,title='LE-Xa R_0001 PM25 2008-05-01_03',dpi=300,grid=True)
               




#LETKF_aod_Xa.graph_map(variable='aod_563nm',level=1,biascorr=1,vmin=0,vmax=0.5,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='LE-Xa AOD565 All 2008-05',dpi=300,grid=True)
#LETKF_aod_Xb.graph_map(variable='aod_563nm',level=1,biascorr=1,vmin=0,vmax=0.5,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='LE-Xb AOD565 All 2008-05',dpi=300,grid=True)

 
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

gridded_Xb=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_polder_200805_xb.nc')
gridded_Xa=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_polder_200805_xa.nc')
gridded_Xb_post=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_post_m_polder_200805_xb.nc')

LETKF_PM10_Xa=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm10_200805_xa.nc')
LETKF_PM10_Xb=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm10_200805_xb.nc')

LETKF_PM25_Xa=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm25_200805_xa.nc')
LETKF_PM25_Xb=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_tpm25_200805_xb.nc')

LETKF_aod_Xa=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LE_m_aod563_200805_xa.nc')
LETKF_aod_Xb=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LE_m_aod563_200805_xb.nc')

dc_Xa=DataManager_LE('/Users/santiago/Documents/LE_outputs/LETKF/LE_CSO_V5_dc_20080509_xa.nc')

gridded_Xb.graph_map(variable='ys',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='LE-XB AOD565 2008-05',dpi=300,grid=True)
gridded_Xb.graph_map(variable='yr',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='POLDER AOD565 2008-05',dpi=300,grid=True)

gridded_Xb_post.graph_map(variable='ys',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='LE-XB post AOD565 2008-05',dpi=300,grid=True)
gridded_Xb_post.graph_map(variable='yr',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='POLDER post AOD565 2008-05',dpi=300,grid=True)

gridded_Xa.graph_map(variable='ys',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='LE-Xa AOD565 2008-05',dpi=300,grid=True)
gridded_Xa.graph_map(variable='ys',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='LE-Xa post AOD565 2008-05',dpi=300,grid=True)
gridded_Xa.graph_map(variable='yr',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='POLDER AOD565 2008-05 Xa',dpi=300,grid=True)
        
LETKF_PM10_Xa.graph_map(variable='tpm10',level=1,biascorr=1e9,vmin=0,vmax=30,n_levels_plot=20,cmap=orange,type_graph='mean',date=False,save=True,title='PM10 2008-05 Xa',dpi=300,grid=True)
LETKF_PM10_Xb.graph_map(variable='tpm10',level=1,biascorr=1e9,vmin=0,vmax=30,n_levels_plot=20,cmap=orange,type_graph='mean',date=False,save=True,title='PM10 2008-05 Xb',dpi=300,grid=True)

LETKF_PM25_Xa.graph_map(variable='tpm25',level=1,biascorr=1e9,vmin=0,vmax=20,n_levels_plot=20,cmap=orange,type_graph='mean',date=False,save=True,title='PM25 2008-05 Xa',dpi=300,grid=True)
LETKF_PM25_Xb.graph_map(variable='tpm25',level=1,biascorr=1e9,vmin=0,vmax=20,n_levels_plot=20,cmap=orange,type_graph='mean',date=False,save=True,title='PM25 2008-05 Xb',dpi=300,grid=True)
             
LETKF_aod_Xa.graph_map(variable='aod_563nm',level=1,biascorr=1,vmin=0,vmax=0.5,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='LE-Xa AOD565 All 2008-05',dpi=300,grid=True)
LETKF_aod_Xb.graph_map(variable='aod_563nm',level=1,biascorr=1,vmin=0,vmax=0.5,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='LE-Xb AOD565 All 2008-05',dpi=300,grid=True)

 
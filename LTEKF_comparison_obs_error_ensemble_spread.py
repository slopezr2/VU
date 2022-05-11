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
RdYlBu_r=cm.get_cmap('RdYlBu_r', 256)
#===Create Color bar====
def make_Ramp( ramp_colors ): 
    from colour import Color
    from matplotlib.colors import LinearSegmentedColormap

    color_ramp = LinearSegmentedColormap.from_list( 'my_list', [ Color( c1 ).rgb for c1 in ramp_colors ] )
    return color_ramp


custom_ramp = make_Ramp( ["#ffffff","#aec8e0","#9dafb2","#e2b85f","#f0ae46","#ea9842","#d36f2b","#a04617","#6c423d"] ) 


#====Save===
save=False

#==Period for larger simulations==
ini_period=0 #02-May
#end_period=12 #14-May
end_period=-1 #14-May

#===Without R factor 2008-05-01 --- 2008-07-31====
gridded_Xb=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_polder_200805_07_S_T_all_xb.nc')
gridded_Xa_Just_emiss=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_polder_200805_Just_emiss_V1_xa.nc')
gridded_Sa_Just_emiss=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/LETKF/LETKF_m_polder_200805_Just_emiss_V1_sa.nc')


#====Observation error and Ensemble Spread===
gridded_Xa_Just_emiss.graph_map(variable='obs_error',biascorr=1,vmin=0,vmax=0.2,n_levels_plot=20,cmap=RdYlBu_r,type_graph='mean',date=False,save=save,title='POLDER Obs error 2008-05_07',dpi=300,grid=True,ini_period=ini_period,end_period=end_period)
gridded_Sa_Just_emiss.graph_map(variable='ys',biascorr=1,vmin=0,vmax=0.2,n_levels_plot=20,cmap=RdYlBu_r,type_graph='mean',date=False,save=save,title='LE-Sa Emiss  2008-05_07',dpi=300,grid=True,ini_period=ini_period,end_period=end_period)

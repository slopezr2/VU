import sys

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
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


#===Load user libraries====
sys.path.append('./libraries/')
from datamanager import DataManager_GRIDDED
from datamanager import DataManager_LE
from matplotlib import cm


viridis = cm.get_cmap('viridis', 256)
redgreen = cm.get_cmap('RdYlGn', 256)
redblue = cm.get_cmap('RdYlBu', 256)
redgreen=redgreen.reversed()
redblue=redblue.reversed()

gridded_nh3=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_iasi-nh3_2008.nc')
gridded_so2=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_iasi-so2_2008.nc')


gridded_nh3.graph_map(variable='yr',facecolor=(1,1,1),biascorr=17.03052*1000,vmin=0,vmax=5,n_levels_plot=10,cmap=redgreen,type_graph='mean',date=False,save=True,title='IASI NH$_{3}$ total column',ocean=True)
gridded_nh3.graph_map(variable='ys',facecolor=(1,1,1),biascorr=17.03052*1000,vmin=0,vmax=5,n_levels_plot=10,cmap=redgreen,type_graph='mean',date=False,save=True,title='LE NH$_{3}$ total column',ocean=True)

gridded_so2.graph_map(variable='yr',facecolor=(1,1,1),biascorr=64.0638*1000,vmin=0,vmax=100,n_levels_plot=10,cmap=redblue,type_graph='mean',date=False,save=True,title='IASI SO$_{2}$ total column',ocean=False)
gridded_so2.graph_map(variable='ys',facecolor=(1,1,1),biascorr=64.0638*1000,vmin=0,vmax=100,n_levels_plot=10,cmap=redblue,type_graph='mean',date=False,save=True,title='LE SO$_{2}$ total column',ocean=False)

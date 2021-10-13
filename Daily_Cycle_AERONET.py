import sys

from netCDF4 import Dataset
import pandas as pd
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
import mpl_scatter_density
from matplotlib import dates as d
import datetime as dt
import matplotlib
#===Load user libraries====
sys.path.append('./libraries/')
from datamanager import DataManager_GRIDDED
from matplotlib import cm

gridded_500=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_aeronet_500_2008.nc')

timestamp_2008 = datetime.timestamp(datetime.strptime('2008-1-1 00:00:00', '%Y-%m-%d %H:%M:%S'))

lat=gridded_500.nc.variables['lat'][:]
lon=gridded_500.nc.variables['lon'][:]

fechas={}
fechas['Fecha']=[]
fechas['Aeronet_station_1']=[]
fechas['LE_station_1']=[]
fechas['Aeronet_station_2']=[]
fechas['LE_station_2']=[]
fechas['Aeronet_Minsk']=[]
fechas['LE_Minsk']=[]
for i in range(gridded_500.nc.variables['ys'].shape[0]):
    fechas['Fecha'].append(datetime.fromtimestamp(gridded_500.nc.variables['time'][i]+timestamp_2008))
    aux=gridded_500.nc.variables['ys'][i,:,:]
    aux_pos=np.argwhere(aux)
    flaq_1=0
    flaq_2=0
    flaq_minsk=0
    for j in range(len(aux_pos)):
        if aux_pos[j][0]==6 and aux_pos[j][1]==35:
            flaq_1=1
            fechas['Aeronet_station_1'].append(gridded_500.nc.variables['yr'][i,aux_pos[j][0],aux_pos[j][1]])
            fechas['LE_station_1'].append(gridded_500.nc.variables['ys'][i,aux_pos[j][0],aux_pos[j][1]])
        elif aux_pos[j][0]==10 and aux_pos[j][1]==60:
            flaq_2=1
            fechas['Aeronet_station_2'].append(gridded_500.nc.variables['yr'][i,aux_pos[j][0],aux_pos[j][1]])
            fechas['LE_station_2'].append(gridded_500.nc.variables['ys'][i,aux_pos[j][0],aux_pos[j][1]])
        elif aux_pos[j][0]==75 and aux_pos[j][1]==85 and gridded_500.nc.variables['yr'][i,aux_pos[j][0],aux_pos[j][1]]>0:
            flaq_minsk=1
            fechas['Aeronet_Minsk'].append(gridded_500.nc.variables['yr'][i,aux_pos[j][0],aux_pos[j][1]])
            fechas['LE_Minsk'].append(gridded_500.nc.variables['ys'][i,aux_pos[j][0],aux_pos[j][1]])    
            
    if flaq_1==0:
        fechas['Aeronet_station_1'].append(np.nan)
        fechas['LE_station_1'].append(np.nan)
    if flaq_2==0:
        fechas['Aeronet_station_2'].append(np.nan)
        fechas['LE_station_2'].append(np.nan)
    if flaq_minsk==0:
        fechas['Aeronet_Minsk'].append(np.nan)
        fechas['LE_Minsk'].append(np.nan)        
aeronet_dataframe=pd.DataFrame.from_dict(fechas,dtype=float)
aeronet_dataframe.set_index('Fecha',inplace=True)
aeronet_dataframe['hour'] =aeronet_dataframe.index.hour
aeronet_dataframe['month'] =aeronet_dataframe.index.month
data=aeronet_dataframe.groupby(aeronet_dataframe['hour']).mean()


fig, ax = plt.subplots(1, figsize=(12,6))
plt.plot(data['Aeronet_station_1'],'g',label='AERONET Lamezia',linewidth=3)
plt.plot(data['Aeronet_station_2'],'g--',label='AERONET Blida',linewidth=3)
plt.plot(data['Aeronet_Minsk'],'g-.',label='AERONET Minsk',linewidth=3)
plt.plot(data['LE_station_1'],'b',label='LE Lamezia',linewidth=3)
plt.plot(data['LE_station_2'],'b--',label='LE Blida',linewidth=3)
plt.plot(data['LE_Minsk'],'b-.',label='LE Minsk',linewidth=3)
plt.ylabel('AOD 550nm')
plt.xlabel('Hour')
plt.legend()
plt.savefig('./Figures/Daily_Cycle_Aeronet.png',format='png', dpi=1000,bbox_inches = "tight")

plt.show()

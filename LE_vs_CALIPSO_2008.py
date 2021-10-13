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

# "Viridis-like" colormap with white background
white_viridis = LinearSegmentedColormap.from_list('white_viridis', [
    (0, '#ffffff'),
    (1e-20, '#440053'),
    (0.2, '#404388'),
    (0.4, '#2a788e'),
    (0.6, '#21a784'),
    (0.8, '#78d151'),
    (1, '#fde624'),
], N=256)

#===Load user libraries====
sys.path.append('./libraries/')
from datamanager import DataManager_GRIDDED
from matplotlib import cm


viridis = cm.get_cmap('viridis', 256)

#===Create Color bar====
def make_Ramp( ramp_colors ): 
    from colour import Color
    from matplotlib.colors import LinearSegmentedColormap

    color_ramp = LinearSegmentedColormap.from_list( 'my_list', [ Color( c1 ).rgb for c1 in ramp_colors ] )
    return color_ramp


custom_ramp = make_Ramp( ["#ffffff","#aec8e0","#9dafb2","#e2b85f","#f0ae46","#ea9842","#d36f2b","#a04617","#6c423d"] ) 

gridded_calipso=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_calipso_V6_2008_gridded.nc')

calipso=np.ma.masked_invalid(gridded_calipso.nc.variables['yr'][:,:,:,:-1])
le=np.ma.masked_invalid(gridded_calipso.nc.variables['ys'][:])


calipso_mean=calipso.mean()
le_mean=le.mean()



# gridded_calipso.transversal(variable='ys',biascorr=1, extern=True,  calipso_mean=calipso_mean, le_mean=le_mean, vmin=0,vmax=3e-4,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='Profile extinction LE 1',direction='row',position=10,extend='max',cut=True,normalize=True)
# gridded_calipso.transversal(variable='ys',biascorr=1, extern=True,  calipso_mean=calipso_mean, le_mean=le_mean, vmin=0,vmax=3e-4,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='Profile extinction LE 2',direction='row',position=40,extend='max',cut=True,normalize=True)
# gridded_calipso.transversal(variable='ys',biascorr=1, extern=True,  calipso_mean=calipso_mean, le_mean=le_mean, vmin=0,vmax=3e-4,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='Profile extinction LE 3',direction='row',position=90,extend='max',cut=True,normalize=True)
# gridded_calipso.transversal(variable='yr',biascorr=1, extern=True,  calipso_mean=calipso_mean, le_mean=le_mean, vmin=0,vmax=3e-4,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='Profile extinction CALIPSO 1',direction='row',position=10,extend='max',cut=False)
# gridded_calipso.transversal(variable='yr',biascorr=1, extern=True,  calipso_mean=calipso_mean, le_mean=le_mean, vmin=0,vmax=3e-4,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='Profile extinction CALIPSO 2',direction='row',position=40,extend='max',cut=False)
# gridded_calipso.transversal(variable='yr',biascorr=1, extern=True,  calipso_mean=calipso_mean, le_mean=le_mean, vmin=0,vmax=3e-4,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='Profile extinction CALIPSO 3',direction='row',position=90,extend='max',cut=False)

# gridded_calipso.transversal(variable='ys',biascorr=1,extern=True,  calipso_mean=calipso_mean, le_mean=le_mean, vmin=0,vmax=3e-4,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='Profile extinction LE 4',direction='col',position=30,extend='max',cut=True,normalize=True)
# gridded_calipso.transversal(variable='ys',biascorr=1,extern=True,  calipso_mean=calipso_mean, le_mean=le_mean, vmin=0,vmax=3e-4,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='Profile extinction LE 5',direction='col',position=80,extend='max',cut=True,normalize=True)
# gridded_calipso.transversal(variable='ys',biascorr=1,extern=True,  calipso_mean=calipso_mean, le_mean=le_mean, vmin=0,vmax=3e-4,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='Profile extinction LE 6',direction='col',position=50,extend='max',cut=True,normalize=True)
# gridded_calipso.transversal(variable='yr',biascorr=1,extern=True,  calipso_mean=calipso_mean, le_mean=le_mean, vmin=0,vmax=3e-4,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='Profile extinction CALIPSO 4',direction='col',position=30,extend='max',cut=False)
# gridded_calipso.transversal(variable='yr',biascorr=1,extern=True,  calipso_mean=calipso_mean, le_mean=le_mean, vmin=0,vmax=3e-4,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='Profile extinction CALIPSO 5',direction='col',position=80,extend='max',cut=False)
# gridded_calipso.transversal(variable='yr',biascorr=1,extern=True,  calipso_mean=calipso_mean, le_mean=le_mean, vmin=0,vmax=3e-4,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='Profile extinction CALIPSO 6',direction='col',position=50,extend='max',cut=False)

# gridded_calipso.graph_map(variable='ys',level=0,biascorr=1,  vmin=0,vmax=3e-4,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=False,title='Level 0 extinction LE',extend='max')
# gridded_calipso.graph_map(variable='ys',level=1,biascorr=1, vmin=0,vmax=3e-4,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=False,title='Level 1 extinction LE',extend='max')
# # # gridded_calipso.graph_map(variable='ys',level=2,biascorr=1,extern=True,  calipso_mean=calipso_mean, le_mean=le_mean, vmin=0,vmax=3e-4,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=False,title='Level 2 extinction LE',extend='max')
# # # gridded_calipso.graph_map(variable='ys',level=3,biascorr=1,extern=True,  calipso_mean=calipso_mean, le_mean=le_mean, vmin=0,vmax=3e-4,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=False,title='Level 3 extinction LE',extend='max')
# # # gridded_calipso.graph_map(variable='ys',level=4,biascorr=1,extern=True,  calipso_mean=calipso_mean, le_mean=le_mean, vmin=0,vmax=3e-4,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=False,title='Level 4 extinction LE',extend='max')

# gridded_calipso.graph_map(variable='yr',level=0,biascorr=1,  vmin=0,vmax=3e-4,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='V5 Level 0 extinction CALIPSO',extend='max')
# gridded_calipso.graph_map(variable='yr',level=1,biascorr=1, vmin=0,vmax=3e-4,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='V5 Level 1 extinction CALIPSO',extend='max')
# gridded_calipso.graph_map(variable='yr',level=2,biascorr=1,vmin=0,vmax=3e-4,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='V5 Level 2 extinction CALIPSO',extend='max')
# gridded_calipso.graph_map(variable='yr',level=3,biascorr=1, vmin=0,vmax=3e-4,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='V5 Level 3 extinction CALIPSO',extend='max')
# gridded_calipso.graph_map(variable='yr',level=4,biascorr=1,vmin=0,vmax=3e-4,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=False,save=True,title='V5 Level 4 extinction CALIPSO',extend='max')


path_calipso='/Users/santiago/Documents/CALIPSO/2008/'
onlyfiles = [f for f in listdir(path_calipso) if isfile(join(path_calipso, f))]
onlyfiles.sort()

n=len(onlyfiles[1:200])

calipso_raw=np.zeros([344,n])
for orbit in range(1,n):
    print(orbit)
    calipso=Dataset(path_calipso+onlyfiles[orbit])
    calipso.set_auto_mask(False)


    CAD=np.flip(calipso.variables['CAD_Score'][:,:,0],1)
    CAD=CAD[:,:344]

    qs=np.flip(calipso.variables['Extinction_Coefficient_Uncertainty_532'][:],1)
    qs=qs[:,:344]
    qs[qs<0]=10000
    qc=np.flip(calipso.variables['Extinction_QC_Flag_532'][:],1)
    qc=qc[:,:344,0]

    column2=np.flip(calipso.variables['Extinction_Coefficient_532'][:],1) 
    column2=column2[:,:344]#To take just the troposphere
    column2=np.ma.masked_array(column2, mask=((column2<=0)|(column2>1)| (CAD>-80) | (qs>column2*0.50)| ((qc!=1) & (qc!=0))))
    #column2=np.ma.masked_array(column2, mask=((column2<=0)|(column2>1)| (CAD>-80) | (qs>column2*0.99)))
    #column2=np.ma.masked_array(column2, mask=((column2<=0)|(column2>1.25)))
    calipso_raw[:,orbit]=np.nanmean(column2,axis=0)


calipso_profile_raw=np.zeros([calipso_raw.shape[0]])
for i in range(calipso_raw.shape[0]):
    if (np.count_nonzero(~np.isnan(calipso_raw[i,:]))>100):
        calipso_profile_raw[i]=np.nanmean(calipso_raw[i,:])


orography=Dataset('/Users/santiago/Documents/LE_outputs/LE_Orography_meteo_20080101.nc')
orography=orography.variables['oro'][:]

#calipso_2=calipso.data
#calipso_2[calipso_2>100]=0
calipso_yearly=np.nanmean(calipso_raw,axis=0)

calipso_surface=np.zeros(calipso_yearly.shape)
calipso_surface[calipso_surface==0]=np.nan
nx=calipso_yearly.shape[0]
ny=calipso_yearly.shape[1]
nz=calipso_yearly.shape[2]
layer_meters=1000

for i in range(nx):
    for j in range(ny):
        for l in range(nz):
            altitude=l*layer_meters+layer_meters/2 - orography[i,j]; #Mean altitude layer
            if altitude>-300:
                indice=int(altitude//layer_meters)
                indice=max(0,indice)
                calipso_surface[i,j,indice]=calipso_yearly[i,j,l]


#calipso_surface[np.isnan(calipso_surface)]=0
le=le*calipso_mean/le_mean
le=le*le_mean/calipso_mean
calipso_profile_surface=np.nanmean(np.nanmean(calipso_surface,axis=0),axis=0)
calipso_profile=np.nanmean(np.nanmean(np.nanmean(calipso,axis=0),axis=0),axis=0)
le_profile=np.nanmean(np.nanmean(np.nanmean(le,axis=0),axis=0),axis=0)

y=np.arange(0,20)
y1=np.arange(0,calipso_profile_raw.shape[0])*60/1000-0.5
fig, ax = plt.subplots()
#ax.plot(le_profile*1000,y,'*-',label='LOTOS-EUROS',linewidth=3)
ax.plot(calipso_profile_raw,y1,'s-',label='CALIOP',linewidth=3)
ax.ticklabel_format(axis='x', style='sci', scilimits=(1,4))
ax.legend()
ax.set_xlabel(r'Extinction Coefficient [km$^{-1}$]',fontsize=12)
ax.set_ylabel(r'Altitud [km]',fontsize=12)
ax.set_ylim([0,20])
#ax.set_xlim([0,1.5e-1])
#plt.savefig('./Figures/Average_AEC_profile.png',format='png', dpi=1000,bbox_inches = "tight")
plt.show()
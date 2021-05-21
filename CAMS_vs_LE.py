import sys
#===Load user libraries====
sys.path.append('./libraries/')
from datamanager import DataManager_LE
from datamanager import DataManager_POLDER
from datamanager import DataManager_GRIDDED
from datamanager import DataManager_CAMS
from matplotlib import cm
viridis = cm.get_cmap('viridis', 256)

#===Create Color bar====
def make_Ramp( ramp_colors ): 
    from colour import Color
    from matplotlib.colors import LinearSegmentedColormap

    color_ramp = LinearSegmentedColormap.from_list( 'my_list', [ Color( c1 ).rgb for c1 in ramp_colors ] )
    return color_ramp

custom_ramp = make_Ramp( ["#ffffff","#aec8e0","#9dafb2","#e2b85f","#f0ae46","#ea9842","#d36f2b","#a04617","#6c423d"] ) 

cams=DataManager_CAMS('/Users/santiago/Documents/CAMS_reanalysis/CAMS_AOD_2008.nc')
LE=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_aod_2008.nc')

lon_lim=[LE.nc.variables['lon'][0],LE.nc.variables['lon'][-1]]
lat_lim=[LE.nc.variables['lat'][0],LE.nc.variables['lat'][-1]]
vmin=0
vmax=0.6
n_levels_plot=30
ini_dia=0
end_dia=31
cams.graph_map(variable='aod550',title='CAMS AOD 550nm',cmap=custom_ramp,n_levels_plot=n_levels_plot,vmin=vmin, vmax=vmax,type_graph='mean',save=True, ini_period=ini_dia,end_period=end_dia*8-1,lon_lim=lon_lim,lat_lim=lat_lim)
cams.graph_map(variable='aod670',title='CAMS AOD 670nm',cmap=custom_ramp,n_levels_plot=n_levels_plot,vmin=vmin, vmax=vmax,type_graph='mean',save=True, ini_period=ini_dia,end_period=end_dia*8-1,lon_lim=lon_lim,lat_lim=lat_lim)
cams.graph_map(variable='aod865',title='CAMS AOD 8650nm',cmap=custom_ramp,n_levels_plot=n_levels_plot,vmin=vmin, vmax=vmax,type_graph='mean',save=True, ini_period=ini_dia,end_period=end_dia*8-1,lon_lim=lon_lim,lat_lim=lat_lim)

LE.graph_map(variable='aod_563nm',title='LE AOD 550nm',cmap=custom_ramp,n_levels_plot=n_levels_plot,vmin=vmin, vmax=vmax,type_graph='mean',save=True, ini_period=ini_dia,end_period=end_dia*24-3,step_period=3)
LE.graph_map(variable='aod_670nm',title='LE AOD 670nm',cmap=custom_ramp,n_levels_plot=n_levels_plot,vmin=vmin, vmax=vmax,type_graph='mean',save=True, ini_period=ini_dia,end_period=end_dia*24-1,step_period=3)
LE.graph_map(variable='aod_865nm',title='LE AOD 865nm',cmap=custom_ramp,n_levels_plot=n_levels_plot,vmin=vmin, vmax=vmax,type_graph='mean',save=True, ini_period=ini_dia,end_period=end_dia*24-1,step_period=3)


for i in [2,3,4,5,6,7,8,9,10,11]:
    if i==3: #Marzo, sumo solo 29 dias
        ini_dia=ini_dia+29
        end_dia=end_dia+31
    elif i==2:
        ini_dia=ini_dia+31
        end_dia=end_dia+29
    elif i in [4,6,9,11]:
        ini_dia=ini_dia+31
        end_dia=end_dia+30
    elif i ==8:   
        ini_dia=ini_dia+31
        end_dia=end_dia+31
    else:
        ini_dia=ini_dia+30
        end_dia=end_dia+31
        
        
    cams.graph_map(variable='aod550',title='CAMS AOD 550nm',cmap=custom_ramp,n_levels_plot=n_levels_plot,vmin=vmin, vmax=vmax,type_graph='mean',save=True, ini_period=ini_dia*(8),end_period=end_dia*8-1,lon_lim=lon_lim,lat_lim=lat_lim)
    cams.graph_map(variable='aod670',title='CAMS AOD 670nm',cmap=custom_ramp,n_levels_plot=n_levels_plot,vmin=vmin, vmax=vmax,type_graph='mean',save=True, ini_period=ini_dia*(8),end_period=end_dia*8-1,lon_lim=lon_lim,lat_lim=lat_lim)
    cams.graph_map(variable='aod865',title='CAMS AOD 865nm',cmap=custom_ramp,n_levels_plot=n_levels_plot,vmin=vmin, vmax=vmax,type_graph='mean',save=True, ini_period=ini_dia*(8),end_period=end_dia*8-1,lon_lim=lon_lim,lat_lim=lat_lim)
    LE.graph_map(variable='aod_563nm',title='LE AOD 550nm',cmap=custom_ramp,n_levels_plot=n_levels_plot,vmin=vmin, vmax=vmax,type_graph='mean',save=True, ini_period=ini_dia*24,end_period=end_dia*24-3,step_period=3)
    LE.graph_map(variable='aod_670nm',title='LE AOD 670nm',cmap=custom_ramp,n_levels_plot=n_levels_plot,vmin=vmin, vmax=vmax,type_graph='mean',save=True, ini_period=ini_dia*24,end_period=end_dia*24-1,step_period=3)
    LE.graph_map(variable='aod_865nm',title='LE AOD 865nm',cmap=custom_ramp,n_levels_plot=n_levels_plot,vmin=vmin, vmax=vmax,type_graph='mean',save=True, ini_period=ini_dia*24,end_period=end_dia*24-1,step_period=3)


cams.graph_map(variable='aod550',title='CAMS AOD 550nm',cmap=custom_ramp,n_levels_plot=n_levels_plot,vmin=vmin, vmax=vmax,type_graph='mean',save=True, ini_period=end_dia*8,end_period=-1,lon_lim=lon_lim,lat_lim=lat_lim)
cams.graph_map(variable='aod670',title='CAMS AOD 670nm',cmap=custom_ramp,n_levels_plot=n_levels_plot,vmin=vmin, vmax=vmax,type_graph='mean',save=True, ini_period=end_dia*8,end_period=-1,lon_lim=lon_lim,lat_lim=lat_lim)
cams.graph_map(variable='aod865',title='CAMS AOD 865nm',cmap=custom_ramp,n_levels_plot=n_levels_plot,vmin=vmin, vmax=vmax,type_graph='mean',save=True, ini_period=end_dia*8,end_period=-1,lon_lim=lon_lim,lat_lim=lat_lim)
LE.graph_map(variable='aod_563nm',title='LE AOD 550nm',cmap=custom_ramp,n_levels_plot=n_levels_plot,vmin=vmin, vmax=vmax,type_graph='mean',save=True, ini_period=end_dia*24,end_period=-1,step_period=3)
LE.graph_map(variable='aod_670nm',title='LE AOD 670nm',cmap=custom_ramp,n_levels_plot=n_levels_plot,vmin=vmin, vmax=vmax,type_graph='mean',save=True, ini_period=end_dia*24,end_period=-1,step_period=3)
LE.graph_map(variable='aod_865nm',title='LE AOD 865nm',cmap=custom_ramp,n_levels_plot=n_levels_plot,vmin=vmin, vmax=vmax,type_graph='mean',save=True, ini_period=end_dia*24,end_period=-1,step_period=3)


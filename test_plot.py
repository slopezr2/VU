import sys
#===Load user libraries====
sys.path.append('./libraries/')
from datamanager import DataManager_LE
from datamanager import DataManager_POLDER
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

gridded_LE=DataManager_GRIDDED('/Users/santiago/Documents/POLDER/LE_POLDER/gridded/LE_20080501_1500_gridded.nc')
gridded_LE.graph_map(variable='ys',biascorr=2.5,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='instant',time=0,date=False,save=False,title='LE AOD565 Biascorr=2.5 2008-05-01')
gridded_LE.graph_map(variable='yr',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='instant',time=0,date=False,save=False,title='POLDER AOD565 2008-05-01')


# aod_file=DataManager_LE('/Users/santiago/Documents/POLDER/LE_POLDER/LE_aod.nc')


# mypath='/Users/santiago/Documents/POLDER/2008/05/'
# polder_file=DataManager_POLDER(mypath)


# for le_wave in ['490','563','670','865']:
#     aod_file.graph_map(biascorr=1,variable='aod_'+le_wave +'nm',vmin=0,vmax=1.2,n_levels_plot=15,cmap=custom_ramp,type_graph='mean',date=False,save=True,title_on=True,title='LE AOD'+le_wave +' average 2008-05')
#     aod_file.graph_map(biascorr=2.5,variable='aod_'+le_wave +'nm',vmin=0,vmax=1.2,n_levels_plot=15,cmap=custom_ramp,type_graph='mean',date=False,save=True,title_on=True,title='LE AOD'+le_wave +' average 2008-05 Biascorr=2.5')
# for polder_wave in  ['490','565','670','865']:
#     polder_file.graph_map(variable='AOD'+polder_wave,vmin=0,vmax=1.2,n_levels_plot=15,cmap=custom_ramp,type_graph='mean',save=True,title_on=True,title='POLDER AOD'+polder_wave +' average 2008-05')







# aod_file.graph_map(biascorr=1,variable='angstrom_modis',vmin=0,vmax=1.8,n_levels_plot=50,cmap=viridis,type_graph='mean',date=False,save=True,title_on=True,title='LE Angstrom Modis average 2008-05')
# aod_file.graph_map(biascorr=1,variable='angstrom_aeronet',vmin=0,vmax=1.8,n_levels_plot=50,cmap=viridis,type_graph='mean',date=False,save=True,title_on=True,title='LE Ansgtrom Aeronet average 2008-05')
# polder_file.graph_map(variable='AExp',vmin=0,vmax=1.8,n_levels_plot=50,cmap=viridis,type_graph='mean',save=True,title_on=True,title='POLDER Angstrom average 2008-05')







# for level in range(5):
#     aod_file.graph_map(variable='tau_550nm',vmin=0,vmax=0.1,level=level,n_levels_plot=15,cmap=custom_ramp,type_graph='mean',save=True)
#     aod_file.graph_map(variable='extinction_550nm',vmin=0,vmax=1e-4,level=level,n_levels_plot=15,cmap=custom_ramp,type_graph='mean',save=True)
    


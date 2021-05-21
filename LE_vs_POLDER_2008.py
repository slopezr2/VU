import sys
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



for month in ['01','02','03','04','05','06','07','08','09','10','11','12']:
    file='/Users/santiago/Documents/LE_outputs/2008_POLDER_Complete/gridded/LE_2008'+month+'.nc'
    
    
    gridded_LE=DataManager_GRIDDED(file)
    gridded_LE.graph_map(variable='ys',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=True,save=True,title='LE AOD565')
    gridded_LE.graph_map(variable='yr',biascorr=1,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=True,save=True,title='POLDER AOD565')
    gridded_LE.graph_map(variable='ys',biascorr=2.5,vmin=0,vmax=1,n_levels_plot=20,cmap=custom_ramp,type_graph='mean',date=True,save=True,title='LE AOD565 Biascorr=2.5')

gridded_LE_year=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_POLDER_Complete/gridded/LE_2008.nc')
gridded_LE_year.timeseries(variable=['ys','yr'],biascorr=[2.5,1],label=['LE Biascorr=2.5','POLDER'],type_graph='sum',save=True,title='Domain Sum time series')


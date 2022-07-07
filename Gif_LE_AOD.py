import numpy as np
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.morphology import generate_binary_structure, binary_erosion
import matplotlib.pyplot as plt
import sys
import imageio
from os import listdir
from os.path import isfile, join
#===Load user libraries====
sys.path.append('./libraries/')
from datamanager import DataManager_GRIDDED
from datamanager import DataManager_LE
from matplotlib.colors import ListedColormap
import met_brewer
g=met_brewer.met_brew('Demuth',n=256,brew_type='continuous')
cmap = ListedColormap(g, name='myColorMap', N=len(g))




# gridded_LE_year=DataManager_GRIDDED('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_m_polder-aod-563_2008.nc')

# mdata = np.ma.filled(gridded_LE_year.nc.variables['yr'][:], np.nan)
# quartil_3rd = np.nanpercentile(mdata,90)

# #applying the detection and plotting results
# dias = []
# for i in range(mdata.shape[0]):
#     dia_polder = mdata[i,:,:,0].copy()
#     dia_polder[dia_polder<quartil_3rd] = np.nan
#     if np.sum(~np.isnan(dia_polder))>300:
#         dias.append(i)


LE_year=DataManager_LE('/Users/santiago/Documents/LE_outputs/2008_Complete/LE_aod_2008.nc')
onlyfiles = []
for i in range(0,600,3):
    i_str = f'{i:03d}'
    LE_year.graph_map(variable='aod_563nm',cmap=cmap.reversed(),vmin=0,vmax=0.3,type_graph='instant',time=i,save=True,stock_image=False,save_title='/gif/LE_AOD_figure_'+i_str,extend='max')
    onlyfiles.append('/gif/LE_AOD_figure_'+i_str+'.png')

mypath='./Figures'
onlyfiles.sort()
images=[]
for file in onlyfiles:
    if file[5:11]=='LE_AOD':
        images.append(imageio.imread(mypath+file))

imageio.mimsave('./Figures/LE_AOD_2008.gif', images,duration=0.3)

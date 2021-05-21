from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from cartopy import config
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm
from os import listdir
from os.path import isfile, join



#===Create Color bar====
def make_Ramp( ramp_colors ): 
    from colour import Color
    from matplotlib.colors import LinearSegmentedColormap

    color_ramp = LinearSegmentedColormap.from_list( 'my_list', [ Color( c1 ).rgb for c1 in ramp_colors ] )
    return color_ramp

custom_ramp = make_Ramp( ["#ffffff","#aec8e0","#9dafb2","#e2b85f","#f0ae46","#ea9842","#d36f2b","#a04617","#6c423d"] ) 

mypath='/Users/santiago/Documents/POLDER/2008/05/'
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f)) and f[0:12]=='GRASP_POLDER']

onlyfiles.sort()


#polder=Dataset('/Users/santiago/Documents/POLDER/GRASP_POLDER_L3_20120813.01degree.nc')
polder=Dataset(mypath+onlyfiles[0])
aux=polder.variables['AOD565'][:]
aod_polder=np.empty([aux.shape[0],aux.shape[1],len(onlyfiles)])
i=0
aod_polder[:,:,i]=aux.filled(fill_value=np.nan)

for file_name in onlyfiles[1:]:
    if file_name[0:12]=='GRASP_POLDER':
        i=i+1
        polder=Dataset(mypath+file_name)
        aux=polder.variables['AOD565'][:]
        aod_polder[:,:,i]=aux.filled(fill_value=np.nan)
    
        latitude=polder.variables['Latitude'][:]
        longitude=polder.variables['Longitude'][:]
        
    
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
european_borders=cfeature.NaturalEarthFeature(
                    category='cultural',
                    name='admin_0_countries',
                    scale='50m',
                    facecolor='none')


coastlines=cfeature.NaturalEarthFeature(
                    category='physical',
                    name='coastline',
                    scale='50m',
                    facecolor='none')
ax.add_feature(european_borders,edgecolor='black',linewidth=0.2)
ax.add_feature(coastlines,edgecolor='black',linewidth=1)
n_levels_plot=15
vmin=0
vmax=1.2
levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)


# pick the desired colormap, sensible levels, and define a normalization
# instance which takes data values and translates those into levels.
cmap = custom_ramp
cmap.set_bad([0.7, 0.7, 0.7],1.)
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
# im=ax.contour(longitude[0,1652:2148],latitude[0:350,0],aod_polder[0:350,1652:2148,i],levels=levels,cmap=custom_ramp,extend='both')
plt.title('POLDER average 2008-05')      
im=plt.pcolor(longitude[0,1652:2148],latitude[0:350,0],np.nanmean(aod_polder[0:350,1652:2148,:],axis=2),cmap=cmap,norm=norm,vmin=vmin,vmax=vmax)
cax = fig.add_axes([ax.get_position().x1+0.02,
                    ax.get_position().y0+0.09,
                    0.08,
                    ax.get_position().height-0.2])
 
plt.colorbar(im,cax=cax,extend='both')
#plt.savefig('./Figures/Polder_20120220_'+str(i)+'.png',format='png', dpi=1000,bbox_inches = "tight")
plt.show()


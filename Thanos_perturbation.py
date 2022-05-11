import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.stats as sts
from netCDF4 import Dataset
from matplotlib import cm
sys.path.append('./libraries/')
from datamanager import graph_map

RdBu = cm.get_cmap('RdBu', 256)

def getmeanofneighbors(matrix, i, j,r):
    region = matrix[max(0, i-r) : i+r+1,
                    max(0, j-r) : j+r+1]
    #print(region)
    #return region
    return (np.sum(region) - matrix[i, j])/(region.size-1) # Sum the region and subtract center



nx=100
ny=110
r=2
repeat=5
new_mean=1
new_std=0.65


domain=1+np.random.randn(nx,ny)*10
domain_backup2=domain.copy()
plt.imshow(domain,vmin=0,vmax=2,cmap='plasma')
plt.title('mean='+str(domain.mean())[0:4]+' std='+str(domain.std())[0:4])
plt.colorbar()
plt.show()

for k in range(repeat):
    domain_backup=domain.copy()
    for i in range(domain.shape[0]):
        for j in range(domain.shape[1]):
            domain[i,j]=getmeanofneighbors(domain_backup, i, j,r)
            
    plt.imshow(domain,vmin=0,vmax=2,cmap='plasma')
    plt.title('mean='+str(domain.mean())[0:4]+' std='+str(domain.std())[0:4])
    plt.colorbar()
    plt.show()

domain=np.exp(domain)
plt.imshow(domain,vmin=0,vmax=2,cmap='plasma')
plt.title('mean='+str(domain.mean())[0:4]+' std='+str(domain.std())[0:4])
plt.colorbar()
plt.show()

domain=new_std*((domain-domain.mean())/domain.std())+new_mean
plt.imshow(domain,vmin=0,vmax=2,cmap='plasma')
plt.title('mean='+str(domain.mean())[0:4]+' std='+str(domain.std())[0:4])
plt.colorbar()
plt.show()

samples=domain.ravel()
h,e = np.histogram(samples, bins=100, density=True)
x = np.linspace(e.min(), e.max())

# plot the histogram
plt.figure(figsize=(8,6))
plt.bar(e[:-1], h, width=np.diff(e), ec='k', align='edge', label='histogram')

# plot the real KDE
kde = sts.gaussian_kde(samples)
plt.plot(x, kde.pdf(x), c='r', lw=4, label='KDE')
plt.xlim([0,5])

samples=domain_backup2.ravel()
h,e = np.histogram(samples, bins=100, density=True)
x = np.linspace(e.min(), e.max())

# plot the histogram
plt.figure(figsize=(8,6))
plt.bar(e[:-1], h, width=np.diff(e), ec='k', align='edge', label='histogram')

# plot the real KDE
kde = sts.gaussian_kde(samples)
plt.plot(x, kde.pdf(x), c='r', lw=4, label='KDE')



# dc=Dataset('/Users/santiago/Downloads/LE_Noise_Thanos_V1_dc_20080502_xi01a.nc')
# dc_map=dc.variables['dc'][0,0,0,:,:]
# samples=dc_map.ravel()
# h,e = np.histogram(samples, bins=100, density=True)
# x = np.linspace(e.min(), e.max())

# fig,ax=plt.subplots(figsize=(8,6))
# # plot the histogram

# ax.bar(e[:-1], h, width=np.diff(e), ec='k', align='edge', label='histogram')

# # plot the real KDE
# kde = sts.gaussian_kde(samples)
# ax.plot(x, kde.pdf(x), c='r', lw=4, label='KDE')
# ax.set_xlim([0,5])

# fig.savefig('./Figures/Hist_dc_Thanos.png',dpi=300,bbox_inches='tight')


# graph_map(variable=dc_map,lat=dc.variables['latitude'][:],lon=dc.variables['longitude'][:],cmap=RdBu.reversed(),date=False,title='Perturbation factor',vmax=2, dpi=300,save=True)
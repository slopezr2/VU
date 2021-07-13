


iasi_so2_CSO=Dataset('/Users/santiago/Documents/IASI/SO2/CSO/CSO_IASI_20080510.nc')

lon=iasi_so2_CSO.variables['longitude'][:]
lat=iasi_so2_CSO.variables['latitude'][:]
iasi_so2=iasi_so2_CSO.variables['SO2'][:]
vmin=0
vmax=3
n_levels_plot=10
#cmap=viridis
cmap=RdYlGn.reversed()
biascorr=64.0638*1000

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
levels = MaxNLocator(nbins=n_levels_plot).tick_values(vmin, vmax)
ax.set_facecolor((1, 1, 1))

norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
plt.scatter(lon, lat, 5,c=biascorr*iasi_so2, 
                   transform=ccrs.PlateCarree(),cmap=cmap, norm=norm,vmin=vmin, vmax=vmax)
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
plt.xlim(-15,35)
plt.ylim(35,70)

        

cax = fig.add_axes([ax.get_position().x1+0.02,
         ax.get_position().y0+0.09,
         0.08,
         ax.get_position().height-0.2])

plt.colorbar(cax=cax,extend='both')
plt.show()


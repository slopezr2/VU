#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 11:43:18 2021

@author: fermat
"""

# imports 
import matplotlib
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from os import listdir
from os.path import isfile, join
from wrf import (getvar, to_np, get_cartopy, latlon_coords, vertcross,
                 cartopy_xlim, cartopy_ylim, interpline, CoordPair)
# path="/media/fermat/DD_JOE1/calificatorio_1/data/febrero"
# path="/media/fermat/DD_JOE1/calificatorio_1/data"
path="/media/fermat/DD_JOE1/calificatorio_1/jimenez"

p_level = 950
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
onlyfiles.sort()
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  5 11:47:09 2021

@author: User
"""
#sns.set_theme(style="whitegrid")
# sns.set_theme(style="white")
p_level = 950
level=7
level_win=20

zorder_pressure=10
zorder_mountain=1000
zorder_velocity=0
zorder_contour=0
w_levels = np.arange(0, 20, 2)
#### Info del corte
latmin=5; latmax=5
lonmin=-74; lonmax=-66
    
cross_start = CoordPair(lat=latmin, lon=lonmin); cross_end   = CoordPair(lat=latmax, lon=lonmax)
 
i=0
#lst=2  #0300, 0700, 1100, 1500, 1900, y 2300 ### Hora local de interes
lst=0

print('vamos en la hora ' + str(lst))
hour=lst#+5 # hora UTC
step=1 # cada cuanto leo archivos.
#lst=hour-5
name='crossSectionJimenez'#+ str(lst)
# name='plot_average'
for file in onlyfiles[hour:-1:step]: #onlyfiles[hour:-1:step]: #onlyfiles[0:-1:24] 

 if file[0:5]=='wrfou':
        i=i+1
        print(file)
        wrf_file = Dataset(path+"/"+file)
        print(path+"/"+file)
        p = getvar(wrf_file, "pressure")
        w = getvar(wrf_file, "wspd_wdir", units="m s-1")[0,:]
        print('here')
        if i==1:
           p_mean=p
           w_mean=w
        else:
            p_mean=p+p_mean
            w_mean=w+w_mean
print('vamos en la hora ' + str(lst))
p_mean=p_mean/i
w_mean=w_mean/i

 # wrf_file = Dataset("/home/fermat/Dropbox/EAFIT/Calificatorio_1/wrfout_d01_2014-01-30_00")
 # wrf_file = Dataset("D:/Dropbox/EAFIT/Calificatorio_1/wrfout_d01_2014-01-30_00")

# get the lat, lon grid
lats, lons = latlon_coords(p_mean)
ht = getvar(wrf_file, "z",units="km")
ht = ht[0:level_win,:,:]
ter = getvar(wrf_file, "ter",units="km")

w = w_mean[0:level_win,:,:]
W = 10**(w/10.) # Use linear Z for interpolation
w_cross = vertcross(W, ht, wrfin=wrf_file,
                    start_point=cross_start,
                    end_point=cross_end,
                    latlon=True, meta=True)
ter_line = interpline(ter, wrfin=wrf_file, start_point=cross_start,
                      end_point=cross_end)
xs = np.arange(0, w_cross.shape[-1], 1)
ys = to_np(w_cross.coords["vertical"])

w_cross_filled = np.ma.copy(to_np(w_cross))
for i in range(w_cross_filled.shape[-1]):
    column_vals = w_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -20).nonzero())[0])
    w_cross_filled[0:first_idx, i] = w_cross_filled[first_idx, i]
y=p[0:level,:,:]
modelLevels=None

labelFontSize = "small"
#     # create figure and axes
# Axes=plt.axes()
fig = plt.figure()

# ax3 = plt.Axes(fig)#, rect, kwargs) axes()
# ax1 = fig.add_subplot(111)
ax1 = plt.axes()

    # according to model levels on the right y axis
ax1.set_ylabel("Pressure [hPa]")
ax1.set_yscale('log')
ax1.set_ylim(10.*np.ceil(y.max()/10.), y.min()) # avoid truncation of 1000 hPa
subs = [1,2,15]
if y.max()/y.min() < 30.:
    subs = [1,2,3,4,5,6,7,8,9]
y1loc = matplotlib.ticker.LogLocator(base=10., subs=subs)
ax1.yaxis.set_major_locator(y1loc)
fmt = matplotlib.ticker.FormatStrFormatter("%g")
ax1.yaxis.set_major_formatter(fmt)
for t in ax1.get_yticklabels():
    t.set_fontsize(labelFontSize)
# calculate altitudes from pressure values (use fixed scale height)
z0 = 8.400    # scale height for pressure_to_altitude conversion [km]
altitude = z0 * np.log(1015.23/y)
# add second y axis for altitude scale 
ax2 = ax1.twinx()
# change values and font size of x labels
xloc = matplotlib.ticker.FixedLocator(np.arange(-74.,-64.,2.))
ax1.xaxis.set_major_locator(xloc)
# for t in ax1.get_xticklabels():
#     t.set_fontsize(labelFontSize)
# draw horizontal lines to the right to indicate model levels
if not modelLevels is None:
    pos = ax1.get_position()
    axm = fig.add_axes([pos.x1,pos.y0,0.02,pos.height], sharey=ax2)
    axm.set_xlim(0., 1.)
    axm.xaxis.set_visisble(False)
    modelLev = axm.hlines(altitude, 0., 1., color='0.5')
    axr = axm     # specify y axis for right tick marks and labels
    # turn off tick labels of ax2
    for t in ax2.get_yticklabels():
        t.set_visible(False)
    label_xcoor = 3.7
else:
    axr = ax2
    label_xcoor = 1.05
axr.set_ylabel("Altitude [km]")
axr.yaxis.set_label_coords
axr.yaxis.set_label_coords(label_xcoor+0.05, 0.5)
axr.set_ylim(altitude.min(), altitude.max())
yrloc = matplotlib.ticker.MaxNLocator(steps=[1,2,5,10])
axr.yaxis.set_major_locator(yrloc)
axr.yaxis.tick_right()
for t in axr.yaxis.get_majorticklines():
    t.set_visible(False)
for t in axr.get_yticklabels():
    t.set_fontsize(labelFontSize)
#########################################3
plot_windspd= ax2.pcolormesh(xs, ys, w_cross_filled,
                         cmap= plt.get_cmap('rainbow'),
                         zorder=zorder_velocity, alpha = 0.7, shading='auto')#,vmin=0,vmax=14)  

ht_fill = ax2.fill_between(xs, 0, to_np(ter_line),
                                facecolor="saddlebrown",zorder=zorder_mountain)
contours = plt.contour(xs, ys, to_np(w_cross_filled),
                       alpha=0.7,levels=w_levels, linestyles='dashed',
                                 colors='black',zorder=zorder_contour) 
plt.clabel(contours,colors = "black", inline=1, fontsize=12,fmt="%i")
x=np.arange(600,1000,50)
plt.ylim((0, 5))

cbar = fig.colorbar(plot_windspd, ax=ax1,
                    orientation="horizontal", pad=.05,extend="max")
cbar.set_label(r'ms$^{-1}$',fontsize=12)
# a=cbar.get_ticks
for p in x:
    ax1.axhline(p,color='black',zorder=zorder_pressure,alpha=1) 

# plt.xticks(np.arange(lonmin,lonmax))

# ax1.get_xticklabels()
a=ax1.get_xticks()
a
# ax1.set_xticks(a)
xaxis=np.linspace(-75, -67, num=5)
ax1.set_xticklabels([i for i in xaxis])

# ax1.set_xticklabels([i for i in [-74,-72,-71,-70,-69,-68,-67,-66]])
ax1.set_xticks# set_xticks
# show plot
plt.savefig(name, dpi=1000,bbox_inches = "tight")
plt.show()


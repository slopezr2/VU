#! /usr/bin/env python



import mdf
import go
import datetime
import numpy
import matplotlib.pyplot as plt





runids = ['oldLU','LUOR']

t = datetime.datetime(2016,4,1,20,0)

res = {}

for runid in runids :

    fname = 'LE_Dcol_%s_meteo_%s.nc' % (runid,t.strftime('%Y%m%d'))
    lons = mdf.get_var(fname,'longitude')
    lats = mdf.get_var(fname,'latitude')
    oro = mdf.get_var(fname,'oro')
    h = mdf.get_var(fname,'h')
    temper = mdf.get_var(fname,'temper')

    nt,nz,ny,nx = temper.shape

#Medellin
    iy=164


    ix0,ix1 = 20,150
    irec = t.hour

    zz = numpy.zeros((nz+3,nx),float)
    zz[0,:] = 0.0
    zz[1,:] = oro[iy,:]
    zz[2:nz+2,:] = oro[iy,:] + h[irec,:,iy,:]
    zz[nz+2,:] = 10e3
    
    xx = numpy.zeros((nz+3,nx))
    for k in range(nz+3) : xx[k,:] = lons
    
    # aver over lons:
    zz2 = (zz[:,0:nx-1] + zz[:,1:nx])/2.0
    xx2 = (xx[:,0:nx-1] + xx[:,1:nx])/2.0

    pat = numpy.ma.array( data=numpy.zeros((nz+2,nx),float), 
                          mask=numpy.ones((nz+2,nx),bool) )
    
    
    res[runid]=pat
    
    pat[1:nz+1,:] = temper[irec,:,iy,:]

    title = '%s - %s (UTC) at %10.1fN' % (runid,t.strftime('%Y-%m-%d %H:%M'),lats[iy])
    label = 'temperature [oC]'
    
    fig = go.plot2.QuickPat( pat[:,ix0:ix1-1]-273, 
                             xx=xx2[:,ix0:ix1], yy=zz2[:,ix0:ix1],
                              vmin=-10, vmax=30,
                              cbar=dict(label=label),
                              figsize=(9,5) )
    fig.ax.set_xlabel('longitude [degrees east]')
    fig.ax.set_ylabel('altitude [m]')
    fig.ax.set_ylim([0,8000])
    fig.ax.set_title( title )
    
    fig.Export( 'temper-xz-%s.png' % runid )
    
    #res[runid] = {}
    #res[runid]['xx'] = xx
    #res[runid]['zz'] = zz
    #res[runid]['temper'] = pat
    

#endfor

# difference figures in the dictionary
dpat = res[runids[1]] - res[runids[0]]
title = '%s - %s (UTC) at %10.1fN' % (runid,t.strftime('%Y-%m-%d %H:%M'),lats[iy])
label = 'Rest temperature between updated and old orography [oC]'
fig = go.plot2.QuickPat( dpat[:,ix0:ix1-1], 
                             xx=xx2[:,ix0:ix1], yy=zz2[:,ix0:ix1],
                              vmin=-5, vmax=5, cmap=dict(colors='pwb'),
                              cbar=dict(label=label),
                              figsize=(9,5) )
fig.ax.set_xlabel('longitude [degrees east]')
fig.ax.set_ylabel('altitude [m]')
fig.ax.set_ylim([0,8000])
fig.ax.set_title( title )
fig.Export( 'temper-xz-%s.png' % runid )
  
plt.show()


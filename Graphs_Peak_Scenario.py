#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 11:55:32 2022

@author: santiago
"""
subdomain='Benelux+Germany'
if subdomain=='Netherlands':
    lat_ini=60
    lat_end=76
    lon_ini=37
    lon_end=45
elif subdomain=='All':
    lat_ini=0
    lat_end=-1
    lon_ini=0
    lon_end=-1    
elif subdomain=='Benelux+Germany':
    lat_ini=50
    lat_end=85
    lon_ini=33
    lon_end=60
elif subdomain=='Scandinavia':
    lat_ini=85
    lat_end=140
    lon_ini=37
    lon_end=100
elif subdomain=='Iberia':
    lat_ini=2
    lat_end=45
    lon_ini=10
    lon_end=45
elif subdomain=='East Europe':
    lat_ini=20
    lat_end=80
    lon_ini=55
    lon_end=100

gridded_LE_year.graph_map(variable='yr',cmap=cmap.reversed(),vmin=0,vmax=0.7,type_graph='mean',time=i,ini_period=0,end_period=-1,save=True,
                          stock_image=True,save_title='Benelux_yearly_POLDER',date=True,grid=True,subdomain=True,lat_lim=[lat_ini,lat_end],lon_lim=[lon_ini,lon_end],title='POLDER-AOD565')
gridded_LE_year.graph_map(variable='yr',cmap=cmap.reversed(),vmin=0,vmax=0.7,type_graph='mean',time=i,ini_period=87,end_period=116,save=True,
                          stock_image=True,save_title='Benelux_April_POLDER',date=True,grid=True,subdomain=True,lat_lim=[lat_ini,lat_end],lon_lim=[lon_ini,lon_end],title='POLDER-AOD565')
gridded_LE_year.graph_map(variable='yr',cmap=cmap.reversed(),vmin=0,vmax=0.7,type_graph='mean',time=i,ini_period=87-30,end_period=116-30,save=True,
                          stock_image=True,save_title='Benelux_March_POLDER',date=True,grid=True,subdomain=True,lat_lim=[lat_ini,lat_end],lon_lim=[lon_ini,lon_end],title='POLDER-AOD565')
gridded_LE_year.graph_map(variable='yr',cmap=cmap.reversed(),vmin=0,vmax=0.7,type_graph='mean',time=i,ini_period=87+30,end_period=116+30,save=True,
                          stock_image=True,save_title='Benelux_May_POLDER',date=True,grid=True,subdomain=True,lat_lim=[lat_ini,lat_end],lon_lim=[lon_ini,lon_end],title='POLDER-AOD565')

gridded_Aeronet.graph_map(variable='yr',vmin=0,vmax=0.4,cmap=cmap.reversed(),type_graph='mean',ini_period=0,end_period=-1,date=True,save=True,title='AERONET-AOD550',stations=True,size_stations=60,grid=True,
                          stock_image=True,subdomain=True,lat_lim=[lat_ini,lat_end],lon_lim=[lon_ini,lon_end],save_title='Benelux_yearly_AERONET')
gridded_Aeronet.graph_map(variable='yr',vmin=0,vmax=0.4,cmap=cmap.reversed(),type_graph='mean',ini_period=1110,end_period=1540,date=True,save=True,title='AERONET-AOD550',stations=True,size_stations=60,grid=True,
                          stock_image=True,subdomain=True,lat_lim=[lat_ini,lat_end],lon_lim=[lon_ini,lon_end],save_title='Benelux_April_AERONET')
gridded_Aeronet.graph_map(variable='yr',vmin=0,vmax=0.4,cmap=cmap.reversed(),type_graph='mean',ini_period=690,end_period=1100,date=True,save=True,title='AERONET-AOD550',stations=True,size_stations=60,grid=True,
                          stock_image=True,subdomain=True,lat_lim=[lat_ini,lat_end],lon_lim=[lon_ini,lon_end],save_title='Benelux_March_AERONET')
gridded_Aeronet.graph_map(variable='yr',vmin=0,vmax=0.4,cmap=cmap.reversed(),type_graph='mean',ini_period=1550,end_period=2050,date=True,save=True,title='AERONET-AOD550',stations=True,size_stations=60,grid=True,
                          stock_image=True,subdomain=True,lat_lim=[lat_ini,lat_end],lon_lim=[lon_ini,lon_end],save_title='Benelux_May_AERONET')



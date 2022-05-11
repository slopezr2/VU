#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  2 13:03:52 2022

@author: ayarcebotero
"""

#===Archivo de ejemplo para Yarce===

from netCDF4 import Dataset  #La libreria para leer archivos .nc
#==Estas librerias siempre las cargo por si necesito usar algo, son basicas mejor dicho===


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm #Para los mapas de colores
import pandas as pd
import os


#===Esto es para importar las clases de la libreria que hice, que estan en otra carpeta====
import sys
sys.path.append('./libraries/')
from datamanager import DataManager_GRIDDED # Esta es para una salida de CSO, es decir, que tenga el modelo y el satelite ya en una grilla
from datamanager import DataManager_LE # Esta es para una salida cualquiera de LOTOS-EUROS


#====Algunos mapas de colores====
viridis = cm.get_cmap('viridis', 256)
RdYlGn=cm.get_cmap('Spectral', 256)
oranges=cm.get_cmap('Oranges', 256)
bwr=cm.get_cmap('bwr', 256)

#=====Salida normal de LE======
LE_file='/Users/santiago/Downloads/mean_nested_MACC.nc'
#LE=DataManager_LE(LE_file)
LE=DataManager_GRIDDED(LE_file)


#===Si queres ver las variables dentro del NC====
LE.nc.variables.keys()
print(LE.nc.variables.keys())

#==Si queres ver la informacion de una variable===
LE.nc.variables['ys']
print(LE.nc.variables['ys'])

#== Si queres acceder a los datos de una variable==
LE.nc.variables['ys'][:]
#print(LE.nc.variables['temper'][:])

# ===Para graficar===
LE.graph_map(variable='ys',biascorr=1,vmin=0,vmax=0.00024998753,n_levels_plot=10,level=0,cmap=cmo.haline,type_graph='instant',time=0,date=True,save=False,title='ys',orientation='horizontal')
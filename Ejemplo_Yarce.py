#===Archivo de ejemplo para Yarce===

from netCDF4 import Dataset  #La libreria para leer archivos .nc
#==Estas librerias siempre las cargo por si necesito usar algo, son basicas mejor dicho===
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm #Para los mapas de colores
import pandas as pd

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
LE_file='/Users/santiago/Documents/LE_outputs/LE_Orography_meteo_20080101.nc'
LE=DataManager_LE(LE_file)


#===Si queres ver las variables dentro del NC====
LE.nc.variables.keys()
print(LE.nc.variables.keys())

#==Si queres ver la informacion de una variable===
LE.nc.variables['temper']
print(LE.nc.variables['temper'])

#== Si queres acceder a los datos de una variable==
LE.nc.variables['temper'][:]
#print(LE.nc.variables['temper'][:])

# ===Para graficar===
LE.graph_map(variable='temper',biascorr=1,vmin=260,vmax=280,n_levels_plot=20,level=0,cmap=oranges,type_graph='mean',ini_period=0,end_period=-1,date=True,save=False,title='Temperatura',orientation='horizontal')
"""=Parametros importantes que no son muy obvios==
type_graph, si es 'mean' te grafica la media entre ini_period y end_period (-1 significa el ultimo). si es 'instant' te grafica unicamente el tiempo definido en la variable time
n_levels_plot, es el numero de niveles en el colorbar
cmap, es el colormap
date, si es True te ponde la fecha en el titulo, si es False, no la pone
biascorr, es un valor para multiplicar la variable

Igualmente, poder mirar cuales son los parametros disponibles, la mayoria es bastante claro que hacen solo con el nombre, o es de probar

"""

"""

Para leer la salida de CSO se usa es DataManager_GRIDDED, funciona exactamente igual, solo que en esos archivos siempre hay solo dos variables, 'ys' y 'yr'

"""



"""
Si queres cargar un NC sin usar esas librerias, podes usar el obejto Dataset de la lbreria netcdf4. Cuando vos usas alguna de las clases de mi libreria, se crea una variable
que se llama .nc, y es exactamente un Dataset

"""

LE_nc=Dataset(LE_file)
#===Si queres ver las variables dentro del NC====
LE_nc.variables.keys()
print(LE_nc.variables.keys())

#==Si queres ver la informacion de una variable===
LE_nc.variables['temper']
print(LE_nc.variables['temper'])

#== Si queres acceder a los datos de una variable==
LE_nc.variables['temper'][:]
#print(LE.nc.variables['temper'][:])

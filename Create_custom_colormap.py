#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 14:16:17 2022

@author: santiago
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib import cm


a=np.random.randn(100,100)
upper = cm.RdYlGn_r(np.arange(256))
lower = np.ones((1,4))
#lower = cm.viridis(np.arange(256))
value_union=1
for i in range(3):
  lower[-value_union:,i] = np.linspace(lower[-value_union,i], upper[0,i], value_union)


# combine parts of colormap
cmap = np.vstack(( lower, upper ))

# convert to matplotlib colormap
cmap = ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

RdYlGn=cm.get_cmap('Spectral', 256)

plt.imshow(a,cmap=cmap)
plt.colorbar()
from os import listdir
from os.path import isfile, join
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime,timedelta
from datetime import date as date_datetime
from netCDF4 import Dataset


path_calipso='/Users/santiago/Documents/CALIPSO/2008/'
onlyfiles = [f for f in listdir(path_calipso) if isfile(join(path_calipso, f))]
onlyfiles.sort()

y1=np.arange(0,344)*60/1000-0.5

obs_dates=pd.DataFrame()
for i in range(344):
    obs_dates[i]=np.full((366),np.nan)


obs_dates.index=pd.date_range(start='1/1/2008',end='31/12/2008',freq='D')



for file in onlyfiles[2:-1]:
    date=file[35:45]
    print(date)
    calipso=Dataset(path_calipso+file)
    calipso.set_auto_mask(False)
    CAD=np.flip(calipso.variables['CAD_Score'][:,:,0],1)
    CAD=CAD[:,:344]
    
    qs=np.flip(calipso.variables['Extinction_Coefficient_Uncertainty_532'][:],1)
    qs=qs[:,:344]
    qs[qs<0]=10000
    qc=np.flip(calipso.variables['Extinction_QC_Flag_532'][:],1)
    qc=qc[:,:344,0]
    
    column2=np.flip(calipso.variables['Extinction_Coefficient_532'][:],1) 
    column2=column2[:,:344]#To take just the troposphere
    column2=np.ma.masked_array(column2, mask=((column2<=0)|(column2>1)| (CAD>-80)| (CAD<-100) | (qs>column2*0.99)))
    column2[column2.mask==True]=np.nan
    for i in range(column2.shape[1]):
        n=np.count_nonzero(~np.isnan(column2[:,i]))
        if (n>0) and (np.isnan(obs_dates.loc[date][i])):
            obs_dates.loc[date][i]=y1[i]
            
columns=np.arange(344)
fig, ax = plt.subplots()
obs_dates[columns].plot(ax=ax,linestyle='None',marker='s',markersize=0.2,color='b',legend=False)
ax.set_ylabel('Altitude [km]',fontsize=12)
ax.set_ylim([-0.5,20])
plt.savefig('./Figures/Calipso_dates_observations.png',format='png', dpi=300,bbox_inches = "tight")

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
import pandas as pd
#Create listing table
template=pd.read_csv('/Users/santiago/scratch/CAMS-61/WP6130/CSO-data/S5p/listing-NO2-CAMS.csv',delimiter=';')
list_CSO={name: [] for name in list(template.columns)}
valid_files=[]

#AOD565_SD=[1.04,0.5,0.34,0.32,0.29,0.34]
AOD565_SD=888*np.ones(6)

AOD865_SD=[1.07,0.42,0.44,0.36,0.38,0.39]
AOD443_SD=[1.04,0.5,0.34,0.32,0.29,0.34]
AAOD565_SD=[2.14,1.36,0.96,0.74,0.64,0.63]
AAOD865_SD=[1.04,0.5,0.34,0.32,0.29,0.34]
AAOD443_SD=[1.04,0.5,0.34,0.32,0.29,0.34]
SSA565_SD=[1.04,0.5,0.34,0.32,0.29,0.34]
SSA865_SD=[1.04,0.5,0.34,0.32,0.29,0.34]
SSA443_SD=[1.04,0.5,0.34,0.32,0.29,0.34]
AExp_SD=[1.04,0.5,0.34,0.32,0.29,0.34]

#months=['01','02','03','04','05','06','07','08','09','10','11','12']
months=['01']
for i in months:
    mypath='/Users/santiago/Documents/POLDER/Model/'+i+'/'
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    onlyfiles.sort()
    cso_path=mypath+'/CSO/'
    
    
    for file_name in onlyfiles[1:2]:
        if file_name[0:12]=='GRASP_POLDER':
            # Reading POLDER data
            print(file_name)
            polder=Dataset(mypath+file_name)
            aod_polder=polder.variables['AOD565'][0:350,1652:2148]
            aod_polder_aux=aod_polder.copy()
            aod_polder_443=polder.variables['AOD443'][0:350,1652:2148]
            aod_polder_865=polder.variables['AOD865'][0:350,1652:2148]
            
            aaod_polder=polder.variables['AAOD565'][0:350,1652:2148]
            aaod_polder_443=polder.variables['AAOD443'][0:350,1652:2148]
            aaod_polder_865=polder.variables['AAOD865'][0:350,1652:2148]
            
            ssa_polder=polder.variables['SSA565'][0:350,1652:2148]
            ssa_polder_443=polder.variables['SSA443'][0:350,1652:2148]
            ssa_polder_865=polder.variables['SSA865'][0:350,1652:2148]
            
            ae_polder=polder.variables['AExp'][0:350,1652:2148]
            
            latitude_polder=polder.variables['Latitude'][0:350,1652:2148]
            longitude_polder=polder.variables['Longitude'][0:350,1652:2148]
            latitude_polder=latitude_polder.reshape((latitude_polder.shape[0]*latitude_polder.shape[1]))
            longitude_polder=longitude_polder.reshape((longitude_polder.shape[0]*longitude_polder.shape[1]))
            
            aod_polder=aod_polder.reshape((aod_polder.shape[0]*aod_polder.shape[1]))
            aod_polder_443=aod_polder_443.reshape((aod_polder_443.shape[0]*aod_polder_443.shape[1]))
            aod_polder_865=aod_polder_865.reshape((aod_polder_865.shape[0]*aod_polder_865.shape[1]))
            
            aaod_polder=aaod_polder.reshape((aaod_polder.shape[0]*aaod_polder.shape[1]))
            aaod_polder_443=aaod_polder_443.reshape((aaod_polder_443.shape[0]*aaod_polder_443.shape[1]))
            aaod_polder_865=aaod_polder_865.reshape((aaod_polder_865.shape[0]*aaod_polder_865.shape[1]))
            
            ssa_polder=ssa_polder.reshape((ssa_polder.shape[0]*ssa_polder.shape[1]))
            ssa_polder_443=ssa_polder_443.reshape((ssa_polder_443.shape[0]*ssa_polder_443.shape[1]))
            ssa_polder_865=ssa_polder_865.reshape((ssa_polder_865.shape[0]*ssa_polder_865.shape[1]))
            
            ae_polder=ae_polder.reshape((ae_polder.shape[0]*ae_polder.shape[1]))
            
            
            
            # dimensions
            npixel = latitude_polder[np.logical_not(aod_polder.mask)].size
            print(npixel)
            if npixel<10:
                continue
            #Creating netcdf4
            
            cso_file='CSO_POLDER_'+file_name[16:24]+'.nc'
            valid_files.append(cso_file)
            print(cso_file)
            root_grp = Dataset(cso_path+cso_file, 'w', format='NETCDF4')
            root_grp.description = 'POLDER observation following CSO estructure'
            
            
            
            root_grp.createDimension('pixel', npixel)
            root_grp.createDimension('retr', 1)
            root_grp.createDimension('corner', 4)
            # variables
            pixel = root_grp.createVariable('pixel', 'i4', ('pixel',),fill_value=-2147483647)
            retr = root_grp.createVariable('retr', 'i4', ('retr',))
            longitude = root_grp.createVariable('longitude', 'f4', ('pixel',),fill_value=3.4028235e+38)
            longitude_bounds = root_grp.createVariable('longitude_bounds', 'f4', ('pixel','corner',),fill_value=3.4028235e+38)
            latitude = root_grp.createVariable('latitude', 'f4', ('pixel',),fill_value=3.4028235e+38)
            latitude_bounds = root_grp.createVariable('latitude_bounds', 'f4', ('pixel','corner',),fill_value=3.4028235e+38)
            AOD565 = root_grp.createVariable('AOD565', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
            AOD865 = root_grp.createVariable('AOD865', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
            AOD443 = root_grp.createVariable('AOD443', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
            AAOD565 = root_grp.createVariable('AAOD565', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
            AAOD865 = root_grp.createVariable('AAOD865', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
            AAOD443 = root_grp.createVariable('AAOD443', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
            SSA565 = root_grp.createVariable('SSA565', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
            SSA865 = root_grp.createVariable('SSA865', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
            SSA443 = root_grp.createVariable('SSA443', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
            AExp = root_grp.createVariable('AExp', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
            
            obs_error_AOD565 = root_grp.createVariable('obs_error_AOD565', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
            obs_error_AOD865 = root_grp.createVariable('obs_error_AOD865', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
            obs_error_AOD443 = root_grp.createVariable('obs_error_AOD443', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
            obs_error_AAOD565 = root_grp.createVariable('obs_error_AAOD565', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
            obs_error_AAOD865 = root_grp.createVariable('obs_error_AAOD865', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
            obs_error_AAOD443 = root_grp.createVariable('obs_error_AAOD443', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
            obs_error_SSA565 = root_grp.createVariable('obs_error_SSA565', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
            obs_error_SSA865 = root_grp.createVariable('obs_error_SSA865', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
            obs_error_SSA443 = root_grp.createVariable('obs_error_SSA443', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
            obs_error_AExp = root_grp.createVariable('obs_error_AExp', 'f4', ('pixel', 'retr',),fill_value=3.4028235e+38)
            
            
            #Attributes
            pixel.long_name= 'pixel index'
            pixel.units= '1'
           
            longitude.long_name= 'pixel center longitude'
            longitude.units= 'degrees_east'
            longitude.standard_name= 'longitude'
            longitude.valid_min: -180.0
            longitude.valid_max: 180.0
            
            longitude_bounds.long_name= 'longitude_bounds'
            longitude_bounds.units= 'degrees_east'
            
            latitude.long_name= 'pixel center latitude'
            latitude.units= 'degrees_north'
            latitude.standard_name= 'latitude'
            latitude.valid_min: -90.0
            latitude.valid_max: 90.0
            
            latitude_bounds.long_name= 'latitude_bounds'
            latitude_bounds.units= 'degrees_north'
            
            AOD565.long_name= 'aod 565 nm'
            AOD565.units= '1'
            AOD565.standard_name= 'aod_565'
            
            AOD443.long_name= 'aod 443 nm'
            AOD443.units= '1'
            AOD443.standard_name= 'aod_443'
            
            AOD865.long_name= 'aod 865 nm'
            AOD865.units= '1'
            AOD865.standard_name= 'aod_865'
            
            AAOD565.long_name= 'aod 565 nm'
            AAOD565.units= '1'
            AAOD565.standard_name= 'aod_565'
            
            AAOD443.long_name= 'aod 443 nm'
            AAOD443.units= '1'
            AAOD443.standard_name= 'aod_443'
            
            AAOD865.long_name= 'aod 865 nm'
            AAOD865.units= '1'
            AAOD865.standard_name= 'aod_865'
            
            SSA565.long_name= 'ssa 565 nm'
            SSA565.units= '1'
            SSA565.standard_name= 'ssa_565'
            
            SSA443.long_name= 'ssa 443 nm'
            SSA443.units= '1'
            SSA443.standard_name= 'ssa_443'
            
            SSA865.long_name= 'ssa 865 nm'
            SSA865.units= '1'
            SSA865.standard_name= 'ssa_865'
            
            
            AExp.long_name= 'Angstrom Exponent 443-865 nm'
            AExp.units= '1'
            AExp.standard_name= 'aexp'
            
            #===Observation errors==
            
            obs_error_AOD565.long_name= 'Observation Error aod 565 nm'
            obs_error_AOD565.units= '1'
            obs_error_AOD565.standard_name= 'obs_error_aod_565'
            
            obs_error_AOD443.long_name= 'Observation Error aod 443 nm'
            obs_error_AOD443.units= '1'
            obs_error_AOD443.standard_name= 'obs_error_aod_443'
            
            obs_error_AOD865.long_name= 'Observation Error aod 865 nm'
            obs_error_AOD865.units= '1'
            obs_error_AOD865.standard_name= 'obs_error_aod_865'
            
            obs_error_AAOD565.long_name= 'Observation Error aod 565 nm'
            obs_error_AAOD565.units= '1'
            obs_error_AAOD565.standard_name= 'obs_error_aod_565'
            
            obs_error_AAOD443.long_name= 'Observation Error aod 443 nm'
            obs_error_AAOD443.units= '1'
            obs_error_AAOD443.standard_name= 'obs_error_aod_443'
            
            obs_error_AAOD865.long_name= 'Observation Error aod 865 nm'
            obs_error_AAOD865.units= '1'
            obs_error_AAOD865.standard_name= 'obs_error_aod_865'
            
            obs_error_SSA565.long_name= 'Observation Error ssa 565 nm'
            obs_error_SSA565.units= '1'
            obs_error_SSA565.standard_name= 'obs_error_ssa_565'
            
            obs_error_SSA443.long_name= 'Observation Error ssa 443 nm'
            obs_error_SSA443.units= '1'
            obs_error_SSA443.standard_name= 'obs_error_ssa_443'
            
            obs_error_SSA865.long_name= 'Observation Error ssa 865 nm'
            obs_error_SSA865.units= '1'
            obs_error_SSA865.standard_name= 'obs_error_ssa_865'
            
            
            obs_error_AExp.long_name= 'Observation Error Angstrom Exponent 443-865 nm'
            obs_error_AExp.units= '1'
            obs_error_AExp.standard_name= 'obs_error_aexp'
                        
            #values
            pixel[:]=list(range(1,npixel+1))
            retr[:]=1
            longitude[:]=longitude_polder[np.logical_not(aod_polder.mask)]
            longitude_bounds_aux=np.zeros((npixel,4))
            longitude_bounds_aux[:,0]=longitude_polder[np.logical_not(aod_polder.mask)] -0.005
            longitude_bounds_aux[:,1]=longitude_polder[np.logical_not(aod_polder.mask)] +0.005
            longitude_bounds_aux[:,2]=longitude_polder[np.logical_not(aod_polder.mask)] +0.005
            longitude_bounds_aux[:,3]=longitude_polder[np.logical_not(aod_polder.mask)] -0.005
            longitude_bounds[:]=longitude_bounds_aux
            latitude[:]=latitude_polder[np.logical_not(aod_polder.mask)]
            latitude_bounds_aux=np.zeros((npixel,4))
            latitude_bounds_aux[:,0]=latitude_polder[np.logical_not(aod_polder.mask)] -0.005
            latitude_bounds_aux[:,1]=latitude_polder[np.logical_not(aod_polder.mask)] -0.005
            latitude_bounds_aux[:,2]=latitude_polder[np.logical_not(aod_polder.mask)] +0.005
            latitude_bounds_aux[:,3]=latitude_polder[np.logical_not(aod_polder.mask)] +0.005
            latitude_bounds[:]=latitude_bounds_aux
            
            AOD565[:]=aod_polder[np.logical_not(aod_polder.mask)]
            AOD443[:]=aod_polder_443[np.logical_not(aod_polder.mask)]
            AOD865[:]=aod_polder_865[np.logical_not(aod_polder.mask)]
            
            AAOD565[:]=aaod_polder[np.logical_not(aod_polder.mask)]
            AAOD443[:]=aaod_polder_443[np.logical_not(aod_polder.mask)]
            AAOD865[:]=aaod_polder_865[np.logical_not(aod_polder.mask)]
            
            SSA565[:]=ssa_polder[np.logical_not(aod_polder.mask)]
            SSA443[:]=ssa_polder_443[np.logical_not(aod_polder.mask)]
            SSA865[:]=ssa_polder_865[np.logical_not(aod_polder.mask)]
            
            AExp[:]=ae_polder[np.logical_not(aod_polder.mask)]
            
            
            for i in range(len(AOD565)):
                if AOD565[i]<=0.05:
                  aux_sd=AOD565_SD[0]
                elif 0.05<AOD565[i]<=0.1:
                  aux_sd=AOD565_SD[1]  
                elif 0.1<AOD565[i]<=0.2:
                  aux_sd=AOD565_SD[2]  
                elif 0.2<AOD565[i]<=0.3:
                  aux_sd=AOD565_SD[3]  
                elif 0.3<AOD565[i]<=0.5:
                  aux_sd=AOD565_SD[4]  
                elif 0.5<AOD565[i]:
                  aux_sd=AOD565_SD[5]  
                obs_error_AOD565[i]= AOD565[i]*aux_sd
            
            for i in range(len(AOD865)):
                if AOD865[i]<=0.05:
                  aux_sd=AOD865_SD[0]
                elif 0.05<AOD865[i]<=0.1:
                  aux_sd=AOD865_SD[1]  
                elif 0.1<AOD865[i]<=0.2:
                  aux_sd=AOD865_SD[2]  
                elif 0.2<AOD865[i]<=0.3:
                  aux_sd=AOD865_SD[3]  
                elif 0.3<AOD865[i]<=0.5:
                  aux_sd=AOD865_SD[4]  
                elif 0.5<AOD865[i]:
                  aux_sd=AOD865_SD[5]  
                obs_error_AOD865[i]= AOD865[i]*aux_sd
            
            for i in range(len(AOD443)):
                if AOD443[i]<=0.05:
                  aux_sd=AOD443_SD[0]
                elif 0.05<AOD443[i]<=0.1:
                  aux_sd=AOD443_SD[1]  
                elif 0.1<AOD443[i]<=0.2:
                  aux_sd=AOD443_SD[2]  
                elif 0.2<AOD443[i]<=0.3:
                  aux_sd=AOD443_SD[3]  
                elif 0.3<AOD443[i]<=0.5:
                  aux_sd=AOD443_SD[4]  
                elif 0.5<AOD443[i]:
                  aux_sd=AOD443_SD[5]  
                obs_error_AOD443[i]= AOD443[i]*aux_sd
            
            for i in range(len(AAOD565)):
                if AAOD565[i]<=0.05:
                  aux_sd=AAOD565_SD[0]
                elif 0.05<AAOD565[i]<=0.1:
                  aux_sd=AAOD565_SD[1]  
                elif 0.1<AAOD565[i]<=0.2:
                  aux_sd=AAOD565_SD[2]  
                elif 0.2<AAOD565[i]<=0.3:
                  aux_sd=AAOD565_SD[3]  
                elif 0.3<AAOD565[i]<=0.5:
                  aux_sd=AAOD565_SD[4]  
                elif 0.5<AAOD565[i]:
                  aux_sd=AAOD565_SD[5]  
                obs_error_AAOD565[i]= AAOD565[i]*aux_sd
            
            
            for i in range(len(AAOD865)):
                if AAOD865[i]<=0.05:
                  aux_sd=AAOD865_SD[0]
                elif 0.05<AAOD865[i]<=0.1:
                  aux_sd=AAOD865_SD[1]  
                elif 0.1<AAOD865[i]<=0.2:
                  aux_sd=AAOD865_SD[2]  
                elif 0.2<AAOD865[i]<=0.3:
                  aux_sd=AAOD865_SD[3]  
                elif 0.3<AAOD865[i]<=0.5:
                  aux_sd=AAOD865_SD[4]  
                elif 0.5<AAOD865[i]:
                  aux_sd=AAOD865_SD[5]  
                obs_error_AAOD865[i]= AAOD865[i]*aux_sd
            
            for i in range(len(AAOD443)):
                if AAOD443[i]<=0.05:
                  aux_sd=AAOD443_SD[0]
                elif 0.05<AAOD443[i]<=0.1:
                  aux_sd=AAOD443_SD[1]  
                elif 0.1<AAOD443[i]<=0.2:
                  aux_sd=AAOD443_SD[2]  
                elif 0.2<AAOD443[i]<=0.3:
                  aux_sd=AAOD443_SD[3]  
                elif 0.3<AAOD443[i]<=0.5:
                  aux_sd=AAOD443_SD[4]  
                elif 0.5<AAOD443[i]:
                  aux_sd=AAOD443_SD[5]  
                obs_error_AAOD443[i]= AAOD443[i]*aux_sd
            
            for i in range(len(SSA565)):
                if SSA565[i]<=0.05:
                  aux_sd=SSA565_SD[0]
                elif 0.05<SSA565[i]<=0.1:
                  aux_sd=SSA565_SD[1]  
                elif 0.1<SSA565[i]<=0.2:
                  aux_sd=SSA565_SD[2]  
                elif 0.2<SSA565[i]<=0.3:
                  aux_sd=SSA565_SD[3]  
                elif 0.3<SSA565[i]<=0.5:
                  aux_sd=SSA565_SD[4]  
                elif 0.5<SSA565[i]:
                  aux_sd=SSA565_SD[5]  
                obs_error_SSA565[i]= SSA565[i]*aux_sd
               
            for i in range(len(SSA865)):
                if SSA865[i]<=0.05:
                  aux_sd=SSA865_SD[0]
                elif 0.05<SSA865[i]<=0.1:
                  aux_sd=SSA865_SD[1]  
                elif 0.1<SSA865[i]<=0.2:
                  aux_sd=SSA865_SD[2]  
                elif 0.2<SSA865[i]<=0.3:
                  aux_sd=SSA865_SD[3]  
                elif 0.3<SSA865[i]<=0.5:
                  aux_sd=SSA865_SD[4]  
                elif 0.5<SSA865[i]:
                  aux_sd=SSA865_SD[5]  
                obs_error_SSA865[i]= SSA865[i]*aux_sd
            
            for i in range(len(SSA443)):
                if SSA443[i]<=0.05:
                  aux_sd=SSA443_SD[0]
                elif 0.05<SSA443[i]<=0.1:
                  aux_sd=SSA443_SD[1]  
                elif 0.1<SSA443[i]<=0.2:
                  aux_sd=SSA443_SD[2]  
                elif 0.2<SSA443[i]<=0.3:
                  aux_sd=SSA443_SD[3]  
                elif 0.3<SSA443[i]<=0.5:
                  aux_sd=SSA443_SD[4]  
                elif 0.5<SSA443[i]:
                  aux_sd=SSA443_SD[5]  
                obs_error_SSA443[i]= SSA443[i]*aux_sd
                
            for i in range(len(AExp)):
                if AExp[i]<=0.05:
                  aux_sd=AExp_SD[0]
                elif 0.05<AExp[i]<=0.1:
                  aux_sd=AExp_SD[1]  
                elif 0.1<AExp[i]<=0.2:
                  aux_sd=AExp_SD[2]  
                elif 0.2<AExp[i]<=0.3:
                  aux_sd=AExp_SD[3]  
                elif 0.3<AExp[i]<=0.5:
                  aux_sd=AExp_SD[4]  
                elif 0.5<AExp[i]:
                  aux_sd=AExp_SD[5]  
                obs_error_AExp[i]= AExp[i]*aux_sd
            
            root_grp.close()
            polder.close()
        
    
    
    for filename in valid_files:
        if filename[-3:]=='.nc':
            list_CSO['filename'].extend(['CSO/'+filename])
            list_CSO['start_time'].extend([filename[11:15]+'-'+filename[15:17]+'-'+filename[17:19]+'T13:30:00.00'])
            list_CSO['end_time'].extend([filename[11:15]+'-'+filename[15:17]+'-'+filename[17:19]+'T15:30:00.00'])
            list_CSO['orbit'].extend([filename[11:19]])
pd.DataFrame.from_dict(list_CSO).to_csv('/Users/santiago/Documents/POLDER/2008/listing_POLDER.csv',sep=';',index=False)
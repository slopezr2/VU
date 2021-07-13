import pandas as pd
import numpy as np
from os import listdir
from os.path import isfile, join
from datetime import datetime,timedelta
from netCDF4 import Dataset

def read_EMEP(filename):
    with open(filename) as f:
        lines = f.readlines()
    
    station={}
    for i in lines:
        if i[0:16]=='Station latitude':
            latitude=[]
            for j in range(2,len(i)):
                if i[-j]==' ':
                    break
                latitude.append(i[-j])
    latitude.reverse()
    aux=''
    latitude=aux.join(latitude)
    latitude=float(latitude)
    station['latitude']=latitude
    
    for i in lines:
        if i[0:17]=='Station longitude':
            longitude=[]
            for j in range(2,len(i)):
                if i[-j]==' ':
                    break
                longitude.append(i[-j])
    longitude.reverse()
    aux=''
    longitude=aux.join(longitude)
    longitude=float(longitude)
    station['longitude']=longitude
    
    for i in lines:
        if i[0:16]=='Station altitude':
            altitude=[]
            for j in range(4,len(i)):
                if i[-j]==' ':
                    break
                altitude.append(i[-j])
    altitude.reverse()
    aux=''
    altitude=aux.join(altitude)
    altitude=float(altitude)
    station['altitude']=altitude
    
    
    for i in lines:
        if i[0:12]=='Station code':
            code=[]
            for j in range(2,len(i)):
                if i[-j]==' ':
                    break
                code.append(i[-j])
    code.reverse()
    aux=''
    code=aux.join(code)
    station['code']=code
    
    j=1
    for i in lines:
        if i[0:9]=='starttime':
            break
        j=j+1    
    
    #start_time=[]
    #PM10=[]
    #PM25=[]

    emep = pd.read_csv(filename, skiprows=j-1,delimiter=r"\s+")
    
    
   # emep['starttime']=datetime.strptime('2008-1-1 00:00:00', '%Y-%m-%d %H:%M:%S')+timedelta(hours=emep['starttime'])
    emep['starttime']=[datetime.strptime('2008-1-1 00:00:00', '%Y-%m-%d %H:%M:%S')+timedelta(days=np.floor(fs)) for fs in emep['starttime']]
    emep.set_index('starttime',inplace=True)
    emep.drop('endtime',inplace=True, axis=1)
    station['values']=emep
    
    
    
    

    return station


path_EMEP='/Users/santiago/Documents/EMEP/PM/'
cso_path='/Users/santiago/Documents/EMEP/CSO/'

onlyfiles = [f for f in listdir(path_EMEP) if isfile(join(path_EMEP, f))]
onlyfiles.sort()
EMEP={}
for file in onlyfiles:
    
    if file[8:12]=='2008':
        
        aux=read_EMEP(path_EMEP+'/'+file)
        EMEP[file[0:7]]=aux

template=pd.read_csv('/Users/santiago/scratch/CAMS-61/WP6130/CSO-data/S5p/listing-NO2-CAMS.csv',delimiter=';')
list_CSO={name: [] for name in list(template.columns)}
datelist = pd.date_range(datetime.strptime('2008-01-01 00:00:00','%Y-%m-%d %H:%M:%S'), periods=366,freq='D').tolist()
valid_files=[]
for fecha in datelist:
    npixel=0
    sites_fecha=[]
    aux_lon=[]
    aux_lat=[]
    aux_pm10=[]
    for site in EMEP.keys():
        if ((fecha in EMEP[site]['values'].index) and ('PM10' in EMEP[site]['values'].columns)):
            npixel=npixel+1
            sites_fecha.append(site)
            aux=EMEP[site]['longitude']
            aux_lon.append(aux)
            aux=EMEP[site]['latitude']
            aux_lat.append(aux)
            aux=EMEP[site]['values']['PM10'][EMEP[site]['values'].index==fecha].values
            aux_pm10.append(aux[0])
    if np.mean(aux_pm10)==999.0:
        continue
    if npixel>0:
            cso_file='CSO_EMEP_PM10_'+fecha.strftime('%Y')+fecha.strftime('%m')+fecha.strftime('%d')+'.nc'
            
            valid_files.append(cso_file)
            root_grp = Dataset(cso_path+cso_file, 'w', format='NETCDF4')
            root_grp.description = 'EMEP observation following CSO estructure'
            
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
            pm10 = root_grp.createVariable('PM10', 'f4', ('pixel', 'retr',),fill_value=999.0)
            
            
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
            
            pm10.long_name= 'PM10 daily average'
            pm10.units= 'ug/m3'
            pm10.standard_name= 'pm10'
            
            #values
            pixel[:]=list(range(1,npixel+1))
            retr[:]=1
            resolution=0.01
            longitude[:]=np.array(aux_lon)
            longitude_bounds_aux=np.zeros((npixel,4))
            aux=np.array(aux_lon)
            longitude_bounds_aux[:,0]=aux -resolution/2
            longitude_bounds_aux[:,1]=aux +resolution/2
            longitude_bounds_aux[:,2]=aux +resolution/2
            longitude_bounds_aux[:,3]=aux -resolution/2
            longitude_bounds[:]=longitude_bounds_aux
            latitude[:]=np.array(aux_lat)
            latitude_bounds_aux=np.zeros((npixel,4))
            aux=np.array(aux_lat)
            latitude_bounds_aux[:,0]=aux -resolution/2
            latitude_bounds_aux[:,1]=aux -resolution/2
            latitude_bounds_aux[:,2]=aux +resolution/2
            latitude_bounds_aux[:,3]=aux +resolution/2
            latitude_bounds[:]=latitude_bounds_aux
            
            pm10[:]=np.array(aux_pm10)
            print(cso_file)
            
            root_grp.close()

        
    
    
for filename in valid_files:
    if filename[-3:]=='.nc':
        list_CSO['filename'].extend(['CSO/'+filename])
        list_CSO['start_time'].extend([filename[14:18]+'-'+filename[18:20]+'-'+filename[20:22]+'T'+'00:00:00.00'])
        list_CSO['end_time'].extend([filename[14:18]+'-'+filename[18:20]+'-'+filename[20:22]+'T'+'23:59:00.00'])
        list_CSO['orbit'].extend([filename[14:22]])
pd.DataFrame.from_dict(list_CSO).to_csv('/Users/santiago/Documents/EMEP/listing_EMEP_PM10.csv',sep=';',index=False)        




valid_files=[]
for fecha in datelist:
    npixel=0
    sites_fecha=[]
    aux_lon=[]
    aux_lat=[]
    aux_pm25=[]
    for site in EMEP.keys():
        if ((fecha in EMEP[site]['values'].index) and ('PM2.5' in EMEP[site]['values'].columns)):
            npixel=npixel+1
            sites_fecha.append(site)
            aux=EMEP[site]['longitude']
            aux_lon.append(aux)
            aux=EMEP[site]['latitude']
            aux_lat.append(aux)
            aux=EMEP[site]['values']['PM2.5'][EMEP[site]['values'].index==fecha].values
            aux_pm25.append(aux[0])
    if np.mean(aux_pm25)==999.0:
        continue
    if npixel>0:
            cso_file='CSO_EMEP_PM25_'+fecha.strftime('%Y')+fecha.strftime('%m')+fecha.strftime('%d')+'.nc'
            
            valid_files.append(cso_file)
            root_grp = Dataset(cso_path+cso_file, 'w', format='NETCDF4')
            root_grp.description = 'EMEP observation following CSO estructure'
            
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
            pm25 = root_grp.createVariable('PM25', 'f4', ('pixel', 'retr',),fill_value=999.0)
            
            
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
            
            pm25.long_name= 'PM25 daily average'
            pm25.units= 'ug/m3'
            pm25.standard_name= 'PM25'
            
            #values
            pixel[:]=list(range(1,npixel+1))
            retr[:]=1
            resolution=0.01
            longitude[:]=np.array(aux_lon)
            longitude_bounds_aux=np.zeros((npixel,4))
            aux=np.array(aux_lon)
            longitude_bounds_aux[:,0]=aux -resolution/2
            longitude_bounds_aux[:,1]=aux +resolution/2
            longitude_bounds_aux[:,2]=aux +resolution/2
            longitude_bounds_aux[:,3]=aux -resolution/2
            longitude_bounds[:]=longitude_bounds_aux
            latitude[:]=np.array(aux_lat)
            latitude_bounds_aux=np.zeros((npixel,4))
            aux=np.array(aux_lat)
            latitude_bounds_aux[:,0]=aux -resolution/2
            latitude_bounds_aux[:,1]=aux -resolution/2
            latitude_bounds_aux[:,2]=aux +resolution/2
            latitude_bounds_aux[:,3]=aux +resolution/2
            latitude_bounds[:]=latitude_bounds_aux
            
            pm25[:]=np.array(aux_pm25)
            print(cso_file)
            
            root_grp.close()

        
    
list_CSO={name: [] for name in list(template.columns)}    
for filename in valid_files:
    if filename[-3:]=='.nc':
        list_CSO['filename'].extend(['CSO/'+filename])
        list_CSO['start_time'].extend([filename[14:18]+'-'+filename[18:20]+'-'+filename[20:22]+'T'+'00:00:00.00'])
        list_CSO['end_time'].extend([filename[14:18]+'-'+filename[18:20]+'-'+filename[20:22]+'T'+'23:59:00.00'])
        list_CSO['orbit'].extend([filename[14:22]])
pd.DataFrame.from_dict(list_CSO).to_csv('/Users/santiago/Documents/EMEP/listing_EMEP_PM25.csv',sep=';',index=False)        
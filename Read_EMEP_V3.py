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
    
    emep['starttime']=[datetime.strptime('2008-1-1 00:00:00', '%Y-%m-%d %H:%M:%S')+timedelta(days=np.floor(fs)) for fs in emep['starttime']]
    emep.set_index('starttime',inplace=True)
    emep.drop('endtime',inplace=True, axis=1)
    station['values']=emep
    station['values'] = station['values'].loc[station['values'][station['values'].keys()[0]]<900].groupby(['starttime']).mean()
    
    
    

    return station


path_EMEP_General='/Users/santiago/Documents/EMEP/'


species=['NH3','SO2','NO2']
for specie in species:
    cso_path ='/Users/santiago/Documents/EMEP/CSO/'+specie
    path_EMEP = '/Users/santiago/Documents/EMEP/'+specie
    onlyfiles = [f for f in listdir(path_EMEP) if isfile(join(path_EMEP, f))]
    onlyfiles.sort()
    EMEP = {}
    for file in onlyfiles:
        
        if file[8:12]=='2008':
            aux = read_EMEP(path_EMEP+'/'+file)
            EMEP[file[0:7]] = aux
    
    template = pd.read_csv('/Users/santiago/scratch/CAMS-61/WP6130/CSO-data/S5p/listing-NO2-CAMS.csv',delimiter=';')
    list_CSO = {name: [] for name in list(template.columns)}
    datelist = pd.date_range(datetime.strptime('2008-01-01 00:00:00','%Y-%m-%d %H:%M:%S'), periods=366,freq='D').tolist()
    valid_files = []
    for fecha in datelist:
        npixel = 0
        sites_fecha = []
        aux_lon = []
        aux_lat = []
        aux_emep = []
        for site in EMEP.keys():
            if ((fecha in EMEP[site]['values'].index) and (specie in EMEP[site]['values'].columns)):
                npixel=npixel+1
                sites_fecha.append(site)
                aux=EMEP[site]['longitude']
                aux_lon.append(aux)
                aux=EMEP[site]['latitude']
                aux_lat.append(aux)
                aux=EMEP[site]['values'][specie][EMEP[site]['values'].index==fecha].values.copy()
                aux_emep.append(aux[0])
                #EMEP[site]['values'][specie][EMEP[site]['values'][specie]>900]=np.NaN
        aux_emep = np.array(aux_emep)
        aux_emep[aux_emep>900] = 999.0
        if np.mean(aux_emep)==999.0:
            continue
        if npixel>0:
                cso_file='CSO_EMEP_'+specie+'_'+fecha.strftime('%Y')+fecha.strftime('%m')+fecha.strftime('%d')+'.nc'
                
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
                pm10 = root_grp.createVariable(specie, 'f4', ('pixel', 'retr',),fill_value=999.0)
                
                
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
                
                pm10.long_name= specie + ' daily average'
                pm10.units= 'ug/m3'
                pm10.standard_name= specie
                
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
                
                pm10[:]=np.array(aux_emep)
                print(cso_file)
                
                root_grp.close()
    
            
        
        
    list_CSO={name: [] for name in list(template.columns)}    
    for filename in valid_files:
        if filename[-3:]=='.nc':
            hours_day=['00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23']
            for j in range(23):
                list_CSO['filename'].extend(['CSO/'+filename])
                list_CSO['start_time'].extend([filename[14:18]+'-'+filename[18:20]+'-'+filename[20:22]+'T'+hours_day[j]+':00:00.00'])
                list_CSO['end_time'].extend([filename[14:18]+'-'+filename[18:20]+'-'+filename[20:22]+'T'+hours_day[j+1]+':00:00.00'])
                list_CSO['orbit'].extend([filename[14:22]+hours_day[j]])
    pd.DataFrame.from_dict(list_CSO).to_csv('/Users/santiago/Documents/EMEP/listing_EMEP_'+specie+'.csv',sep=';',index=False)        
    

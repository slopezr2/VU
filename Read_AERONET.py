import pandas as pd
import numpy as np
from os import listdir
from os.path import isfile, join
from datetime import datetime
from netCDF4 import Dataset

def read_aeronet(filename):
    """Read a given AERONET AOT data file, and return it as a dataframe.
    
    This returns a DataFrame containing the AERONET data, with the index
    set to the timestamp of the AERONET observations. Rows or columns
    consisting entirely of missing data are removed. All other columns
    are left as-is.
    """
    dateparse = lambda x: pd.datetime.strptime(x, "%d:%m:%Y %H:%M:%S")
    aeronet = pd.read_csv(filename, skiprows=6, na_values=['N/A'],
                          parse_dates={'times':[0,1]},
                          date_parser=dateparse)

    aeronet = aeronet.set_index('times')
    #del aeronet['Julian_Day']
    
    # Drop any rows that are all NaN and any cols that are all NaN
    # & then sort by the index
    an = (aeronet.dropna(axis=1, how='all')
                .dropna(axis=0, how='all')
                .rename(columns={'Last_Processing_Date(dd/mm/yyyy)': 'Last_Processing_Date'})
                .sort_index())

    return an


path_aeronet='/Users/santiago/Documents/AERONET/AOD_Level20_All_Points_V3/AOD/AOD20/ALL_POINTS'
cso_path='/Users/santiago/Documents/AERONET/CSO/'
sites_aeronet=pd.read_csv('/Users/santiago/Documents/AERONET/List_AERONET_sites.csv')
sites_domain=sites_aeronet[(sites_aeronet['Longitude(decimal_degrees)']>-15) & (sites_aeronet['Longitude(decimal_degrees)']<35) & (sites_aeronet['Latitude(decimal_degrees)']>35) & (sites_aeronet['Latitude(decimal_degrees)']<70  )]
sites_domain_2008=['ATHENS-NOA', 'Arcachon', 'Baneasa', 'Cabo_Raso', 'Ersa', 'Helsinki', 'Hyytiala', 'Kuopio', 'Malaga', 'Palgrunden', 'Saint_Mandrier', 'Valladolid_Sci', 'Wytham_Woods', 'Xanthi']

onlyfiles = [f for f in listdir(path_aeronet) if isfile(join(path_aeronet, f))]
onlyfiles.sort()
aeronet={}
for file in onlyfiles:
    
    if file[18:-6] in sites_domain['Site_Name'].values:
        
        aux=read_aeronet(path_aeronet+'/'+file)
        #if (aux.index[0]>=datetime.strptime('2008-01-01 00:00:00','%Y-%m-%d %H:%M:%S')) and (aux.index[0]<=datetime.strptime('2008-12-31 23:59:59','%Y-%m-%d %H:%M:%S')):
        aux=aux[(aux.index>=datetime.strptime('2008-01-01 00:00:00','%Y-%m-%d %H:%M:%S')) & (aux.index<=datetime.strptime('2008-12-31 23:59:59','%Y-%m-%d %H:%M:%S'))]
        if len(aux.index)==0:
            continue
        aux = aux.resample('60T', label='right').mean()
        aux=aux.dropna(axis=1, how='all').dropna(axis=0, how='all')
        aeronet[file[18:-6]]=aux


template=pd.read_csv('/Users/santiago/scratch/CAMS-61/WP6130/CSO-data/S5p/listing-NO2-CAMS.csv',delimiter=';')
list_CSO={name: [] for name in list(template.columns)}
datelist = pd.date_range(datetime.strptime('2008-01-01 00:00:00','%Y-%m-%d %H:%M:%S'), periods=366*24,freq='60min').tolist()
valid_files=[]
for fecha in datelist:
    npixel=0
    sites_fecha=[]
    aux_lon=[]
    aux_lat=[]
    aux_aod500=[]
    aux_aod550=[]
    aux_aod440=[]
    aux_aod870=[]
    aux_aex=[]
    for site in aeronet.keys():
        if fecha in aeronet[site].index:
            npixel=npixel+1
            sites_fecha.append(site)
            aux=sites_domain['Longitude(decimal_degrees)'][sites_domain['Site_Name']==site].values
            aux_lon.append(aux[0])
            aux=sites_domain['Latitude(decimal_degrees)'][sites_domain['Site_Name']==site].values
            aux_lat.append(aux[0])
            aux=aeronet[site]['AOD_500nm'][aeronet[site].index==fecha].values
            if aux[0]<0:
                aux[0]=-999.0
            aux_aod500.append(aux[0])
            #aux=aeronet[site]['AOD_550nm'][aeronet[site].index==fecha].values
            #aux_aod550.append(aux[0])
            aux=aeronet[site]['AOD_440nm'][aeronet[site].index==fecha].values
            if aux[0]<0:
                aux[0]=-999.0
            aux_aod440.append(aux[0])
            aux=aeronet[site]['AOD_870nm'][aeronet[site].index==fecha].values
            if aux[0]<0:
                aux[0]=-999.0
            aux_aod870.append(aux[0])
            aux=aeronet[site]['440-870_Angstrom_Exponent'][aeronet[site].index==fecha].values
            if aux[0]<0:
                aux[0]=-999.0
            aux_aex.append(aux[0])
    if np.mean(aux_aod500)==-999.0:
        continue
    if npixel>0:
        
            cso_file='CSO_AERONET_'+fecha.strftime('%Y')+fecha.strftime('%m')+fecha.strftime('%d')+'_'+fecha.strftime('%H')+'.nc'
            
            valid_files.append(cso_file)
            root_grp = Dataset(cso_path+cso_file, 'w', format='NETCDF4')
            root_grp.description = 'AERONET observation following CSO estructure'
            
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
            AOD500 = root_grp.createVariable('AOD500', 'f4', ('pixel', 'retr',),fill_value=-999.0)
            #AOD550 = root_grp.createVariable('AOD550', 'f4', ('pixel', 'retr',),fill_value=-999.0)
            AOD440 = root_grp.createVariable('AOD440', 'f4', ('pixel', 'retr',),fill_value=-999.0)
            AOD870 = root_grp.createVariable('AOD870', 'f4', ('pixel', 'retr',),fill_value=-999.0)
            AExp = root_grp.createVariable('AExp', 'f4', ('pixel', 'retr',),fill_value=-999.0)
            
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
            
            AOD500.long_name= 'aod 500 nm'
            AOD500.units= '1'
            AOD500.standard_name= 'aod_500'
            
            #AOD550.long_name= 'aod 500 nm'
            #AOD550.units= '1'
            #AOD550.standard_name= 'aod_500'
            
            AOD440.long_name= 'aod 440 nm'
            AOD440.units= '1'
            AOD440.standard_name= 'aod_440'
            
            AOD870.long_name= 'aod 870 nm'
            AOD870.units= '1'
            AOD870.standard_name= 'aod_870'
            
           
            
            
            AExp.long_name= 'Angstrom Exponent 440-870 nm'
            AExp.units= '1'
            AExp.standard_name= 'aexp'
           
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
            
            AOD500[:]=np.array(aux_aod500)
            print(cso_file)
            print(aux_aod500)
            #AOD550[:]=aux_aod550
            AOD440[:]=np.array(aux_aod440)
            AOD870[:]=np.array(aux_aod870)
            AExp[:]=np.array(aux_aex)
            
            
            
            
            root_grp.close()

        
    
    
for filename in valid_files:
    if filename[-3:]=='.nc':
        list_CSO['filename'].extend(['CSO/'+filename])
        list_CSO['start_time'].extend([filename[12:16]+'-'+filename[16:18]+'-'+filename[18:20]+'T'+filename[21:23]+':00:00.00'])
        list_CSO['end_time'].extend([filename[12:16]+'-'+filename[16:18]+'-'+filename[18:20]+'T'+filename[21:23]+':59:00.00'])
        list_CSO['orbit'].extend([filename[12:23]])
pd.DataFrame.from_dict(list_CSO).to_csv('/Users/santiago/Documents/AERONET/listing_AERONET.csv',sep=';',index=False)        
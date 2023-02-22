#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main script to run the LTEKF-4D (Local Transform Ensemble Smoother) in LOTOS-EUROS
Author: Santiago Lopez-Restrepo, VU Amsterdam/TNO
"""
import sys
from datetime import datetime
from datetime import timedelta
import numpy as np
import pandas as pd
sys.path.append('./LE_communication/')
from LE_communication import LE_model
from LE_communication import IF_NOT_OK
from LETKF_4D import LETKF_4D
status = 0
error_message = ''


#==Create object LE_model to interact with the LEKF simulation
LE = LE_model()
#==Define path to run the LE simulation===
LE.run_path = '/scratch-shared/slr/projects/Smoother/'

#==Definition of the assimilation window==
Tb = 8		 #Background step. The total days of simulation.
Ts = 2 		 #Assimilation step. The last days of the background step where observations will be collected. This value guide the number of days using a constant dc factor.
if Tb%Ts!=0:
	status = 1
	error_message = 'The background step Tb is not divisible by the assimilation step Ts. Use a proper definition to avoid smoother issues.'
IF_NOT_OK(status,error_message)
Ta = Tb-Ts		#Analysis step. The initial days where the dc factor will be estimated.

#==Define id simulation name and dates of the==
LE.ID = 'Test_V9'
start_date = '2008-05-01 00:00:00'
Ncycles = 4			#Number of assimilation cycles.

#==Number of ensembles and localization radius==
LE.Nens = 32
radius = 1 #Degree
#==CSO instrument outputs to be assimilated==
CSO_instruments = ['polder-aod']

#==Using spatial correlation for perturbation factor===
spatial_corr = True

#============================== End of user definitions parameters ======================================================

#==Be sure that CSO_instruments is a list==
status = 0
error_message = ''
if type(CSO_instruments)==str:
	CSO_instruments = [CSO_instruments]
elif type(CSO_instruments)==list:
	pass
else:
	status = 1
	error_message = 'Incorrect format for the CSO instruments to be assimilated. CSO_instruments have to be a list of strings.'

#==Conver start date from string to datetime==
start_date_datetime = datetime.strptime(start_date, '%Y-%m-%d %H:%M:%S')

#==Write ID,and number of ensembles to the LE configuration files===
LE.write_ID()
LE.write_N_ensembles()

#===Init assimilation cycles===
for cycle in range(Ncycles):
	#==Init a new Smoother instance==
	LETKF = LETKF_4D(Nens=LE.Nens,Tb=Tb,Ts=Ts,Ta=Ta)
	flaq_obs = 0		#Detect whether there are observations in the cycle or not
	
	#==Run the new cycle==
	#Check if it is the first cycle, to generate random initial dc
	if cycle==0:
		#==Run the model one day to generate the template dc files==
		LE.write_read_dc_file(read_dc_file='F',path='${my.project.dir}/dc_smoother',model='LE_${run.id}_smoother')		#Option to LEKF reads the dc from files
		LE.write_restart(restart='F')
		new_start_date =  start_date_datetime+timedelta(days=Ts*cycle)
		LE.start_date = new_start_date.strftime('%Y-%m-%d %H:%M:%S')
		end_date_datetime = new_start_date+timedelta(days=1)-timedelta(hours=1)
		LE.end_date = end_date_datetime.strftime('%Y-%m-%d %H:%M:%S')
		LE.write_date()
		print('Running Initial one day simulation to get the DC files shape')
		LE.compile()
		#=Chance the .out and .err for each cycle
		LE.run_command = 'srun  --kill-on-bad-exit=1 --output=le_cycle_'+str(cycle)+'.run.out.%t  --error=le_cycle_'+str(cycle)+'.run.err.%t  ./lekf.x lekf.rc'
		LE.run()
		
		#==Generate the random initial dc files==
		#=Get the dates for the next new cycle =
		dates_cycle_new_cycle = pd.date_range(new_start_date,new_start_date+timedelta(days=Tb),freq='1d',inclusive='left')	#Here it is neccesary to take each day, to create the new DC files
		dates_cycle_new_cycle = [f.strftime('%Y%m%d') for f in dates_cycle_new_cycle] 	#Convert the dates to string format of LE files
		words = [dates_cycle_new_cycle[0]]
		words.extend(['dc','xa'])
		status = 0
		err_message = ''
		file_name = LE.list_output_files(words)
		if file_name==[] or len(file_name)>1:
			status = 1
			err_message = 'Problems geting the list of dc factors for initial cycle: '+str(cycle)+', ensemble: '+ens_str+', at day: '+date 
		else:
			file_name = file_name[0]
		IF_NOT_OK(status,error_message)
		#==Read a dc file to get the dimensions==
		aux_dim = LE.read_one_output(file_name=file_name,variable='dc',vector=False)[None,:,:,:]
		dc_dim = [24,1] #24 hours, and 1 hist
		dc_dim.extend(aux_dim.shape[1:]) #(/hours, hist, noise, lat, lon / )
		aux_dc = np.ones(dc_dim)
		for ens in range(1,LE.Nens+1):
			i_day_LETKF_bef = -1
			for idate,date in enumerate(dates_cycle_new_cycle):
				i_day_LETKF = idate // Ts
				ens_str = f'{ens:02d}'
				if i_day_LETKF_bef!=i_day_LETKF: #Actualize just for a new Ts day
					aux_dc[:,0,:,:,:] = LETKF.create_rand_ensemble_dc(new_mean=1,new_std=1,r=radius+1,auto_shape=False,nnoise=dc_dim[2],nlat=dc_dim[3],nlon=dc_dim[4],spatial_corr=spatial_corr) #The same value for all the hours of the day
					print('Writting a new initial background DC file for the ensemble: '+ens_str+' at day '+date)
				else:
					print('Writting the previous background DC file for the ensemble: '+ens_str+' at day '+date)
				filename = 'LE_'+LE.ID+'_smoother_xi_'+ens_str+'_dc_'+date+'_xa.nc'
				LE.write_one_dc_file(filename=filename,date=date,data=aux_dc )
				filename = 'LE_'+LE.ID+'_smoother_xi_'+ens_str+'_dc_'+date+'_xb_cycle_'+str(cycle)+'_.nc'
				LE.write_one_dc_file(filename=filename,date=date,data=aux_dc )
				i_day_LETKF_bef = i_day_LETKF
		LE.write_read_dc_file(read_dc_file='T',path='${my.project.dir}/dc_smoother',model='LE_${run.id}_smoother')		#Option to LEKF reads the dc from files
	else:
		LE.write_read_dc_file(read_dc_file='T',path='${my.project.dir}/dc_smoother',model='LE_${run.id}_smoother')		#Option to LEKF reads the dc from files
		LE.write_restart(restart='T') 	#Restart the simulation from the last concentration field
	
	#==Create a new cycle of simulations
	new_start_date =  start_date_datetime+timedelta(days=Ts*cycle)
	LE.start_date = new_start_date.strftime('%Y-%m-%d %H:%M:%S')
	end_date_datetime = new_start_date+timedelta(days=Tb)-timedelta(hours=1)
	LE.end_date=end_date_datetime.strftime('%Y-%m-%d %H:%M:%S')
	LE.write_date()
	print('Running Cycle '+str(cycle)+', from: '+LE.start_date+' to: '+LE.end_date)
	
	LE.compile()
	#=Chance the .out and .err for each cycle
	LE.run_command = 'srun  --kill-on-bad-exit=1 --output=le_cycle_'+str(cycle)+'.run.out.%t  --error=le_cycle_'+str(cycle)+'.run.err.%t  ./lekf.x lekf.rc'
	LE.run()
	

	#=Get the dates between the assimilation start and end of period=
	dates_cycle_dc = pd.date_range(new_start_date+timedelta(days=Ts),end_date_datetime,freq=str(Ts)+'d',inclusive='left')	#Is nor necessary read all days, because the dc factors chance each Ts days
	dates_cycle_dc = [f.strftime('%Y%m%d') for f in dates_cycle_dc] 	#Convert the dates to string format of LE files
	
	#=Get the dates of the observation period=
	dates_cycle_obs = pd.date_range(new_start_date+timedelta(days=Ta),end_date_datetime,freq='1d',inclusive='left')	#Read observations of each day
	dates_cycle_obs = [f.strftime('%Y%m%d') for f in dates_cycle_obs] 	#Convert the dates to string format of LE files
	
	for ens in range(1,LE.Nens+1):
		#==Read the dc, dc_lat, and dc_lon of each ensemble
		ens_str = f'{ens:02d}'
		for idate,date in enumerate(dates_cycle_dc):
			words = [date]
			words.extend(['dc','xi'+ens_str+'a'])
			status = 0
			err_message = ''
			file_name = LE.list_output_files(words)
			if file_name==[] or len(file_name)>1:
				status = 1
				err_message = 'Problems geting the list of dc factors for cycle: '+str(cycle)+', ensemble: '+ens_str+', at day: '+date 
			else:
				file_name = file_name[0]
			IF_NOT_OK(status,error_message)
			#==Read dc,lat, and lon values==
			if idate==0: #First day of the cycle
				LETKF.Ens[ens] = LE.read_one_output(file_name=file_name,variable='dc',vector=False)[None,:,:,:]
				LETKF.Ens_lat = LE.read_one_output(file_name=file_name,variable='latitude',vector=True)
				LETKF.Ens_lon= LE.read_one_output(file_name=file_name,variable='longitude',vector=True)
			else:
				LETKF.Ens[ens] = np.concatenate((LETKF.Ens[ens],LE.read_one_output(file_name=file_name,variable='dc',vector=False)[None,:,:,:]),axis=0)
			
			print('Readed DC cycle: '+str(cycle)+', ensemble: '+ens_str+', at day: '+date)
		for idate,date in enumerate(dates_cycle_obs):
			#==Read Y_ens,==
			
			for instrument in CSO_instruments:
				words = [date]
				words.extend([instrument,'xi'+ens_str+'a'])
				file_name = LE.list_output_files(words)
				if file_name==[]: #The are not observations for this day
					continue
				else:
					flaq_obs = 1
				if idate==0: #First day of the cycle
					LETKF.Y_Ens[ens] = LE.read_multiple_outputs(files=file_name,variable='ys',vector=True)
				else:
					LETKF.Y_Ens[ens] = np.concatenate((LETKF.Y_Ens[ens],LE.read_multiple_outputs(files=file_name,variable='ys',vector=True)))
				
				#==Read Y,Y_lat,Y_lon==
				if ens==1:		#Just read one time
					words = [date]
					words.extend([instrument,'data'])
					file_name = LE.list_output_files(words)
					if idate==0: #First day of the cycle
						LETKF.Y = LE.read_multiple_outputs(files=file_name,variable='yr',vector=True)
						LETKF.Y_error = LE.read_multiple_outputs(files=file_name,variable='obs_error',vector=True)*0.001
						LETKF.Y_lat = LE.read_multiple_outputs(files=file_name,variable='latitude',vector=True)
						LETKF.Y_lon = LE.read_multiple_outputs(files=file_name,variable='longitude',vector=True)
					else:
						LETKF.Y = np.concatenate((LETKF.Y,LE.read_multiple_outputs(files=file_name,variable='yr',vector=True)))
						LETKF.Y_error = np.concatenate((LETKF.Y_error,LE.read_multiple_outputs(files=file_name,variable='obs_error',vector=True)*0.001))
						LETKF.Y_lat = np.concatenate((LETKF.Y_lat,LE.read_multiple_outputs(files=file_name,variable='latitude',vector=True)))
						LETKF.Y_lon = np.concatenate((LETKF.Y_lon,LE.read_multiple_outputs(files=file_name,variable='longitude',vector=True)))
				
				print('Readed Observations from: ' +instrument+ ' in the cycle: '+str(cycle)+', ensemble: '+ens_str+', at day: '+date)
		
		#==Initialize the analysis ensemble===
		LETKF.Ens_analy[ens] = np.zeros(LETKF.Ens[ens].shape)
	#===Check wheter flaq_obs is equal to one to continue, in other case, there is nor observations available for the cicle
	if flaq_obs==1:
		for i in range(len(LETKF.Ens_lat)):
			for j in range(len(LETKF.Ens_lon)):
				LETKF.get_subdomain(lat_point=i,lon_point=j,radius=radius)
				if len(LETKF.sub_Y)==0: #Check wheter there is observations in the subdomain
					for ens in range(1,LE.Nens+1):
						LETKF.Ens_analy[ens][:,:,i,j] = LETKF.Ens[ens][:,:,i,j]
				else:
					LETKF.analysis_step()
					for ens in range(1,LE.Nens+1):
						LETKF.Ens_analy[ens][:,:,i,j] = LETKF.map_local_global(ens)
	else:
		#=Get the dates for the next new cycle =
		dates_cycle_new_cycle = pd.date_range(new_start_date+timedelta(days=Ts),end_date_datetime+timedelta(days=Ts),freq='1d',inclusive='left')	#Here it is neccesary to take each day, to create the new DC files
		dates_cycle_new_cycle = [f.strftime('%Y%m%d') for f in dates_cycle_new_cycle] 	#Convert the dates to string format of LE files
		dc_dim = [24,1] #24 hours, and 1 hist
		dc_dim.extend(LETKF.Ens[1].shape[1:])
		aux_dc = np.ones(dc_dim)
		LETKF.ensemble_mean_std() #Calculate the Ensemble mean and std
		#==Following A. Tsikerdekis et al.: Estimating aerosol emission from SPEXone under OSSEs, 2022===
		if LETKF.Ens_mean<=1.1:
			Ens_mean = LETKF.Ens_mean
			Ens_std = 1
		else:
			Ens_mean = LETKF.Ens_mean
			Ens_std = Ens_mean*0.9
		for ens in range(1,LE.Nens+1):
			i_day_LETKF_bef = -1
			for idate,date in enumerate(dates_cycle_new_cycle):
				i_day_LETKF = idate // Ts
				ens_str = f'{ens:02d}'
				if i_day_LETKF_bef!=i_day_LETKF: #Actualize just for a new Ts day
					aux_dc[:,0,:,:,:]= LETKF.create_rand_ensemble_dc(new_mean=Ens_mean,new_std=Ens_std,r=radius+1,spatial_corr=spatial_corr) #The same value for all the hours of the day
				print('No observations in this cycle. Writting the background DC file for the ensemble: '+ens_str+' at day '+date)
				filename = 'LE_'+LE.ID+'_smoother_xi_'+ens_str+'_dc_'+date+'_xa.nc'
				LE.write_one_dc_file(filename=filename,date=date,data=aux_dc)
				filename = 'LE_'+LE.ID+'_smoother_xi_'+ens_str+'_dc_'+date+'_xa_cycle_'+str(cycle)+'_.nc'
				LE.write_one_dc_file(filename=filename,date=date,data=aux_dc )
				i_day_LETKF_bef = i_day_LETKF
		continue
		
	#=Get the dates between the assimilation start and end of period =
	dates_cycle_new_dc = pd.date_range(new_start_date+timedelta(days=Ts),end_date_datetime,freq='1d',inclusive='left')	#Here it is neccesary to take each day, to create the new DC files
	dates_cycle_new_dc = [f.strftime('%Y%m%d') for f in dates_cycle_new_dc] 	#Convert the dates to string format of LE files
	
	dc_dim = [24,1] #24 hours, and 1 hist
	dc_dim.extend(LETKF.Ens_analy[1].shape[1:])
	aux_dc = np.ones(dc_dim)
	dc_xa_dim = [LE.Nens]
	dc_xa_dim.extend(dc_dim)
	aux_dc_xa = np.ones(dc_xa_dim)
	for idate,date in enumerate(dates_cycle_new_dc):
		i_day_LETKF = idate // Ts
		for ens in range(1,LE.Nens+1):
			ens_str = f'{ens:02d}'
			aux_dc[:,0,:,:,:] = LETKF.Ens_analy[ens][i_day_LETKF,:,:,:] #The same value for all the hours of the day
			print('Writting the Analysis DC file for the ensemble: '+ens_str+' at day '+date)
			filename = 'LE_'+LE.ID+'_smoother_xi_'+ens_str+'_dc_'+date+'_xa.nc'
			LE.write_one_dc_file(filename=filename,date=date,data=aux_dc)
			filename = 'LE_'+LE.ID+'_smoother_xi_'+ens_str+'_dc_'+date+'_xa_cycle_'+str(cycle)+'_.nc'
			LE.write_one_dc_file(filename=filename,date=date,data=aux_dc )
			aux_dc_xa[ens-1,:] = aux_dc
			
		print('Writting the Mean Analysis DC file at day '+date)
		filename = 'LE_'+LE.ID+'_smoother_xa_dc_'+date+'_xa.nc'
		LE.write_one_dc_file(filename=filename,date=date,data=aux_dc_xa.mean(axis=0))		
	#=Get the dates for the observation period of next cycle =
	dates_cycle_new_obs = pd.date_range(end_date_datetime+timedelta(days=1),end_date_datetime+timedelta(days=Ts+1),freq='1d',inclusive='left')	#Here it is neccesary to take each day, to create the new DC files
	dates_cycle_new_obs = [f.strftime('%Y%m%d') for f in dates_cycle_new_obs] 	#Convert the dates to string format of LE files

	aux_dc = np.ones(dc_dim)
	LETKF.ensemble_mean_std() #Calculate the Ensemble mean and std
	#==Following A. Tsikerdekis et al.: Estimating aerosol emission from SPEXone under OSSEs, 2022===
	if LETKF.Ens_mean<=1.1:
		Ens_mean = LETKF.Ens_mean
		Ens_std = 1
	else:
		Ens_mean = LETKF.Ens_mean
		Ens_std = Ens_mean*0.9
	
	for ens in range(1,LE.Nens+1):
		i_day_LETKF_bef = -1
		for idate,date in enumerate(dates_cycle_new_obs):
			i_day_LETKF = idate // Ts
			ens_str = f'{ens:02d}'
			if i_day_LETKF_bef!=i_day_LETKF: #Actualize just for a new Ts day
				aux_dc[:,0,:,:,:]= LETKF.create_rand_ensemble_dc(new_mean=Ens_mean,new_std=Ens_std,r=radius+1,spatial_corr=spatial_corr) #The same value for all the hours of the day
				print('Writting a new background DC file for the ensemble: '+ens_str+' at day '+date)
			else:
				print('Writting the previous background DC file for the ensemble: '+ens_str+' at day '+date)
			filename = 'LE_'+LE.ID+'_smoother_xi_'+ens_str+'_dc_'+date+'_xa.nc'
			LE.write_one_dc_file(filename=filename,date=date,data=aux_dc)
			filename = 'LE_'+LE.ID+'_smoother_xi_'+ens_str+'_dc_'+date+'_xb_cycle_'+str(cycle)+'_.nc'
			LE.write_one_dc_file(filename=filename,date=date,data=aux_dc )
			i_day_LETKF_bef = i_day_LETKF
	








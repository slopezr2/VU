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
LE.run_path = '/scratch-shared/slr/projects/Smoother_V2/'

#==Definition of the assimilation window==
Tb = 8		 #Background step. The total days of simulation.
Ts = 2 		 #Assimilation step. The last days of the background step where observations will be collected. This value guide the number of days using a constant dc factor.
if Tb%Ts!=0:
	status = 1
	error_message = 'The background step Tb is not divisible by the assimilation step Ts. Use a proper definition to avoid smoother issues.'
IF_NOT_OK(status,error_message)
Ta = Tb-Ts		#Analysis step. The initial days where the dc factor will be estimated.

#==Define id simulation name and dates of the==
LE.ID = 'Test_New_V2'
start_date = '2008-05-01 00:00:00'
Ncycles = 3			#Number of assimilation cycles.

#==Number of ensembles and localization radius==
LE.Nens = 3
radius = 1 #Degree
#==CSO instrument outputs to be assimilated==
CSO_instruments = ['polder-aod', 'polder-ae', 'polder-ssa-563']

#==Using spatial correlation for perturbation factor===
spatial_corr = True

#==Using initial simulation to create DC factors template files==
initial_cycle = True

#==If you need to restart the simulation from another cycle
first_cycle = 0
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
for cycle in range(first_cycle,Ncycles):
	#==Init a new Smoother instance==
	LETKF = LETKF_4D(Nens=LE.Nens,Tb=Tb,Ts=Ts,Ta=Ta)
	flaq_obs = 0		#Detect whether there are observations in the cycle or not
	
	#==Run the new cycle==
	#Check if it is the first cycle, to generate random initial dc
	if cycle==0:
		new_start_date =  start_date_datetime+timedelta(days=Ts*cycle)
		if initial_cycle==True:
			#==Run the model one day to generate the template dc files==
			LE.write_read_dc_file(read_dc_file='F',path='${my.project.dir}/dc_smoother',model='LE_${run.id}_smoother')		#Option to LEKF reads the dc from files
			LE.write_restart(restart='F')
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
	
	








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
LE=LE_model()
#==Define path to run the LE simulation===
LE.run_path='/scratch-shared/slr/projects/Smoother/'

#==Definition of the assimilation window==
Tb = 6		 #Background step. The total days of simulation.
Ts = 2 		 #Assimilation step. The last days of the background step where observations will be collected. This value guide the number of days using a constant dc factor.
if Tb%Ts!=0:
	status = 1
	error_message = 'The background step Tb is not divisible by the assimilation step Ts. Use a proper definition to avoid smoother issues.'
IF_NOT_OK(status,error_message)
Ta = Tb-Ts		#Analysis step. The initial days where the dc factor will be estimated.

#==Define id simulation name and dates of the==
LE.ID = 'Test_V6'
start_date = '2008-05-01 00:00:00'
Ncycles = 2			#Number of assimilation cycles.

#==Number of ensembles==
LE.Nens = 4
#============================== End of user definitions parameters ======================================================

start_date_datetime = datetime.strptime(start_date, '%Y-%m-%d %H:%M:%S')


#===Init assimilation cycles===
for cycle in range(Ncycles):
	#==Create a new cycle of simulations
	new_start_date =  start_date_datetime+timedelta(days=Ts*cycle)
	LE.start_date = new_start_date.strftime('%Y-%m-%d %H:%M:%S')
	end_date_datetime = new_start_date+timedelta(days=Tb)
	LE.end_date=end_date_datetime.strftime('%Y-%m-%d %H:%M:%S')
	LE.write_date()
	#==Run the new cycle==
	#Check if it is the first cycle, to generate random initial dc
# 	if cycle==0:
# 		LE.write_read_dc_file(read_dc_file='F',path='${my.project.dir}/dc_smoother',model='LE_${run.id}_smoother')		#Option to LEKF reads the dc from files
# 	else:
# 		LE.write_read_dc_file(read_dc_file='T',path='${my.project.dir}/dc_smoother',model='LE_${run.id}_smoother')		#Option to LEKF reads the dc from files
# 	LE.compile()
# 	#=Chance the .out and .err for each cycle
# 	LE.run_command = 'srun  --kill-on-bad-exit=1 --output=le_cycle_'+str(cycle)+'.run.out.%t  --error=le_cycle_'+str(cycle)+'.run.err.%t  ./lekf.x lekf.rc'
# 	LE.run()
	#==Init a new Smoother instance==
	LETKF = LETKF_4D(Nens=LE.Nens,Tb=Tb,Ts=Ts,Ta=Ta)
	for ens in range(1,LE.Nens+1):
		#==Read the dc, dc_lat, and dc_lon of each ensemble
		#=Get the dates between the assimilation start and end of period=
		dates_cycle = pd.date_range(new_start_date+timedelta(days=Ts),end_date_datetime,freq=str(Ts)+'d',inclusive='left')
		dates_cycle = [f.strftime('%Y%m%d') for f in dates_cycle] 	#Convert the dates to string format of LE files
		ens_str = f'{ens:02d}'
		for idate,date in enumerate(dates_cycle):
			words = [date]
			print(date)
			words.extend(['dc','xi'+ens_str+'a'])
			status = 0
			err_message = ''
			file_name = LE.list_output_files(words)
			print(type(file_name))
			if file_name==[] or len(file_name)>1:
				status = 1
				err_message = 'Problems geting the list of dc factors for cycle: '+str(cycle)+', ensemble: '+ens_str+', at day: '+date 
			else:
				file_name=file_name[0]
			IF_NOT_OK(status,error_message)
			#==Read dc values
			if idate==0: #First day of the cycle
				LETKF.Ens[ens]=LE.read_one_output(file_name=file_name,variable='dc',vector=False)
			else:
				LETKF.Ens[ens]=np.stack((LETKF.Ens[ens],LE.read_one_output(file_name=file_name,variable='dc',vector=False)))
		print(LETKF.Ens[ens].shape)


# for j in range(3):
# 	new_date_datetime = start_date_datetime+timedelta(days=j)
# 	new_date = new_date_datetime.strftime('%Y%m%d')
# 	for i in range(1,5):
# 		v =np.ones((24, 1, 6, 140, 100))*(i+j)
# 		ens = f'{i:02d}'
# 		LE. write_one_dc_file(filename='LE_'+LE.ID+'_smoother_'+'xi_'+ens+'_dc_'+new_date+'_xa.nc',date='20/05/2008',data=v)
# 	
# LE.compile()
# LE.run_command='srun  --kill-on-bad-exit=1 --output=le_reading.run.out.%t  --error=le_reading.run.err.%t  ./lekf.x lekf.rc'
# LE.run()

# #==Read de dc factors===
# new_date_file = (start_date_datetime+timedelta(days=1)).strftime('%Y%m%d')
# file_name = LE.list_output_files([new_date_file,'dc','xi03a'])
# dc_aux = LE.read_one_output(file_name=file_name[0],variable='dc',vector=False)
# 
# file_name = LE.list_output_files([new_date_file,'polder','data'])
# dc_aux = LE.read_one_output(file_name=file_name[0],variable='yr',vector=False)
# 










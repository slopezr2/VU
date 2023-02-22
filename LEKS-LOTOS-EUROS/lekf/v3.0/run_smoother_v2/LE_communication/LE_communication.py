import os
import sys
import subprocess
from datetime import datetime
from datetime import timedelta
from os import listdir
from os.path import isfile, join
from netCDF4 import Dataset
import numpy as np
import fileinput
import pathlib

def IF_NOT_OK(status,message):
	if status!=0:
		sys.exit('[ERROR SMOOTHER]             '+message)

def look_files(words,files):
	"""
	Select the files that contains the word from a list of files
	"""
	if type(words)==str:
		words_aux = words.split()
	elif type(words)==list:
		words_aux = words 
    
	files_with_word = files.copy()    
	for word in words_aux:
		filter = [word in file for file in files_with_word]
		files_with_word  =  [i for (i, v) in zip(files_with_word, filter) if v]

	return files_with_word

class LE_model:
	"""
	The class :py:mod:`LE_model` allows the interaction with a LEKF simulation.
	With this class is possible to compile, run, modify ID and simulations dates, read model outputs and write restart values of a LEKF simulation.
	"""
	def __init__(self):
		"""
		Create LE_model instance.
		"""
		self.ID = ''
		self.start_date = ''
		self.end_date = ''
		self.compile_command = './base/005/bin/setup-lekf proj/smoother_V2/002/rc/lekf_v2.rc'
		self.run_command = 'srun  --kill-on-bad-exit = 1 --output = le.run.out.%t  --error = le.run.err.%t  ./lekf.x lekf_v2.rc'
		self.run_path = '/scratch-shared/slr/projects/Smoother_V2/'
		self.lekf_rc_path = '../proj/smoother_V2/002/rc/'
		self.smoother_path = os.getcwd()
		self.compiled = False
		self.Nens = 20
	
	def write_ID(self):
		"""
		Write the ID value in the lekf_id.rc file, to create a new LEKF simulation.
		
		"""		
		write_id = 'run.id             :  '+self.ID
		
		file = self.lekf_rc_path+'lekf_v2.rc'
		status = 1
		try:
			for line in fileinput.input(file, inplace  =  1): 
				if line.startswith("run.id"):
					print(line.replace(line,write_id.rstrip()))
					status = 0
				else:
					print(line.replace(line,line.rstrip()))
			fileinput.close()
			if status==1:
				error_message = 'The file does not contain the run.id definition'
			else:
				error_message = ''
		except Exception as e:
			status = 1
			error_message = str(e)
		
		IF_NOT_OK(status,error_message)

	def write_date(self):
		"""
		Write the dates in the lotos-euros-dates.rc file, to define the LEKF simulation period.
		
		"""
		file = self.lekf_rc_path+'lotos-euros-v2.2.002.rc'		
		write_start_date = 'timerange.start            :  '+self.start_date
		write_end_date = 'timerange.end               :  '+self.end_date
		status = 1
		flaq_start = 0
		flaq_end = 0
		try:
			for line in fileinput.input(file, inplace  =  1): 
				if line.startswith('timerange.start'):
					print(line.replace(line,write_start_date.rstrip()))
					flaq_start = 1
				elif line.startswith('timerange.end'):
					print(line.replace(line,write_end_date.rstrip()))
					flaq_end = 1
				else:
					print(line.replace(line,line.rstrip()))
			fileinput.close()
			if flaq_start==1 and flaq_end==1:
				status = 0
			if status==1:
					error_message = 'The file does not contain the timerange.start or timerange.end definitions'
			else:
					error_message = ''
		except Exception as e:
			status = 1
			error_message = str(e)
		
		IF_NOT_OK(status,error_message)
	
	def write_N_ensembles(self):
		"""
		Write the number of ensembles in the lekf_v2.rc file.
		
		"""
		file = self.lekf_rc_path+'lekf_v2.rc'
		write_Nens = 'kf.nmodes             :   '+str(self.Nens)
		status = 1
		flaq_Nens = 0
		try:
			for line in fileinput.input(file, inplace  =  1): 
				if line.startswith('kf.nmodes'):
					print(line.replace(line,write_Nens.rstrip()))
					flaq_Nens = 1
				else:
					print(line.replace(line,line.rstrip()))
			fileinput.close()
			if flaq_Nens==1 :
				status = 0
			if status==1:
					error_message = 'The file does not contain the kf.nmodes  definition'
			else:
					error_message = ''
		except Exception as e:
			status = 1
			error_message = str(e)
		
		IF_NOT_OK(status,error_message)
	
	
	def write_restart(self,restart='F'):
		"""
		Write the LE option restart
		"""
		file = self.lekf_rc_path+'lekf_v2.rc'		
		write_restart = 'kf.restart  :  ' + restart
		status = 1
		flaq_restart = 0
		try:
			for line in fileinput.input(file, inplace  =  1): 
				if line.startswith('kf.restart  :'):
					print(line.replace(line,write_restart.rstrip()))
					flaq_restart = 1
				else:
					print(line.replace(line,line.rstrip()))
			fileinput.close()
			if flaq_restart==1 :
				status = 0
			if status==1:
					error_message = 'The file does not contain the kf.restart definition'
			else:
					error_message = ''
		except Exception as e:
			status = 1
			error_message = str(e)
		
		IF_NOT_OK(status,error_message)
	
	def write_read_dc_file(self,read_dc_file='T',path='',model=''):
		"""
		Write the read DC file option and the needed information. This option allows LE to read the DC values of 
		each ensemble from a file
		"""
		file = self.lekf_rc_path+'lekf_v2.rc'	
		write_read_dc_file = 'kf.read_dc_file  :  ' + read_dc_file
		write_read_dc_file_path = 'kf.read_dc_file.path  :  '+path
		write_read_dc_file_model= 'kf.read_dc_file.model  :  '+model
		status = 1
		flaq_read_dc_file = 0
		flaq_read_dc_file_path = 0
		flaq_read_dc_file_model = 0
		try:
			for line in fileinput.input(file, inplace  =  1): 
				if line.startswith('kf.read_dc_file.path'):
					print(line.replace(line,write_read_dc_file_path.rstrip()))
					flaq_read_dc_file_path = 1 
				elif line.startswith('kf.read_dc_file.model'):
					print(line.replace(line,write_read_dc_file_model.rstrip()))
					flaq_read_dc_file_model = 1
				elif line.startswith('kf.read_dc_file'):
					print(line.replace(line,write_read_dc_file.rstrip()))
					flaq_read_dc_file = 1
				else:
					print(line.replace(line,line.rstrip()))
			fileinput.close()
			if flaq_read_dc_file==1 and flaq_read_dc_file_path==1 and flaq_read_dc_file_model==1:
				status = 0
			if status==1:
					error_message = 'The file does not contain the kf.read_dc_file definitions of path and model'
			else:
					error_message = ''
		except Exception as e:
			status = 1
			error_message = str(e)
		IF_NOT_OK(status,error_message)
	
	
	
	def compile(self,option = ''):
		"""
		Compile the LEKF simulation.
		
		Parameters
		
		* ``option`` : Compilation options. -n for a new compilation.
		
		"""		
		compile_command = self.compile_command+' '+option
		os.chdir('..')
		status = subprocess.run(compile_command, shell = True)
		IF_NOT_OK(status.returncode,'Error compiling the LEKF simulation')
		self.compiled = True
		os.chdir(self.smoother_path)
	
	def run(self):
		"""
		Run the LEKF simulation.
		
		"""		
		LE_run_path = self.run_path+self.ID+'/run/'
		
		if self.compiled:
			os.chdir(LE_run_path)
			status = subprocess.run(self.run_command, shell = True)
			IF_NOT_OK(status.returncode,'Error running the LEKF simulation')
			os.chdir(self.smoother_path)
		else:
			self.compile()
			os.chdir(LE_run_path)
			status = subprocess.run(self.run_command, shell = True)
			IF_NOT_OK(status.returncode,'Error running the LEKF simulation')
			os.chdir(self.smoother_path)			
	
	def list_output_files(self,words = None):
		"""
		List the output files of the LEKF simulation that contains specific words. if words is empty, none or '', all the files will be listed.
		
		"""		
		LE_output_path = self.run_path+self.ID+'/output/'
		try:
			output_files  =  [f for f in listdir(LE_output_path) if isfile(join(LE_output_path, f))]
			output_files.sort()
			status = 0
			error_message = ''
		except Exception as e:
			status = 1
			error_message = str(e)
		IF_NOT_OK(status,error_message)
		
		if (words==None or words=='' or not(words)):
			output_files = output_files
		else:
			output_files = 	look_files(words,output_files)
		return output_files
		
	def read_one_output(self,file_name,variable,vector=True,squeeze=True):
		"""
		Read one single LEKF output file and create a vector with the values or a variable with the output dimensions.
		
		"""		
		LE_output_path = self.run_path+self.ID+'/output/'
		file = LE_output_path+file_name
		try: #Try to open the required file
			nc_file = Dataset(file)
			status = 0
			error_message = ''
		except Exception as e:
			status = 1
			error_message = str(e)
		IF_NOT_OK(status,error_message)
		
		if variable in nc_file.variables: #check if the required variable exist in the nc file
			if squeeze==True:
				nc = nc_file.variables[variable][:].squeeze() #Remove dimensions on length one
			else:
				nc = nc_file.variables[variable][:] #Remove dimensions on length one
			status = 0
			error_message = ''
		else:
			status = 1
			error_message = 'The required variable '+variable+'  is not in the nc file.'+file
		IF_NOT_OK(status,error_message)
		
		if vector==True:
			if nc.ndim==1: #If ndim is 1, the variable is already a vector likely the observations vector
				nc = nc
			elif nc.ndim==3 or nc.ndim==2: #The dimensions are lat,lon (just one noise), or noise,lat,lon
				nc = nc.flatten()
			else: #If ndim is greater than 3, we have to take a time (first dimension) value and flatten it. Likely is the state vector (dc files)
				#If there are more than one noise, the flatten vector will be organice as noise1x1y1,noise1x1y2,..noise1x2y1,noise1x2y2,...,noise2x1y1,noise2x1y2,......
				nc = nc[0,:].flatten()
		else:
			if nc.ndim>3:
				nc = nc[0,:] #Taking just the first time of the file.
		return nc
	
	def read_multiple_outputs(self,files,variable,vector=True,squeeze=True):
		"""
		Read the same variable of multiple files, and concatenate the values in a single vector.
		"""		
		if type(files)==list:
			status = 0
			error_message = ''
		else:
			status = 1
			error_message = 'Files has to be a list of file names'
		IF_NOT_OK(status,error_message)	
		nc = self.read_one_output(files[0],variable)
		for file in files[1:]:
			nc_aux = self.read_one_output(file,variable,vector,squeeze)
			if vector==True:
				nc = np.concatenate([nc,nc_aux])
			else:
				nc=np.stack([nc,nc_aux])
		
		return nc

	def write_one_dc_file(self,filename,date,data):
		"""
		Write a dc-like file with the information of the variable data. The function takes a existing dc file from LE to copy all the dimensions,
		variables and attributes, and just change de DC data.
		"""
		
		LE_output_path = self.run_path+self.ID+'/output/'
		LE_dc_path = self.run_path+self.ID+'/dc_smoother/'
		
		pathlib.Path(LE_dc_path).mkdir(parents=True, exist_ok=True)
		status = 0
		error_message = ''
		files = self.list_output_files(words = ['_dc_','xa','.nc']) #List the dc files
		if files==[]:
			status = 1
			error_message = 'There are no DC files available to use as template'
		IF_NOT_OK(status,error_message)
		
		#Create de source file and target file
		src  = Dataset(LE_output_path+files[1])#Read the second dc file as dummy file 
		trg  =  Dataset(LE_dc_path+filename, mode = 'w',format='NETCDF3_CLASSIC')
		#Verify that the source file is a DC file
		if not('dc' in src.variables.keys()):
			status = 1
			error_message = 'The dummy DC file does not containt the dc variable'
		IF_NOT_OK(status,error_message)
		
		# Create the dimensions of the file   
		for name, dim in src.dimensions.items():
			trg.createDimension(name, len(dim) if not dim.isunlimited() else None)
			# Copy the global attributes
			trg.setncatts({a:src.getncattr(a) for a in src.ncattrs()})
		
		# Create the variables in the file
		aux_var_names = ['dc','time']
		time_2008 = datetime.strptime('31/12/2007 22:0:0' , '%d/%m/%Y %H:%M:%S').timestamp()
		date_timestamp = datetime.strptime(date, '%Y%m%d').timestamp()
		aux_time = [f*3600-time_2008+date_timestamp for f in range(24)]
		
		for name, var in src.variables.items():
			trg.createVariable(name, var.dtype, var.dimensions)
			# Copy the variable attributes
			trg.variables[name].setncatts({a:var.getncattr(a) for a in var.ncattrs()})
			if name in aux_var_names:
				if name=='dc':
					try:
						dc_data = np.reshape(data,trg['dc'][:].shape)
						dc_data=dc_data.astype(dtype="float32") #Ensure that the data is float32
						trg.variables[name][:]  =  dc_data	
					except Exception as e:
						status = 1
						error_message = 'Problems with the DC dimensions: '+str(e)
				else:
					trg.variables[name][:]  =  aux_time
				IF_NOT_OK(status,error_message)
			else:
				trg.variables[name][:]  =  src.variables[name][:]	
		
		# Save the file
		trg.close()




















        
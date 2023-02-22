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


def getmeanofneighbors(matrix, i, j,r):
	region = matrix[max(0, i-r) : i+r+1,max(0, j-r) : j+r+1]
	return (np.sum(region) - matrix[i, j])/(region.size-1) # Sum the region and subtract center


class LETKF_4D:
	"""
	The class :py:mod:`LETKF_4D` contains routines to implement the LETKF analysis using observations and states for different time steps.
	Additionally, create the new ensemble members for next anaysis step.
	"""
	def __init__(self,Nens,Tb,Ts,Ta):
		"""
		Creates LETKF_4D instance.
		"""
		self.Ens = dict.fromkeys(range(1,Nens+1))
		self.Ens_analy = dict.fromkeys(range(1,Nens+1))
		self.Y_Ens = dict.fromkeys(range(1,Nens+1))
		self.Nens = Nens

	def get_subdomain(self,lat_point,lon_point,radius):
		"""
		Obtains de values corresponding with the subdomain centered in (lat_poin,lon_point) with radius=radius. Create new variables with the state vector and observation vector.
		"""
		self.lat_point = lat_point
		self.lon_point = lon_point
		subdomain_lat_dc = ((self.Ens_lat[lat_point]-radius)<self.Ens_lat) & (self.Ens_lat<=(self.Ens_lat[lat_point]+radius))
		subdomain_lon_dc = ((self.Ens_lon[lon_point]-radius)<self.Ens_lon) & (self.Ens_lon<=(self.Ens_lon[lon_point]+radius))
		
		subdomain_yr = (((self.Ens_lat[lat_point]-radius)<self.Y_lat) & (self.Y_lat<=(self.Ens_lat[lat_point]+radius))) & (((self.Ens_lon[lon_point]-radius)<self.Y_lon) & (self.Y_lon<=(self.Ens_lon[lon_point]+radius)))
		self.sub_Y = self.Y[subdomain_yr]
		self.sub_Y_error = self.Y_error[subdomain_yr]
		
		for ens in range(1,self.Nens+1):
			subdomain_dc = self.Ens[ens][:,:,subdomain_lat_dc.nonzero()[0][:,np.newaxis],subdomain_lon_dc.nonzero()[0]]
			#==save the position of the global point in the subdomain==
			if ens==1:
				self.index_sub = np.where(subdomain_dc[0,0,:,:]==self.Ens[ens][0,0,lat_point,lon_point])
				self.sub_Ens_shape = subdomain_dc.shape
			subdomain_dc = subdomain_dc.flatten()
			subdomain_ys = self.Y_Ens[ens][subdomain_yr]
			if ens==1:
				self.sub_Ens = subdomain_dc[:,np.newaxis]
				self.sub_Y_Ens = subdomain_ys[:,np.newaxis]
			else:
				self.sub_Ens = np.concatenate([self.sub_Ens,subdomain_dc[:,np.newaxis]],axis=1)	#Create the ensemble matrix with dimensions n*N
				self.sub_Y_Ens = np.concatenate([self.sub_Y_Ens,subdomain_ys[:,np.newaxis]],axis=1)	#Create the ensemble outputs matrix with dimensions m*N


	def analysis_step(self):
		"""
		Applys the LETKF analysis to each subdomain.
		"""	
		#==Calculate the state ensemble mean and observation ensemble mean===
		x_mean = self.sub_Ens.mean(axis=1)
		y_mean = self.sub_Y_Ens.mean(axis=1)
		
		#==Calculate deviation matrices===
		#X = [ ..., xi - <xi>, ..]
		#Y = [ ..., Hxi - <Hxi>, ..]
		X = np.add(self.sub_Ens, -x_mean[:,None])
		Y = np.add(self.sub_Y_Ens, -y_mean[:,None])
		
		#==Calculate matrix C===
		#C = Y^T R^{-1}
		C = np.zeros(Y.T.shape)
		for i in range(Y.shape[0]):
			C[:,i] = [y / self.sub_Y_error[i]**2 for y in (Y[i,:])]
		
		#==Calculate Pa^-{1}==
		#Pa^{-1} = [ C Y + (m-1) I ]
		inv_Pa = np.dot(C,Y) + np.diag(np.ones(C.shape[0])*(self.Nens-1))
		
		#==Calculate eigendecomposition of Pa^{-1}
		#Pa^{-1} = Q Lambda Q^T
		try: 
			Lambda,Q = np.linalg.eigh(inv_Pa)
			status = 0
			error_message = ''
		except Exception as e:
			status = 1
			error_message = str(e)
		IF_NOT_OK(status,error_message)
		
		#==Calculate Pa ===
		#Pa = Q Lambda^{-1} Q^T
		Pa = np.zeros(inv_Pa.shape)
		for i in range(Pa.shape[0]):
			Pa[:,i] = Q[:,i] / Lambda[i] 
		Pa = np.dot(Pa,Q.T)
		
		#==Calculate the symmetric factorition==
		# Symetric square roots W of (m-1)Pa :
		#  (m-1) Pa = (m-1) Q Lambda^{-1} Q^T 
		#   W W^T = (m-1) Q Lambda^{-1/2} Q^T Q Lambda^{-1/2} Q^T
		#    W = sqrt(m-1) Q Lambda^{-1/2} Q^T
		Wa = np.zeros(Pa.shape)
		for i in range(Wa.shape[0]):
			Wa[:,i] = np.sqrt(self.Nens - 1) * Q[:,i] / np.sqrt(Lambda[i])
		Wa = np.dot(Wa,Q.T)
		
		#==Calculate the analysis weight in the ensemble space==
		# w = Pa C ( y - ybar )
		wa = np.dot(Pa,np.dot(C,self.sub_Y - y_mean))
		
		#==Compute the analysis mean==
		xa_mean = x_mean + np.dot(X,wa)
		xa_mean = xa_mean*1.5
		#==Create the analysis ensemble==
		self.sub_Ens_analy = np.zeros(self.sub_Ens.shape)
		for i in range(self.Nens):
			self.sub_Ens_analy[:,i] = xa_mean + np.dot(X,Wa[:,i])
	
		#Test
		self.sub_Ens_analy = np.zeros(self.sub_Ens.shape)
		for i in range(self.Nens):
			self.sub_Ens_analy[:,i] = x_mean
	
	def map_local_global(self,ens):
		"""
		Maps the local domain analysis into the global analysis.
		"""
		#==Reshape each ensemble analysis state==
		#new shape (\ t,noise,lat,lon\)

		subdomain_dc = self.sub_Ens_analy[:,ens-1].reshape(self.sub_Ens_shape)
		
		subdomain_dc[:,:,self.index_sub[0][0],self.index_sub[1][0]] = 59
		return subdomain_dc[:,:,self.index_sub[0][0],self.index_sub[1][0]]
	
	def ensemble_mean_std(self):
		"""
		Calculate the ensemble mean and standard deviation of all the domain.
		"""
		
		for ens in range(1,self.Nens+1):
			dc = self.Ens[ens][:,:,:]
			#==save the position of the global point in the subdomain==
			dc = dc.flatten()
			if ens==1:
				Ens_dc = dc[:,np.newaxis]
			else:
				Ens_dc = np.concatenate([Ens_dc,dc[:,np.newaxis]],axis=1)	#Create the ensemble matrix with dimensions n*N
			
		self.Ens_mean = 	Ens_dc.mean()
		self.Ens_std = Ens_dc.std()
	
	def create_rand_domain(self,repeat=5,new_mean=1,new_std=1,r=2,nx=140,ny=100,spatial_corr=True):
		"""
		Create a random DC with spatial correlation, following A. Tsikerdekis et al.: Assimilating aerosol optical properties related to size and absorption,2021
		"""
		
		if spatial_corr==True:
			domain = 1+np.random.randn(nx,ny)*10
			for k in range(repeat):
				domain_backup = domain.copy()
				for i in range(domain.shape[0]):
					for j in range(domain.shape[1]):
						domain[i,j] = getmeanofneighbors(domain_backup, i, j,r)
			domain = np.exp(domain)
			domain = new_std*((domain-domain.mean())/domain.std())+new_mean
			domain[domain<0] = 0
		else:
			domain = np.ones((nx,ny))*(new_mean+np.random.randn(1,1)*new_std) #All the correction factors with the same value
		return domain
	
	def create_rand_ensemble_dc(self,repeat=5,r=2,new_mean=1,new_std=1,nnoise=6,nlat=110,nlon=140,auto_shape=True,spatial_corr=True):
		"""
		Create an ensemble of random DC spatial factor for each noise.
		"""
		if auto_shape==True:
			dim_dc = self.Ens[1].shape[1:] #Shape (/noise,lat,lon / )
		else:
			dim_dc = [nnoise,nlat,nlon]
		dc_aux = np.zeros(dim_dc)
		for i in range(dim_dc[0]):
			dc_aux[i,:,:] = self.create_rand_domain(repeat=repeat,new_mean=new_mean,new_std=new_std,r=r,nx=dim_dc[1],ny=dim_dc[2],spatial_corr=spatial_corr)
		
		return dc_aux
		
		
		
		
		
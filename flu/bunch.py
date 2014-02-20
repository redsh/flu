import numpy as np
from stencils.cyl import *
from fields import Sim
from matplotlib.pyplot import *

_z,_r,   _uz,_ur,   _q,_m = [0,1,2,3,4,5]

class kinetic_bunch(object):

	def __init__(self, N, mean, sigma, gamma, emit):

		data = np.zeros( (6, N) )

		data[_z, :] = np.random.normal(loc=mean[0],scale=sigma[0],size=(N))
		data[_r, :] = np.random.normal(loc=mean[0],scale=sigma[0],size=(N))
		data[_r, :] = np.abs(data[_r, :])

		data[_uz, :] = gamma
		#TODO mettere emittanza bene
		data[_ur, :] = np.random.normal(loc=0, scale=gamma*emit,size=(N))

		data[_q, :] = 1.*(data[_r, :]+dz*0.5)
		data[_m, :] = 1.

		self.data = data
		self.shape_order = 1

	def coordinates(self):
		zz,rr = self.data[_z,:],self.data[_r,:]
		
		idz = 1.0/dz
		idr = 1.0/dr

		nz = (zz-zmin) * idz
		nr = (rr) * idr 

		if self.shape_order == 0:
			i = np.floor( nz + 0.5 ).astype(int)
			j = np.floor( nr + 0.5 ).astype(int)

		elif self.shape_order == 1:
			nz = (zz-zmin) * idz
			nr = (rr) * idr 

			i = np.floor( nz ).astype(int)
			j = np.floor( nr ).astype(int)

			wz = nz - i
			wr = nr - j

		return i,j,nz,nr

	def interpolate(self, fields):
		ret = []

		i,j,nz,nr = self.coordinates()

		if self.shape_order == 0:
			for f in fields:
				ret += [ f[j,i] ]

		elif self.shape_order == 1:
			wz = nz - i
			wr = nr - j

			for f in fields:
				ret += [  f[j+1,i+1]*(wr)*(wz) \
						+ f[j  ,i+1]*(wr)*(1.0-wz) \
						+ f[j+1,i  ]*(1.0-wr)*(wz)  \
						+ f[j  ,i  ]*(1.0-wr)*(1.0-wz) 
						]

		self.debug = dict(i=i,j=j,ez=ret[0])

		return ret

	def deposition(self, s, values):
		ret = []

		#for histograms
		lz = np.linspace(0,s.Ez.shape[1],s.Ez.shape[1]+1)
		lr = np.linspace(0,s.Ez.shape[0],s.Ez.shape[0]+1)

		print lz.min(),lz.max(),lz.shape

		zz,rr = self.data[_z,:],self.data[_r,:]
		i,j,nz,nr = self.coordinates()

		if self.shape_order == 0:
			for v in values:
				J = np.zeros(s.Ez.shape)
				J,xx,yy = np.histogram2d(j,i,bins=(lr,lz),weights=v)
				ret += [J]
			error('not implemented')

		elif self.shape_order == 1:
			wz = nz - i
			wr = nr - j

			wz = 1.0-wz
			wr = 1.0-wr

			inv_j_0 = j
			inv_j_0[inv_j_0 == 0] = 0.25
			inv_j_0 = 1.0/inv_j_0

			inv_j_1 = 1.0/(j+1.0)

			for v in values:

				J ,xx,yy = np.histogram2d(j,i,bins=(lr,lz),weights=v*(wr)*(wz)*inv_j_0)
				
				J1,xx,yy = np.histogram2d(j,i,bins=(lr,lz),weights=v*(wr)*(1.0-wz)*inv_j_0)
				J += J1

				J1,xx,yy = np.histogram2d(j,i,bins=(lr,lz),weights=v*(1.0-wr)*(wz)*inv_j_1)
				J += J1
				
				J1,xx,yy = np.histogram2d(j,i,bins=(lr,lz),weights=v*(1.0-wr)*(1.0-wz)*inv_j_1)
				J += J1
				
				ret += [J]

		self.debug['J'] = ret

		return ret

	def deriv(self, s):

		data = self.data
		N = data.shape[1]

		zz,rr,ur,uz,q = self.data[_z,:],self.data[_r,:],self.data[_ur,:],self.data[_uz,:],self.data[_q,:]

		Ez_sampled,Er_sampled,Bp_sampled = self.interpolate([s0.Ez,s0.Er,s0.Bp])

		gamma = np.sqrt(1.+uz*uz+ur*ur)
		beta_z = uz/gamma
		beta_r = ur/gamma

		Jz, Jr = self.deposition(s,[q*beta_z,q*beta_r])
		#indici su griglia

		ds_dt = np.zeros( (6, N) )
		
		ds_dt[_z, :] = beta_z - 1.
		ds_dt[_r, :] = beta_r

		ds_dt[_uz, :] = -Ez_sampled - beta_r*Bp_sampled
		ds_dt[_ur, :] = -Er_sampled + beta_z*Bp_sampled

		return ds_dt,Jz,Jr



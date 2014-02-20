import numpy as np
from stencils.cyl import *
from fields import Sim
from matplotlib.pyplot import *

_z,_r,   _uz,_ur,   _q,_m = [0,1,2,3,4,5]
shape_order = 1
debug = {}

def init_kinetic_bunch(N, mean, sigma, gamma, emit, q0):
	data = np.zeros((6,N))
	data[_z, :] = np.random.normal(loc=mean[1],scale=sigma[1],size=(N))
	
	data[_r, :] = np.random.normal(loc=mean[0],scale=sigma[0],size=(N))
	data[_r, :] = np.abs(data[_r, :])

	data[_uz, :] = gamma
	#TODO mettere emittanza bene
	data[_ur, :] = np.random.normal(loc=0, scale=gamma*emit,size=(N))

	data[_q, :] = (data[_r, :]/dr+0.5)*q0
	data[_m, :] = 1.

	return data

def coordinates(data):
	zz,rr = data[_z,:],data[_r,:]
	
	idz = 1.0/dz
	idr = 1.0/dr

	nz = (zz-zmin) * idz
	nr = (rr) * idr 

	if shape_order == 0:
		i = np.floor( nz + 0.5 ).astype(int)
		j = np.floor( nr + 0.5 ).astype(int)

	elif shape_order == 1:
		nz = (zz-zmin) * idz
		nr = (rr) * idr 

		i = np.floor( nz ).astype(int)
		j = np.floor( nr ).astype(int)

		wz = nz - i
		wr = nr - j

	return i,j,nz,nr

def interpolate(data, fields):
	ret = []

	i,j,nz,nr = coordinates(data)

	if shape_order == 0:
		for f in fields:
			ret += [ f[j,i] ]

	elif shape_order == 1:
		wz = nz - i
		wr = nr - j

		for f in fields:
			ret += [  f[j+1,i+1]*(wr)*(wz) \
					+ f[j  ,i+1]*(wr)*(1.0-wz) \
					+ f[j+1,i  ]*(1.0-wr)*(wz)  \
					+ f[j  ,i  ]*(1.0-wr)*(1.0-wz) 
					]

	global debug
	debug = dict(i=i,j=j,ez=ret[0])

	return ret

def deposition(data, s, values):
	ret = []

	#for histograms
	lz = np.linspace(0,s.Ez.shape[1],s.Ez.shape[1]+1)
	lr = np.linspace(0,s.Ez.shape[0],s.Ez.shape[0]+1)

	#print lz.min(),lz.max(),lz.shape

	zz,rr = data[_z,:],data[_r,:]
	i,j,nz,nr = coordinates(data)

	if shape_order == 0:
		for v in values:
			J = np.zeros(s.Ez.shape)
			J,xx,yy = np.histogram2d(j,i,bins=(lr,lz),weights=v)
			ret += [J]
		error('not implemented')

	elif shape_order == 1:
		wz = nz - i
		wr = nr - j

		wz = 1.0-wz
		wr = 1.0-wr

		inv_j_0 = j.astype(float)
		inv_j_0[inv_j_0 <= 0] = 1.0/6.
		#print 'inv_j_0',inv_j_0.min(),inv_j_0.max(),'m0',(inv_j_0 <= 0).sum()
		#print 'inv_j_0_2',inv_j_0.min(),inv_j_0.max(),'m0',(inv_j_0 <= 0).sum()
		
		inv_j_0 = 1.0/inv_j_0

		inv_j_1 = 1.0/(j+1.0)

		for v in values:

			J ,xx,yy = np.histogram2d(j+1,i+1,bins=(lr,lz),weights=v*(wr)*(wz)*inv_j_0)
			
			J1,xx,yy = np.histogram2d(j,i+1,bins=(lr,lz),weights=v*(wr)*(1.0-wz)*inv_j_0)
			J += J1

			J1,xx,yy = np.histogram2d(j+1,i,bins=(lr,lz),weights=v*(1.0-wr)*(wz)*inv_j_1)
			J += J1
			
			J1,xx,yy = np.histogram2d(j,i,bins=(lr,lz),weights=v*(1.0-wr)*(1.0-wz)*inv_j_1)
			J += J1
			
			ret += [J]

	global debug
	debug['J'] = ret

	return ret

def deriv(data, s):

	N = data.shape[1]

	zz,rr,ur,uz,q = data[_z,:],data[_r,:],data[_ur,:],data[_uz,:],data[_q,:]

	Ez_sampled,Er_sampled,Bp_sampled = interpolate(data,[s0.Ez,s0.Er,s0.Bp])

	gamma = np.sqrt(1.+uz*uz+ur*ur)
	beta_z = uz/gamma
	beta_r = ur/gamma

	Jz, Jr = deposition(data, s,[q*beta_z,q*beta_r])
	#indici su griglia

	ds_dt = np.zeros( (6, N) )
	
	ds_dt[_z, :] = beta_z - 1.
	ds_dt[_r, :] = beta_r

	ds_dt[_uz, :] = -Ez_sampled - beta_r*Bp_sampled
	ds_dt[_ur, :] = -Er_sampled + beta_z*Bp_sampled

	#ds_dt = np.zeros( (6, N) )
	ds_dt *= cfg.bunch_adiabatic_T*0

	return ds_dt,Jz,Jr



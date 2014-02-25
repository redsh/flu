import numpy as np
from stencils.cyl import *
from fields import Sim
from plot import *

_z,_r,   _uz,_ur,   _q,_m = [0,1,2,3,4,5]
shape_order = 0

def init_kinetic_bunch(s, N, mean, sigma, gamma, emit, q0):
	data = np.zeros((6,N))
	
	data[_z, :] = np.random.normal(loc=mean[1],scale=sigma[1],size=(N))

	y = np.random.normal(loc=mean[0],scale=sigma[0],size=(N))
	z = np.random.normal(loc=mean[0],scale=sigma[0],size=(N))
	data[_r, :] = np.sqrt(y*y+z*z)
	data[_r, :] = np.random.normal(loc=mean[0],scale=sigma[0],size=(N))

	data[_uz, :] = gamma*np.sqrt(1.0-1.0/(gamma*gamma))
	#TODO mettere emittanza bene
	if emit > 0:
		data[_ur, :] = np.random.normal(loc=0, scale=gamma*emit,size=(N))

	data[_q, :] = np.abs(data[_r, :])/dr+0.5
	data[_m, :] = 1.

	#deposition

	rhob   = deposition(data,s,[ data[_q,:] ])[0]
	center = (mean[0]/dr,(mean[1]-zmin)/dz)
	delta  = (sigma[0]/dr*0.5,sigma[1]/dz*0.5)

	rhob = rhob[center[0]:center[0]+delta[0], center[1]-delta[1]:center[1]+delta[1]]
	data[_q, :] *= q0/rhob.mean()

	imshow(rhob*(q0/rhob.mean()))
	colorbar()
	show()
	#fads
	#fasdfds

	#r0=data.rmin+j*data.dr;
	#r0=r0+(k+0.5)*dr;
	#particles.q[n]=density*r0*idr/Nppc;

	return data

def coordinates(data):
	zz,rr = data[_z,:],data[_r,:]
	
	idz = 1.0/dz
	idr = 1.0/dr

	nz = (zz-zmin) * idz
	nr = np.abs(rr * idr)

	if shape_order == 0:
		i = np.floor( nz + 0.5 ).astype(int)
		j = np.floor( nr + 0.5 ).astype(int)

	elif shape_order == 1:
		nz = (zz-zmin) * idz
		nr = np.abs(rr) * idr 

		i = np.floor( nz ).astype(int)
		j = np.floor( nr ).astype(int)

		wz = nz - i
		wr = nr - j

	return i,j,nz,nr

def interpolate(data, fields):
	ret = []

	i,j,nz,nr = coordinates(data)
	i = np.clip(i,0,fields[0].shape[1]-1)
	j = np.clip(j,0,fields[0].shape[0]-1)

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


	return ret

def deposition(data, s, values):
	ret = []

	#for histograms
	#lz = np.linspace(0,s.Ez.shape[1],s.Ez.shape[1]+1)
	#lr = np.linspace(0,s.Ez.shape[0],s.Ez.shape[0]+1)
	kw = dict(bins=s.Ez.shape,range=((0,s.Ez.shape[0]),(0,s.Ez.shape[1])))

	#print lz.min(),lz.max(),lz.shape

	i,j,nz,nr = coordinates(data)

	if shape_order == 0:
		for v in values:
			J = np.zeros(s.Ez.shape)

			inv_j_0 = j.astype(float)
			inv_j_0[inv_j_0 <= 0] = -.125
			inv_j_0 = 1.0/(inv_j_0+0.5)

			J ,xx,yy = np.histogram2d(j,i,weights=v*inv_j_0, **kw)
			ret += [J]

	elif shape_order == 1:
		wz = nz - i
		wr = nr - j

		wz = 1.0-wz
		wr = 1.0-wr

		inv_j_0 = j.astype(float)
		inv_j_0[inv_j_0 <= 0] = 1.0/6.
		inv_j_0 = 1.0/inv_j_0

		inv_j_1 = 1.0/(j+1.0)

		for v in values:

			J ,xx,yy = np.histogram2d(j+1,i+1,weights=v*(wr)*(wz)*inv_j_1, **kw)
			
			J1,xx,yy = np.histogram2d(j,i+1,weights=v*(wr)*(1.0-wz)*inv_j_0, **kw)
			J += J1

			J1,xx,yy = np.histogram2d(j+1,i,weights=v*(1.0-wr)*(wz)*inv_j_1, **kw)
			J += J1
			
			J1,xx,yy = np.histogram2d(j,i,weights=v*(1.0-wr)*(1.0-wz)*inv_j_0, **kw)
			J += J1
			


			ret += [J]

	
	return ret

def deriv(data, s):

	N = data.shape[1]

	rr,ur,uz,q = data[_r,:],data[_ur,:],data[_uz,:],data[_q,:]
	sgn = (rr > 0)*2.0-1.0

	Ez_sampled,Er_sampled,Bp_sampled = interpolate(data,[s0.Ez,s0.Er,s0.Bp])

	gamma = np.sqrt(1.+uz*uz+ur*ur)
	beta_z = uz/gamma
	beta_r = ur/gamma

	Jz, Jr = deposition(data, s,[-q*beta_z, -q*beta_r*sgn])
	#indici su griglia

	ds_dt = np.zeros( (6, N) )
	
	ds_dt[_z, :] = beta_z - 1.
	ds_dt[_r, :] = beta_r

	ds_dt[_uz, :] = -Ez_sampled - beta_r*Bp_sampled
	ds_dt[_ur, :] = -Er_sampled + beta_z*Bp_sampled

	#if cfg.bunch_adiabatic_T < 0.999:
	if 1:
		ds_dt *= 0.

	return ds_dt,Jz,Jr


## field initialization
import poisson

def rho(data,s):
	return deposition(data,s,[data[_q,:]])[0]

def initial_fields_rho(rho_,gamma,s,convergence=1e-4):
	solver = poisson.PoissonSolver2DCylGamma(gamma,dr,dz)
	phi = solver.cycle(None,rho_,callback=None,max_iters=10000,convergence=convergence)
	
	figure(figsize=(18,5))
	cols = 7
	plot_field(z,r,phi,'phi',x=1,columns=cols)
	

	for f in [phi]:
		f[:,-1] = 0.
		f[:,-2] = 0.

	beta = np.sqrt(1.0-1.0/(gamma*gamma))

	Ez = -diff_z_1(phi)*(1.0/(gamma*gamma))
	Er = -diff_r_1(phi)
	Bphi = Er*beta

	divE = diff_z_1(Ez)+div_r(Er,False)

	plot_field(z,r,Ez,'Ez',x=2,columns=cols)

	plot_field(z,r,Er,'Er',x=3,columns=cols)

	plot_field(z,r,Bphi,'Bp',x=4,columns=cols)	

	plot_field(z,r,divE+rho_,'divE-rho',x=5,columns=cols)

	plot_field(z,r,divE,'divE',x=6,columns=cols)

	plot_field(z,r,rho_,'rho',x=7,columns=cols)

	show()


	#fdasfds

	return Ez,Er,Bphi
	

def initial_fields(b,s):
	return initial_fields_rho(rho(b,s),s)






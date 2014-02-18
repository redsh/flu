import copy,math
import numpy as np

#if current_config.geometry == 'cyl':

class Sim(object):
	def __init__(self):
		pass

	def _unpack(self,glob):
		for p in dir(self):
			if p[0] != '_':
				glob[p] = getattr(self,p)
				#print p 

	#def pack(self, locs): #TODO can be static
	#	r = Sim()
	#	for p in locs.keys():
	#		if p[0] != '_':
	#			setattr(r,p,locs[p])
	#			print 'packing'+p
	#	return r

	def __mul__(self,b):
		ret = Sim()
			
		if isinstance(b,Sim):
			for k in self.__dict__:
				if hasattr(b,k):
					setattr(ret,k,  getattr(self,k)*getattr(b,k) )
		else:
			for k in self.__dict__:
				setattr(ret,k,  getattr(self,k)*b )

		return ret

	def __add__(self,b):
		ret = Sim()
			
		if isinstance(b,Sim):
			for k in self.__dict__:
				if hasattr(b,k):
					setattr(ret,k,  getattr(self,k)+getattr(b,k) )
		else:
			for k in self.__dict__:
				setattr(ret,k,  getattr(self,k)+b )

		return ret


def init_fields(current_config, globals):
	ax_z,ax_r = (1,0)

	cfg = current_config
	z,r = cfg.zr()
	dz,dr = cfg.dz,cfg.dr
	dt = cfg.dt

	one__r = r.copy()
	print 'fields dr',dr

	for i in range(one__r.shape[0]):
		if i != 0:
			one__r[i,:] = 1.0/(dr*i)

	dd  = (cfg.dr,cfg.dz)
	shape = z.shape

	#status
	#e.m.

	s0 = Sim()

	s0.Ez = np.zeros(shape)
	s0.Er = np.zeros(shape)
	s0.Bp = np.zeros(shape)

	#plasma: fluid
	#s0.rho = np.zeros(shape)
	s0.ruz = np.zeros(shape)
	s0.rur = np.zeros(shape)

	#envelope
	s0.Ar = np.zeros(shape)
	s0.Ai = np.zeros(shape)
	#s0.A2 = np.zeros(shape)

	#init: envelope
	zl0 = 0

	#cfg.z,cfg.r = z,r

	L = copy.copy(locals())
	for l in L:
		if l != 'current_config' and l != 'globals':
			globals[l] = locals()[l]

def gaussian_envelope(t):
	gamma_g=cfg.k0_kp;
	delta_g=math.sqrt(1-1/(gamma_g*gamma_g))-1;
	
	print 'gaussian envelope @',t,'delta_g=',delta_g
	return cfg.a0*np.exp(-np.power(r/cfg.W,2)-np.power((z-(zl0+delta_g*t))/cfg.L,2));

def cos2_envelope(t):
	gamma_g=cfg.k0_kp;
	delta_g=math.sqrt(1-1/(gamma_g*gamma_g))-1;
	
	buf = (z-(zl0+delta_g*t))/cfg.L;
 
	ret = cfg.a0*exp(-np.power(r/cfg.W,2))*np.power(np.cos(0.5*math.pi*buf),2)

	ret[ np.abs(buf) > .99995 ] = 0.0

	return ret





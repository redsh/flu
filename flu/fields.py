import numpy as np
from configs import current_config
from setup import *

#if current_config.geometry == 'cyl':
ax_z,ax_r = (1,0)

cfg = current_config
z,r = current_config.zr()
dz,dr = cfg.dz,cfg.dr
dt = cfg.dt

one__r = r.copy()
print 'fields dr',dr

for i in range(one__r.shape[0]):
	if i != 0:
		one__r[i,:] = 1.0/(dr*i)

dd  = (current_config.dr,current_config.dz)
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
	#retur

	#s0.Ar = cfg.a0*np.exp(-((z-zl0)*(z-zl0))/(cfg.L*cfg.L))*np.exp(-r*r/(cfg.W*cfg.W))

#from stencils.cyl import *

from plot import *




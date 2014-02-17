import math
import numpy as np

current_config = 123

class Config(object):
	def __init__(self):
		pass

	def zr(self):
		lz = np.linspace(-self.Lz*0.5, self.Lz*0.5, self.nz)
		lr = np.linspace(0, self.Lr, self.nr)
		return np.meshgrid(lz,lr)

	def set(self):
		
		self.nz=int((self.Lz)/self.dz+0.5)
		self.dz=(self.Lz)/(self.nz-1)

		self.nr=int((self.Lr)/self.dr+0.5)
		self.dr=(self.Lr)/(self.nr-1)

		self.dt = 0.2*self.dz
		#self.nz,self.nr=self.Lz/self.dz,self.Lr/self.dz
		global current_config
		current_config = self
		
def sample_config():
	r = Config()
	r.geometry = 'cyl'
	r.method   = 'fluid'

	r.l0=1e-4

	r.Lz,r.Lr = 16.,6.
	r.dz,r.dr = 1./20,1./20
	r.L = 2.
	r.W = 2.
	r.k0 = 2.0*math.pi/r.l0
	r.a0 = 1.0
	r.dt = r.dz*0.25

	return r

#l0=1e-4
#Lz,Lr = 120e-4,80e-4
#dz,dr = 0.5e-4,0.5e-4
#L = 10e-4
#W = 4e-4
#k0 = 2.0*math.pi/l0

def envelope_config():
	r = Config()
	r.geometry = 'cyl'
	r.method   = 'fluid'

	r.l0=1e-4

	r.Lz,r.Lr = 60e-4,30e-4
	r.dz,r.dr = 0.5e-4,0.5e-4
	r.L = 10e-4
	r.W = 4e-4
	r.k0 = 2.0*math.pi/r.l0
	r.a0 = 1.0
	r.dt = r.dz*0.2

	r.nz,r.nr=r.Lz/r.dz,r.Lr/r.dz

	return r



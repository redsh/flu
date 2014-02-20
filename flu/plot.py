import math
import numpy as np
from matplotlib.pyplot import *

def plot_field(z,r,a,name='',x=0,columns=1):

	subplot(2,columns,x)
	title(name)

	if len(a.shape) == 1:
		plot(a,'-.')

		return
	#print 'ffff',2,columns,x,2*columns+x
	sgn = 1.0
	if name == 'Er' or name == 'Bp':
		sgn = -1.0

	if len(z.shape) == 2:
		plot(z[1,:],a[1,:]) #r,z

		subplot(2,columns,columns+x)
	
		nr = a.shape[0]
		af = np.zeros((nr*2-1,a.shape[1]))
		af[0:nr,:]     = a[::-1,:]*sgn
		af[nr:nr*2-1,:]= a[1:,:]
		#af[nr,:]=0
		
		imshow(af, aspect='auto', extent=[z.min(), z.max(), r.max(), -r.max()])
		colorbar()
	else:
		plot(z, a[a.shape[0]/2+1,:])

		subplot(2,columns,columns+x)
	
		imshow(a, aspect='auto', extent=[z.min(), z.max(), r.max(), -r.max()])
		colorbar()

def plots(z,r,fields, names=[], supt=''):
	figure(figsize=(20,7))
	suptitle(supt)

	for i,f in enumerate(fields):
		name = ''
		if i < len(names): name = names[i]
		plot_field(z,r,f, name ,i+1,len(fields))

def plots_explicit(z_bounds,r_bounds,fields,names=[],supt=''):
	z = np.arange(z_bounds[0],z_bounds[1],(z_bounds[1]-z_bounds[0])/fields[0].shape[1])
	r = np.arange(r_bounds[0],r_bounds[1],(r_bounds[1]-r_bounds[0])/fields[0].shape[0])
	#print z,r
	plots(z,r,fields,names,supt)

def plot_sim(z,r,t,s):
	figure(figsize=(20,7))

	plots(z,r, \
		[np.sqrt(s.Ar*s.Ar + s.Ai*s.Ai), s.Ez, s.Er, s.Bp, s.rho], \
		['A','Ez','Er','Bp','rho'],str(t) )

	show()

def phys_imshow(z,r,a):
	imshow(a, aspect='auto', extent=[z.min(), z.max(), r.max(), 0])
		
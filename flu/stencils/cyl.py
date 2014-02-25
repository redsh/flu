import numpy as np
# z derivs


def diff_z_1(a):
	ap = np.roll(a, -1, ax_z)
	am = np.roll(a,  1, ax_z)

	ap[:,-1] = 0
	am[:, 0] = 0

	return (ap-am)/(2.0*dz)

#@np.vectorize	
def diff_z_1_upwind(a):
	
	app = np.roll(a, -2, ax_z)
	ap  = np.roll(a, -1, ax_z)

	app[:,-1] = 0
	app[:,-2] = 0
	ap [:,-1] = 0

	return (-app+4.0*ap-3.0*a)/(2.0*dz)

def diff_z_2(a):
	#print dz,ax_z
	ap = np.roll(a, -1, ax_z)
	am = np.roll(a,  1, ax_z)
	##print 'dz1=',dz
    
	return (ap - 2.0*a + am)/(dz*dz)


# r derivs

#@np.vectorize(excluded=['even'])
def diff_r_1(a, even=True):
	ap = np.roll(a, -1, ax_r)
	am = np.roll(a,  1, ax_r)

	if even:
		am[0,:] =  a[1,:] 
	else:
		am[0,:] = -a[1,:] 

	d = (ap-am)/(2.0*dr)
	d[-1,:] = 0
	if even:
		d[0,:] = 0 #??

	return d

#@np.vectorize
def diff_r_2(a,even=True):
	ap = np.roll(a, -1, ax_r)
	am = np.roll(a,  1, ax_r)

	if even:
		am[0,:] =  a[1,:] 
	else:
		am[0,:] = -a[1,:] 

	ap[-1:,:] = 0

	d = (ap - 2.0*a + am)/(dr*dr)
	d[-1,:] = 0

	return d

'''
rr=data.dr*j + data.rmin;
d_Ez[i][0]  +=  2*Bphi[i][1]/data.dr-(Jz+Jz_i);
d_Ez[i][j]  +=  0.5*(Bphi[i][j-1]+Bphi[i][j+1])/rr+(Bphi[i][j+1]-Bphi[i][j-1])*(2dr)
'''

#@np.vectorize
def div_r(a,even=True):
	ap = np.roll(a, -1, ax_r)
	am = np.roll(a,  1, ax_r)
	ap[-1,:] = 0
	am[ 0,:] = 0 #not used

	d = (ap-am)/(2.0*dr) + one__r*a

	if even:
		d[0,:]  = 0
	else:
		d[0,:] = 2.0*a[1,:]/dr

	d[-1,:] = 0

	return d

def div_r_filt(a,even=True):
	#return div_r(a,even)
	#print dr
	ap = np.roll(a, -1, ax_r)
	am = np.roll(a,  1, ax_r)
	ap[-1,:] = 0
	am[ 0,:] = 0 #not used

	d = (ap-am)/(2.0*dr) + one__r*(ap+am)*0.5

	if even:
		d[0,:]  = 0
	else:
		d[0,:] = 2.0*a[1,:]/dr

	d[-1,:] = 0

	return d

#@np.vectorize
def lapl_r(a, even=True):
	d1 = diff_r_1(a,even)
	d = diff_r_2(a,even)+one__r*d1

	d[0,:] = 2.0*d1[1,:]/dr
	d[-1,:]= 0

	return d
	#return div_r(    diff_r_1(r*a, even) , even  )


import numpy as np

def rk4(yn, t, h, f):

	h__2 = h/2

	k1 = f(t     ,yn        )
	k2 = f(t+h__2,yn+k1*h__2)
	k3 = f(t+h__2,yn+k2*h__2)
	k4 = f(t+h   ,yn+k3*h   )
	
	#print dir(yn)
	#print dir(k1)
	#print dir(k1 + k2*2.0)
	#print dir(k1 + k2*2.0 + k3*2.0)
	#print dir(k1 + k2*2.0 + k3*2.0 + k4)

	return yn + (k1 + k2*2.0 + k3*2.0 + k4) * ((1.0/6.0)*h)

def rk2(yn, t, h, f):
	b = yn + (f(t,yn)*(0.5*h))
	return yn + (f(t+0.5*h,b)*h)

def euler(yn, t, h, f):

	return yn + f(t,yn)*h


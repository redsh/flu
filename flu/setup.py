from configs import current_config

import numpy as np
import math


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



def _unpack(obj,glob):
	for p in dir(self):
		if p[0] != '_':
			glob[p] = getattr(obj,p)
			#print p 



import numpy as np
import os,sys,json
if os.getcwd() not in sys.path: sys.path.append(os.getcwd())


from fields import *
from stencils.cyl import *
from integrators import *



import numpy as np
import math

def _unpack(obj,glob):
	for p in dir(obj):
		if p[0] != '_':
			glob[p] = getattr(obj,p)
			#print p 

def _unpack_dict(obj,glob):
	for p in obj:
		if p[0] != '_':
			glob[p] = obj[p]
			#print p 

def setup_module(mod,root):
	for k in root:
		setattr(mod,k,root[k])

import numpy as np
import os,sys,json
if os.getcwd() not in sys.path: sys.path.append(os.getcwd())

from fields import init_fields,gaussian_envelope
import plasma, envelope

from stencils.cyl import *
from integrators import *
from plot import *
from envelope import S,MI

def init(current_config, globals):
	root = {}
	current_config.set()

	init_fields(current_config,root)

	import fields,stencils,stencils.cyl
	setup_module(fields,root)
	setup_module(stencils,root)
	setup_module(stencils.cyl,root)
	setup_module(envelope,root)
	setup_module(plasma,root)

	_unpack(current_config, globals)
	_unpack_dict(root,globals) 


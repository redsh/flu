import numpy as np
from stencils.cyl import *

def Norm(a,b):
    return np.sqrt( np.max( (a-b)*(a-b)  ) )

def Error(x,b):
    return Norm(cfg.k0_kp*cfg.k0_kp*x + diff_z_2(x), b)

def S(A,xi):
    return -0.5*lapl_r(A)-0.5*xi*A

def Jacobi(x,b):
    #print cfg.k0_kp*cfg.k0_kp,2.0/(dz*dz)
    D  = cfg.k0_kp*cfg.k0_kp - 2.0/(dz*dz)
    Rx = np.roll(x,-1,ax_z)/(dz*dz) + np.roll(x,1,ax_z)/(dz*dz)
    return (b-Rx)/D

def MI(guess,b):
    x = guess.copy()
    e = []
    #print dz
    for i in range(200):
        x = Jacobi(x,b)
        e += [ np.log10(Error(x,b)) ]
        
        if len(e) > 5:
            if e[-1]/e[-2] < 1.+1e-3:
                break
    
    #plot(e)
    #show()
    return x

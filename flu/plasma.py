import numpy as np
from stencils.cyl import *
from fields import Sim

def deriv(t,s):
    #boundary: r
    if 1:
        s.ruz[:,-2:] = 0
        s.rur[:,-2:] = 0
        s.Bp[:,-2:] = 0
        s.Ez[:,-2:] = 0
        s.Er[:,-2:] = 0
      
        s.ruz[-2:,:] = 0
        s.rur[-2:,:] = 0
        s.Bp[-2:,:] = 0
        s.Ez[-2:,:] = 0
        s.Er[-2:,:] = 0
        
        s.Bp [0,:]=0
        s.Er [0,:]=0
        s.rur[0,:]=0

    s.rho = ltprofile(t)-diff_z_1_upwind(s.Ez)-div_r_filt(s.Er,False)
   
    s._unpack(globals())
 
    #inizializzazione adiabatica all'inizio
    epsilon = 1.0
    if cfg.T0>0:
        epsilon = 1.-(t/cfg.T0) 
        epsilon*=2
        if epsilon > 1.0: epsilon = 1.0
    #
    
    A2    = (Ar*Ar + Ai*Ai) * (epsilon * epsilon) # envelope
    uz,ur = ruz/rho,rur/rho
    gamma = np.sqrt( 1 + A2*0.5 + uz*uz + ur*ur )
  
    s.gamma = gamma
    
    bz,br  = uz/gamma, ur/gamma
    jz,jr = -ruz/gamma, -rur/gamma
    
    for P in s.plasmas:
        ds_dt,pjz,pjr = P.deriv(s)
        jz += pjz
        jr += pjr

    ds = Sim()
    
    ds.Ez =  diff_z_1_upwind(Ez   ) + div_r_filt(Bp,False) - jz
    ds.Er =  diff_z_1_upwind(Er-Bp)                        - jr
    ds.Bp = -diff_z_1_upwind(Er-Bp) + diff_r_1(Ez)
    
    ds.ruz = diff_z_1_upwind(ruz*(1.0-bz))  - div_r_filt  ( br*ruz,False ) + rho*( - (Ez+Bp*br) - (0.25/(gamma)) * diff_z_1_upwind(A2)  )
    ds.rur = diff_z_1_upwind(rur*(1.0-bz))  - div_r_filt  ( br*rur,True  ) + rho*( - (Er-Bp*bz) - (0.25/(gamma)) * diff_r_1       (A2)  )

    ds.Ar  = 0
    ds.Ai  = 0
    ds.rho = 0
    ds.gamma = 0
    #ds.rho = np.sqrt(1 + A2*0.5)
    
    return ds

#profilo densita'
def ltprofile(t):
    return 1.0

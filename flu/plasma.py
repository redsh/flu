import numpy as np
from stencils.cyl import *
from fields import Sim
import bunch, plot

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

        #for P in s.plasmas:
        #    neg = P[bunch._r,:] < 0
        #    P[bunch._r,neg]  = -P[bunch._r,neg]
        #    P[bunch._ur,neg] = -P[bunch._ur,neg]

    s.rho = ltprofile(t)-diff_z_1_upwind(s.Ez)-div_r_filt(s.Er,False)
   
    s._unpack(globals())
 
    #inizializzazione adiabatica all'inizio
    epsilon = 1.0
    if cfg.T0 < 0:
        # at 0: 1
        # at T1: 0

        T1 = cfg.T0
        if hasattr(cfg,'bunch_adiabatic_T'):
            T1 += cfg.bunch_adiabatic_T

        epsilon = (T1-t)/T1
        #print 't',t,'t0',cfg.T0,'t1=',T1,'dt=',T1-t,'eps',epsilon
        #epsilon*=2
        if epsilon < 0.0: epsilon = 0.0
        if epsilon > 1.0: epsilon = 1.0
        cfg.epsilon = epsilon

    
    A2    = (Ar*Ar + Ai*Ai) * (epsilon * epsilon) # envelope
    uz,ur = ruz/rho,rur/rho
    gamma = np.sqrt( 1 + A2*0.5 + uz*uz + ur*ur )
  
    s.gamma = gamma
    
    bz,br  = uz/gamma, ur/gamma
    jz,jr = -ruz/gamma, -rur/gamma
    
    ds = Sim()
   
    if hasattr(cfg,'bunch_current_factor'):
        if cfg.bunch_current_factor < 0.9999:
            jz *= epsilon
            jr *= epsilon
    #    #print epsilon
    
    if hasattr(cfg,'bunch_corrective_jz'):
        jz = cfg.bunch_corrective_jz
        jr = cfg.bunch_corrective_jr

        ds.plasmas = [0. for P in s.plasmas]
    else:
        ds.plasmas = []

        for P in s.plasmas:
            ds_dt,pjz,pjr = bunch.deriv(P,s)

            #if hasattr(cfg,'bunch_current_factor'):
            #    pjz *= cfg.bunch_current_factor
            #    pjr *= cfg.bunch_current_factor
                #print cfg.bunch_current_factor

            jz += pjz
            jr += pjr

            ds.plasmas += [ds_dt]
        
            cfg.debug.fields.z0_Jzb = jz
            cfg.debug.fields.z1_Jrb = jr
            cfg.debug.fields.z2_rhob = bunch.deposition(P,s,[P[bunch._q,:]])[0]

            cfg.debug.fields.z3_divE = diff_z_1_upwind(s.Ez)+div_r_filt(s.Er,False)
            cfg.debug.fields.z4_charge = cfg.debug.fields.z3_divE-cfg.debug.fields.z2_rhob

            cfg.debug.lineouts['charge_err'] += [ np.abs(cfg.debug.fields.z4_charge).max() ]
            cfg.debug.lineouts['sigma'] += [ np.std(P[bunch._r,:])+np.std(P[bunch._z,:]) ]
            #plot.plot_field(z,r,pjz)
            #plot.show()
            #plot.plot_field(z,r,pjr)
            #plot.show()

    
    ds.Ez =  diff_z_1_upwind(Ez   ) + div_r_filt(Bp,False) - jz
    ds.Er =  diff_z_1_upwind(Er-Bp)                        - jr
    ds.Bp = -diff_z_1_upwind(Er-Bp) + diff_r_1(Ez)
    
    ds.ruz = diff_z_1_upwind(ruz*(1.0-bz))  - div_r_filt  ( br*ruz,False ) + rho*( - (Ez+Bp*br) - (0.25/(gamma)) * diff_z_1_upwind(A2)  )
    ds.rur = diff_z_1_upwind(rur*(1.0-bz))  - div_r_filt  ( br*rur,True  ) + rho*( - (Er-Bp*bz) - (0.25/(gamma)) * diff_r_1       (A2)  )

    ds.Ar  = 0
    ds.Ai  = 0
    ds.rho = 0
    ds.gamma = 0

    if hasattr(cfg,'bunch_current_factor'):
        ds.ruz *= cfg.epsilon
        ds.rur *= cfg.epsilon
    #ds.rho = np.sqrt(1 + A2*0.5)
    
    return ds

#profilo densita'
def ltprofile(t):
    return 1.0

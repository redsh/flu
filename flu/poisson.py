import copy
import numpy as np
import jacobi
#from stencils.cyl import *

def diff_z_2(dz,a):
    #print dz,ax_z
    ap = np.roll(a, -1, 1)
    am = np.roll(a,  1, 1)
    ##print 'dz1=',dz
    d = (ap - 2.0*a + am)/(dz*dz)

    d[0,:] = 0.
    d[-1,:] = 0.
    
    d[:,0] = 0.
    d[:,-1] = 0.

    return d
    
def diff_r_2(dr,a):
    ap = np.roll(a, -1, 0)
    am = np.roll(a,  1, 0)

    d = (ap - 2.0*a + am)/(dr*dr)

    d[0,:] = 0.
    d[-1,:] = 0.
    
    d[:,0] = 0.
    d[:,-1] = 0.

    return d




class PoissonSolver1DGamma(jacobi.JacobiSolver):
    
    def __init__(self,gamma,dr,dz):
        self.gamma = gamma
        self.dr = dr
        self.dz = dz

    def F(self, a):
        dz = self.dz

        ap = np.roll(a, -1, 1)
        am = np.roll(a,  1, 1)
        d = (ap - 2.0*a + am)/(dz*dz)
        
        d[:,0] = 0.
        d[:,-1] = 0.

        ig2 = (1.0/(self.gamma*self.gamma))
        return d*ig2

    def D(self):
        #print 'dz=',dz
        ig2 = (1.0/(self.gamma*self.gamma))
        return -2.0/(self.dz*self.dz)*ig2

    def R(self,rho):
        d = self.F(rho)-self.D()*rho

        d[:,0] = 0.
        d[:,-1] = 0.

        return d

    def down(self,mg):
        r = copy.copy(self)
        r.dz *= 2.0
        r.dr *= 2.0
        return r

    def __init__(self, gamma, dr, dz):
        self.gamma = gamma
        self.dz = dz
        self.dr = dr

        super(PoissonSolver1DGamma,self).__init__()


class PoissonSolver2DGamma(jacobi.JacobiSolver):

    def __init__(self,gamma,dr,dz):
        self.gamma = gamma
        self.dr = dr
        self.dz = dz

    def F(self, rho):
        ig2 = (1.0/(self.gamma*self.gamma))

        d =  diff_z_2(self.dz,rho)*ig2+diff_r_2(self.dr,rho) #lapl_r(rho) + 

        d[0,:] = 0.
        d[-1,:] = 0.
        
        d[:,0] = 0.
        d[:,-1] = 0.
        
        return d

    def D(self):
        #print 'dz=',dz
        ig2 = (1.0/(self.gamma*self.gamma))

        return -2.0/(self.dz*self.dz)*ig2-2.0/(self.dr*self.dr) #-3.0/(2.0*dz)

    def R(self,rho):
        d = self.F(rho)-self.D()*rho

        d[0,:] = 0.
        d[-1,:] = 0.
        
        d[:,0] = 0.
        d[:,-1] = 0.

        return d

    def down(self,mg):
        r = copy.copy(self)
        r.dz *= 2.0
        r.dr *= 2.0
        return r

    def __init__(self, gamma, dr, dz):
        self.gamma = gamma
        self.dr = dr
        self.dz = dz
        
        super(PoissonSolver2DGamma,self).__init__()

def diff_r_1(dr, a, even=True):
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

def lapl_r(one__r, dr, a, even=True):
    d1 = diff_r_1(a,even)
    d = diff_r_2(a,even)+one__r*d1

    d[0,:] = 2.0*d1[1,:]/dr
    d[-1,:]= 0

    return d

class PoissonSolver2DCylGamma(jacobi.JacobiSolver):

    def __init__(self,gamma,dr,dz):
        self.gamma = gamma
        self.dr = dr
        self.dz = dz
        self.one__r = one__r.copy()
        
        super(PoissonSolver2DCylGamma,self).__init__()

    def F(self, rho):
        ig2 = (1.0/(self.gamma*self.gamma))

        return diff_z_2(self.dz,rho)*ig2+lapl_r(self.one__r,dr,self.dr,rho) #lapl_r(rho) + 

    def D(self):
        #print 'dz=',dz
        gamma = self.gamma
        ig2 = (1.0/(gamma*gamma))

        return -2.0/(self.dz*self.dz)*ig2-2.0/(self.dr*self.dr) #-3.0/(2.0*dz)

    def R(self,rho):
        d = self.F(rho)-self.D()*rho

        #d[0,:] = 0.
        #d[-1,:] = 0.
        
        #d[:,0] = 0.
        #d[:,-1] = 0.

        return d

    def down(self, mg):
        r = copy.copy(self)
        r.dz *= 2.0
        r.dr *= 2.0
        r.one__r = mg.restrict(self.one__r)
        return r

        

'''

    def __init__(self,gamma):
        def _F(rho):
            ig2 = (1.0/(gamma*gamma))
            R = diff_z_2(rho)*ig2+lapl_r(rho) #lapl_r(rho) + 
            
            R[:,0] = 0.
            R[:,-1] = 0.
            R[-1,:] = 0.

            return R

        def _F_D():
            #print 'dz=',dz
            ig2 = (1.0/(gamma*gamma))
            return -2.0/(dz*dz)*ig2-2.0/(dr*dr) #-3.0/(2.0*dz)

        def _F_R(rho):
            return _F(rho)-_F_D()*rho

        super(PoissonSolver2DCylindricalGamma,self).__init__(_F, _F_R, _F_D)


'''

#def lapl_r(a, even=True):
#    d1 = diff_r_1(a,even)
#    d = diff_r_2(a,even)+one__r*d1
#
#    d[0,:] = 2.0*d1[1,:]/dr
#    d[-1,:]= 0
#    return d
    #return div_r(    diff_r_1(r*a, even) , even  )



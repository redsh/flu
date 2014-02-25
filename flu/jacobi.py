import numpy as np
from stencils.cyl import *
from plot import *

class IterativeSolver(object):
    def cycle(self,guess,b,max_iters=200, min_iters=10, callback=None, convergence=1e-4):
        if not guess:
            x = np.zeros(b.shape)
        else:
            x = guess.copy()
        e = []

        for i in range(max_iters):
            x,err = self.smooth(x,b)
            
            #e += [err]
            e += [ np.log10(err) ]

            if i > min_iters:
                if e[-1]/e[-2] < 1.+convergence and convergence > 0.:
                    if callback:
                        callback(i,-1,x,e)
            
                    break

        if callback:
            callback(-1,max_iters,x,e)

        return x

class JacobiSolver(IterativeSolver):

    def __init__(self):
        #self.F = F
        #self.D_inv = 1.0/D()
        #self.R = R
        pass

    def norm(self,a,b):
        return np.sqrt( np.max( (a-b)*(a-b)  ) )

    def error(self,x,b):
        return self.norm(self.F(x),b)

    def residual(self,x,b):
        return self.F(x)-b

    def smooth(self,x,b):
        x1 = (b-self.R(x))/self.D()
        return x1, self.error(x,b)

    

class VCycle(IterativeSolver):

    def __init__(self, smoother, max_depth=4):
        self.max_depth = max_depth
        self.pre_smooth = 8
        self.post_smooth = 8
        self.final_smooth = 8
        self.smoothers = [ smoother ]
        self.verbose = False

        for depth in range(1,self.max_depth+1):
            self.smoothers.append(self.smoothers[-1].down(self))

    def restrict(self,x):
        b = 1
        c = b*2
        shape = [(nx-2*b)/2+2*b for nx in x.shape]
        out = np.zeros(shape)

        in_row0 = x[c:-c:2,:]
        in_rowp = x[c+1:-c+1:2,:]
        in_rowm = x[c-1:-c-1:2,:]

        out[b:-b,b:-b] = \
            0.2500*in_row0[:,c+0:-c+0:2]+\
            0.1250*in_row0[:,c+1:-c+1:2]+\
            0.1250*in_row0[:,c-1:-c-1:2]+\
            0.1250*in_rowp[:,c+0:-c+0:2]+\
            0.1250*in_rowm[:,c+0:-c+0:2]+\
            0.0625*in_rowp[:,c-1:-c-1:2]+\
            0.0625*in_rowp[:,c+1:-c+1:2]+\
            0.0625*in_rowm[:,c-1:-c-1:2]+\
            0.0625*in_rowm[:,c+1:-c+1:2]

        return out#[b:-b,b:-b]

    def prolungate(self,x):
        b = 1
        c = b*2
        shape = [(nx-2*b)*2+1+2*b for nx in x.shape]
        out = np.zeros(shape)

        in_row0 = x[b:-b,:]
        in_row1 = x[b+1:,:]
        
        #out_row0 = out[c:-c:2,:]
        #out_row1 = out[c+1:-c+1:2,:]

        x00 = in_row0[:,b:-b]
        x10 = in_row0[:,b+1:]
        x01 = in_row1[:,b:-b]
        x11 = in_row1[:,b+1:]
        
        #print 'ss',out_row0[:,c+0:-c+0:2].shape,x00.shape,x01.shape,x10.shape[1]+1

        out[c:-c:2,c+0:-c+0:2] = x00
        out[c:-c:2,c+1:-c+1:2] = (x10+x00)*0.5
        out[c+1:-c+1:2,c+0:-c+0:2] = (x01+x00)*0.5
        out[c+1:-c+1:2,c+1:-c+1:2] = (x00+x10+x01+x11)*0.25

        return out#[c:-c,c:-c]

    def recurse_odd(self,depth,xi,b):
        #print depth,xi.shape
        smoother = self.smoothers[depth]

        for i in range(self.pre_smooth):
            xi,err = smoother.smooth(xi, b)

        if depth == self.max_depth:
            for i in range(self.final_smooth):
                xi,err = smoother.smooth(xi, b)
            return xi

        residual = smoother.residual(xi,b)

        residual_small = self.restrict(residual)

        d_small = self.recurse(depth+1, None, residual_small)

        d = self.prolungate(d_small)
        
        if self.verbose:
            print depth,smoother.dz
            figure(figsize=(20,4))
            suptitle('depth %d'%depth)
            
            fields = [('xi',xi),('res',residual),('residual_small',residual_small),\
                      ('d_small',d_small),('d',d),('xi-d',xi-d)]
            for i,(n,f) in enumerate(fields):
                subplot(2,len(fields),1+i)
                title(n)
                plt.plot(f[f.shape[0]/2,:])
                subplot(2,len(fields),1+i+len(fields))
                imshow(f)
            show()

        xi = xi - d

        for i in range(self.post_smooth):
            xi,err = smoother.smooth(xi, b)

        return xi

    def recurse(self,depth,xi,b):
        if xi is None:
            xi = np.zeros(b.shape)

        oldshape = xi.shape
        newshape = []
        for k in oldshape:
            newshape += [(k % 2 == 0) + k]

        if list(newshape) == list(oldshape):
            return self.recurse_odd(depth,xi,b)
        else:
            #print 'resizing',depth,oldshape,newshape
            xi_odd = np.zeros(newshape)
            b_odd = np.zeros(newshape)

            xi_odd[0:oldshape[0],0:oldshape[1]] = xi
            b_odd[0:oldshape[0],0:oldshape[1]] = b

            xi_odd = self.recurse_odd(depth,xi_odd,b_odd)

            return xi_odd[0:oldshape[0],0:oldshape[1]]

    def smooth(self,guess,b):
        x = self.recurse(0,guess,b)
        err = self.smoothers[0].error(x,b)
        return x,err




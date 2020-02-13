R"""
classes for integrating second-order differential equations of the 
form u''(r) = f(r)u(r) on different radial grids. For the equidistant
grid the very efficient Numerov method is used.
"""
from numpy import ediff1d, hstack, zeros

def is_equidistant(r):
    from numpy import mean
    h = ediff1d(r)
    if abs(sum(h-mean(h))) <= 1.0e-13:
        return True
    else:
        return False

class RadialIntegrator:
    def __init__(self):
        pass
    def integrate_outward(self,r, f, u0out, u1out, transformBack=True):
        # overload this method
        pass
    def integrate_inward(self,r, f, u1in, u0in, transformBack=True):
        # overload this method
        pass
    def integrate_match(self, R, F, u0out, u1out, u1in, u0in, imatch, transformBack=True):
        """
        Input:
        ======
        R: grid, which can be any dimension, has to be equidistant
        F: can be either a numpy array with values of f on the grid
          or a function F(r) that computes theses values
        u0out, u1out: starting values for outward integration
        u1in, u0in: starting values for inward integration, u1in and u0in 
           are in opposite order because the iteration proceeds from right
           to left on the grid.
        imatch: index of radial point rmatch=r[imatch] where to match the
           values and first derivatives from the outward and the inward integration
        transformBack: boolean
           For finding the eigenvalues the wave function has to be known
           only around the matching point. When looking for eigenvalues
           the inverse transformation w --> u does not have to be performed.
   
        Output:
        =======
        g: mismatch at rmatch
        node_count: number of nodes of the solution
        U: values of u on the grid (if transformBack is True, otherwise None)
        """
        r = R.ravel()
        assert is_equidistant(r)
        indx = r.argsort()
        r = r[indx]
        rout = r[:imatch+1]
        rin =  r[imatch:]

        f = F.ravel()[indx]
        fout = f[:imatch+1]
        fin =  f[imatch:]

        uout,nodes_out = self.integrate_outward(rout, fout, u0out, u1out, transformBack=transformBack)
        uin,nodes_in = self.integrate_inward(rin, fin, u1in, u0in, transformBack=transformBack)
        
        node_count = nodes_out + nodes_in
        """match the solutions"""
        #print "Match solutions at r=%s, imatch=%s" % (r[imatch], imatch)
        scale_fac = uin[0]/uout[-1]
# The following two lines caused the bug, that only bound solutions with even node counts were found!
#        if abs(scale_fac < 1.0e-15):
#            scale_fac = 1.0e-15
#        print "scale_fac = %s" % scale_fac
#        print "uin[0] = %s, uout[-1] = %s" % (uin[0], uout[-1])
        g = scale_fac*(uout[-1]-uout[-2]) - (uin[1]-uin[0])

        if transformBack:
            u = hstack( (scale_fac*uout[:-1], uin) )
            U = (u[indx]).reshape(R.shape)
        else:
            U = None
        return g,node_count,U

class NumerovIntegrator(RadialIntegrator):
    """for equidistant grids"""
    def integrate_outward(self, r, f, u0out, u1out, transformBack=True):
        h = ediff1d(r, to_end=r[-1]-r[-2])
        w = zeros(r.shape)
        node_count = 0
        w[0] = (1.0 - 1.0/12.0*pow(h[0],2) * f[0])*u0out
        w[1] = (1.0 - 1.0/12.0*pow(h[1],2) * f[1])*u1out    
        # r,f and u are flat arrays
        for i,ri in enumerate(r):
            if i < 2:
                # first two values given as initial conditions
                continue
            w[i] = 2.0 * (1.0 + 5.0/12.0*pow(h[i],2)*f[i-1])/(1.0 - 1.0/12.0*pow(h[i],2)*f[i-1]) * w[i-1] - w[i-2]
            if w[i]*w[i-1] < 0.0:
                node_count += 1
        if transformBack:
            """for matching the inward and outward solutions only the values at u[0]
            and u[1] are needed"""
            u = 1.0/(1.0 - pow(h,2)/12.0 * f) * w
        else:
            """only transform the points back which are needed for matching"""
            u = zeros(r.shape)
            u[-2:] = 1.0/(1.0 - pow(h[-2:],2)/12.0 * f[-2:]) * w[-2:]
        u[0] = u0out
        u[1] = u1out
        return u, node_count
    
    def integrate_inward(self, r, f, u1in, u0in, transformBack=True):
        h = ediff1d(r, to_end=r[-1]-r[-2])
        w = zeros(r.shape)
        node_count = 0
        w[-1] = (1.0 - 1.0/12.0*pow(h[0],2) * f[-1])*u0in
        w[-2] = (1.0 - 1.0/12.0*pow(h[1],2) * f[-2])*u1in    
        for j in range(0, len(r)):
            if j < 2:
                # first two values given as initial conditions
                continue
            i = len(r)-j-1
            w[i] = 2.0 * (1.0 + 5.0/12.0*pow(h[i],2)*f[i+1])/(1.0 - 1.0/12.0*pow(h[i],2)*f[i+1]) * w[i+1] - w[i+2]
            if w[i]*w[i+1] < 0.0:
                node_count += 1
        if transformBack:
            u = 1.0/(1.0 - pow(h,2)/12.0 * f) * w
        else:
            """only transform the points back which are needed for matching"""
            u = zeros(r.shape)
            u[:2] = 1.0/(1.0 - pow(h[:2],2)/12.0 * f[:2]) * w[:2]
        u[-1] = u0in
        u[-2] = u1in
        return u, node_count

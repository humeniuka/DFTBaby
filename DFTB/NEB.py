"""
Nudged elastic band (NEB) method for finding minimum energy paths (MEP) and saddle points.
Implementation based on 
     "Improved tangent estimate in the nudged elastic band method for finding minimum energy paths and saddle points" by Henkelman, G.; Jonsson, H. J.Chem.Phys. 113, 9978 (2000) 
"""
from DFTB import XYZ
from numpy import zeros, cos, sin, pi, linspace, array, dot, vstack, cumsum, argmin, frompyfunc, sign
import numpy as np
from numpy.linalg import norm
import numpy.linalg as la
from scipy import special

class NEB:
    def __init__(self, force_constant=1.0, force_constant_surface_switch=5.0, mass=1.0, nr_processors=1):
        self.force_constant = 1.0
        self.force_constant_surface_switch = force_constant_surface_switch
        self.mass = 1.0
        self.nr_processors = nr_processors
        self.V = []
        self.istep = 0
    def setImages(self, images, states):
        """
        Parameters:
        ===========
        images: list of numpy arrays with the coordinates of each image. If no intermediate states
           along the reaction coordinate are known, this has to be at least [initial_position, final_positions]
        states: list of indeces indicating on which electronic state the images reside
        """
        #self.R = [array(img) for img in images] # make sure we have numpy arrays
        self.R = images # vectors of atom positions
        self.states = states
    def addImagesLinearly(self, nimg=10):
        """
        between segments interpolate linearly. Add nimg additional images
        to each segments. If there are initially only two segments, the products
        and the educts, this creates a path of lengths nimg+2.
        """
        path = []
        st = []
        for i,Ri in enumerate(self.R[:-1]):
            for j,a in enumerate(linspace(0.0, 1.0, nimg)):
                Rinterp = (1.0-a)*self.R[i] + a*self.R[i+1]
                # 
                if j < nimg/2.0:
                    st.append( self.states[i] )
                else:
                    st.append( self.states[i+1] )
                if j < nimg-1: # do not append the last image
                    path.append(Rinterp)
        path.append(self.R[-1])
        st.append(self.states[-1])
        print "initial guess for path contains %s images" % len(path)
        self.R = path
        self.states = st
    def setEnergyCalculator(self, ecalc):
        """
        determine by which method energies and gradients are calculated

        Parameters:
        ==========
        calc: object with method
                getEnergiesAndGradient(R, state)
        """
        self.ecalc = ecalc
    def _getPES(self):
        """calculate energies and gradients at the image positions"""
        self.V = [None for Ri in self.R]  # potential energies of each image
        self.F = [None for Ri in self.R]  # true forces acting on each image
        for i,Ri in enumerate(self.R):
            print "image %d" % i
            energies, grad = self.ecalc.getEnergiesAndGradient(Ri, self.states[i])
            self.V[i] = energies[self.states[i]]
            self.F[i] = - grad
    def _getTangents(self):
        """tangents along the path at the image positions"""
        self.tangents = [None for Ri in self.R]
        for i in xrange(1,len(self.R)-1):
            if self.V[i-1] <= self.V[i] <= self.V[i+1]:
                self.tangents[i] = self.R[i+1] - self.R[i]
            elif self.V[i+1] < self.V[i] < self.V[i-1]:
                self.tangents[i] = self.R[i] - self.R[i-1]
            else:
                dVmax = max(abs(self.V[i+1] - self.V[i]),abs(self.V[i-1] - self.V[i]))
                dVmin = min(abs(self.V[i+1] - self.V[i]),abs(self.V[i-1] - self.V[i]))
                taup = self.R[i+1] - self.R[i]
                taum = self.R[i] - self.R[i-1]
                if self.V[i+1] > self.V[i-1]:
                    self.tangents[i] = taup*dVmax + taum*dVmin
                elif self.V[i+1] < self.V[i-1]:
                    self.tangents[i] = taup*dVmin + taum*dVmax
                else:
                    self.tangents[i] = self.R[i+1] - self.R[i-1] # ?
            # normalize tangents
            self.tangents[i] /= norm(self.tangents[i])
    def _getEffectiveForces(self):
        # effective total force
        self.effF = [None for Ri in self.R]
        for i in xrange(1,len(self.R)-1):
            # spring force parallel to tangents
            km = self.force_constant
            if self.states[i+1] != self.states[i]:
                kp = self.force_constant_surface_switch
                dR = self.R[i+1] - self.R[i]
                dE = self.V[i+1] - self.V[i]
                dF = self.F[i+1] - self.F[i]

                F1 = dR + (dE - np.dot(self.F[i+1], dR)) * self.F[i+1]
                F2   = -dR + (-dE + np.dot(self.F[i], dR)) * self.F[i]
                print "dE = %s" % dE
                print "erf(dE) = %s" % special.erf(dE)
                print "Fi+1 = %s" % self.F[i+1]
                print "Fi   = %s" % self.F[i]
                Fspring = kp * (F1 + F2)
                print "Fspring = %s" % Fspring
            else:
                kp = self.force_constant
                Fspring = kp * norm(self.R[i+1] - self.R[i]) * self.tangents[i] # new implementation by Henkelman/Jonsson
                                
            if self.states[i] != self.states[i-1]:
                km = self.force_constant_surface_switch
                Fspring -= dot(km*(self.R[i] - self.R[i-1]), self.tangents[i]) * self.tangents[i] # from original implementation of NEB
            else:
                km = self.force_constant
                Fspring -= km*norm(self.R[i] - self.R[i-1]) * self.tangents[i] # new implementation by Henkelman/Jonsson
            #Fspring = dot(self.force_constant*((self.R[i+1] - self.R[i]) - (self.R[i] - self.R[i-1])), self.tangents[i]) * self.tangents[i] # from original implementation of NEB
            # perpendicular component of true forces
            Fnudge = self.F[i] - dot(self.F[i],self.tangents[i])*self.tangents[i]
            self.effF[i] = Fspring + Fnudge
        if self.optimize_endpoints == True:
            # initial and final weights move toward minima
            self.effF[0] = self.F[0]
            self.effF[-1] = self.F[-1]
        else:
            # supress force on ends so that they stay put
            self.effF[0] = zeros(self.F[0].shape)
            self.effF[-1] = zeros(self.F[-1].shape)
    def _converged(self, tolerance):
        """Check if average forces have dropped below certain threshold"""
        self.avgForce = 0.0
        for i in xrange(1,len(self.R)-1):
            self.avgForce += norm(self.effF[i])
        if self.optimize_endpoints == True:
            # force on enpoints should only add to the convergence measure
            # is they can be optimized
            self.avgForce += norm(self.effF[0])
            self.avgForce += norm(self.effF[-1])
        self.avgForce /= len(self.R)
        print "average force = %2.5f (tolerance = %2.5f)" % (self.avgForce, tolerance)
        if self.avgForce < tolerance:
            return True
        else:
            return False
    def findMEP(self, tolerance=0.001, nsteps=1000, dt=0.01, friction=0.05, optimize_endpoints=True):
        """
        find minimum energy path by allowing the band to move with some friction to
        damp out oscillations. After the optimization the optimally distributed images
        can be retrieved with .getImages().

        Parameters:
        ===========
        tolerance: If the average force on the band drops below this value, the optimization is stopped
        nsteps: run damped dynamics for nsteps.
        dt: time step
        friction: damping coefficient between 0.0 (no damping) and 1.0 (do not use such high damping)
        optimize_endpoints: if True, end points of band move towards energy minima.
        """
        self.optimize_endpoints = optimize_endpoints
        Rlast = [Ri for Ri in self.R] # R(t-dt), R[0] and R[-1] stay always the same
        Rnext = [Ri for Ri in self.R] # R(t+dt)
        for self.istep in xrange(0, nsteps):
            self._getPES()
            self._getTangents()
            self._getEffectiveForces()
            # optimized positions of intermediate images
            # and minimize positions of ends
            for i in xrange(0, len(self.R)):
                if i in [0, len(self.R)-1] and (self.optimize_endpoints == False):
                    # as effF[0] = 0 and effF[-1] = 0 this line should not be neccessary ????
                    continue
                if self.istep == 0:
                    # Euler step, without initial velocity
                    Rnext[i] = self.R[i] + 0.5*self.effF[i]/self.mass*pow(dt,2)
                    Rlast[i] = self.R[i]
                    self.R[i] = Rnext[i]
                else:
                    # damped Verlet algorithm
                    Rnext[i] = (2.0-friction)*self.R[i] - (1.0-friction)*Rlast[i] + self.effF[i]/self.mass*pow(dt,2)
                    Rlast[i] = self.R[i]
                    self.R[i] = Rnext[i]
            if self._converged(tolerance):
                break
            self._writeIteration()
            self.plot()
        #
        else:
            raise Warning("Could not find minimum energy path in %s iterations (average force = %2.5f > tolerance = %2.5f)." % (self.istep+1, self.avgForce, tolerance))
    def _writeIteration(self):
        import sys
        print "Iteration = %s  " % self.istep
        sys.stdout.flush()
    def plot(self):
        if self.istep % 10 == 0:
            if hasattr(self.ecalc, "plot"):
                self.ecalc.plot(images=self.getImages(), energies=self.V, istep=self.istep)
#
            if hasattr(self.ecalc, "plot3d"):
                self.ecalc.plot3d(images=self.getImages(), states=self.getImageStates())

    def getImages(self):
        return self.R
    def getImageStates(self):
        return self.states
    def splineMEProfile(self):
        """
        interpolate the energy along MEP between images with a cubic spline.
        
        Returns:
        ========
        me: callable function, that returns the energy as a function of the reaction coordinate
          which runs from 0 (educts) to 1 (product)
        """
        n = len(self.R)
        x = zeros(n)
        f = zeros(n)
        f1 = zeros(n) # first derivative of f
        for i in range(0,n):
            if i == 0:
                t = self.R[1] - self.R[0]
                x[i] = 0.0
            else:
                t = self.R[i] - self.R[i-1]
                x[i] = x[i-1] + norm(t)
            f[i] = self.V[i]
            f1[i] = dot(self.F[i], t/norm(t))
        def MEfunc(rxc):
            """
            potential energy along minimum energy path.

            Parameters:
            ===========
            rxc: reaction coordinate (between 0.0 and 1.0)
            """
            assert 0.0 <= rxc <= 1.0
            s = rxc*x[-1]
            ic = argmin(abs(x - s))
            if s < x[ic]:
                ic -= 1
            if s == x[ic]:
                return f[ic]
            dx = norm(x[ic+1]-x[ic])
            df = f[ic+1]-f[ic]
            a = (s-x[ic])/dx
            assert 0.0 <= a <= 1.0
            fa = (1-a)*f[ic] + a*f[ic+1] + a*(1-a)*((1-a)*(-f1[ic]*dx-df) + a*(+f1[ic+1]*dx + df))
            return fa
        return frompyfunc(MEfunc,1,1)

    def splineMEPath(self):
        """
        interpolate the minimum energy path with a cubic spline along the images.
        
        Returns:
        ========
        mep: callable function, that returns the geometry as a function of the reaction coordinate
          which runs from 0 (educts) to 1 (product)
        """
        n = len(self.R)
        x = zeros(n)
        images = self.getImages()
        gradients =  [-f for f in self.F]
        for i in range(0,n):
            if i == 0:
                t = self.R[1] - self.R[0]
                x[i] = 0.0
            else:
                t = self.R[i] - self.R[i-1]
                x[i] = x[i-1] + norm(t)
        def MEPfunc(rxc):
            """
            linearly interpolated geometries along the minimum energy path.

            Parameters:
            ===========
            rxc: reaction coordinate (between 0.0 and 1.0)
            """
            assert 0.0 <= rxc <= 1.0
            s = rxc*x[-1]
            ic = argmin(abs(x - s))
            if s < x[ic]:
                ic -= 1
            if s == x[ic]:
                return images[ic]
            dx = norm(x[ic+1]-x[ic])
            assert x[ic] <= s <= x[ic+1]
            a = (s-x[ic])/dx
            assert 0.0 <= a <= 1.0
            geom = (1-a)*images[ic] + a*images[ic+1]
            return geom
        return MEPfunc

surf_colours = ["blue", "green"]
beads_colours = ["black", "red"]

class PEScalculator:
    def __init__(self, ax,ay,axy):
        self.ax = ax
        self.ay = ay
        self.axy = axy
    def getEnergiesAndGradient(self, x, gradient_state):
        self._setGeometry(x)
        en = self._getEnergy(gradient_state)
        grad = self._getGradient(gradient_state)
        energies = [en]
        return energies, grad
    def _setGeometry(self, R):
        self.R = R
    def _getEnergy(self, state):
        x,y = self.R[0],self.R[1]
        return self._potential(x,y)[state]
    def _getGradient(self, state):
        x,y = self.R[0],self.R[1]
        F = zeros(2)
        F[0] = -2*pi*sin(2*pi*x) * (1 + sin(2*pi*y))
        F[1] = -2*pi*sin(2*pi*y) + 2*pi*cos(2*pi*x)*cos(2*pi*y)
        return F
    ###
    def _potential(self,x,y):
        V = cos(2*pi*x) + cos(2*pi*y) + cos(2*pi*x)*sin(2*pi*y)
        return [V]        
    def plot(self, images=[]):
        from numpy import outer, ones, linspace
        import matplotlib.pyplot as plt
        plt.ion()
        plt.clf()
        N = 100
        x = outer(linspace(-1.0, 1.0, N), ones(N))
        y = outer(ones(N), linspace(-1.0, 1.0, N))
        z = self._potential(x,y)[0]
        plt.contour(x,y,z, zdir='z', levels=linspace(z.min(), z.max(), 20))
        # plot images individually
        for R in images:
            x,y = R[0],R[1]
            plt.plot(x,y, "o", color="red")
        if len(images) > 1:
            # segments between images
            xy = vstack(images).transpose()
            plt.plot(xy[0], xy[1], ls="-", color="black")
#        plt.colorbar()
        plt.draw()
    def plot3d(self, images=[], states=[]):
        from numpy import outer, ones, linspace
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import axes3d
        plt.ion()
        if not hasattr(self, "axes"):
            fig = plt.figure()
            self.axes = fig.add_subplot(111, projection='3d')
        print self.axes
        self.axes.cla()
        N = 100
        lim = 1.0
        x = outer(linspace(-lim, lim, N), ones(N))
        y = outer(ones(N), linspace(-lim, lim, N))
        zs = self._potential(x,y)
        for i,z in enumerate(zs):
            self.axes.plot_wireframe(x,y,z, color=surf_colours[i])
            # plot images individually
        for R,st in zip(images, states):
            x,y = R[0],R[1]
            z = self._potential(x,y)[st]
            self.axes.plot([x],[y],[z], "o", color=beads_colours[st])
        if len(images) > 1:
            # segments between images
            xs, ys, zs = [], [], []
            for im,st in zip(images,states):
                x,y = im
                z = self._potential(x,y)[st]
                xs.append(x)
                ys.append(y)
                zs.append(z)
            self.axes.plot(xs, ys, zs, ls="-", color="black")
#        plt.colorbar()
        plt.draw()


class Morse2D(PEScalculator):
    def __init__(self, x0,y0, a, D, V0):
        self.x0 = x0
        self.y0 = y0
        self.a = a
        self.D = D
        self.V0 = V0
    def _potential(self,x,y):
        R = np.sqrt((x-self.x0)**2 + (y-self.y0)**2)
        Veduct = self.V0 + self.D * pow(1.0 - np.exp(-self.a * R),2)
        return Veduct

class HarmonicOscillator2D(PEScalculator):
    def __init__(self, x0, y0, a):
        self.x0 = x0
        self.y0 = y0
        self.a = a
    def _potential(self,x,y):
        V = self.a*( pow(x-self.x0,2) + pow(y-self.y0,2) )
        return V
        
class ValenceBondPotential2D(PEScalculator):
    def __init__(self, PESproduct, PESeduct):
        self.PESproduct = PESproduct
        self.PESeduct = PESeduct
    def _potential(self,x,y):
        Vp = self.PESproduct._potential(x,y)
        Ve = self.PESeduct._potential(x,y)
        # coupling is strong where energy levels are close
        xTS, yTS = 0.0, 0.0
        VeTS = self.PESeduct._potential(xTS,yTS)
        VpTS = self.PESproduct._potential(xTS,yTS)
        VTS = 0.5*(VeTS + VpTS)
        A = (Vp - VTS)*(Ve - VTS)
        coupling = 10.0*x*y * np.exp(-0.5*(pow(x-xTS,2) + pow(y-yTS,2)))
        # valence bond hamiltonian

        def diag(V11,V12,V21,V22):
            Hvb = np.array([[V11,  V12],
                            [V21,  V22]])
        
            w,v = la.eigh(Hvb)
            sort_indx = np.argsort(w)
            w = w[sort_indx]
            v = v[:,sort_indx]
            return w[0], w[1]

        udiag = np.frompyfunc(diag, 4, 2)
        
        Vadiab = udiag(Vp,coupling,coupling,Ve)
        return Vadiab

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    task = "TS"
    if task == "TS":
        # find transition state
        pes = PEScalculator(1,1,1)
        eps = 0.1
        educt = array([-0.5-eps, 0.375+eps])
    #    product = array([0.5+eps, 0.375-eps])
        product = array([0.5, -0.63])
        
        neb = NEB()
        neb.setEnergyCalculator(pes)
        neb.setImages([educt,product], states=[0,0])
        neb.addImagesLinearly(15)
        neb.findMEP(tolerance=0.02)

        me = neb.splineMEProfile()
        mep = neb.splineMEPath()
        rxc = linspace(0.0, 1.0, 30)
        
        plt.clf()
    #    plot(rxc, me(rxc))
#        pes.plot(images=[mep(x) for x in rxc])
#
        pes.plot3d(images=[mep(x) for x in rxc], states=[0 for x in rxc])
#
        plt.ioff()
        plt.show()
    elif task == "CI":
#        eductPES = Morse2D(-0.5, -0.5, 5.0, 5.0, 0.0)
#        productPES = Morse2D(+0.5, +0.5, 5.0, 5.0, 0.4)
        eductPES = HarmonicOscillator2D(-0.5,-0.5, 3.0)
        productPES = HarmonicOscillator2D(+0.5,+0.5, 3.0)
        vb = ValenceBondPotential2D(eductPES, productPES)

        educt = np.array([-0.5, -0.5])
        product = np.array([+0.5, +0.5])

        neb = NEB()
        neb.setEnergyCalculator(vb)
        neb.setImages([educt,product], states=[1,0])
        neb.addImagesLinearly(15)
        neb.findMEP(tolerance=0.2)

#        vb.plot3d()

        plt.ioff()
        plt.show()


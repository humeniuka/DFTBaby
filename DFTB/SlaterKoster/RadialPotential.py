"""
The module defines classes that allow to solve the radial Schroedinger equation 

  [ 1/2 d^2/dr^2 + (E - V(r) - l(l+1)/(2 r^2)) ] u_l =  0

for different potentials V(r) by inward and outward integration and matching 
to find the eigenenergies (shooting method)
(see Thijssen, "Computational Physics" chapter 2)

Note that for some elements such as gold the relativistic
Schroedinger equation has to be used which cannot be treated 
with Numerov's method because they it contains d/dr u_l(r).
"""

from numpy import *
from scipy.special import jn, yn
import mpmath
import RadialIntegrator
import sys
from DFTB.AtomicData import nuclear_charge

def bisect(f, x1, x2):
    """between x1 and x2 f(x) changes sign."""
    f1 = f(x1)
    f2 = f(x2)
    if f1*f2 > 0.0:
        raise ValueError("f does not change sign between x1 = %s and x2 = %s, f(x1) = %s and f(x2) = %s" % (x1, x2, f1,f2))
    if f1 < 0.0:
        a = x1
        b = x2
    else:
        a = x2
        b = x1
    while True:
        midpoint = (a+b)/2.0
        fm = f(midpoint)
        if fm < 0.0:
            a = midpoint
        else:
            b = midpoint
        yield fm,(a,b)

class RadialPotential(object):
    eps = 1.0e-7 # to remove singularity in l*(l+1)/r^2
    integrator = RadialIntegrator.NumerovIntegrator()
    #### virtual funcitions which have to be overloaded by derived classes ###
    def V(self, r):
        pass
    def F(self, r, E, l):
        f = 2.0*(self.Vr - E) + l*(l+1)/pow(r+self.eps,2)
        return f
    def zero_limit(self, r, E, l):
        """analytical solution for small r for initializing outward integration"""
        pass
    def asymptotic_decaying(self, r, E, l):
        """analytical solution to radial Schroedinger eq. for large r for
        initializing inward iteration, for bound states"""
        pass
    def asymptotic_oscillating(self, r, E, l):
        """asymptotic solution for positive energies. The general
        solution is a combination of two independent solutions with
        coefficients fixed by the scattering phase shift."""
        pass
    ##### access methods #####
    def getRadialGrid(self):
        return self.r
    def getVolumeElement(self):
        """
        gives the differential volume element dr^3
        for integration. If the grid was transformed
        this might not be constant, e.g. for exponential grid
        dr^3 = 4.0*pi*exp(r) * dr
        """
        return 4.0*pi*pow(self.r,2)
    ##### methods common to most potentials ####
    def setIntegrator(self, integrator):
        """
        determine the method by which the radial Schroedinger eq. is solved.
        Parameters:
        ===========
        integrator: instance of NumerovIntegrator or RungeKuttaIntegrator
        """
        self.integrator = integrator
    def setRadialGrid(self, r):
        """
        set the radial grid used for solving the radial Schroedinger equation 
        by integrating from rmin=r[0] outward to rmax=r[-1]. The potential
        is assumed to be zero after rmax.

        Parameters:
        ===========
        r: 1D numpy array, uniform grid
        """
        self.r = r
        self.h = ediff1d(self.r, to_end=self.r[-1]-self.r[-2])
    def getRadialPotential(self):
        self.Vr = self.V(self.r)
        return self.Vr
    def getAngularPotential(self, l):
        return 0.5 * l*(l+1)/pow(self.r+self.eps,2)
    def getEffectivePotential(self, l):
        Veff = self.Vr + self.getAngularPotential(l) 
        return Veff
    def classical_turning_points(self, E, l):
        tpts = []
        Veff = self.getEffectivePotential(l)
        """find minimum of potential, the outer classical turning point
        has to lie to the right of it"""
        imin = Veff.argmin()
        for i in range(imin, len(Veff)-2):
            if (Veff[i-1]-E)*(Veff[i+1]-E) < 0.0:
                """sign change reveals root of V(r)-E"""
                tpts.append(i)
        return tpts
    def bound_state(self, E, l, normalize=False, **opts):
        """
        E: energy of state
        l: angular momentum
        """
        r = self.r
        # The solutions from outward (r=0 --> r=oo) and inward  (r=0 <-- r=oo) integration
        # are matched at the r_match = r[imatch]. If the energy is an eigenvalue, the solution and its
        # r-derivative can be matched continuously.

        # r_match can be provided as an option ... 
        rmatch = opts.get("rmatch", None)
        # ... or it is chosen as the outmost turning point, where V(r) = E
        if rmatch is None:
            tpts = self.classical_turning_points(E, l)
            if len(tpts) == 0 or tpts == [0]:
                # For d-orbitals it turns out, that bound states seem to be possible
                # even if there is no classical turning point???
                tpts = [len(r)/2]  # match in the middle
                #return None
            # select outermost turning point
            imatch = tpts[-1]
        else:
            assert 0 < rmatch < r.max()
            imatch = argmin(abs(r - rmatch))
            del opts["rmatch"]
            
        f = self.F(r, E, l)
        u0out = self.zero_limit(r[0], E, l)
        u1out = self.zero_limit(r[1], E, l)
        u0in = self.asymptotic_decaying(r[-1], E, l)
        u1in = self.asymptotic_decaying(r[-2], E, l)
        gE,node_count,U = self.integrator.integrate_match(r, f, u0out, u1out, u1in, u0in, imatch, **opts)
        gE /= max(U)  # normalize mismatch, otherwise we cannot 
                      # define universal convergence thresholds
                      
        if opts.get("transformBack", True):
            if normalize == True:
                # normalize wave function to 1
                Umax = U.max()
                U_rescaled = U/Umax
                norm = sqrt(sum(U_rescaled**2 * self.h))
                U = U_rescaled / norm
            return gE,node_count, U
        else:
            return gE,node_count, None               
    def match_scattering_solution(self, E, l, r1, r2, ul1, ul2):
#        print "r1 = %s, r2 = %s" % (r1,r2)
        solA1, solB1 = self.asymptotic_oscillating(r1, E, l)
        solA2, solB2 = self.asymptotic_oscillating(r2, E, l)
        Det = solA1*solB2 - solA2*solB1
        a = (solB2 * ul1/r1 - solB1 * ul2/r2)/Det
        b = (-solA2* ul1/r1 + solA1 * ul2/r2)/Det
        """scattering phase"""
        phase_shift = arctan2(-b,a)
        """scale factor"""
        scale_fac = 1.0/sqrt(a*a + b*b)
        return phase_shift, scale_fac
    def scattering_state(self, E, l):
        """
        Find scattering state of this potential.

        Parameters:
        ===========
        E: energy of continuum state
        l: angular momentum

        Returns:
        ========
        delta_l: phase shift
        u_l: radial wave function u_l(r)=r*R_l(r) on the grid
        """
        r = self.r
        f = self.F(r, E, l)

        u0out = self.zero_limit(r[0], E, l)
        u1out = self.zero_limit(r[1], E, l)

        u,node_count = self.integrator.integrate_outward(r, f, u0out, u1out)

        k = sqrt(2.0*E)
        deBroglieLambda = 2.0*pi/k
#        print "deBroglie-wavelength lambda = %s." % deBroglieLambda
#        print "Make sure that the radial grid covers at least 2*lamda."

        """evaluate wave function at two distinct points larger than rmax"""
        n = len(r)
        while u[n-1] == 0.0:
            n -= 1

        delta_l, scale_fac = self.match_scattering_solution(E, l, r[n-1], r[n-2], u[n-1], u[n-2])
#        print "scale_fac = %s" % scale_fac
        return delta_l, scale_fac*u
    def discrete_spectrum_angmom(self, l, energy_range, eps_conv=1.0e-8):
        """
        Find all radial wavefunctions u(r) = r*R(r) with angular momentum l and eigen
        energies within energy_range.

        Parameters:
        ===========
        l: angular momentum
        energy_range: numpy array of energies, the intervals
           between subsequent energies are searched
           for matching conditions by bisection.
        eps_conv: if the mismatch of derivatives from outward
           and inward integration at the turning points is smaller 
           than this constant stop bisection.

        Returns:
        ========
        energy_scan: numpy array of all energies that were scanned
        penalty: numpy array of values of mismatch at scanned energies
        spectrum: numpy array of eigenenergies that were found
           spectrum[i] contains the i-th eigenvalue within energy_range
        nodes: nodes[i] contains number of sign changes in i-th wavefunction,
           this helps checking if you missed some eigenenergies in between.
        wavefunctions: list of radial wavefunctions u(r), wavefunctions[i] contains
           i-th radial wavefunction u[i] = r*R[i] belonging to eigen energy spectrum[i]   
        """
        self.node_count = 0
        spectrum = []
        nodes = []
        wavefunctions = []
        def f(E):
            ret = self.bound_state(E, l)#, transformBack=False)
            if ret == None:
                """no bound state possible"""
                gE = inf
            else:
                gE, self.node_count, U = ret
                """
                from matplotlib.pyplot import ion, cla, plot, draw
                tp = self.classical_turning_points(E, l)[-1]
                ion()
                cla()
                plot(self.r, U)                # linear grid
                plot([self.r[tp]], [U[tp]], "o", markersize=20)
#                plot(exp(self.r), sqrt(exp(self.r))*U) # exponential grid
#                plot([exp(self.r)[tp]], [(sqrt(exp(self.r))*U)[tp]], "o", markersize=20)
                draw()
                """
#            print "%s %s" % (E, gE)
            return gE
        elast = energy_range[0]
        energy_scan = [elast]
        flast = f(elast)
        penalty = [flast]
        for i in range(1, len(energy_range)):
            e1 = elast
            f1 = flast
            e2 = energy_range[i]
            f2 = f(e2)
            elast = e2
            flast = f2
            energy_scan.append(e2)
            penalty.append(f2)

            if f1*f2 < 0.0:
                print "sign change in interval (%s,%s) from %s to %s" % (e1,e2, f1,f2)
                for fm,(a,b) in bisect(f,e1,e2):
                    energy_scan.append((a+b)/2.0)
                    penalty.append(fm)
                    if abs(fm) < eps_conv:
                        print "eigenvalue found at E = %s" % ((a+b)/2.0)
                        print "node_count = %s" % self.node_count
                        if self.node_count != len(spectrum):
                            if self.node_count == len(spectrum)-1:
                                print "Found eigenstate with the same node count twice !?"
                                break
                            else:
                                msg = "Incorrect node count.\nSome eigenvalues have probably been missed.\nTry Increasing range and resolution of energy range."
                                print msg
                                #raise Exception(msg)
                        E = (a+b)/2.0
                        spectrum.append(E)
                        nodes.append(self.node_count)
                        dummy,dummy,u = self.bound_state(E,l,normalize=True)
                        wavefunctions.append(u)
                        break
                    if abs(a-b) < eps_conv:
                        import time
                        print "No eigenvalue, despite sign change???"
                        #time.sleep(10)
                        break
            sys.stdout.flush()
            sys.stderr.flush()
        return energy_scan, penalty, spectrum, nodes, wavefunctions
    def discrete_spectrum(self, energy_range, lmax, nmax=100, eps_conv=1.0e-8):
        """
        Find all eigen states within energy_range with angular momentum between l=0 and l=lmax.

        Parameters:
        ===========
        energy_range: numpy array of energies, each interval
           is searched for eigen energies
        lmax: max. angular momentum
        nmax: check that the first nmax solutions for each angular
           momentum have the correct node counts. 

        Returns: lists of energies and wave functions sorted by angular momenta
        ========
        spectrum: 2D list, spectrum[l][i] contains eigen energy
           of i-th wave function with angular momentum l
        wavefuncs: 2D list, wavefunc[l][i] contains i-th wave function
           with angular momentum l, u_l(r)
        """
        spectrum = []
        wavefuncs = []
        nodes = []
        for l in range(0, lmax+1):
            print "******** l = %s ***********" % l
            energy_scan_l, penalty_l, spectrum_l, nodes_l, wavefuncs_l = self.discrete_spectrum_angmom(l, energy_range, eps_conv=eps_conv)
            spectrum.append(spectrum_l)
            wavefuncs.append(wavefuncs_l)
            nodes.append(nodes_l)
        # check consistency
        for l in range(0, lmax+1):
            for nml,nrnode in enumerate(nodes[l]):
                n = nml + l
                if n >= nmax:
                    # Shooting method can generate unoccupied orbitals, but
                    # we do not care if they have correct node count
                    continue
                if nrnode != nml:
                    msg = "Node counts incorrect. n-th wave function with angular momentum l must have (n-1)-l nodes.\n"
                    msg +="\ngot nodes = %s" % nodes
                    raise Exception(msg)
        return spectrum, wavefuncs
    def tot_xsec(self, E, lmax):
        """
        compute total cross section for scattering of an electron
        with energy E=k^2/2 from this potential into scattered waves with 
        angular momentum up to lmax.
        """
        xsec = 0.0
        k = sqrt(2.0*E)
        for l in range(0, lmax+1):
            delta_l, ul = self.scattering_state(E, l)
            xsec += (2*l+1)*pow(sin(delta_l), 2)
        return xsec#*4.0*pi/pow(k,2)

                
class LennardJones(RadialPotential):
    """
    V(r) = Vm * ( (rm/r)^12 - 2 (rm/r)^6 )
    """
    def __init__(self, Vm, rm, eps=1e-7):
        self.eps = eps # remove singularity at r = 0
        self.Vm = Vm
        self.rm = rm
    def V(self, r):
        return self.Vm*( pow(self.rm/(r+self.eps), 12) - 2.0*pow(self.rm/(r+self.eps), 6) )
    def zero_limit(self, r, E, l):
        C = 1.0/5.0 * sqrt(2.0*self.Vm) * pow(self.rm,6)
        return exp(-C/pow(r+self.eps, 5))
    def asymptotic_decaying(self, r, E, l):
        return exp(-sqrt(2.0*abs(E))*r)
    def asymptotic_oscillating(self, r, E, l):
        k = sqrt(2.0*E)
        solution1 = jn(l,k*r)
        solution2 = yn(l,k*r)
        return solution1, solution2

class Coulomb(RadialPotential):
    """
    single electron in a Coulomb potential generated by a nucleus with positive
    charge Z

    V(r) = -Z/r
    """
    def __init__(self, Z, eps=1e-7):
        self.eps = eps # remove singularity at r=0  1/r -> 1/(r+eps)
        self.Z = Z
        self.nuc_charge = nuclear_charge(self.Z)
        # In the presence of electrons, the nuclear charge can be screened, so that
        # asymptotically only the monopole term is felt
        self.unscreened_charge = self.nuc_charge
    def V(self, r):
        return -self.nuc_charge/(r+self.eps)
    def zero_limit(self, r, E, l):
        return pow(r,l+1)
    def asymptotic_decaying(self, r, E, l):
        if hasattr(self, "has_confinement") and self.has_confinement \
                and hasattr(self, "r0") and 1.0/self.r0 != 0.0:
            """
            If some derived class defines a parameter self.r0, add a confinement potential
            Vconf(r) = (r/r0)^2 to the radial potential. Asymptotically the wave function
            become the eigenstates of the 3D harmonic oscillator with frequency 1/r0^2.
            """
            return pow(r/self.r0,l)*exp(-0.5*pow(r/self.r0, 2))
        else:
            return exp(-sqrt(2.0*abs(E))*r)
    def asymptotic_oscillating(self, r, E, l):
        k = sqrt(2.0*E)
        solution1 = float(mpmath.coulombf(l, -self.unscreened_charge/k, r*k))/r
        solution2 = float(mpmath.coulombg(l, -self.unscreened_charge/k, r*k))/r
        return solution1, solution2
    
class CoulombExpGrid(Coulomb):
    """
    solve radial Schroedinger equation on an exponential grid
    rho = log(r)
    """
    def getVolumeElement(self):
        rho = self.getRadialGrid()
        # dr^3 = 4*pi r^2 dr = 4*pi exp(rho)^2 exp(rho) drho
        return 4.0*pi*exp(3*rho)
    def zero_limit(self, rho, E, l):
        return exp(-rho/2.0) * super(CoulombExpGrid, self).zero_limit(exp(rho), E, l)
    def asymptotic_decaying(self, rho, E, l):
        return exp(-rho/2.0) * super(CoulombExpGrid, self).asymptotic_decaying(exp(rho), E, l)
    def V(self, rho):
        return super(CoulombExpGrid, self).V(exp(rho))
    def getAngularPotential(self, l):
        return 0.5 * l*(l+1)/pow(exp(self.r)+self.eps,2)
    def F(self, rho, E, l):
        """
        using rho=log(r) as radial coordinate and solving for v_l(rho) = r^{-1/2} u_l(r)
        leads to the transformed SE:
          d^2/drho^2 v_l(r) = [ 1/4 + l*(l+1) - 2*exp(2*rho)*(E - V(exp(rho))) ] v_l(rho)
                            = F(rho,E,l) v_l(rho)
        """
        f = 0.25 + l*(l+1) - 2.0*exp(2*rho)*(E-self.Vr)
        return f
    def bound_state(self, E, l, **opts):
        """
        wrap around RadialPotential.bound_state to
        transform wave functions v_l(rho) -> u_l(r)
        u_l(r) = exp(1/2 rho)*v_l(rho)
        """
        if not (opts.get("rmatch") is None):
            # convert matching point to logarithmic scale
            opts["rmatch"] = log(opts["rmatch"])
        ret = super(CoulombExpGrid, self).bound_state(E, l, **opts)
        if type(ret) == type(None) or type(ret[2]) == type(None):
            # no bound state found or the wave function wasn't transformed back
            return ret
        gE, node_count, v = ret
        u = exp(0.5*self.r)*v
        """Renormalize since normalization bound_state 
        is not correct for an exponential grid.
                ul^2(r) dr = ul^2(exp(rho)) * exp(rho) drho """
        # on an exponential grid the unnormalized wave function can be very large
        # (1.0E200). To avoid nans in the normalization, rescale it before squaring
        umax = abs(u.max())
        u_rescaled = u/umax
        norm = sqrt(sum(u_rescaled**2 * exp(self.r) * self.h))
        u = u_rescaled / norm
        """
        print "v.max() = %s" % v.max()
        print "u.max() = %s" % u.max()
        print "u_rescaled.max() = %s" % u_rescaled.max()
        print "norm = %s" % norm
        print "u_normalized.max() = %s" % u.max()
        """
        return gE, node_count, u
    def match_scattering_solution(self, E, l, r1, r2, ul1, ul2):
        ul1 *= exp(0.5*r1)
        ul2 *= exp(0.5*r2)
        delta_l, scale_fac =  super(CoulombExpGrid, self).match_scattering_solution(E,l,exp(r1),exp(r2),ul1,ul2)
        return delta_l, scale_fac*exp(0.5*self.r)

######### TESTING ########################################################

def test_spectrum_Coulomb():
    Coulpot = Coulomb(1.0)
    h,rmin,rmax = 0.03, 0.0, 30.0
    Npts = round(rmax/h)
    r = linspace(rmin, rmax, Npts)

    Coulpot.setRadialGrid(r)
    Coulpot.getRadialPotential()
    energy_range = linspace(-2.0, -0.01, 200)
    lmax = 3

    spectrum, wavefunctions = Coulpot.discrete_spectrum(energy_range, lmax)
    print "Eigenenergies: \n%s" % spectrum
    for l in range(0, lmax+1):
        for n, u in enumerate(wavefunctions[l]):
            plot(Coulpot.getRadialGrid(), u, label="n=%s l=%s" % (n,l))
    legend()
    savefig("test_spectrum_Coulomb.png")

def test_spectrum_Coulomb_on_expgrid():
    Coulpot = CoulombExpGrid(1.0)
    rmin = 0.003
    rmax = 30.0
    rho = linspace(log(rmin), log(rmax), 1000)

    Coulpot.setRadialGrid(rho)
    Coulpot.getRadialPotential()
    energy_range = linspace(-2.0, -0.001, 200)
    lmax = 3

    spectrum, wavefunctions = Coulpot.discrete_spectrum(energy_range, lmax)
    print "Eigenenergies: \n%s" % spectrum
    for l in range(0, lmax+1):
        for n, u in enumerate(wavefunctions[l]):
            plot(exp(Coulpot.getRadialGrid()), u, ls="-.", lw=3, label="n=%s l=%s" % (n,l))
    legend()
    savefig("test_spectrum_Coulomb_on_expgrid.png")


def test_scattering_Coulomb(l, E):
    from CoulombWave import CoulombWave
    Coulpot = Coulomb(1.0)
    h,rmin,rmax = 0.03, 0.0, 30.0
    Npts = round(rmax/h)
    r = linspace(rmin, rmax, Npts)
    Coulpot.setRadialGrid(r)
    Coulpot.getRadialPotential()
    delta_l, u = Coulpot.scattering_state(E, l)
    r = Coulpot.getRadialGrid()
    ion()
    plot(r, u, label="E=%s a.u., l=%s" % (E,l))
    k = sqrt(2.0*E)
    plot(r, CoulombWave(r, 1.0, k, l)*k*r, ls="-.", label="Coulomb function")
    legend()
    draw()
    savefig("scatteringCoulomb.png")

def test_scattering_Coulomb_on_expgrid(l, E):
    from CoulombWave import CoulombWave
    Coulpot = CoulombExpGrid(1.0)
    h,rmin,rmax = 0.03, 0.0001, 30.0
    Npts = round(rmax/h)
    rho = linspace(log(rmin), log(rmax), Npts)
    Coulpot.setRadialGrid(rho)
    Coulpot.getRadialPotential()
    delta_l, u = Coulpot.scattering_state(E, l)
    rho = Coulpot.getRadialGrid()
    r = exp(rho)
    ion()
    plot(r, u, ls="-.", lw=2, label="E=%s a.u., l=%s" % (E,l))
    k = sqrt(2.0*E)
    plot(r, CoulombWave(r, 1.0, k, l)*k*r, ls="-.", label="Coulomb function")
    legend()
    draw()
    savefig("scatteringCoulomb_on_expgrid.png")


def test_phase_shifts_Coulomb():
    Coulpot = Coulomb(1.0)
    h,rmin,rmax = 0.03, 0.0, 300.0
    Npts = round(rmax/h)
    r = linspace(rmin, rmax, Npts)
    Coulpot.setRadialGrid(r)

    Coulpot.getRadialPotential()
    energy_range = linspace(0.01, 500.0, 300)

    for l in range(0, 3):
        phases = []
        for E in energy_range:
            delta_l, u = Coulpot.scattering_state(E, l)
            phases.append(delta_l)
        plot(energy_range, array(phases), label="l = %s" % l)
    xlabel("E")
    ylabel("$\\delta_l$")
    legend()
    savefig("Coulomb_phase_shifts.png")

if __name__ == "__main__":
    from matplotlib.pyplot import plot, ion, legend, savefig, draw, xlabel, ylabel

    test_spectrum_Coulomb()
    test_spectrum_Coulomb_on_expgrid()
    test_scattering_Coulomb(2, 1.74651)
    test_scattering_Coulomb_on_expgrid(2, 1.74651)

    #test_phase_shifts_Coulomb()

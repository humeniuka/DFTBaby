from numpy import array, linspace, arctan2, sqrt, log, exp, pi, mgrid, sum, reshape, where, zeros
import numpy as np
from scipy import interpolate
from PolarTwoCenterGrid import ptcgrid
import slako_transformations as T
import slako_transformations_dipole as Tdip
from DFTB import AtomicData
from DFTB import utils
from hotbit_format import parseHotbitParameters
import os.path

############# CREATION OF SLATER-KOSTER TABLES ###########################################

def integrands_tau(x, h,     e1, l1,m1,R1_spl, Veff1_spl,r1_0,  e2, l2,m2,R2_spl, Veff2_spl,r2_0):
    """
    Find the integrands for calculating the overlap and hamiltonian matrix elements
    between an orbital psi_1(r1,th2,ph2) = R1(r1) Yreal_(l1,m1)(th1,ph1)
    located at z=+h and an orbital psi_2(r2,th2,ph2) located at z=-h.

    Integrating these functions over z=x[0] from  -inf to inf 
    and rho=[1] from 0 to inf gives the overlap S_12(h) and the
    hamiltonian matrix element H_12(h) as a function of the the 
    distance 2*h between the centers.

    Parameters:
    ===========
    x: tuple of coordinates (z,rho)
    h: 1/2 distance between orbital 1 and orbital 2
    e1,e2: energies of confined pseudo atomic orbitals
    l1,m1: angular quantum numbers of orbital 1
    l2,m2: angular quantum numbers of orbital 2
    R1_spl, R2_spl: callable functions that return
       the radial part of the wavefunctions 1 and 2 at the position
       given as the argument.
    Veff1_spl, Veff2_spl: callable functions that return
       the effective radial potential of atoms 1 and 2 (without confinement and
       without l*(l+1)/(2 r^2) angular potential)
    r1_0,r2_0: confinement radius used in calculating pseudo atoms 1 and 2

    Returns:
    ========
    s,h: 
       - s(z,rho) is the integrand such that S_12(tau,h) = integral dz drho Is(z,rho)
         s(z,rho) = rho * R1(r1)*R2(r2)*phi_tau(th1,th2)
          see eq. 142 in Koskinen, Maekinen 2009
       - h(z,rho) 
         h(z,rho) = rho * R1(r1)*R2(r2)*phi_tau(th1,th2) * ( Veff2(r2) - Vconf1(r1) )
    """
    z = x[0]
    rho = x[1]
    th1 = arctan2(rho, z-h)
    th2 = arctan2(rho, z+h)
    r1 = sqrt(rho**2 + (z-h)**2)
    r2 = sqrt(rho**2 + (z+h)**2)
    tau = (l1,m1,l2,m2)

    """
    # plot wave function
    from matplotlib.pyplot import plot, show, cla
    cla()
    sort_indx = np.argsort(r1)
    plot(r1[sort_indx], R1_spl(r1[sort_indx]), ls="-.", color="green")
    show()
    #
    """
    R1r1 = R1_spl(r1)
    R2r2 = R2_spl(r2)
    s = rho*R1r1*R2r2*T.angular_phi[tau](th1,th2)
    # h = s*((effective potential of atom 2) - (confinement of atom 1))
    # THERE MIGHT BE SOMETHING WRONG HERE !!!!
    h = s*(Veff2_spl(r2) - pow(r1/r1_0,2)) # this differs from formula (149) in Koskinen/Maekinen
                           # [Veff2 - Vconf2 - Vconf1] because in the splined effective potential
                           # the confinement potential is already subtracted out.

    # orbital probabilities, check normalization
    # r^2*sin(th) = rho*r
    p1 = rho*pow(R1r1,2)*T.angular_phi[(l1,m1,l1,m1)](th1,th1)
    p2 = rho*pow(R2r2,2)*T.angular_phi[(l2,m2,l2,m2)](th2,th2)

    """
    if l1 == 0 and l2 == 1:
        from mpl_toolkits.mplot3d import Axes3D
        from matplotlib.pyplot import figure, show
        fig = figure()
        ax = fig.gca(projection='3d')
        print "l1 = %s m1 = %s" % (l1,m1)
        print "max(R1) = %s at rho=%s, z = %s" % (max(abs(R1_spl(r1))), rho[(abs(R1_spl(r1))).argmax()], z[abs(R1_spl(r1)).argmax()])
        print "max(r1) = %s at rho=%s, z = %s" % (max(r1), rho[r1.argmax()], z[r1.argmax()])
        print "max(rho) = %s at rho=%s, z = %s" % (max(rho), rho[rho.argmax()], z[rho.argmax()])
        ax.plot(rho, z, s)
        try:
            show()
        except AttributeError as e:
            print e
    """
    return s,h, p1,p2


def integrands_tau_dipole(x, h, l1,m1,R1_spl, lM,mM, l2,m2,R2_spl):
    """
    Find the integrands for calculating the dipole matrix elements
    between an orbital psi_1(r1,th2,ph2) = R1(r1) Yreal_(l1,m1)(th1,ph1)
    located at z=+h and an orbital psi_2(r2,th2,ph2) located at z=-h.

    Integrating these functions over z=x[0] from  -inf to inf 
    and rho=[1] from 0 to inf gives the dipole vector (Dx,Dy,Dz) a function of the the 
    distance 2*h between the centers.

    Parameters:
    ===========
    x: tuple of coordinates (z,rho)
    h: 1/2 distance between orbital 1 and orbital 2
    l1,m1: angular quantum numbers of orbital 1
    l2,m2: angular quantum numbers of orbital 2
    R1_spl, R2_spl: callable functions that return
       the radial part of the wavefunctions 1 and 2 at the position
       given as the argument.
    lM,mM: angular quantum numbers of multipole operator, (1,1) -> px, (1,-1) -> py, (1,0) -> pz

    Returns:
    ========
    dip  multipole integrand (actually only x,y,z so far)
    """
    z = x[0]
    rho = x[1]
    th1 = arctan2(rho, z-h)
    th2 = arctan2(rho, z+h)
    r1 = sqrt(rho**2 + (z-h)**2)
    r2 = sqrt(rho**2 + (z+h)**2)

    R1r1 = R1_spl(r1)
    R2r2 = R2_spl(r2)
    # The indeces for dipoles are enumerated differently, find the unique ones
    tau = (l1,m1,lM,mM,l2,m2)
    dip = r1 * rho*R1r1*R2r2*Tdip.angular_phi[tau](th1,th2)  

    return dip

def spline_wavefunction(rg,ug,ug_asymptotic=lambda r: 0*r):
    """
    Create a spline of the radial part R(r) from a 
    a discretized version ug of u(r)=r*R(r) on the radial grid rg.

    Parameters:
    ===========
    rg: numpy array with r-points
    ug: numpy array with values of u on rg

    Optional:
    =========
    ug_asymptotic: function u(r) for large r, 
         outside the interpolation range this asymptotic form of the wavefunction
         is used, for bound wavefunctions the function should give 0, while 
         for scattering wavefunctions the oscillating solution with the correct phase shift 
         should be returned.

    Returns:
    ========
    a callable function R(r)
    """
    """
    # plot wave function
    from matplotlib.pyplot import plot, show, cla
    cla()
    plot(rg, ug/rg)
    show()
    #
    """
    tck = interpolate.splrep(rg,ug, s=0)
    rmin, rmax = min(rg), max(rg)

    def R(r):
        rflat = r.ravel()
        """
        # old code
        u = reshape(where((rflat<=rmax) & (rmin<=rflat), \
                          interpolate.splev(rflat,tck,der=0), ug_asymptotic(rflat)), r.shape)
        """
        # interpolate wavefunction in the range [0,rmax]
        u = 0.0*rflat
        u[rflat <= rmax] = interpolate.splev(rflat[rflat <= rmax], tck, der=0)
        # use asymptotic form for rmax < r
        u[rmax < rflat] = ug_asymptotic(rflat[rmax < rflat])

        u = reshape(u, r.shape)
        
        return u/r

    return R

def spline_effective_potential(rg,Veffg,r0):
    """
    Create spline of the effective potential without the confinement potential.

    Parameters:
    ===========
    rg: numpy array with r-points
    Veffg: self-consistent effective potential (including confinement (r/r0)^2 on rg 
    r0: confinement radius

    Returns:
    =======
    callable function with approximation to true crystal potential Veff(r)
    """
    # WARNING: When rg includes the origin (r[0] == 0), the spline
    # will oscillate wildly. 
    tck = interpolate.splrep(rg,Veffg-pow(rg/r0,2), s=0)
    rmax = max(rg)
    def Veff(r):
        rflat = r.ravel()
        # outside rmax, the effective potential should be zero
        return reshape(where(rflat<rmax, interpolate.splev(rflat,tck,der=0), 0.0), r.shape)
    return Veff

def spline_radial_density(rg,rhog):
    """
    Create spline of the total radial electron density

    Parameters:
    ===========
    rg: numpy array with r-points
    rhog: self-consistent radial electron density

    Returns:
    =======
    callable function rho(r)
    """
    tck = interpolate.splrep(rg,rhog, s=0)
    rmax = max(rg)
    def rho(r):
        rflat = r.ravel()
        # outside rmax, the density should be zero
        return reshape(where(rflat<rmax, interpolate.splev(rflat,tck,der=0), 0.0), r.shape)
    return rho

def test_valence_orbitals():
    from confined_pseudo_atoms import h,he,li,be,b,c,n,o,f,ne,na,si

    assert h.valence_orbitals == [0]
    assert he.valence_orbitals == [0]
    assert c.valence_orbitals == [1,3]
    assert n.valence_orbitals == [1,3]
    assert o.valence_orbitals == [1,3]
    # check atoms with higher Z, too !

class Atom:
    def __init__(self, atom_data):
        self.atom = atom_data
    def splineValenceOrbitals(self):
        """
        produce spline representation of the radial parts of the valence orbitals
        """
        self.radial_val = [] # list of splines of the radial parts of the valence orbitals
        self.energies_val = [] # list of energies of the valence orbitals
        for hb_indx,indx in enumerate(self.atom.valence_orbitals):
            R_spl = spline_wavefunction(self.atom.r, self.atom.radial_wavefunctions[indx])
            #
            """
            from matplotlib.pyplot import plot, legend, show
            from hotbit_format import parseHotbitParameters
            hbdata = parseHotbitParameters("DFTB/test_parameters/C.elm")
            indx2name = {1: "2s", 2: "2p"}
            orbname = "radial_wavefunction_%s" % indx2name[hb_indx+1]
            plot(hbdata[orbname][:,0], hbdata[orbname][:,1], ls="-.", label=orbname)

            r = linspace(0.0, 30.0, 1000)
            plot(r,R_spl(r)*r)
            plot(self.atom.r, self.atom.radial_wavefunctions[indx], lw=2, ls="-.", label="indx = %s" % indx)
            legend()
            show()
            """
            #
            self.radial_val.append(R_spl)
            self.energies_val.append(self.atom.energies[indx])
    def splineEffectivePotential(self):
        self.Veff_spl = spline_effective_potential(self.atom.r, self.atom.effective_potential, self.atom.r0)
    def getValenceOrbitals(self):
        """gives dictonary with tuples (<orbital energy>, <spline of orbital>) of the 
        valence orbitals indexed by angular momentum."""
        R_valence = {}
        for i,indx in enumerate(self.atom.valence_orbitals):
            n, l = self.atom.nshell[indx], self.atom.angular_momenta[indx]
            assert not l in R_valence.keys() # assume that all valence orbitals have different l
                                         # e.g. 4s and 5s cannot be valence orbitals at the 
                                         # same time
            # check normalization
            r = linspace(0.00001, 20.0, 300)
            norm = sqrt(sum(pow(self.radial_val[i](r),2)*pow(r,2)*(r[1]-r[0]))*2*pi)
#            print "norm integral r^2*R(r)d3r = %s" % norm

            R_valence[l] = (self.energies_val[i], self.radial_val[i])
        return R_valence
    def getEffectivePotential(self):
        """returns spline of the effective SCF potential without confining potential"""
        return self.Veff_spl

class AtomPair:
    def __init__(self, atom_data1, atom_data2):
        """
        atom_data1 and atom_data2 are precomputed modules in the
        directory pseudo_atoms. 
        e.g. from pseudo_atoms import li as atom_data1
        """
        # Because Slater-Koster files for A-B and B-A contain the
        # same information, sort atoms by atomic number so that 
        # we only need to write one table but still know which atom 
        # is first and which is the second.
        if atom_data1.Z > atom_data2.Z:
            tmp = atom_data1
            atom_data1 = atom_data2
            atom_data2 = tmp
        self.A1 = Atom(atom_data1)
        self.A2 = Atom(atom_data2)
        self.A1.splineValenceOrbitals()
        self.A1.splineEffectivePotential()
        self.A2.splineValenceOrbitals()
        self.A2.splineEffectivePotential()
    def getSKIntegrals(self, d, grid):
        """
        compute Slater-Koster integrals for overlaps S and hamiltonian
        matrix elements H between valence orbitals as functions of the
        distance d between the atomic centers.

        Parameters:
        ===========
        d: 1D numpy array with distances
        grid: tuple (rhos,zs,areas) as returned by ptcgrid(...)
           -rhos[i] contains a 1D numpy array with rho-coordinates
            of the polar grid where the two centers lie d[i] apart.
           -zs[i] contains a 1D numpy array with z-coordinates
           -areas[i] contains the areas associated with the grid points

        Returns:
        ========
        S,H: S[(l1,l2,i)], H[(l1,l2,i)] contain the tau(i) integrals for
           matrix elements between the valence orbital with angular momentum l1
           on atom 1 and the valence orbital with angular momentum l2 on atom 2.
        """
        R_valence1 = self.A1.getValenceOrbitals()
        R_valence2 = self.A2.getValenceOrbitals()
        Veff1_spl = self.A1.getEffectivePotential()
        Veff2_spl = self.A2.getEffectivePotential()

        self.S = {}
        self.H = {}
        # distance
        self.d = d
        # overlaps and hamiltonian matrix elements
        for i in T.index2tau.keys():
            l1,m1,l2,m2 = T.index2tau[i]
            if not (R_valence1.has_key(l1) and R_valence2.has_key(l2)):
                continue
            en1, R1_spl = R_valence1[l1]
            en2, R2_spl = R_valence2[l2]
            """
            # integration on cartesian grid
            # dmin, dmax and Nd are global variables
            d, zs,rhos = mgrid[dmin:dmax:Nd*1j,-25.0:25.0:500j,0.0:10.0:300j] 
            s,h, p1,p2 = integrands_tau((zs,rhos),d/2, \
                        en1,l1,m1,R1_spl,Veff1_spl,self.A1.atom.r0, \
                        en2,l2,m2,R2_spl,Veff2_spl,self.A2.atom.r0)
            area = (zs[0,1,0]-zs[0,0,0])*(rhos[0,0,1]-rhos[0,0,0])
            olap = sum(sum(s*area, axis=2), axis=1)
            Hpart = sum(sum(h*area, axis=2), axis=1)
            #
            """
            # integration on two center polar grid
            olap = []
            Hpart = []
            for k,dk in enumerate(self.d):
                rhos,zs,areas = grid[0][k], grid[1][k], grid[2][k]
                #
                """
                if dk > 3.0:
                    from matplotlib.pyplot import plot,show
                    plot(rhos, zs, "o")
                    show()
                """
                #
                s,h,p1,p2 = integrands_tau((zs,rhos),dk/2.0, \
                        en1,l1,m1,R1_spl,Veff1_spl,self.A1.atom.r0, \
                        en2,l2,m2,R2_spl,Veff2_spl,self.A2.atom.r0)
                norm1 = sum(p1*areas)
                norm2 = sum(p2*areas)
#                print "norm1 = %s" % norm1
#                print "norm2 = %s" % norm2
                olap.append( sum(s*areas) )
                Hpart.append( sum(h*areas) )
            olap = array(olap)
            Hpart = array(Hpart)
            self.S[(l1,l2,i)] = olap
#            self.H[(l1,l2,i)] = en2*olap + Hpart
            self.H[(l1,l2,i)] = en1*olap + Hpart
        # dipoles
        self.Dipole = {}

        for i in Tdip.index2tau.keys():
            l1,m1,lM,mM,l2,m2 = Tdip.index2tau[i]
            if not (R_valence1.has_key(l1) and R_valence2.has_key(l2)):
                continue
            en1, R1_spl = R_valence1[l1]
            en2, R2_spl = R_valence2[l2]
       
            # integration on two center polar grid
            Dippart = []
            for k,dk in enumerate(self.d):
                rhos,zs,areas = grid[0][k], grid[1][k], grid[2][k]
                #
                """
                if dk > 3.0:
                    from matplotlib.pyplot import plot,show
                    plot(rhos, zs, "o")
                    show()
                """
                #
                dip = integrands_tau_dipole((zs,rhos),dk/2.0, \
                        l1,m1,R1_spl, lM,mM, l2,m2,R2_spl)
                Dippart.append( sum(dip*areas) )
            Dippart = array(Dippart)
            self.Dipole[(l1,l2,i)] = Dippart

        return self.S,self.H,(self.Dipole)
    def plotSKIntegrals(self, img_dir):
        from matplotlib.pyplot import cla, xlabel, ylabel, legend, plot, savefig, title, show
    
        atname1 = AtomicData.atom_names[self.A1.atom.Z-1]
        atname2 = AtomicData.atom_names[self.A2.atom.Z-1]
        imgfileS = os.path.join(img_dir, "overlap_%s_%s.png" % (atname1, atname2))
        imgfileH = os.path.join(img_dir, "hamiltonian_%s_%s.png" % (atname1, atname2))
        imgfileDip = os.path.join(img_dir, "dipoles_%s_%s.png" % (atname1,atname2))

        cla()
        xlabel("distance $d$ between centers / bohr")
        ylabel("overlap")
        title("Overlaps between $Z_1=%s$ and $Z_2=%s$" % (self.A1.atom.Z,self.A2.atom.Z))
        for (l1,l2,i),olap in self.S.iteritems():
            plot(self.d, olap, label="$S_{%s,%s}(%s)(d)$" % \
                     (l1,l2,T.tau2symbol[T.index2tau[i]]), lw=2)
        legend()
        savefig(imgfileS)
#        show()

        cla()
        xlabel("distance $d$ between centers / bohr")
        ylabel("H integrals")
        title("Hamiltonian integrals between $Z_1=%s$ and $Z_2=%s$" % (self.A1.atom.Z,self.A2.atom.Z))
        for (l1,l2,i),H in self.H.iteritems():
            plot(self.d, H, label="$H_{%s,%s}(%s)(d)$" % \
                     (l1,l2,T.tau2symbol[T.index2tau[i]]), lw=2)
        legend()
        savefig(imgfileH)
#        show()

        cla()
        xlabel("distance $d$ between centers / bohr")
        ylabel("Dipole integrals")
        title("Dipole integrals between $Z_1=%s$ and $Z_2=%s$" % (self.A1.atom.Z,self.A2.atom.Z))
        angmom_to_xyz = {(0,0): "s", (1,-1): "px", (1,1): "py", (1,0): "pz"}
        for (l1,l2,i),D in self.Dipole.iteritems():
            lo1,mo1,lM,mM,lo2,mo2 = Tdip.index2tau[i]
            plot(self.d, D, label="$Dipole_{%s,%s}(%s)(d)$" % \
                     (angmom_to_xyz[(lo1,mo1)],angmom_to_xyz[(lo2,mo2)],angmom_to_xyz[(lM,mM)]), lw=2)
        legend()
        savefig(imgfileDip)
#        show()                      
    def write_py_slako(self, slako_dir):
        """
        write Slater-Koster tables to a python module that can be loaded.

        Parameters:
        ===========
        slako_dir: directory where the Slator Koster modules are located.
           Dictionaries for overlaps and hamiltonian matrix elements
           for the atom pair will be stored in a  module file called 
           <slako_dir>/<atom 1>_<atom 2>.py.
        """
        atname1 = AtomicData.atom_names[self.A1.atom.Z-1]
        atname2 = AtomicData.atom_names[self.A2.atom.Z-1]
        slako_file = os.path.join(slako_dir, "%s_%s.py" % (atname1, atname2))

        import pprint
        import sys
        np.set_printoptions(threshold=sys.maxint)
        pp = pprint.PrettyPrinter(depth=10)
        fh = open(slako_file, "w")
        print>>fh, "# This file has been generated automatically by %s" % sys.argv[0]
        print>>fh, "# from %s and %s." % (self.A1.atom.__file__, self.A2.atom.__file__)
        print>>fh, "from numpy import array"
        print>>fh, "Z1 = %s" % self.A1.atom.Z
        print>>fh, "Z2 = %s" % self.A2.atom.Z
        print>>fh, "# overlaps S[(l1,l2,i)] and hamilton matrix elements H[(l1,l2,i)]"
        print>>fh, "# l1 and l2 are the angular quantum numbers of valence orbitals"
        print>>fh, "# on atom1 and atom2 respectively."
        print>>fh, "# i enumerates the Slater-Koster integrals:"
        print>>fh, "index2symbol = \\\n%s" % dict([(i, T.tau2symbol[T.index2tau[i]]) for (l1,l2,i) in self.S.keys()])
        print>>fh, "# grid for distance d between atomic centers"
        print>>fh, "d = \\\n%s" % pp.pformat(self.d)
        print>>fh, "# overlaps"
        print>>fh, "S = \\\n%s" % pp.pformat(self.S)
        print>>fh, "# hamiltonian matrix elements"
        print>>fh, "H = \\\n%s" % pp.pformat(self.H)
        print>>fh, "# dipoles"
        print>>fh, "Dipole = \\\n%s" % pp.pformat(self.Dipole)
        if hasattr(self, "PhaseShifts"):
            print>>fh, "# phase shifts"
            print>>fh, "PhaseShifts = \\\n%s" % pp.pformat(self.PhaseShifts)
        fh.close()


################## READING SLATER-KOSTER FILES IN DIFFERENT FORMATS ###################

def process_slako_line(l):
    """
    convert a line into a list of column values respecting the
    strange format conventions used in DFTB+ Slater-Koster files.

    In Slater-Koster files used by DFTB+ zero columns 
    are not written: e.g. 4*0.0 has to be replaced
    by four columns with zeros 0.0 0.0 0.0 0.0.
    """
    parts = l.replace(",", " ").strip().split()
    columns = []
    for p in parts:
        if "*" in p:
            f,v = p.split("*")
            for c in range(0, int(f)):
                columns.append(float(v))
        else:
            columns.append(float(p))
    return columns

class SlakoModule:
    pass

def read_slakoformat_skf(skfile, sk=None, atom_order="AB"):
    """
    read table for S and H in the format used by DFTB+ 
    see http://www.dftb.org/fileadmin/DFTB/public/misc/slakoformat.pdf

    Parameters:
    ===========
    skfile: path to <atom 1>-<atom 2>.skf
    sk: slako_module loaded from a different file, to which the data from current table
      is appended.
    atom_order: Which of the atoms in the SK table is to be treated as the first atom?
      "AB": The first orbital belongs to the first atom
      "BA": The first orbital belongs to the second atom

    
    Returns:
    ========
    sk

    The Slater Koster data is stored in attributes of sk:
    sk.Z1, sk.Z2, sk.d, sk.S, sk.H    
    """
    from os.path import basename
    at1, at2 = basename(skfile).replace(".skf", "").split("-")
    fh = open(skfile, "r")
    # line 1
    parts =  process_slako_line(fh.readline())
    gridDist, nGridPoints = float(parts[0]), int(parts[1])
    if at1 == at2:
        # line 2
        # homonuclear case
        Ed, Ep, Es, SPE, Ud, Up, Us, fd, fp, fs = map(float, process_slako_line(fh.readline()))
    # line 2 or 3
    mass, c2, c3, c4, c5, c6, c7, c8, c9, rcut, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10 = \
        map(float, process_slako_line(fh.readline()))

    d = linspace(0.0, gridDist*(nGridPoints-1), nGridPoints)
    if sk == None:
        sk = SlakoModule()
        sk.Z1 = AtomicData.atomic_number(at1)
        sk.Z2 = AtomicData.atomic_number(at2)
        sk.d = d
        sk.S = {}
        sk.H = {}
    S = sk.S
    H = sk.H

    if atom_order == "BA":
        tausymbols = T.tausymbols_BA
    else:
        tausymbols = T.tausymbols_AB

    lines = fh.readlines()
    parts = process_slako_line(lines[0])
    n = len(parts)/2
    assert n == 10 # for orbitals up to d functions there are 10 Slater Koster integrals
    for tausym in tausymbols[-n:]:
        try:
            tau = T.symbol2tau[tausym]
        except KeyError:
            continue
#        print "found %s integrals" % (T.tau2symbol[tau])
        H[(tau[0], tau[2], T.tau2index[tau])] = array([0.0 for i in range(0, int(nGridPoints)) ])
        S[(tau[0], tau[2], T.tau2index[tau])] = array([0.0 for i in range(0, int(nGridPoints)) ])
    # line 4 to (4 + nGridPoints -1)
    for i in range(0, nGridPoints):
        parts = process_slako_line(lines[i])
#Hdd0 Hdd1 Hdd2 Hpd0 Hpd1 Hpp0 Hpp1 Hsd0 Hsp0 Hss0 Sdd0 Sdd1 Sdd2 Spd0 Spd1 Spp0 Spp1 Ssd0 Ssp0 Sss0
        for pos, tausym in enumerate(tausymbols[-n:]):
            try:
                tau = T.symbol2tau[tausym]
            except KeyError:
                continue
            l1,l2 = tau[0], tau[2]
            if atom_order == "BA": 
                orbital_parity = pow(-1,l1+l2)
            else:
                orbital_parity = 1                

            H[(l1, l2, T.tau2index[tau])][i] = orbital_parity * float(parts[pos])
            S[(l1, l2, T.tau2index[tau])][i] = orbital_parity * float(parts[len(tausymbols)+pos])

    return sk

def read_slakoformat_par(parfile, atom_order="AB"):
    """
    read table for S and H in the format used by hotbit (https://trac.cc.jyu.fi/projects/hotbit)

    Paramters:
    ==========
    parfile: path to <atom1>_<atom2>.par

    Returns:
    ========
    sk

    The Slater Koster data is stored in attributes of sk:
    sk.Z1, sk.Z2, sk.d, sk.S, sk.H    
    """

    from os.path import basename
    at1, at2 = basename(parfile).replace(".par", "").split("_")
    if atom_order == "BA":
        tmp = at1
        at1 = at2
        at2 = tmp
    data_blocks = parseHotbitParameters(parfile)
    data_AB = data_blocks["slako_integrals_%s_%s" % (at1,at2)]
    data_BA = data_blocks["slako_integrals_%s_%s" % (at2,at1)]
    m,n = data_AB.shape
    assert data_AB.shape == data_BA.shape

    d = data_AB[:,0]
    assert np.all(data_AB[:,0] == data_BA[:,0])

    sk = SlakoModule()
    sk.Z1 = AtomicData.atomic_number(at1)
    sk.Z2 = AtomicData.atomic_number(at2)
    sk.d = d
    sk.S = {}
    sk.H = {}

    S = sk.S
    H = sk.H
    for pos, (tausym_AB, tausym_BA) in enumerate(zip(T.tausymbols_AB[-n:], T.tausymbols_BA[-n:])):
        try:
            tau_AB = T.symbol2tau[tausym_AB]
            tau_BA = T.symbol2tau[tausym_BA]
        except KeyError:
            continue
        # AB
        l1,l2 = tau_AB[0], tau_AB[2]

        H[(l1, l2, T.tau2index[tau_AB])] = data_AB[:,pos+1]
        S[(l1, l2, T.tau2index[tau_AB])] = data_AB[:,len(T.tausymbols_AB)+pos+1]
        # BA
        l1,l2 = tau_BA[0], tau_BA[2]

##        if sk.Z1 == sk.Z2: # I think this is the right condition
        if sk.Z1 > 1 and sk.Z2 > 1: # but with this illogical condition I can reproduce Roland's results for HCNO compounds
            orbital_parity = pow(-1,l1+l2)
        else:
            orbital_parity = 1
        H[(l1, l2, T.tau2index[tau_BA])] = orbital_parity * data_BA[:,pos+1]
        S[(l1, l2, T.tau2index[tau_BA])] = orbital_parity * data_BA[:,len(T.tausymbols_BA)+pos+1]

    return sk

def read_slakoformat(filename):
    """
    read Slater Koster tables from different file formats, which
    are recognized by the file extension.
    
    The following extensions are supported:
    *.skf    used by DFTB+
    *.par    used by hotbit
    *.py     

    Paramters:
    ==========
    filename: path to Slater Koster file

    Returns:
    ========
    sk

    The Slater Koster data is stored in attributes of sk:
    sk.Z1, sk.Z2, sk.d, sk.S, sk.H    
    """
    from os.path import exists, basename, dirname, join

    if ".skf" in filename:
        filename_AB = filename
        # load Slater Koster integrals for A atom
        sk = read_slakoformat_skf(filename_AB, sk=None, atom_order="AB")
        at1, at2 = basename(filename).replace(".skf", "").split("-")
        # append Slater Koster integrals for B atom
        filename_BA = join(dirname(filename), "%s-%s.skf" % (at2,at1))
        sk = read_slakoformat_skf(filename_BA, sk=sk, atom_order="BA")
        return sk
    elif ".par" in filename:
        atom_order = "AB"
        if not exists(filename):
            # file file A_B.par does not exist, look for file B_A.par
            at1, at2 = basename(filename).replace(".par", "").split("_")
            filename_BA = join(dirname(filename), "%s_%s.par" % (at2,at1))
            print "parameter file %s does not exist -> try to load %s" % (filename, filename_BA)
            filename = filename_BA
            atom_order = "BA"
            assert exists(filename)
        sk = read_slakoformat_par(filename, atom_order=atom_order)
        return sk
    elif ".py" in filename:
        # access dictionary by .-notation
        mod = utils.dotdic()
        execfile(filename, mod)
#        mod = weird_sign_change(mod)
        return mod

#def weird_sign_change(slako_mod):
#    """
#    change sign of all integrals for atom pairs H-X (X is not a hydrogen).
#    This does not have an impact on final energies or gradients 
#    but makes it easier to compare with hotbit Slater Koster integrals.
#    """
#    if slako_mod.Z1 == 1 and slako_mod.Z2 != 1:
#        print "WEIRD SIGN CHANGE"
#        for k,v in slako_mod.S.iteritems():
#            slako_mod.S[k] = -v
#        for k,v in slako_mod.H.iteritems():
#            slako_mod.H[k] = -v
#    return slako_mod
        

##################################################################

def test_hh():
    from matplotlib.pyplot import plot, show, savefig
    from confined_pseudo_atoms import h,li
    test_valence_orbitals()
    # d: distance between centers
    # zs,rhos: integration variables
    R1_spl = spline_wavefunction(h.r, h.radial_wavefunctions[0])
    R2_spl = spline_wavefunction(li.r, li.radial_wavefunctions[1])
    I = overlap_tau((zs,rhos), d/2, 0,0, R1_spl, 0,0, R2_spl)

    area = (zs[0,1,0]-zs[0,0,0])*(rhos[0,0,1]-rhos[0,0,0])
    S_12 = sum(sum(I*area, axis=2), axis=1)
    print len(d[:,0,0])
    print len(S_12)
    
    plot(d[:,0,0],R1_spl(d[:,0,0]))
    plot(d[:,0,0],S_12)
    savefig("h_h_S_12.png")

def test_cc():
    from confined_pseudo_atoms import c

    script_dir = os.path.dirname(os.path.realpath(__file__))
    slako_dir = os.path.join(script_dir, "slako_tables/")

    dimer = AtomPair(c,c)

    dmin, dmax, Nd = 0.0, 10.0, 2#150
    rmin, rmax, Nr = 0.00001, 25.0,1000
    Na = 150
    angles = linspace(0.0, pi, Na)
    r = linspace(max(1.0e-10, rmin), rmax, Nr)#exp(linspace(log(rmin), log(rmax), Nr))
    d = linspace(dmin, dmax, Nd)
    grid = ptcgrid(d/2.0,r,angles)

    dimer.getSKIntegrals(d, grid)
    dimer.write_py_slako(slako_dir)
    dimer.plotSKIntegrals(slako_dir)
    

########################## USAGE OF SLATER-KOSTER TABLES #################

def spline_slako_integrals(dg, SorHg):
    """
    Create a spline of a Slator Koster integral as a function of the
    separation between the two atomic centers.

    Parameters:
    ===========
    dg: numpy array with distances
    SorHg: hash with numpy arrays containing integrals (overlaps S or Hamiltonian matrix elements H)
      on the grid dg.

    Returns:
    ========
    - a callable function S(i,d,der=0) or H(i,d,der=0) where i refers to the i-th Slater-Koster integral and d
      is the distance between the atoms. With der=1,2,.. the 1st, 2nd etc derivative of S or H is returned.
    - a dictionary with tuples containing the spline data, (knots, coefs, degree), for all Slater-Koster integrals.
    """
    dmax = max(dg)
    splines_tck = {}
    for (l1,l2,i),Ig in SorHg.iteritems():
        splines_tck[i] = interpolate.splrep(dg,Ig, s=None) #s=None)
        # Without smoothing (s=0) numerical and analytical gradients can deviate
        # for molecules such as formonitrile because of kinks in the curves
        # But with too much smoothing (s=1) we get the wrong band structure energy for benzene, because
        # the spline deviates more from the original data points.
        # When using analytical gradients smoothing does not matter and should be low to ensure the correct energy.
    def SorH(i,d,der=0):
        dflat = d.ravel()
        # set to zero outside [0.0,dmax] instead of interpolating which gives nonsense
        interp = interpolate.splev(dflat,splines_tck[i],der=der)
        SorH_interp = reshape(where((dflat<=dmax) & (0.0<=dflat), \
                          interp, 0.0), d.shape)
        return SorH_interp

    return SorH, splines_tck

class NoSKDipolesException(Exception):
    pass

class SlakoTransformations:
    """
    Slater Koster table with transformation for an atom pair
    """
    def __init__(self, slako_module):
        """
        slako_module: python module which was loaded from slako_tables
        """
        self.slako = slako_module
        self.S_spl, self.S_tck = spline_slako_integrals(self.slako.d, self.slako.S)
        self.H_spl, self.H_tck = spline_slako_integrals(self.slako.d, self.slako.H)
        if self.haveDipoles():
            self.Dipole_spl, self.Dipole_tck = spline_slako_integrals(self.slako.d, self.slako.Dipole)
    def _directional_cosines(self, pos1,pos2):
        """
        compute directional cosines for the vector going from
        pos1 to pos2

        Returns:
        ========
        r,x,y,z   r length of vector, x,y,z directional cosines
        """
        xc = pos2[0]-pos1[0]
        yc = pos2[1]-pos1[1]
        zc = pos2[2]-pos1[2]
        r = sqrt(xc**2+yc**2+zc**2)
        # directional cosines
        if r > 0.0:
            x,y,z = xc/r,yc/r,zc/r
        else:
            x,y,z = 0,0,1
        return r,x,y,z
    def getOverlap(self, l1,m1,pos1,  l2,m2,pos2, deriv=0):
        """
        compute overlap (and its derivatives with respect to the r=pos2-pos1) 
        between valence orbital on atom1 and valence orbital on atom2 
        using tabulated Slater Koster integrals and transformation rules

        Paramters:
        ==========
        l1,m1: angular momentum and magnetic qnumber of valence orbital on atom 1
        l2,m2: angular momentum and magnetic qnumber of valence orbital on atom 2
        pos1: (x1,y1,z1) cartesian coordinates of atom 1
        pos2: (x2,y2,z2) cartesian coordinates of atom 2
        deriv: 0 -> overlap
               1 -> gradient of overlap
               2 -> hessian of overlap

        Returns:
        ========
        overlap S or gradient of S or hessian of S (depending on deriv) 
        """
        r,x,y,z = self._directional_cosines(pos1,pos2)
        if deriv == 0:
            S = T.slako_transformations[(l1,m1,l2,m2)](r,x,y,z, self.S_spl)
            """
            # test new Fortran code
            indeces = [i for (l1dummy,l2dummy,i) in self.slako.S.keys()]
            S_spl_arr = np.zeros(max(indeces)+1)
            for i in indeces:
                S_spl_arr[i] = self.S_spl(i,r)
            #print "S_spl_arr = %s" % S_spl_arr
            from DFTB.extensions import slako
            S_f90 = slako.slako.slako_transformation(r,x,y,z, S_spl_arr, l1,m1,l2,m2)
            # compare with python code
            assert abs(S-S_f90) < 1.0e-8, "S_py = %s   S_f90 = %s for (%s,%s,%s,%s)" % (S, S_f90, l1,m1,l2,m2)
            #
            """
            return S
        elif deriv == 1:
            gradS = list(zeros(3))
            for xyz in range(0,3):
                gradS[xyz] = T.slako_transformations_gradient[xyz][(l1,m1,l2,m2)](r,x,y,z, self.S_spl)
            """
            # test new Fortran code
            indeces = [i for (l1dummy,l2dummy,i) in self.slako.S.keys()]
            S_spl_arr = np.zeros(max(indeces)+1)   # SK integrals for S
            Sd_spl_arr = np.zeros(max(indeces)+1)  # derivatives of SK integrals for S
            for i in indeces:
                S_spl_arr[i] = self.S_spl(i,r)
                Sd_spl_arr[i] = self.S_spl(i,r,der=1)
            from DFTB.extensions import slako
            gradS_f90 = slako.slako.slako_transformation_gradient(r,x,y,z, S_spl_arr, Sd_spl_arr, l1,m1,l2,m2)
            # compare with python code
            assert np.sum(abs(gradS-gradS_f90)) < 1.0e-8, "grad S_py = %s   grad S_f90 = %s for (%s,%s,%s,%s)" % (gradS, gradS_f90, l1,m1,l2,m2)
            #
            """
            return gradS
        elif deriv == 2:
            hessS = list(zeros((3,3)))
            for xyz1 in range(0,3):
                for xyz2 in range(0,3):
                    #print xyz1,xyz2
                    hessS[xyz1][xyz2] = T.slako_transformations_hessian[xyz1][xyz2][(l1,m1,l2,m2)](r,x,y,z, self.S_spl)
            return hessS
        else:
            raise Exception("Higher than 2nd derivatives of S or H not implemented")
    def getHamiltonian0(self, l1,m1,pos1, l2,m2,pos2, deriv=0):
        """
        compute Hamiltonain matrix element (and its derivatives with respect to the r=pos2-pos1) 
        between valence orbital on atom1 and valence orbital on atom2 
        using tabulated Slater Koster integrals and transformation rules

        Paramters:
        ==========
        l1,m1: angular momentum and magnetic qnumber of valence orbital on atom 1
        l2,m2: angular momentum and magnetic qnumber of valence orbital on atom 2
        pos1: (x1,y1,z1) cartesian coordinates of atom 1
        pos2: (x2,y2,z2) cartesian coordinates of atom 2
        deriv: 0 -> overlap
               1 -> gradient of overlap
               2 -> hessian of overlap

        Returns:
        ========
        (Hamiltonian matrix element H) or (gradient of H) or (hessian of H)   (depending on deriv)
        """
        r,x,y,z = self._directional_cosines(pos1,pos2)
        if deriv == 0:
            H = T.slako_transformations[(l1,m1,l2,m2)](r,x,y,z, self.H_spl)
            """
            # test new Fortran code
            indeces = [i for (l1dummy,l2dummy,i) in self.slako.H.keys()]
            H_spl_arr = np.zeros(max(indeces)+1)
            for i in indeces:
                H_spl_arr[i] = self.H_spl(i,r)
            from DFTB.extensions import slako
            H_f90 = slako.slako.slako_transformation(r,x,y,z, H_spl_arr, l1,m1,l2,m2)
            # compare with python code
            assert abs(H-H_f90) < 1.0e-8, "H_py = %s   H_f90 = %s for (%s,%s,%s,%s)" % (H, H_f90, l1,m1,l2,m2)
            #
            """
            return H
        elif deriv == 1:
            gradH = list(zeros(3))
            for xyz in range(0,3):
                gradH[xyz] = T.slako_transformations_gradient[xyz][(l1,m1,l2,m2)](r,x,y,z, self.H_spl)
            """
            # test new Fortran code
            indeces = [i for (l1dummy,l2dummy,i) in self.slako.S.keys()]
            H_spl_arr = np.zeros(max(indeces)+1)   # SK integrals for H
            Hd_spl_arr = np.zeros(max(indeces)+1)  # derivatives of SK integrals for H
            for i in indeces:
                H_spl_arr[i] = self.H_spl(i,r)
                Hd_spl_arr[i] = self.H_spl(i,r,der=1)
            from DFTB.extensions import slako
            gradH_f90 = slako.slako.slako_transformation_gradient(r,x,y,z, H_spl_arr, Hd_spl_arr, l1,m1,l2,m2)
            # compare with python code
            assert np.sum(abs(gradH-gradH_f90)) < 1.0e-8, "grad H_py = %s   grad H_f90 = %s for (%s,%s,%s,%s)" % (gradH, gradH_f90, l1,m1,l2,m2)
            #
            """
            return gradH
        elif deriv == 2:
            hessH = list(zeros((3,3)))
            for xyz1 in range(0,3):
                for xyz2 in range(0,3):
                    hessH[xyz1][xyz2] = T.slako_transformations_hessian[xyz1][xyz2][(l1,m1,l2,m2)](r,x,y,z, self.H_spl)
            return hessH
        else:
            raise Exception("Higher than 2nd derivatives of S or H not implemented")
    def haveDipoles(self):
        """
        test if Slater Koster tables for dipoles are available
        """
        return hasattr(self.slako, 'Dipole')
    def havePhaseShifts(self):
        """
        test if Slater Koster tables for scattering states has phase shifts
        """
        return hasattr(self.slako, 'PhaseShifts')
    def getPhaseShift(self, l2):
        if self.havePhaseShifts():
            return self.slako.PhaseShifts[l2]
        else:
            return 1.0
    def getDipole(self, l1,m1,pos1, l2,m2,pos2, deriv=0):
        """
        compute dipole matrix elements between valence orbital on atom1 and valence orbitals on atom2
        using tabulated Slater Koster integrals and transformation rules

        Paramters:
        ==========
        l1,m1: angular momentum and magnetic qnumber of valence orbital on atom 1
        l2,m2: angular momentum and magnetic qnumber of valence orbital on atom 2
        pos1: (x1,y1,z1) cartesian coordinates of atom 1
        pos2: (x2,y2,z2) cartesian coordinates of atom 2

        Returns:
        ========
        dipole vector [Dx,Dy,Dz]
        """
        if l1 > Tdip.lmax or l2 > Tdip.lmax:
            #print "Warning: Slater Koster files for dipole matrix elements between orbitals with l > %d not implemented yet!" % Tdip.lmax
            return [0.0, 0.0, 0.0]
        r,x,y,z = self._directional_cosines(pos1,pos2)
        if not self.haveDipoles():
            raise NoSKDipolesException("No Slater Koster tables for dipoles available for these parameters (atom pair %s-%s)" \
              % (AtomicData.atom_names[self.slako.Z1-1], AtomicData.atom_names[self.slako.Z2-1]))
        ## HACK
        # The Slater-Koster rules for dipole matrix elements between d-orbitals
        # contain divisions by (x**2+y**2), which may raise a division by zero error.
        # To avoid this, x and y are shifted slightly away from 0 
        if l1 > 1 or l2 > 1:
            if (x**2+y**2 < 1.0e-20):
                x += 1.0e-20
                y += 1.0e-20
        ##
        if deriv == 0:
            Dx = Tdip.slako_transformations[(l1,m1,1, 1,l2,m2)](r,x,y,z, self.Dipole_spl)
            Dy = Tdip.slako_transformations[(l1,m1,1,-1,l2,m2)](r,x,y,z, self.Dipole_spl)
            Dz = Tdip.slako_transformations[(l1,m1,1, 0,l2,m2)](r,x,y,z, self.Dipole_spl)
            """
            # test new Fortran code
            indeces = [i for (l1dummy,l2dummy,i) in self.slako.Dipole.keys()]
            Dipole_spl_arr = np.zeros(max(indeces)+1)
            for i in indeces:
                Dipole_spl_arr[i] = self.Dipole_spl(i,r)
            from DFTB.extensions import slako
            Dx_f90 = slako.slako.slako_transformation_dipole(r,x,y,z, Dipole_spl_arr, l1,m1, 1, l2,m2)
            Dy_f90 = slako.slako.slako_transformation_dipole(r,x,y,z, Dipole_spl_arr, l1,m1,-1, l2,m2)
            Dz_f90 = slako.slako.slako_transformation_dipole(r,x,y,z, Dipole_spl_arr, l1,m1, 0, l2,m2)
            # compare with python code
            assert abs(Dx-Dx_f90) < 1.0e-8, "Dx_py = %s   Dx_f90 = %s for (%s,%s,%s,%s)" % (Dx, Dx_f90, l1,m1,l2,m2)
            assert abs(Dy-Dy_f90) < 1.0e-8, "Dy_py = %s   Dy_f90 = %s for (%s,%s,%s,%s)" % (Dy, Dy_f90, l1,m1,l2,m2)
            assert abs(Dz-Dz_f90) < 1.0e-8, "Dz_py = %s   Dz_f90 = %s for (%s,%s,%s,%s)" % (Dz, Dz_f90, l1,m1,l2,m2)
            #
            """
            #
            # <nu(r-R1)|r|mu(r-R2)> = <nu(r)|r|mu(r-(R2-R1))> + R1*<nu(r)|mu(r-(R2-R1))> = D + R1*S
            S = self.getOverlap(l1,m1,pos1, l2,m2,pos2, deriv=0)
            # apparently there is some confusion in the formula above,
            # if pos1 and pos2 are exchanged we get the expected result
            # for the dipole matrix elements for h-h
            """
            Dx += pos1[0]*S
            Dy += pos1[1]*S
            Dz += pos1[2]*S
            """
            Dx += pos2[0]*S
            Dy += pos2[1]*S
            Dz += pos2[2]*S
            return [Dx,Dy,Dz]
        else:
            raise Exception("Derivative of the dipole vector are not implemented")

# GLUE FOR FORTRAN CODE
def combine_slako_tables_f90(SKT):
    """
    Fortran cannot deal with python hashes. Therefore the dictionary SKT which contains
    the splines of the Slater-Koster tables for each pair (Z1,Z2) have to be combine into giant
    numpy array that can be passed to a Fortran function.

    The atoms are assigned to atom types 0,...,Ntype. Then Slater-Koster tables which are labelled
    by pairs of atomic numbers (Zi,Zj), are relabled by a single index running from 0,1,...,Npair
    where Npair is the number of unique pairs with Zi <= Zj:
     Npair = Ntype*(Ntype+1)/2
    The rules for translating the indeces are:
     pair of atomic numbers (Zi,Zj) -> pair of atom type (i,j) ->  pair index ipair = j*(j+1)/2 + i

    Nsk is the maximum number of Slater-Koster tables for any pair. For instance, the pair h-h has only one table
    for s-s-sigma type interaction, whereas the pair h-c has two tables, for ss-sigma and for sp-sigma
    type interactions, the pair zn-zn has 6 tables, ss-sigma, sd-sigma, ds-sigma, dd-delta, dd-pi and dd-sigma.
    Unfortunately tables are not enumerated consecutively, which is a little bit weird. In effect, the
    two tables for the h-c pair have indeces 0 and 2, while the table for index 1 is empty; the 6 tables for zn-zn
    have indeces 0,3,9,12,13 and 14 with the tables in between empty. 
    Therefore Nsk would be 15 for the atoms h,c,n,o and zn.
    NskD has a similar meaning for Slater-Koster tables for dipoles.

    Parameters:
    ===========
    SKT: dictionary, SKT[(Z1,Z2)] is an instance of SlakoTransformations, only the member variables
       S_tck and H_tck (and if available Dipole_tck) are used, which contain the knots and coefficients of the splines.

    Returns:
    ========
    atom_type_dic: dictionary that translates atomic numbers into atom types
    spline_deg: degree of B-splines
    tab_filled_SH0: shape (Npair,Nsk), if the table t contains data for atom pair ipair, 
       tab_filled_SH0[ipair,t] = 1 and otherwise 0
    tab_filled_D: shape (Npair,NskD), if the table t contains data for atom pair ipair, 
       tab_filled_D[ipair,t] = 1 and otherwise 0
    S_knots, S_coefs, H_knots, H_coefs, D_knots, D_coefs: numpy arrays with shape (Npair, Nsk, Npts) for S and H0 and (Npair, NskD, Npts) for dipoles
      that contain the knot positions and B-spline coefficients for the overlap S, hamiltonian H0 and dipole integrals D
    """
    # unique atomic numbers
    Zs = list(set([Z1 for (Z1,Z2) in SKT.keys()]))
    Zs.sort()
    Ntype = len(Zs)
    # The atom with the lowest atomic numbers becomes type 0, the next one becomes type 1, and so on
    atom_type_dic = {}
    for itype,Zi in enumerate(Zs):
        atom_type_dic[Zi] = itype
    # number of unique pairs
    Npair=(Ntype*(Ntype+1))/2
    # maximum number of tables for S and H0
    Nsk=1
    # maximum number of tables for D if available
    NskD=1
    # number of knots in the B-spline
    Npts=None
    # degree of B-splines
    spline_deg = None
    for (Zi,Zj),sktab in SKT.iteritems():
        # S and H0
        tab_indeces = sktab.S_tck.keys()
        Nsk=max(Nsk, max(tab_indeces))+1
        knots,coefs,deg = sktab.S_tck[tab_indeces[0]]
        # The knots have to be spaced uniformly
        dx = np.ediff1d(knots)[deg+1:-(deg+1)]
        assert np.std(dx) < 1.0e-7, "ERROR: Grid for Slater Koster table of atom pair (%d,%d) is not uniformly spaced, deviation = %s" % (Zi,Zj, np.std(dx))
        if Npts == None:
            Npts = len(knots)
            spline_deg = deg
        else:
            # All tables should have been calculated for the same grid
            assert Npts == len(knots)
            assert spline_deg == deg
        # dipoles
        if sktab.haveDipoles():
            tab_indeces = sktab.Dipole_tck.keys()
            NskD=max(NskD, max(tab_indeces))+1
    # fill numpy arrays with the SK tables
    sh=(Npair,Nsk,Npts)
    S_knots = np.zeros(sh)
    S_coefs = np.zeros(sh)
    H_knots = np.zeros(sh)
    H_coefs = np.zeros(sh)
    # dipoles
    shD=(Npair,NskD,Npts)
    D_knots = np.zeros(shD)
    D_coefs = np.zeros(shD)
    tab_filled_SH0 = np.zeros((Npair,Nsk),dtype=int)
    tab_filled_D = np.zeros((Npair,NskD),dtype=int)
    for (Zi,Zj),sktab in SKT.iteritems():
        i,j = atom_type_dic[Zi], atom_type_dic[Zj]
        ipair = (j*(j+1))/2+i
        # mark filled tables
        tab_indeces = sktab.S_tck.keys()
        tab_filled_SH0[ipair,tab_indeces] = 1
        if sktab.haveDipoles():
            tab_indeces = sktab.Dipole_tck.keys()
            tab_filled_D[ipair,tab_indeces] = 1
        # overlap S
        for isk,(knots,coefs,deg) in sktab.S_tck.iteritems():
            S_knots[ipair,isk,:] = knots
            S_coefs[ipair,isk,:] = coefs
        # hamiltonian H0
        for isk,(knots,coefs,deg) in sktab.H_tck.iteritems():
            H_knots[ipair,isk,:] = knots
            H_coefs[ipair,isk,:] = coefs
        # dipoles
        if sktab.haveDipoles():
            for isk,(knots,coefs,deg) in sktab.Dipole_tck.iteritems():
                D_knots[ipair,isk,:] = knots
                D_coefs[ipair,isk,:] = coefs
            
    return atom_type_dic, spline_deg, tab_filled_SH0, tab_filled_D, (S_knots, S_coefs, H_knots, H_coefs, D_knots, D_coefs)
        
    
    
def test_slako_transformations():
    from matplotlib.pyplot import plot, show
    from slako_tables import c_c
    SK = SlakoTransformations(c_c)
    x = linspace(0.0, 10.0, 100)
    y = linspace(0.0, 10.0, 100)
    z = linspace(0.0, 10.0, 100)
    r = sqrt(x*x+y*y+z*z)
    print SK.getOverlap(1,-1,(0,0,0),0,0,(x[10],y[10],z[10]))
    print SK.getOverlap(1,0,(0,0,0),0,0,(1,0,0),deriv=1)
    print SK.getOverlap(1,0,(0,0,0),0,0,(1,0,0),deriv=2)
    print SK.getHamiltonian0(1,0,(0,0,0),0,0,(1,0,0),deriv=1)
    print SK.getHamiltonian0(1,0,(0,0,0),0,0,(1,0,0),deriv=2)

    show()

if __name__ == "__main__":
    test_valence_orbitals()
#    test_slako_transformations()

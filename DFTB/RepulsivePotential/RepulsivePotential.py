"""

"""
from numpy import array, where, reshape, zeros, linspace, hstack, exp, argmax
import numpy as np
from scipy import interpolate
import os.path
import imp

from DFTB.SlaterKoster import hotbit_format
from DFTB import AtomicData
from DFTB import utils

############# USAGE OF REPULSIVE POTENTIAL TABLES ######################
from collections import defaultdict

class RepulsivePotential:
    def __init__(self, reppot_module, scaling_factors=defaultdict(lambda: 1.0,{}), smooth_decay=False):
        """
        Parameters:
        ===========
        reppot_module should be an object with the following members
        
        Z1,Z2: atomic numbers of atom pair
        d: 1D numpy array, distance between atoms
        Vrep: repulsive potential on the grid d

        smooth_decay controls whether Vrep and its derivatives are set abruptly to
        0 after the cutoff radius or whether a smooting function is added.
        WARNING: the smooting function can change the nuclear repulsion energy between
        atoms that are far apart. Therefore you should check visually the tails of
        the repulsive potential and check that the additional energy is not negligible.
        """
        # the core potential can be made more or less repulsive by scaling the
        # distances with a factor > 1 or < 1
        Z1,Z2 = reppot_module.Z1, reppot_module.Z2
        Z1,Z2 = min(Z1,Z2), max(Z1,Z2)
        scaling=scaling_factors[(Z1,Z2)]
        reppot_module.d *= scaling
        # cut-off radius
        self.Rcut = reppot_module.d[-1]
        if smooth_decay:
            # To avoid discontinuities we attach a smooth function 
            # with matched 1st and 2nd derivative at dmax that goes slowly to zero
            # The tail has the functional form g(d) = Vrep(dmax) * exp(-a1*(d-dmax) - a2*(d-dmax)^2)
            # a1 and a2 are determined such that
            #     g(dmax) = Vrep(dmax),  g'(dmax) = Vrep'(dmax) and g''(dmax) = Vrep''(dmax)

            # find largest distance where Vrep is not zero
            dmax_indx = argmax(reppot_module.d[abs(reppot_module.Vrep) >= 0.1])
            self.dmin = min(reppot_module.d)
            self.dmax = reppot_module.d[dmax_indx]

            self.tck = interpolate.splrep(reppot_module.d[:dmax_indx], reppot_module.Vrep[:dmax_indx], s=0)

            Vrep0 = interpolate.splev(self.dmax, self.tck, der=0)
            Vrep1 = interpolate.splev(self.dmax, self.tck, der=1)
            Vrep2 = interpolate.splev(self.dmax, self.tck, der=2)
            Vrep3 = interpolate.splev(self.dmax, self.tck, der=3)

            #assert Vrep1 < 0.0
            a1 = -Vrep1/Vrep0
            a2 = 0.5*(pow(a1,2) - Vrep2/Vrep0)
            a3 = abs(a2/a1) #0.1 # make sure cubic term in g(x) is positive so that g(x)
            # decays to 0 independently of a1 or a2
            assert a3 >= 0.0
            def smooth_tail(d,deriv=0):
#                print "dmax = %s deriv = %s a1=%s a2=%s a3=%s" % (self.dmax, deriv, a1,a2,a3)
                g = Vrep0 * exp(-a1*(d-self.dmax) -a2*pow(d-self.dmax,2) - a3*pow(d-self.dmax,3))
                if deriv == 0:
                    return g
                elif deriv == 1:
                    gp = g*(-a1 - 2.0*a2*(d-self.dmax) - 3.0*a3*pow(d-self.dmax,2))
                    return gp
                elif deriv == 2:
                    gpp = g*(pow(-a1 - 2.0*a2*(d-self.dmax) - 3.0*a3*pow(d-self.dmax,2),2) - 2.0*a2 - 6.0*a3*(d-self.dmax))
                    return gpp
                else:
                    raise Exception("If higher derivatives are needed, add a3,a4,...and match appropriately")
        else:
            self.dmin = min(reppot_module.d)
            self.dmax = max(reppot_module.d)
            self.tck = interpolate.splrep(reppot_module.d, reppot_module.Vrep,s=0)
            def smooth_tail(d,deriv=0):
                return 0.0
        self.smooth_tail = smooth_tail
    def getVrep(self, d, deriv=0):
        """
        interpolate repulsive potential (or its derivatives)
        """
        if d < self.dmin:
            print "WARNING: Repulsive potential is tabulated for %s <= r, got r = %s, probably atoms are too close!" % (self.dmin, d)
        if d >= self.dmax:
            v = self.smooth_tail(d,deriv=deriv)
        else:
            v = interpolate.splev(d,self.tck,der=deriv)
        return v
    def plotVrep(self, deriv=0):
        from numpy import frompyfunc
        from matplotlib.pyplot import plot, show, xlim
        d = linspace(0.0, 10*self.dmax, 1000)
        xlim((0.0,d.max()))
        Vrep_ufunc = frompyfunc(lambda x: self.getVrep(x,deriv=deriv), 1,1)
        Vrep = Vrep_ufunc(d)
        plot(d, Vrep, label=label)
        show()


class ReppotModule():
    pass

def dummy_reppot_module():
    """
    create a module that has the same interface as a ReppotModule but which is empty
    """
    print "DUMMY MODULE !!!!!!!!!!!!!"
    mod = ReppotModule()
    mod.Z1,mod.Z2 = 0,0
    mod.d = linspace(0.0, 10.0, 100)
    mod.Vrep = 50.0*exp(-mod.d)
    return mod

def read_repulsive_potential(filename):
    """
    read repulsive potential in the format used by Hotbit, DFTBaby or DFTB+
    """
    from os.path import exists, basename, dirname, join
    if ".par" in filename:
        # Hotbit's format
        at1, at2 = basename(filename).replace(".par", "").split("_")
        if not exists(filename):
            # file file A_B.par does not exist, look for file B_A.par
            filename_BA = join(dirname(filename), "%s_%s.par" % (at2,at1))
            print "parameter file %s does not exist -> try to load %s" % (filename, filename_BA)
            filename = filename_BA
        data_blocks = hotbit_format.parseHotbitParameters(filename)
        mod = ReppotModule()
        mod.d = data_blocks["repulsive_potential"][:,0]
        mod.Vrep = data_blocks["repulsive_potential"][:,1]
        mod.Z1 = AtomicData.atomic_number(at1)
        mod.Z2 = AtomicData.atomic_number(at2)
        return mod
    elif ".py" in filename:
        # DFTBaby format
        # access dictionary by .-notation
        mod = utils.dotdic()
        execfile(filename, mod)
        return mod
    elif ".skf" in filename:
        # DFTB+ format, only the Splines part is read
        at1, at2 = basename(filename).replace(".skf", "").split("-")
        if not exists(filename):
            # file file A-B.skf does not exist, look for file B-A.skf
            filename_BA = join(dirname(filename), "%s-%s.skf" % (at2,at1))
            print "parameter file %s does not exist -> try to load %s" % (filename, filename_BA)
            filename = filename_BA
        fh = open(filename)
        lines = fh.readlines()
        fh.close()
        # find section with repulsive potential
        for i,l in enumerate(lines):
            if l.strip() == "Spline":
                break
        else:
            raise Exception("No 'Spline' section found in parameter file %s" % filename)
        # Line 2:
        nInt, cutoff = lines[i+1].strip().split()
        nInt, cutoff = int(nInt), float(cutoff)
        # Line 3: V(r < r0) = exp(-a1*r+a2) + a3   is r too small to be covered by the spline
        a1, a2, a3 = map(float, lines[i+2].strip().split())
        # Line 4 to 4+nInt-2
        rs = np.zeros(nInt)
        cs = np.zeros((4,nInt))
        for j in range(0, nInt):
            # spline for the range [rj=start,r_(j+1)=end]
            # V(r_j <= r < r_(j+1)) = c0 + c1*(r-r0) + c2*(r-r0)^2 + c3*(r-r0)^3
            start, end, c0, c1, c2, c3 = map(float, lines[i+3+j].strip().split()[:6])
            rs[j] = start
            cs[:,j] = np.array([c0,c1,c2,c3])
        # V(r_nInt < r) = 0.0
        assert end == cutoff
        # Now we evaluate the spline on a equidistant grid
        Npts = 100
        d = np.linspace(0.0, cutoff, Npts)
        Vrep = np.zeros(Npts)
        j = 0
        for i,di in enumerate(d):
            if (di < rs[0]):
                Vrep[i] = np.exp(-a1*di+a2) + a3
            else:
                # find interval such that r[j] <= di < r[j+1]
                while (di >= rs[j+1]) and j < nInt-2:
                    j += 1
                if j < nInt-2:
                    assert rs[j] <= di < rs[j+1]
                    c0,c1,c2,c3 = cs[:,j]
                    Vrep[i] = c0 + c1*(di-rs[j]) + c2*(di-rs[j])**2 + c3*(di-rs[j])**3
                else:
                    Vrep[i] = 0.0
        # create python module
        mod = ReppotModule()
        mod.d = d
        mod.Vrep = Vrep
        mod.Z1 = AtomicData.atomic_number(at1)
        mod.Z2 = AtomicData.atomic_number(at2)
        return mod
    else:
        raise Exception("Format of %s not understood" % filename)

    
def write_repulsive_potential(mod, reppot_file, title=""):
    # NOTE: a similar function for writing repulsive potentials
    # exists in FitForces.py
    fh = open(reppot_file, "w")
    atname1 = AtomicData.atom_names[mod.Z1-1]
    atname2 = AtomicData.atom_names[mod.Z2-1]

    import pprint
    import sys
    import numpy
    numpy.set_printoptions(threshold=sys.maxint)
    pp = pprint.PrettyPrinter(depth=10)
    fh = open(reppot_file, "w")
    print>>fh, "# %s" % title
    print>>fh, "from numpy import array"
    print>>fh, "Z1 = %s" % mod.Z1
    print>>fh, "Z2 = %s" % mod.Z2
    print>>fh, "# grid for distance d between atomic centers in bohr"
    print>>fh, "d = \\\n%s" % pp.pformat(mod.d)
    print>>fh, "# repulsive potential in hartree/bohr"
    print>>fh, "Vrep = \\\n%s" % pp.pformat(mod.Vrep)
        
    fh.close()


def load_repulsive_potentials(atpairs, missing_reppots="error", reppot_paths=[]):
    """
    Import repulsive potentials for all atom pairs (A,B).  
    Folders are searched for the repulsive potential modules `a_b.py`
    in the order
     1) any folder specified in `reppot_paths`
     2) The default location for repulsive potentials, 'DFTB/RepulsivePotentials/reppot_tables/'

    Parameters
    ----------
    atpairs         : list of tuples (Z1,Z2) with all unique combinations of atomic numbers present
                      in the molecule
    
    Optional
    --------
    missing_reppots : How should missing repulsive potentials be handled? 
                      "error" - stop with an error message or
                      "dummy" - import a dummy potential

    Returns
    -------
    reppot_tables   : module with a repulsive potentials in the member variable, `reppot_tables.atompairs`
    """
    atompairs = {}
    for (Zi,Zj) in atpairs:
        atI = AtomicData.atom_names[Zi-1]
        atJ = AtomicData.atom_names[Zj-1]
        # First we look for repulsive potentials in non-standard locations
        for path in reppot_paths:
            try:
                name = "%s_%s" % (atI,atJ)
                module_path = os.path.relpath(os.path.join(path, name+".py"))
                atompairs[(Zi,Zj)] = imp.load_source(name, module_path)
            except IOError as e:
                continue
            break
        else:
            # We could not find the repulsive potential anywhere else, so we look in the
            # standard location.
            try:
                atompairs[(Zi,Zj)] = __import__("DFTB.RepulsivePotential.reppot_tables.%s_%s" % (atI,atJ), fromlist=['Z1','Z2'])
            except ImportError as e:
                if missing_reppots == "dummy":
                    print "WARNING: Repulsive potential for atom pair %s-%s not found!" % (atI,atJ)
                    print "         A dummy potential with Vrep=0 is loaded instead. Optimizations and"
                    print "         dynamics simulations with dummy potentials are bound to produce nonsense."
                    atompairs[(Zi,Zj)] = __import__("DFTB.RepulsivePotential.reppot_tables.dummy", fromlist=['Z1','Z2'])
                    # set correct atomic numbers
                    atompairs[(Zi,Zj)].Z1 = Zi
                    atompairs[(Zi,Zj)].Z1 = Zj
                else:
                    raise e
            continue
        print "repulsive potential for %s-%s loaded from %s" % (atI.upper(), atJ.upper(), module_path)
    reppot_tables = imp.new_module("reppot_tables")
    reppot_tables.atompairs = atompairs

    return reppot_tables
        
        
if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print "Usage: %s <.par file with repulsive potential>" % sys.argv[0]
        print "Plot repulsive potential"
        exit(-1)
    filename = sys.argv[1]
    reppot_mod = read_repulsive_potential(filename)
    RepPot = RepulsivePotential(reppot_mod)
    RepPot.plotVrep(deriv=0)

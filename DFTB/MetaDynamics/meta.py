#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
from sympy.utilities.lambdify import lambdify
import multiprocessing as mp

import cvfunctions as cvf
if os.path.exists("cvcustom.py"):
    import cvcustom as cvc
import mulliken as mul

au2ang = 0.52917721067
kb = (1.38064852e-23)/(4.35974465e-18)

class metadynamics:
    def __init__(self, symbols=None, infile="meta-config.py", restart=False):
        """Create metadynamics object
        
        Parameters
        ----------
        symbols : list
            element symbols needed for calculation of coordination numbers (default: None)
        inputfile : string
            name of the parameter input file (default: meta-config.py)
        restart : boolean
            read gaussian centers from vg_centers.dat to restart previous run (default: False)
        
        """
        self.symbols = symbols
        self.infile = infile
        self.restart = restart
        self.cvs = []
        self.welltemp = False

    def read_input(self):
        """Parse input from the input file"""
        with open(self.infile) as f:
            self.cfg = eval("".join(f.readlines()))
        self.ncvs = len(self.cfg["collective variables"])
        for idx, cv in enumerate(self.cfg["collective variables"]):
            if cv["type"] == "bond":
                self.cvs.append(cv_bond(idx, cv, self.symbols))
            elif cv["type"] == "angle":
                self.cvs.append(cv_angle(idx, cv, self.symbols))
            elif cv["type"] == "torsion":
                self.cvs.append(cv_torsion(idx, cv, self.symbols))
            elif cv["type"] == "cn":
                self.cvs.append(cv_cn(idx, cv, self.symbols))
            elif cv["type"] == "mullcharge":
                self.cvs.append(cv_mullcharge(idx, cv, self.symbols))
            elif cv["type"] == "custom":
                self.cvs.append(cv_custom(idx, cv, self.symbols))
            if "plot" in cv.keys():
                self.cvs[-1].set_plot(cv["plot"])
        if "well-tempered" in self.cfg.keys():
            self.welltemp = True
            self.deltat = self.cfg["well-tempered"]["deltaT"]
            self.heights = []

    def read_cfg(self):
        lines = []
        with open(self.infile) as f:
#            lines = f.readlines()
            for line in f:
                lines += line.split("#")[0]
        print lines
        

    def set_pes_object(self, pes):
        """
        pes is an instance of DFTB.PotentialEnergySurface
        """
        for cv in self.cvs:
            if cv.type == "mullcharge":
                cv.set_pes_object(pes)

    def initialize_vg(self):
        """Initialize sympy expressions and output file for the metadynamics potential vg"""
        # create list of sympy variables xi and sympy symbol for vg
        if self.ncvs > 1:
            self.x = sym.symbols("x"+" x".join([str(i) for i in range(self.ncvs)]))
        else:
            self.x = [sym.symbols("x")]
        # initialize vg and lamdvg
        self.vg = 0.0
        self.set_lamdvg()
        if self.restart:
            self.reload_vg()
        else:
        # file in which the centers of the gaussians (s) are stored, needed for reconstruction
            with open("vg_centers.dat", "w") as outf:
                outf.write("# gaussian centers of vg, for parameters see "+str(self.infile)+"\n")
        # file in which gaussian heights are stored if well-tempered simulation
            if self.welltemp:
                with open("vg_heights.dat", "w") as outf:
                    outf.write("# gaussian heights of vg\n")

    def reload_vg(self):
        with open("vg_centers.dat") as f:
            f.next()
            for line in f:
                self.s = [float(s) for s in line.split()]
                self.set_vg()

    def set_s_and_ds(self):
        """Calculate values of CVs s and derivatives wrt coordinates ds_dr"""
        for cv in self.cvs:  # delete old tmp files if they exist
            cv.delete_tmpfile()
        self.s = np.zeros(self.ncvs)
        self.ds_dr = np.zeros((self.ncvs, len(self.coordinates), 3))
        for i, cv in enumerate(self.cvs):
            self.s[i], self.ds_dr[i] = cv.get_s_and_ds(self.coordinates)

    def get_height(self):
        if self.welltemp:
            lamvg = lambdify(self.x, self.vg)
            height = self.cfg["height"]*np.exp(-lamvg(*self.s)/(kb*self.deltat))
            self.heights.append(height)
            return height
        else:
            return self.cfg["height"]

    def set_vg(self):
        """Add new gaussian to current meta potential vg"""
        def sym_gauss(i, cv):
            return sym.exp(-(self.x[i]-self.s[i])**2/2/cv.width**2)
        g = [sym_gauss(i, cv) for i, cv in enumerate(self.cvs)]
        self.vg += self.get_height()*np.prod(g)

    def set_lamdvg(self):
        """Differentiate vg wrt s and lambdify derivative"""
        self.lamdvg = [lambdify(self.x, sym.diff(self.vg, self.x[i])) for i in range(self.ncvs)]

    def get_forces(self, coords, step):
        """
        Add new gaussian to vg every tau_g steps and return metadynamics forces
        
        Parameters
        ----------
        coords : numpy array
            coordinates (shape: (N, 3))
        step : number
            step number of the dynamics simulation
        
        Returns
        -------
        forces : numpy array
            forces obtained from the current metadynamics potential
        
        """
        # update coordinates, values of the CVs s and derivatives ds_dr
        self.coordinates = coords
        self.set_s_and_ds()
        # if metadynamics step, update lamdvg and write s to vg_centers.dat
        if not step % self.cfg["tau_g"]:
            fmt = self.ncvs*"%12.6f"+"\n"
            self.set_vg()
            self.set_lamdvg()
            with open("vg_centers.dat", "a") as outf:
                s_converted = [cv.convert_units(self.s[i]) for i, cv in enumerate(self.cvs)]
                outf.write(fmt % tuple(s_converted))
            if self.welltemp:
                with open("vg_heights.dat", "a") as outf:
                    outf.write("%12.6f\n" % (self.heights[-1]))
        # get metadynamics forces
        forces = np.zeros(self.coordinates.shape)
        for i, cv in enumerate(self.cvs):
            if cv.type == "torsion":
                # torsion angles are defined in the range from -180 to 180 degrees. in order
                # to account for periodicity, dvg_ds is additionally calculated in the ranges
                # from -540 to -180 and from 180 to 540 degrees and summed up with the values
                # in the original range.
                tmp = np.array([self.s, self.s, self.s])
                tmp[:, i] += np.array([-1, 0, 1])*360
                dvg_ds = np.sum([self.lamdvg[i](*t) for t in tmp])
            else:
                dvg_ds = self.lamdvg[i](*self.s)
            forces -= dvg_ds*self.ds_dr[i]
        return forces


class cv():
    def __init__(self, idx, cfg, symb):
        self.idx = idx
        self.tmp = "tmp.npy"
        self.name = cfg["name"]
        self.type = cfg["type"]
        self.width = cfg["width"]
        self.symbols = symb

    def delete_tmpfile(self):
        if os.path.exists(self.tmp):
            os.remove(self.tmp)

    def set_plot(self, plot):
        self.x = np.linspace(plot["min"], plot["max"], plot["npoints"])
        self.ticks = plot["ticks"]


class cv_bond(cv):
    def __init__(self, idx, cfg, symb):
        cv.__init__(self, idx, cfg, symb)
        self.atoms = cfg["atoms"]
        self.x = np.linspace(0.5, 2.5, 200)
        self.ticks = None

    def get_s_and_ds(self, coords):
        coord = [coords[i] for i in self.atoms]
        bond = cvf.bond(coord)
        dbond = np.zeros(coords.shape)
        dbond[self.atoms] = cvf.dbond(coord)
        return bond, dbond

    def convert_units(self, s):
        return s*au2ang


class cv_angle(cv):
    def __init__(self, idx, cfg, symb):
        cv.__init__(self, idx, cfg, symb)
        self.atoms = cfg["atoms"]
        self.x = np.linspace(0., 180., 200)
        self.ticks = 30.

    def get_s_and_ds(self, coords):
        coord = [coords[i] for i in self.atoms]
        angle = cvf.angle(coord)
        dangle = np.zeros(coords.shape)
        dangle[self.atoms] = cvf.dangle(coord)
        return angle, dangle

    def convert_units(self, s):
        return s


class cv_torsion(cv):
    def __init__(self, idx, cfg, symb):
        cv.__init__(self, idx, cfg, symb)
        self.atoms = cfg["atoms"]
        self.x = np.linspace(-180., 180., 200)
        self.ticks = 45.

    def get_s_and_ds(self, coords):
        coord = [coords[i] for i in self.atoms]
        torsion = cvf.torsion(coord)
        dtorsion = np.zeros(coords.shape)
        dtorsion[self.atoms] = cvf.dtorsion(coord)
        return torsion, dtorsion

    def convert_units(self, s):
        return s


class cv_cn(cv):
    def __init__(self, idx, cfg, symb):
        cv.__init__(self, idx, cfg, symb)
        self.atom = cfg["atom"]
        if "reference" in cfg.keys():  # if reference atoms are defined, get indices
            self.ref = cfg["reference"]
        else:  # if not, list of reference atoms has to be every other atom
            self.ref = range(self.atom)+range(self.atom+1, len(self.symbols))
        self.n = cfg["n"]
        self.m = cfg["m"]
        if cfg["d"] == "auto":  # if set to auto, determin d for every atom pair from cov. radii
            self.d = [cvf.get_d(self.symbols, self.atom, atom2) for atom2 in self.ref]
        elif type(cfg["d"]) == type([]):  # if list is parsed, take it as it is
            self.d = list(np.array(cfg["d"])/au2ang)
        elif type(cfg["d"]) == type(0.0):  # if float is parsed, use it for all atom pairs
            self.d = len(self.ref)*[cfg["d"]/au2ang]
        print "d: ", self.d
        self.x = np.linspace(0., 1., 200)
        self.ticks = None

    def get_s_and_ds(self, coords):
        return cvf.cn(self.symbols, coords, self.atom, self.n, self.m, self.d, self.ref)

    def convert_units(self, s):
        return s

    def check_parameters(self):
        if self.d == "auto":
            atom1 = raw_input("Please enter the first atom type: ").capitalize()
            atom2 = raw_input("Please enter the second atom type: ").capitalize()
            d = cvf.get_d([atom1, atom2], 0, 1)
            plt.title(r"%s-%s ($d = %.2f \:\AA$)" % (atom1, atom2, d*au2ang))
        else:
            d = self.d/au2ang
            plt.title(r"$d = %.2f \:\AA$" % (d*au2ang))
        x = np.linspace(0., 4., 100)
        cn = cvf.cn_i(x/au2ang, d, self.n, self.m)
        plt.xlabel(r"bond length / $\AA$")
        plt.ylabel(r"coordination number")
        plt.axvline(x=d*au2ang, color="grey", ls="dashed")
        plt.plot(x, cn)
        plt.tight_layout(pad=0.1)
        plt.savefig(self.name.replace(" ", "-")+".png", dpi=600)
        plt.show()


class cv_mullcharge(cv):
    def __init__(self, idx, cfg, symb):
        cv.__init__(self, idx, cfg, symb)
        if type(cfg["atoms"]) == type([]):
            self.atoms = cfg["atoms"]
        else:
            self.atoms = [cfg["atoms"]]
        self.x = np.linspace(-1.0, 1.0, 200)
        self.ticks = None
        self.method = cfg["method"]
        self.tmp = "tmp_charges.npy"
        if "nproc" in cfg.keys():
            self.nproc = cfg["nproc"]
        else:
            self.nproc = 1

    def set_pes_object(self, pes):
        self.pes = pes

    def get_s_and_ds(self, coords):
        nat = len(self.symbols)
        charge = np.sum(self.pes.tddftb.dftb2.dq[self.atoms])
        if self.method == "analytic":
            dcharges = self.pes.grads.getChargeGradients()
            dcharge = np.sum(dcharges[:, self.atoms], axis=1).reshape(nat, 3)
            return charge, dcharge
        # numerical charge gradients:
        h = 1e-3
        dcharge = np.zeros(coords.shape)
        if os.path.exists(self.tmp):
            charges_h = np.load(self.tmp)
            for i in range(3*nat):
                dcharge[i/3, i%3] = (np.sum(charges_h[i/3, i%3, self.atoms])-charge)/h
        else:
            pool = mp.Pool(processes=self.nproc)
            coords_h = []
            for i in range(3*nat):
                coords_h.append(np.copy(coords))
                coords_h[i][i/3, i%3] += h
            i_sym_coords_h = [(i, self.symbols, coord) for i, coord in enumerate(coords_h)]
            charges_h = np.array(pool.map(mul.get_new_charges_dftb, i_sym_coords_h))
            pool.close()
            pool.join()
            dcharge = (np.sum(charges_h[:, self.atoms], axis=1).reshape(nat, 3)-charge)/h
            np.save(self.tmp, charges_h)
            print "charge: ", charge
            print "dcharge:\n", dcharge
        return charge, dcharge

    def convert_units(self, s):
        return s


class cv_custom(cv):
    def __init__(self, idx, cfg, symb):
        cv.__init__(self, idx, cfg, symb)
        self.prm = cfg["params"]
        self.x = np.linspace(0., 1.0, 200)
        self.ticks = None

    def get_s_and_ds(self, coords):
        custom = cvc.customcv(self.symbols, coords, self.prm)
        dcustom = cvc.dcustomcv(self.symbols, coords, self.prm)
        return custom, dcustom

    def convert_units(self, s):
        return s

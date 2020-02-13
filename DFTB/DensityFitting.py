"""

"""
import numpy as np
import numpy.linalg as la

from DFTB.BasisSets import AtomicBasisSet, AuxiliaryBasisSet
from DFTB import AtomicData
from DFTB import XYZ

class DensityFitting:
    def __init__(self, atomlist, hubbard_U):
        """
        
        Parameters:
        ===========
        bfs: list of numeric basis functions
        bfs_aux: list of auxiliary basis functions
        """
        self.atomlist = atomlist
        self.bfs = AtomicBasisSet(atomlist).bfs
        self.bfs_aux = AuxiliaryBasisSet(atomlist, hubbard_U).bfs
        # create a list of points that should have the symmetry
        # of the molecule on which the true density should agree
        # with the model density
        Xs,Ys,Zs = [], [], []
        # atomic centers
        for Zi,posi in atomlist:
            x,y,z = posi
            Xs.append(x)
            Ys.append(y)
            Zs.append(z)
        # points between atoms
        for i,(Zi,posi) in enumerate(atomlist):
            for j,(Zj,posj) in enumerate(atomlist):
                if i == j:
                    continue
                Rij = np.array(posj)-np.array(posi)
                x,y,z = np.array(posi) + 0.33*Rij
                Xs.append(x)
                Ys.append(y)
                Zs.append(z)
                x,y,z = np.array(posi) - 0.33*Rij
                Xs.append(x)
                Ys.append(y)
                Zs.append(z)


        # points perpendicular to plane between 3 atoms
        for i,(Zi,posi) in enumerate(atomlist):
            for j,(Zj,posj) in enumerate(atomlist):
                if i == j:
                    continue
                for k,(Zk,posk) in enumerate(atomlist):
                    if i == k:
                        continue
                    if i == j:
                        continue
                    Rij = np.array(posj) - np.array(posi)
                    Rik = np.array(posk) - np.array(posi)
                    x,y,z = np.array(posi) + 0.25*np.cross(Rij, Rik)
                    Xs.append(x)
                    Ys.append(y)
                    Zs.append(z)

        Xs,Ys,Zs = np.mgrid[-6.0:6.0:30j,-6.0:6.0:30j,-6.0:6.0:30j]
        Xs = Xs.ravel()
        Ys = Ys.ravel()
        Zs = Zs.ravel()
        grid = np.array(Xs), np.array(Ys), np.array(Zs)

        # number of fit parameters
        self.Nfit = len(self.bfs_aux)
        assert self.Nfit % 4 == 0
        # number of basis functions
        self.Nbfs = len(self.bfs)
        # number of sample points
        self.Npts = len(Xs)
        assert self.Npts > self.Nfit
        # save grid to xyz file
        gridlist = atomlist[:]
        for i in range(0, self.Npts):
            gridlist.append( (1, (Xs[i], Ys[i], Zs[i])) )
        XYZ.write_xyz("/tmp/fitgrid.xyz", [gridlist])
        # evaluate all basis functions on the grid
        self.bf_grid = np.zeros((self.Npts, self.Nbfs))
        for m,bfm in enumerate(self.bfs):
            self.bf_grid[:,m] = bfm.amp(*grid)        
        # evaluate fit functions on the grid
        self.bf_aux_grid = np.zeros((self.Npts, self.Nfit))
        for m,bfm_aux in enumerate(self.bfs_aux):
            self.bf_aux_grid[:,m] = bfm_aux.amp(*grid)

        #
        M = self.bf_aux_grid
        self.M = M
         # constrained   C.x - Q = 0
        Id4 = np.eye(4)
        C = np.hstack([Id4 for i in range(0, self.Nfit/4)])
        Ct = C.transpose()
        #
        Mt = M.transpose()
        MtM = np.dot(Mt, M)
        MtMinv = la.inv(MtM)
        #
        self.A1 = np.dot(C,np.dot(MtMinv, Mt))
        self.A2 = np.dot(C,np.dot(MtMinv, Ct))
        self.A3 = np.dot(MtMinv, Ct)
        self.A4 = np.dot(MtMinv,Mt)
    def fit(self, P, Qtot, Dtot):
        """
        Parameters:
        ===========
        P: (partial) density matrix
        Qtot: total charge
        Dtot: total dipole moment
        """
        # compute true density
        y = np.zeros(self.Npts)
        for m in range(0, self.Nbfs):
            for n in range(0, self.Nbfs):
                y += self.bf_grid[:,m] * P[m,n] * self.bf_grid[:,n]
        # linear least square fit
        Q = np.array([Qtot, Dtot[0], Dtot[1], Dtot[2]])
        print "Q = %s" % Q
        # solve for Lagrange multiplier
        c = Q - np.dot(self.A1,y)
        lagr = la.solve(self.A2,c)
        print "Lagrange multiplier = %s" % lagr
        #
        x = np.dot(self.A4,y) + np.dot(self.A3,lagr)
        print "  Density Fitting"
        print "  ==============="
        print "  Atom     dQ             partial dipole"
        for i,(Zi,posi) in enumerate(self.atomlist):
            print "  %s%d   %+5.3f  [%+5.3f %+5.3f %+5.3f]" \
                % (AtomicData.atom_names[Zi-1], i, x[4*i+0], x[4*i+1], x[4*i+2], x[4*i+3])

        # compare fit density with true density at the fit points
        yfit = np.dot(self.M,x)
        err = np.sum(abs(y-yfit))
        """
        print "RHO"
        print y
        print "RHO fit"
        print yfit
        """
        print "|y-yfit| = %s" % err
        print "|y-0   | = %s" % np.sum(abs(y))
        return x
    

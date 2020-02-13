"""
The symmetry label of an excited state can be inferred from the sign
changes of the transition density under the application of the symmetry
operations. For non-degenerate irreps the transition density at most changes
sign. The sign changes allow to map each state to a row in the character table
of the symmetry group.
Degenerate irreps have to be assigned by visual inspection
"""
import numpy as np
from DFTB.Analyse.Cube import electron_density
from DFTB.BasisSets import AtomicBasisSet
import time

class SymmetryAssignmentEx:
    """
    assign symmetry labels to excited states
    """
    def __init__(self, tddftb, symmetry_group):
        self.tddftb = tddftb
        self.dftb = tddftb.dftb2
        if self.dftb != None:
            self.bs = AtomicBasisSet(self.dftb.atomlist)
        self.symmetry_group = symmetry_group
        # How do I choose the test points?
        npts = 10
        # random points in a box [-3,3]x[-3,3]x[-3,3] around the first atom
        pos_at1 = np.array(self.dftb.atomlist[0][1]) # position of 1st atom
        self.test_pts = pos_at1 + 3.0 * 2.0*(np.random.rand(npts, 3) - 0.5)
    def assign_symmetry(self, state):
        Ptrans = self.tddftb.TransitionDensityMatrix(state)
        # transform test points according to the symmetry operations
        # of the group
        xs = []  
        npts = len(self.test_pts)
        # Put all test points into one long list so that we
        # can compute the values of the transition density
        # for all of the in one step
        for x in self.test_pts:
            # image of x under all symmetry transformations
            # e.g.  xts = {Id(x), sigma(x), ....}
            xts = self.symmetry_group.symmetry_equivalents(x)
            xs += xts 
        xs = np.array(xs)
        # The coordinates of the transformed vectors belonging
        # to the the i-th test point can be indexed as 
        # xs[i*ng:(i+1)*ng,:]

        # size of group - number of symmetry elements
        ng = self.symmetry_group.size()
        # evaluate electron density on all transformed points
        fs = electron_density(xs.transpose(), self.bs.bfs, Ptrans)
        """
        def f(pos):
            # compute the transition density at the points pos
            x,y,z = pos
            fx = electron_density([np.array([x]),np.array([y]),np.array([z])], \
                                      self.bs.bfs, Ptrans)
            return fx
        """
        for i in range(0, npts):
            xts = xs[i*ng:(i+1)*ng,:]
            fxs = fs[i*ng:(i+1)*ng]
            irrep = self.symmetry_group.assign_irrep(xts, fxs)
            if irrep != "?":
                # an irrep could be assigned 
                break
        return irrep


class SymmetryAssignmentEx_new:
    """
    assign symmetry labels to excited states
    """
    def __init__(self, tddftb, symmetry_group):
        self.tddftb = tddftb
        self.dftb2 = tddftb.dftb2
        if self.dftb2 != None:
            self.bs = AtomicBasisSet(self.dftb2.atomlist)
        self.symmetry_group = symmetry_group
        # How do I choose the test points?
        npts = 5
        # random points in a box [-3,3]x[-3,3]x[-3,3] around the first atom
        pos_at1 = np.array(self.dftb2.atomlist[0][1]) # position of 1st atom
        test_pts1 = pos_at1 + 3.0 * 2.0*(np.random.rand(npts, 3) - 0.5)
        pos_atN = np.array(self.dftb2.atomlist[-1][1]) # position of last atom
        test_ptsN = pos_atN + 3.0 * 2.0*(np.random.rand(npts, 3) - 0.5)
        self.test_pts = np.vstack((test_pts1,test_ptsN))
        # transform test points according to the symmetry operations
        # of the group
        xs = []  
        self.npts = len(self.test_pts)
        # Put all test points into one long list so that we
        # can compute the values of the transition density
        # for all of the in one step
        for x in self.test_pts:
            # image of x under all symmetry transformations
            # e.g.  xts = {Id(x), sigma(x), ....}
            xts = self.symmetry_group.symmetry_equivalents(x)
            xs += xts 
        self.xs = np.array(xs)
        # The coordinates of the transformed vectors belonging
        # to the the i-th test point can be indexed as 
        # xs[i*ng:(i+1)*ng,:]

        # size of group - number of symmetry elements
        self.ng = self.symmetry_group.size()
        # number of basis functions
        self.nbf = len(self.bs.bfs)
        # evaluate basis functions on the grid of symmetry equivalent points
        grid = self.xs.transpose()
        self.bf_grid = []
        for bf in self.bs.bfs:
            self.bf_grid.append( bf.amp(*grid) )
        self.bf_grid = np.array(self.bf_grid)
        ## classify molecular orbitals by symmetry
        #self.orbs_irreps = self.assign_orbital_symmetries(self.dftb2.orbs)
        # 
        self.selected_irreps = ["B1U"]
        #
    def transition_densities(self, Ptrans):
        """
        Parameters:
        ===========
        Ptrans[I,a,b] is the a,b component (a,b AOs) of the transition density matrix belonging to state I

        Returns:
        ========
        rho: transition densities evaluated at the symmetry equivalent points
          rho[I,:] are the values for state I
        """

        nst,nao,nao = Ptrans.shape
        rho = np.zeros((nst, len(self.xs)))
        for m in range(0, self.nbf):
            for n in range(0, self.nbf):
                rho += Ptrans[:,m,n] * self.bf_grid[m].conjugate() * self.bf_grid[n] # NxN grid*grid multiplications

        """
        # fast numpy code does the same
        rho = np.sum(np.dot(Ptrans, self.bf_grid) * self.bf_grid, axis=0)
        """
        return rho
    def assign_symmetry(self, state):
        PtransI = self.tddftb.TransitionDensityMatrix(state)
        # evaluate electron density on all transformed points
        nao,nao = PtransI.shape
        Ptrans = np.reshape(PtransI, (1,nao,nao))
        fs = self.transition_densities(Ptrans)[0,:]

        for i in range(0, self.npts):
            xts = self.xs[i*self.ng:(i+1)*self.ng,:]
            fxs = fs[i*self.ng:(i+1)*self.ng]
            irrep = self.symmetry_group.assign_irrep(xts, fxs)
            if irrep != "?":
                # an irrep could be assigned 
                break
        return irrep
    def assign_symmetries(self, T):
        time1 = time.time()
        nocc,nvirt,nst = T.shape
        irreps = []
        # active occupied and virtual orbitals
        orbs_occ = self.dftb2.orbs[:,self.tddftb.active_occupied_orbs]
        orbs_virt = self.dftb2.orbs[:,self.tddftb.active_virtual_orbs]

        for n in range(0, nst):
            # transition density matrix
            Ptrans = np.dot(orbs_occ, np.dot(T[:,:,n], orbs_virt.transpose()))
            # transition density evaluated at the symmetry equivalent points
            """
            # very slow python code
            rho = np.zeros(len(self.xs))
            for a in range(0, self.nbf):
                for b in range(0, self.nbf):
                    rho += Ptrans[a,b] * self.bf_grid[a].conjugate() * self.bf_grid[b] # grid*grid multiplications
            #
            """
            """
            # little bit faster, but still slow
            rho = np.zeros(len(self.xs))
            for i in range(0, len(self.xs)):
                rho[i] = np.dot(self.bf_grid[:,i], np.dot(Ptrans, self.bf_grid[:,i]))
            """

            # fast numpy code does the same
            rho = np.sum(np.dot(Ptrans, self.bf_grid) * self.bf_grid, axis=0)

            # try to assign irrep based on the transformation of the transition density
            for i in range(0, self.npts):
                xts = self.xs[i*self.ng:(i+1)*self.ng,:]
                fxs = rho[i*self.ng:(i+1)*self.ng]
                irrep = self.symmetry_group.assign_irrep(xts, fxs)
                if irrep != "?":
                    # an irrep could be assigned 
                    break
            irreps.append( irrep )
        irreps = np.array(irreps) # array of strings
        time2 = time.time()
        print "assigning symmetries took %s seconds" % (time2-time1)
        return irreps
    
    def select_irreps(self, selected_irreps):
        self.selected_irreps = selected_irreps
    def select_vectors(self, Ts):
        """
        select expansion vectors with desired symmetry
        """
        irreps = self.assign_symmetries(Ts)
        # indeces of desired symmetry
        good = []
        # other indeces
        bad = []
        for i,irrep in enumerate(irreps):
            if irrep in self.selected_irreps:
                good.append(i)
            else:
                bad.append(i)
        good = np.array(good)
        bad = np.array(bad)
        return good, bad, irreps
    ############ ORBITAL SYMMETRIES #########
    def assign_orbital_symmetries(self, orbs):
        """
        classify molecular orbitals by symmetry
        """
        nao,nmo = orbs.shape
        amp = np.dot(orbs.transpose(), self.bf_grid)
        irreps = []
        for m in range(0, nmo):
            #print "Trying to determine symmetry of orbital %d" % m
            # try to assign irrep based on the transformation of the orbital amplitude
            for i in range(0, self.npts):
                xts = self.xs[i*self.ng:(i+1)*self.ng,:]
                fxs = amp[m, i*self.ng:(i+1)*self.ng]
                irrep = self.symmetry_group.assign_irrep(xts, fxs)
                if irrep != "?":
                    # an irrep could be assigned 
                    break
            irreps.append( irrep )
        irreps = np.array(irreps) # array of strings
        #
        print "Symmetries of Molecular Orbitals"
        print "================================"
        for m in range(0, nmo):
            print "%3d      %.7f         %s" % (m,self.dftb2.orbe[m], irreps[m])
        #
        return irreps

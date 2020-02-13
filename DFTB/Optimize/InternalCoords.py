# -*- coding: utf-8 -*-
"""
redundant internal coordinates as explained in Refs. [1] and [4]

References
----------
[1]  P. Pulay, G. Fogerasi,
     "Geometry optimization in redundant internal coordinates",
     J. Chem. Phys. 96(4), 1992
[2]  V. Bakken, T. Helgaker,
     "The efficient optimization of molecular geometries using redundant internal coordinates",
     J. Chem. Phys. 117, 9160, 2002
[3]  J. Baker, A. Kessi, B. Delley
     "The generation and use of delocalized internal coordinates in geometry optimization"
     J. Chem. Phys. 105(1), 1996
[4]  C. Peng, P. Ayala, H. Schlegel, M. Frisch
     "Using Redundant Internal Coordinates to Optimize Equilibrium Geometries and Transition States",
     J. Comput. Chem., Vol. 17, No. 1, 49-56 (1996)
"""

from DFTB.ForceField.PeriodicForceField import PeriodicForceField, read_force_field
from DFTB import XYZ
from DFTB import AtomicData
from DFTB import utils
from DFTB.Modeling import MolecularCoords as MolCo
from DFTB.Analyse import MolecularGraph

import numpy as np
import numpy.linalg as la
import scipy.linalg as sla

class InternalValenceCoords:
    def __init__(self, atomlist, freeze=[], explicit_bonds=[], verbose=0):
        """
        setup system of internal coordinates using
        valence bonds, angles and dihedrals

        Parameters
        ----------
        atomlist   :  list of tuples (Z,[x,y,z]) with molecular
                      geometry, connectivity defines the valence
                      coordinates

        Optional
        --------
        freeze          :  list of tuples of atom indices (starting at 0) corresponding
                           to internal coordinates that should be frozen
        explicit_bonds :   list of pairs of atom indices (starting at 0) between which artificial
                           bonds should be inserted, i.e. [(0,1), (10,20)]. 
                           This allows to connect separate fragments.
        verbose        :   write out additional information if > 0
        """
        self.verbose=verbose
        self.atomlist = atomlist
        
        self.masses = AtomicData.atomlist2masses(self.atomlist)
        # Bonds, angles and torsions are constructed by the force field.
        # Atom types, partial charges and lattice vectors
        # all don't matter, so we assign atom type 6 (C_R, carbon in resonance)
        # to all atoms.
        atomtypes = [6 for atom in atomlist]
        partial_charges = [0.0 for atom in atomlist]
        lattice_vectors = [ [0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0],
                            [0.0, 0.0, 0.0] ]
        # But since the covalent radii are wrong, we have to provide
        # the connectivity matrix
        conmat = XYZ.connectivity_matrix(atomlist, hydrogen_bonds=True)

        # insert artificial bonds
        for (I,J) in explicit_bonds:
            print "explicit bond between atoms %d-%d" % (I+1,J+1)
            conmat[I,J] = 1

        # Internal coordinates only work if the molecule does not
        # contain disconnected fragments, since there is no way how the
        # interfragment distance could be expressed in terms of internal coordinates.
        # We need to check that there is only a single fragment.
        fragment_graphs = MolecularGraph.atomlist2graph(self.atomlist, conmat=conmat)
        nr_fragments = len(fragment_graphs)
        error_msg  = "The molecule consists of %d disconnected fragments.\n" % nr_fragments
        error_msg += "Internal coordinates only work if all atoms in the molecular graph are connected.\n"
        error_msg += "Disconnected fragments may be joined via an artificial bond using the\n"
        error_msg += "`explicit_bonds` option.\n"
        assert nr_fragments == 1, error_msg
     
        # Frozen degrees of freedom do not necessarily correspond to physical bonds
        # or angles. For instance we can freeze the H-H distance in water although there
        # is no bond between the hydrogens. To allow the definition of such 'unphysical'
        # internal coordinates, we have to modify the connectivity matrix and introduce
        # artificial bonds.
        for IJKL in freeze:
            if len(IJKL) == 2:
                I,J = IJKL
                # create artificial bond between atoms I and J
                conmat[I,J] = 1
            elif len(IJKL) == 3:
                I,J,K = IJKL
                # create artifical bonds I-J and J-K so that the valence angle I-J-K exists
                conmat[I,J] = 1
                conmat[J,K] = 1
            elif len(IJKL) == 4:
                I,J,K,L = IJKL
                # create artifical bonds I-J, J-K and K-L so that the dihedral angle I-J-K-L
                # exists
                conmat[I,J] = 1
                conmat[J,K] = 1
                conmat[K,L] = 1

        # cutoff for small singular values when solving the 
        # linear system of equations B.dx = dq in a least square
        # sense.
        self.cond_threshold = 1.0e-10
                
        self.force_field = PeriodicForceField(atomlist, atomtypes, partial_charges, lattice_vectors, [],
                                              connectivity_matrix=conmat, verbose=1)        
        x0 = XYZ.atomlist2vector(atomlist)
        # shift molecule to center of mass
        self.x0 = MolCo.shift_to_com(x0, self.masses)

        self._selectActiveInternals(freeze=freeze)
        
    def _selectActiveInternals(self, freeze=[], use_internal_types=['B','A','D']):
        """
        Here a redundant set of active internal coordinates is selected. The
        coordinates are included if they are of one of the types specified
        in `use_internal_types` (B - bond, A - bending angle, D - dihedral, I - inversion)

        With 
            use_internal_types=['B','A','D']
        all internal coordinates are used except to inversions. 
        """
        nat = len(self.atomlist)

        # compute internal coordinates
        #   qprim - list all internal coordinates (bonds, angles, dihedrals, inversions) in bohr or radians
        #   Bprim - Wilson's B-matrix, derivatives of all internal coordinates w/r/t cartesian coordinates
        #         Bprim[i,j] = dqprim(i)/dx(j)    i = 1,...,nred    j=1,...,3*Nat
        qprim, Bprim = self.force_field.getRedundantInternalCoordinates(self.x0, 0)
        # number of redundant internal coordinates
        nred = len(qprim)

        # Perform singular value decomposition of B.B^T. There should be exactly 3*Nat-6 or 3*Nat-5
        # singular values > 0, and 6 or 5 singular values = 0
        G = np.dot(Bprim, Bprim.transpose())
        U,s,Vh = la.svd(G)

        nint = len(s[s > 1.0e-8])
        error_msg = "There should be 3*Nat-6= %d or 3*Nat-5= %d non-zero eigenvalues of matrix G = B.B^T, but got %d ones!" % (3*nat-6, 3*nat-5, nint)
        #assert (nint == 3*nat-6) or (nint == 3*nat-5), error_msg

        self.coord_types, self.atom_indices = self.force_field.getInternalCoordinateDefinitions()
        
        # Here we select a subset of the internal coordinates
        self.active_internals = []
        for k,ktype in enumerate(self.coord_types):
            # select internal coordinate k if it has the right type
            if not (ktype in use_internal_types):
                continue
            # Internal coordinate k has passed tests, so we add it to the active set.
            self.active_internals.append( k )
                
                
        # Indices into set of active coordinate which are frozen. The
        # gradient along these coordinates is set to 0 when transforming
        # the cartesian gradient to internal coordinates.
        self.frozen_internals = []
        for IJKL in freeze:
            self.freeze(IJKL)
                
        # In the following we prepare the columns and row labels for printing
        # internal coordinates nicely
        
        # labels for internal coordinates and cartesian coordinates
        self.internal_labels = []
        max_width = 0
        for iact,k in enumerate(self.active_internals):
            # shift indices by 1
            atoms = [I+1 for I in self.atom_indices[k]]
            label = self.coord_types[k] + "(" + "-".join(map(lambda i: str(i), atoms)) + ")"
            self.internal_labels.append(label)
            max_width = max(max_width, len(label))
        # all labels should have the same width
        for i in range(0, len(self.internal_labels)):
            self.internal_labels[i] = self.internal_labels[i].ljust(max_width, " ")

        # labels for cartesian positions, e.g. C1x, C1y, C1z, ...
        self.cartesian_labels = []
        for i,(Z,pos) in enumerate(self.atomlist):
            at = AtomicData.atom_names[Z-1].upper()
            xyz_labels = [at+str(i+1)+"x",at+str(i+1)+"y",at+str(i+1)+"z"]
            self.cartesian_labels += xyz_labels
        
        # values of active coordinates and Wilson's B-matrix restricted to active set
        q = qprim[self.active_internals]
        B = Bprim[self.active_internals,:]

        if self.verbose > 0:
            # table with internal coordinates, atom indices and current values
            print self._table_internal_coords(q)
            print self._table_Bmatrix(B)

    def _internal_from_indices(self, IJKL):
        """
        find the internal coordinate corresponding to the atom
        indices IJKL.
        """
        for k,ktype in enumerate(self.coord_types):
            # compare atom indices in any order
            if set(self.atom_indices[k]) == set(IJKL):
                break
        else:
            raise ValueError("The tuple of atom indices %s does not correspond to an existing bond, angle or dihedral!" % [I+1 for I in IJKL])

        if self.verbose > 1:
            print "Atom indices %s correspond to the %d-th internal coordinate, which is of type %s." % (IJKL, k+1, ktype)
            
        return k
            
    def freeze(self, IJKL):
        """
        freezes the internal coordinate defined by the atom indices
        IJKL. The atom indices do not necessarily have to correspond to a "physical" bond, angle or dihedral
        that is actually present in the molecule. So, for instance, you can also freeze the 
        distance between two atoms that are not bonded. 

        Parameters
        ----------
        IJKL:      tuple of 2, 3 or 4 atom indices (starting at 0)
                   (I,J)     -   bond between atoms I and J
                   (I,J,K)   -   valence angle I-J-K
                   (I,J,K,L) -   dihedral angle between the bonds I-J, J-K and K-L
        """
        k = self._internal_from_indices(IJKL)

        if not k in self.active_internals:
            print "Internal coordinate %d does not belong to the active set"
            return

        kact = self.active_internals.index(k)
        if not kact in self.frozen_internals:
            self.frozen_internals.append(kact)

    def coordinate_value(self, x, IJKL):
        """
        current value of internal coordinate defined by atom indices IJKL

        Parameters
        ----------
        x           : vector of cartesian coordinates
        IJKL        : tuple of 2,3 or 4 atom indices (starting from 0)

        Returns
        -------
        val         : value of selected internal coordinate,
                      bond lengths in Angstrom, angles in degrees
        """
        k = self._internal_from_indices(IJKL)
        # determine internal coordinates q ~ x
        qprim, Bprim = self.force_field.getRedundantInternalCoordinates(x, 0)
        val = qprim[k]
        
        # convert units
        if len(IJKL) == 2:
            # bond
            val *= AtomicData.bohr_to_angs
        elif len(IJKL) in [3,4]:
            # angles
            val *= 180.0/np.pi

        return val
        
    def cartesian2internal(self, x):
        """ 
        convert cartesian to redundant internal coordinates 
        
        Parameters
        ----------
        x          : cartesian coordinates, array of length 3*Nat

        Returns
        -------
        q          : redundant internal coordinates of length n >= 3*nat-6
        """
        if self.verbose > 1:
            print "cartesian -> internal"
        x = MolCo.shift_to_com(x, self.masses)

        # values of internal coordinates at cartesian position x
        qprim, Bprim = self.force_field.getRedundantInternalCoordinates(x, 0)

        #test_wilson_bmatrix(self.force_field, x)
        
        # 
        q = qprim[self.active_internals]
        B = Bprim[self.active_internals,:]

        # At the cartesian position x0 we know the internal coordinates q0 exactly.
        # We need this information later to deduce the new cartesian coordinates x (close to x0)
        # which corresponds to the new internal coordinates q (also close to q0).
        self.q0 = q
        self.x0 = x
        self.B0 = B
        # projector 
        self.P0 = self.feasibility_projector(B)
        
        return q

    def internal2cartesian(self, q, max_iter=1000):
        """
        transform internal coordinates back to cartesians.
        
        Since the internal coordinates are curvilinear the transformation
        has to be done iteratively and depends on having a closeby point q0
        for which we know the cartesian coordinates x0. If the displacement
          dq = q-q0 
        is too large, the iteration will not converge.

        Given the initial point
          x0 ~ q0
        we wish to find the cartesian coordinate x that corresponds to q
          x ~ q      q = q0 + dq

        Parameters
        ----------
        q          :  redundant internal coordinates, should not be too far
                      from the coordinates obtained by the last
                      call to `cartesian2delocalized(...)`

        Returns
        -------
        x          :  cartesian coordinates corresponding to q,
                      x ~ q

        Optional
        --------
        max_iter   :  maximum number of iterative refinements
                      of cartesian position
        """
        # If the geometry did not change, there is nothing to do,
        # just return the old cartesian coordinates.
        if la.norm(q-self.q0) < 1.0e-15:
            return self.x0

        if self.verbose > 0:
            print "internal -> cartesian"
            print " The internal coordinates are iteratively transformed"
            print " transformed back to cartesian coordinates."

        # bending and torsion angles should be in the range [0,2*pi] and
        # inversions in the range [-pi/2,pi/2]
        q = self._wrap_angles(q)

        if self.verbose > 1:
            print "previous internal coordinates q0"
            print self._table_internal_coords(self.q0)
            print "current internal coordinates q"
            print self._table_internal_coords(q)

        # By how much did the internal coordinates change relative
        # to the reference geometry, where we know the transformation
        # x0 ~ q0 ?
        dq = q-self.q0
        # solve  B(q0).dx = dq    for qx
        ret = sla.lstsq(self.B0, dq, cond=self.cond_threshold)
        dx = ret[0]
        # Initial quess for updated cartesian coordinates
        xi = self.x0+dx

        if self.verbose > 2:
            err_dx = la.norm( np.dot(self.B0, dx) - dq )
            print "error |B0.dx - (q-q0)|= %e        |dq|= %e" % (err_dx, la.norm(dq))

        if self.verbose > 0:
            print "  Iteration= %4.1d    |dx|= %e   |q-q0|= %e" % (0, la.norm(dx), la.norm(dq))

        # The initial guess xi is refined iteratively        
        for i in range(0, max_iter):
            qi_prim, Bi_prim = self.force_field.getRedundantInternalCoordinates(xi, 0)
            qi = qi_prim[self.active_internals]
            Bi = Bi_prim[self.active_internals,:]

            if self.verbose > 2:
                print "intermediate internal coordinates q(%d)" % (i+1)
                print self._table_internal_coords(qi)
                print self._table_Bmatrix(Bi)

            ddq = q-qi
            # solve B(xi).dx = q-qi for refinement dx
            ret = sla.lstsq(Bi, ddq, cond=self.cond_threshold)
            dx = ret[0]
            xi += dx

            # error of least square solution
            err_dx = la.norm( np.dot(Bi, dx) - ddq )

            if self.verbose > 2:
                print "changes in internal coordinates q-q(%d)" % (i+1)
                print self._table_internal_coords(ddq)
            if self.verbose > 0:
                print "  Iteration= %4.1d   |B.dx-dq|= %e   |dx|= %e   |q-q(%d)|= %e" % (i+1, err_dx, la.norm(dx), i+1, la.norm(ddq))

            # If qi has converged to q or the desired accuracy cannot
            # be reached, because the solution of B.dx = q-qi has a
            # large error anyway, we accept the final cartesian coordinates xi.
            # The factor 1.01 is to avoid an endless loop if |dq| == err_dx.
            if la.norm(ddq) < 1.01*max(err_dx, 1.0e-10):
                break
            
            # Abort if new cartesian coordinates make no sense.
            dx_norm = la.norm(dx)
            if (dx_norm > 1.0e6):
                raise RuntimeError("change in cartesian coordinates too large |dx|= %e !" % dx_norm)

        else:
            dq_norm = la.norm(q-qi)
            if dq_norm < 0.5:
                print "WARNING: Internal->cartesian transformation did not converge!"
                print "         But |dq|= %e is not too large, so we try to continue anyway." % la.norm(q-qi)
            else:
                raise NotConvergedError("ERROR: internal->cartesian transformation did not converge! |q-qi|= %e" % dq_norm)

        # We have succesfully found the cartesian coordinates x ~ q
        # and use them as the new reference geometry.
        self.q0 = np.copy(qi)
        self.x0 = np.copy(xi)
        self.B0 = np.copy(Bi)
        
        return xi

    def transform_gradient(self, x, g_cart, max_iter=200):
        """
        transform cartesian gradient g_cart = dE/dx to redundant internal coordinates according to

              B^T g_intern = g_cart

        The underdetermined system of linear equations is solved
        in a least square sense.
        """
        if self.verbose > 0:
            print "transform gradient to internal coordinates"
            print "  dE/dx -> dE/dq"
        # Since the energy of a molecule depends only on the internal
        # coordinates q it should be possible to transform the cartesian
        # gradient exactly into internal coordinates according to
        #   dE/dx(i) = sum_j dE/dq(j) dq(j)/dx(i)
        #            = sum_j (B^T)_{i,j} dE/dq(j)
        #
        # cartesian force
        f_cart = -g_cart
        # cartesian accelerations
        a_cart = f_cart / self.masses
        # cartesian velocity
        dt = 1.0
        vel = a_cart * dt
        # shift position to center of mass
        x = MolCo.shift_to_com(x, self.masses)
        # eliminate rotation and translation from a_cart
        I = MolCo.inertial_tensor(self.masses, x)
        L = MolCo.angular_momentum(self.masses, x, vel)
        P = MolCo.linear_momentum(self.masses, vel)
        omega = MolCo.angular_velocity(L,I)
        vel = MolCo.eliminate_translation(self.masses, vel, P)
        vel = MolCo.eliminate_rotation(vel, omega, x)
        # Check that the total linear momentum and angular momentum
        # indeed vanish.
        assert la.norm( MolCo.linear_momentum(self.masses, vel) ) < 1.0e-14
        assert la.norm( MolCo.angular_momentum(self.masses, x, vel) ) < 1.0e-14
        # go back from velocities to gradient
        g_cart = - vel / dt * self.masses 
        
        # solve  B^T . g_intern = g_cart  by linear least squares.
        ret = sla.lstsq(self.B0.transpose(), g_cart,
                        cond=self.cond_threshold)
        g_intern = ret[0]

        if self.verbose > 1:
            print self._table_gradients(g_cart, g_intern)

        # check solution
        err = la.norm( np.dot(self.B0.transpose(), g_intern) - g_cart )
        assert err < 1.0e-10, "|B^T.g_intern - g_cart|= %e" % err
            
        # Abort if gradients are not reasonable
        gnorm = la.norm(g_intern)
        if gnorm > 1.0e5:
            raise RuntimeError("ERROR: Internal gradient is too large |grad|= %e" % gnorm)
        
        # Components of gradient belonging to frozen internal
        # coordinates are zeroed to avoid changing them.
        g_intern[self.frozen_internals] = 0.0

        # Simply setting some components of the gradient to zero
        # will most likely lead to inconsistencies if the coordinates
        # are coupled (Maybe it's impossible to change internal coordinate X
        # without changing coordinate Y at the same time). These
        # inconsistencies are removed by applying the projection
        # operator repeatedly until the components of the frozen
        # internal coordinates have converged to 0.
        if self.verbose > 0:
            print " apply projector P=G.G- to internal gradient"

        # The projector is updated everytime cartesian2internal(...)
        # is called.
        proj = self.P0
        
        for i in range(0, max_iter):
            g_intern = np.dot(proj, g_intern)
            # gradient along frozen coordinates should be zero
            gnorm_frozen = la.norm(g_intern[self.frozen_internals])
            
            if self.verbose > 0:
                print "  Iteration= %4.1d   |grad(frozen)|= %s" % (i, gnorm_frozen)
                
            if gnorm_frozen < 1.0e-10:
                break
            else:
                g_intern[self.frozen_internals] = 0.0
        else:
            if gnorm_frozen < 1.0e-5:
                print "WARNING: Projection of gradient vector in internal coordinates did not converge!"
                print "         But |grad(frozen)|= %e is not too large, so let's continue anyway." % gnorm_frozen
            else:
                raise RuntimeError("ERROR: Projection of gradient vector in internal coordinates did not converge! |grad(frozen)|= %e" % gnorm_frozen)
            
        if self.verbose > 1:
            print "gradients after applying projector (only internal gradient changes)"
            print self._table_gradients(g_cart, g_intern)
        
        return g_intern

    def feasibility_projector(self, B):
        """
        compute the projector P into the space of feasible displacement
        vectors 

        Not all displacements in internal coordinates are feasible,
        since internal coordinates often depend on each other so that changing
        one coordinate requires changing another coordinate at the same time.
        A displacement that would change only one of them, would be impossible
        to realize. Given a displacement vector `g`, the projection matrix
        P resolves these dependencies, so that `g_prj = P.g` is a feasible displacement
        vector.

        Parameters
        ----------
        B          : Wilson's B-matrix for active set of 
                     internal coordinates

        Returns
        -------
        P          : projector onto feasible directions
        """
        # n is the number of active internal coordinates
        n,ncart = B.shape
        assert n == len(self.active_internals)
        # G = B.B^T
        G = np.dot(B, B.transpose())
        # compute pseudo-inverse G^-
        Gm = la.pinv(G)
        # projector due to coupling of internal coordinates
        P = np.dot(G,Gm)
        
        return P
    
    def internal_step(self, x0, IJKL, incr, max_iter=200):
        """
        take a step of size `dq` along the internal coordinate defined
        by the atom indices IJKL starting from the cartesian coordinates `x0`.

        Parameters
        ----------
        x0      :  initial cartesian coordinates
        IJKL    :  tuple of 2, 3 or 4 atom indices (starting at 0)
                   (I,J)     -   bond between atoms I and J
                   (I,J,K)   -   valence angle I-J-K
                   (I,J,K,L) -   dihedral angle between the bonds I-J, J-K and K-L
        incr    :  displacement, in bohr for bond lengths and in radians for angles

        Optional
        --------
        max_iter : maximum number of interators for projecting
                   the step into the feasible subspace

        Returns
        -------
        x1      :  cartesian coordinates corresponding to the displaced
                   geometry x0 ~ q0, x1 ~ q0+incr*e_IJKL
        """
        if self.verbose > 0:
            print "take internal step along coordinate %s" % map(lambda I: I+1, IJKL)
        # Which internal coordinate should be changed?
        k = self._internal_from_indices(IJKL)
        # 1) transform cartesian to internal coordinates, x0 ~ q0
        q0 = self.cartesian2internal(x0)
    
        # 2) take a step along the internal coordinate k
        dq = np.zeros(len(q0))
        dq[k] = incr

        # 3) Since internal coordinates are coupled it's not
        # always possible to change only one coordinate without
        # changing others. Therefore the displacement dq has to be
        # projected into the space of displacements that are compatible
        # with the dependencies between internal coordinates.
        # For instance the angles A-B-C, A-B-D and C-B-D are coupled. 
        # If we increase the angle A-B-C by 2 degrees we have
        # to decrease the angles A-B-D and C-B-D each by 1 degree. This is
        # accomplished by the projection.
        
        if self.verbose > 0:
            print " apply projector P=G.G- to step dq in internal coordinates"

        # projector onto feasible displacements
        proj = self.P0
        
        for i in range(0, max_iter):
            dq = np.dot(proj, dq)
            ddq = dq[k] - incr
            if self.verbose > 0:
                print "  Iteration= %4.1d   |dq[k]-incr|= %s" % (i, abs(ddq))
            if abs(ddq) < 1.0e-10:
                break
            else:
                dq[k] = incr
        else:
            if abs(ddq) < 1.0e-6:
                print "WARNING: Projection of displacement vector in internal coordinates did not converge! |dq[k]-incr|= %e" % abs(ddq)
                print "         But the deviation is not too large, so let's try to continue anyway."
            else:
                raise RuntimeError("ERROR: Projection of displacement vector in internal coordinates did not converge! |dq[k]-incr|= %e" % abs(ddq))

        q1 = q0 + dq

        if self.verbose > 0:
            print "initial coordinates"
            print self._table_internal_coords(q0)
            print "coordinates after step along %d-th coordinate" % (k+1)
            print self._table_internal_coords(q1)
        
        # 4) transform back x1 ~ q1
        x1 = self.internal2cartesian(q1)
        
        # verify that x really corresponds to internal coordinate q
        q_test = self.cartesian2internal(x1)
        err = la.norm(q1 - q_test)
#        assert err < 1.0e-10, "|q(x(q)) - q|= %e" % err

        if self.verbose > 0:
            print "coordinates after step (determined from cartesians)"
            print self._table_internal_coords(q_test)
        
        return x1
        
    def _wrap_angles(self, q):
        """
        Bending angles and dihedral angles have to be in the
        range [0,pi], while inversion angles have to be in the range [-pi/2,pi/2].
        Angles outside these ranges are wrapped back to the equivalent
        angle inside the range.

        Parameters
        ----------
        q          :  values of active internal coordinates with angles outside
                      the valid ranges

        Returns
        -------
        q_wrap     :  values of internal coordinates with all angles inside valid ranges
        """
        # types of the active coordinates
        act_coord_types = self.coord_types[self.active_internals]
        # select internal coordinates which are angles
        bending_indices   = np.where(act_coord_types == 'A')[0]
        torsion_indices   = np.where(act_coord_types == 'D')[0]
        inversion_indices = np.where(act_coord_types == 'I')[0]

        bending   = q[bending_indices]
        torsion   = q[torsion_indices]
        inversion = q[inversion_indices]

        # Shifting angles by multiples of 2*pi has no effect
        # for any type of angle.
        bending = bending % (2*np.pi)
        torsion = torsion % (2*np.pi)
        inversion = inversion % (2*np.pi)

        # wrap bending and torsion angles back into the range [0,pi]
        bending[bending > np.pi] = 2*np.pi - bending[bending > np.pi]
        torsion[torsion > np.pi] = 2*np.pi - torsion[torsion > np.pi]
        
        # wrap inversion angle back into the range [-pi/2,pi/2]
        # shift interval [pi,2*pi] to [-pi,0]
        inversion[inversion > np.pi] = inversion[inversion > np.pi] - 2*np.pi
        # wrap [-pi,-pi/2] to [-pi/2,0]
        # and  [pi/2,pi] to [0, pi/2]
        inversion[inversion > np.pi/2] = np.pi - inversion[inversion > np.pi/2]
        inversion[inversion < -np.pi/2] = -np.pi - inversion[inversion < -np.pi/2]

        # update internal coordinates with wrapped angles
        q[bending_indices]   = bending
        q[torsion_indices]   = torsion
        q[inversion_indices] = inversion
        
        return q
        
    def _table_internal_coords(self, q):
        """
        make table with internal coordinates, atom indices and current values
        """
        txt  = "\n"
        txt += " ============================================================================= \n"
        txt += "                 Selected Active Internal Coordinates                         \n"
        txt += "                                                                              \n"
        txt += " Internal Coordinate           Atom Indices         Current Value      Frozen \n"
        txt += "    Nr.    Type              (starting at 1)                                  \n"
        txt += " ============================================================================= \n"
        for iact,k in enumerate(self.active_internals):
            # If this coordinate frozen?
            if iact in self.frozen_internals:
                frozen_flag = "Y"
            else:
                frozen_flag = "N"
                
            ktype = self.coord_types[k]
            if ktype == 'B':
                I,J = self.atom_indices[k]
                txt += "  %4.1d     BOND            %3.1d  %3.1d                   %+8.5f Ang    " % (iact+1,I+1,J+1, q[iact]*AtomicData.bohr_to_angs)
            elif ktype == 'A':
                I,J,K = self.atom_indices[k]
                txt += "  %4.1d     ANGLE           %3.1d  %3.1d  %3.1d          %+12.5f degrees" % (iact+1,I+1,J+1,K+1, q[iact]*180.0/np.pi)
            elif ktype == 'D':
                I,J,K,L = self.atom_indices[k]
                txt += "  %4.1d     DIHEDRAL        %3.1d  %3.1d  %3.1d  %3.1d     %+12.5f degrees" % (iact+1,I+1,J+1,K+1,L+1, q[iact]*180.0/np.pi)
            elif ktype == 'I':
                I,J,K,L = self.atom_indices[k]
                txt += "  %4.1d     INVERSION       %3.1d  %3.1d  %3.1d  %3.1d     %+12.5f degrees" % (iact+1,I+1,J+1,K+1,L+1, q[iact]*180.0/np.pi)
            txt += "    %s \n" % frozen_flag
        txt += "\n"
        
        return txt

    def _table_gradients(self, g_cart, g_intern):
        """
        make tables with gradient vectors in cartesian and internal coordinates
        """
        txt =  ""
        txt += " =========================== \n"
        txt += " Cartesian Gradient dE/dx(i) \n"
        txt += " =========================== \n"
        txt += utils.annotated_matrix(np.reshape(g_cart, (len(g_cart),1)),
                                      self.cartesian_labels, ["grad E"])
        txt += " cartesian gradient norm : %e \n" % la.norm(g_cart)
        txt += "\n"

        txt += " ========================== \n"
        txt += " Internal Gradient dE/dq(i)  \n"
        txt += " ========================== \n"
        txt += utils.annotated_matrix(np.reshape(g_intern, (len(g_intern),1)),
                                      self.internal_labels, ["grad E"])
        txt += " internal gradient norm : %e \n" % la.norm(g_intern)
        txt += "\n"

        return txt
        
        
    def _table_Bmatrix(self, B):
        """
        make table with rows of Wilson's B-matrix belonging to the
        active internal coordinates
        """

        txt  = " ======================================================================\n"
        txt += " Wilson's B-matrix B(i,j) = dq(i)/dx(j) for active internal coordinates\n"
        txt += "     B - bond,   A - valence angle,   D - dihedral,   I - inversion    \n"
        txt += " ======================================================================\n"
 
        txt += utils.annotated_matrix(B, self.internal_labels, self.cartesian_labels)
        txt += "\n"

        return txt
        
# This exception is thrown if the iterative transformation from
# internal to cartesian coordinates fails to converge
class NotConvergedError(Exception):
    pass


############## Testing ###########################

def test_wilson_bmatrix(force_field, x0):
    """
    check that Wilson B matrix obtained from numerical differentiation
    of internal coordinates w/r/t cartesian coordinates gives the same
    result as the analytical Fortran code
    """
    # redundant internal coordinates
    q0, B = force_field.getRedundantInternalCoordinates(x0, 0)
    nred = len(q0)
    # compare B-matrix B_ij = dqi/dxj from numerical differentiation
    # with analytical one
    for i in range(0, nred):
        def f(x):
            q0i,Bi = force_field.getRedundantInternalCoordinates(x,0)
            return q0i[i]
        dqi_dx = utils.numerical_gradient(f,x0)
        err = la.norm(B[i,:] - dqi_dx)
        print "Internal coordinate %d" % i
        print " |dqi/dx (numerical) - dqi/dx (analytical)|= %e" % err
        assert err < 1.0e-8
    

if __name__ == "__main__":
    import sys
    
    atomlist0 = XYZ.read_xyz(sys.argv[1])[-1]
    x0 = XYZ.atomlist2vector(atomlist0)
    
    IC = InternalValenceCoords(atomlist0, verbose=3)

    test_wilson_bmatrix(IC.force_field, x0)
    
    # 1) transform cartesian to internal coordinates, x0 ~ q0
    q0 = IC.cartesian2internal(x0)

    # 2) Take a step along the last internal coordinate,
    # this will only work if the coordinate is not coupled
    # to any other.
    dq = np.zeros(len(q0))
    dq[-1] = 0.1
    q = q0+dq
    
    # 3) and transform back x ~ q
    x = IC.internal2cartesian(q)

    # veriy that really x corresponds to internal coordinate q
    q_test = IC.cartesian2internal(x)
    err = la.norm(q - q_test)
#    assert err < 1.0e-10, "|q(x(q)) - q|= %e" % err
    
    # 4) save initial and displaced geometries
    atomlist = XYZ.vector2atomlist(x, atomlist0)
    XYZ.write_xyz("/tmp/test_internal_coords.xyz", [atomlist0, atomlist])
    

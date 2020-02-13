"""
implicit solvent model

For the electrostatic embedding COSMO (see [1]) is implemented.

References
----------
[1] A.Klamt, G. Schuurmann,
   "COSMO: A New Approach to Dielectric Screening in Solvents with Explicit
    Expressions for the Screening Energy and its Gradient"
    Journal of the Chemical Society, Perkin Transactions 2 5 (1993): 799-805.
"""
from DFTB import XYZ, AtomicData
from DFTB.MolecularIntegrals.LebedevQuadrature import get_lebedev_grid, Lebedev_L2max, Lebedev_Npts

from DFTB.Timer import GlobalTimer as T

# faster Fortran implementation for constructing gamma_solvent
# and its gradients
from DFTB.extensions import cosmo

import numpy as np
import numpy.linalg as la

class SolventCavity:
    def __init__(self, implicit_solvent=0, solvent_radius=1.4,
                 permittivity=80.0, delta_sc=0.9, points_per_sphere=50,
                 trim_solvent="1"):
        """
        Solvation.implicit_solvent: turn the implicit solvent model on (1) or off (0), COSMO is used for the electrostatic screening
        Solvation.solvent_radius: radius (in Ang) of solvent probe (e.g. 1.4 for water)
        Solvation.permittivity: dielectric constant eps of solvent (e.g. 80.0 for water), the screening is computed for a perfect conductor and is then scaled by the function f(eps) = (eps-1)/(eps+x)
        Solvation.delta_sc: separation between charge on solvent molecule and its center (in Ang)
        Solvation.points_per_sphere: number of points on each vdW sphere that make up the solvent accessible surface (SAS), since points inside the solute volume are removed, the number of points exposed to the solvent is smaller
        Solvation.trim_solvent: Remove surface points that do not satisfy a selection condition cond(x,y,z) given in the form of the body of a lambda expression, `lambda x,y,z: cond(x,y,z)`. Example: To model a molecule at the interface (perpendicular to the z-axis) between vacuum (z >= 0) and a liquid (z < 0), only surface points in liquid should be retained. trim_solvent="z < 0.0" includes only those surface points that lie in the lower half-space. By default all points are included.
        """
        self.implicit_solvent = implicit_solvent
        self.solvent_radius = solvent_radius / AtomicData.bohr_to_angs
        self.permittivity = permittivity
        self.delta_sc = delta_sc / AtomicData.bohr_to_angs
        assert 0.0 <= self.delta_sc <= self.solvent_radius
        self.points_per_sphere = points_per_sphere
        # find Lebedev grid with the number of requested points
        ileb = np.argmin( abs(np.array(Lebedev_Npts) - points_per_sphere) )
        # order of request grid
        self.lebedev_order = Lebedev_L2max[ileb]

        # scaling of the dielect screening energy as a function of the permittivity
        # epsilon relative to the conductor-like case
        x = 0.5 # parameter x should be in the range [0,2]
        self.f = (self.permittivity - 1.0)/(self.permittivity + x)

        # selection condition for including surface points
        self.selection_condition = eval("lambda x,y,z: %s" % trim_solvent)
        
    @T.timer
    def constructSAS(self, atomlist):
        """
        construct solvent accessible surface (SAS)
        """
        self.atomlist = atomlist
        # number of atoms
        nats = len(atomlist)
        # coordinates of centers of spheres
        centers = []
        # radii of the spheres
        radii = []
        # On each atom a sphere is placed and points are distributed
        # evenly on the surface of each sphere.
        points = []
        # Each sphere is subdivided into little patches around each
        # point. The weight of a point is the area of the patch
        # its sits on.
        areas = []
        # indeces of spheres to which each point belongs
        parent_spheres = []
        # The points are distributed on each sphere on a Lebedev grid.
        # The number of points is determined by the `order` of the grid. 
        ths,phs,ws = get_lebedev_grid(self.lebedev_order)
        # convert from spherical coordinates on the unit sphere
        # to cartesian coordinates
        xs = np.sin(ths)*np.cos(phs)
        ys = np.sin(ths)*np.sin(phs)
        zs = np.cos(ths)
        
        # loop over atoms
        for i,(Zi,pos) in enumerate(atomlist):
            center = np.array(pos)
            centers.append(center)
            # radius of sphere on atom i + solvent radius
            radius = AtomicData.vdw_radii[AtomicData.atom_names[Zi-1]] + self.solvent_radius
            radii.append(radius)
            # put a mesh of points on the i-th sphere
            for k in range(0, len(ws)):
                # x,y and z-coordinates of the k-th point on the i-th sphere
                vec = np.array([xs[k],ys[k],zs[k]])
                point = center + radius * vec
                area = 4.0 * np.pi * radius**2 * ws[k]
                points.append(point)
                areas.append(area)
                # point k belongs to sphere i
                parent_spheres.append(i)
        
        # Select points which lie on the surface. These are those points,
        # which are not contained in any other sphere
        surface_points = []
        surface_areas = []
        # indeces of atoms to which the surface points are attached
        parent_atoms = []
        # total number of points
        npts = len(points)
        for k in range(0, npts):
            i = parent_spheres[k]
            # find distance of points k to other centers
            for j in range(0, nats):
                if i == j:
                    # exclude sphere to which the point belongs
                    continue
                dist = la.norm(points[k]-centers[j])
                if dist < radii[j]:
                    # point k lies both in sphere i and j
                    break
            else:
                # Point k has passed all test, it really lies on the surface

                # Now we contract the remaining points to a sphere of radius
                #   Rvdw + Rsolv - delta_sc
                point = (points[k] - centers[i]) * (radii[i] - self.delta_sc)/radii[i] + centers[i]
                # The area contracts accordingly.
                area = areas[k] * (radii[i] - self.delta_sc)**2 / radii[i]**2

                # Should the point be included based on the selection condition?
                if self.selection_condition(*point):
                    surface_points.append(point)
                    surface_areas.append(area)
                    parent_atoms.append(i)
        
        self.surface_points = np.array(surface_points)
        self.surface_areas = np.array(surface_areas)
        self.centers = np.array(centers)
        self.parent_atoms = np.array(parent_atoms)

    @T.timer
    def constructCOSMO(self):
        """
        construct the additional gamma matrix

                                     -1  T
           g         = -f(eps) * B. A  .B
            solvent 

        due to the conductor like screening model.

        Returns
        -------
        gamma_solvent   :  numpy array with gamma matrix of size Nat x Nat
        """
        # slow python implementation
        #A, B = self._cosmo_ab_matrices()
        # faster Fortran implementation
        A, B = cosmo.cosmo_ab_matrices(self.centers.transpose(),
                                       self.surface_points.transpose(),
                                       self.surface_areas)
        #
        BinvA = np.dot(B, la.inv(A))
        gamma_solvent = -self.f * np.dot(B, BinvA.transpose())

        # save quantities that might be needed later
        self.BinvA = BinvA
        self.gamma_solvent = gamma_solvent
        
        return gamma_solvent

    def _cosmo_ab_matrices(self):
        """
        construct matrices A and B that give the Coulomb interaction
        between surface charges among themselves and the Coulomb
        interaction between charges inside the cavity and the surface charges, respectively

        Returns
        -------
        A       : numpy array of size Npts x Npts
        B       : numpy array of size Nat x Npts
        """
        nats = len(self.atomlist)
        npts = len(self.surface_points)
        # Coulomb interactions between induced charges on surface
        A = np.zeros((npts,npts))
        for k in range(0, npts):
            # diagonal corresponds to self-interaction of charges on the
            # same patch, see eqn. 7b) in Ref. [1]
            A[k,k] = 3.8 / self.surface_areas[k]
            for l in range(k+1,npts):
                # off-diagonal corresponds to Coulomb interaction of
                # induced charges at different patches k != l
                A[k,l] = 1.0/la.norm(self.surface_points[k,:] - self.surface_points[l,:])
                A[l,k] = A[k,l]
                
        # Coulomb interactions between charges inside cavity (of solute) and
        # induced charges on the surface of the cavity
        B = np.zeros((nats ,npts))
        for i in range(0, nats):
            for k in range(0, npts):
                B[i,k] = 1.0/la.norm(self.centers[i,:] - self.surface_points[k,:])
        
        return A,B
                
    @T.timer
    def constructCOSMOgradient(self):
        """
        construct the gradient of the gamma matrix due to screening w/r/t the
        nuclear coordinates

                                            -1 T       -1          T       -1 T              -1 T
         grad g         = -f { (grad B).(B.A  )  + (B.A  ).(grad B)  - (B.A  ) .(grad A).(B.A  )  }
               solvent 

        Returns
        -------
        g        :  gamma matrix due to solvent
        dg       :  dg[3*k+xyz,i,j] is the derivative of the matrix element gamma_solvent_{i,j}
                    with respect to the xyz-th coordinate of atom k
        """
        # slow python implementation
        #grad_gamma_solvent = self._constructCOSMOgradient()
        # Fortran implementation
        grad_gamma_solvent = cosmo.cosmo_gamma_gradient(\
            self.centers.transpose(),
            self.surface_points.transpose(),
            self.parent_atoms,
            self.BinvA, self.f)
        return self.gamma_solvent, grad_gamma_solvent

    def _constructCOSMOgradient(self):
        """
        python implementation
        """
        nats = len(self.atomlist)
        npts = len(self.surface_points)
        # 
        grad_gamma_solvent = np.zeros((3*nats,nats,nats))
        # build gradient for each atom position
        for k in range(0, nats):
            # construct grad_k A
            # The surface points tm and tn belong to atoms i and j, respectively.
            #
            #                     tm - tn
            #   grad   A    =  - --------- (delta_{k,i} - delta_{k,j})
            #       rk  m,n      |tm-tn|^3
            #
            gradkA = np.zeros((3,npts,npts))
            for m in range(0, npts):
                tm = self.surface_points[m,:]
                i = self.parent_atoms[m]
                for n in range(0, npts):
                    tn = self.surface_points[n,:]
                    j = self.parent_atoms[n]
                    if i == j:
                        # surface points belonging to the same atom exert no force
                        # this atom
                        continue
                    vec = - (tm-tn)/la.norm(tm-tn)**3
                    if i == k:
                        gradkA[:,m,n] += vec
                    if j == k:
                        gradkA[:,m,n] -= vec                        
            # construct grad_k B
            # Atom i is at position ri, surface point tm belongs to atom j
            #
            #                    tm - ri
            #   grad   B     =  --------- (delta_{k,i} - delta_{k,j})
            #       rk   i,m    |tm-ri|^3
            #
            gradkB = np.zeros((3,nats,npts))
            for i in range(0, nats):
                ri = self.centers[i,:]
                for m in range(0, npts):
                    tm = self.surface_points[m,:]
                    j = self.parent_atoms[m]
                    vec = (tm-ri)/la.norm(tm-ri)**3
                    if i == j:
                        # surface point m does not exert force onto its parent atom
                        continue
                    if i == k:
                        gradkB[:,i,m] += vec
                    if j == k:
                        gradkB[:,i,m] -= vec

            # assemble gradient of gamma_solvent from gradients gradkA and gradkB
            #
            for xyz in [0,1,2]:
                grad_gamma_solvent[3*k+xyz,:,:] = \
                     np.dot(gradkB[xyz,:,:], self.BinvA.transpose()) \
                    +np.dot(self.BinvA, gradkB[xyz,:,:].transpose()) \
                    -np.dot(self.BinvA, np.dot(gradkA[xyz,:,:], self.BinvA.transpose()))

            
        # scale to finite dielectic constant
        grad_gamma_solvent *= -self.f

        return grad_gamma_solvent
            
    def getGammaSolvent(self):
        # assumes a previous call to 'constructCOSMO()'
        return self.grad_gamma_solvent
    
    def getInducedCharges(self, Q):
        """
        compute the charges induced on the cavity surface by the
        solute point charges Q.

        Parameters
        ----------
        Q     : numpy array with Mulliken charges

        Returns
        -------
        q_ind : numpy array with induced surface charges in the
                same order as the surface points
        """
        q_ind = -np.dot(self.BinvA.transpose(), Q)

        return q_ind

    def getScreeningEnergy(self, Q):
        """
        compute the screening energy for Mulliken charges Q

                     T
          dE = 1/2 Q . g        . Q
                         solvent

        dE is negative because of the minus sign in the definition of g_solvent
        """
        screening_energy = 0.5 * np.dot(Q, np.dot(self.gamma_solvent, Q))
        return screening_energy
    
    
    def getSurfaceArea(self):
        """ 
        surface area in bohr^2 

        Warning: The surface area is only an approximation because
        points in the reentrant surfaces (close to where two vdW spheres intersect)
        were removed.
        """
        return np.sum(self.surface_areas)
    
    def getSurfacePoints(self):
        """
        retrieve points on the solvent accessible surface at which
        the surface charges are induced

        Returns
        -------
        surface_points  :  numpy array with with positions of points,
                           surface_points[i,:] are the cartesian coordinates 
                           of the i-th point
        """
        return self.surface_points
    
#############################################
# Testing
#############################################

def check_analytic_gradients(atomlist0):
    """
    compare analytical gradients for gamma_solvent with numerical ones
    """
    from DFTB import utils
    from DFTB import XYZ
    
    cavity = SolventCavity(points_per_sphere=6)
    nats = len(atomlist0)
    x0 = XYZ.atomlist2vector(atomlist0)

    # ... analytical gradients
    cavity.constructSAS(atomlist0)
    cavity.constructCOSMO()
    gamma_solvent, grad_ana = cavity.constructCOSMOgradient()
    # ... numerical gradients
    grad_num = np.zeros((3*nats,nats,nats))
    for i in range(0, nats):
        for j in range(0, nats):
            # define function for computing gamma_solvent[i,j]
            def f(x):
                atomlist = XYZ.vector2atomlist(x,atomlist0)
                cavity.constructSAS(atomlist)
                gamma_solvent = cavity.constructCOSMO()
                return gamma_solvent[i,j]
            grad_num[:,i,j] = utils.numerical_gradient(f,x0)

    print "analytical gradient"
    print np.round(grad_ana, 6)
    print "numerical gradient"
    print np.round(grad_num, 6)
    
    err = la.norm(grad_ana-grad_num)
    print "|(grad g_solvent)_num - (grad g_solvent)_ana|= %e" % err
    assert err < 1.0e-5
    
            
def check_surface_area(atomlist):
    # check convergence of surface area as a function of the
    # number of points
    print "# order      points        area    "
    print "#          per sphere      Ang^2   "
    print "# ---------------------------------"
    for io,order in enumerate(Lebedev_L2max):
        # number of points per sphere
        nr_points = Lebedev_Npts[io]
        cavity = SolventCavity(points_per_sphere=nr_points)
        cavity.constructSAS(atomlist)
        area = cavity.getSurfaceArea()
        cavity.constructCOSMO()
        cavity.constructCOSMOgradient()
        print "    %4.1d      %4.1d       %8.4f    " % (order, nr_points, area * AtomicData.bohr_to_angs**2)

    print ""

###########################################################################
# Surface area of intersecting spheres
###########################################################################

def vdW_surface_area(atomlist, nr_slices=50):
    """
    surface area of intersecting van der Waals spheres

    Parameters
    ----------
    atomlist   :   list of tuples (Zi,[x,y,z]) with atoms

    Optional
    --------
    nr_slices  :   number of horizontal slices used to approximate the
                   surface area of a sphere

    Returns
    -------
    area       :   float, exposed area
    """
    # A sphere with vdW radius is placed around each atom
    nat = len(atomlist)
    # radii of spheres
    radii   = np.zeros(nat)
    # positions of spheres
    centers = np.zeros((3,nat))
    
    # bring data into the form expected by `surface_areas_spheres(...)`
    for i,(Zi,posi) in enumerate(atomlist):
        vdw_radius = AtomicData.vdw_radii[AtomicData.atom_names[Zi-1]]
        radii[i] = vdw_radius
        centers[:,i] = posi

    areas = surface_areas_spheres(centers, radii, nr_slices=nr_slices)
    
    return sum(areas)
        
def surface_areas_spheres(centers, radii, nr_slices=10, debug=0):
    """
    surface area of intersecting spheres.

    The general idea for computing the area consists in cutting the volume into slices (along the xy-plane).
    Then the circumference of each slice is obtained by considering the set of intersecting
    disks. The circumferences of all slices are then added weighted by their height. 

    Parameters
    ----------
    centers    :  numpy array of shape (n,3), centers of spheres in 3d space, 
                  centers[:,i] is the position (xi,yi,zi) of the i-th sphere
    radii      :  numpy array of length n, radii of spheres

    Optional
    --------
    nr_slices  :  number slices along z-axis, for each sphere
    debug      :  controls amount of printed messages (0 - no output, 1 - debug)
    
    Returns
    -------
    areas      :  numpy array of length n, exposed surface areas
                  areas[i] is the area of sphere i that lies on the surface
    """
    n = len(radii)
    # azimuthal angles
    # If the center of the sphere is located at z=0, the sphere is cut
    # at z_i = r*cos(thetas_slice[i])
 
    thetas_slice = np.linspace(0.0, np.pi, nr_slices)
    cos_slice = np.cos(thetas_slice)
    sin_slice = np.sin(thetas_slice)
    # differentials
    dthetas = np.ediff1d(thetas_slice, to_end=thetas_slice[-1]-thetas_slice[0])

    areas = np.zeros(n)
    for i in range(0, n):
        xi,yi,zi = centers[:,i]
        # surface area (circumference * height) of k-th slice
        for k in range(0, nr_slices):
            if debug:
                print " ==> cut plane k=%d <== " % k
            # radius of circle i (intersection of i-th sphere with k-th cut plane)
            radius_ik = radii[i]*sin_slice[k]
            # global z-coordinate of k-th cut plane
            z_cut_k = zi + radii[i]*cos_slice[k]

            # Find arcs of circle i occluded by circle j
            # An arc of a circle is defined by the angle of the first point
            # with the x-axis (phi_state) and the angle of the last point (phi_end).
            occluded_arcs = []
            for j in range(0, n):
                if i == j:
                    continue
                
                xj,yj,zj = centers[:,j]
                if abs(zj - z_cut_k) > radii[j]:
                    # sphere j does not intersect k-th cut plane
                    continue
                # radius of circle j (intersection of j-th sphere with k-th cut plane)
                radius_jk = np.sqrt(radii[j]**2 - (zj - z_cut_k)**2)

                # distance between circles i and j
                # vector between centers of circles in the xy-plane
                dx = xj-xi
                dy = yj-yi
                # distance
                d2 = dx**2 + dy**2
                d = np.sqrt(d2)

                if (d + radius_jk < radius_ik):
                    if debug:
                        print "case 1: circle %d lies inside circle %d" % (j,i)
                    # case 1:
                    #   circle j is fully contained in circle i => all area is accessible
                    phi_start = 0.0
                    phi_end   = 0.0
                elif (d + radius_ik < radius_jk):
                    if debug:
                        print "case 2: circle %d lies inside circle %d" % (i,j)
                    # case 2:
                    #   circle i is fully contained in circle j => no accessible surface area
                    phi_start = 0.0
                    # We have to subtract an infinitesimal angle, since otherwise 0 rad = 2*pi rad
                    phi_end   = 2.0*np.pi - 1.0e-12
                elif (radius_ik + radius_jk < d):
                    if debug:
                        print "case 3: circle %d and %d do not intersect" % (i,j)
                    # case 3:
                    #   circles do not intersect => all area of i is accessible
                    phi_start = 0.0
                    phi_end   = 0.0
                else:
                    if debug:
                        print "case 4: circles %d and %d intersect" % (i,j)
                    # case 4:
                    #   circles intersect 
                    # angle between vector from center i to center j and the x-axis
                    phi0 = np.arctan2(dy,dx)
                    # invoke law of cosines to get angles at which circle j intersects
                    # circle i
                    cosphi = (d2 + radius_ik**2 - radius_jk**2)/(2*d*radius_ik)
                    dphi = np.arccos(cosphi)

                    phi_start = phi0 - dphi
                    phi_end   = phi0 + dphi
                    
                if debug:
                    print "angles of crossing points of sphere %d and sphere %d (relative to center of sphere %d)" % (j,i,i)
                    print "  phi_start          = %s" % phi_start
                    print "  phi_end            = %s" % phi_end

                # Multiplies of 2*pi can be added to angles without changing them.
                # Make sure all angles are > 0
                phi_start += 2*np.pi
                phi_end   += 2*np.pi
                assert phi_start >= 0.0
                assert phi_end   >= 0.0
                occluded_arcs.append( (phi_start, phi_end) )

            # Intervals may wrap around at 2*pi, which makes the identification of overlapping
            # intervals harder. Therefore we split intervals containing the boundary 2*pi into
            # two.
            wrapped_occluded_arcs = wrap_intervals_2pi(occluded_arcs)
            # merge all occluded arcs
            inside_merged_arcs = merge_overlapping_intervals(wrapped_occluded_arcs)
            # find complementary arcs, which lie on the outside
            vmin = min([phi_start for (phi_start, phi_end) in inside_merged_arcs], 0.0)
            outside_arcs = complementary_intervals(inside_merged_arcs, vmin=vmin, vmax=vmin+2.0*np.pi)

            # compute length of exposed circumference of circle i
            # sum of exposed angles
            total_angle = 0.0
            for (phi_start, phi_end) in outside_arcs:
                assert phi_start <= phi_end
                total_angle += phi_end - phi_start
            if debug:
                print "Exposed angle of sphere i=%d in cut plane k=%d  :  %s rad" % (i,k,total_angle)

            # The area of the exposed surface of this slice is computed as
            # d(Area) =          d(Height)     x    d(Length)
            #         =    radius_i * dtheta   x  radius_ik * d(angle phi)
            total_area_ik = radii[i]*dthetas[k]  *  radius_ik*total_angle

            areas[i] += total_area_ik
            
    return areas
    
def unique_angle(phi):
    """
    Angles are only defined up to multiples of 2*pi. This function
    transforms angles into the unique range [0, 2*pi]
    """
    phi = phi % (2*np.pi)
    if (phi < 0.0):
        phi += 2*np.pi
    assert 0.0 <= phi <= 2*np.pi
    
    return phi

def wrap_intervals_2pi(intervals):
    """
    split an interval (phi_start, phi_end) into two non-overlapping intervals
    whenever it contains the angle 2*pi. The split is made at 2*pi

          (phi_start, phi_end)  ->  (phi_start, 2*pi), (0, phi_end)

    Parameters
    ----------
    intervals  :   list of tuples (phi_start, phi_end)

    Returns
    -------
    split_intervals  :  list of tuples (phi_start', phi_end')
                        with phi_start' < phi_end'
    """
    split_intervals = []
    for (phi_start, phi_end) in intervals:
        phi_start = unique_angle(phi_start)
        phi_end   = unique_angle(phi_end)
        if phi_start > phi_end:
            # interval wraps around at 2*pi
            # split it into two intervals at 2*pi ~ 0.0
            split_intervals += [(phi_start, 2*np.pi), (0.0, phi_end)]
        else:
            split_intervals.append( (phi_start, phi_end) )
            
    return split_intervals

def merge_overlapping_intervals(intervals):
    """
    given a set of intervals 
         
        I1,I2,...,In 

    merge the overlapping intervals, so as to obtain the minimal number of intervals 

        I1',I2',...,Im'  with m < n

    which cover the same range, i.e.

       I1' U I2' U ... U Im' = I1 U I2 U ... U In

    but are non-overlapping, i.e.

       Ii' intersection with Ij' = 0

    Parameters
    ----------
    intervals : list of tuples with start and endpoints of intervals
                [(start1,end1), (start2, end2), ..., (startN,endN)]
                 
    Parameters
    ----------
    merged    :  list of tuples with start end endpoints of non-overlapping,
                 intervals [(start1',end1'), ... (startM',endM')]
                 sorted such that end_I' < start_(I+1)'
    """
    if len(intervals) == 0:
        return intervals
    # sort intervals by starting point
    sorted_intervals = sorted(intervals, key=lambda (start,end): start)
    #
    merged_intervals = []
    # current interval
    startC,endC = sorted_intervals[0]
    for (startI,endI) in sorted_intervals[1:]:
        assert startC <= startI
        if endI <= endC:
            # case 1: interval I is fully contained in the current interval 
            #  |---------------|    current C
            #     |------|            I
            # skip it
            pass
        elif (startI <= endC):
            # case 2: interval I overlaps partially with the current interval
            #  |---------------|     current C
            #       |-------------|   I
            assert endC < endI
            # merge intervals
            #  |------------------|  merged interval
            endC = endI
        else:
            # interval I is separate from current interval
            #  |---------------|                      current
            #                     |-----------|         I
            assert endC < startI
            # add previous finished interval ...
            merged_intervals.append( (startC, endC) )
            # ... and create a new interval
            startC = startI
            endC = endI
    # append last interval
    merged_intervals.append( (startC, endC) )

    return merged_intervals

def complementary_intervals(intervals, vmin=0.0, vmax=2.0*np.pi):
    """
    given a set of non-overlapping, sorted intervals 

         I1, I2, ..., In

    contained in the range D = [vmin, vmax], find the complementary
    set of non-overlapping intervals

         D - (I1 U I2 U ... U In)

    Parameters
    ----------
    intervals : list of tuples with start and endpoints of sorted, non-overlapping intervals
                [(start1,end1), (start2, end2), ..., (startN,endN)]

    Returns
    -------
    complementary :  list of tuples with start and endpoints of complementary intervals
    """
    #
    #  vmin                                                                 vmax
    #   |        I1               I2                 I3                      |
    #   |    |---------|     |---------|   |------------------|              |
    #   |                                                                    |
    #
    # complementary intervals
    #
    #     I1'            I2'            I3'                       I4'               
    #   |----|         |-----|         |---|                  |--------------|      
    #
    n = len(intervals)
    if n == 0:
        return [(vmin, vmax)]
    
    complementary = []
    # Is there a gap between the first interval and the initial point vmin?
    # Then we need to add an interval at the beginning (I1' in the sketch above)
    start0, end0 = intervals[0]
    if start0 > vmin:
        startC = vmin
        endC = start0
        complementary.append( (startC, endC) )
    assert start0 >= vmin, "First interval [%s,%s] lies outside range [%s,%s] " % (start0, end0, vmin,vmax)
        
    # complementary interval in the inner region (I2' and I3')
    for i in range(0,n-1):
        startI, endI = intervals[i]
        startIp1, endIp1 = intervals[i+1]
        complementary.append( (endI, startIp1) )

        assert endIp1 <= vmax, "Interval [%s,%s] lies outside range [%s,%s] " % (startIp1, endIp1, vmin,vmax)
    
    # Is there a gap between the last interval and the enpoint vmax?
    # The we have to add an interval at the end (I4')
    startN, endN = intervals[-1]
    if endN < vmax:
        complementary.append( (endN, vmax) )

    return complementary

def exposed_angle(intervals):
    """
    compute the angle not covered by any of the n overlapping intervals I1,I2,...,In
    using the inclusion-exclusion principle (see M. Hall, "Combinatorial Theory", chapter 2):

       exposed angle = 2*pi - len(I1 U I2 U ... U In)

                     = 2*pi - sum len(I  ) +  sum  len(I  /\ I  )  -   sum    len(I  /\ I  /\  )
                               i1      i1    i1<i2      i1    i2     i1<i2<i3      i1    I2  I3

                                s                                             n
                      +...+ (-1)    sum    len(I   /\ ... /\ I  ) + ... + (-1)  len(I  /\ ... /\ I  )
                                 i1<...<is      i1            is                     i1           in

    """
    pass

########################################################
# testing
########################################################

def test_interval_merging():
    intervals = [(1,2), (-5,1.5), (2,5), (-2,3), (10,12), (11,13)]
    print intervals
    merged = merge_overlapping_intervals(intervals)
    print merged
    complementary = complementary_intervals(merged, vmin=-20, vmax=20)
    print complementary
    
def circumference_2d():
    """
    circumference of intersecting disks in 2d
    """
    import matplotlib.pyplot as plt
    from matplotlib import patches

    """
    centers = np.array([
        [0.0, 0.0],
        [0.4, 1.5],
        [1.0, 0.8],
        [2.0, 1.5]])
    radii = np.array([1.0, 0.9,0.8, 0.5])
    """
    """
    centers = np.array([
        [0.0, 0.0],
        [0.0, 0.0]])
    radii = np.array([1.0, 0.9])
    """

    centers = np.array([
        [0.0, 0.0],
        [0.0, 0.01]])
    radii = np.array([1.0, 1.0])


    n = len(radii)

    # plot circles and label them
    ax = plt.gca()
    for i in range(0, n):
        xi,yi = centers[i,:]
        ri = radii[i]
        circle = patches.Arc((xi,yi), ri/0.5, ri/0.5)
        ax.add_patch(circle)
        ax.text(xi,yi, "%d" % i, horizontalalignment='center', verticalalignment='center')

        
    # exposed circumference of sphere i
    exposed_lengths = []
    for i in range(0, n):
        xi,yi = centers[i,:]
        print "sphere %d at (%s,%s)" % (i,xi,yi)
        # find length of arcs of i occluded by sphere j
        # An arc of a circle is defined by the angle of the first point
        # with the x-axis (phi_state) and the angle of the last point (phi_end).
        occluded_arcs = []
        for j in range(0, n):
            if i == j:
                continue
            # distance between spheres i and j
            xj,yj = centers[j,:]
            print "sphere j at (%s,%s)" % (xj,yj)
            # vector between centers of circles in the xy-plane
            dx = xj-xi
            dy = yj-yi
            # distance
            d2 = dx**2 + dy**2
            d = np.sqrt(d2)

            if (d + radii[j] < radii[i]):
                print "case 1: sphere %d lies inside sphere %d" % (j,i)
                # case 1:
                #  sphere j is fully contained in sphere i => all area is accessible
                phi_start = 0.0
                phi_end   = 0.0
            elif (d + radii[i] < radii[j]):
                print "case 2: sphere %d lies inside sphere %d" % (i,j)
                # case 2:
                #  sphere i is fully contained in sphere j => no accessible surface area
                phi_start = 0.0
                # We have to subtract an infinitesimal angle, since otherwise 0 rad = 2*pi rad
                phi_end   = 2.0*np.pi - 1.0e-12
            elif (d > radii[i] + radii[j]):
                print "case 3: spheres %d and %d do not intersect" % (i,j)
                # case 3:
                # spheres do not intersect => all area of i is accessible
                phi_start = 0.0
                phi_end   = 0.0
            else:
                print "case 4: sphere %d and %d intersect" % (i,j)
                # case 4:
                # spheres intersect 
                # angle between vector from center i to center j and the x-axis
                phi0 = np.arctan2(dy,dx)
                # invoke law of cosines to get angles at which sphere j intersects
                # sphere i
                cosphi = (d2 + radii[i]**2 - radii[j]**2)/(2*d*radii[i])
                dphi = np.arccos(cosphi)

                phi_start = phi0 - dphi
                phi_end   = phi0 + dphi


            print "angles of crossing points of sphere %d and sphere %d (relative to center of sphere %d)" % (j,i,i)
            print "  phi_start          = %s" % phi_start
            print "  phi_end            = %s" % phi_end

            #assert phi_start <= phi_end
            # Multiplies of 2*pi can be added to angles without changing them.
            # Make sure all angles are > 0
            phi_start += 2*np.pi
            phi_end   += 2*np.pi
            assert phi_start >= 0.0
            assert phi_end   >= 0.0
            occluded_arcs.append( (phi_start, phi_end) )

            # plot section of circumference of sphere i occluded by sphere j
            occluded_arc = patches.Arc((xi,yi), radii[i]/0.5,radii[i]/0.5, theta1=phi_start*180.0/np.pi, theta2=phi_end*180.0/np.pi, color="red", lw=3)
            ax.add_patch(occluded_arc)

        # Intervals may wrap around at 2*pi, which makes the identification of overlapping
        # intervals harder. Therefore we split intervals containing the boundary 2*pi into
        # two.
        wrapped_occluded_arcs = wrap_intervals_2pi(occluded_arcs)
        # merge all occluded arcs
        inside_merged_arcs = merge_overlapping_intervals(wrapped_occluded_arcs)
        # find complementary arcs, which lie on the outside
        vmin = min([phi_start for (phi_start, phi_end) in inside_merged_arcs], 0.0)
        outside_arcs = complementary_intervals(inside_merged_arcs, vmin=vmin, vmax=vmin+2.0*np.pi)

        print " occluded_arcs         = %s" % occluded_arcs
        print " wrapped_occluded arcs = %s" % wrapped_occluded_arcs
        print " inside_merged_arcs    = %s" % inside_merged_arcs
        print " outside_arcs          = %s" % outside_arcs
        
        # plot exposed surface of sphere i
        for (phi_start, phi_end) in outside_arcs:
            exposed_arc = patches.Arc((xi,yi), radii[i]/0.5,radii[i]/0.5, theta1=phi_start*180.0/np.pi, theta2=phi_end*180.0/np.pi, color="green", lw=3)
            ax.add_patch(exposed_arc)

        # compute length of exposed circumference of sphere i
        # sum of exposed angles
        total_angle = 0.0
        for (phi_start, phi_end) in outside_arcs:
            assert phi_start <= phi_end
            total_angle += phi_end - phi_start
        # total length of exposed arcs
        #  L = radius * angle
        total_length_i = radii[i]*total_angle

        print "length of exposed circumference of sphere %d  :  %s" % (i, total_length_i)

        exposed_lengths.append( total_length_i )

    print "circumference of intersecting disks: %s" % np.sum(exposed_lengths)
    
    ax.set_aspect('equal')
    plt.xlim((-2.0, 4.0))
    plt.ylim((-2.0, 4.0))
    red_patch = patches.Patch(color='red', label='occluded surfaces')
    green_patch = patches.Patch(color='green', label='exposed surfaces')
    plt.legend(handles=[red_patch, green_patch])
    plt.legend()
    # hide axes
    plt.axis('off')

    plt.show()

    return exposed_lengths

def test_surface_areas_spheres():
    """
    check convergence of surface area as a function of the number of 
    horizontal slices
    """

    centers = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0001]
        ])
    radii = np.array([1.0, 1.0])

    """
    centers = np.array([
        [0.0, 0.0, 10.0]])
    radii = np.array([1.0])
    """
    
    print "exposed surface areas"
    print "  nr. slices            areas "
    for nr_slices in range(5, 50):
        areas = surface_areas_spheres(centers.transpose(), radii, nr_slices=nr_slices, debug=0)
        print "   %4.d        %s" % (nr_slices, areas)

    print "4 pi r**2 = %s" % (4.0 * np.pi * radii**2)


def test_vdw_surface_area():
    import sys
    import os.path
    if len(sys.argv) < 2:
        print "Usage: %s input.xyz" % os.path.basename(sys.argv[0])
        print "  compute area of van der Waals surface of a molecule"
        exit(-1)

    # load geometry
    xyz_file = sys.argv[1]
    atomlist = XYZ.read_xyz(xyz_file)[0]
    # compute vdW area
    area = vdW_surface_area(atomlist)

    print "van der Waals area: %s bohr^2  %s Ang^2" % (area, area*AtomicData.bohr_to_angs**2)
    
if __name__ == "__main__":
    import sys
    #atomlist = XYZ.read_xyz(sys.argv[1])[0]
    #check_analytic_gradients(atomlist)    
    #check_surface_area(atomlist)

    #circumference_2d()

    #test_interval_merging()    
    #test_surface_areas_spheres()
    test_vdw_surface_area()

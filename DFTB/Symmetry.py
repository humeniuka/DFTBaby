"""
determine the molecular symmetry for some groups automatically
"""
import numpy as np
import numpy.linalg as la

from DFTB.Modeling import MolecularCoords as MolCo
from DFTB import XYZ
from DFTB import AtomicData

# SYMMETRY ELEMENTS

class SymOp:
    def __call__(self, x):
        return self.apply(x)
    def transform(self, atomlist):
        """
        transform a molecular geometry by applying the symmetry operations
        to all its atoms
        """
        atomlist_trans = []
        for (Zi,posi) in atomlist:
            posi_trans = self.apply(posi)
            atomlist_trans.append( (Zi, posi_trans) )
        return atomlist_trans

class Id(SymOp):
    """
    identity 
    """
    def apply(self, x):
        return x
    def name(self):
        return "Identity"

class Inv(SymOp):
    """
    inversion about origin
    """
    def apply(self, x):
        return -x
    def name(self):
        return "Inversion"

class Ref(SymOp):
    """
    reflection through a plane 
    """
    def __init__(self, axis):
        """
        axis is the normal to the plane of reflection
        """
        self.axis = axis
        self.S = MolCo.reflection_matrix(axis)
    def apply(self, x):
        return np.dot(self.S, x)
    def name(self):
        return "Reflection perpendicular to axis %s" % (self.axis)

class Rot(SymOp):
    """
    rotation around an axis
    """
    def __init__(self, axis, angle):
        self.axis = axis
        self.angle = angle
        self.R = MolCo.rotation_matrix(axis, angle)
    def apply(self, x):
        return np.dot(self.R, x)
    def name(self):
        return "Rotation around axis %s by %2.4f deg" % (self.axis, self.angle*180.0/np.pi)
    def order(self):
        # order of the rotation axis
        # R^n = Id
        n = int(round(2.0*np.pi / self.angle))
        return n
    
class SRot(SymOp):
    """
    An improper rotation is rotation followed by a reflection perpendicular
    to the axis of rotation

      S_n = sigma C_n
    """
    def __init__(self, axis, angle):
        self.axis = axis
        self.angle = angle
        self.R = MolCo.rotation_matrix(axis, angle)
        self.S = MolCo.reflection_matrix(axis)
    def apply(self, x):
        return np.dot(self.S, np.dot(self.R, x))
    def name(self):
        return "Improper Rotation around axis %s by %2.4f deg" % (self.axis, self.angle*180.0/np.pi)

# SYMMETRY GROUPS AND CHARACTER TABLES
# the character tables are taken from Phil R. Bunker "Molecular Symmetry and Spectroscopy"

class SymGroup(object):
    # axes relative to which the symmetry operations are defined
    a = np.array([1,0,0])
    b = np.array([0,1,0])
    c = np.array([0,0,1])
    axes = [a,b,c]
    # rotational symmetry number
    sigma_rot = 1
    def __init__(self, axis_order=None):
        if not (axis_order is None):
            # order the axes so that the principal axes (the one with
            # the highest rotation order) becomes the x-axis
            ia,ib,ic = axis_order
            self.axes = self.axes[ia], self.axes[ib], self.axes[ic]
    def size(self):
        """
        number of symmetry elements in the group
        """
        return len(self.elems)
    def check_symmetry(self, atomlist):
        """
        test whether all symmetry operations preserve the molecular
        geometry (up to reordering of atoms)
        """
        issym = 1
        for e in self.elems:
            atomlist_e = e.transform(atomlist)
            if symmetry_related(atomlist, atomlist_e) == False:
                issym *= 0
                print "%s changes the geometry from" % e.name()
                print write_atomlist(atomlist)
                print "to"
                print write_atomlist(atomlist_e)
        if issym == 1:
            print ">>> Molecule HAS %s symmetry <<<" % self.name()
        else:
            print ">>> Molecule DOES NOT HAVE %s symmetry <<<" % self.name()

        return issym
    # 
    def symmetry_equivalents(self, x):
        """
        The symmetry equivalent points to a 3d vector x result by applying
        all symmetry operations of a group g to x

            x' = {e(x) for e elements of g}

        Parameters:
        ===========
        x: point in 3D space
        
        Returns:
        ========
        xts: list of tranformed points {e1(x),e2(x),...,eN(x)}
        """
        xts = []
        for i,e in enumerate(self.elems):
            xts.append( e(x) )
        return xts
    def assign_irrep(self, xts, fxs, tol=1.0e-2):
        """
        determine according to which irrep the function f transforms
        by observing the sign changes of the function.

        Parameters:
        ===========
        xts: list of symmetry equivalent points, for each element in the group
           there should be one transformed point
        fxs: list values of the function f on the points xts

        Returns:
        ========
        irrep: name of the irreducible representation or '?' is the assignment
           was not successful
        """
        f0 = fxs[0]   # identity element comes first
        charF = np.zeros(len(self.elems))
        for i,e in enumerate(self.elems):
            fi = fxs[i]
            sigi = fi/f0
            if abs(sigi - 1.0) < tol:
                sigi = 1
            elif abs(sigi + 1.0) < tol:
                sigi = -1
            else:
                #print "Cannot determine to which irrep the function belongs fi/f0 = %s" % sigi
                pass
            charF[i] = sigi
        # find the row in the character table that matches charF
#        print "Characters for f"
#        print charF
#        print "Character table"
#        print self.chars
        nirrep,nelem = self.chars.shape
        for i in range(0, nirrep):
            if np.sum(abs(charF - self.chars[i,:])) < 1.0e-10:
                irrep = self.irreps[i]
                break
        else:
            irrep = "?"
#        print " => Irrep = %s" % irrep
        return irrep

    def name(self):
        return self.__class__.__name__
    def rotational_symmetry_number(self):
        """
        the number of equivalent orientations of the molecule that can be achieved by
        rotations.
        """
        return self.sigma_rot
    
class C1(SymGroup):
    sigma_rot = 1
    def __init__(self, axes_order=None):
        super(C1,self).__init__(axes_order)
        #
        self.elems = [Id()]
        self.chars = np.array([
                [1]])
        self.irreps = ["A"]

class Cs(SymGroup):
    sigma_rot = 1
    def __init__(self, axes_order=None):
        super(Cs,self).__init__(axes_order)
        #
        a,b,c = self.axes
        # symmetry elements
        self.elems = [Id(), Ref(a)]
        # character table
        self.chars = np.array([
                [1, 1],
                [1,-1]])
        # names of irreps
        self.irreps = ["A'", "A\""]
        

class C2v(SymGroup):
    sigma_rot = 2
    def __init__(self, axes_order=None):
        super(C2v, self).__init__(axes_order)
        #
        a,b,c = self.axes
       # symmetry elements
#        self.elems = [Id(), Rot(b,np.pi), Ref(c), Ref(a)]
        self.elems = [Id(), Rot(c,np.pi), Ref(a), Ref(b)]
        # character table
        self.chars = np.array([
                [1, 1, 1, 1],
                [1, 1,-1,-1],
                [1,-1,-1, 1],
                [1,-1, 1,-1]])
        # names of irreps
        self.irreps = ["A1", "A2", "B1", "B2"]

class C2h(SymGroup):
    sigma_rot = 2
    def __init__(self, axes_order=None):
        super(C2h,self).__init__(axes_order)
        #
        a,b,c = self.axes
        # symmetry elements
        self.elems = [Id(), Rot(c,np.pi), Ref(c), Inv()]
        # character table
        self.chars = np.array([
                [1, 1, 1, 1],
                [1, 1,-1,-1],
                [1,-1,-1, 1],
                [1,-1, 1,-1]])
        # names of irreps
        self.irreps = ["AG", "AU", "BG", "BU"]

class D2h(SymGroup):
    sigma_rot = 4
    def __init__(self, axes_order=None):
        super(D2h,self).__init__(axes_order)
        #
        c,b,a = self.axes # this order is correct for ethene and tetrazine, but does it work in general?
        # symmetry elements
        self.elems = [Id(), Rot(a,np.pi), Rot(b,np.pi), Rot(c,np.pi), Ref(c), Ref(b), Ref(a), Inv()]
        # character table for non-degenerate irreps
        self.chars = np.array([
                [1, 1, 1, 1, 1, 1, 1, 1],
                [1, 1, 1, 1,-1,-1,-1,-1],
                [1, 1,-1,-1,-1,-1, 1, 1],
                [1, 1,-1,-1, 1, 1,-1,-1],
                [1,-1, 1,-1,-1, 1,-1, 1],
                [1,-1, 1,-1, 1,-1, 1,-1],
                [1,-1,-1, 1, 1,-1,-1, 1],
                [1,-1,-1, 1,-1, 1, 1,-1]])
        # names of irreps
        self.irreps = ["AG", "AU", "B1G", "B1U", "B2G", "B2U", "B3G", "B3U"]
        # degenerate irreps are not so simple

class D3h(SymGroup):
    sigma_rot = 6
    def __init__(self, axes_order=None):
        super(D3h, self).__init__(axes_order)
        #
        a,b,c = self.axes
        # symmetry elements
        self.elems = [Id(), Rot(c,np.pi*4.0/3.0), Rot(b,np.pi), Ref(c), SRot(c,np.pi*4.0/3.0), Ref(a)]
        # character table for non-degenerate irreps
        self.chars = np.array([
                [ 1, 1, 1, 1, 1, 1],
                [ 1, 1,-1, 1, 1,-1],
                [ 1, 1, 1,-1,-1,-1],
                [ 1, 1,-1,-1,-1, 1]])
        # names of irreps
        self.irreps = ["A1'", "A2'", "A1''", "A2''"]
        # degenerate irreps are not so simple
        # E' and E'' cannot be detected so easily

class D4h(SymGroup): # NOT TESTED
    sigma_rot = 4
    def __init__(self, axes_order=None):
        super(D4h, self).__init__(axes_order)
        #
        a,b,c = self.axes
        
        # symmetry elements
        self.elems = [Id(), Rot(c,np.pi/2.0), Rot(c,np.pi), Rot(a,np.pi), Rot(b,np.pi), Inv(), SRot(c,np.pi/2.0), Ref(c), Ref(a), Ref(b)]
        # character table for non-degenerate irreps
        self.chars = np.array([
                [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                [ 1, 1, 1,-1,-1, 1, 1, 1,-1,-1],
                [ 1,-1, 1, 1,-1, 1,-1, 1, 1,-1],
                [ 1,-1, 1,-1, 1, 1,-1, 1,-1, 1],
                [ 1, 1, 1, 1, 1,-1,-1,-1,-1,-1],
                [ 1, 1, 1,-1,-1,-1,-1,-1, 1, 1],
                [ 1,-1, 1, 1,-1,-1, 1,-1,-1, 1],
                [ 1,-1, 1,-1, 1,-1, 1,-1, 1,-1]])
        # names of irreps
        self.irreps = ["A1G", "A2G", "B1G", "B2G", "A1U", "A2U", "B1U", "B2U"]
        # degenerate irreps are not so simple
        # Eg and Eu have to be detected by looking at the orbitals or transition densities
    

class D6h(SymGroup):
    sigma_rot = 12
    def __init__(self, axes_order=None):
        super(D6h, self).__init__(axes_order)
        #
        a,b,c = self.axes
        # symmetry elements
        self.elems = [Id(), Rot(c,np.pi/3.0), Rot(c,np.pi*2.0/3.0), Rot(c, np.pi), \
                            Rot(a,np.pi), Rot(b, np.pi), Inv(), SRot(c,np.pi*2.0/3.0), SRot(c,np.pi/3.0), \
                            Ref(c), Ref(a), Ref(b)]
        # character table for non-degenerate irreps
        self.chars = np.array([
                [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                [ 1, 1, 1,-1,-1, 1, 1, 1, 1, 1,-1,-1],
                [ 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1],
                [ 1,-1, 1,-1,-1, 1, 1,-1, 1,-1,-1, 1],
                [ 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1],
                [ 1, 1, 1, 1,-1,-1,-1,-1,-1,-1, 1, 1],
                [ 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1],
                [ 1,-1, 1,-1,-1, 1,-1, 1,-1, 1, 1,-1]])
        # names of irreps
        self.irreps = ["A1g", "A2g", "B1g", "B2g", "A1u", "A2u", "B1u", "B2u"]
        # degenerate irreps are not so simple
        # E1g,E2g and E1u and E2u cannot be detected so easily

class D2d(SymGroup):
    sigma_rot = 2
    def __init__(self, axes_order=None):
        super(D2d, self).__init__(axes_order)
        # 
        a,b,c = self.axes
        d = (a+b)/np.sqrt(2.0) # C2 rotation axis perpendicular to S4 axis
        # symmetry elements
        self.elems = [Id(), SRot(c, np.pi/2.0), Rot(c, np.pi), Rot(d, np.pi),  Ref(a)]
        # character table for non-degenerate irreps
        self.chars = np.array([
                [1, 1, 1, 1, 1],
                [1, 1, 1,-1,-1],
                [1,-1, 1, 1,-1],
                [1,-1, 1,-1, 1]])
        # names of irreps
        self.irreps = ["A1", "A2", "B1", "B2"]
        # degenerate irreps are not so simple
    

# If you add new symmetries you also have to add them to the list 'point_groups'
# in the function 'detect_symmetry_brute_force' below.        
    
##############

def symmetry_related(atomlist1, atomlist2, tol=1.0e-3):
    """
    check whether atomlist1 and atomlist2 contain the same atoms, maybe reordered
    """
    Nat1 = len(atomlist1)
    Nat2 = len(atomlist2)
    assert Nat1 == Nat2
    # permutation matrix
    P = np.zeros((Nat1,Nat2))
    for i,(Zi,posi) in enumerate(atomlist1):
        for j,(Zj,posj) in enumerate(atomlist2):
            if Zi == Zj:
                dist = la.norm(np.array(posi)-np.array(posj))
                if dist < tol:
                    # atom i has been mapped to atom j
                    P[i,j] = 1
    # The determinant of a permutation matrix has to be +1 or -1
    if abs(la.det(P)) == 1.0:
        return True
    else:
        return False

def detect_symmetry(atomlist):
    return detect_symmetry_brute_force(atomlist)

def detect_symmetry_brute_force(atomlist):
    """
    Find symmetry by checking all implemented point groups
    """
    # List of the point groups that can be detected
    # ordered by increasing complexity
    o = axes_ordering(atomlist)
    point_groups = [C1(o), Cs(o), C2v(o), C2h(o), D2h(o), D2d(o), D3h(o), D4h(o), D6h(o)]
    hassym = np.zeros(len(point_groups))
    for i,g in enumerate(point_groups):
        hassym[i] = g.check_symmetry(atomlist)
    # use the most complex symmetry group
    isym = np.where(hassym == 1)[0][-1]
    g = point_groups[isym]
    print " => use the most complex of the detected symmetry groups, %s" % g.name()
    return g

def axes_ordering(atomlist_std):
    """
    reorder x,y and z axis so that the axis with the highest rotation
    order becomes the z-axis (c-axis)

    Parameters:
    ===========
    atomlist_std: list of atoms that are oriented in such a way
      that the axes of inertia concide with the x,y and z-axes    

    Returns:
    ========
    o: list of indeces for permutation, o[0] is the index of the principal axis
    """
    order = np.array([1,1,1])
    order_improp = np.array([0,0,0])
    axes = [np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1])]
    for iax,axis in enumerate(axes):
        # try 2-fold, 3-fold, 4-fold, ... rotations around axis
        for n in range(2,7):
            angle = 2.0*np.pi / float(n)
            # proper rotation
            R = Rot(axis, angle)
            atomlist_e = R.transform(atomlist_std)
            if symmetry_related(atomlist_std, atomlist_e) == True:
                order[iax] = n
            # improper rotation
            SR = SRot(axis, angle)
            atomlist_e = SR.transform(atomlist_std)
            if symmetry_related(atomlist_std, atomlist_e) == True:
                order_improp[iax] = n        
    print "Rotation orders of x,y and z-axes: %s" % order
    if len(np.unique(order)) == 1:
        # all axes have same proper rotation order
        # use order of improper rotations
        order = order_improp
    # 
    if len(np.unique(order)) == 1:
        # all axes have the same proper and improper rotation order
        # count the number of atoms that lie on each axis
        atoms_on_axis = np.array([0,0,0])
        for iax, axis in enumerate(axes):
            for (Zi, posi) in atomlist_std:
                c = np.dot(posi, axis)/( la.norm(posi)*la.norm(axis) )
                if abs(abs(c)-1.0) < 1.0e-8:
                    # position vector is parallel to axis
                    atoms_on_axis[iax] += 1
        print "Number of atoms on x,y and z-axes: %s" % atoms_on_axis
        # increase order of axis with many atoms
        order += atoms_on_axis
    # index of principle axis comes last
    o = np.argsort(order)
    print "Axis ordering: %s" % o
    return o


#### IDENTIFY POINT GROUP ####################################
def is_linear(atomlist):
    pass

def has_inversion(atomlist):
    pass

def has_two_or_more_Cn(atomlist):
    pass

def has_Cn(atomlist):
    pass

def has_reflection(atomlist):
    pass

########################################

def write_atomlist(atomlist):
    Nat = len(atomlist)
    pos = XYZ.atomlist2vector(atomlist)
    txt  = "                in Angstrom:\n"
    txt +=     "Atom          X           Y           Z\n"
    for i in range(0, Nat):
        txt += ("  %s  " % i).rjust(3)
        txt += "%+4.7f  %+4.7f  %+4.7f\n" % tuple(pos[3*i:3*(i+1)]*AtomicData.bohr_to_angs)
    txt += "\n"
    return txt

########################################


if __name__ == "__main__":
    import sys
    atomlist = XYZ.read_xyz(sys.argv[1])[0]
    atomlist_std = MolCo.standard_orientation(atomlist)
    group = detect_symmetry_brute_force(atomlist_std)
#    group = C2v()
    group.check_symmetry(atomlist_std)
#    atomlist_std_trans = group.elems[1].transform(atomlist_std)
#    print atomlist_std
#    print atomlist_std_trans
#    print symmetry_related(atomlist_std, atomlist_std_trans)

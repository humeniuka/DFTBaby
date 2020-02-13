#!/usr/bin/env python
"""
systematically create all kind of 2D Zn-porphyrin oligomers
"""
from DFTB import XYZ
from DFTB import Periodic
from DFTB import AtomicData
from DFTB.DFTB2 import DFTB2
from DFTB import utils
from DFTB.Molden import MoldenExporter

import numpy as np
import numpy.linalg as la
import scipy.linalg as sla
import os
from os.path import expandvars, expanduser

class Flake:
    def __init__(self, grid, orientation=None, weight=1, parent=None):
        """
        Parameters:
        ===========
        grid: 2D array filled with 0s or 1s
             at the position of the 1s monomer units are placed
        orientation: 2D list filled with indeces which give
             the orientation of a flake. If the position [i,j] is not occupied
             orientation[i][j] should be None
        weight: How many symmetry equivalent flakes do exist?
        parent: from which flake can this flake be generated?
        """
        self.grid = grid
        if type(orientation) == type(None):
            self.orientation = np.zeros(self.grid.shape, dtype=int)
        else:
            self.orientation = orientation
        self.weight = weight
        self.parent = parent
        # moves are: left, right, up, down
        self.moves = [(-1,0), (1,0), (0,1), (0,-1)]
        # rotation out of plane
        twist_angle = 5.0 * np.pi/180.0
        rot_angle = 90.0 * np.pi/180.0
        x_axis = np.array([1.0, 0.0, 0.0])
        y_axis = np.array([0.0, 1.0, 0.0])
        diag1 = np.array([1.0,1.0,0.0])/np.sqrt(2.0)
        diag2 = np.array([-1.0,1.0,0.0])/np.sqrt(2.0)
        self.rotations = [(x_axis, 0.0), 
                          (x_axis, -twist_angle), (y_axis, -twist_angle), (x_axis, twist_angle), (y_axis, twist_angle),
                          (diag1, rot_angle), (diag2, rot_angle)]
        # which moves are compatible with which rotations
        self.compatible_rotations =\
            {"meso_beta_wire":
                                    {((-1,0), 0): 0,
                                     (( 1,0), 0): 0,
                                     (( 0,1), 0): None,
                                     ((0,-1), 0): None},                 
             "meso_meso_wire":
                                    {((-1,0), 0): 5,
                                     (( 1,0), 0): 5,
                                     (( 0,1), 0): None,
                                     ((0,-1), 0): None,

                                     ((-1,0), 5): 0,
                                     (( 1,0), 5): 0,
                                     (( 0,1), 6): None,
                                     ((0,-1), 6): None,
                                     
                                     ((-1,0), 6): None,
                                     (( 1,0), 6): None,
                                     (( 0,1), 5): None,
                                     ((0,-1), 5): None},
             "meso_beta": 
                                    {((-1,0), 0): 0,
                                     (( 1,0), 0): 0,
                                     (( 0,1), 0): 0,
                                     ((0,-1), 0): 0},
             "meso_meso":
                                    {((-1,0), 0): 5,
                                     (( 1,0), 0): 5,
                                     (( 0,1), 0): 6,
                                     ((0,-1), 0): 6,

                                     ((-1,0), 5): 0,
                                     (( 1,0), 5): 0,
                                     (( 0,1), 6): 0,
                                     ((0,-1), 6): 0,
                                     
                                     ((-1,0), 6): None,
                                     (( 1,0), 6): None,
                                     (( 0,1), 5): None,
                                     ((0,-1), 5): None,

                                     ((-1,0), 1): 4,
                                     (( 1,0), 1): 4,
                                     (( 0,1), 1): 2,
                                     ((0,-1), 1): 2,
                                     
                                     ((-1,0), 2): 3,
                                     (( 1,0), 2): 3,
                                     (( 0,1), 2): 1,
                                     ((0,-1), 2): 1,
                                     
                                     ((-1,0), 3): 2,
                                     (( 1,0), 3): 2,
                                     (( 0,1), 3): 4,
                                     ((0,-1), 3): 4,

                                     ((-1,0), 4): 1,
                                     (( 1,0), 4): 1,
                                     (( 0,1), 4): 3,
                                     ((0,-1), 4): 3}
             }

    def grow_meso_beta(self):
        """
        grow flake by attaching one new porphyrin monomer at all possible sites.
        The new monomer is connected in the meso- and in both beta-positions.
        """
        n,m = self.grid.shape
        grown_flakes = []
        for i in range(0,n):
            for j in range(0,m):
                if self.grid[i,j] == 1:
                    # fuse an additional monomer to the existing flake
                    for lr,tb in self.moves:
                        i_new = i+lr
                        j_new = j+tb
                        if 0 <= i_new < n and 0 <= j_new < m \
                                and self.grid[i_new,j_new] == 0:
                            grid_new = np.copy(self.grid)
                            grid_new[i_new,j_new] = 1
                            grown_flakes.append( Flake(grid_new, weight=self.weight, parent=self) )
        return grown_flakes
    def grow_connectivity(self, connectivity):
        """
        grow flake by adding monomers as specified by connectivity.

        connectivity can be:
          "meso_beta" - 2D planar grid of meso and beta connected porphyrins
          "meso_meso" - monomers are only connected at the meso position
                        and are rotated by 90 deg
          "meso_beta_wire" - 1D chain of meso and beta connected porphyrins
        """
        assert connectivity in self.compatible_rotations.keys()
        n,m = self.grid.shape
        grown_flakes = []
        for i in range(0,n):
            for j in range(0,m):
                if self.grid[i,j] == 1:
                    # fuse an additional monomer to the existing flake
                    for lr,tb in self.moves:
                        i_new = i+lr
                        j_new = j+tb
                        if 0 <= i_new < n and 0 <= j_new < m \
                                and self.grid[i_new,j_new] == 0:
                            grid_new = np.copy(self.grid)
                            grid_new[i_new,j_new] = 1
                            orientation_new = np.copy(self.orientation)
                            orient = self.compatible_rotations[connectivity][((lr,tb),self.orientation[i,j])]
                            if orient == None:
                                continue
                            orientation_new[i_new,j_new] = orient
                            grown_flakes.append( Flake(grid_new, orientation=orientation_new, weight=self.weight, parent=self) )
        return grown_flakes        
    def __str__(self):
        xmin, ymin, xmax, ymax = grid_bounding_box(self.grid)
        return grid_to_txt(self.grid[xmin:xmax+1,:][:,ymin:ymax+1])
    def __cmp__(self, other):
        """
        check if two flakes are related by symmetry operations
        """
        if other == None:
            return +1
        # create all symmetry equivalent flakes
        x = shift_to_origin(self.grid)
#        print "Compare two grids"
#        print "================="
#        print "Grid 1"
#        print grid_to_txt(x)
        for T in symmetry_operations:
            symequiv = T(other.grid)
            Ty = shift_to_origin(symequiv)
#            print "Grid 2"
#            print grid_to_txt(Ty)
            if (Ty == x).all():
#                print "ARE THE SAME"
                return 0
        if np.sum(self.grid) < np.sum(other.grid):
            return -1
        else:
            return +1
    def build_xyz(self, a1, a2, a3, unitcell, meso, beta):
        n,m = self.grid.shape
        center = (n/2, m/2)
        flake_atomlist = []
        unitcell_vec = XYZ.atomlist2vector(unitcell)
        for i in range(0, n):
            for j in range(0, m):
                if self.grid[i,j] == 1:
                    if self.orientation[i,j] != 0:
                        axis, angle = self.rotations[self.orientation[i,j]]
                        Rot = axis_angle2rotation(axis, angle)
                    else:
                        Rot = np.identity(3)
                    R = (i-n/2)*a1 + (j-m/2)*a2 + 0.0*a3
                    shifted_unitcell_vec = np.zeros(unitcell_vec.shape)
                    for at in range(0, len(unitcell)):
                        # rotate and translate
                        shifted_unitcell_vec[3*at:3*(at+1)] = \
                            np.dot(Rot, unitcell_vec[3*at:3*(at+1)]) + R
                    flake_atomlist += XYZ.vector2atomlist(shifted_unitcell_vec, unitcell)
                    # hydrogenize or add meso/beta substituents
                    for (lr,tb) in self.moves:
                        i_hyd = i+lr
                        j_hyd = j+tb
                        direction = lr*a1 + tb*a2
                        direction /= la.norm(direction)
                        if 0 <= i_hyd < n and 0 <= j_hyd < m:
                            if self.grid[i_hyd,j_hyd] != 1:
                                # no neighbouring porphine in that direction => add both meso and beta hydrogens
                                hydrogens = meso + beta # hydrogens can also contain protecting groups, so it's not only H-atoms!
                            elif self.grid[i_hyd,j_hyd] == 1 \
                                    and (((self.orientation[i,j] == 0) and (self.orientation[i_hyd,j_hyd] in [5,6])) \
                                             or ((self.orientation[i_hyd,j_hyd] == 0) and (self.orientation[i,j] in [5,6]))):
                                hydrogens = beta
                            else:
                                hydrogens = []
                            # add hydrogens
                            hydadd = []
                            for at in range(0, len(hydrogens)):
                                (Zat, pos) = hydrogens[at]
                                hyddir = np.array(pos)
                                hyddir /= la.norm(hyddir)
                                angle = np.arccos(np.dot(hyddir, direction))
                                if abs(angle) < np.pi/4.0:
                                    hydadd += [(Zat, np.dot(Rot, np.array(pos)) + R)]
                            flake_atomlist += hydadd
                # meso-meso
                
                    
        return flake_atomlist
    def hueckel(self, orbe0, Sh, Sv, H0h, H0v):
        """
        perform a Hueckel calculation using the matrix elements passed as arguments

        Parameters:
        ===========
        orbe0: energies of basis orbitals
        Sh: overlap between orbitals on neighbouring horizontally fused porphyrins
        Sv: overlap between orbitals on neighbouring vertically fused porphyrins
        H0h: matrix elements between orbitals on neighbouring horizontally fused porphyrins
        H0v: matrix elements between orbitals on neighbouring vertically fused porphyrins

        Results:
        ========
        en_tot, HLgap: total pi-energy and HOMO-LUMO gap in Hartree
        orbe, orbs: orbe[i] is the Hueckel energy belonging to the orbital with the coefficients
           orbs[:,i]
        """
        # arrange the cells of the polyomino (or flake, we'll use the terms interchangeably)
        # on a linear grid. This defines a mapping (i,j) -> k
        n,m = self.grid.shape
        k = 0
        grid2chain = {}  # mapping from 2D grid to linear chain
        chain2grid = {}  # inverse mapping
        for i in range(0, n):
            for j in range(0, m):
                if self.grid[i,j] == 1:
                    # occupied
                    grid2chain[(i,j)] = k
                    chain2grid[k] = (i,j)
                    k += 1
        # number of monomers in the flake == number of 1s in the grid
        N = np.sum(self.grid[self.grid == 1])
        # each monomer contributes 4 valence electrons in the 2 orbitals a1u and a2u
        nelec = 4*N
        assert k == np.sum(self.grid[self.grid == 1])
        # construct matrix elements for the whole flake from the matrix elements between two porphyrins
        norb = Sh.shape[0]  # number of orbitals per porphyrin
        Zero = np.zeros((norb,norb))
        # list of block matrices
        S  = [[Zero for l in range(0, N)] for k in range(0, N)]
        H0 = [[Zero for l in range(0, N)] for k in range(0, N)]
        for k in range(0, N):
            (i1,j1) = chain2grid[k]
            for l in range(0, N):
                (i2,j2) = chain2grid[l]
                if k == l:
                    # diagonal elements are orbital energies
                    S[k][l] = np.eye(norb)
                    H0[k][l] = np.diag(orbe0)
                else:
                    # matrix elements between orbitals on porphyrin (i1,j1) and (i2,j2)
                    # are non-zero only if they are nearest neighbours
                    if i2-i1 == 1 and j2-j1 == 0:
                        # horizontal nearest neighbours: #-#
                        S[k][l]  = Sh
                        H0[k][l] = H0h
                    elif i2-i1 == 0 and j2-j1 == 1:
                        # vertical nearest neighbours: #
                        #                              |
                        #                              #
                        S[k][l]  = Sv
                        H0[k][l] = H0v
        #
        S = np.bmat(S)
        H0 = np.bmat(H0)
        # make matrices hermitian by adding symmetric matrix elements, that were skipped above
        S  = 0.5*(S+S.transpose())
        H0 = 0.5*(H0+H0.transpose())
        # Hueckel energies and coefficients
        orbe, orbs = sla.eigh(H0, S)
        # fill the orbitals of lowest energy with 2 electrons each according to the Aufbau principle
        HOMO = nelec/2-1
        print "number of porphyrin monomers N = %d" % N
        print "number of electron pairs = %d" % (nelec/2)
        print "index of HOMO = %d" % HOMO
        LUMO = HOMO+1
        en_tot = 2.0 * np.sum(orbe[:HOMO+1]) # sum over energies of doubly occupied orbitals
        HLgap = orbe[LUMO]-orbe[HOMO]
        print "--------------------------------------------------"
        print self
        print "Hueckel total energy: %8.6f eV" % (en_tot*AtomicData.hartree_to_eV)
        print "Hueckel total energy/site: %8.6f" % ((en_tot/float(N))*AtomicData.hartree_to_eV)
        print "Hueckel HOMO: %8.6f eV" % (orbe[HOMO]*AtomicData.hartree_to_eV)
        print "Hueckel LUMO: %8.6f eV" % (orbe[LUMO]*AtomicData.hartree_to_eV)
        print "Hueckel HOMO-LUMO gap: %8.6f eV" % ( HLgap*AtomicData.hartree_to_eV )
        print "--------------------------------------------------"
        # save additional data as member variables for later use
        # mapping from k to (i,j)
        self.chain2grid = chain2grid
        # indexes of HOMO and LUMO
        self.HOMO = HOMO
        self.LUMO = LUMO
        # number of orbitals per site
        self.norb = norb
        # number of sites
        self.N = N
        
        return en_tot, HLgap, orbe, orbs
    def count_new_bonds(self):
        """
        count the number of new conjugated C~C bonds that were formed by fusing
        """
        new_bonds = 0
        n,m = self.grid.shape
        for i1 in range(0, n):
            for j1 in range(0, m):
                if self.grid[i1,j1] == 1:
                    for i2 in range(0, n):
                        for j2 in range(0, m):
                            if self.grid[i2,j2] == 1:
                                if abs(i1-i2)+abs(j1-j2) == 1:
                                    new_bonds += 3  # 3 C-C bonds per fusion
        new_bonds /= 2   # a-b and b-a are calculated twice
        print "new conjugated C-C bonds: %d (%d x 3)" % (new_bonds, new_bonds/3)
        return new_bonds
        
def grid_to_txt(grid):
    n,m = grid.shape
    txt = ""
    for i in range(0,n):
        for j in range(0,m):
            if grid[i,j] == 1:
                txt+= "#"
            elif grid[i,j] == -1:
                # internal holes
                txt+= "H"
            elif 1 < grid[i,j] < 10:
                # show color
                txt+= "%d" % grid[i,j]
            else:
                txt+= "*"
        txt += "\n"
    return txt

def grid_bounding_box(grid):
    n,m = grid.shape
    xmin,xmax = n,0
    ymin,ymax = m,0
    for i in range(0, n):
        for j in range(0, m):
            if grid[i,j] == 1:
                xmin = min(xmin, i)
                ymin = min(ymin, j)
                xmax = max(xmax, i)
                ymax = max(ymax, j)
    return (xmin,ymin, xmax,ymax)

def grid_center(grid):
    """center of mass of polyomino"""
    n,m = grid.shape
    # center of grid is obtained as the average position
    x0, y0 = 0, 0
    # size of polyomino, number of 1's in grid
    size = 0
    for i in range(0,n):
        for j in range(0,m):
            if grid[i,j] == 1:
                x0 += i
                y0 += j
                size += 1
    x0 /= float(size)
    y0 /= float(size)
    return (x0,y0)
                
def plot_polyomino(flake, ax, origin=(0,0), s=1.0):
    """
    draw a porphyrin flake as a polyomino on a matplotlib axis
    s is the size of a single box
    """
    # find bounding box of polyomino
    (xmin,ymin, xmax,ymax) = grid_bounding_box(flake.grid)
    # find center of mass of polyomino
    x0, y0 = grid_center(flake.grid)
    # width and height
    w = xmax-xmin
    h = ymax-ymin
    for i in range(0,w+1):
        for j in range(0, h+1):
            if flake.grid[i+xmin,j+ymin] == 1:
                #
                d = 0.5
                # The (0,0) origin of the grid is in the upper left corner
                # and the y-axis points from top to bottom, while the y-axis
                # used for plotting points upwards.
                xpos = + (np.array([i-d,i+d,i+d,i-d,i-d]) - (x0-xmin)) - origin[0]
                ypos = - (np.array([j-d,j-d,j+d,j+d,j-d]) - (y0-ymin)) - origin[1]
                ax.plot(s*xpos,s*ypos, ls="-", lw=3, color="black")
                ax.fill(s*xpos,s*ypos, color="grey")

###############################
# only meso-meso connected

def axis_angle2rotation(axis, angle):
    """
    Parameters:
    ===========
    axis: 3d unit vector
    angle: angle in radians

    Returns:
    ========
    R: rotation matrix, which performs active rotation by angle around axis
    """
    # cross product matrix such that (axis x v) = K.v
    K = np.array([[0.0,  -axis[2],  axis[1]],
                  [axis[2], 0.0 ,  -axis[0]],
                  [-axis[1], axis[0], 0.0]])
    R = np.identity(3) + np.sin(angle)*K + (1.0-np.cos(angle))*np.dot(K,K)
    return R

###############################

def shift_to_origin(grid):
    row, col = np.nonzero(grid)
    shifted_grid = np.zeros(grid.shape)
    shifted_grid[(row-row.min(), col-col.min())] = grid[(row,col)]
    return shifted_grid

# Symmetry operations
def Id(grid):
    return np.copy( grid )

# from these 3 operations all other symmetry operations can be generated
def mirror_vertical(grid):
    return np.copy( grid[:,::-1] )

def mirror_horizontal(grid):
    return np.copy( grid[::-1,:] )

def mirror_diag_1(grid):
    return grid.transpose()

#
def mirror_diag_2(grid):
    return rotate_180_deg(mirror_diag_1(grid))

def rotate_90_deg(grid):
    return mirror_vertical(mirror_diag_1(grid))

def rotate_180_deg(grid):
    return mirror_vertical(mirror_horizontal(grid))

def rotate_270_deg(grid):
    return rotate_180_deg(rotate_90_deg(grid))


symmetry_operations = [Id, mirror_vertical, mirror_horizontal, \
                           mirror_diag_1, mirror_diag_2, \
                           rotate_90_deg, rotate_180_deg, rotate_270_deg]

#
def group_by_equivalence(A):
    """
    group elements in A together if they are equivalent

    Parameters:
    ===========
    A: list of elements with __cmp__ method

    Returns:
    ========
    [E1,E2,...,EN], where Ei contains elements that are equivalent
    """
    equiv_groups = []
    removed = [0 for a in A] # elements that have been dealt with get a 1
    for i,a in enumerate(A):
        if removed[i] == 1:
            continue
        equiv_groups.append([a])
        removed[i] = 1
        for j,b in enumerate(A):
            if removed[j] == 1:
                continue
            if a == b:
                equiv_groups[-1].append(b)
                removed[j] = 1

    assert sum([len(eqg) for eqg in equiv_groups]) == len(A)
    return equiv_groups

def grow_flakes(max_generations, connectivity="meso_beta", thresh=0.01):
    if connectivity in ["meso_meso_wire", "meso_beta_wire"]:
        size = (max_generations+1)*2
    else:
        size = max((max_generations+1)*2,11)
    gr = np.zeros((size,size), dtype=int)
    orient = np.zeros((size,size), dtype=int)
    # place a monomer in the center
    gr[size/2,size/2] = 1
    orient[size/2,size/2] = 0

    fl = Flake(gr, orient)
    grown_flakes = {0: [fl]}
    for gen in range(1, max_generations+1):
        flakes = []
        for fl in grown_flakes[gen-1]:
#            flakes += fl.grow_meso_beta()
            flakes += fl.grow_connectivity(connectivity)
        equiv_groups = group_by_equivalence(flakes)
        # find unique flakes
        grown_flakes[gen] = []
        for eqg in equiv_groups:
            # select one representative of equivalence group
            flake_repr = eqg[0]
            weight = sum([fl.weight for fl in eqg])
            flake_repr.weight = weight
            grown_flakes[gen].append(flake_repr)
        # normalize weights in each generation
        weight_sum = 0.0
        for fl in grown_flakes[gen]:
            weight_sum += fl.weight
        for fl in grown_flakes[gen]:
            fl.weight = float(fl.weight) / float(weight_sum)

        if thresh > 0.0:
            # remove isomers with low probability
            grown_flakes[gen] = [fl for fl in grown_flakes[gen] if fl.weight > thresh]
            # normalize again weights in each generation
            weight_sum = 0
            for fl in grown_flakes[gen]:
                weight_sum += fl.weight
            for fl in grown_flakes[gen]:
                fl.weight = float(fl.weight) / float(weight_sum)

        weight_sum = 0
        for fl in grown_flakes[gen]:
            weight_sum += fl.weight
        print "WEIGHT SUM = %s" % weight_sum
        assert abs(weight_sum - 1.0) < 1.0e-8

    # count flakes
    nr_flakes = sum([len(v) for v in grown_flakes.values()])
    print "%s flakes were grown" % nr_flakes
    """
    # Show flakes
    for gen in grown_flakes.keys():
        print "Generation %s (%s flakes)" % (gen, len(grown_flakes[gen]))
        print "========================="
        for fl in grown_flakes[gen]:
            print "weight = %s" % fl.weight
            print fl
    """
    return grown_flakes


def grow_square_flakes(max_generations):
    """
    Create the monomer, 2x2, 3x3, ..., NxN square oligomers
    """
    grown_flakes = {}
    size = max(max_generations*2+1,11)
    for gen in range(1, max_generations+1):
        grid_square = np.zeros((size,size), dtype=int)
        grid_square[size/2-gen/2:size/2+(gen+1)/2,size/2-gen/2:size/2+(gen+1)/2] = 1
        square_flake = Flake(grid_square, weight=1.0, parent=None)        
        grown_flakes[gen] = [square_flake]
    return grown_flakes

##### Pattern Functions ############
def sierpinsky_pattern(F):
    if F == None:
        return np.array([[1]], dtype=int)
    n,m = F.shape
    F_new= np.zeros((3*n,3*n))
    for i in range(0, n):
        for j in range(0, n):
            # repeat F in all blocks except for the middle
            for k in range(0, 3):
                for l in range(0, 3):
                    if k == 1 and l == 1:
                        continue
                    F_new[k*n+i,l*n+j] = F[i,j]
    return F_new

def edge_pattern(F):
    #   1111
    #   1010
    if F == None:
        return np.array([[1,1],[1,0]], dtype=int)
    n,m = F.shape
    F_new= np.zeros((2*n,2*n))
    for i in range(0, n):
        for j in range(0, n):
            # repeat F in all blocks except for the middle
            for k in range(0, 2):
                for l in range(0, 2):
                    F_new[k*n+i,l*n+j] = F[i,j]
    return F_new

##############################

def grow_self_similar_flakes(max_generations, pattern_func):
    """
    creates fractal flake
    """
    grown_flakes = {}
    def embed_flake(F):
        # add margin with 0s
        n,m = F.shape
        grid = np.zeros((n+2,n+2))
        for i in range(0, n):
            for j in range(0, n):
                grid[i+1,j+1] = F[i,j]
        return grid

    F = pattern_func(None)
    grid = embed_flake(F)
    grown_flakes[0] = [Flake(grid, weight=1.0, parent=None)]
    for gen in range(1, max_generations+1):
        F = pattern_func(F)
        # embed 
        flake = Flake(embed_flake(F), weight=1.0, parent=None)
        grown_flakes[gen] =[flake]
    return grown_flakes


def grow_square_flakes(max_generations):
    """
    Create the monomer, 2x2, 3x3, ..., NxN square oligomers
    """
    grown_flakes = {}
    size = max(max_generations*2+1,11)
    for gen in range(1, max_generations+1):
        grid_square = np.zeros((size,size), dtype=int)
        grid_square[size/2-gen/2:size/2+(gen+1)/2,size/2-gen/2:size/2+(gen+1)/2] = 1
        square_flake = Flake(grid_square, weight=1.0, parent=None)        
        grown_flakes[gen] = [square_flake]
    return grown_flakes


############################ RATE EQUATIONS ###################################
def rate_equations_matrix(flakes, c):
    """
    change of concentrations in one time step

    all reaction are of the type:
       F^i_n + Monomer -> F_(n+1)
    """
    nFl = len(flakes)
    
    dc = np.zeros(nFl)
    M = np.zeros((nFl,nFl))
    for i,flakeI in enumerate(flakes):
        # How does the population of flake I change depending on the
        # population of flakes J?
        for j,flakeJ in enumerate(flakes):
            # I -> J
            # flake I is the parent of flake J
            # F_I + Monomer -> F_J
            if flakeJ.parent == flakeI:
                    dc[i] -= flakeI.weight * c[i] * c[0]
                    M[i,i] -= flakeI.weight
            # flake J is the parent of flake I
            # J -> I
            # F_J + Monomer -> F_I
            if flakeI.parent == flakeJ:
                dc[i] += flakeJ.weight * c[j] * c[0]
                M[i,j] += flakeJ.weight
    return M

from scipy.integrate import ode

def integrate_rate_equations(M, c0, t0, t1, dt):
    def f(t,c):
        dcdt = np.dot(M, c) * c[0]
        return dcdt
    integrator = ode(f)
    integrator.set_integrator('vode', method='bdf')
    integrator.set_initial_value(c0,t0)
    ts = [t0]
    cs = [c0]
    while integrator.successful() and integrator.t < t1:
        integrator.integrate(integrator.t + dt)
        cs.append(integrator.y)
        ts.append(integrator.t)
        
    return np.array(ts),np.vstack(cs)

def plot_concentrations(flakes, ts, cs):
    from matplotlib import pyplot as plt
    tot = np.sum(cs, axis=1)
    for n in range(0, len(flakes)):
        plt.plot(ts, cs[:,n], label="%d" % n)
    plt.plot(ts,tot, label="Total")
    plt.xlabel("reaction time")
    plt.ylabel("rel. concentrations")
    plt.legend()
    plt.show()

def solve_rate_equations(grown_flakes, out_concentrations, out_weights):
    """
    solve rate equations and save the time dependent concentrations
    to file
    """
    generations = grown_flakes.keys()
    generations.sort()

    flakes = reduce(lambda l1,l2: l1+l2, [grown_flakes[gen] for gen in generations])

    c0 = np.zeros(len(flakes))
    c0[0] = 1.0
    M = rate_equations_matrix(flakes, c0)
#    print "M"
#    print M
    ts,cs = integrate_rate_equations(M, c0, 0.0, 10.0, 0.5)
#    plot_concentrations(flakes, ts, cs)

    # save time-dependent concentrations
    nt, nfl = cs.shape
    time_conc = np.zeros((nt, nfl+1))
    time_conc[:,1:] = cs
    time_conc[:,0] = ts
    if out_concentrations != "":
        np.savetxt(expandvars(expanduser(out_concentrations)), time_conc)
    # How many monomers are contained in each flake?
    # the absorption needs to be divided by this number to make
    # spectra with different number of monomers comparable
    nr_monomers = []
    for gen in generations:
        nr_monomers += [float(gen+1) for f in grown_flakes[gen]]
    nr_monomers = np.array(nr_monomers, dtype=float)
    # Dividing the concentrations by the number of monomers gives
    # the weights each absorption spectrum should have in the average
    print "TOTAL NUMBER OF FLAKES = %s" % len(nr_monomers)
    time_weights = np.zeros((nt,nfl+1))
    time_weights[:,1:] = cs/nr_monomers
    time_weights[:,0] = ts
    if out_weights != "":
        np.savetxt(expandvars(expanduser(out_weights)), time_weights)

############################ GEOMETRY DATA ####################################
# B3lyp optimized Zn-Porphyrin
def zn_porphene_lattice(cc1=2.56, substitutions=["meso-phenyl"]):
    d = 13.133 # diameter of porphyrine in bohr
    a = d + cc1 

    # B3-Lyp/TZVP optimized geometry of Zn-Porphyrin (with d4h symmetry)
    unitcell = [
     (30, (0.0, 0.0, 0.0)),
     (6, (-8.024201735853538, -1.2849475318240677, 0.0)),
     (6, (-5.415109975750299, -2.0829259044732664, 0.0)),
     (7, (-3.8805825531626366, 0.0, 0.0)),
     (6, (-5.415109975750299, 2.0829259044732664, 0.0)),
     (6, (-8.024201735853538, 1.2849475318240677, 0.0)),
     (6, (-4.5813685387672445, -4.5813685387672445, 0.0)),
     (6, (-2.0829259044732664, -5.415109975750299, 0.0)),
     (7, (0.0, -3.8805825531626366, 0.0)),
     (6, (2.0829259044732664, -5.415109975750299, 0.0)),
     (6, (1.2849475318240677, -8.024201735853538, 0.0)),
     (6, (-1.2849475318240677, -8.024201735853538, 0.0)),
     (6, (4.5813685387672445, -4.5813685387672445, 0.0)),
     (6, (5.415109975750299, -2.0829259044732664, 0.0)),
     (7, (3.8805825531626366, 0.0, 0.0)),
     (6, (5.415109975750299, 2.0829259044732664, 0.0)),
     (6, (8.024201735853538, 1.2849475318240677, 0.0)),
     (6, (8.024201735853538, -1.2849475318240677, 0.0)),
     (6, (4.5813685387672445, 4.5813685387672445, 0.0)),
     (6, (2.0829259044732664, 5.415109975750299, 0.0)),
     (6, (1.2849475318240677, 8.024201735853538, 0.0)),
     (6, (-1.2849475318240677, 8.024201735853538, 0.0)),
     (6, (-2.0829259044732664, 5.415109975750299, 0.0)),
     (7, (0.0, 3.8805825531626366, 0.0)),
     (6, (-4.5813685387672445, 4.5813685387672445, 0.0))]
    hydrogens_meso = [
     (1, (6.02929737820229, 6.02929737820229, 0.0)),
     (1, (6.02929737820229, -6.02929737820229, 0.0)),
     (1, (-6.02929737820229, -6.02929737820229, 0.0)),
     (1, (-6.02929737820229, 6.02929737820229, 0.0))]
    CH_meso_bondlength = 2.04768 # in a.u.
    hydrogens_beta = [
     (1, (9.632056195976029, 2.5403907717884526, 0.0)),
     (1, (9.632056195976029, -2.5403907717884526, 0.0)),
     (1, (2.5403907717884526, -9.632056195976029, 0.0)),
     (1, (-2.5403907717884526, -9.632056195976029, 0.0)),
     (1, (-9.632056195976029, -2.5403907717884526, 0.0)),
     (1, (-9.632056195976029, 2.5403907717884526, 0.0)),
     (1, (-2.5403907717884526, 9.632056195976029, 0.0)),
     (1, (2.5403907717884526, 9.632056195976029, 0.0))]

    # substitute meso-hydrogens by protecting Phenyl groups
    # Phenyl lies in the x-y plane, the carbon atom without a hydrogen lies roughly at (0,0,0)
    phenyl = [
     (1, (0.33741057526076673, -0.09488314188654774, 4.08698598453918)),
     (1, (5.0465321497599005, -0.08717305985314573, 4.067521806856818)),
     (1, (7.3840287113326, 0.007710082033402008, -0.01970984206087817)),
     (1, (5.012120239507878, 0.09488314188654774, -4.08698598453918)),
     (1, (0.30297976774885876, 0.08717305985314573, -4.067521806856818)),
     (6, (1.3626058213247185, -0.05327137561803985, 2.294373014513328)),
     (6, (4.006294684826861, -0.048943903104194116, 2.2833937065196848)),
     (6, (5.318444822256522, 0.004327472513845734, -0.010979307993643545)),
     (6, (3.9869060961840415, 0.05327137561803985, -2.294373014513328)),
     (6, (1.3432172326818987, 0.048943903104194116, -2.2833937065196848)),
     (6, (0.031067095252237498, -0.004327472513845734, 0.010979307993643545))]

    protection_meso = phenyl

    substituents_meso = []
    xaxis = np.array([1.0, 0.0, 0.0])
    CC_bondlength = 2.7779  # bond length between protecting group and porphyrin in a.u.
                            # both C's are sp^2 hybridized => bond length of 147 pm
    zaxis = np.array([0.0, 0.0, 1.0])
    angles = np.pi * np.array([1.0/4.0, -1.0/4.0, 5.0/4.0, 3.0/4.0])
    for ih,(Zh, pos_meso) in enumerate(hydrogens_meso):
        for (Z, pos_prot) in protection_meso:
            RotZ = axis_angle2rotation(zaxis, angles[ih])
            RotX = axis_angle2rotation(xaxis, 0.0) #(90-61.5) * np.pi/180.0) # Phenyl rings not orthogonal
            Rot = np.dot(RotZ, RotX)
            substituents_meso.append( (Z, np.array(pos_meso) + np.dot(Rot, np.array(pos_prot + (CC_bondlength - CH_meso_bondlength)*xaxis))) )

    protection_beta = phenyl
    substituents_beta = []
    angles = np.pi * np.array([1.0/4.0, -1.0/4.0, -1.0/4.0, 5.0/4.0, 5.0/4.0, 3.0/4.0, 3.0/4.0, 1.0/4.0])
    for ih,(Zh, pos_beta) in enumerate(hydrogens_beta):
        for (Z, pos_prot) in protection_beta:
            RotZ = axis_angle2rotation(zaxis, angles[ih])
            RotX = axis_angle2rotation(xaxis, 0.0) #(90-61.5) * np.pi/180.0) # Phenyl rings not orthogonal
            Rot = np.dot(RotZ, RotX)
            substituents_beta.append( (Z, np.array(pos_beta) + np.dot(Rot, np.array(pos_prot + (CC_bondlength - CH_meso_bondlength)*xaxis))) )

    if "meso-phenyl" in substitutions:
        meso = substituents_meso
    else:
        meso = hydrogens_meso
    if "beta-phenyl" in substitutions:
        beta = substituents_beta
    else:
        beta = hydrogens_beta

    cm = center_of_mass(unitcell + meso + beta)
    unitcell = shift_atomlist(unitcell, cm)
    meso = shift_atomlist(meso, cm)
    beta = shift_atomlist(beta, cm)

    lat = Periodic.PorpheneTetragonal(a)    
    return lat, unitcell, meso, beta

def center_of_mass(atomlist):
    """
    shift atomlist to center of mass
    """
    CM = np.zeros(3)
    M = 0.0
    for Zi,posi in atomlist:
        mi = AtomicData.atom_masses[AtomicData.atom_names[Zi-1]]
        CM += mi*np.array(posi)
        M += mi
    CM /= M
    return CM

def shift_atomlist(atomlist, vec):
    shifted = []
    for Zi,posi in atomlist:
        shifted.append( (Zi, (posi[0]-vec[0], posi[1]-vec[1], posi[2]-vec[2])) )
    return shifted

def flakes_to_xyz(grown_flakes, out_xyz, cc1=2.56, substitutions=[]):
    generations = grown_flakes.keys()
    generations.sort()

    lat, unitcell, meso, beta = zn_porphene_lattice(cc1=cc1, substitutions=substitutions)
    a1,a2,a3 = lat.getLatticeVectors("real")

    if out_xyz[-4:] == ".xyz":
        try:
            os.remove(opts.out_xyz)
        except OSError as e:
            print e
    for gen in generations:
        print "Generation %s (%s flakes)" % (gen, len(grown_flakes[gen]))
        print "========================="
        for fl in grown_flakes[gen]:
            print "weight = %s" % fl.weight
            print fl
            flake_atomlist = fl.build_xyz(a1,a2,a3,unitcell, meso, beta)
            # reorder atomlist so that all hydrogens are placed at the end of the file
            flake_atomlist = XYZ.hydrogens_to_end(flake_atomlist)
            XYZ.write_xyz(expandvars(expanduser(out_xyz)), [flake_atomlist], title="weight=%d" % fl.weight, mode="a")

##################### Hueckel calculations ###################
"""
Hueckel calculation for porphyrin polyominos using only the four Gouterman frontier orbitals 
as a minimal basis for each lattice point.
"""

def Gouterman_matelems():
    """
    computes overlap and hamiltonian matrix elements between the 4 frontier orbitals (a1u,a2u,eg,eg)
    of neighbouring porphyrine squares in a porphyrine flake (polyomino).
    The matrix elements between vertically and horizontally fused porphyrins will be related
    by a rotation of 90 degrees and will be different.
    
    Returns:
    ========:
    orbe: energies of the four orbitals
    Sh,Sv: two 4x4 matrices with overlaps between frontier orbitals of 
        horizontally and vertically fused porphyrins.
    H0h,H0v: two 4x4 matrices with matrix elements of the 0th-order DFTB hamiltonian 
        (neglecting interaction between partial charges) 
    """
    # perform DFTB calculation on the monomer to obtain the MO coefficients for the frontier
    # orbitals.
    lat, unitcell, meso, beta = zn_porphene_lattice(substitutions=[])
    monomer = unitcell+meso+beta

    dftb = DFTB2(monomer, missing_reppots="dummy")
    dftb.setGeometry(monomer)
    dftb.getEnergy()    
    orbs = dftb.getKSCoefficients()
    orbe = dftb.getKSEnergies()
    HOMO,LUMO = dftb.getFrontierOrbitals()

    # save DFTB MOs so that we can identify the symmetries of the orbitals.
    molden = MoldenExporter(dftb, title="Gouterman orbitals")
    molden.export(molden_file="/tmp/Gouterman_orbitals.molden")

    # In DFTB the orbital symmetries are:
    #   H-2  a1u
    #   H-1  rubbish
    #   H    a2u
    #   L    eg
    #   L+1  eg
    # The H-1 belongs to the sigma-framework and has the wrong energetic position in DFTB.
    Gouterman_orbitals = np.array([HOMO-2,HOMO,LUMO,LUMO+1], dtype=int) # indeces of G. orbitals
    # remove orbitals belonging to hydrogen
    valorbs = dftb.getValorbs()
    mu = 0  # iterates over all orbitals
    hydrogen_orbs = []  # indeces of orbitals belonging to hydrogen that should be removed
    for i,(Zi,posi) in enumerate(monomer):
        for (ni,li,mi) in valorbs[Zi]:
            if (Zi == 1):
                hydrogen_orbs.append(mu)
            mu += 1
    # remove rows or columns belonging to hydrogen
    orbs = np.delete(orbs, hydrogen_orbs, 0)
    orbs = np.delete(orbs, hydrogen_orbs, 1)
    nmo,nmo = orbs.shape
    
    a1,a2,a3 = lat.getLatticeVectors()
    C = orbs[:,Gouterman_orbitals]
    # energies of Gouterman orbitals
    orbe = orbe[Gouterman_orbitals]
    # horizontally fused
    # matrix elements between AOs
    H0h,Sh = dftb._constructH0andS_AB(nmo, nmo, unitcell, shift_atomlist(unitcell, a1))
    # transform to the basis of frontier orbitals
    Sh = np.dot(C.transpose(), np.dot(Sh,C))
    H0h = np.dot(C.transpose(), np.dot(H0h,C))
    # vertically fused
    # matrix elements between AOs
    H0v,Sv = dftb._constructH0andS_AB(nmo, nmo, unitcell, shift_atomlist(unitcell, a2))
    # transform to the basis of frontier orbitals
    Sv = np.dot(C.transpose(), np.dot(Sv,C))
    H0v = np.dot(C.transpose(), np.dot(H0v,C))

    #
    orb_names = ["a1u", "a2u", "eg", "eg"]

    print "Lattice vectors"
    print "   a1 = %s bohr" % a1
    print "   a2 = %s bohr" % a2
    print "orbital energies:"
    for name,en in zip(orb_names, orbe*AtomicData.hartree_to_eV):
        print "  energy(%3s) = %8.4f eV" % (name, en)
    print "matrix elements with neighbouring porphyrin"
    print "along lattice vector a1 = %s" % a1
    print "overlap Sh"
    print utils.annotated_matrix(Sh, orb_names, orb_names)
    print "hamiltonian H0h"
    print utils.annotated_matrix(H0h, orb_names, orb_names)
    print "along lattice vector a2 = %s" % a2
    print "overlap Sv"
    print utils.annotated_matrix(Sv, orb_names, orb_names)
    print "hamiltonian H0v"
    print utils.annotated_matrix(H0v, orb_names, orb_names)

    return orbe, Sh,Sv, H0h,H0v


if __name__ == "__main__":
    import sys
    import optparse
    usage = "Usage: python %s" % sys.argv[0]

    parser = optparse.OptionParser(usage)
    parser.add_option("--max_generations", dest="max_generations", help="generate all flakes with at most (max_generations+1) monomer units [default: %default]", default=2, type=int)
    parser.add_option("--out_xyz", dest="out_xyz", help="save the flake geometries to this xyz-file [default: %default]", default="/tmp/flakes.xyz") 
    parser.add_option("--out_concentrations", dest="out_concentrations", help="save time dependent concentations to this file [default: %default]", default="/tmp/concentrations.dat")
    parser.add_option("--out_weights", dest="out_weights", help="save time dependent concentrations divided by the number of monomers in each isomer to this file [default: %default]. This has nothing to do with the statistical weights of each isomer due to symmetry.", default="/tmp/weights.dat")
    parser.add_option("--connectivity", dest="connectivity", help="Specify in which ways monomers can be linked together in 2D: 'meso_beta', 'meso_meso', 'meso_beta_wire' or 'meso_meso_wire', 'square' [default: %default]", default="meso_beta")
    parser.add_option("--weight_threshold", dest="weight_threshold", help="Flake isomers with a relative weight lower than this threshold are not grown further [default: %default]", type=float, default=0.0)
    parser.add_option("--fused_cc_length", dest="fused_cc_length", help="Length of a fused meso-meso or beta-beta bond in Angstrom [default: %default]", type=float, default=2.56*AtomicData.bohr_to_angs)
    parser.add_option("--substitutions", dest="substitutions", help="Replace hydrogens by protecting groups. List of substitutions patterns, possible values 'meso-phenyl' => replace meso hydrogens by phenyl groups, 'beta-phenyl' => replace beta hydrogens by phenyl groups [default: %default]", default="[]")
    parser.add_option("--hueckel", dest="hueckel", help="Perform a Hueckel calculation with the 4 Gouterman orbitals as a basis for porphyrin in the flake. Energies of HOMO, LUMO and and total energies are written to the console for each flake that is generated. [default: %default]", type=int, default=0)
    parser.add_option("--polyomino_figures", dest="polyomino_figures", help="If set to 1, the porphyrin flakes are sketched as little connected boxes (polyominos). The resulting png-files are saved to polyomino_####.svg [default: %default]", default=0, type=int)
    
    (opts,args) = parser.parse_args()

    # Grow flakes by adding one monomer at a time
    if opts.connectivity == "square":
        grown_flakes = grow_square_flakes(opts.max_generations)
    elif opts.connectivity == "sierpinsky":
        grown_flakes = grow_self_similar_flakes(opts.max_generations, sierpinsky_pattern)
    elif opts.connectivity == "edge":
        grown_flakes = grow_self_similar_flakes(opts.max_generations, edge_pattern)
    else:
        grown_flakes = grow_flakes(opts.max_generations, opts.connectivity, opts.weight_threshold)

    # Build geometries
    if opts.out_xyz != "":
        flakes_to_xyz(grown_flakes, opts.out_xyz, \
                          opts.fused_cc_length / AtomicData.bohr_to_angs, \
                          substitutions=eval(opts.substitutions))
    #
    if opts.polyomino_figures == 1:
        # extent of largest polyomino
        width = 1
        height = 1
        for gen,flake_generation in grown_flakes.iteritems():
            for flake in flake_generation:
                (xmin, ymin, xmax, ymax) = grid_bounding_box(flake.grid)
                width = max(xmax-xmin, width)
                height = max(ymax-ymin, height)
        # plot squares
        import matplotlib.pyplot as plt
        i = 0
        for gen,flake_generation in grown_flakes.iteritems():
            for flake in flake_generation:
                #
                ax = plt.gca()
                ax.clear()
                ax.set_xlim((-width,width))
                ax.set_ylim((-height,height))
                ax.axis('off')
                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)
                ax.set_aspect('equal', 'datalim')
                plot_polyomino(flake, ax)
                plt.savefig("polyomino_%.4d.svg" % (i+1), transparent=True)
                #
                i += 1

    # find time dependent concentrations
    solve_rate_equations(grown_flakes, opts.out_concentrations, opts.out_weights)

    if opts.hueckel == 1:
        print "HUECKEL CALCULATION with Gouterman orbitals"
        # matrix elements between Gouterman orbitals fused horizontally or vertically
        orbe, Sh,Sv, H0h,H0v = Gouterman_matelems()
        i = 0
        for gen,flake_generation in grown_flakes.iteritems():
            print "Generation %d" % gen
            hueckel_out = "polyomino_hueckel_gen_%d.dat" % (gen+1)
            fh = open(hueckel_out, "w")
            print>>fh, "# ID   N    EN(TOT)    HOMO-LUMO gap    (-1)*(NR. NEW BONDS)"
            for flake in flake_generation:
                print "Flake %d" % i
                en_tot, HLgap, orbeHueckel, orbsHueckel = flake.hueckel(orbe, Sh,Sv, H0h,H0v)
                new_bonds = flake.count_new_bonds()
                print>>fh,"%d  %d  %10.6f %10.6f  %d" % (i,gen+1,en_tot, HLgap, -new_bonds)
                i += 1
            fh.close()
            print "Hueckel results for generation %d written to file %s" % (gen, hueckel_out)

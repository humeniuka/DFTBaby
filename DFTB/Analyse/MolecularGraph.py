"""
This module accomplishes two tasks:
  - it identifies disconnected molecular fragments
  - and brings atoms into canonical order, so that molecules can be compared
    
"""

from DFTB import XYZ, AtomicData

import numpy as np

##### DETECTION OF FRAGMENTS ###########################

class AtomNode(dict):
    def __init__(self, i, Z):
        self.i = i  # atom label
        self.Z = Z  # atomic number
        # the child nodes are stored in the dictionary
        # The carbene H-C-H could be translated into the following tree:
        #  v = AtomNode(0,6)
        #  v[1] = AtomNode(1,1)
        #  v[2] = AtomNode(2,1)
        
def atomlist2graph(atomlist, hydrogen_bonds=False, conmat=None):
    """
    Based on the adjacency matrix a molecular graph is constructed. 
    """
    Nat = len(atomlist)
    if conmat is None:
        # adjacency matrix, conmat[i,j] = 1 if there is a bond between atoms i and j
        conmat = XYZ.connectivity_matrix(atomlist, thresh=1.3, hydrogen_bonds=hydrogen_bonds)
    else:
        # check shape of adjacency matrix 
        assert conmat.shape == (Nat,Nat)
    # visited flags the atoms that have been added to the graph to avoid
    # repetitions and to identify disconnected parts
    visited = [0 for i in range(0, Nat)]
    # atomic numbers
    Zs = [Zi for (Zi, posi) in atomlist]
    # This recursive function starts at one atoms and adds all the connected
    # atoms as child nodes
    def make_graph(i):
        graph = AtomNode(i,Zs[i])
        visited[i] = 1
        connected = np.where(conmat[i,:] == 1)[0]
        for j in connected:
            if visited[j] == 1:
                continue
            graph[j] = make_graph(j)
            visited[j] = 1
        return graph
    #
    fragment_graphs = []
    while True:
        # find the first atom that has not been visited yet
        try:
            start = visited.index(0)
        except ValueError:
            # all atoms have been visited
            break
        graph = make_graph(start)
        fragment_graphs.append(graph)
    #print "number of disconnected fragments: %d" % len(fragment_graphs)
    return fragment_graphs

def graph2atomlist(graph, atomlist_ref):
    """
    converts a molecular graph back to the original geometry
    """
    atomlist = [atomlist_ref[graph.i]]
    for node in graph.values():
        atomlist += graph2atomlist(node, atomlist_ref)
    return atomlist

def graph2indeces(graph):
    """
    converts a molecular graph to a list of atom indeces into the list of atoms
    """
    indeces = [graph.i]
    for node in graph.values():
        indeces += graph2indeces(node)
    return indeces

def disconnected_fragments(atomlist, hydrogen_bonds=False):
    """
    separate disconnected molecular fragments into separate atomlists

    Parameters:
    ===========
    atomlist: list of tuples (Z,[x,y,z])
    
    Optional:
    =========
    hydrogen_bonds: flag that controls whether fragments connected by hydrogen bonds
       should be considered connected or not.
       

    Returns:
    ========
    fragments: list of fragments
    """
    fragment_graphs = atomlist2graph(atomlist, hydrogen_bonds=hydrogen_bonds)
    fragments = []
    for g in fragment_graphs:
        fragments.append( graph2atomlist(g, atomlist) )
    return fragments
    
########## CANONICAL ORDERING ###############################

def morgan_ordering(atomlist, hydrogen_bonds=False):
    """
    order atoms canonically using a variant of Morgan's algorithm according
    to 
     Morgan, H.L., The Generation of a Unique Machine Description for Chemical Structures
    and 
     "Morgan Revisited" by Figueras,J.  J.Chem.Inf.Comput.Sci. 1993, 33, 717-718

    Optional:
    =========
    hydrogen_bonds: hydrogen bonds are included in the connectivity matrix

    Returns:
    ========
    atomlist_can: canonically reordered atoms
    A_can: adjacency matrix for canonical order of atoms
    """
    Nat = len(atomlist)
    # adjacency matrix, A[i,j] = 1 if there is a bond between atoms i and j
    A = XYZ.connectivity_matrix(atomlist, thresh=1.3, hydrogen_bonds=hydrogen_bonds)
    # initial invariants are the atom types
    v = [Zi for (Zi, posi) in atomlist]
    k = len(np.unique(v))
    #
    for i in range(0, 100):
    #while True:
        vnext = np.dot(A, v)
        # count the number of unique invariants
        knext = len(np.unique(vnext))
        v = vnext
        if knext <= k:
            # 
            break
        k = knext
    else:
        raise Exception("Morgan algorithm did not converge in %d iterations!" % i)    
    # Symmetry equivalent atoms have same labels in v, while symmetry-inequivalent
    # ones should have different labels
    #print "Labels"
    #print v
    
    # Now the molecular graph is traversed starting with the atom with the lowest
    # label. Nearest neighbours are added in the order of their invariants.

    # visited flags the atoms that have been added to the graph to avoid
    # repetitions and to identify disconnected parts
    visited = [0 for i in range(0, Nat)]

    def traverse_ordered(i, ordering=[]):
        # atoms are added as they are encountered while traversing the graph
        ordering += [i]
        # mark atom i as visited
        visited[i] = 1
        # sort connected atoms by values of invariants
        connected = np.where(A[i,:] == 1)[0]
        sort_indx = np.argsort(v[connected])
        # The child nodes are added in depth-first order. Nodes on the same level
        # are sorted by v.
        for j in connected[sort_indx]:
            if visited[j] == 1:
                continue
            ordering = traverse_ordered(j, ordering)
        return ordering
    # start with the atom with the lowest label
    istart = np.argmin(v)
    ordering = traverse_ordered(istart)

    assert len(ordering) == Nat, "It seems that the molecule contains disconnected fragments. Try to compute the canonical ordering for each fragment separately."

    #print "Canonical Ordering"
    #print ordering

    # Atoms are reorderd according to the canonical ordering
    atomlist_can = []
    for i in ordering:
        atomlist_can.append( atomlist[i] )
    # The rows of the adjacency matrix are also reordered
    A_can = A[ordering,:][:,ordering]
    # check
    """
    print "A reordered"
    print A_can
    print "A calculated"
    print XYZ.connectivity_matrix(atomlist_can)
    """
    return atomlist_can, A_can

def identifier(atomlist_can, A_can):
    """
    This function finds a string for a molecule in canonical order. The identifier
    is composed of the atom characters in canonical order and is not necessarily unique.
    E.g. for methane => 'CHHHH'
    """
    idstr = ""
    for (Zi, posi) in atomlist_can:
        idstr += AtomicData.atom_names[Zi-1]
    return idstr

def fragment_trajectory(geometries):
    """
    Parameters:
    ===========
    geometries: list of atomlists along a trajectory

    Each geometry is decomposed into separate fragments and each fragment is assigned a string
    representation.

    Returns:
    ========
    fragtraj: list of fragment identifiers for each geometry
    fragdic: dictionary that maps fragment identifiers to geometries
    """
    fragtraj = []
    fragdic = {} # dictionary that translates fragment identifiers into molecular geometries
    for atomlist in geometries:
        fragments = disconnected_fragments(atomlist)
        fragment_labels = [identifier( *morgan_ordering(atomlist) ) for atomlist in fragments]
        #
        for ifrag,fl in enumerate(fragment_labels):
            if not fragdic.has_key(fl):
                fragdic[fl] = fragments[ifrag]
        #
        fragtraj.append( fragment_labels )
    return fragtraj, fragdic
        
if __name__ == "__main__":
    import sys
    atomlist = XYZ.read_xyz( sys.argv[1] )[0]
    fragments = disconnected_fragments(atomlist)

    atomlists_ordered = [morgan_ordering(atomlist)[0] for atomlist in fragments]
        
    XYZ.write_xyz("/tmp/fragments.xyz", atomlists_ordered) #atomlists_frag)
    
    print fragment_trajectory( XYZ.read_xyz( sys.argv[1] ) )

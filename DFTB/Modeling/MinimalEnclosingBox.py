"""
This module finds the box with minimal volume that encloses a molecule.

References:

   [1] Freeman,H.;\ \ Shapira,~R.
       Determining the Minimum-Area Encasing Rectangle for an Arbitrary Closed Curve.
       Communications of the ACM, 1975, 18, 409-413

   [2] O'Rourke,~J.
       Finding Minimal Enclosing Boxes.
       International Journal of Computer and Information Sciences, 1985, 14, 183-199
"""
from DFTB import XYZ, AtomicData
import DFTB.Modeling.MolecularCoords as MolCo

import numpy as np
import numpy.linalg as la
from scipy.spatial import ConvexHull

class MinimalEnclosingRectangle(object):
    def __init__(self, hull2D):
        self.hull2D = hull2D
        # list of possible bounding rectangles
        # and their areas
        rectangles = []
        areas = []
        # We iterate over edges of the convex hull as we know from a the proof in [1]
        # that the minimum-area rectangle is flush with one of the edges of the convex
        # polygon.
        for edge in self.hull2D.simplices:
            # coordinates of vertices in the edge
            pos0 = self.hull2D.points[edge[0]]
            pos1 = self.hull2D.points[edge[1]]
            # find two orthogonal vectors u,v which are parallel and perpendicular to the edge
            # u is parallel to the edge
            u = pos1-pos0
            u /= la.norm(u)
            (a,b) = u
            # v is perpendicular perpendicular
            v = np.array([b,-a])
            # consistency checks
            assert abs(np.dot(u,u) - 1.0) < 1.0e-10
            assert abs(np.dot(v,v) - 1.0) < 1.0e-10
            assert abs(np.dot(u,v)) < 1.0e-10
            # Now we project all vertex positions onto the axes
            # and find the minimal and maximum distances
            xs = []  # vertex positions the coordinate system whose origin lies at pos0
            ys = []  # and whose axes are aligned with u and v
            for pos in self.hull2D.points:
                dpos = pos-pos0
                x = np.dot(dpos, u)
                y = np.dot(dpos, v)
                xs.append(x)
                ys.append(y)
            # extent of bounding rectangle aligned with this edge
            xmin = min(xs)
            xmax = max(xs)
            ymin = min(ys)
            ymax = max(ys)
            # area
            area = (xmax-xmin)*(ymax-ymin)
            #
            rectangles.append( (pos0, u,v, [xmin,xmax,ymin,ymax]) )
            areas.append(area)
        # select rectangle with minimum area
        i = np.argmin(areas)
        (pos0, u,v, [xmin,xmax,ymin,ymax]) = rectangles[i]
        #
        self.area = areas[i]
        # compute the vertices of the minimal enclosing rectangle
        self.vertices = np.array(
            [pos0 + xmin*u + ymin*v,
             pos0 + xmax*u + ymin*v,
             pos0 + xmax*u + ymax*v,
             pos0 + xmin*u + ymax*v])
        self.edges = [[0,1],
                      [1,2],
                      [2,3],
                      [3,0]]
        
class MinimalEnclosingBox(object):
    def __init__(self, hull3D):
        self.hull3D = hull3D
        # list of possible oriented bounding boxes
        # and their volumes
        boxes = []
        volumes = []
        # find unique edges
        edges = []
        for simplex in self.hull3D.simplices:
            a,b,c = simplex
            for e in ( [a,b],[b,c],[c,a] ):
                
                if not (e in edges or e[::-1] in edges):
                    edges.append( e )
        # iterate over the edges of the convex polyhedron
        for edge in edges:
            # coordinates of vertices in the edges
            pos0 = self.hull3D.points[edge[0]]
            pos1 = self.hull3D.points[edge[1]]
            # w is the unit vector passing through the edge
            w = pos1-pos0
            w /= la.norm(w)
            # find 2 additional unit vectors u,v so that u,v,w forms a right-handed coordinate system
            utilde = np.random.rand(3)
            # ensure that utilde is not parallel to w
            while abs(np.dot(utilde, w)) < 1.0e-10:
                utilde = np.random.rand(3)
            # subtract projection onto w
            u = utilde - np.dot(w, utilde) * w
            u /= la.norm(u)
            v = -np.cross(u,w)
            assert abs(np.dot(u,w)) < 1.0e-10
            assert abs(np.dot(u,v)) < 1.0e-10
            assert abs(np.dot(v,w)) < 1.0e-10
            # All points are projected onto the plane perpendicular to w.
            xs = []
            ys = []
            zs = []
            for pos in self.hull3D.points:
                dpos = pos-pos0
                x = np.dot(u,dpos)
                y = np.dot(v,dpos)
                z = np.dot(w,dpos)
                xs.append(x)
                ys.append(y)
                zs.append(z)
            points2d = np.vstack([xs,ys]).transpose()
            # compute the convex hull of the projected points ...
            hull2D = ConvexHull(points2d, qhull_options="QbB Qt")
            # ... and bound the projected points by a rectangle in the plane perpendicular to w
            rect2d = MinimalEnclosingRectangle(hull2D)
            
            # extent of bounding box along w-axis
            zmin = min(zs)
            zmax = max(zs)
            # volume of the 3d box
            volume = (zmax-zmin)*rect2d.area
            #
            boxes.append( (pos0, u,v,w, rect2d,[zmin,zmax]) )
            volumes.append(volume)
        # select box with minimal volume
        i = np.argmin(volumes)
        (pos0, u,v,w, rect2d,[zmin,zmax]) = boxes[i]
        #
        self.volume = volumes[i]
        # compute the vertices of the minimal enclosing box
        self.vertices = []
        for zw in [zmin,zmax]:
            for corner2d in rect2d.vertices:
                (xu,yv) = corner2d
                corner3d = pos0 + xu*u + yv*v + zw*w
                self.vertices.append( corner3d )
        self.vertices = np.array(self.vertices)
        #
        self.edges = [ [0,1],[1,2],[2,3],[3,0],   # edges around lower face
                       [4,5],[5,6],[6,7],[7,4],   # edges around upper face
                       [0,4],[1,5],[2,6],[3,7] ]  # edges between vertices in lower and upper faces
        # compute center of the box
        self.center = np.sum(self.vertices, axis=0)/len(self.vertices)
        # Next we define a coordinate system whose axes are aligned with the box endges.
        # The x-axis is parallel to the longest edge, the y-axis to the second longest
        # and the z-axis to the shortest
        # Which of the edges a=0->1, b=0->3, or c=0->4 is the longest one?
        a = self.vertices[1]-self.vertices[0]
        b = self.vertices[3]-self.vertices[0]
        c = self.vertices[4]-self.vertices[0]
        vectors = [a,b,c]
        lengths = [la.norm(a), la.norm(b), la.norm(c)]
        iz,iy,ix = np.argsort(lengths)
        xaxis = vectors[ix]/lengths[ix]
        yaxis = vectors[iy]/lengths[iy]
        zaxis = vectors[iz]/lengths[iz]
        # The axes should be right-handed
        signum = np.sign( np.dot(zaxis, np.cross(xaxis,yaxis)) )
        self.axes = np.array([xaxis,yaxis,signum*zaxis])
        
class MoleculeBox(MinimalEnclosingBox):
    def __init__(self, atomlist):
        # convert the list of atoms into a list of 3d points
        nat = len(atomlist)
        pos = XYZ.atomlist2vector(atomlist)
        points = []
        for i in range(0, nat):
            points.append( pos[3*i:3*(i+1)] )
        # QHull needs at least 4 point to construct the initial
        # simplex. 
        if nat < 4:
            # If there are less than 4 atoms,
            # we add the center of mass as an additional point
            masses = AtomicData.atomlist2masses(atomlist)
            com = MolCo.center_of_mass(masses, pos)
            points.append( com )
        if nat < 3:
            # If there are less than 3 atoms we add 
            # an arbitrary point on the x-axis
            points.append( com + np.array([0.0005,0.0,0.0]) )
        if nat < 2:
            # If there is only one atom, we add another arbitrary
            # point on the y-axis
            points.append( com + np.array([0.0,0.0005,0.0]) )
        # We add small random numbers to the input coordinates, so that
        # we get a 3D convex hull even if the molecule is planar
        points = np.array(points) + 0.0001*np.random.rand(len(points),3)
        # find the convex hull using the qhull code
        hull = ConvexHull(points, qhull_options="QbB Qt")

        # call the constructor of the parent class (MinimalEnclosingBox)
        super(MoleculeBox, self).__init__(hull)
        
#### TEST FUNCTIONS ########
        
def test_minimal_enclosing_rectangle():
    points = np.random.rand(20,2)
    
    hull = ConvexHull(points)

    rect = MinimalEnclosingRectangle(hull)
    
    import matplotlib.pyplot as plt
    # plot points
    plt.plot(points[:,0], points[:,1], "o")
    # convex hull
    for simplex in hull.simplices:
        print "simplex = %s" % simplex
        plt.plot(points[simplex,0], points[simplex,1], "k-")
    # minimal enclosing rectangle
    print "rectangle vertices:"
    print rect.vertices
    for edge in rect.edges:
        plt.plot(rect.vertices[edge,0], rect.vertices[edge,1], "r-", lw=2)

    plt.axes().set_aspect('equal', 'datalim')
    plt.show()

def test_minimal_enclosing_box():
    points = np.random.rand(20,3)

    hull = ConvexHull(points)
    box = MinimalEnclosingBox(hull)

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    # plot points
    plt.plot(points[:,0], points[:,1], points[:,2], "o")
    # convex hull
    for simplex in hull.simplices:
        print "simplex = %s" % simplex
        plt.plot(points[simplex,0], points[simplex,1], points[simplex,2], "k-")
    # minimal enclosing box
    print box.vertices
    for edge in box.edges:
        plt.plot(box.vertices[edge,0], box.vertices[edge,1], box.vertices[edge,2], "r-", lw=2)
    # show the center of the box
    plt.plot(box.center[0:1], box.center[1:2], box.center[2:], "x", lw=5, color="green")
    
    
#    plt.axes().set_aspect('equal', 'datalim')
    plt.show()

def pyrene_dimer_analysis(xyz_traj_file, out_file, qmmm_partitioning=None):
    """
    This function analyses the relative displacement between two pyrene units in a dimer.
    'xyz_traj_file' is a xyz-file that contains the geometries for each time step.
    The horizontal (Rx,Ry) and vertical (Rz) displacements are written to a table in 'out_file'.
    If the file contains additional molecules (e.g. MM molecules), the atom indeces belonging
    to the dimer should be specified as a list in qmmm_partioning.
    """
    fh = open(out_file, "w")
    print>>fh, "# TSTEP     R_X / Angstrom    R_Y / Angstrom   R_Z / Angstrom"
    # process one geometry after the other
    for i,atomlist in enumerate(XYZ.read_xyz_it(xyz_traj_file)):
        if qmmm_partitioning != None:
            # select the atoms belonging to the dimer pair
            atomlist = [atomlist[j] for j in qmmm_partitioning]
        nat = len(atomlist)
        assert nat == 52, "Pyrene dimer should contain 52 atoms!"
        # 
        monomer1 = atomlist[:nat/2]
        monomer2 = atomlist[nat/2:]
        # compute the enclosing box for each pyrene unit
        box1 = MoleculeBox(monomer1)
        box2 = MoleculeBox(monomer2)

        # distance between the centers of the enclosing boxes
        vec = box2.center - box1.center
        # decompose distance between centers into a part that is perpendicular to the
        # plane of the 1st monomer and another one that is horizontal to it by projecting
        # onto the z-axis of the box belonging to the 1st monomer
        vec_vertical = np.dot(box1.axes[2], vec) * box1.axes[2]
        ## The x-axis is parallel to the vector passing through atoms C12 and C10
        #xaxis = np.array(atomlist[9][1]) - np.array(atomlist[11][1])
        # In the geometries with periodic boundary conditions the x-axis passes
        # through atoms C3 and C23
        xaxis = np.array(atomlist[2][1]) - np.array(atomlist[22][1])
        # remove z-component from x-axis
        xaxis -= np.dot(box1.axes[2], xaxis) * box1.axes[2]
        xaxis /= la.norm(xaxis)
        ## The y-axis is parallel to the vector passing through atoms C11 and C13
        #yaxis = np.array(atomlist[10][1]) - np.array(atomlist[12][1])
        # In the geometries with periodic boundary conditions the y-axis passes
        # through atoms C14 and C4
        yaxis = np.array(atomlist[13][1]) - np.array(atomlist[3][1])
        yaxis -= np.dot(box1.axes[2], yaxis) * box1.axes[2]
        yaxis /= la.norm(yaxis)

        vec_horizontal = vec - vec_vertical
        Rx = np.dot(xaxis, vec_horizontal) * AtomicData.bohr_to_angs
        Ry = np.dot(yaxis, vec_horizontal) * AtomicData.bohr_to_angs
        Rz = la.norm(vec_vertical) * AtomicData.bohr_to_angs

        print>>fh, "%d  %8.6f  %8.6f  %8.6f" % (i,Rx,Ry,Rz)
    fh.close()
    
if __name__ == "__main__":
    import sys
    test_minimal_enclosing_rectangle()
    test_minimal_enclosing_box()
    #pyrene_dimer_analysis(sys.argv[1], sys.argv[2])
    
    

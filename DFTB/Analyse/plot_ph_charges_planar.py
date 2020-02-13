#!/usr/bin/env python
"""
For planar molecules plot the projection of the geometry into the xy-plane 
and show the distribution of the particle-hole charges.
"""
from DFTB import XYZ
from DFTB import AtomicData
from DFTB.Modeling.porphyrin_flakes import axis_angle2rotation

import numpy as np
import numpy.linalg as la
from matplotlib import pyplot as plt

def rotate_molecule(atomlist, angle, axis):
    Rot = axis_angle2rotation(axis, angle)
    atomlist_rot = []
    for ZA,posA in atomlist:
        posArot = np.dot(Rot, np.array(posA)).tolist()
        atomlist_rot.append((ZA,posArot))
    return atomlist_rot

def shift_molecule(atomlist, shift_vec):
    atomlist_shifted = []
    for ZA,posA in atomlist:
        posAshifted = (np.array(posA) + shift_vec).tolist()
        atomlist_shifted.append((ZA,posAshifted))
    return atomlist_shifted
    
def get_2d_bbox(atomlist):
    pos = XYZ.atomlist2vector(atomlist)
    x, y, z = pos[0::3], pos[1::3], pos[2::3]
    return x.min(), y.min(), x.max(), y.max()

def get_label_position(atomlist):
    xmin,ymin,xmax,ymax = get_2d_bbox(atomlist)
    right_center = (xmax*1.1, (ymin+ymax)/2.0)
    return right_center

def plot_planar_geometry(ax, atomlist, zaxis=2):
    Con = XYZ.connectivity_matrix(atomlist)
    ax.set_aspect('equal', 'datalim')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    # atomic positions
    for A,(ZA,posA) in enumerate(atomlist):
        x,y = posA[0:2]
        atname = AtomicData.atom_names[ZA-1]
        if atname in ["h", "c"]:
            # do not draw the carbon skeleton
            ax.plot([x],[y], "o", color="black", markersize=AtomicData.covalent_radii[atname]*3.0)
        else:
            ax.text(x,y,"%s" % atname.capitalize(), ha='center', va='center')
    # bonds
    Nat = len(atomlist)
    for A in range(0, Nat):
        ZA,posA = atomlist[A]
        xA,yA = posA[0:2]
        connected_atoms = np.where(Con[A,:] > 0)[0]
        for Bcon in connected_atoms:
            ZB,posB = atomlist[Bcon]
            xB,yB,zB = posB[0:3]
            #
            if zB > 0.3:
                ax.annotate("", xy=(xA,yA), xytext=(xB,yB), arrowprops=dict(arrowstyle="wedge", connectionstyle="arc3", color="gray", lw=4))
            elif zB < -0.3:
                ax.plot([xA,xB],[yA,yB], ls="-.", lw=4, color="gray")
            else:
                ax.plot([xA,xB],[yA,yB], ls="-", lw=4, color="gray")

def plot_charges(ax, atomlist, charges, scale=2.0):
    for A,((ZA,posA), qA) in enumerate(zip(atomlist, charges)):
        x,y = posA[0:2]
        if qA > 0.0:
            colour = "blue"
        else:
            colour = "red"
        circ = plt.Circle((x,y), abs(qA)*scale, color=colour, alpha=0.5)
        ax.add_artist(circ)

def plot_ph_distance(ax, atomlist, particle_charges, hole_charges):
    # distance between particle and hole
    p_pos = np.zeros(3)
    h_pos = np.zeros(3)
    for A,(ZA,posA) in enumerate(atomlist):
        posA = np.array(posA)
        p_pos += posA * particle_charges[A]
        h_pos += posA * hole_charges[A]
    p_pos /= float(np.sum(particle_charges))
    h_pos /= float(np.sum(hole_charges))
    ph_dist = la.norm(p_pos - h_pos)
    
    p_x,p_y = p_pos[0:2]
    h_x,h_y = h_pos[0:2]
    ax.text(p_x,p_y,"-", fontsize=30, color="darkblue", fontweight="bold", va="bottom", ha="center")
    ax.text(h_x,h_y,"+", fontsize=30, color="darkred", fontweight="bold", va="bottom", ha="center")
    # distance between + and -
    ax.annotate("", xy=(p_x,p_y), xytext=(h_x,h_y), arrowprops=dict(arrowstyle="|-|", connectionstyle="arc3"))
    angle = np.arctan((h_y-p_y)/(h_x-p_x)) * 180.0/np.pi
    ax.text(0.5*(p_x+h_x),0.5*(p_y+h_y), "$d_{e-h} = %2.1f $ bohr" % ph_dist, fontsize=20, rotation=angle, va="bottom", ha="center")

if __name__ == "__main__":
    import sys
    from os.path import expandvars, expanduser
    from optparse import OptionParser
    usage = "Usage: python %s <xyz file> <file with particle-hole charges>" % sys.argv[0]
    parser = OptionParser(usage)
    parser.add_option("--save_png", dest="save_png", help="save plot to this file [default: %default]", default="")
    parser.add_option("--scale", dest="scale", help="Radii of circles are proportional to the partial charge on the atomic center, <scale> is the proportionality factor. [default: %default]", default=3.0, type=float)
    parser.add_option("--rotate", dest="rotate", help="Rotate molecule around z-axis by this angle [default: %default]", default=0.0, type=float)
    parser.add_option("--densities", dest="densities", help="Select difference densities that should be plotted, list can contain 'hole', 'particle' or 'total' [default: %default]", default="['particle', 'hole', 'total']")

    (opts, args) = parser.parse_args()
    if len(args) < 2:
        print usage
        exit(-1)

    xyz_file = args[0]
    ph_charges_file = args[1]

    atomlist = XYZ.read_xyz(xyz_file)[0]
    atomlist = rotate_molecule(atomlist, opts.rotate*np.pi/180.0, np.array([0.0, 0.0, 1.0]))
#    charges = np.loadtxt(charge_file)
    ph_charges = np.loadtxt(ph_charges_file)

    fig = plt.figure()
    ax = fig.gca()
#    plot_charges(ax, atomlist, charges)
    scale = 1.0/abs(ph_charges).max() * opts.scale
    print "scale = %s" % scale
    xmin,ymin,xmax,ymax = get_2d_bbox(atomlist)
    sep = 1.2*(ymax - ymin)
    nsep = 0
    # distribution of particle
    if "particle" in opts.densities:
        atomlist_p = atomlist
        xl,yl = get_label_position(atomlist_p)
        ax.text(xl, yl, "$\Delta \\rho_p$", fontsize=20, va="center")
        plot_planar_geometry(ax, atomlist_p)
        plot_charges(ax, atomlist_p, ph_charges[:,0], scale=scale)
        nsep += 1
    # distribution of hole
    if "hole" in opts.densities:
        atomlist_h = shift_molecule(atomlist, np.array([0.0, -nsep*sep, 0.0]))
        xl,yl = get_label_position(atomlist_h)
        ax.text(xl, yl, "$\Delta \\rho_h$", fontsize=20, va="center")
        plot_planar_geometry(ax, atomlist_h)
        plot_charges(ax, atomlist_h, ph_charges[:,1], scale=scale)
        nsep += 1
    if "total" in opts.densities:
        # total charges
        atomlist_ph = shift_molecule(atomlist, np.array([0.0, -nsep*sep, 0.0]))
        xl,yl = get_label_position(atomlist_ph)
        ax.text(xl, yl, "$\Delta \\rho$", fontsize=20, va="center")
        plot_planar_geometry(ax, atomlist_ph)
        plot_charges(ax, atomlist_ph, ph_charges[:,0] + ph_charges[:,1], scale=scale)
    
#    plot_ph_distance(ax, atomlist, ph_charges[:,0], ph_charges[:,1])

    if opts.save_png != "":
        plt.savefig(expandvars(expanduser(opts.save_png)))
    plt.show()


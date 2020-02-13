#!/usr/bin/env python

import numpy as np
from numpy.linalg import norm


a0 = 0.52917721067
rcov = {"H": 0.32, "C": 0.77, "N": 0.71, "O": 0.73, "Ru": 1.26}


def bond(coord):
    ij = coord[0]-coord[1]
    return norm(ij)


def dbond(coord):
    ij = coord[0]-coord[1]
    d_di = ij/norm(ij)
    d_dj = -d_di
    return np.array([d_di, d_dj])


def angle(coord):
    ij = coord[0]-coord[1]
    kj = coord[2]-coord[1]
    return np.arccos(np.dot(ij, kj)/(norm(ij)*norm(kj)))*180/np.pi


def dangle(coord):
    a = angle(coord)*np.pi/180
    ij = coord[0]-coord[1]
    kj = coord[2]-coord[1]
    ij0 = ij/norm(ij)
    kj0 = kj/norm(kj)
    d_di = 1/np.sqrt(1-np.cos(a)**2)*(1/norm(ij))*(ij0*np.cos(a)-kj0)*180/np.pi
    d_dk = 1/np.sqrt(1-np.cos(a)**2)*(1/norm(kj))*(kj0*np.cos(a)-ij0)*180/np.pi
    d_dj = -d_di-d_dk
    return np.array([d_di, d_dj, d_dk])


def torsion(coord):
    ji = (coord[1]-coord[0])/norm(coord[1]-coord[0])
    jk = (coord[1]-coord[2])/norm(coord[1]-coord[2])
    lk = (coord[3]-coord[2])/norm(coord[3]-coord[2])
    cross1 = np.cross(ji, jk)
    cross2 = np.cross(jk, lk)
    cross3 = np.cross(cross1, jk)
    dot1 = np.dot(cross1, cross2)
    dot2 = np.dot(cross3, cross2)
    return np.arctan2(dot2, dot1)*180/np.pi


def dtorsion(coord):
    """
    Implementation according to van Schalk et al, J. Mol. Biol. 234, 751 /
    https://salilab.org/modeller/9v6/manual/node436.html
    """
    ij = coord[0]-coord[1]
    kj = coord[2]-coord[1]
    kl = coord[2]-coord[3]
    mj = np.cross(ij, kj)
    nk = np.cross(kj, kl)
    dot1 = np.dot(ij, kj)
    dot2 = np.dot(kl, kj)
    d_di = (norm(kj)/norm(mj)**2)*mj*180/np.pi
    d_dl = -(norm(kj)/norm(nk)**2)*nk*180/np.pi
    d_dj = ((dot1/norm(kj)**2)-1)*d_di-(dot2/norm(kj)**2)*d_dl
    d_dk = ((dot2/norm(kj)**2)-1)*d_dl-(dot1/norm(kj)**2)*d_di
    return np.array([d_di, d_dj, d_dk, d_dl])


def cn_i(r, d, n, m):
    return (1-(r/d)**n)/(1-(r/d)**m)


def dcn_i(r, d, n, m):
    t1 = m*(1-(r/d)**n)*(r/d)**m/(r*(1-(r/d)**m)**2)
    t2 = n*(r/d)**n/(r*(1-(r/d)**m))
    return t1-t2


def get_d(symb, i, j):
    tol = 1.0
    return (rcov[symb[i]]+rcov[symb[j]])*tol/a0


def cn(symb, coord, ndx1, n, m, d, ndx2=[]):
    cn = 0.0
    ds_dr = np.zeros((len(coord), 3))
    for i, coord_i in enumerate(coord):
        # if no set of reference atoms is defined or if i is in the list of reference atoms:
        if (i != ndx1 and not ndx2) or (i != ndx1 and i in ndx2):
            r = bond([coord[ndx1], coord_i])
            # determine d for the current couple of atom types
            if ndx2:  # if reference atoms are defined, d contains only values for atoms in ndx2
                d_i = d[ndx2.index(i)]
            else:  # if no reference atoms are defined, d contains values for all atoms
                d_i = d[i]
            cn += cn_i(r, d_i, n, m)
            dcn = dcn_i(r, d_i, n, m)
            dr = dbond([coord[ndx1], coord_i])
            ds_dr[ndx1] += dcn*dr[0]
            ds_dr[i] = dcn*dr[1]
    print cn
    return cn, ds_dr

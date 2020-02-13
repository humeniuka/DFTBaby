# -*- coding: utf-8 -*-
"""
conversion between cartesian and spherical coordinates
"""

import numpy as np

def arctan3(y,x):
    """translate angles from the range [-pi, pi] to the range [0, 2.0*pi]"""
    phi = np.arctan2(y,x)
    return np.where(phi >= 0.0, phi, 2.0*np.pi+phi)

def cartesian2spherical(C):
    x,y,z = C
    r = np.sqrt(x*x+y*y+z*z)
    #th = np.where(r != 0.0, np.arccos(z/r), 0.0)
    th = 0.0*z
    th[r>0.0] = np.arccos(z[r>0.0]/r[r>0.0])
    """for r=0 the angle th is not well defined, choose th = 0.0"""
    phi = np.where(x != 0.0, \
                    np.where(y != 0.0, arctan3(y,x), (1.0-np.sign(x))/2.0*np.pi), \
                np.pi/2.0+(1.0-np.sign(y))/2.0*np.pi)
    """if x == 0, phi = pi/2 if y is positive and phi=3/2*pi if y is negative
    and if y == 0, phi=0 for x positive and phi=pi for x negative
    """
    return (r,th,phi) 

def spherical2cartesian(C):
    r,th,phi = C
    st = np.sin(th)
    x = r*st*np.cos(phi)
    y = r*st*np.sin(phi)
    z = r*np.cos(th)
    return (x,y,z)

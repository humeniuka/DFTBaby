#!/usr/bin/env python

from numpy import *
from scipy import interpolate
from matplotlib.pyplot import *

from pseudo_atoms import c

r = linspace(0.0, 2.0, 1000)
tck = interpolate.splrep(c.r,c.radial_density,s=0)

rdens_new = interpolate.splev(r,tck,der=0)

xlim((0.0, 2.0))
plot(c.r, c.radial_density, "o")
plot(r, rdens_new, "x")
show()

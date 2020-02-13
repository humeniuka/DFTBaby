#!/usr/bin/env python
"""
plot and compare overlaps and hamilton matrix elements from 
different Slater-Koster tables.
"""
from DFTB import AtomicData
from SKIntegrals import read_slakoformat
import slako_transformations as T

import matplotlib as mpl
label_size = 12
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size

from scipy import interpolate
from matplotlib.pyplot import *
from numpy import *

if len(sys.argv) < 2: 
    print "Usage: %s slako-file1 slako-file2 ..." % sys.argv[0]
    print "slako-files with the following extensions can be read: .skf, .par and .py"
    exit(-1)

def plot_interpolation(dg, Ig, d):
    dmax = max(dg)
    splines_tck = interpolate.splrep(dg,Ig, s=0) 
    dflat = d.ravel()
    # set to zero outside [0.0,dmax] instead of interpolating which gives nonsense
    interp = interpolate.splev(dflat,splines_tck,der=0)
    SorH_interp = reshape(where((dflat<=dmax) & (0.0<=dflat), \
                                    interp, 0.0), d.shape)
    plot(dflat, SorH_interp, ls="-.")

sks = []
for slakofile in sys.argv[1:]:
    sks.append( read_slakoformat(slakofile) )

linestyles = ['-','-.','--',':']
cla()
xlabel(r"distance $d$ between centers / bohr", fontsize=15)
ylabel("overlap", fontsize=15)
title("Overlaps from %s" % slakofile)

for fnr,slako_module in enumerate(sks):
    for (l1,l2,i),olap in slako_module.S.iteritems():
        # do not plot overlaps that are zero
        if all(abs(olap) < 1.0e-10):
            continue
        plot(slako_module.d, olap,
             ls=linestyles[fnr], label="$S_{%s,%s}(%s)(d)$" % \
                 (l1,l2,T.tau2symbol[T.index2tau[i]]), lw=2)
        plot_interpolation(slako_module.d, olap, linspace(0.0, 30.0, 1000))
legend()
show()
#savefig(slakofile + ".S.png")

cla()
xlabel(r"distance $d$ between centers / bohr", fontsize=15)
ylabel("hamiltonian", fontsize=15)
title("Hamilton integrals from %s" % slakofile)
for fnr,slako_module in enumerate(sks):
    for (l1,l2,i),H in slako_module.H.iteritems():
        # do not plot hamilton matrix elements that are zero
        if all(abs(H) < 1.0e-10):
            continue
        plot(slako_module.d, H,
             ls=linestyles[fnr], label="$H_{%s,%s}(%s)(d)$" % \
                 (l1,l2,T.tau2symbol[T.index2tau[i]]), lw=2)
        plot_interpolation(slako_module.d, H, linspace(0.0, 30.0, 1000))
legend()
show()
#savefig(slakofile + ".H.png")

import slako_transformations_dipole as Tdip
angmom_to_xyz = {(0,0): "s",
                 (1,-1): "px", (1,1): "py", (1,0): "pz",
                 (2,-2): "dxy", (2,-1): "dyz", (2,0): "dz2", (2,1): "dzx", (2,2): "dx2y2"}

cla()
xlabel(r"distance $d$ between centers / bohr")
ylabel("Dipole integrals")
title("Dipole integrals %s" % slakofile)

for fnr,slako_module in enumerate(sks):
        if not hasattr(slako_module, "Dipole"):
            continue
        for (l1,l2,i),D in slako_module.Dipole.iteritems():
            lo1,mo1,lM,mM,lo2,mo2 = Tdip.index2tau[i]
            plot(slako_module.d, D, label="$Dipole_{%s,%s}(%s)(d)$" % \
                     (angmom_to_xyz[(lo1,mo1)],angmom_to_xyz[(lo2,mo2)],angmom_to_xyz[(lM,mM)]), lw=2)
            plot_interpolation(slako_module.d, D, linspace(0.0, 30.0, 1000))


legend()
savefig(slakofile + ".D.png")
show()                      


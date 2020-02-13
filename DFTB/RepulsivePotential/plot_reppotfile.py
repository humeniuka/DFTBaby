#!/usr/bin/env python
"""
plot and compare repulsive potentials from
different parametrizations.
"""
from DFTB import AtomicData
from DFTB import utils
from RepulsivePotential import RepulsivePotential, read_repulsive_potential
from matplotlib.pyplot import *
from numpy import *

if len(sys.argv) < 2: 
    print "Usage: %s reppot-file1 reppot-file2 ..." % sys.argv[0]
    print "Plot and compare repulsive potentials."
    print "reppot-files with the following extensions can be read: .par and .py"
    exit(-1)

files = sys.argv[1:]
reppots = []
for reppotfile in files:
    reppots.append( read_repulsive_potential(reppotfile) )

linestyles = ['-','-.','--',':']
cla()
xlabel("distance $d$ between centers [$\AA$]")
ylabel("$V_{rep}(R)$ in [eV]")
title("Repulsive Potential")

for fnr,reppot_module in enumerate(reppots):
    RepPot = RepulsivePotential(reppot_module)
    d = linspace(min(reppot_module.d), 3.0*RepPot.dmax, 1000)
    xlim((0.0,d.max()))
    Vrep_ufunc = frompyfunc(lambda x: RepPot.getVrep(x,deriv=0), 1,1)
    Vrep = Vrep_ufunc(d)
    curve = plot(d*AtomicData.bohr_to_angs, Vrep*AtomicData.hartree_to_eV, \
             label=files[fnr], ls=linestyles[fnr])[0]
    # show where the cut-off radius lies
    annotate("$R_{cut}$", (RepPot.Rcut*AtomicData.bohr_to_angs, 0.0), \
                xytext=(-50, 30+10*fnr), textcoords='offset points',\
                color=curve.get_color(), \
                arrowprops=dict(arrowstyle="->"))
legend()
show()



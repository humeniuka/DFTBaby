#!/usr/bin/env python
"""
Create a tcl script that can be run in vmd as

  vmd -e vmd_input.tcl

to visualize a trajectory. The dftbaby.cfg configuration file is read to determine which 
atoms belong to the QM and MM parts.
"""

from DFTB import utils
from DFTB.DFTB2 import DFTB2
from DFTB.Dynamics.SurfaceHopping import MolecularDynamics

if __name__ == "__main__":
    import sys
    import os.path
    
    usage = "Usage: %s\n" % os.path.basename(sys.argv[0])
    usage += "  creates a tcl script called 'vmd_input.tcl' that allows to visualize the geometries of\n"
    usage += "  a DFTBaby dynamics simulation. The script should be run inside the folder where\n"
    usage += "  the 'dftbaby.cfg' and the 'dynamics.xyz' files are.\n"
    usage += "  To visualize the trajectory, run\n"
    usage += "  \n  vmd -e vmd_input.tcl\n\n"

    # read dftbaby.cfg

    # This wrapper makes the optional parameters of the python function __init__ visible
    # as optional command line argument.
    parser = utils.OptionParserFuncWrapper([DFTB2.__init__, MolecularDynamics.__init__],
                                           usage, section_headers=["SurfaceHopping", "DFTBaby"],
                                           ignore_unknown_options=True)
    # extract optional parameters from command line
    (options,args) = parser.parse_args()

    # QM/MM partitioning
    if options.has_key("qmmm_partitioning"):
        qm_indeces = eval(options["qmmm_partitioning"])
        qm_selection = "{" + "or".join(map(str, [" index %d " % (i-1) for i in qm_indeces])) + "}"
    else:
        qm_selection = "all"
        
    #

    vmd_commands=r"""
# Execute this script as
#   vmd -e vmd_input.tcl

# In the tcl console, type
# > renderMovie
# to produce a sequence of rgb-file snap.####.rgb
#
# This script loads all trajectories with the names dynamics*.xyz into vmd and animates them
foreach dynfile [glob dynamics*.xyz] {
    mol new $dynfile waitfor all
    # VIEW
    axes location off
    display depthcue off
    display projection orthographic
    display resize 1000 500
    display resetview
    # COLORS
    color Display Background white
    # of atoms
    color Name H silver
    color Name C gray
    #
    # REPRESENTATIONS
    # select the QM atoms
    mol selection %s
    # bonds should change dynamically
    mol representation DynamicBonds 1.6 0.2 10
    mol addrep top
    mol representation VDW 0.3 12
    mol addrep top
    
    # hydrogen bonds
    mol selection all
    mol representation HBonds 3.0 20 1
    mol addrep top
}

proc renderMovie {} {
    animate goto start
    set skip 100
    set num [molinfo top get numframes]
    for {set i 0} {$i < $num} {incr i $skip} {
	puts $i
	set filename snap.[format "%%06d" $i].rgb
	render snapshot $filename
	animate goto $i
    }
}

#renderMovie
#exit
    """ % (qm_selection)

    vmd_input = "vmd_input.tcl"
    fh = open(vmd_input, "w")
    fh.write(vmd_commands)
    fh.close()

    print "VMD input script written to %s" % vmd_input

# Execute this script as
#   vmd -e show_trajectories.tcl
#
# This script loads all trajectories into vmd and animates them
foreach dynfile [glob dynamics_*.xyz] {
    mol new $dynfile waitfor all
    # REPRESENTATIONS
#    mol delrep 0 top
#    # bonds should change dynamically
#    mol representation DynamicBonds 1.6 0.2 10
#    mol addrep top
#    mol representation VDW 0.3 12
#    mol addrep top

}

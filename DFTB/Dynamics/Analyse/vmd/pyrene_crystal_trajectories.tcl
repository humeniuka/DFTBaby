# Execute this script as
#   vmd -e pyrene_crystal_trajectories.tcl

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
    # select the 2 central pyrene molecules
    mol selection {index 51 or index 87 or index 235 or index 244 or index 247 or index 250 or index 253 or index 262 or index 283 or index 289 or index 292 or index 345 or index 390 or index 473 or index 476 or index 506 or index 509 or index 654 or index 672 or index 687 or index 690 or index 720 or index 797 or index 800 or index 827 or index 830 or index 14 or index 15 or index 24 or index 25 or index 105 or index 111 or index 116 or index 117 or index 127 or index 322 or index 323 or index 333 or index 334 or index 418 or index 433 or index 520 or index 523 or index 524 or index 525 or index 526 or index 529 or index 536 or index 538 or index 539 or index 736 or index 748}
    # bonds should change dynamically
    mol representation DynamicBonds 1.6 0.2 10
    mol addrep top
    mol representation VDW 0.3 12
    mol addrep top
}

proc renderMovie {} {
    animate goto start
    set skip 100
    set num [molinfo top get numframes]
    for {set i 0} {$i < $num} {incr i $skip} {
	puts $i
	set filename snap.[format "%06d" $i].rgb
	render snapshot $filename
	animate goto $i
    }
}

#renderMovie
#exit

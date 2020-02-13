# Execute this script as
#   vmd -e show_unit_cell.tcl
#

mol new tube_unit_cell.ff type xyz


# VIEW
axes location off
display projection orthographic
# COLORS
color Display Background black
# REPRESENTATIONS
# Zinc atoms
mol selection "name Zn"
mol representation VDW 1.0 12
mol addrep top
# replica atoms
mol selection "all"
mol representation Points
mol addrep top

# PERIODIC BOUNDARY CONDITIONS
pbc set {0.0 0.0 13.9975371   90.0 90.0 90.0}

mol showperiodic top 1 z
mol numperiodic top 1 5
mol showperiodic top 2 z
mol numperiodic top 2 5

# PLOT SPIRALS
set radius 26.0
set height 28.0
set pi [expr {atan(1) * 4}]
set num_points 4000
set num_windings 4
# center of molecule
set x0 0.0
set y0 -26.0
set z0 -2.92
set phi_starts { -30 60 150 240 }
set colors {blue orange red green}
set j 0
foreach phi_start $phi_starts {
    graphics top color [lindex $colors $j]
    set phi0 [expr {$phi_start * $pi/180.0}]

    for {set i 0} {$i < $num_points} {incr i} {
	# polar angle
	set phi [expr {- $num_windings * 2 * $pi / $num_points * $i}]
	# cartesian coordinates
	set x [expr { $radius * cos($phi - $phi0) + $x0 }]
	set y [expr { $radius * sin($phi - $phi0) + $y0 }]
	set z [expr { - $height / (2 * $pi) * $phi + $z0 }]
	
	graphics top point "$x $y $z"
    }
    incr j
}

# diameter
label add Bonds "0/310" "0/1830"
# height of one winding
color white
graphics top line {15.25 -47.75 -2.92} {15.25 -47.75 25.08} width 4 style dashed
graphics top text {15.25 -47.75 11.08} "height = 28 Ang" 

# C=O --- H-O  hydrogen bond distance
label add Bonds "0/358" "0/599"
# C=O --- H-O  hydrogen bond angle
label add Angles "0/358" "0/600" "0/599"
# Zn - O distance
label add Bonds "0/310" "0/447"
#
#                                                This Model       Holzwarth Model
#                                                =================================
#
#  number of stack per closed ring                  20             20 +- 2
#  Zn-Zn distance in one stack                     7.0 Ang        6.8 +- 0.2 Ang
#  tube diameter (Zn to opposite Zn distance)    53.13 Ang           54 Ang              
#  H-bond bond length (donor-acceptor)            3.67 Ang            ?
#  H-bond length (hydrogen-acceptor)              2.03 Ang           2.2 Ang
#  H-bond angle                                   159.3            139-153 deg
#

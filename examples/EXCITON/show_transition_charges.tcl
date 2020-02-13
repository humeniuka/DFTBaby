#
# Execute this script as
#    vmd -e show_transition_charges.tcl  -args  aggregate.ff  chromophores.chromo
#
# Purpose:
# --------
# The atoms are coloured by their transition charges and the electric and magnetic transition
# dipoles for each chromophore (each geometry in chromophores.chromo) are computed and displayed
# as a blue and red vectors.
# WARNING: If the magnetic dipoles are defined in the local molecular frame and were simply
# copied between different chromophores, they point into the wrong direction. You need to
# use the option --rotate_dipoles=1 of `copy_transition_charges.py` to rotate each dipole
# according to the orientation of each chromophore.
#
# Arguments:
# ----------
# aggregate.ff        :   .xyz or .ff file with geometry of aggregate
# chromophores.chromo :   .chromo file with transition charges for each chromophore
#

# 
set ff_file [lindex $argv 0]
set chromo_file [lindex $argv 1]
# excited (monomer) state
set state 1

mol new $ff_file type xyz

# VIEW
axes location off
display projection orthographic
# COLORS
color Display Background black

proc load_transition_charges {filename {state 1}} {
    # read transition charges from file
    set fh [open $filename r]
    # enumerate chromophores
    set nchromo 0  
    while {[gets $fh nat] >= 0} {
	# compute center of mass
	set center {0.0 0.0 0.0}
	set total_mass 0.0
	# and electric transition dipole
	set elec_dipole {0.0 0.0 0.0}
	# and magnetic transition dipole
	set magn_dipole {0.0 0.0 0.0}
	
	# 1st line: number of atoms nat
	# 2nd line: comment line with excitation energies
	set comment [gets $fh]
	# lines 3-nat+3: atom type, coordinates, index into main list of atoms and transition charges  q mx my mz
	for {set i 0} {$i < $nat} {incr i} {
	    set line [gets $fh]
	    set index [lindex $line 4]
	    set q [lindex $line [expr "4 + 4 * ( $state - 1) + 1"]]
	    set mx [lindex $line [expr "4 + 4 * ( $state - 1) + 2"]]
	    set my [lindex $line [expr "4 + 4 * ( $state - 1) + 3"]]
	    set mz [lindex $line [expr "4 + 4 * ( $state - 1) + 4"]]
	    puts "index=$index charge=$q  mx=$mx my=$my mz=$mz"
	    # write the transition charge to the charge field of the atom
	    set atom [atomselect top "index $index"]
	    $atom set charge $q
	    set pos "[$atom get x] [$atom get y] [$atom get z]"
	    set mass [$atom get mass]
	    # weight position by mass
	    set mass_pos [vecscale $mass $pos]
	    set total_mass [expr "$total_mass + $mass"]
	    set center [vecadd $center $mass_pos]
	    # weight position by charge
	    set charge_pos [vecscale $q $pos]
	    set elec_dipole [vecadd $elec_dipole $charge_pos]
	    #
	    set magn_dipole [vecadd $magn_dipole "$mx $my $mz"]
     
	}
	# convert dipoles from angstrom to atomic units
	set b2a 0.529177249
	set elec_dipole [vecscale [expr "1.0 / $b2a"] $elec_dipole]
	set magn_dipole [vecscale [expr "1.0 / $b2a"] $magn_dipole]
	
	# display the transition dipole moment for each chromophore
	set center [vecscale [expr "1.0/$total_mass"] $center]
	puts "center of mass of chromophore $nchromo (in Ang) = $center"
	# draw vector for electric transition dipole
	puts "electric dipole (in au): $elec_dipole"
	set molid 0
	graphics $molid color blue
	vmd_draw_vector $molid $center $elec_dipole scale 4 radius 0.3
	# draw vector for magnetic transition dipole
	#puts "magnetic dipole (in au): $magn_dipole"
	graphics $molid color red
	vmd_draw_vector $molid $center $magn_dipole scale 4 radius 0.3
	
	incr nchromo
	# 
    }
    close $fh
}

#############################################################################
# The code for the function `vmd_draw_vector` is taken from vmd_draw_lib.tcl
#
# vmd draw extension procedures
#
# version 1.0
# Copyright (c) 2006 Axel Kohlmeyer <akohlmey@cmm.chem.upenn.edu>
#
#
# define a vector drawing function. contrary to the arrow example
# from the users guide, the vector is defined by it's origin,
# a diretion at this point, and a scaling factor.
# the function returns a list of the graphics ids for easy deletion.
proc vmd_draw_vector {args} {
    set usage {"vmd_draw_vector molid {x1 y1 z1} {x2 y2 z2} [scale <s>] [resolution <res>] [radius <r>] [filled <yes/no>]"}
    # defaults
    set scale 0.8
    set res 6
    set radius 0.2
    set filled yes

    if {[llength $args] < 3} {
	error "wrong # args: should be $usage"
    }
    set mol    [lindex $args 0]
    set center [lindex $args 1]
    set vector [lindex $args 2]
    if {[llength $center] != 3 || [llength $vector] != 3} {
	error "wrong type of args: should be $usage"
    }

    foreach {flag value} [lrange $args 3 end] {
	switch -glob $flag {
	    scale  {set scale  $value}
	    res*   {set res    $value}
	    rad*   {set radius $value}
	    fill*  {set filled $value}
	    default {error "unknown option '$flag': should be $usage" }
	}
    }

    set vechalf [vecscale [expr $scale * 0.5] $vector]
    return [list \
		[graphics $mol cylinder [vecsub $center $vechalf] \
		     [vecadd $center [vecscale 0.7 $vechalf]] \
		     radius $radius resolution $res filled $filled] \
		[graphics $mol cone [vecadd $center [vecscale 0.7 $vechalf]] \
		     [vecadd $center $vechalf] radius [expr $radius * 1.7] \
		     resolution $res]]
}
#################################################################################

puts "load transition charges"
load_transition_charges $chromo_file 1

# REPRESENTATIONS
# color atoms by transition charges
mol selection "all"
mol color Charge
mol representation VDW 0.2 12
mol addrep top

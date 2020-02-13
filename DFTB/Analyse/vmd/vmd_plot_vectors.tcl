#
# run this script in VMD as
#
#  vmd -e vmd_plot_vectors.tcl -args geometry.xyz  vectors.dat  scale
#
# to plot non-adiabatic coupling vectors as arrows centered at each atom.
#
# Arguments:
#   - `geometry.xyz` should contain the molecular geometry in XYZ-format
#   - `vectors.dat` should contain at matrix of shape (nr.atoms,3) with the
#     vectors for each atom
#   - `scale` is a floating point number by which the length of all vectors is
#     scaled

proc align_axes {sel} {
    #
    # aligns the principal axes of inertia of the molecule
    # with the global axes
    #
    # Returns
    # =======
    #   m4     :  the 4x4 transformation matrix
    
    # measure tensor of inertia (a 3x3 matrix)
    set m3 [lindex [measure inertia $sel] 1]
    # build 4x4 transformation matrix
    set v0 "[lindex $m3 0] 0"
    set v1 "[lindex $m3 1] 0"
    set v2 "[lindex $m3 2] 0"
    set v3 {0 0 0  1}
    set m4 [list $v0 $v1 $v2 $v3]
    puts "Transformation matrix"
    puts $m4
    # rotate molecule so that the principal axes of inertia
    # coincide with the x,y and z axes
    $sel move $m4

    return $m4
}


axes location off
display projection orthographic
display depthcue off
color Display Background white

# input files
set xyz_file [lindex $argv 0]
set vector_file [lindex $argv 1]
# scaling of vectors
set scale [lindex $argv 2]

# load geometry
mol new $xyz_file

# choose representation
mol delrep 0 top
mol representation CPK 0.2 0.1 12 12
mol addrep top


set atoms [atomselect top "all"]

# rotate molecule, so that it is aligned with the global axes
set m4 [align_axes $atoms]
# top view on molecule
rotate y by -90

set vector_fh [open $vector_file]
set i 0
foreach coord [$atoms get {x y z}] name [$atoms get {name}] {
    # read next line from file which is not a comment
    set line "#"
    while {[regexp "#" $line] == 1} {
	# skip comment
	set line [gets $vector_fh]
    }
    set vec $line
    set pos0 {}
    set pos1 {}
    set tip {}
    # length**2 of vector
    set len2 0.0
    puts $vec

    # rotate vectors, too
    set vec [vectrans $m4 $vec]
    # tail of arrow
    foreach x $coord y $vec {
	lappend pos0 [expr {$x - 0.5 * $scale * $y}]
	lappend pos1 [expr {$x + 0.2 * $scale * $y}]
	set len2 [expr {$len2 + $y * $y}]
    }
    set len [expr { sqrt( $len2 ) }]

    # Tip of arrow
    foreach x $coord y $vec {
	lappend tip  [expr {$x + ( 0.5 * $scale * $y )}]
    }
    
    set cone_radius [expr {0.2 * sqrt( $len2 ) * $scale }]
    set cylinder_radius [expr {0.5 * $cone_radius}]

    # atom label
#    label add Atoms 0/$i
#    label textformat Atoms $i {%e}
    
    # draw arrow made of a cylinder and a cone as the tip
    graphics 0 color red
    graphics 0 cylinder $pos0 $pos1 radius $cylinder_radius
    graphics 0 cone $pos1 $tip  radius $cone_radius resolution  10

    incr i
}

# show scale factor in the upper left corner
set max [lindex [measure minmax $atoms] 1]
set ymax [lindex $max 1]
set ymax [expr { $ymax * 1.3 }]
set zmax [lindex $max 2]

graphics 0 color black
graphics 0 text "0 $ymax $zmax" "x$scale" size 2 thickness 3

## render image and save it to a file
#render TachyonInternal vectors.tga

#exit


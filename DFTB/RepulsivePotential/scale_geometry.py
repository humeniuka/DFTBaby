#!/usr/bin/env python
"""
generate different geometries by scaling the position vector of each atom by a constant 
factor.
"""
from numpy import linspace
from DFTB import utils, XYZ


class ScaledGeometries:
    def __init__(self, xyz_in_file, xyz_scaled_file):
        self.xyz_in_file = xyz_in_file
        self.xyz_scaled_file = xyz_scaled_file
    def generate_scaled_geometries(self, scale_min=0.1, scale_max=1.9, Nscale=100):
        """
        Parameters:
        ===========
        scale_min: scale factor for smallest structure
        scale_max: scale factor for largest structure
        Nscale: number of points in the interval [scale_min, scale_max]
        """
        geometry = XYZ.read_xyz(self.xyz_in_file)[0]
        scaled_geoms = []
        for s in linspace(scale_min, scale_max, Nscale):
            scaled_atomlist = []
            for (Zi,(xi,yi,zi)) in geometry:
                scaled_atomlist.append( (Zi, (xi*s,yi*s,zi*s)) )
            scaled_geoms.append(scaled_atomlist)
        XYZ.write_xyz(self.xyz_scaled_file, scaled_geoms, title="scaled geometries")

if __name__ == "__main__":
    import sys
    usage = "python %s <geometry .xyz-file> <scaled geometries .xyz-file>\n\n" % sys.argv[0]
    usage += "generate different geometries by scaling the position vector of each atom by a constant factor."
    if len(sys.argv) < 3:
        print usage
        exit(-1)
    SG = ScaledGeometries(sys.argv[1], sys.argv[2])
    parser = utils.OptionParserFuncWrapper(SG.generate_scaled_geometries,usage)
    (options, args) = parser.parse_args()
    SG.generate_scaled_geometries(**options)

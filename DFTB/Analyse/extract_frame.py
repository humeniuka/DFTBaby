from DFTB import XYZ
import sys

frame=433
geometries = XYZ.read_xyz("dynamics.xyz")

XYZ.write_xyz("frame_%d.xyz" % frame, [geometries[433]])

# Blender script for plotting xyz data
# requires Matthias' pylib

import sys, copy
sys.path = ["/usr/lib/python3/dist-packages"] + sys.path

import pylib.blender.path as bpath
import pylib.blender.molecule as bmol
import pylib.molecule as molecule
import pylib.blender.materials as bmat
import imp, glob

# USER INTERFACE
import bpy

# ImportHelper is a helper class, defines filename and
# invoke() function which calls the file selector.
from bpy_extras.io_utils import ImportHelper
from bpy.props import StringProperty, IntProperty
from bpy.types import Operator


class ImportXYZData(Operator, ImportHelper):
    """import and plot geometries in .xyz format"""
    bl_idname = "import_xyz.xyz_data"  
    bl_label = "Import XYZ structure"

    # ImportHelper mixin class uses this
    filename_ext = ".xyz"

    filter_glob = StringProperty(
            default="*.xyz",
            options={'HIDDEN'},
            )

    def execute(self, context):
        xyz = molecule.read_xyz(self.filepath)
        print("xyz file has been read")
        bpy.context.scene.frame_set(0)
        mol = bmol.Molecule(xyz)
        print("mol is ready")
        mol.generate().create()
        print("molecule has been created")
        return {'FINISHED'}

def menu_func_import(self, context):
    self.layout.operator(ImportXYZData.bl_idname, text="XYZ Import Operator")


def register():
    bpy.utils.register_class(ImportXYZData)
    bpy.types.INFO_MT_file_import.append(menu_func_import)


def unregister():
    bpy.utils.unregister_class(ImportXYZData)
    bpy.types.INFO_MT_file_import.remove(menu_func_import)



if __name__ == "__main__":
    register()

    bpy.ops.import_xyz.xyz_data('INVOKE_DEFAULT')


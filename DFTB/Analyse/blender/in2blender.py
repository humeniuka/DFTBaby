#!/usr/bin/env python
"""
This blender script visualizes the initial geometry of a trajectory 
in the *.in format that is used by Jens' program and plots the initial
velocities as arrows. The script only works with blender-2.49b. 
You habe to load it in a text-window in blender and type Alt-P to start it.
"""

from Blender import *
from Blender import Mathutils
from Blender.Mathutils import Vector
from Blender.BGL import *
from Blender.Draw import *
import bpy
import math

# Global data for Graphical User Interface

# Global data for molecule visualization
isothreshold=0.10
spheresubdivisions=5
covalentradii={"c":0.77,"h":0.37,"n":0.77,"o":0.77,"zn": 1.50,"au":1.50}
atomicradii={"h":0.4,"c":0.7,"n":0.7,"o":0.7,"zn":1.50,"au":2.00}
stickradius=0.20

# Materials for atoms, bonds and orbital (density) lobes
matC = Material.New('c')                                
matC.rgbCol = [0.0, 1.0, 0.0]          
matC.setAlpha(0.8)                     
matC.emit = 0.2
matC.setHardness(400)
matC.setSpec(2.0)

matH = Material.New('h')                                
matH.rgbCol = [0.6, 0.6, 0.6]          
matH.setAlpha(0.8)                     
matH.emit = 0.2 
matH.setHardness(400)
matH.setSpec(2.0)

matN = Material.New('n')                                
matN.rgbCol = [0.2, 0.1, 0.9]          
matN.setAlpha(0.8)                     
matN.emit = 0.2  
matN.setHardness(400)
matN.setSpec(2.0)

matO = Material.New('o')                                
matO.rgbCol = [1.0, 0.2, 0.1]          
matO.setAlpha(0.8)                     
matO.emit = 0.2 
matO.setHardness(400)
matO.setSpec(2.0)

matS = Material.New('s')                                
matS.rgbCol = [0.9, 0.9, 0.0]          
matS.setAlpha(0.8)                     
matS.emit = 0.2 
matS.setHardness(400)
matS.setSpec(2.0)

matZn = Material.New('zn')                                
matZn.rgbCol = [0.1,0.1, 0.1]          
matZn.setAlpha(0.8)                     
matZn.emit = 0.2 
matZn.setHardness(400)
matZn.setSpec(2.0)

matAu = Material.New('au')                                
matAu.rgbCol = [228.0/255.0, 202.0/255.0, 66.0/255.0]          
matAu.setAlpha(0.8)                     
matAu.emit = 0.2 
matAu.setHardness(400)
matAu.setSpec(2.0)

matstick = Material.New('stick')                                
matstick.rgbCol = [0.8, 0.2, 0.2]          
matstick.setAlpha(0.8)                     
matstick.emit = 0.2                 
matstick.setHardness(400)
matstick.setSpec(2.0)

matpos = Material.New('positivelobe')                                
matpos.rgbCol = [1.0, 0.0, 0.0]          
matpos.setAlpha(0.8)                     
matpos.emit = 0.2                 
matpos.setHardness(400)
matpos.setSpec(2.0)

matneg = Material.New('positivelobe')                                
matneg.rgbCol = [0.0, 0.0, 1.0]          
matneg.setAlpha(0.8)                     
matneg.emit = 0.2                 
matneg.setHardness(400)
matneg.setSpec(2.0)

matArrow = Material.New("arrow")
matArrow.rgbCol = [0.5, 0.0, 1.0]          
matArrow.setAlpha(0.8)                     
matArrow.emit = 0.2                 
matArrow.setHardness(400)
matArrow.setSpec(2.0)

arrowradius=0.1

materials={"h":matH,"c":matC,"n":matN,"o":matO, "s": matS, "zn":matZn, "au":matAu,"stick":matstick,"positivelobe":matpos,"negativelobe":matneg, "arrow": matArrow}
class xyz:
	def __init__(self,filename):
		self.filename=filename
		self.coord=[]
                self.vectors=[]
		self.atomtypes=[]
	def readxyz(self,):
		fh=open(self.filename,"r")
		while 1:
			line = fh.readline().strip()
			if line == "":
                                break
                        nat = int(line.strip())
#			title = fh.readline()
			atoms = []
                        # read coordinates
			for i in range(0,nat):
				line = fh.readline()
				words = line.split()
                                # convert to Angstrom
				x,y,z = map(lambda f: 0.529177249*float(f),words[1:])
				atname = words[0].lower()
				self.atomtypes.append(atname)
				self.coord.append([x,y,z])
                        # read vectors
			for i in range(0,nat):
				line = fh.readline()
				words = line.split()
                                # convert to Angstrom
				vx,vy,vz = map(lambda f: 0.529177249*float(f),words)
				self.vectors.append([vx,vy,vz])
		fh.close()


class cube2blender:
	def __init__(self,xyz):
		self.xyz=xyz
		self.scene=Scene.GetCurrent()

	def blenderstructure(self,):
                nat = len(self.xyz.atomtypes)
                # create atoms
	        for i in range(0,nat):
                        loc=Mathutils.Vector(self.xyz.coord[i])
	                me = Mesh.Primitives.Icosphere(spheresubdivisions,atomicradii[self.xyz.atomtypes[i]])
			me.materials=[materials[self.xyz.atomtypes[i]]]
		        for face in me.faces:
				face.smooth=True
	                obj=self.scene.objects.new(me,'Atom')
	                obj.setLocation(loc)
                # form bonds between atoms
	        for i in range(0,nat):
	                for j in range(i+1,nat):
	                        vec1=Mathutils.Vector(self.xyz.coord[i])
	                        vec2=Mathutils.Vector(self.xyz.coord[j])
	                        vec=vec2-vec1
	                        distcovalent=covalentradii[self.xyz.atomtypes[i]]+covalentradii[self.xyz.atomtypes[j]]
                                print "vec.length = %s  distcovalent = %s" % (vec.length, distcovalent)
	                        if (vec.length-distcovalent)<=0.10*distcovalent:
	                                me=Mesh.Primitives.Tube(32,stickradius,vec.length)
		        		for face in me.faces:
						face.smooth=True
	                                obj=self.scene.objects.new(me,'Bond')
	                                axis=Mathutils.CrossVecs(Vector([0,0,1]),vec)
	                                angle=Mathutils.AngleBetweenVecs(Vector([0,0,1]),vec)
	                                rotmat=Mathutils.RotationMatrix(angle,4,"R",axis)
	                                obj.setMatrix(obj.matrix*rotmat)
	                                obj.setLocation((vec1+vec2)*0.5)
		# vectors
                scale = 1.0 #1000.0 #1.0 #1000.0
                for i in range(0, nat):
                        loc=Mathutils.Vector(self.xyz.coord[i])
                        vec=Mathutils.Vector(self.xyz.vectors[i])*scale
                        # arrow tail
                        me=Mesh.Primitives.Tube(32,arrowradius,vec.length)
			me.materials=[materials["arrow"]]
                        for face in me.faces:
                                face.smooth=True
                        obj=self.scene.objects.new(me,"Arrow-Tail")
                        axis=Mathutils.CrossVecs(Vector([0,0,1]),vec)
                        angle=Mathutils.AngleBetweenVecs(Vector([0,0,1]),vec)
                        rotmat=Mathutils.RotationMatrix(angle,4,"R",axis)
                        obj.setMatrix(obj.matrix*rotmat)
                        obj.setLocation(loc+0.5*vec)
                        # arrow head
                        me=Mesh.Primitives.Cone(32,2*arrowradius,0.5)
			me.materials=[materials["arrow"]]
                        for face in me.faces:
                                face.smooth=True
                        obj=self.scene.objects.new(me,"Arrow-Head")
                        axis=Mathutils.CrossVecs(Vector([0,0,1]),vec)
                        angle=Mathutils.AngleBetweenVecs(Vector([0,0,1]),vec)
                        rotmat=Mathutils.RotationMatrix(angle+180.0,4,"R",axis)
                        obj.setMatrix(obj.matrix*rotmat)
                        obj.setLocation(loc+vec)
                        
EVENT_none=0
EVENT_Button_Import=1
EVENT_Button_Cancel=2
EVENT_Button_Browse=3
#inputfile=Create("/tmp/dyn.in")
inputfile=Create("/tmp/mulliken_dipoles.dat")

class GUI:
	def __init__(self,):
		pass		
	
	def getFilename_callback(self,filename):
		inputfile.val=filename
	

	def gui_draw(self):
		global EVENT_Button_Import
		global EVENT_Button_Cancel
		global EVENT_Button_Browse
		global EVENT_Button_Browse2
		global isovalue
		global sliderMax
		global sliderMin

		glClearColor(0.5, 0.5, 0.5, 0.0)
		glClear(BGL.GL_COLOR_BUFFER_BIT)
		glColor3f(0.0, 0.0, 0.0)

		glColor3f(1.0, 1.0, 1.0)
		glRasterPos2i(180, 340)
		Text(".IN File Importer", "large")


		String("File path: ", EVENT_none,40, 300, 360, 20,inputfile.val,399, ".in file path")

		Button("Import", EVENT_Button_Import, 110, 260, 70, 24, "")
		Button("Cancel", EVENT_Button_Cancel, 320, 260, 70, 24, "")
		Button("Browse", EVENT_Button_Browse, 400, 300, 70, 20, "")
	
	def event(self,event, value):
		if (event == ESCKEY or event == QKEY) and not value:
			Exit()

	def b_event(self,event):
		global EVENT_Button_Import
		global EVENT_Button_Cancel
		global EVENT_Button_Browse
		if event == 0: 
                        pass
		elif event == EVENT_Button_Import:
			self.xyzobject=xyz(inputfile.val)
			self.xyzobject.readxyz()
			self.cube2blenderobject=cube2blender(self.xyzobject)
			self.cube2blenderobject.blenderstructure()
		elif event == EVENT_Button_Cancel:
			Exit()
		elif event == EVENT_Button_Browse:
			Window.FileSelector(self.getFilename_callback, "Select xyz file")
		Draw()
		
		
guiobject=GUI()
Register(guiobject.gui_draw, guiobject.event, guiobject.b_event)

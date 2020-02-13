#!/usr/bin/env python

from Blender import *
from Blender import Mathutils
from Blender.Mathutils import Vector
from Blender.BGL import *
from Blender.Draw import *
from Blender.Window import *
import bpy
import math
import os
import vtk

pwd = os.getcwd()

# Global data for Graphical User Interface
width = int(680)
height = int(3*width/4.)

EVENT_none=0
EVENT_Button_Import=1
EVENT_Button_Cancel=2
EVENT_Button_Browse=3
EVENT_Button_Plot=4
EVENT_Slider_Isovalue=5
inputfile=Create(pwd + "/DENSE/ed_0000.cube")
sliderMin=Create(0.0)
sliderMax=Create(1.0)
isovalue=Create(0.0)

# Global data for molecule visualization
isothreshold=0.10
spheresubdivisions=5
covalentradii={"6":0.77/0.529,"1":0.37/0.529,"7":0.77/0.529,"8":0.77/0.529, "16": 1.98, "30":1.35/0.529 ,"79":1.50/0.529}
atomicradii={"1":0.7,"6":1.50,"7":1.50,"8":1.50, "16": 1.50, "30":2.00, "79":2.00}
stickradius=0.30

# Materials for atoms, bonds and orbital (density) lobes
matC = Material.New('C')                                
matC.rgbCol = [0.0, 1.0, 0.0]          
matC.setAlpha(0.8)                     
matC.emit = 0.2
matC.setHardness(400)
matC.setSpec(2.0)

matH = Material.New('H')                                
matH.rgbCol = [0.6, 0.6, 0.6]          
matH.setAlpha(0.8)                     
matH.emit = 0.2 
matH.setHardness(400)
matH.setSpec(2.0)

matN = Material.New('N')                                
matN.rgbCol = [0.2, 0.1, 0.9]          
matN.setAlpha(0.8)                     
matN.emit = 0.2  
matN.setHardness(400)
matN.setSpec(2.0)

matO = Material.New('O')                                
matO.rgbCol = [1.0, 0.2, 0.1]          
matO.setAlpha(0.8)                     
matO.emit = 0.2 
matO.setHardness(400)
matO.setSpec(2.0)

matS = Material.New('S')                                
matS.rgbCol = [0.9, 0.9, 0.0]          
matS.setAlpha(0.8)                     
matS.emit = 0.2 
matS.setHardness(400)
matS.setSpec(2.0)

matZn = Material.New('Zn')                                
matZn.rgbCol = [0.1,0.1, 0.1]          
matZn.setAlpha(0.8)                     
matZn.emit = 0.2 
matZn.setHardness(400)
matZn.setSpec(2.0)

matAu = Material.New('Au')                                
matAu.rgbCol = [228.0/255.0, 202.0/255.0, 66.0/255.0]          
matAu.setAlpha(0.8)                     
matAu.emit = 0.2 
matAu.setHardness(400)
matAu.setSpec(2.0)

matDefault = Material.New('Default')
matDefault.rgbCol = [228.0/255.0, 202.0/255.0, 66.0/255.0]          
matDefault.setAlpha(0.8)                     
matDefault.emit = 0.2 
matDefault.setHardness(400)
matDefault.setSpec(2.0)

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

materials={"1":matH,"6":matC,"7":matN,"8":matO, "16": matS, "30":matZn, "79":matAu,"stick":matstick,"positivelobe":matpos,"negativelobe":matneg, "default": matDefault}

class cube:
	def __init__(self,filename):
		self.filename=filename
		self.coord=[]
		self.atomtypes=[]
	def readCube(self,):
		file=open(self.filename,"r")
		self.header=file.readline()
		self.header=self.header+" "+file.readline()
		line=file.readline().split()
	        self.natom=abs(int(line[0]))                                            # number of atoms
        	self.origin=[float(line[1]),float(line[2]),float(line[3])]        # origin of the cube file
        	line=file.readline().split()
        	self.nx=int(line[0])                                                # number of points in the x-direction
        	self.ivec=[float(line[1]),float(line[2]),float(line[3])]
		self.dx=math.sqrt(self.ivec[0]**2+self.ivec[1]**2+self.ivec[2]**2)
        	line=file.readline().split()
        	self.ny=int(line[0])                                                # number of points in the y-direction
        	self.jvec=[float(line[1]),float(line[2]),float(line[3])]
		self.dy=math.sqrt(self.jvec[0]**2+self.jvec[1]**2+self.jvec[2]**2)
        	line=file.readline().split()
        	self.nz=int(line[0])                                                # number of points in the z-direction
        	self.kvec=[float(line[1]),float(line[2]),float(line[3])]
		self.dz=math.sqrt(self.kvec[0]**2+self.kvec[1]**2+self.kvec[2]**2)
        	for i in range(self.natom):
                	line=file.readline()
                	splitline=line.split()
                        self.atomtypes.append(splitline[0])
                        self.coord.append([float(splitline[2]),float(splitline[3]),float(splitline[4])])
        	self.gridpoints=[]
        	for i in range(self.nx):
                	for j in range(self.ny):
                        	for k in range(self.nz):
                                	x=self.origin[0]+i*self.ivec[0]+j*self.jvec[0]+k*self.kvec[0]
                	             	y=self.origin[1]+i*self.ivec[1]+j*self.jvec[1]+k*self.kvec[1]
                                	z=self.origin[2]+i*self.ivec[2]+j*self.jvec[2]+k*self.kvec[2]
                                	self.gridpoints.append([x,y,z])

		self.data1d=[]
        	for line in file:
                	splitline=line.split()
                	self.data1d.extend(map(float,splitline))
		file.close()
		self.data3d=[[[0 for i in range(self.nx)] for j in range(self.ny)] for k in range(self.nz)]
		ind=0
		for i in range(self.nx):
			for j in range(self.ny):
				for k in range(self.nz):
					self.data3d[k][j][i]=self.data1d[ind]
					ind+=1

	def cube2vtk(self,filename):
		# create temporarily VTK file equivalent to the .cube file
		vtkfile=open(filename,"w")
		vtkfile.write("# vtk DataFile Version 2.0\n")
		vtkfile.write("vtk representation of a gaussian cube file\n")
		vtkfile.write("ASCII\n")
		vtkfile.write("\n")
		vtkfile.write("DATASET STRUCTURED_POINTS\n")
		vtkfile.write("DIMENSIONS "+str(self.nx)+" "+str(self.ny)+" "+str(self.nz)+"\n")
		vtkfile.write("ORIGIN "+str(self.origin[0])+" "+str(self.origin[1])+" "+str(self.origin[2])+"\n")
		vtkfile.write("\n")
		vtkfile.write("SPACING "+str(self.dx)+" "+str(self.dy)+" "+str(self.dz)+"\n")
		vtkfile.write("POINT_DATA "+str(self.nx*self.ny*self.nz)+"\n")
		vtkfile.write("SCALARS scalars float\n")
		vtkfile.write("LOOKUP_TABLE defaults\n")
		for i in range(len(self.data)):
		       vtkfile.write(str(self.data[i])+"\n")
		vtkfile.close()

class cube2blender:
	def __init__(self,cube):
		self.cube=cube
		self.scene=Scene.GetCurrent()
		for ob in self.scene.objects:
			if ob.name=='Plane':
				self.scene.unlink(ob)
		# insert the molecule at the current cursor position
		self.cursorXYZ = GetCursorPos()
		#

	def blenderstructure(self,):
	        for i in range(len(self.cube.atomtypes)):
	                me = Mesh.Primitives.Icosphere(spheresubdivisions,atomicradii.get(self.cube.atomtypes[i], 2.0))
			me.materials=[materials.get(self.cube.atomtypes[i], materials["default"])]
		        for face in me.faces:
				face.smooth=True
	                obj=self.scene.objects.new(me,'Mesh')
	                obj.setLocation(self.cube.coord[i][0]+self.cursorXYZ[0],self.cube.coord[i][1]+self.cursorXYZ[1],self.cube.coord[i][2]+self.cursorXYZ[2])
	        for i in range(len(self.cube.atomtypes)):
	                for j in range(i+1,len(self.cube.atomtypes)):
	                        vec1=Mathutils.Vector(self.cube.coord[i])
	                        vec2=Mathutils.Vector(self.cube.coord[j])
	                        vec=vec2-vec1
	                        distcovalent=covalentradii.get(self.cube.atomtypes[i],2.0)+covalentradii.get(self.cube.atomtypes[j],2.0)
	                        if (vec.length-distcovalent)<=0.10*distcovalent:
	                                me=Mesh.Primitives.Tube(32,stickradius,vec.length)
		        		for face in me.faces:
						face.smooth=True
	                                obj=self.scene.objects.new(me,'Cylinder')
	                                axis=Mathutils.CrossVecs(Vector([0,0,1]),vec)
	                                angle=Mathutils.AngleBetweenVecs(Vector([0,0,1]),vec)
	                                rotmat=Mathutils.RotationMatrix(angle,4,"R",axis)
	                                obj.setMatrix(obj.matrix*rotmat)
					loc = (vec1+vec2)*0.5 + Mathutils.Vector(self.cursorXYZ)
	                                obj.setLocation(loc)
		
	def isosurface(self,cubeobject,isovalue):
        	me = Mesh.New()
        	me2=Mesh.New()
        	faces=[]
        	vertices = []
		scalars = vtk.vtkFloatArray()
		for k in range(cubeobject.nz):
			for j in range(cubeobject.ny):
				for i in range(cubeobject.nx):
					scalars.InsertNextValue(cubeobject.data3d[k][j][i])
		grid=vtk.vtkStructuredPoints()
                info=vtk.vtkInformation()
		grid.SetOrigin(cubeobject.origin[0],cubeobject.origin[1],cubeobject.origin[2])
		grid.SetDimensions(cubeobject.nx,cubeobject.ny,cubeobject.nz)
		grid.SetSpacing(cubeobject.dx,cubeobject.dy,cubeobject.dz)
		grid.SetNumberOfScalarComponents(cubeobject.nx*cubeobject.ny*cubeobject.nz, info)
		grid.GetPointData().SetScalars(scalars)
        	Marching = vtk.vtkContourFilter()
        	Marching.SetInputData(grid)
        	Marching.SetValue(0, isovalue)
        	Marching.Update()
        	contoursurface=Marching.GetOutput()
        	for i in range(contoursurface.GetNumberOfPoints()):
       		        point = contoursurface.GetPoint(i)
        	        vertices.append([point[0]+self.cursorXYZ[0],point[1]+self.cursorXYZ[1],point[2]+self.cursorXYZ[2]])
        	me.verts.extend(vertices)
		me.materials=[materials["positivelobe"]]
        	for i in range(contoursurface.GetNumberOfCells()):
        	        cell = contoursurface.GetCell(i)
        	        n1 = cell.GetPointId(0)
        	        n2 = cell.GetPointId(1)
        	        n3 = cell.GetPointId(2)
			faces.append([me.verts[n1], me.verts[n2], me.verts[n3]])
        	me.faces.extend(faces)
		for face in me.faces:
			face.smooth=True
        	ob = self.scene.objects.new(me, 'Positive Lobe')
        	Marching.SetValue(0, -isovalue)
        	Marching.Update()
        	contoursurface=Marching.GetOutput()
        	vertices = []
        	for i in range(contoursurface.GetNumberOfPoints()):
        	        point = contoursurface.GetPoint(i)
        	        vertices.append([point[0]+self.cursorXYZ[0],point[1]+self.cursorXYZ[1],point[2]+self.cursorXYZ[2]])
        	me2.verts.extend(vertices)
		faces=[]
        	for i in range(contoursurface.GetNumberOfCells()):
        	        cell = contoursurface.GetCell(i)
        	        n1 = cell.GetPointId(0)
        	        n2 = cell.GetPointId(1)
       		        n3 = cell.GetPointId(2)
        		faces.append([me2.verts[n1], me2.verts[n2], me2.verts[n3]])
		me2.faces.extend(faces)
		me2.materials=[materials["negativelobe"]]
   		for face in me2.faces:
			face.smooth=True    
		ob = self.scene.objects.new(me2, 'Negative Lobe')
	
	def animatedensity(self,filenamelist,renderingpath,isovalue):
		context=self.scene.getRenderingContext()
		context.extensions=True
		context.renderPath=renderingpath
		#context.imageType=Render.JPEG
		context.sFrame=1
		context.eFrame=1
		context.sizeX=width
		context.sizeY=height
		
		cubeobject=cube(filenamelist[0])
		cubeobject.readCube()
                atoms=[]
                ipokeytypeloc=Object.IpoKeyTypes.LOC
                ipokeytyperot=Object.IpoKeyTypes.ROT
                for i in range(len(cubeobject.atomtypes)):
                        me = Mesh.Primitives.Icosphere(spheresubdivisions,atomicradii.get(cubeobject.atomtypes[i], 2.0))
                        me.materials=[materials.get(cubeobject.atomtypes[i], materials["default"])]
                        for face in me.faces:
                                face.smooth=True
                        obj=self.scene.objects.new(me,'Atom')
                        obj.setLocation(cubeobject.coord[i][0],cubeobject.coord[i][1],cubeobject.coord[i][2])
                        atoms.append(obj)

                bonds=[]
                for i in range(len(cubeobject.atomtypes)):
                        for j in range(i+1,len(cubeobject.atomtypes)):
                                vec1=Mathutils.Vector(cubeobject.coord[i])
                                vec2=Mathutils.Vector(cubeobject.coord[j])
                                vec=vec2-vec1
                                distcovalent=covalentradii.get(cubeobject.atomtypes[i],2.0)+covalentradii.get(cubeobject.atomtypes[j],2.0)
                                if (vec.length-distcovalent)<=0.10*distcovalent:
                                        me=Mesh.Primitives.Tube(32,stickradius,vec.length)
                                        for face in me.faces:
                                                face.smooth=True
                                        obj=self.scene.objects.new(me,'Cylinder')
                                        axis=Mathutils.CrossVecs(Vector([0,0,1]),vec)
                                        angle=Mathutils.AngleBetweenVecs(Vector([0,0,1]),vec)
                                        rotmat=Mathutils.RotationMatrix(angle,4,"R",axis)
                                        obj.setMatrix(obj.matrix*rotmat)
                                        obj.setLocation((vec1+vec2)*0.5)
                                        bonds.append([i,j,vec,obj,vec.length])
		self.isosurface(cubeobject,isovalue)
		context.render()
		filename="RENDER/image_0000.jpg"
		context.saveRenderedImage(filename)


                for j in range(1,len(filenamelist)):
			print "Rendering:", j
			filename="RENDER/image_%04d.jpg" % (j)
			for ob in self.scene.objects:
				if "Lobe" in ob.name:
					self.scene.unlink(ob)
			cubeobject=cube(filenamelist[j])
			cubeobject.readCube() 
			for i in range(len(atoms)): 
				atoms[i].setLocation(cubeobject.coord[i][0],cubeobject.coord[i][1],cubeobject.coord[i][2]) 
				atoms[i].insertIpoKey(ipokeytypeloc)
                        for i in range(len(bonds)):
                                vec1=Mathutils.Vector(cubeobject.coord[bonds[i][0]])
                                vec2=Mathutils.Vector(cubeobject.coord[bonds[i][1]])
                                vec=vec2-vec1
                                dist=vec.length
                                axis=Mathutils.CrossVecs(bonds[i][2],vec)
                                angle=Mathutils.AngleBetweenVecs(bonds[i][2],vec)
                                rotmat=Mathutils.RotationMatrix(angle,4,"R",axis)
                                bonds[i][3].setMatrix(bonds[i][3].matrix*rotmat)
                                bonds[i][3].setLocation((vec1+vec2)*0.5)
                                bonds[i][3].setSize(1.0,1.0,dist/bonds[i][4])
                                bonds[i][3].insertIpoKey(ipokeytypeloc)
                                bonds[i][3].insertIpoKey(ipokeytyperot)
                                bonds[i][2]=vec
                                bonds[i][4]=dist
			self.isosurface(cubeobject,isovalue)
			context.render()
			context.saveRenderedImage(filename)

		


EVENT_none=0
EVENT_Button_Import=1
EVENT_Button_Cancel=2
EVENT_Button_Browse=3
EVENT_Button_Plot=4
EVENT_Slider_Isovalue=5
EVENT_Button_Animate=6
inputfile=Create(pwd + "/Dyson_0001.cube")
renderingpath=Create(pwd +"/")
sliderMin=Create(0.0)
sliderMax=Create(1.0)
isovalue=Create(0.001)

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
		Text("Animate Electron Densities", "large")


		String("File path: ", EVENT_none,40, 300, 360, 20,inputfile.val,399, "cube file path")

		Button("Import (to cursor pos.)", EVENT_Button_Import, 110, 260, 140, 24, "")
		Button("Cancel", EVENT_Button_Cancel, 320, 260, 70, 24, "")
		Button("Browse", EVENT_Button_Browse, 400, 300, 70, 20, "")

		String("Rendering path: ", EVENT_none,40, 140, 360, 20,renderingpath.val,399, "Rendering directory")
		Button("Browse", EVENT_Button_Browse, 400, 140, 70, 20, "")
	
		glRasterPos2i(190, 220)
		Text("Plot Isosurface", "large")
	
		isovalue=Slider("Isovalue:",EVENT_Slider_Isovalue,110,180,210,24,isovalue.val,sliderMin.val,sliderMax.val,1,"isovalue")

		Button("Plot", EVENT_Button_Plot, 320, 180, 70, 24, "")
		Button("Animate", EVENT_Button_Animate, 200, 100, 70, 24, "")

	def event(self,event, value):
		if (event == ESCKEY or event == QKEY) and not value:
			Exit()

	def b_event(self,event):
		global EVENT_Button_Import
		global EVENT_Button_Cancel
		global EVENT_Button_Browse
		global EVENT_Button_Plot
		global sliderMin
		global sliderMax
		if event == 0: pass
		elif event == EVENT_Button_Import:
			self.cubeobject=cube(inputfile.val)
			self.cubeobject.readCube()
			self.cube2blenderobject=cube2blender(self.cubeobject)
			self.cube2blenderobject.blenderstructure()
			sliderMin.val=min(self.cubeobject.data1d)
			sliderMax.val=max(self.cubeobject.data1d)
			
		

		elif event == EVENT_Button_Animate:
			filenamelist1=os.listdir(os.path.dirname(inputfile.val))
			indexlist=[]
                        name0 = os.path.basename(inputfile.val)
			for name in filenamelist1:
                                print "File: %s" % name
                                if name[-5:] == ".cube" and name[:-9] == name0[:-9]:
                                        idx = int(name[-9:-5])
                                        indexlist.append(idx)
                                        print "added with index %d" % idx
			filenamelist=[]
			indexlist.sort()
			for i in indexlist:
        			path2=os.path.dirname(inputfile.val)+"/%s%04d.cube" % (name0[:-9],i)
				filenamelist.append(path2)
			self.cube2blenderobject.animatedensity(filenamelist,renderingpath.val,isovalue.val)

		elif event == EVENT_Button_Plot:
			self.cube2blenderobject.isosurface(self.cubeobject,isovalue.val)
		elif event == EVENT_Button_Cancel:
			Exit()
		elif event == EVENT_Button_Browse:
			Window.FileSelector(self.getFilename_callback, "Select cube file")
		Draw()
		
		
guiobject=GUI()
Register(guiobject.gui_draw, guiobject.event, guiobject.b_event)

#!/usr/bin/env python
"""
a graphical tool for analysing trajectories
"""
from pyface.qt import QtGui, QtCore
from PySide import QtUiTools

from DFTB.Analyse.mayavi.CubeViewerWidget import QCubeViewerWidget

import matplotlib
matplotlib.rc('xtick', labelsize=16)
matplotlib.rc('ytick', labelsize=16)

from matplotlib.figure import Figure
from matplotlib import patches
from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)

#import Avogadro

import numpy as np
from DFTB import XYZ
import os.path
import glob
import sys

##########
from DFTB.AtomicData import atomic_number, bohr_to_angs

def read_xyz_datafields_it(filename, units="Angstrom", fragment_id=atomic_number):
    """
    The columns after the x-,y- and z-coordinates in an XYZ-file can be used
    to store additional data, such as charges. This function reads both
    the geometry and the datafield and stores it as a numpy array

    Parameters:
    ===========
    filename: path to xyz-file
    units: specify the units of coordinates in xyz-file, "Angstrom" or "bohr"

    Returns:
    ========
    iterator to individual structures and datafields
    """
    assert units in ["bohr", "Angstrom", "hartree/bohr", ""]
    fh = open(filename)
    igeo = 1
    while 1:
        line = fh.readline()
        if not line: break
        try:
            nat = int(line.split()[0])
        except ValueError as e:
            print e
            raise Exception("Probably wrong number of atoms in xyz-file '%s'" % filename)
        title = fh.readline()
        atoms = []
        datafield = []
        for i in xrange(nat):
            line = fh.readline()
            words = line.split()
            atno = fragment_id(words[0])
            x,y,z = map(float,words[1:4])
            datafield.append( map(float,words[4:]) )
            if units == "Angstrom":
                x,y,z = map(lambda c: c/bohr_to_angs, [x,y,z])
            atoms.append((atno,(x,y,z)))
        igeo += 1
        datafield = np.array(datafield)
        yield atoms, datafield

def read_xyz_datafields(filename, units="Angstrom"):
    geometries = []
    datafields = []
    for atoms,datafield in read_xyz_datafields_it(filename, units=units):
        geometries.append(atoms)
        datafields.append(datafield)
    return geometries, datafields

#####

def loadUiWidget(uifilename, parent=None):
    loader = QtUiTools.QUiLoader()
    uifile = QtCore.QFile(uifilename)
    uifile.open(QtCore.QFile.ReadOnly)
    ui = loader.load(uifile, parent)
    uifile.close()
    return ui

class SaveMovieDialog(QtGui.QFileDialog):
    def __init__(self, parent = None):
        super(SaveMovieDialog, self).__init__(parent, "Save Movie")
        layout = self.layout()
        # select start and end frame and step
        timeFrame = QtGui.QFrame()
        layout.addWidget(timeFrame)
        timeLayout = QtGui.QHBoxLayout(timeFrame)
        startLabel = QtGui.QLabel("start frame:")
        self.startInput = QtGui.QLineEdit("0")
        self.startInput.setValidator(QtGui.QIntValidator())
        endLabel = QtGui.QLabel("end frame:")
        self.endInput = QtGui.QLineEdit("-1")
        self.endInput.setValidator(QtGui.QIntValidator())
        stepLabel = QtGui.QLabel("step:")
        self.stepInput = QtGui.QLineEdit("1")
        self.stepInput.setValidator(QtGui.QIntValidator())
        timeLayout.addWidget(startLabel)
        timeLayout.addWidget(self.startInput)
        timeLayout.addWidget(endLabel)
        timeLayout.addWidget(self.endInput)
        timeLayout.addWidget(stepLabel)
        timeLayout.addWidget(self.stepInput)

    @staticmethod
    def getDirnameAndFrames(parent = None):
        dialog = SaveMovieDialog(parent)
        result = dialog.exec_()
        start = int(dialog.startInput.text())
        end   = int(dialog.endInput.text())
        step  = int(dialog.stepInput.text())
        path = str(dialog.selectedFiles()[0])
        return (result == QtGui.QDialog.Accepted, path, start, end, step)
        
class Main(QtGui.QMainWindow):
    def __init__(self, ):
        super(Main, self).__init__()
        self.ui = loadUiWidget(os.path.join(os.path.dirname(os.path.abspath(__file__)), "window2.ui"))
        self.setCentralWidget(self.ui)
        self.canvases = []
        self.trajdirs = []

        self.dt = 0.0
        self.geometries = []
        self.particle_hole_charges = []
        
        self.ui.trajSelection.itemSelectionChanged.connect(self.changeTrajectory)
        
        # load Trajectories
        dyndir = "./"
        if len(sys.argv) > 1:
            dyndir = sys.argv[1]
        self.ui.dirBrowser.setText(dyndir)
        self.ui.dirBrowser.returnPressed.connect(self.fetchURL)
        self.ui.openButton.clicked.connect(self.openDirectory)

        # save movie
        self.ui.movieButton.clicked.connect(self.renderMovie)
        self.stop = False
        # plots
        import matplotlib.pyplot as plt
        cmap = plt.get_cmap('gnuplot')
        self.colors = ['black', 'red', 'green', 'blue', 'orange', 'm'] + [cmap(f) for f in np.linspace(0.4,1.0,10)]

        self.energiesAxis = self.addTab(self.ui.energiesLayout)
        self.coefficientsAxis = self.addTab(self.ui.coefficientsLayout)
        self.couplingsAxis = self.addTab(self.ui.couplingsLayout)
        self.currEnergyAxis = self.addTab(self.ui.currEnergyLayout)
        self.currStateAxis = self.addTab(self.ui.currStateLayout)
        self.populationsAxis = self.addTab(self.ui.populationsLayout)
        self.axes = [self.energiesAxis, self.coefficientsAxis, self.couplingsAxis, self.currEnergyAxis, self.currStateAxis]
        
        self.objects = []

        # Window with geometries
        self.geometryViewer = QCubeViewerWidget()
        self.geometryViewer.setAnimation(True)
        #self.geometryViewer.hideIsoControls()
        self.geometryViewer.selectShowOptions(options=["charges"])
        self.ui.viewerWidgetLayout.addWidget(self.geometryViewer)

        """
        self.glWidget = self.newGLWidget()
        self.mol = Avogadro.molecules.addMolecule()
        self.objects.append( self.mol )
        self.glWidget.molecule = self.mol
        self.ui.viewerWidgetLayout.addWidget(Avogadro.toPyQt(self.glWidget))
        """
        #
        self.ui.animationSlider.valueChanged.connect(self.showMolecules)

        self.loadTrajectories(dyndir)

        #
        self.pictureAorDiab = QtGui.QComboBox()
        self.pictureAorDiab.addItems(["adiabatic", "local diabatic"])
        self.pictureAorDiab.setToolTip("Choose whether energies and couplings should be shown in the adiabatic or the local diabatic picture.")
        self.pictureAorDiab.currentIndexChanged.connect(self.switchPicture)
        self.picture = self.pictureAorDiab.currentText()
        self.ui.gridLayout.addWidget(self.pictureAorDiab, 1,0)

    def switchPicture(self):
        self.picture = self.pictureAorDiab.currentText()
        if self.picture == "adiabatic":
            self.ui.tabWidget.setTabText(0, "Adiab. Energies")
            self.ui.tabWidget.setTabText(2, "Nonadiab. Couplings")
        else:
            self.ui.tabWidget.setTabText(0, "Local Diab. Energies")
            self.ui.tabWidget.setTabText(2, "Local Diab. Coupl.")        
        self.changeTrajectory()
        
    def openDirectory(self):
        dyndir = QtGui.QFileDialog.getExistingDirectory()
        self.ui.dirBrowser.setText(dyndir)
        self.fetchURL()
        
    def fetchURL(self):
        dyndir = str( self.ui.dirBrowser.text() )
        print "TOP DIRECTORY: %s" % dyndir
        self.loadTrajectories(dyndir)
        
    def loadTrajectories(self, dyndir):
        dirs = glob.glob(os.path.join(dyndir, "*TRAJ*/"))
        if len(dirs) == 0:
            # open a single trajectory
            if os.path.isfile(os.path.join(dyndir, "dynamics.xyz")):
                self.ui.trajSelection.addItem("%d" % len(self.trajdirs))
                self.trajdirs.append( dyndir )
        else:
            dirs.sort()
            # remove previous trajectories
            self.trajdirs = []
            for i,trajdir in enumerate(dirs):
                if os.path.isfile(os.path.join(trajdir, "dynamics.xyz")):
                    print "Load trajectory form %s" % trajdir
                    self.ui.trajSelection.addItem("%d" % i)
                    self.trajdirs.append( trajdir )
                    
            self.plotPopulations()
            
            for ax in self.axes:
                ax.clear()
                ax.set_xlim((-1.0,1.0))
                ax.set_ylim((-1.0,1.0))
                ax.text(-0.75, 0.2, "Click on a trajectory", fontsize=30)
                ax.text(-0.75,-0.2, "    on the left"      , fontsize=30)
            for canvas in self.canvases:
                canvas.draw()    

                
    def plotPopulations(self):
        tmin=0.0
        tmax=0.0
        Nst=0
        Nt=0
        data_list = []
        for i,trajdir in enumerate(self.trajdirs):
            state_file = os.path.join(trajdir, "state.dat")
            try:
                data = np.loadtxt(state_file)
            except IOError as e:
                print e
                continue
            Nst = max(data[:,1].max()+1, Nst)
            tmin = min(data[:,0].min(), tmin)
            tmax = max(data[:,0].max(), tmax)
            Nt = max(len(data[:,0]), Nt)
            data_list.append(data) 
        Nst = int(Nst)
        Nt = int(Nt)
        print "%d electronic states" % Nst
        print "%d time steps between %s and %s" % (Nt, tmin, tmax)
        pop = np.zeros((Nt,Nst+1))
        pop[:,0] = np.linspace(tmin,tmax,Nt) # time axis in fs
        Ntraj = [0 for t in range(0, Nt)]  # count trajectories available at each time step
        for i,data in enumerate(data_list):
            # only consider trajectories that finished nicely
            Nt_i = data.shape[0]
            if Nt_i != Nt:
                print "Trajectory %d has only %d time steps" % (i,Nt_i)
            for t in range(0, Nt_i):
                st = data[t,1]
                pop[t,st+1] += 1
                Ntraj[t] += 1.0
        print "%s trajectories" % Ntraj[0]
        # divide populations by the number of trajectories
        for t in range(0, Nt):
            pop[t,1:] /= float(Ntraj[t])

        self.populationsAxis.clear()
        self.populationsAxis.set_xlabel("Time / fs", fontsize=15)
        self.populationsAxis.set_ylabel("Populations", fontsize=15)
        for i in range(0, Nst):
            self.populationsAxis.plot(pop[:,0], pop[:,1+i], lw=2, color=self.colors[i], label="State %d" % i)
                
    def addTab(self, layout):
        fig = Figure()
        ax = fig.add_subplot(111)
        ax.set_xlim((-1.0,1.0))
        ax.set_ylim((-1.0,1.0))
        ax.text(-0.75, 0.2, "Click on a trajectory", fontsize=30)
        ax.text(-0.75,-0.2, "    on the left"      , fontsize=30)

        canvas = FigureCanvas(fig)
        self.canvases.append(canvas)
        layout.addWidget(canvas)
        canvas.draw()
        toolbar = NavigationToolbar(canvas,
                                    canvas, coordinates=True)
        return ax

    def newGLWidget(self):
        glWidget = Avogadro.GLWidget()
        glWidget.loadDefaultEngines()
        glWidget.quality = 4
        toolGroup = Avogadro.ToolGroup()
        tool = Avogadro.PluginManager.instance.tool('Navigate', None)
        if tool:
            toolGroup.append(tool)
            self.objects.append(tool)
        glWidget.toolGroup = toolGroup
        self.objects.append(glWidget)
        self.objects.append(toolGroup)
        return glWidget
        
    def changeTrajectory(self):
        selected = self.ui.trajSelection.selectedItems()
        trajs = map(int, [s.text() for s in selected])
        self.geometries = []
        self.nsteps = 0
        print "Selected: %s" % map(str, trajs)
        print "Plot trajectories %s" % trajs
        for ax in self.axes:
            ax.clear()
        for t in trajs:
            trajdir = self.trajdirs[t]
            print trajdir
            try:
                self.plotTrajectory(trajdir, t)
                if os.path.isfile(os.path.join(trajdir, "particle_hole_charges.xyz")):
                    geoms, ph_charges = read_xyz_datafields(os.path.join(trajdir, "particle_hole_charges.xyz"))
                else:
                    geoms = XYZ.read_xyz(os.path.join(trajdir, "dynamics.xyz"))
                    ph_charges = None                
                self.geometries.append( geoms )
                self.particle_hole_charges.append( ph_charges )
                self.nsteps = max(len(self.geometries[-1]), self.nsteps)
            except (ValueError, IOError) as e:
                print "Unable to load %s" % trajdir
                print "ERROR: %s" % e
        for canvas in self.canvases:
            canvas.draw()    

        self.ui.animationSlider.setMinimum(0)
        self.ui.animationSlider.setMaximum(self.nsteps)
        self.ui.animationTime.setText("Time: %7.2f fs" % 0.0)
        self.geometryViewer.clear()
        self.showMolecules()
        
    def showMolecules(self):
        tstep = self.ui.animationSlider.value()
        self.ui.animationTime.setText("Time: %7.2f fs" % (tstep * self.dt))
        molecules = []
        charge_lists = []
        for geoms,ph_charges in zip(self.geometries, self.particle_hole_charges):
            geom_time = geoms[min(tstep, len(geoms)-1)]
            if type(ph_charges) == type(None):
                ph_charges_times = np.zeros((len(geom_time),2))
            else:
                ph_charges_times = ph_charges[min(tstep, len(ph_charges)-1)]
                charge_lists.append( -ph_charges_times)
            molecules.append(geom_time)
            # change sign of charges, internally charges are represented in units of e-!
        if len(charge_lists) == 0:
            self.geometryViewer.setGeometries(molecules)
        else:
            self.geometryViewer.setGeometriesAndCharges(molecules, charge_lists)

    def abortCurrentAction(self):
        self.stop = True
        #raise AbortException() 
    def renderMovie(self):
        if len(self.geometries) == 0:
            # no trajectories
            return
        # ask for filename and start and end frames
        (ok, path, start, end, step) = SaveMovieDialog.getDirnameAndFrames()
        if ok == False:
            return
        steps = range(0, self.nsteps)
        self.ui.movieButton.setText("Stop")
        self.ui.movieButton.clicked.disconnect()
        self.ui.movieButton.clicked.connect(self.abortCurrentAction)

        dirname = os.path.dirname(path)
        basename = os.path.basename(path)
        if basename == '':
            basename = "movie.png"
        name, suffix = os.path.splitext(basename)
        print "Frames will be stored in %s" % os.path.join(dirname, "%s_######.png" % name)
        for i in steps[start:end:step]:
            print "render frame %d" % i
            self.ui.animationSlider.setValue(i)
            image_file = os.path.join(dirname, "%s_%06d.png" % (name, i))
            self.geometryViewer.savefig(image_file)
            if self.stop == True:
                print "stop rendering!"
                self.stop = False
                break
        
        self.ui.movieButton.setText("Movie...")
        self.ui.movieButton.clicked.disconnect()
        self.ui.movieButton.clicked.connect(self.renderMovie)
        
    def plotTrajectory(self, trajdir, trajID):
        # load adiabatic energies
        energy_files = glob.glob(os.path.join(trajdir, "energy_*.dat"))
        states = []
        energies_list = []
        for f in energy_files:
            states.append( int(f[-15:-4].split("_")[1]) )
            data = np.loadtxt(f)
            times = data[:,0]
            energies_list.append( data[:,1])
        Nst = len(states)
        Nt = len(times)
        energies = np.zeros((Nt, Nst))
        indx = np.argsort(states)
        for i in range(0, Nst):
            energies[:,i] = energies_list[indx[i]]
        # plot levels
        if self.picture == "adiabatic":
            print "adiabatic levels ...",
            self.energiesAxis.set_xlabel("Time / fs", fontsize=17)
            self.energiesAxis.set_ylabel("Adiab. Energy / eV", fontsize=17)
            for i in range(0, Nst):
                self.energiesAxis.plot(times, energies[:,i]*27.211, lw=2, color=self.colors[i], label="State %d" % i)
            print "...done"
        else:
            print "local diabatic levels ..."
            # local diabatic picture
            ld_energy_file = os.path.join(trajdir, "local_diabatic_energies.dat")
            self.energiesAxis.set_xlabel("Time / fs", fontsize=17)
            self.energiesAxis.set_ylabel("Local Diab. Energy / eV", fontsize=17)
            try:
                data = np.loadtxt(ld_energy_file)
                for i in range(0, Nst):
                    self.energiesAxis.plot(data[:,0], data[:,1+i]*27.211, lw=2, color=self.colors[i], label="State %d" % i)
                print "...done"
            except IOError as e:
                print "No diabatic energies"
                print e
        # plot current state
        print "current states ...",
        curr_states = np.loadtxt( os.path.join(trajdir, "state.dat") )
        curr_energy = []
        state_changes = [0]
        st_last = curr_states[0,1]
        for i in range(0, Nt):
            st = curr_states[i,1]
            if st_last != st:
                state_changes.append( i )
                st_last = st
            curr_energy.append( energies[i,st] )
        state_changes.append(i) 
        curr_energy = np.array(curr_energy)
        # show current energy
        # ... as a dashed brown line
        self.energiesAxis.plot(times, curr_energy*27.211, ls="-.", lw=3, color="brown", label="Curr. State")
        # ... or as circles
        #self.energiesAxis.plot(times[::15], (curr_energy*27.211)[::15], 'o', color='brown', markersize=13, mfc='none', label="Curr. State")
        print "...done"
        # time-step in fs
        self.dt = times[1]-times[0]
        # couplings
        if self.picture == "adiabatic":
            # non-adiabatic couplings
            print "non-adiabatic couplings...",
            self.couplingsAxis.set_xlabel("Time / fs", fontsize=15)
            self.couplingsAxis.set_ylabel("Scalar Non-adiab. Coupling", fontsize=15)
            for i in range(0, Nst):
                for j in range(i+1,Nst):
                    try:
                        data = np.loadtxt(os.path.join(trajdir, "nonadiabatic"+str(i)+str(j)+".dat"))
                        self.couplingsAxis.plot(data[:,0], data[:,1], lw=2, label="Non-adiab. coupling %d-%d" % (i,j))
                    except IOError as e:
                        print e
            print "...done"
        else:
            # local diabatic couplings
            print "local diabatic coupling...",
            ld_coupling_file = os.path.join(trajdir, "local_diabatic_couplings.dat")
            self.couplingsAxis.set_xlabel("Time / fs", fontsize=15)
            self.couplingsAxis.set_ylabel("Scalar Local Diab. Coupling", fontsize=15)
            try:
                data = np.loadtxt(ld_coupling_file)
                c = 1
                for i in range(0, Nst):
                    for j in range(i+1,Nst):
                        self.couplingsAxis.plot(data[:,0], data[:,c], lw=2, label="Local diab. coupling %d-%d" % (i,j))
                        c += 1
                print "...done"
            except IOError as e:
                print "No diabatic couplings"
                print e
            
        # coefficients
        print "coefficients ...",
        self.coefficientsAxis.set_xlabel("Time / fs", fontsize=15)
        self.coefficientsAxis.set_ylabel("Electronic Coefficients $\\vert C_i \\vert^2$", fontsize=15)
        for i in range(0, Nst):
            try:
                data = np.loadtxt(os.path.join(trajdir, "coeff_"+str(i)+".dat"))
                self.coefficientsAxis.plot(data[:,0], data[:,1], lw=2, color=self.colors[i], label="$\\vert C_%d \\vert^2$" % i)
            except IOError as e:
                print e
        print "...done"
        # current energy
        print "current energy ...",
        self.currEnergyAxis.set_xlabel("Time / fs", fontsize=15)
        self.currEnergyAxis.set_ylabel("Current Energy / eV", fontsize=15)
        try:
            data = np.loadtxt(os.path.join(trajdir, "curr_energy.dat"))
            self.currEnergyAxis.plot(data[:,0], data[:,1]*27.211, lw=2, color="black", label="kinetic")
            self.currEnergyAxis.plot(data[:,0], data[:,2]*27.211, lw=2, color="red"  , label="potential")
            if data.shape[1] > 3:
                self.currEnergyAxis.plot(data[:,0], data[:,3]*27.211, lw=2, color="green", label="total")
            legend = self.currEnergyAxis.get_legend()
            if legend == None:
                handles, labels = self.currEnergyAxis.get_legend_handles_labels()
                self.currEnergyAxis.legend(handles, labels, loc=1)
        except IOError as e:
            print e
        print "...done"
        # current states as blocks
        print "current state blocks ...",
        self.currStateAxis.set_yticks(range(0, len(self.trajdirs)))
        self.currStateAxis.set_xlabel("Time / fs", fontsize=15)
        self.currStateAxis.set_ylabel("Trajectory", fontsize=15)
        state = curr_states[0,1]
        for i in range(1, len(state_changes)):
            xi = times[state_changes[i-1]]
            yi = trajID-0.45
            width = times[state_changes[i]] - times[state_changes[i-1]]
            st = int(curr_states[state_changes[i-1],1])
            self.currStateAxis.add_patch(
                patches.Rectangle( (xi,yi), width, 0.90 , color=self.colors[st])
                )
            self.currStateAxis.text(xi+0.5*width, yi+0.25, "%d" % st)
        xmin,xmax = self.currStateAxis.get_xlim()
        self.currStateAxis.set_xlim((min(times[0], xmin), max(times[-1], xmax)))
        ymin,ymax = self.currStateAxis.get_ylim()
        self.currStateAxis.set_ylim((min(ymin, trajID-1), max(ymax,trajID+1)))
        print "...done"
            
if __name__ == "__main__":
    #app = QtGui.QApplication(sys.argv)
    app = QtGui.QApplication.instance()
    main = Main()
    main.show()
    sys.exit(app.exec_())

#!/usr/bin/env python
"""
a graphical tool for analysing trajectories
"""
from PyQt4.uic import loadUiType
from PyQt4 import QtGui, QtCore

from matplotlib.figure import Figure
from matplotlib import patches
from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)

import Avogadro

import numpy as np
from DFTB import XYZ
import os.path
import glob
import sys

Ui_MainWindow, QMainWindow = loadUiType(os.path.join(os.path.dirname(os.path.abspath(__file__)), "window.ui"))

class Main(QMainWindow, Ui_MainWindow):
    def __init__(self, ):
        super(Main, self).__init__()
        self.setupUi(self)
        self.canvases = []
        self.trajdirs = []

        self.trajSelection.itemSelectionChanged.connect(self.changeTrajectory)
        
        # load Trajectories
        dyndir = "./"
        if len(sys.argv) > 1:
            dyndir = sys.argv[1]
        self.dirBrowser.setText(dyndir)
        self.dirBrowser.returnPressed.connect(self.fetchURL)
        self.openButton.clicked.connect(self.openDirectory)

        # plots
        import matplotlib.pyplot as plt
        cmap = plt.get_cmap('gnuplot')
        self.colors = ['black', 'red', 'green', 'blue', 'k', 'm'] + [cmap(f) for f in np.linspace(0.4,1.0,10)]

        self.energiesAxis = self.addTab(self.energiesLayout)
        self.coefficientsAxis = self.addTab(self.coefficientsLayout)
        self.couplingsAxis = self.addTab(self.couplingsLayout)
        self.currEnergyAxis = self.addTab(self.currEnergyLayout)
        self.currStateAxis = self.addTab(self.currStateLayout)
        self.populationsAxis = self.addTab(self.populationsLayout)
        self.axes = [self.energiesAxis, self.coefficientsAxis, self.couplingsAxis, self.currEnergyAxis, self.currStateAxis]
        
        self.objects = []

        # Window with geometries
        self.glWidget = self.newGLWidget()
        self.mol = Avogadro.molecules.addMolecule()
        self.objects.append( self.mol )
        self.glWidget.molecule = self.mol
        self.viewerWidgetLayout.addWidget(Avogadro.toPyQt(self.glWidget))

        #
        self.animationSlider.valueChanged.connect(self.showMolecules)

        self.loadTrajectories(dyndir)

        #
        self.pictureAorDiab = QtGui.QComboBox()
        self.pictureAorDiab.addItems(["adiabatic", "local diabatic"])
        self.pictureAorDiab.setToolTip("Choose whether energies and couplings should be shown in the adiabatic or the local diabatic picture.")
        self.pictureAorDiab.currentIndexChanged.connect(self.switchPicture)
        self.picture = self.pictureAorDiab.currentText()
        self.gridLayout.addWidget(self.pictureAorDiab, 1,0)
        
    def switchPicture(self):
        self.picture = self.pictureAorDiab.currentText()
        if self.picture == "adiabatic":
            self.tabWidget.setTabText(0, "Adiab. Energies")
            self.tabWidget.setTabText(2, "Nonadiab. Couplings")
        else:
            self.tabWidget.setTabText(0, "Local Diab. Energies")
            self.tabWidget.setTabText(2, "Local Diab. Coupl.")        
        self.changeTrajectory()
        
    def openDirectory(self):
        dyndir = QtGui.QFileDialog.getExistingDirectory()
        self.dirBrowser.setText(dyndir)
        self.fetchURL()
        
    def fetchURL(self):
        dyndir = str( self.dirBrowser.text() )
        print "TOP DIRECTORY: %s" % dyndir
        self.loadTrajectories(dyndir)
        
    def loadTrajectories(self, dyndir):
        dirs = glob.glob(os.path.join(dyndir, "*TRAJ*/"))
        if len(dirs) == 0:
            # open a single trajectory
            if os.path.isfile(os.path.join(dyndir, "dynamics.xyz")):
                self.trajSelection.addItem("%d" % len(self.trajdirs))
                self.trajdirs.append( dyndir )
        else:
            dirs.sort()
            # remove previous trajectories
            self.trajdirs = []
            for i,trajdir in enumerate(dirs):
                if os.path.isfile(os.path.join(trajdir, "dynamics.xyz")):
                    print "Load trajectory form %s" % trajdir
                    self.trajSelection.addItem("%d" % i)
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
        selected = self.trajSelection.selectedItems()
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
                self.geometries.append( XYZ.read_xyz(os.path.join(trajdir, "dynamics.xyz")) )
                self.nsteps = max(len(self.geometries[-1]), self.nsteps)
            except (ValueError, IOError) as e:
                print "Unable to load %s" % trajdir
                print "ERROR: %s" % e
        for canvas in self.canvases:
            canvas.draw()    

        self.animationSlider.setMinimum(0)
        self.animationSlider.setMaximum(self.nsteps)
        self.animationTime.setText("Time: %7.2f fs" % 0.0)
        self.showMolecules()
        
    def showMolecules(self):
        tstep = self.animationSlider.value()
        self.animationTime.setText("Time: %7.2f fs" % (tstep * self.dt))
        self.mol.clear()
        for geoms in self.geometries:
            geom_time = geoms[min(tstep, len(geoms)-1)]
            # atoms
            atoms = []
            for (Z,pos) in geom_time:
                a = self.mol.addAtom()
                a.pos = np.array(pos, dtype=float)
                a.atomicNumber = Z
                atoms.append(a)
            # bond
            C = XYZ.connectivity_matrix(geom_time)
            Nat = len(geom_time)
            for i in range(0, Nat):
                for j in range(i+1,Nat):
                    if C[i,j] == 1:
                        bond = self.mol.addBond()
                        bond.setBegin(atoms[i])
                        bond.setEnd(atoms[j])
                        bond.order = 1
                        
        self.glWidget.mol = self.mol
        self.glWidget.updateGeometry()
        Avogadro.toPyQt(self.glWidget).update()
        
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
            self.energiesAxis.set_xlabel("Time / fs", fontsize=15)
            self.energiesAxis.set_ylabel("Adiab. Energy / eV", fontsize=15)
            for i in range(0, Nst):
                self.energiesAxis.plot(times, energies[:,i]*27.211, lw=2, color=self.colors[i], label="State %d" % i)
            print "...done"
        else:
            print "local diabatic levels ..."
            # local diabatic picture
            ld_energy_file = os.path.join(trajdir, "local_diabatic_energies.dat")
            self.energiesAxis.set_xlabel("Time / fs", fontsize=15)
            self.energiesAxis.set_ylabel("Local Diab. Energy / eV", fontsize=15)
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
        self.energiesAxis.plot(times, curr_energy*27.211, ls="-.", lw=3, color="brown", label="Curr. State")
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
    app = QtGui.QApplication(sys.argv)
    main = Main()
    main.show()
    sys.exit(app.exec_())

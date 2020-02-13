"""
graphical user interface for displaying the results of a TD-DFTB calculation
"""
from pyface.qt import QtGui, QtCore
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)

import string
import numpy as np

from DFTB import AtomicData
from DFTB.Analyse.absorption_spectrum import broadened_spectrum, convert_energy
from DFTB.Analyse.mayavi.CubeViewerWidget import QCubeViewerWidget, CubeData
from DFTB.BasisSets import AtomicBasisSet
from DFTB.Analyse import Cube
from DFTB.GammaApproximation import AuxiliaryBasisSet

class QSpectrumWidget(QtGui.QWidget):
    def __init__(self, tddftb, parent=None):
        self.tddftb = tddftb
        
        QtGui.QWidget.__init__(self, parent)
        layout = QtGui.QVBoxLayout(self)
        #
        upperFrame = QtGui.QFrame()
        layout.addWidget(upperFrame)
        upperLayout = QtGui.QHBoxLayout(upperFrame)
        # table with excitation energies and states
        tableFrame = QtGui.QFrame()
        upperLayout.addWidget(tableFrame)
        tableLayout = QtGui.QVBoxLayout(tableFrame)
        tableLayout.addWidget(QtGui.QLabel("Excited States:"))
        self.tableStates = QtGui.QTableWidget()
        self.tableStates.setToolTip("Select a row to highlight the respective state in the absorption spectrum")
        self.tableStates.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
        self.tableStates.itemSelectionChanged.connect(self.plotSpectrum)
        tableLayout.addWidget(self.tableStates)
        cubeFrame = QtGui.QFrame()
        tableLayout.addWidget(cubeFrame)
        cubeLayout = QtGui.QHBoxLayout(cubeFrame)
        self.ppbGrid = QtGui.QComboBox()
        self.ppbGrid.addItems(["crude: 1.0", "coarse: 2.0", "medium: 3.0", "fine: 4.0"])
        self.ppbGrid.setCurrentIndex(1)
        self.ppbGrid.setToolTip("The number of points per bohr determines the resolution of the grid on which the cubes are calculated")
        self.ppbGrid.currentIndexChanged.connect(self.recalculateCubes)
        cubeLayout.addWidget(self.ppbGrid)
                             
        calcCubesButton = QtGui.QPushButton("Calc. cubes")
        calcCubesButton.setToolTip("Calculate cubes with the transition densities for the selected states")
        calcCubesButton.clicked.connect(self.calculateCubes)
        cubeLayout.addWidget(calcCubesButton)
        # transition densities
        tdenseFrame = QtGui.QFrame()
        upperLayout.addWidget(tdenseFrame)
        tdenseLayout = QtGui.QVBoxLayout(tdenseFrame)
        self.tdenseLabel = QtGui.QLabel("Excitated State:")
        tdenseLayout.addWidget(self.tdenseLabel)
        tdenseTab = QtGui.QTabWidget()
        tdenseLayout.addWidget(tdenseTab)
        # 3D cubes - transition density
        self.bs = AtomicBasisSet(self.tddftb.dftb2.atomlist)
        self.transitionDensityViewer = QCubeViewerWidget()
        tdenseTab.addTab(self.transitionDensityViewer, "3D Transition Density")
        # 3D cubes - density difference
        self.bs_aux = AuxiliaryBasisSet(self.tddftb.dftb2.atomlist, self.tddftb.dftb2.hubbard_U)
        self.differenceDensityViewer = QCubeViewerWidget()
        tdenseTab.addTab(self.differenceDensityViewer, "3D Difference Density")
        # 2D occ-virt plane
        exvec2dFrame = QtGui.QFrame()
        exvec2dLayout = QtGui.QVBoxLayout(exvec2dFrame)
        self.figExvec = Figure()
        self.canvasExvec = FigureCanvas(self.figExvec)
        exvec2dLayout.addWidget(self.canvasExvec)
        tdenseTab.addTab(exvec2dFrame, "2D Excitation Vector")
        NavigationToolbar(self.canvasExvec, exvec2dFrame, coordinates=True)
        # atomic charges in the excited state
        self.chargesViewer = QCubeViewerWidget()
        tdenseTab.addTab(self.chargesViewer, "Partial Charges")
        atomlist = self.tddftb.dftb2.getGeometry()
#        partial_charges = -self.tddftb.dftb2.getPartialCharges()
#        self.chargesViewer.setGeometriesAndCharges([atomlist], [partial_charges])
        self.chargesViewer.selectShowOptions(options=["charges"])
        # transition charges between ground and excited state
        self.trans_chargesViewer = QCubeViewerWidget()
        tdenseTab.addTab(self.trans_chargesViewer, "Transition Charges")
        atomlist = self.tddftb.dftb2.getGeometry()
#        # for S0->S0 transition, the transition charges are equal to the partial charges
#        transition_charges = -self.tddftb.dftb2.getPartialCharges()
#        self.trans_chargesViewer.setGeometriesAndCharges([atomlist], [transition_charges])
        self.trans_chargesViewer.selectShowOptions(options=["charges", "transition charges"])
        # figure with spectrum
        spectrumFrame = QtGui.QFrame()
        layout.addWidget(spectrumFrame)
        spectrumLayout = QtGui.QVBoxLayout(spectrumFrame)
        self.figSpectrum = Figure()
        self.canvasSpectrum = FigureCanvas(self.figSpectrum)
        spectrumLayout.addWidget(self.canvasSpectrum)
        NavigationToolbar(self.canvasSpectrum, spectrumFrame, coordinates=True)
        # controls
        controlFrame = QtGui.QFrame()
        spectrumLayout.addWidget(controlFrame)
        controlLayout = QtGui.QHBoxLayout(controlFrame)
        self.energyUnits = QtGui.QComboBox()
        
        self.energyUnits.addItems(["Hartree", "eV", "nm", "cm-1"])
        self.energyUnits.currentIndexChanged.connect(self.plotSpectrum)
        controlLayout.addWidget(QtGui.QLabel("units:"))
        controlLayout.addWidget(self.energyUnits)
        controlLayout.addWidget(QtGui.QLabel("broadening:"))
        self.broadening = QtGui.QLineEdit()
        self.broadening.setText("0.0")
        self.broadening.editingFinished.connect(self.plotSpectrum)
        self.broadening.setToolTip("The stick spectrum is convolved with a Gaussian to simulated a temperature broadened spectrum")
        controlLayout.addWidget(self.broadening)
        
        # load spectrum
        nr_dominant_ex = 2
        tab = self.tableStates
        tab.setRowCount(len(tddftb.Omega))
        headers = ["N", "Spin", "Sym", "exc. en. / hartree", "exc. en. / eV", "exc. en. / nm", "osc. strength", "Lambda diagn.", "Cube"]
        tab.setColumnCount(len(headers))
        tab.setHorizontalHeaderLabels(headers)
        for I in range(0, len(tddftb.Omega)):
            row = [string.rjust(str(I+1),5), \
                               tddftb.multiplicity, \
                               string.rjust(tddftb.Irreps[I], 3), \
                               string.rjust("%.7f" % tddftb.Omega[I],20), \
                               string.rjust("%.7f" % (tddftb.Omega[I]*AtomicData.hartree_to_eV), 17), \
                               string.rjust("%.7f" % (AtomicData.hartree_to_nm / tddftb.Omega[I]), 17), \
                               string.rjust("%.7f" % tddftb.oscillator_strength[I], 12), \
                               string.rjust("%.4f" % tddftb.Lambda2[I], 7)]
            for j,r in enumerate(row):
                tab.setItem(I,j, QtGui.QTableWidgetItem("%s" % r))
        tab.resizeColumnsToContents()

        self.plotSpectrum()
        # cubes with transition densities and difference densities for each state if calculated
        self.tdense_cubes = {}
        self.difdense_cubes = {}
        self.partial_charges = {}  # partial charges on excited states
    def plotSpectrum(self):
        #print "plot spectrum"
        self.figSpectrum.clf()
        ax = self.figSpectrum.add_subplot(111)
        
        units = self.energyUnits.itemText(self.energyUnits.currentIndex())
        broadening = abs(float(self.broadening.text()))
        ax.set_title("Absorption Spectrum")
        ax.set_xlabel("Excitation Energy / %s" % units, fontsize=15)
        ax.set_ylabel("Oscillator strength", fontsize=15)

        ens = self.tddftb.Omega
        osz = self.tddftb.oscillator_strength
        # stick spectrum
        ens_out = convert_energy(ens, "Hartree", units)
        ax.vlines(ens_out, np.zeros(ens_out.shape), osz, lw=2, color="black", picker=5)

        # selected states are highlighted
        selected = self.tableStates.selectionModel().selectedRows()
        for s in selected:
            I = s.row()
            ax.vlines(ens_out[I], [0.0], osz[I], lw=3, color="green")
            ax.text(ens_out[I], osz[I], "%s (%s$_{%d}$)" % (self.tddftb.Irreps[I], self.tddftb.multiplicity, I+1))

        def onpick(event):
            self.tableStates.setCurrentCell(event.ind[0], 0)
        self.figSpectrum.canvas.mpl_connect('pick_event', onpick)
            
        if broadening > 0.0:
            ens_broadened, spec = broadened_spectrum(ens, osz, broadening)
            ens_broadened_out = convert_energy(ens_broadened, "Hartree", units)
            scale = spec.max()/osz.max()/1.1
            ax.plot(ens_broadened_out, spec/scale, lw=1, color="blue")
                    
        self.canvasSpectrum.draw()

        if len(selected) > 0:
            I = selected[0].row()
            self.tdenseLabel.setText("Excited State: %s (%s%d)" % (self.tddftb.Irreps[I],self.tddftb.multiplicity, I+1))
            self.showPartialCharges(I)
            if self.tdense_cubes.has_key(I):
                tdense_cube = self.tdense_cubes[I]
                self.transitionDensityViewer.setCubes([tdense_cube])
                difdense_cube = self.difdense_cubes[I]
                self.differenceDensityViewer.setCubes([difdense_cube])
            else:
                print "No cube for transition and difference density of state %d" % (I+1)
            # transition densities in occ-virt plane
            self.plotExcitationCoefficients2D(I)
        self.canvasExvec.draw()
    def plotExcitationCoefficients2D(self, I):
        self.figExvec.clf()
        ax = self.figExvec.add_subplot(111)
        ax.set_title("$C^{%d}_{(o \\to v)}$" % (I+1))
        ax.set_xlabel("occupied orbital")
        ax.set_ylabel("virtual orbital")
        active_occupied_orbs, active_virtual_orbs = self.tddftb.getActiveOrbitals()
        xlabel_positions = np.arange(0, len(active_occupied_orbs))
        xlabels = [str(o+1) for o in active_occupied_orbs]
        xlabels[-1] = "HOMO %s" % xlabels[-1]
        ylabel_positions = np.arange(0, len(active_virtual_orbs))
        ylabels = [str(v+1) for v in active_virtual_orbs]
        ylabels[0] = "LUMO %s" % ylabels[0]
        ax.set_xticks(xlabel_positions)
        ax.set_xticklabels(xlabels, rotation=45, fontsize=10)
        ax.set_yticks(ylabel_positions)
        ax.set_yticklabels(ylabels, fontsize=10)
        self.Cexvec = self.tddftb._getExcitationCoefficients()[I,:,:]
        image = ax.imshow(self.Cexvec.transpose(), interpolation='none', origin='lower', picker=True, vmin=-1.0, vmax=1.0)
        self.figExvec.colorbar(image, ax=ax)

        def onpick(event):
            x,y = event.mouseevent.xdata, event.mouseevent.ydata
            o,v = int(np.floor(x+0.5)), int(np.floor(y+0.5))
            if self.Cexvec[o,v] < -0.5:
                color = "white"
            else:
                color = "black"
            print "%s -> %s    %s" % (active_occupied_orbs[o]+1,active_virtual_orbs[v]+1, self.Cexvec[o,v])
            ax.text(o-0.4,v-0.25,"$%d \\to %d$\n%+3.3f" % (active_occupied_orbs[o]+1,active_virtual_orbs[v]+1, self.Cexvec[o,v]), color=color, fontsize=12)
            self.canvasExvec.draw()
            
        self.figExvec.canvas.mpl_connect('pick_event', onpick)
    def showPartialCharges(self, I):
        atomlist = self.tddftb.dftb2.getGeometry()
        # particle-hole charges
        print "compute particle-hole charges for excited state %d" % (I+1)
        particle_charges, hole_charges = self.tddftb.ParticleHoleCharges(I)
        # partial charges on ground state
        dq0 = self.tddftb.dftb2.getPartialCharges()
        # partial charges on excited state
        dqI = particle_charges + hole_charges + dq0
        self.partial_charges[I] = dqI
        # internally charges are measured in units of e,
        charges = -dqI
        self.chargesViewer.setGeometriesAndCharges([atomlist], [charges])
        print "compute transition charges between ground and excited state %d" % (I+1)
        transition_charges = self.tddftb.TransitionChargesState(I)
        self.trans_chargesViewer.setGeometriesAndCharges([atomlist], [transition_charges])
    def calculateCubes(self):
        selected = self.tableStates.selectionModel().selectedRows()
        for s in selected:
            I = s.row()
            if not (self.tdense_cubes.has_key(I)):
                self.showPartialCharges(I)
                # need to calculated cube
                print "compute transition and difference density for state %d ..." % (I+1),
                tdense_cube, difdense_cube = self._computeTransitionAndDifferenceDensity(I)
                print "...done"
                self.tdense_cubes[I] = tdense_cube
                self.difdense_cubes[I] = difdense_cube
                
                self.tableStates.setItem(I, self.tableStates.columnCount()-1, QtGui.QTableWidgetItem("Y"))
                self.tableStates.item(I, self.tableStates.columnCount()-1).setBackground(QtGui.QColor(255,0,0))
        self.plotSpectrum()
    def recalculateCubes(self):
        # delete cubes
        for I in self.tdense_cubes.keys():
                self.tableStates.setItem(I, self.tableStates.columnCount()-1, QtGui.QTableWidgetItem(""))
                self.tableStates.item(I, self.tableStates.columnCount()-1).setBackground(QtGui.QColor(255,255,255))
        self.tdense_cubes = {}
        self.difdense_cubes = {}
        self.calculateCubes()
    def _computeTransitionAndDifferenceDensity(self, I):
        Ptrans = self.tddftb.TransitionDensityMatrix(I)
        P0,PI = self.tddftb.ExcitedDensityMatrix(I)
        Pdif = PI-P0
        
        cube = CubeData()
        atomlist = self.tddftb.dftb2.getGeometry()
        
        (xmin,xmax),(ymin,ymax),(zmin,zmax) = Cube.get_bbox(atomlist, dbuff=5.0)
        dx,dy,dz = xmax-xmin,ymax-ymin,zmax-zmin

        ppb_text = self.ppbGrid.itemText(self.ppbGrid.currentIndex())
        ppb = float(ppb_text.split()[1]) # points per bohr
        #print "points per bohr: %s" % ppb

        nx,ny,nz = int(dx*ppb),int(dy*ppb),int(dz*ppb)
        x,y,z = np.mgrid[xmin:xmax:nx*1j, ymin:ymax:ny*1j, zmin:zmax:nz*1j]
        grid = (x,y,z)
        # transition density
        tdenseGrid = Cube.electron_density(grid, self.bs.bfs, Ptrans)

        tdense_cube = CubeData()
        tdense_cube.data = tdenseGrid
        tdense_cube.grid = grid
        tdense_cube.atomlist = atomlist
        ## approximate difference density based on particle and hole charges
        #dqI = self.partial_charges[I]
        #difdenseGrid = self.bs_aux.partial_density(dqI, 0.0*self.tddftb.dftb2.ddip, x,y,z)
        # exact difference density
        difdenseGrid = Cube.electron_density(grid, self.bs.bfs, Pdif)
        
        difdense_cube = CubeData()
        difdense_cube.data = difdenseGrid
        difdense_cube.grid = grid
        difdense_cube.atomlist = atomlist
        
        return tdense_cube, difdense_cube

######################### MOLECULAR ORBITALS AND KOHN-SHAME ENERGIES #################
class QMolecularOrbitalWidget(QtGui.QWidget):
    def __init__(self, tddftb, parent=None):
        self.tddftb = tddftb
        
        QtGui.QWidget.__init__(self, parent)
        layout = QtGui.QVBoxLayout(self)
        #
        upperFrame = QtGui.QFrame()
        layout.addWidget(upperFrame)
        upperLayout = QtGui.QHBoxLayout(upperFrame)
        # table with MOs and Kohn-Sham energies
        tableFrame = QtGui.QFrame()
        upperLayout.addWidget(tableFrame)
        tableLayout = QtGui.QVBoxLayout(tableFrame)
        tableLayout.addWidget(QtGui.QLabel("Molecular Orbitals:"))
        self.tableMOs = QtGui.QTableWidget()
        self.tableMOs.setToolTip("Select a row to highlight the respective orbital in the density of states")
        self.tableMOs.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
        self.tableMOs.itemSelectionChanged.connect(self.plotDOS)
        tableLayout.addWidget(self.tableMOs)
        cubeFrame = QtGui.QFrame()
        tableLayout.addWidget(cubeFrame)
        cubeLayout = QtGui.QHBoxLayout(cubeFrame)
        self.ppbGrid = QtGui.QComboBox()
        self.ppbGrid.addItems(["crude: 1.0", "coarse: 2.0", "medium: 3.0", "fine: 4.0"])
        self.ppbGrid.setCurrentIndex(1)
        self.ppbGrid.setToolTip("The number of points per bohr determines the resolution of the grid on which the cubes are calculated")
        self.ppbGrid.currentIndexChanged.connect(self.recalculateCubes)
        cubeLayout.addWidget(self.ppbGrid)
                             
        calcCubesButton = QtGui.QPushButton("Calc. cubes")
        calcCubesButton.setToolTip("Calculate cubes for the selected orbitals")
        calcCubesButton.clicked.connect(self.calculateCubes)
        cubeLayout.addWidget(calcCubesButton)
        # MOs
        moFrame = QtGui.QFrame()
        upperLayout.addWidget(moFrame)
        moLayout = QtGui.QVBoxLayout(moFrame)
        self.moLabel = QtGui.QLabel("Ground State:")
        moLayout.addWidget(self.moLabel)
        moTab = QtGui.QTabWidget()
        moLayout.addWidget(moTab)
        # 3D cubes - MOs
        self.bs = AtomicBasisSet(self.tddftb.dftb2.atomlist)
        self.molecularOrbitalViewer = QCubeViewerWidget()
        self.molecularOrbitalViewer.setIsoValue(0.02)
        moTab.addTab(self.molecularOrbitalViewer, "3D Molecular Orbital")
        # 3D cubes - MOs
        # Mulliken Charges
        self.chargesViewer = QCubeViewerWidget()
        moTab.addTab(self.chargesViewer, "Mulliken Charges")
        atomlist = self.tddftb.dftb2.getGeometry()
        partial_charges = -self.tddftb.dftb2.getPartialCharges()
        self.chargesViewer.setGeometriesAndCharges([atomlist], [partial_charges])
        self.chargesViewer.selectShowOptions(options=["charges"])
        # figure with DOS
        dosFrame = QtGui.QFrame()
        layout.addWidget(dosFrame)
        dosLayout = QtGui.QVBoxLayout(dosFrame)
        self.figDOS = Figure()
        self.canvasDOS = FigureCanvas(self.figDOS)
        dosLayout.addWidget(self.canvasDOS)
        NavigationToolbar(self.canvasDOS, dosFrame, coordinates=True)
        # controls
        controlFrame = QtGui.QFrame()
        dosLayout.addWidget(controlFrame)
        controlLayout = QtGui.QHBoxLayout(controlFrame)
        self.energyUnits = QtGui.QComboBox()
        
        self.energyUnits.addItems(["Hartree", "eV", "cm-1"])
        self.energyUnits.currentIndexChanged.connect(self.plotDOS)
        controlLayout.addWidget(QtGui.QLabel("units:"))
        controlLayout.addWidget(self.energyUnits)
        controlLayout.addWidget(QtGui.QLabel("broadening:"))
        self.broadening = QtGui.QLineEdit()
        self.broadening.setText("0.01")
        self.broadening.editingFinished.connect(self.plotDOS)
        self.broadening.setToolTip("The sticks are convolved with a Gaussian to simulated a temperature broadened DOS")
        controlLayout.addWidget(self.broadening)
        
        # load KS orbitals
        tab = self.tableMOs
        dftb = self.tddftb.dftb2
        orbe = dftb.getKSEnergies()
        f = dftb.getOccupation()
        self.HOMO,self.LUMO = dftb.getFrontierOrbitals()
        tab.setRowCount(len(orbe))
        headers = ["N", "name", "occ.", "orb. en. / hartree", "orb. en. / eV", "Cube"]
        tab.setColumnCount(len(headers))
        tab.setHorizontalHeaderLabels(headers)
        self.mo_names = []
        for i in range(0, len(orbe)):
            if i == self.HOMO:
                name = "HOMO"
            elif i == self.LUMO:
                name = "LUMO"
            else:
                name = ""
            self.mo_names.append(name)
            row = [string.rjust(str(i+1),5), \
                   string.rjust(name,5), 
                   string.rjust(str(f[i]),5),
                   string.rjust("%.7f" % orbe[i],20), \
                   string.rjust("%.7f" % (orbe[i]*AtomicData.hartree_to_eV), 17) ]
            for j,r in enumerate(row):
                tab.setItem(i,j, QtGui.QTableWidgetItem("%s" % r))
        tab.resizeColumnsToContents()
        self.plotDOS()
        # cubes with MOs
        self.mo_cubes = {}
        # select HOMO
        self.tableMOs.setCurrentCell(self.HOMO, 0)
    def plotDOS(self):
        #print "plot density of states"
        self.figDOS.clf()
        ax = self.figDOS.add_subplot(111)
        
        units = self.energyUnits.itemText(self.energyUnits.currentIndex())
        broadening = abs(float(self.broadening.text()))
        ax.set_title("Density of States (DOS)")
        ax.set_xlabel("Kohn-Sham Energy / %s" % units, fontsize=15)
        ax.set_ylabel("DOS / arbitrary units", fontsize=15)

        ens = self.tddftb.dftb2.getKSEnergies()
        dos = 1.0 * np.ones(ens.shape)
        # stick spectrum
        ens_out = convert_energy(ens, "Hartree", units)
        ax.vlines(ens_out, np.zeros(ens_out.shape), dos, lw=2, color="black", picker=5)

        # selected states are highlighted
        selected = self.tableMOs.selectionModel().selectedRows()
        for s in selected:
            i = s.row()
            ax.vlines(ens_out[i], [0.0], dos[i], lw=3, color="green")
            ax.text(ens_out[i], dos[i], "%s %s" % (i+1, self.mo_names[i]))

        def onpick(event):
            self.tableMOs.setCurrentCell(event.ind[0], 0)
        self.figDOS.canvas.mpl_connect('pick_event', onpick)
            
        if broadening > 0.0:
            ens_broadened, spec = broadened_spectrum(ens, dos, broadening)
            ens_broadened_out = convert_energy(ens_broadened, "Hartree", units)
            scale = spec.max()/dos.max()/1.1
            x, y = ens_broadened_out, spec/scale
            ax.plot(x, y, lw=1, color="blue")

            x_occ = x[x <= ens_out[self.HOMO]]
            y_occ = y[x <= ens_out[self.HOMO]]
            x_virt = x[x >= ens_out[self.LUMO]]
            y_virt = y[x >= ens_out[self.LUMO]]
            
            ax.fill_between(x_occ,  0, y_occ, alpha=0.5, color="blue")
            ax.fill_between(x_virt, 0, y_virt, alpha=0.5, color="red")

            ax.annotate(
                '', xy=(ens_out[self.HOMO], 1.1), xycoords='data',
                xytext=(ens_out[self.LUMO], 1.1), textcoords='data',
                arrowprops={'arrowstyle': '<->'})
            ax.annotate(
                'gap = %4.3f' % (ens_out[self.LUMO]-ens_out[self.HOMO]),
                xy=((ens_out[self.HOMO]+ens_out[self.LUMO])/2.0, 1.1), xycoords='data',
                xytext=(0, 5), textcoords='offset points')
            
        self.canvasDOS.draw()

        self.showMullikenCharges()
        if len(selected) > 0:
            i = selected[0].row()
            self.moLabel.setText("Molecular Orbital: %s %s" % (i+1, self.mo_names[i]))
            if self.mo_cubes.has_key(i):
                mo_cube = self.mo_cubes[i]
                self.molecularOrbitalViewer.setCubes([mo_cube])
            else:
                print "No cube for molecular orbital %d" % (i+1)
    def showMullikenCharges(self):
        atomlist = self.tddftb.dftb2.getGeometry()
        dq0 = self.tddftb.dftb2.getPartialCharges()
        # internally charges are measured in units of e,
        charges = -dq0
        self.chargesViewer.setGeometriesAndCharges([atomlist], [charges])
    def calculateCubes(self):
        selected = self.tableMOs.selectionModel().selectedRows()
        for s in selected:
            i = s.row()
            if not (self.mo_cubes.has_key(i)):
                # need to calculated cube
                print "compute cube for molecular orbital %d ..." % (i+1),
                mo_cube = self._computeMO(i)
                print "...done"
                self.mo_cubes[i] = mo_cube

                self.tableMOs.setItem(i, self.tableMOs.columnCount()-1, QtGui.QTableWidgetItem("Y"))
                self.tableMOs.item(i, self.tableMOs.columnCount()-1).setBackground(QtGui.QColor(255,0,0))
        self.plotDOS()
    def recalculateCubes(self):
        # delete cubes
        for i in self.mo_cubes.keys():
                self.tableMOs.setItem(i, self.tableMOs.columnCount()-1, QtGui.QTableWidgetItem(""))
                self.tableMOs.item(i, self.tableMOs.columnCount()-1).setBackground(QtGui.QColor(255,255,255))
        self.mo_cubes = {}
        delattr(Cube.orbital_amplitude, "cached_grid") # cached grid has wrong resolution
        self.calculateCubes()
    def _computeMO(self, i):
        cube = CubeData()
        atomlist = self.tddftb.dftb2.getGeometry()
        
        (xmin,xmax),(ymin,ymax),(zmin,zmax) = Cube.get_bbox(atomlist, dbuff=5.0)
        dx,dy,dz = xmax-xmin,ymax-ymin,zmax-zmin

        ppb_text = self.ppbGrid.itemText(self.ppbGrid.currentIndex())
        ppb = float(ppb_text.split()[1]) # points per bohr
        #print "points per bohr: %s" % ppb

        nx,ny,nz = int(dx*ppb),int(dy*ppb),int(dz*ppb)
        x,y,z = np.mgrid[xmin:xmax:nx*1j, ymin:ymax:ny*1j, zmin:zmax:nz*1j]
        grid = (x,y,z)
        #
        orbs = self.tddftb.dftb2.getKSCoefficients()
        moGrid = Cube.orbital_amplitude(grid, self.bs.bfs, orbs[:,i], cache=True)
        
        mo_cube = CubeData()
        mo_cube.data = moGrid.real
        mo_cube.grid = grid
        mo_cube.atomlist = atomlist
        
        return mo_cube

    
class QdftbGUI(QtGui.QMainWindow):
    def __init__(self, tddftb, parent=None):
        QtGui.QMainWindow.__init__(self, parent)
        self.tddftb = tddftb

        main = QtGui.QTabWidget()
        moWidget = QMolecularOrbitalWidget(tddftb)
        main.addTab(moWidget, "Ground State (DFTB)")
        specWidget = QSpectrumWidget(tddftb)
        main.addTab(specWidget, "Excited States (TD-DFTB)")
        main.setCurrentIndex(1)
        self.setCentralWidget(main)
        
def start_graphical_analysis(tddftb):
    import sys
    #app = QtGui.QApplication(sys.argv)
    app = QtGui.QApplication.instance()
    window = QdftbGUI(tddftb)
    window.show()
    app.exec_()


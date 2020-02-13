#!/usr/bin/env python

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import meta

class reconstruct_metadynamics(meta.metadynamics):
    def reconstruct_vg(self, x0=None, ngauss=None, plot=True, animate=False, stride=1):
        """
        Reconstruct and plot potential from data in vg_centers.dat and meta-config.py. Grid and
        VG are saved as numpy arrays to vg_grid.npy and vg_values.npy
        
        Parameters
        ----------
        x0 : list of numpy arrays
            manually specified x-values for each collective variable (default: None)
        ngauss : number
            first ngauss gaussian functions of the simulation will be evaluated (default: all)
        plot : boolean
            plot single image of vg with the specified number of gaussians (default: True)
        animate : boolean
            animate evolution of vg with number of gaussians (default: False)
        stride : number
            number of gaussians to be added in each frame of the animation (default: 1)
        
        """
        self.centers = np.loadtxt("vg_centers.dat", ndmin=1)
        if not x0:  # if x0 not set manually, use standard ranges attributed to the CVs
            x0 = [cv.x for cv in self.cvs]
        if not ngauss:  # if ngauss not set manually, use all gaussians present in vg_centers.dat
            ngauss = len(self.centers)
        h, w = self.cfg["height"], [cv.convert_units(cv.width) for cv in self.cvs]
        if self.ncvs == 1:
            x, w, cv = x0[0], w[0], self.cvs[0]
            vg = self.get_1dvg(x, h, w, ngauss)
            if plot:
                self.plot_1dvg(vg, x, cv, ngauss)
            if animate:
                self.animate_1dvg(vg, x, h, w, cv, ngauss, stride)
        else:
            x = np.meshgrid(*x0)
            vg = self.get_ndvg(x, h, w, ngauss)
            if self.ncvs == 2:
                if plot:
                    self.plot_2dvg(vg, x, self.cvs, ngauss)
                if animate:
                    self.animate_2dvg(vg, x, h, w, self.cvs, ngauss, stride)
            else:
                print "Plotting the metadynamics potential is only supported for up to 2 CVs"
        np.save("vg_grid.npy", x)
        print "Grid written to vg_grid.npy"
        np.save("vg_values.npy", vg)
        print "Metadynamics potential written to vg_values.npy"

    def gauss(self, x, g, w):
        return np.exp(-(x-g)**2/2/w**2)

    def get_1dvg(self, x, h, w, ngauss):
        """
        Reconstruct vg on the 1D-grid x
        
        Parameters
        ----------
        x : numpy array
            x-values for the collective variable
        h : number
            height of the added gaussians
        w : number
            width of the added gaussians
        ngauss : number
            first ngauss gaussian functions of the simulation will be evaluated
        
        Returns
        -------
        vg : numpy array
            metadynamics potential on the grid x
        
        """
        cv = self.cvs[0]
        vg = np.zeros(x.size)
        for g in self.centers[:ngauss]:
            vg += h*self.gauss(x, g, w)
            if cv.type == "torsion":  # account for periodic angles
                vg += h*self.gauss(x+360, g, w)
                vg += h*self.gauss(x-360, g, w)
        return vg

    def get_ndvg(self, x, h, w, ngauss):
        """
        Reconstruct vg on the nD-grid x
        
        Parameters
        ----------
        x : list of numpy arrays
            x-values for each collective variable
        h : number
            height of the added gaussians
        w : list of numbers
            widths of the added gaussians (one per CV)
        ngauss : number
            first ngauss gaussian functions of the simulation will be evaluated
        
        Returns
        -------
        vg : numpy array
            metadynamics potential on the grid x
        
        """
        vg = np.zeros(x[0].shape)
        for g in self.centers[:ngauss]:  # add one gaussian in each step
            glist = np.array([self.gauss(x[i], g[i], w[i]) for i in range(self.ncvs)])
            vg += h*np.prod(glist, axis=0)  # combine single gaussians
            for j, cv in enumerate(self.cvs):
                if cv.type == "torsion":  # account for periodic angles
                    for shift in [-360., 360.]:
                        xs = [x_i.copy() for x_i in x]
                        xs[j] += shift
                        glist = np.array([self.gauss(xs[i], g[i], w[i]) for i in range(self.ncvs)])
                        vg += h*np.prod(glist, axis=0)  # combine single gaussians
        return vg

    def plot_1dvg(self, vg, x, cv, ngauss):
        """
        Plot 1D metadynamics potential and save to file vg_ngauss-g.png
        
        Parameters
        ----------
        vg : 1d numpy array
            metadynamics potential on the grid x
        x : 1d numpy array
            grid for the metadynamics potential
        cv : cv object
            parameters of the collective variable read from input file
        ngauss: number
            number of gaussians contained in vg
        
        """
        if cv.ticks:
            plt.xticks(np.arange(x.min(), x.max()+cv.ticks, cv.ticks))
        plt.xlabel("CV: "+cv.name)
        plt.title("Metadynamics Potential ("+str(ngauss)+" gaussians)")
        plt.plot(x, vg)
        plt.tight_layout(pad=0.1)
        plt.savefig("vg_"+str(ngauss)+"-g.png", dpi=600)
        plt.show()

    def plot_2dvg(self, vg, x, cvs, ngauss):
        """
        Plot 2D metadynamics potential and save to file vg_ngauss-g.png
        
        Parameters
        ----------
        vg: 2d numpy array
            metadynamics potential on the grid x
        x : 2d numpy meshgrid
            grid for the metadynamics potential
        cvs : list of cv objects
            list of parameters of the collective variables read from input file
        ngauss: number
            number of gaussians contained in vg
        
        """
        xmin, xmax, ymin, ymax = x[0].min(), x[0].max(), x[1].min(), x[1].max()
        asp = (xmax-xmin)/(ymax-ymin)
        plt.imshow(vg, origin='lower', aspect=asp, extent=(xmin, xmax, ymin, ymax))
        plt.colorbar()
        plt.title("Metadynamics Potential ("+str(ngauss)+" gaussians)")
        type1, type2 = [cv.type for cv in cvs]
        if cvs[0].ticks:
            plt.xticks(np.arange(xmin, xmax+cvs[0].ticks, cvs[0].ticks))  # CV 1
        if cvs[1].ticks:
            plt.yticks(np.arange(ymin, ymax+cvs[1].ticks, cvs[1].ticks))  # CV 2
        plt.xlabel("CV 1: "+cvs[0].name)
        plt.ylabel("CV 2: "+cvs[1].name)
        plt.savefig("vg_"+str(ngauss)+"-g.png", dpi=600)
        plt.show()

    def animate_1dvg(self, vg, x, h, w, cv, ngauss, stride=1):
        """
        Animate 1D metadynamics potential and save frames to directory animation
        
        Parameters
        ----------
        vg : 1d numpy array
            metadynamics potential on the grid x
        x : 1d numpy array
            grid for the metadynamics potential
        h : number
            height of the added gaussians
        w : number
            width of the added gaussians
        cv : cv object
            parameters of the collective variable read from input file
        ngauss: number
            number of gaussians contained in vg
        stride : number
            number of gaussians to be added in each frame of the animation (default: 1)
        
        """
        # animation function
        def animate(ng):
            line.set_ydata(self.get_1dvg(x, h, w, ng))
            plt.title("Metadynamics Potential (%i gaussians)" % (ng))
            plt.savefig("animation/%04i.png" % (ng))
            return line,
        # configure plot
        fig, ax = plt.subplots()
        if cv.ticks:
            plt.xticks(np.arange(x.min(), x.max()+cv.ticks, cv.ticks))
        plt.xlabel("CV: "+cv.name)
        # plot and animate
        line, = ax.plot(x, vg)
        if not os.path.isdir("animation"):
            os.mkdir("animation")
        ani = animation.FuncAnimation(fig, animate, np.arange(0, ngauss+1, stride), interval=25,
                                      repeat=False)
        plt.show()

    def animate_2dvg(self, vg, x, h, w, cvs, ngauss, stride=1):
        """
        Animate 2D metadynamics potential and save frames to directory animation
        
        Parameters
        ----------
        vg : 2d numpy array
            metadynamics potential on the grid x
        x : 2d numpy array
            grid for the metadynamics potential
        h : number
            height of the added gaussians
        w : list of numbers
            widths of the added gaussians (1 per CV)
        cvs : list of cv objects
            parameters of the collective variables read from input file
        ngauss: number
            number of gaussians contained in vg
        stride : number
            number of gaussians to be added in each frame of the animation (default: 1)
        
        """
        # animation function
        def animate(ng):
            im.set_array(self.get_ndvg(x, h, w, ng))
            plt.title("Metadynamics Potential (%i gaussians)" % (ng))
            plt.savefig("animation/%04i.png" % (ng))
            return im,
        # configure plot
        fig, ax = plt.subplots()
        xmin, xmax, ymin, ymax = x[0].min(), x[0].max(), x[1].min(), x[1].max()
        asp = (xmax-xmin)/(ymax-ymin)
        type1, type2 = [cv.type for cv in cvs]
        if cvs[0].ticks:
            plt.xticks(np.arange(xmin, xmax+cvs[0].ticks, cvs[0].ticks))  # CV 1
        if cvs[1].ticks:
            plt.yticks(np.arange(ymin, ymax+cvs[1].ticks, cvs[1].ticks))  # CV 2
        plt.xlabel("CV 1: "+cvs[0].name)
        plt.ylabel("CV 2: "+cvs[1].name)
        # plot and animate
        im = plt.imshow(vg, origin='lower', aspect=asp, extent=(xmin, xmax, ymin, ymax))
        plt.colorbar()
        if not os.path.isdir("animation"):
            os.mkdir("animation")
        ani = animation.FuncAnimation(fig, animate, np.arange(0, ngauss+1, stride), interval=25,
                                      repeat=False)
        plt.show()


if __name__ == "__main__":
    import argparse
    # initialize parser
    parser = argparse.ArgumentParser(description="Reconstruct the metadynamics potential using \
                                     data in vg_centers.dat and parameters in meta-config.py")
    parser.add_argument("--ngauss", type=int, default=None,
                        help="use only the first ngauss gaussians (default: all)")
    parser.add_argument("--plot", type=bool, default=True,
                        help="plot the potential and save png (default: True)")
    parser.add_argument("--animate", type=bool, default=False,
                        help="animate the potential evolution (default: False)")
    parser.add_argument("--stride", type=int, default=1,
                        help="number of gaussians per animation frame (default: 1)")
    args = parser.parse_args()
    ngauss = vars(args)["ngauss"]
    doplot = vars(args)["plot"]
    doanimate = vars(args)["animate"]
    stride = vars(args)["stride"]
    # do reconstruction
    metad = reconstruct_metadynamics()
    metad.read_input()
    metad.reconstruct_vg(ngauss=ngauss, plot=doplot, animate=doanimate, stride=stride)

#!/usr/bin/env python
"""
broaden a stick spectrum with a line shape.
"""
from DFTB import AtomicData
import os.path

import numpy as np
from matplotlib import pyplot as plt

def line_shape(nu, nu0, FWHM):
    # Lorentzian distribution
    return ( 1.0/np.pi * (FWHM/2.0) / ((nu-nu0)**2 + (FWHM/2.0)**2) )

def convert_energy(en_in, in_units, out_units):
    en = np.copy(en_in)
    # first convert en to Hartree
    if in_units == "Hartree":
        pass
    elif in_units == "eV":
        en /= AtomicData.hartree_to_eV
    elif in_units == "nm":
        en = AtomicData.hartree_to_nm / en
    elif in_units == "cm-1":
        en /= AtomicData.hartree_to_wavenumbers 
    else:
        raise NotImplemented("Unknown energy units: %s" % in_units)
    # convert to output units
    if out_units == "Hartree":
        pass
    elif out_units == "eV":
        en *= AtomicData.hartree_to_eV
    elif out_units == "nm":
        en = AtomicData.hartree_to_nm / en
    elif out_units == "cm-1":
        en *= AtomicData.hartree_to_wavenumbers
    else:
        raise NotImplemented("Unknown energy units: %s" % out_units)
    return en

def broadened_spectrum(energies, oscillator_strengths, broadening):
    """
    energies: numpy array with excitation energies in a.u.
    oscillator_strengths:
    broadening: line broadening in Hartree

    Returns:
    ========
    x,y: x are the energies in Hartree, eV or nm depending on the keyword units
    """
    if broadening > 0.0:
        emin = energies.min() #max(energies.min() - 0.01, 1.0e-10)
        emax = energies.max() + 0.01
        energy_range = np.linspace(emin, emax, 5000)
        dE = np.ediff1d(energy_range, to_end=energy_range[-1]-energy_range[-2])

        absorption_spectrum = 0.0*energy_range
        for en,f in zip(energies, oscillator_strengths):
            peak = line_shape(energy_range, en, broadening)
#            integ_intensity = sum(peak*dE)
            absorption_spectrum += f * peak * dE
    else:
        energy_range = energies
        absorption_spectrum = oscillator_strengths

    return energy_range, absorption_spectrum

if __name__ == "__main__":
    import sys
    from os.path import expandvars, expanduser, basename
    from optparse import OptionParser
    usage = "Usage: python %s <list of files with tabulated absorption spectrum>\n" % basename(sys.argv[0])
    usage += "  use --help to see all options"
    
    parser = OptionParser(usage)
    parser.add_option("--save_png", dest="save_png", help="save plot to this file [default: %default]", default="")
    parser.add_option("--broadening", dest="broadening", help="convolve stick spectrum with a Lorentzian of FWHM=broadening (in Hartree) [default: %default]", default=0.005)
    parser.add_option("--in_units", dest="in_units", default="Hartree", help="input units, 'Hartree', 'eV', 'nm', 'cm-1' [default: %default]")
    parser.add_option("--out_units", dest="out_units", default="eV", help="output units, 'Hartree', 'eV', 'nm', 'cm-1' [default: %default]")
    parser.add_option("--plot_range", dest="plot_range", default="", help="plot range in output units, e.g. (0.0, 6.0) to plot from 0 to 6 eV [default: %default]")
    parser.add_option("--plot_labels", dest="plot_labels", default=1, help="show file names as plot labels [default: %default]")
    parser.add_option("--plot_sticks", dest="plot_sticks", type=int, default=1, help="plot stick spectrum, too [default: %default]")    
    parser.add_option("--vertical", dest="vertical", default=0, type=int, help="use separate subplot for each spectrum [default: %default]")
    parser.add_option("--invert_xaxis", dest="invert_xaxis", default=0, type=int, help="Invert energy-axis so that it points from right to left [default: %default")
    parser.add_option("--save_converted", dest="save_converted", default="", help="convert the spectrum to out_units and save it to this file [default: %default]")
    parser.add_option("--save_convolved", dest="save_convolved", default="", help="convert the spectrum to out_units, convole it and save it to this file [default: %default]")
    parser.add_option("--transient_weights", dest="transient_weights", default="", help="Produce a transient absorption spectrum. Each spectrum is weighted by the time-dependent concentration of the intermediate to which it belongs. The concentrations are read from the data file <transient_weight>, the first column in the file is the reaction time followed by as many columns as there are spectra.")
    parser.add_option("--transient_difference", dest="transient_difference", default=0, type=int, help="For transient spectra: subtract spectrum at time t=0 from all subsequent spectra.")
    parser.add_option("--labels", dest="labels", default="", type=str, help="List of labels to be used for each absorption spectrum [default: %default]")
    parser.add_option("--log_plot", dest="log_plot", default=0, type=int, help="To sho states with very low oscillator strengths plot the y-axis in logarithmic scale [default: %default]")
    parser.add_option("--normalize", dest="normalize", action="store_true", help="Normalize spectrum so that maximum has intensity of 1 [default: %default]")
    (opts, args) = parser.parse_args()
    if len(args) < 1:
        print usage
        exit(-1)

    spectra = []
    labels = []
    # load all spectra
    for i,spec_file in enumerate(args):
        data = np.loadtxt(spec_file).transpose()
        # If can happend that there is only a single line in the file
        # adjust shape
        if len(data.shape) == 1:
            data = np.reshape(data, (len(data),1))
        spectra.append(data)
    # labels
    if opts.labels == "":
        for i,spec_file in enumerate(args):
            labels.append( os.path.basename(spec_file) )
    else:
        labels = eval(opts.labels)
    if opts.transient_weights != "":
        # transient absorption spectrum
        concentrations = np.loadtxt(expandvars(expanduser(opts.transient_weights)))
        ts = concentrations[:,0]
        cs = concentrations[:,1:]
        assert len(spectra) == cs.shape[1], "The number of spectrum files has to match the number of concentrations in the file <transient weights>! #spectra=%s != #weight=%s" % (len(spectra),cs.shape[1])
        transient_spectra = []
        transient_labels = []
        for itime,t in enumerate(ts):
            print "Time step: %s" % itime
            i_en = []
            i_osz = []
            i_label = "%2.3f" % t
            for iweight,weight in enumerate(cs[itime,:]):
                i_en += list(spectra[iweight][0])
                i_osz += list(weight * spectra[iweight][1])
            data = np.array([i_en, i_osz])
            print data.shape
            transient_spectra.append(data)
            transient_labels.append(i_label)
        spectra = transient_spectra
        labels = transient_labels

    for i,data in enumerate(spectra):
        if opts.vertical == 1:
            if i == len(spectra)-1:
                shared_axes = plt.subplot(len(spectra),1,i+1)
                shared_axes.yaxis.set_ticklabels([])
            else:
                axes = plt.subplot(len(spectra),1,i+1)
                axes.xaxis.set_visible(False)
                axes.yaxis.set_ticklabels([])
        print "spectrum %s" % i
        en_in, osz = data[0], data[1]
        # sort energies
        sort_indx = np.argsort(en_in)
        en_in = en_in[sort_indx]
        osz = osz[sort_indx]
        #
        en = convert_energy(en_in, opts.in_units, "Hartree")
        if opts.save_converted != "":
            print "Save converted spectrum to %s" % opts.save_converted
            en_conv = convert_energy(en, "Hartree", opts.out_units)
            np.savetxt(expandvars(expanduser(opts.save_converted)), np.vstack((en_conv, osz)).transpose())
            exit(0)
        # plot broadened spectrum
        ens, spec = broadened_spectrum(en, osz, float(opts.broadening))
        ens_out = convert_energy(ens, "Hartree", opts.out_units)
        # normalize
        if opts.normalize == True:
            if opts.plot_range != "":
                emin, emax = eval(opts.plot_range)
                max_intensity = spec[(ens_out > emin) & (ens_out < emax)].max()
            else:
                max_intensity = spec.max()
            spec /= max_intensity
        if opts.save_convolved != "":
            print "Save convolved spectrum to %s" % opts.save_convolved
            np.savetxt(expandvars(expanduser(opts.save_convolved)), np.vstack((ens_out, spec)).transpose())
            exit(0)

        if int(opts.plot_labels) == 1:
            label=labels[i]
        else:
            label = ""
        if float(opts.broadening) > 0.0:
            l, = plt.plot(ens_out, spec, lw=2, label="%s" % label )
        else:
            l, = plt.plot(ens_out, spec, "o", label="%s" % label )
        # plot stick spectrum
        if opts.plot_sticks == 1:
            scale = spec.max()/osz.max()/1.3
            en_out = convert_energy(en, "Hartree", opts.out_units)
            if len(args) > 1:
                clr = plt.getp(l, "color")
            else:
                clr = "black"
            plt.vlines(en_out, np.zeros(en_out.shape), osz*scale, color=clr)
        if opts.vertical == 0:
            ax = plt.gca()
            loglabel = ""
            if opts.log_plot == 1:
                ax.set_yscale('log')
                loglabel = " (log.)"
            # make labels larger if there is only one plot
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(15) 
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(15) 
            plt.ylabel("oscillator strengths %s" % loglabel, fontsize=18)
        else:
            plt.ylabel("", fontsize=12)            
        plt.xlabel(opts.out_units, fontsize=18)

        if opts.plot_range != "":
            plt.xlim(eval(opts.plot_range))
        if opts.plot_labels == 1:
            plt.legend()
    if opts.vertical == 1:
        plt.subplots_adjust(wspace=0.0, hspace=0.0)
    if opts.invert_xaxis == 1:
        ax = plt.gca()
        ax.invert_xaxis()


    if opts.save_png != "":
        plt.savefig(expandvars(expanduser(opts.save_png)))
    else:
        plt.show()
    

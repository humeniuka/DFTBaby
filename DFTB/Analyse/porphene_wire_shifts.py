"""
In porphene wires of increasing length the shift of the B- and Q-bands with increasing length can be explained by two mechanisms:

1. The shift of the lowest excited state in a wire of meso-meso,beta-beta,beta-beta fused
porphyrins can be explained with the Hueckel model for linear conjugated polyenes.

The energy of the first excited state roughly equals the HOMO-LUMO gap, which for
a linear conjugated polymer consisting of N monomer units is given by:

E(HOMO-LUMO) = -4 beta sin(pi/(2 (N+1) )) ~ -4 beta pi/2 1/(N+1)

where beta = <a|H|a+1> is the interaction between the relevant orbital on adjacent monomers.

2. The splitting of the B-band is due to excitonic coupling
"""
import numpy as np
import numpy.linalg as la
from matplotlib import pyplot as plt


from DFTB import AtomicData

def extract_shifts(spec, enB0, osc_threshold):
    """
    Find the shift of the Q-band and the splitting of the B-band

    Parameters:
    ===========
    spec: 2d numpy array spec[i,:] = (en_i, osc_i)  contains the energy
      and oscillator strength of the i-th state
    enB0: energy of the Soret peak in the monomer
    """
    en = spec[:,0]
    f = spec[:,1]
    # the lowest excitation energy with oscillator strength higher than 
    # osc_threshold is taken to be the lowest Q-band peak
    for i in range(0, len(en)):
        if f[i] > osc_threshold:
            Q_indx = i
            break
    # the brightest peak below enB0 is taken to be the second peak of the split Soret band
    Bsplit_indx = np.argmax(f[en < enB0])
    enQ = en[Q_indx]
    exciton_splitting = en[Bsplit_indx] - enB0
    return enQ, exciton_splitting

def fit_hueckel(enQ, Nmax=12):
    wavelengthQ = AtomicData.hartree_to_nm/enQ

    N = np.array(range(2, len(enQ)+1))
    # fit a line m*x + c through the data
    x = N
    # monomer is not included
    y = wavelengthQ[1:]

    xfit = x[:Nmax]
    yfit = y[:Nmax]
    A = np.vstack([xfit, np.ones(len(xfit))]).transpose()
    m, c = la.lstsq(A, yfit)[0]

    plt.xlabel("N", fontsize=15)
    plt.ylabel("$\lambda_{max}$ / nm", fontsize=15)
    
    plt.plot(N, y, "o", label="Theory (lc-TD-DFTB)")
    plt.plot(xfit, m*(xfit)+c, ls="-", lw=2, label=u"fit to H\u00FCckel Model, $%3.0f N + %3.0f$ (nm)" % (m,c))
    plt.legend(loc="lower center")
    plt.show() 

def fit_exciton_splitting(splitting_energies, Nmax=4):
    N = np.array(range(2, len(splitting_energies)+1))
    # fit a line m*x + c through the data
    x = np.cos(np.pi/(N+1))
    # monomer is not included
    y = np.array(splitting_energies[1:])

    xfit = x[:Nmax]
    yfit = y[:Nmax]
    A = np.vstack([xfit, np.ones(len(xfit))]).transpose()
    m, c = la.lstsq(A, yfit)[0]
    print "Exciton splitting energy Delta E_0 = %3.7f Hartree    %3.7f eV    %3.7f nm     %3.7f cm^(-1)" % (m, m*AtomicData.hartree_to_eV, AtomicData.hartree_to_nm/m, m*AtomicData.hartree_to_wavenumbers)

    plt.xlabel("$cos(\pi/(N+1))$", fontsize=15)
    plt.ylabel("exciton splitting $\Delta E$ / eV", fontsize=15)

    m *= AtomicData.hartree_to_eV
    c *= AtomicData.hartree_to_eV

    plt.xlim((0.45,0.90))
    plt.ylim((-1.5,-0.8))

    plt.plot(xfit, yfit * AtomicData.hartree_to_eV, "o", label="Theory (lc-TD-DFTB)")
    plt.plot(xfit, m*xfit+c, ls="-", lw=2, label="fit to Exciton Model $%3.2f \cos(\pi/(N+1)) %+3.2f$ (eV)" % (m,c))
    plt.legend(loc="upper center", bbox_to_anchor=(0.5,1.05))
    plt.show() 
    

if __name__ == "__main__":
    import sys
    from os.path import expandvars, expanduser
    from optparse import OptionParser
    usage = "Usage: python %s <list of files with tabulated absorption spectrum>" % sys.argv[0]
    parser = OptionParser(usage)
    parser.add_option("--osc_threshold", dest="osc_threshold", help="The Q-band is assigned to the lowest peak having oscillator strength greater than this threshold [default: %default]", type=float, default=0.02)
    parser.add_option("--out_file", dest="out_file", help="Save energies of lowest Q-peak and split Soret peak to this file [default: %default]", default="/tmp/QB_energies.txt")
    
    (opts, args) = parser.parse_args()
    if len(args) < 1:
        print usage
        exit(-1)

    spectra = []
    # load all spectra
    for i,spec_file in enumerate(args):
        data = np.loadtxt(spec_file)
        spectra.append(data)

    enQs = []
    exciton_splittings = []
    for i,spec in enumerate(spectra):
        if i == 0:
            # extract position of Soret band from monomer
            en = spec[:,0]
            f = spec[:,1]
            soret_indx = np.argmax(f[en < (5.0 / AtomicData.hartree_to_eV)])
            soret_peak = en[soret_indx]
            print "Soret peak in the monomer: %3.7f hartree     %3.7f eV" % (soret_peak, soret_peak * AtomicData.hartree_to_eV)
        enQ, dExciton = extract_shifts(spec, soret_peak, opts.osc_threshold)
        enQs.append(enQ)
        exciton_splittings.append(dExciton)
        print "%3.7f   %3.7f" % (enQ, dExciton) 
    enQs = np.array(enQs)
    data = np.vstack((enQs, exciton_splittings)).transpose()

    fit_hueckel(enQs)
    fit_exciton_splitting(exciton_splittings)

    if opts.out_file != "":
        np.savetxt(opts.out_file, data)

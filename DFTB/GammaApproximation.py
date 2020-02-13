import numpy as np
import numpy.linalg as la
from scipy.special import erf 
from scipy import interpolate

from DFTB import AtomicData

def gaussian_decay(hubbard_U, valorbs, Rlr, lc_flag, verbose=0):
    """
    The decay constants for the gaussian charge fluctuations
    are determined from the requirement d^2 E_atomic/d n^2 = U_H.

    In the DFTB approximations with long-range correction one has

     U_H = gamma_AA - 1/2 * 1/(2*l+1) gamma^lr_AA

    where l is the angular momentum of the highest valence orbital

    see "Implementation and benchmark of a long-range corrected functional
         in the DFTB method" by V. Lutsker, B. Aradi and Th. Niehaus

    Here, this equation is solved for sigmaA, the decay constant
    of a gaussian.

    Parameters:
    ===========
    hubbard_U: dictionary with hubbard parameters for each atomic number
    valorbs: dictionary with list of valence orbitals (n,l,m) for each atomic number
    Rlr: long-range radius in bohr
    lc_flag: flag for turning long-range correction on (1) or off (0)

    Returns:
    ========
    sigma: dictionary with decay constants for each atomic number
    """
    sigmas0 = {}
    for Z,U in hubbard_U.iteritems():
        sigmas0[Z] = 1.0/(np.sqrt(np.pi)*U)

    sigmas = {}
    if lc_flag == 0:
        # no long-range correction
        sigmas = sigmas0
    else:

        # sigma is one of the roots of the
        # quartic polynomial
        # p(x) = a*x^4 + b*x^3 + c*x^2 + d*x + e
        for Z,U in hubbard_U.iteritems():
            # angular momentum of highest shell
            lH = max([l for (n,l,m) in valorbs[Z]])
            a = np.pi*U**2
            b = -2*np.sqrt(np.pi)*U
            c = 1+0.25*np.pi*U**2*Rlr**2 - 0.25/(2*lH+1)**2
            d = -0.5*np.sqrt(np.pi)*U*Rlr**2
            e = 0.25*Rlr**2
            x0s = np.roots([a,b,c,d,e])
            x0s_re = x0s[x0s.imag == 0.0].real
            # More than one solution - the minimum is the right one
            sigmas[Z] = min(x0s_re[x0s_re > 0.0])
#            sigmas[Z] = max(x0s_re[x0s_re > 0.0])
            # check solutions
            x = sigmas[Z]
            px = a*x**4+b*x**3+c*x**2+d*x+e
            assert abs(px) < 1.0e-12, "Error in root of 4th order polynomial"
            # 
            err = abs(U - (1.0/np.sqrt(np.pi) * (1.0/sigmas[Z] - 1.0/(2*(2*lH+1)*np.sqrt(sigmas[Z]**2 + 0.25*Rlr**2)))))
            assert err < 1.0e-10, "Error in sigma-calculation: err = %s" % err
            assert sigmas[Z] > 0.0

        # This works better
        sigmas = sigmas0

    if verbose > 0:
        print "Gaussian decay constants"
        print "========================================================="
        print "Atom      U/Hartree      sigma/bohr    sigma/bohr (no lc)"
        for Z,U in hubbard_U.iteritems():
            # for comparison: sigma without lc-correction
            print " %2.s         %2.4f         %2.4f          %2.4f" % (AtomicData.atom_names[Z-1],U,sigmas[Z],sigmas0[Z])   

    return sigmas

def slater_decay(hubbard_U, valorbs, Rlr, lc_flag, verbose=0):
    """
    same for Slater functions
    """
    taus = {}
    if lc_flag == 0:
        for Z,U in hubbard_U.iteritems():
            taus[Z] = 16.0/5.0 * U
    else:
        # TODO: solve eqn. (58) from Lutsker's LC-DFTB article for tau
#        raise NotImplemented("")
        # Meanwhile use the same taus
        for Z,U in hubbard_U.iteritems():
            taus[Z] = 16.0/5.0 * U        
    return taus

class GammaFunction: # base class
    """
    gamma_AB = int F_A(r-RA) * 1/|RA-RB| * F_B(r-RB) d^3r
    """
    def __init__(self, sigmas):
        """
        Parameters:
        ===========
        sigmas: dictionary sigmas[Z] should give the decay parameter belonging
          to an atom with atomic number Z
        """
        self.sigmas = sigmas
    def g(self, R, ZA, ZB):
        """evaluate gamma(R)"""
        pass
    def g_deriv(self, R, ZA, ZB):
        """evaluate the derivative d/dR gamma(R)"""
        pass
    def g_limit0(self, Z):
        """evaluate the limit R->0 of gamma(R)"""
        pass

class GammaGaussian(GammaFunction):
    """
    spherical charge fluctuations are modelled as Gaussians 

       FA(|r-RA|) = 1/(2 pi sigmaA^2)^(3/2) * exp( - (r-RA)^2/(2*sigmaA^2) )

    """
    def __init__(self, sigmas):
        """construct the C_AB matrix"""
        self.sigmas = sigmas
        self.C = {}
        Zs = sigmas.keys()
        for ZA in Zs:
            for ZB in Zs:
                self.C[(ZA,ZB)] = 1.0/np.sqrt(2*(sigmas[ZA]**2 + sigmas[ZB]**2))
    def g(self, R, ZA, ZB):
        assert R > 0.0
        gAB = erf(self.C[(ZA,ZB)]*R)/R
        return gAB
    def g_deriv(self, R, ZA, ZB):
        if abs(R) < 1.0e-5:
            # R -> 0 limit
            gd = 0.0
        else:
            c = self.C[(ZA,ZB)]
            gd = 2.0*c/np.sqrt(np.pi) * np.exp(-(c*R)**2)/R - erf(c*R)/R**2
        return gd
    def g_limit0(self, Z):
        g0 = 1.0/(np.sqrt(np.pi) * self.sigmas[Z])
        return g0

class GammaSlater(GammaFunction):
    """
    spherical charge fluctuation around an atom A are modelled as Slater functions

      FA(|r-RA|) = tA^3/(8 pi)*exp(-tA*|r-RA|)

    """
    def __init__(self, sigmas):
        # for confusion, I call the decay constants tau in the case of Slater functions
        self.tau = sigmas
    def g(self, R, ZA,ZB):
        tA = self.tau[ZA]
        tB = self.tau[ZB]
        if abs(R) < 1.0e-5:
            # R -> 0 limit
            gAB = tA*tB*(tA**2 + 3*tA*tB + tB**2)/(2.0*(tA+tB)**3)
        elif abs(tA-tB) < 1.0e-5:
            # tA == tB limit
            tau = tA
            x = tau*R
            gAB = 1.0/R * (1.0 - np.exp(-tau*R) * (48+33*x+9*x**2+x**3) / 48.0)
        else:
            # general case, R != 0 and tA != tB
            denomAB = tB**4*(tB**2*(2+tA*R)-tA**2*(6+tA*R))
            denomBA = tA**4*(tA**2*(2+tB*R)-tB**2*(6+tB*R))
            num = 2*(tA**2-tB**2)**3
            gAB = 1.0/R * (1.0 + (np.exp(-tA*R)*denomAB - np.exp(-tB*R)*denomBA)/num)
        return gAB
    def g_deriv(self, R, ZA,ZB):
        tA = self.tau[ZA]
        tB = self.tau[ZB]
        if abs(R) < 1.0e-5:
            # R -> 0 limit
            gd = 0.0
        elif abs(tA-tB) < 1.0e-5:
            # tA == tB limit
            x = tA*R
            gd = -1.0/R**2 * (1 \
                    - np.exp(-x)*(1+1.0/48.0*(x*(4+x)*(12+x*(3+x)))))
        else:
            # general case
            tA = self.tau[ZA]
            tB = self.tau[ZB]
            tAR = tA*R
            tBR = tB*R
            tA2 = tA**2
            tB2 = tB**2
            denom = 2*(tA2-tB2)**3
            fB = (2+tBR*(2+tBR))*tA2 - (6+tBR*(6+tBR))*tB2
            fA = (2+tAR*(2+tAR))*tB2 - (6+tAR*(6+tAR))*tA2
            gd = -1.0/R**2 * (1 \
                   - 1.0/denom * (tA2**2*fB*np.exp(-tBR) - tB2**2*fA*np.exp(-tAR)))
        return gd
    def g_limit0(self, Z):
        g0 = 5.0/16.0 * self.tau[Z]
        return g0

class GammaGaussianLC(GammaFunction):
    """
    spherical charge fluctuations are modelled as Gaussians 

       FA(|r-RA|) = 1/(2 pi sigmaA^2)^(3/2) * exp( - (r-RA)^2/(2*sigmaA^2) )

    long-range part of the Coulomb potential
      
       1/r  ->  erf(r/Rlr)/r
    """
    def __init__(self, sigmas, Rlr):
        """construct the C_AB matrix"""
        self.sigmas = sigmas
        self.Rlr = Rlr
        #
        self.C = {}
        Zs = sigmas.keys()
        for ZA in Zs:
            for ZB in Zs:
                self.C[(ZA,ZB)] = 1.0/np.sqrt(2*(self.sigmas[ZA]**2 + self.sigmas[ZB]**2 + 0.5*Rlr**2))
    def g(self, R, ZA, ZB):
        assert R > 0.0
        gAB = erf(self.C[(ZA,ZB)]*R)/R
        return gAB 
    def g_deriv(self, R, ZA, ZB):
        assert R > 0.0
        c = self.C[(ZA,ZB)]
        gd = 2.0*c/np.sqrt(np.pi) * np.exp(-(c*R)**2)/R - erf(c*R)/R**2
        return gd
    def g_limit0(self, Z):
        g0 = 1.0/(np.sqrt(np.pi * (self.sigmas[Z]**2 + 0.25*self.Rlr**2)))
#        # This works better, or maybe not?
#        g0 = 0.0
        return g0

class GammaLC_approx(GammaFunction):
    """ approximate LC integral by taking switching function out of the integral"""
    def __init__(self, gamma_func, switching_func):
        self.gf = gamma_func
        self.sw = switching_func
    
    def g(self, R, ZA, ZB):
        # gamma_lr ~ gamma * f(R)
        gAB = self.gf.g(R,ZA,ZB) * self.sw.f(R)
        return gAB
    def g_deriv(self, R, ZA, ZB):
        gd = self.gf.g_deriv(R,ZA,ZB) * self.sw.f(R) \
            +self.gf.g(R,ZA,ZB) * self.sw.f_deriv(R)
        return gd
    def g_limit0(self, Z):
        return 0.0

class GammaNumerical(GammaFunction):
    """
    The gamma functions are read from a table and interpolated
    """
    def __init__(self, atompairs):
        # import precalculated gamma functions
        #from DFTB.SlaterKoster.confined_pseudo_atoms import gamma_integrals
        from DFTB.SlaterKoster.free_pseudo_atoms import gamma_integrals
        
        # tabulated distances
        rg = gamma_integrals.r

        # dictionary holding splined gamma-functions for atom pairs (Za,Zb)
        self.gamma_splines_dic = {}
        self.gamma_deriv_splines_dic = {}
        
        for (Za,Zb) in atompairs:
            # gamma function on grid
            try:
                gamma_dic = gamma_integrals.gamma_integrals[(Za,Zb)]
            except KeyError as e:
                print "ERROR: No numerical gamma-function found for atom pair %s-%s." % (AtomicData.atom_names[Za-1], AtomicData.atom_names[Zb-1])
                print "You have to add these atoms in the script 'generate_gamma_integrals.py' and recompute."
                raise e
            # There are different gamma-functions for interactions between s- and p-shells, but
            # their values are quite similar. Since we want to have a single gamma function per
            # atom combination we need to average over gamma functions for different shells.
            gamma_g = 0*rg
            for (la,lb),gamma_ab in gamma_dic.iteritems():
                gamma_g += gamma_ab
            gamma_g /= len(gamma_dic)
            
            # interpolate gamma function
            tck = interpolate.splrep(rg, gamma_g)
            rmax = max(rg)

            # functions for evaluating gamma-function and its derivatives for arbitrary distances
            def gamma_AB(r):
                rflat = r.ravel()
                # set gamma(r) to 1/r for r > rmax
                gamma_r = np.reshape(np.where(rflat<=rmax, \
                                              interpolate.splev(rflat,tck,der=0), 1.0/rflat), r.shape)
                return gamma_r

            def gamma_AB_deriv(r):
                # derivative d(gamma)/dr
                rflat = r.ravel()
                # set dg/dr(r) to -1/r^2 for r > rmax
                dgamma_dr = np.reshape(np.where(rflat<=rmax, \
                                                interpolate.splev(rflat,tck,der=1), -1.0/rflat**2), r.shape)
                return dgamma_dr

            # save functions in dictionary
            self.gamma_splines_dic[(Za,Zb)] = gamma_AB
            self.gamma_deriv_splines_dic[(Za,Zb)] = gamma_AB_deriv
            
    def g(self, R, Za, Zb):
        Za,Zb = np.sort([Za,Zb])  # ensure Za <= Zb
        if type(R) == np.ndarray:
            gAB = self.gamma_splines_dic[(Za,Zb)](R)
        else:
            # turn single float into array
            R = np.array([R])
            gAB = self.gamma_splines_dic[(Za,Zb)](R)[0]
        return gAB
    def g_deriv(self, R, Za, Zb):
        Za,Zb = np.sort([Za,Zb])  # ensure Za <= Zb
        if type(R) == np.ndarray:
            gd = self.gamma_deriv_splines_dic[(Za,Zb)](R)
        else:
            # turn single float into array
            R = np.array([R])
            gd = self.gamma_deriv_splines_dic[(Za,Zb)](R)[0]
        return gd
    def g_limit0(self, Z):
        return self.g(0.0, Z, Z)
        
        
class SwitchingFunction: # base class
    def f(self, R):
        pass
    def f_deriv(self, R):
        pass

class ErfSwitching(SwitchingFunction):
    """ f(R) = erf(R/Rlr)"""
    def __init__(self, Rlr):
        self.Rlr = Rlr
    def f(self, R):
        sw = erf(R/self.Rlr)
        return sw
    def f_deriv(self, R):
        dsw = 2.0/(np.sqrt(np.pi)*self.Rlr) * np.exp(-(R/self.Rlr)**2)
        return dsw

class ErfGauSwitching(SwitchingFunction):
    """ f(R) = erf(R/Rlr) - 2/(sqrt(pi)*Rlr) * R * exp(-1/3 * (R/Rlr)**2) """
    def __init__(self, Rlr):
        self.Rlr = Rlr
    def f(self, R):
        if R < 1.0e-8:
            return 0.0
        Rlr = self.Rlr
        sw = erf(R/Rlr)/R \
            - 2.0/(np.sqrt(np.pi)*Rlr)*np.exp(-1.0/3.0 * (R/Rlr)**2)
        return sw
    def f_deriv(self, R):
        if R < 1.0e-8:
            return 0.0
        Rlr = self.Rlr
        r2 = pow(R/Rlr,2)
        sw_deriv = 4.0/(3.0*np.sqrt(np.pi)*Rlr**3) * np.exp(-r2/3.0)*R \
                  +2.0/(np.sqrt(np.pi)*Rlr) * np.exp(-r2)/R \
                  -erf(R/Rlr)/R**2
        return sw_deriv

class NoSwitching(SwitchingFunction):
    """ f(R) = 0"""
    def f(self, R):
        return 0
    def f_deriv(self, R):
        return 0

class GammaMatrix:
    def __init__(self, gamma_func):
        self.gamma_func = gamma_func
    def gamma_atomwise(self, atomlist, distances):
        """build the gamma-matrix and its gradients"""
        Nat = len(atomlist)
        g0 = np.zeros((Nat,Nat))  # gamma
        g1 = np.zeros((Nat,Nat))  # derivative of gamma
        for i,(Zi,posi) in enumerate(atomlist):
            for j,(Zj,posj) in enumerate(atomlist):
                if i == j:
                    g0[i,j] = self.gamma_func.g_limit0(Zi)
                    g1[i,j] = 0.0
                elif i < j:
                    Rij = distances[i,j]
                    g0[i,j] = self.gamma_func.g(Rij, Zi, Zj)
                    g1[i,j] = self.gamma_func.g_deriv(Rij, Zi, Zj)
                else:
                    g0[i,j] = g0[j,i]
                    g1[i,j] = g1[j,i]
        return g0, g1
    def gamma_AOwise(self, atomlist, valorbs, distances, directions):
        """
        form gamma matrix in AO basis
          gamma_ab = sum_A sum_B delta(a in A) delta(b in B) gamma_AB
        """
        g0,g1 = self.gamma_atomwise(atomlist, distances)
        Nat = len(atomlist)
        # count valence orbitals
        Norb = 0
        for i,(Zi,posi) in enumerate(atomlist):
            Norb += len(valorbs[Zi])

        g0_AO = np.zeros((Norb,Norb))
        g1_AO = np.zeros((3*Nat,Norb,Norb))

        # iterate over atoms
        mu = 0
        for i,(Zi,posi) in enumerate(atomlist):
            # iterate over orbitals on center i
            for (ni,li,mi) in valorbs[Zi]:
                # iterate over atoms
                nu = 0
                for j,(Zj,posj) in enumerate(atomlist):
                    # iterate over orbitals on center j
                    for (nj,lj,mj) in valorbs[Zj]:
                        g0_AO[mu,nu] = g0[i,j]
                        g0_AO[nu,mu] = g0[i,j]
                        if i != j:
                            eij = directions[i,j,:]
                            #
                            g1_AO[3*i:3*i+3,mu,nu] += g1[i,j] * eij
                            g1_AO[3*i:3*i+3,nu,mu] += g1[i,j] * eij

                        nu += 1
                mu += 1
        g1_AT = np.zeros((3*Nat,Nat,Nat))
        for i in range(0, Nat):
            for j in range(0, Nat):
                eij = directions[i,j,:]
                g1_AT[3*i:3*i+3,i,j] = g1[i,j] * eij
                
        return g0, g1_AT, g0_AO, g1_AO

def atomwise2AOwise(gamma, grad_gamma, atomlist, valorbs):
    """
    expand gamma matrix g_{i,j}, whose indeces i,j refer to atoms,
    to an enlarged gamma matrix gAO_{mu,mu}, whose indeces mu,nu 
    refer to valence orbitals on the atoms i,j.

    Parameters
    ----------
    gamma     :  gamma[i,j] is the (screened) Coulomb interaction between
                 charges on atoms i and j
    grad_gamma:  gradient of gamma matrix w/r/t atom positions,
                 grad_gamma[3*k:3*(k+1),i,j] = d(gamma[i,j])/dRk
    atomlist  :  atomic geometry, list of tuples (Zi,[x,y,z])
    valorbs   :  dictionary (keys are atomic numbers) with lists of the quantum numbers (n,l,m)
                 of the atomic orbitals for each atom

    Returns
    -------
    gAO       :  gAO[mu,nu] is the (screened) Coulomb interaction between charges
                 in the valence orbitals mu and nu
    grad_gAO  :  gradient of gAO w/r/t nuclear coordinates,
                 grad_gAO[3*k:3*(k+1),mu,nu] = g(gAO[mu,nu])/dRk
    
    """
    Nat = len(atomlist)
    # count valence orbitals
    Norb = 0
    for i,(Zi,posi) in enumerate(atomlist):
        Norb += len(valorbs[Zi])

    # orbital-wise gamma matrix and its gradients
    g0_AO = np.zeros((Norb,Norb))
    g1_AO = np.zeros((3*Nat,Norb,Norb))

    # iterate over atoms
    mu = 0
    for i,(Zi,posi) in enumerate(atomlist):
        # iterate over orbitals on center i
        for (ni,li,mi) in valorbs[Zi]:
            # iterate over atoms
            nu = 0
            for j,(Zj,posj) in enumerate(atomlist):
                # iterate over orbitals on center j
                for (nj,lj,mj) in valorbs[Zj]:
                    # enlarge gamma matrix
                    g0_AO[mu,nu] = gamma[i,j]
                    g1_AO[:,mu,nu] = grad_gamma[:,i,j]
                    
                    nu += 1
            mu += 1

    return g0_AO, g1_AO

    
################### including partial dipoles ################

def gaussian_decay_multipoles(hubbard_U):
    """
    Returns:
    ========
    sigma: dictionary with decay constants for each atomic number
    """
    sigmas = {}
    for Z,U in hubbard_U.iteritems():

        # monopoles only
        sig = 1.0/(np.sqrt(np.pi)*U)

        """
        # with dipoles
        a = 4.0+9.0*np.pi*U**2+3.0*np.sqrt(np.pi)*U*np.sqrt(8.0+9.0*np.pi*U**2)
        sig = 2.0 + 2.0*2.0**(2.0/3.0)/pow(a,1.0/3.0) + pow(2*a,1.0/3.0)
        sig /= 6*np.sqrt(np.pi)*U
        """
        sigmas[Z] = sig
    return sigmas

DIPOLE_SCALE = 1.0
    
class GammaGaussianMultipole(GammaFunction):
    """
    spherical charge fluctuations are modelled as Gaussians 

       FA(|r-RA|) = 1/(2 pi sigmaA^2)^(3/2) * exp( - (r-RA)^2/(2*sigmaA^2) )
       
    """
    def __init__(self, hubbard_U):
        sigmas_s = {}
        sigmas_p = {}
        for Z,U in hubbard_U.iteritems():
            sigmas_s[Z] = 1.0/(np.sqrt(np.pi)*U)
            sigmas_p[Z] = sigmas_s[Z] #sigmas_s[Z] #1.0/(6.0*np.sqrt(np.pi)*U**3) #
        """construct the C_AB matrix"""
        self.C_s = {}
        self.C_p = {}
        Zs = sigmas_s.keys()
        for ZA in Zs:
            for ZB in Zs:
                self.C_s[(ZA,ZB)] = 1.0/np.sqrt(2*(sigmas_s[ZA]**2 + sigmas_s[ZB]**2))
                self.C_p[(ZA,ZB)] = 1.0/np.sqrt(2*(sigmas_p[ZA]**2 + sigmas_p[ZB]**2))
    def g(self, R, rAB, ZA, ZB):
        """
        Parameters:
        ===========
        rAB: 3d unit vector pointing from atom A to atom B
        R: distance between atoms A and B
        """
        assert R > 0.0
        x = self.C_s[(ZA,ZB)]*R
        erfx = erf(x)
        xexpx2 = 2.0/np.sqrt(np.pi)* x*np.exp(-x**2)
        # monopole-monopole 
        gss = erfx/R
        # monopole-dipole
        gsp = -rAB/R**2 * (erfx - xexpx2)
        # dipole-monopole
        gps = -gsp
        # dipole-dipole
        x = self.C_p[(ZA,ZB)]*R
        erfx = erf(x)
        xexpx2 = 2.0/np.sqrt(np.pi)* x*np.exp(-x**2)

        gpp = (np.identity(3) * (erfx - xexpx2) \
               -np.outer(rAB,rAB) * (3*erfx - (3+2*x**2)*xexpx2)) / R**3

###### 
        gsp *= DIPOLE_SCALE
        gps *= DIPOLE_SCALE
        gpp *= DIPOLE_SCALE
######
        
        gAB = np.array([
                [gss,    gsp[0],   gsp[1],   gsp[2]  ],
                [gps[0], gpp[0,0], gpp[0,1], gpp[0,2]],
                [gps[1], gpp[1,0], gpp[1,1], gpp[1,2]],
                [gps[2], gpp[2,0], gpp[2,1], gpp[2,2]]])
        return gAB
    def g_limit0(self, Z):
        C_s = self.C_s[(Z,Z)]
        C_p = self.C_p[(Z,Z)]

        gss = 2.0/np.sqrt(np.pi) * C_s
        gsp = np.zeros(3)
        gps = np.zeros(3)
        gpp = 4.0/(3.0*np.sqrt(np.pi)) * C_p**3 * np.identity(3)


###### 
        gsp *= DIPOLE_SCALE
        gps *= DIPOLE_SCALE
        gpp *= DIPOLE_SCALE
######

        g0 = np.array([
                [gss,    gsp[0],   gsp[1],   gsp[2]  ],
                [gps[0], gpp[0,0], gpp[0,1], gpp[0,2]],
                [gps[1], gpp[1,0], gpp[1,1], gpp[1,2]],
                [gps[2], gpp[2,0], gpp[2,1], gpp[2,2]]])

        return g0


class GammaMatrixMultipoles:
    def __init__(self, gamma_func):
        self.gamma_func = gamma_func
    def gamma_atomwise(self, atomlist, distances, directions):
        """build the gamma-matrix and its gradients"""
        Nat = len(atomlist)
        g0 = np.zeros((4*Nat,4*Nat))  # gamma
        for i,(Zi,posi) in enumerate(atomlist):
            for j,(Zj,posj) in enumerate(atomlist):
                if i == j:
                    # 
                    g0[4*i:4*(i+1),4*j:4*(j+1)] = self.gamma_func.g_limit0(Zi)
                    #
                    """
                    g0_analim = self.gamma_func.g_limit0(Zi)
                    g0_numlim = self.gamma_func.g(1.0e-10, np.zeros(3), Zi,Zi)
                    err = np.sum(abs(g0_analim-g0_numlim))
                    assert err < 1.0e-8, "err = %s\ng0_analim = %s\n g0_numlim = %s" % (err, g0_analim, g0_numlim)
                    """
                else:
                    Rij = distances[i,j]
                    rij = directions[i,j,:]
                    g0[4*i:4*(i+1),4*j:4*(j+1)] = self.gamma_func.g(Rij, rij, Zi, Zj)
        return g0, None

############ Auxiliary Basis ####

class AuxiliaryBasisFunction:
    def __init__(self, Zi, center, sigma_s, sigma_p):
        self.Zi = Zi
        self.center = center
        self.sigma_s = sigma_s
        self.sigma_p = sigma_p
    def drho_dq(self, x,y,z):
        x0,y0,z0 = self.center
        dx,dy,dz = x-x0,y-y0,z-z0
        r2 = dx*dx+dy*dy+dz*dz
        Ns = 1.0/(2.0*np.pi*self.sigma_s**2)**(1.5)
        Np = 1.0/((2.0*np.pi*self.sigma_p**2)**(1.5) * self.sigma_p**2)
        Fs = Ns * np.exp(-r2/(2.0*self.sigma_s**2))
        Fp = Np * np.exp(-r2/(2.0*self.sigma_p**2))
        drho_dq_A = np.array([Fs, dx*Fp, dy*Fp, dz*Fp])
        # (dqA,ddipA_x,ddipA_y,ddipA_z) . (drho_dqA) = drho_A
        return drho_dq_A
    def partial_density(self, dqA, ddipA, x,y,z):
        """
        computes the distribution of the partial density around the atom
        A, which is modelled by a monopole- and a dipole-like Gaussian.
        """
        x0,y0,z0 = self.center
        dx,dy,dz = x-x0,y-y0,z-z0
        r2 = dx*dx+dy*dy+dz*dz
        Ns = 1.0/(2.0*np.pi*self.sigma_s**2)**(1.5)
        Np = 1.0/((2.0*np.pi*self.sigma_p**2)**(1.5) * self.sigma_p**2)
        Fs = Ns * np.exp(-r2/(2.0*self.sigma_s**2))
        Fp = Np * np.exp(-r2/(2.0*self.sigma_p**2))
        drho = dqA * Fs + (ddipA[0]*dx + ddipA[1]*dy + ddipA[2]*dz) * Fp
        return drho
        
class AuxiliaryBasisSet:
    def __init__(self, atomlist, hubbard_U):
        sigmas_s = {}
        sigmas_p = {}
        for Z,U in hubbard_U.iteritems():
            sigmas_s[Z] = 1.0/(np.sqrt(np.pi)*U)
            sigmas_p[Z] = sigmas_s[Z] #1.0/(6.0*np.sqrt(np.pi)*U**3)
        #
        self.bfs = []
        for i,(Zi, posi) in enumerate(atomlist):
            basis_function = AuxiliaryBasisFunction(Zi, posi, sigmas_s[Zi], sigmas_p[Zi])
            self.bfs.append( basis_function )
    def partial_density(self, dq, ddip, x,y,z):
        dense = np.zeros(x.shape)
        for i,bf in enumerate(self.bfs):
            dense += bf.partial_density(dq[i], ddip[i,:], x,y,z)
        return dense
    def drho_dq(self, x,y,z):
        """
        Find the derivative of the partial density with respect to 
        the partial charges and partial dipoles
        """
        deriv = np.zeros((4*len(self.bfs) , x.shape))
        for i,bf in enumerate(self.bfs):
            deriv[4*i:4*(i+1),:] = bf.drho_dq(x,y,z)
        return deriv
    

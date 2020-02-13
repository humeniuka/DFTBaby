"""
implementation of LDA XC-functionals
"""
from numpy import sqrt, exp, log, pi, linspace, arctan

class XC_None:
    def __init__(self):
        pass
    def setGrid(self, x):
        # GGA potentials need to know the grid in order
        # to differentiate
        self.x = x    
    def vxc(self,n):
        return 0.0*n

#### LOCAL DENSITY APPROXIMATION ##############################

class XC_PW92(XC_None):
    def __init__(self):
        """ The Perdew-Wang 1992 LDA exchange-correlation functional. """
        self.small=1E-90
        self.a1 = 0.21370
        self.c0 = 0.031091
        self.c1 = 0.046644
        self.b1 = 1.0/2.0/self.c0*exp(-self.c1/2.0/self.c0)
        self.b2 = 2*self.c0*self.b1**2
        self.b3 = 1.6382
        self.b4 = 0.49294
    def exc(self,n,der=0):
        """ Exchange-correlation with electron density n. """
        return self.e_x(n,der=der)+self.e_corr(n,der=der)

    def e_x(self,n,der=0):
        """ Exchange. """
        if der==0:
            return -3.0/4*(3*n/pi)**(1.0/3)
        elif der==1:
            return -3.0/(4*pi)*(3*n/pi)**(-2.0/3)

    def e_corr(self,n,der=0):
        """ Correlation energy. """
        rs = (3.0/(4*pi*n))**(1.0/3)
        aux=2*self.c0*( self.b1*sqrt(rs)+self.b2*rs+self.b3*rs**(3.0/2)+self.b4*rs**2 )
        if der==0:
            return -2*self.c0*(1+self.a1*rs)*log(1+aux**-1)
        elif der==1:
            return ( -2*self.c0*self.a1*log(1+aux**-1) \
                   -2*self.c0*(1+self.a1*rs)*(1+aux**-1)**-1*(-aux**-2)\
                   *2*self.c0*(self.b1/(2*sqrt(rs))+self.b2+3*self.b3*sqrt(rs)/2+2*self.b4*rs) )*( -(4*pi*n**2*rs**2)**-1 )

    def vxc(self,n):
        """ Exchange-correlation potential (functional derivative of exc). """
        return self.exc(n)+n*self.exc(n,der=1)

class XC_VWN(XC_None): # not sure if this gives correct xc-potential
    """
    Vosko, Wilk and Nusair XC functional as used by NIST's "Atomic Reference Data
    for Electronic Structure Calculations" for spin unpolarized densities

    see http://physics.nist.gov/PhysRefData/DFTdata/chap2.2.html#exchange
    """
    def __init__(self):
        pass
    def __F(self, rs, A, x0, b, c):
        def X(x):
            return x**2+b*x+c
        x = sqrt(rs)
        Q = sqrt(4*c-b**2)
        aux = arctan(Q/(2*x+b))
        F = A*(   log(x**2/X(x))             \
                      + 2*b/Q * aux          \
                - b*x0/X(x0) * ( log( pow(x-x0,2)/X(x) ) + 2*(b+2*x0)/Q * aux) )
        return F
    def __dF(self, rs, A, x0, b, c):
        """
        derivative of F with respect to rs
        """
        def X(x):
            return x**2+b*x+c
        def dX(x):
            return 2*x+b
        x = sqrt(rs)
        Q = sqrt(4*c-b**2)
        dx_drs = 1.0/(2.0*x)
        dF_dx = A*(   2.0/x - dX(x)/X(x)            \
                    - 4.0*b/( pow(2*x+b,2) + Q**2 ) \
                    - b*x0/X(x0)* ( 2.0/(x-x0) - dX(x)/X(x) - 4.0*(b+2.0*x0)/(pow(2.0*x+b,2) + Q**2) ))
        return dF_dx * dx_drs
    def exc(self,n,der=0):
        """ Exchange-correlation with electron density n. """
        return self.e_x(n,der=der)+self.e_corr(n,der=der)

    def e_x(self,n,der=0):
        """ Exchange. """
        if der==0:
            return -3.0/4*(3*n/pi)**(1.0/3) # -3.0/2.0*pow(3*n/pi, 1.0/3) # NIST larger by factor of 2 !?
        elif der==1:
            return -3.0/(4*pi)*(3*n/pi)**(-2.0/3) # -3.0/(2*pi)*pow(3*n/pi, -2.0/3)

    def e_corr(self,n,der=0):
        """ Correlation energy. """
        # paramagnetic
        A_p, x0_p, b_p, c_p = 0.0310907, -0.10498, 3.72744, 12.9352
        # electron gas parameter
        rs = pow(3.0/(4*pi*n), 1.0/3)
        
        if der==0:
            e_c = self.__F(rs, A_p, x0_p, b_p, c_p)
            return e_c
        elif der==1:
            dF_drs = self.__dF(rs, A_p, x0_p, b_p, c_p)
            drs_dn = -1.0/3.0 * pow(3.0/(4.0*pi), 1.0/3.0) * pow(n, -4.0/3.0)
            dec_dn = dF_drs * drs_dn
            return dec_dn 

    def vxc(self,n):
        """ Exchange-correlation potential (functional derivative of exc). """
        return self.exc(n)+n*self.exc(n,der=1)

########### LIB-XC #####################

from pylibxc.pylibxc import libXCFunctional

if __name__ == "__main__":
    from matplotlib.pyplot import plot, show

    rs = linspace(0.1, 5.0, 2000)
    rho = 1.0/ ( 4.0/3.0*pi*pow(rs,3) )
    xcpot_pw92 = XC_PW92()
    plot(rho, xcpot_pw92.vxc(rho))

    xcpot_vwn = XC_VWN()
    plot(rho, xcpot_vwn.vxc(rho), "o")
    
    show()

# -*- coding: utf-8 -*-
"""
routines for integrals between uncontracted cartesian Gaussian basis function
"""

import numpy as np

def basis_overlap(basis):
    """
    compute the overlap matrix between uncontracted basis functions

    Parameters
    ----------
    basis  : instance of UncontractedBasisSet object

    Returns
    -------
    S       : overlap matrix, S[i,j] = <i|j> overlap of the i-th basis function
              from `basis` with the j-th basis function from `basis`
    """
    S = np.zeros((basis.nbfs, basis.nbfs))
    for i in range(0, basis.nbfs):
        alpha1 = basis.exponents[i]
        (l1,m1,n1) = basis.powers[:,i]
        A = basis.centers[:,i]
        for j in range(i, basis.nbfs):
            alpha2 = basis.exponents[j]
            (l2,m2,n2) = basis.powers[:,j]
            B = basis.centers[:,j]
            S[i,j] = overlap(alpha1,(l1,m1,n1),A,alpha2,(l2,m2,n2),B)
            # overlap matrix is symmetric
            S[j,i] = S[i,j]
    return S

def basis_dipoles(basis):
    """
    compute matrix elements of dipole operator between uncontracted basis 
    functions

    Parameters
    ----------
    basis  : instance of UncontractedBasisSet object

    Returns
    -------
    D       : dipole matrix, D[:,i,j] = <i|r|j> dipole matrix element of the 
              i-th basis function from `basis` and the j-th basis function from `basis`
    """
    D = np.zeros((3, basis.nbfs, basis.nbfs), dtype=float)
    for i in range(0, basis.nbfs):
        alpha1 = basis.exponents[i]
        (l1,m1,n1) = basis.powers[:,i]
        A = basis.centers[:,i]
        for j in range(i, basis.nbfs):
            alpha2 = basis.exponents[j]
            (l2,m2,n2) = basis.powers[:,j]
            B = basis.centers[:,j]
            D[:,i,j] = dipole(alpha1,(l1,m1,n1),A,alpha2,(l2,m2,n2),B)
            # dipole matrix is symmetric under exchange of bra and ket
            D[:,j,i] = D[:,i,j]
    return D
    
def basis_angmomL(basis):
    """
    compute matrix elements of angular momentum operator L between uncontracted basis 
    functions

    Parameters
    ----------
    basis  : instance of UncontractedBasisSet object

    Returns
    -------
    L       : angular momentum matrix, L[:,i,j] = <i|L|j> angular momentum matrix element of the 
              i-th basis function from `basis` and the j-th basis function from `basis`
    """
    L = np.zeros((3, basis.nbfs, basis.nbfs), dtype=complex)
    for i in range(0, basis.nbfs):
        alpha1 = basis.exponents[i]
        (l1,m1,n1) = basis.powers[:,i]
        A = basis.centers[:,i]
        for j in range(i, basis.nbfs):
            alpha2 = basis.exponents[j]
            (l2,m2,n2) = basis.powers[:,j]
            B = basis.centers[:,j]
            L[:,i,j] = angmomL(alpha1,(l1,m1,n1),A,alpha2,(l2,m2,n2),B)
            # angular momentum matrix is Hermitian under exchange of bra and ket
            L[:,j,i] = L[:,i,j].conjugate()
    return L

def basis_angmomL2(basis):
    """
    compute matrix elements of total angular momentum operator L^2 = Lx^2 + Ly^2 + Lz^2
    between uncontracted basis functions

    Parameters
    ----------
    basis  : instance of UncontractedBasisSet object

    Returns
    -------
    L2     : total angular momentum, L2[i,j] = <i|L^2|j> total angular momentum matrix element 
             of the i-th basis function from `basis` and the j-th basis function from `basis`
    """
    L2 = np.zeros((basis.nbfs, basis.nbfs), dtype=complex)
    for i in range(0, basis.nbfs):
        alpha1 = basis.exponents[i]
        (l1,m1,n1) = basis.powers[:,i]
        A = basis.centers[:,i]
        for j in range(i, basis.nbfs):
            alpha2 = basis.exponents[j]
            (l2,m2,n2) = basis.powers[:,j]
            B = basis.centers[:,j]
            L2[i,j] = angmomL2(alpha1,(l1,m1,n1),A,alpha2,(l2,m2,n2),B)
            # total angular momentum matrix is Hermitian under exchange of bra and ket
            L2[j,i] = L2[i,j].conjugate()
    return L2


#########################################################################################
#
# INTEGRALS
#
# routines for computing overlaps and dipole matrix elements between primitive Gaussian
# basis functions, stolen from 'pyints.py' in PyQuante
#
# see THO paper
#   H.Taketa, S. Huzinga, K. O-ohata,
#  'Gaussian-Expansion Methods for Molecular Integrals'
#   J. Phys. Soc. Japan, 21, 2313, 1966.
# 
##########################################################################################

def fact(i):
    "Normal factorial"
    val = 1
    while (i>1):
        val = i*val
        i = i-1
    return val

def fact2(i):
    "Double factorial (!!) function = 1*3*5*...*i"
    val = 1
    while (i>0):
        val = i*val
        i = i-2
    return val

def binomial_prefactor(s,ia,ib,xpa,xpb):
    """From Augspurger and Dykstra"""
    sum = 0
    for t in xrange(s+1):
        if s-ia <= t <= ib:
            sum = sum + binomial(ia,s-t)*binomial(ib,t)* \
                  pow(xpa,ia-s+t)*pow(xpb,ib-t)
    return sum

def binomial(a,b):
    """Binomial coefficient"""
    return fact(a)/fact(b)/fact(a-b)

def norm(alpha, (l,m,n)):
    # eqn. (2.2) from THO paper
    nrm = np.sqrt(pow(2,2*(l+m+n)+1.5)*
                 pow(alpha,l+m+n+1.5)/
                 fact2(2*l-1)/fact2(2*m-1)/
                 fact2(2*n-1)/pow(np.pi,1.5))
    return nrm

def overlap_unnorm(alpha1,(l1,m1,n1),A,alpha2,(l2,m2,n2),B):
    """overlap without normalization constants. Taken from THO eq. 2.12"""
    rab2 = dist2(A,B)
    gamma = alpha1+alpha2
    P = gaussian_product_center(alpha1,A,alpha2,B)

    pre = pow(np.pi/gamma,1.5)*np.exp(-alpha1*alpha2*rab2/gamma)
    # intergrals over x,y and z can be performed independently
    wx = overlap_1D(l1,l2,P[0]-A[0],P[0]-B[0],gamma)
    wy = overlap_1D(m1,m2,P[1]-A[1],P[1]-B[1],gamma)
    wz = overlap_1D(n1,n2,P[2]-A[2],P[2]-B[2],gamma)
    return pre*wx*wy*wz

def overlap(alpha1,(l1,m1,n1),A,alpha2,(l2,m2,n2),B):
    """
    overlap <g1|g2> between two normalized cartesian Gaussian functions which
    are centered at A and B

                                  l1        m1       n1   -alpha1 (r-A)^2
    g1(x,y,z) = N(l1,m1,n1) (x-Ax)   (y-Ay)    (z-Az)    e
    
    and a similar expression for g2(x,y,z).
    
    Parameters
    ----------
    alpha1     :  exponent, > 0
    (l1,m1,n1) :  powers, L=l1+m1+n1 is the angular momentum
    A          :  3-dim vector, center of basis function

    and similary for second basis function

    Returns
    -------
    scalar, <g1|g2>
    """
    olap = overlap_unnorm(alpha1,(l1,m1,n1),A,alpha2,(l2,m2,n2),B)
    norm1 = norm(alpha1,(l1,m1,n1))
    norm2 = norm(alpha2,(l2,m2,n2))
    return norm1*norm2*olap

def dist2(A,B):
    return pow(A[0]-B[0],2)+pow(A[1]-B[1],2)+pow(A[2]-B[2],2)

def gaussian_product_center(alpha1,A,alpha2,B):
    """see Szabo & Ostlund"""
    gamma = alpha1+alpha2
    return (alpha1*A[0]+alpha2*B[0])/gamma,\
           (alpha1*A[1]+alpha2*B[1])/gamma,\
           (alpha1*A[2]+alpha2*B[2])/gamma


def overlap_1D(l1,l2,PAx,PBx,gamma):
    """Taken from THO eq. 2.12"""
    sum = 0
    for i in xrange(1+int(np.floor(0.5*(l1+l2)))):
        sum = sum + binomial_prefactor(2*i,l1,l2,PAx,PBx)* \
              fact2(2*i-1)/pow(2*gamma,i)
    return sum

def dipole(alpha1,(l1,m1,n1),A1, alpha2,(l2,m2,n2),A2):
    """dipole matrix element between two Gaussian orbitals"""
    olap = overlap_unnorm(alpha1,(l1,  m1  ,n1  ),A1, alpha2,(l2,m2,n2),A2)
    Dx   = overlap_unnorm(alpha1,(l1+1,m1  ,n1  ),A1, alpha2,(l2,m2,n2),A2) + A1[0]*olap
    Dy   = overlap_unnorm(alpha1,(l1  ,m1+1,n1  ),A1, alpha2,(l2,m2,n2),A2) + A1[1]*olap
    Dz   = overlap_unnorm(alpha1,(l1  ,m1  ,n1+1),A1, alpha2,(l2,m2,n2),A2) + A1[2]*olap
    norm1 = norm(alpha1,(l1,m1,n1))
    norm2 = norm(alpha2,(l2,m2,n2))
    
    return norm1*norm2*np.array([Dx,Dy,Dz])

def angmomL(alphaA, nA, A, alphaB, nB, B):
    """
    matrix element of angular momentum (Lx,Ly,Lz) between two Gaussian orbitals
    """
    L = np.array([0,0,0], dtype=complex)
    # unit vectors
    e = [np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1])]

    for i,j,k in zip([1,2,0], [2,0,1], [0,1,2]):
        #         
        # L_k  = -I hbar (r_i d/dr_j - r_j d/dr_i)
        #
        if nB[i] > 0:
            L[k] += -nB[i]         * overlap_unnorm(alphaA, nA, A, alphaB, nB-e[i]+e[j], B) \
                    -nB[i]*B[j]    * overlap_unnorm(alphaA, nA, A, alphaB, nB-e[i]     , B)
        if nB[j] > 0:
            L[k] +=  nB[j]         * overlap_unnorm(alphaA, nA, A, alphaB, nB+e[i]-e[j], B) \
                    +nB[j]*B[i]    * overlap_unnorm(alphaA, nA, A, alphaB, nB     -e[j], B)
        L[k] +=     +2*alphaB*B[j] * overlap_unnorm(alphaA, nA, A, alphaB, nB+e[i]     , B) \
                    -2*alphaB*B[i] * overlap_unnorm(alphaA, nA, A, alphaB, nB     +e[j], B)

    normA = norm(alphaA,nA)
    normB = norm(alphaB,nB)

    L *= -1.0j * normA * normB
    
    return L

def angmomL2(alphaA, nA, A, alphaB, nB, B):
    """
    matrix element of total angular momentum L^2=Lx^2+Ly^2+Lz^2 between two Gaussian orbitals
    """
    L2 = 0.0+0.0j
    # unit vectors
    e = [np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1])]

    for i,j,k in zip([1,2,0], [2,0,1], [0,1,2]):
        #  2        
        # L  = hbar^2 (r_i d/dr_j - r_j d/dr_i)^2
        #  k

        # <L^2> is computed as <a|L^2|b> = (L|a>)^H (L|b>)
        # Applying L to the bra <a| gives <a|L = sum_i  coeffsA[i] <gaussian(alphaA,powersA[i],A)| 
        powersA = [nA+e[i]-e[j], nA-e[i]+e[j], nA-e[j]   , nA-e[i]    , nA+e[i]      , nA+e[j]       ]
        coeffsA = [nA[j]       , -nA[i]      , nA[j]*A[i], -nA[i]*A[j], 2*alphaA*A[j], -2*alphaA*A[i]]
        # Applying L to the ket |b> gives L|b> = sum_j  coeffsB[j] |gaussian(alphaB,powersB[j],B)>.
        powersB = [nB+e[i]-e[j], nB-e[i]+e[j], nB-e[j]   , nB-e[i]    , nB+e[i]      , nB+e[j]       ]
        coeffsB = [nB[j]       , -nB[i]      , nB[j]*B[i], -nB[i]*B[j], 2*alphaB*B[j], -2*alphaB*B[i]]

        # <L^2> can now be calculated as the overlap between <a|L and L|b>.
        for cA,pA in zip(coeffsA, powersA):
            for cB,pB in zip(coeffsB, powersB):
                c = cA*cB
                if abs(c) > 0.0:
                    assert np.all(pA >= 0)
                    assert np.all(pB >= 0)
                    L2 += c * overlap_unnorm(alphaA, pA, A, alphaB, pB, B)

    # normalization constants
    normA = norm(alphaA,nA)
    normB = norm(alphaB,nB)

    L2 *= normA * normB

    assert abs(L2.imag) < 1.0e-10
    
    return L2.real

############################################################
#
# Electron Repulsion Integrals (ERI)
# stolen from R. Muller's PyQuante (pyints.py)
#
############################################################

def coulomb_repulsion((xa,ya,za),norma,(la,ma,na),alphaa,
                      (xb,yb,zb),normb,(lb,mb,nb),alphab,
                      (xc,yc,zc),normc,(lc,mc,nc),alphac,
                      (xd,yd,zd),normd,(ld,md,nd),alphad):

    rab2 = dist2((xa,ya,za),(xb,yb,zb))
    rcd2 = dist2((xc,yc,zc),(xd,yd,zd))
    xp,yp,zp = gaussian_product_center(alphaa,(xa,ya,za),alphab,(xb,yb,zb))
    xq,yq,zq = gaussian_product_center(alphac,(xc,yc,zc),alphad,(xd,yd,zd))
    rpq2 = dist2((xp,yp,zp),(xq,yq,zq))
    gamma1 = alphaa+alphab
    gamma2 = alphac+alphad
    delta = 0.25*(1/gamma1+1/gamma2)

    Bx = B_array(la,lb,lc,ld,xp,xa,xb,xq,xc,xd,gamma1,gamma2,delta)
    By = B_array(ma,mb,mc,md,yp,ya,yb,yq,yc,yd,gamma1,gamma2,delta)
    Bz = B_array(na,nb,nc,nd,zp,za,zb,zq,zc,zd,gamma1,gamma2,delta)

    sum = 0.
    for I in xrange(la+lb+lc+ld+1):
        for J in xrange(ma+mb+mc+md+1):
            for K in xrange(na+nb+nc+nd+1):
                sum = sum + Bx[I]*By[J]*Bz[K]*Fgamma(I+J+K,0.25*rpq2/delta)

    return 2*pow(np.pi,2.5)/(gamma1*gamma2*np.sqrt(gamma1+gamma2)) \
           *np.exp(-alphaa*alphab*rab2/gamma1) \
           *np.exp(-alphac*alphad*rcd2/gamma2)*sum*norma*normb*normc*normd

def B_term(i1,i2,r1,r2,u,l1,l2,l3,l4,Px,Ax,Bx,Qx,Cx,Dx,gamma1,gamma2,delta):
    "THO eq. 2.22"
    return fB(i1,l1,l2,Px,Ax,Bx,r1,gamma1) \
           *pow(-1,i2)*fB(i2,l3,l4,Qx,Cx,Dx,r2,gamma2) \
           *pow(-1,u)*fact_ratio2(i1+i2-2*(r1+r2),u) \
           *pow(Qx-Px,i1+i2-2*(r1+r2)-2*u) \
           /pow(delta,i1+i2-2*(r1+r2)-u)

def B_array(l1,l2,l3,l4,p,a,b,q,c,d,g1,g2,delta):
    Imax = l1+l2+l3+l4+1
    B = [0]*Imax
    for i1 in xrange(l1+l2+1):
        for i2 in xrange(l3+l4+1):
            for r1 in xrange(i1/2+1):
                for r2 in xrange(i2/2+1):
                    for u in xrange((i1+i2)/2-r1-r2+1):
                        I = i1+i2-2*(r1+r2)-u
                        B[I] = B[I] + B_term(i1,i2,r1,r2,u,l1,l2,l3,l4,
                                             p,a,b,q,c,d,g1,g2,delta)
    return B

def fB(i,l1,l2,P,A,B,r,g):
    return binomial_prefactor(i,l1,l2,P-A,P-B)*B0(i,r,g)

def B0(i,r,g): return fact_ratio2(i,r)*pow(4*g,r-i)

def fact_ratio2(a,b): return fact(a)/fact(b)/fact(a-2*b)

def Fgamma(m,x):
    "Incomplete gamma function"
    SMALL=0.00000001
    x = max(abs(x),SMALL)
    val = gamm_inc(m+0.5,x)
    return 0.5*pow(x,-m-0.5)*val;

def gammln(x):
    "Numerical recipes, section 6.1"
    cof = [76.18009172947146,-86.50532032941677,
           24.01409824083091,-1.231739572450155,
           0.1208650973866179e-2,-0.5395239384953e-5]
    y=x
    tmp=x+5.5
    tmp = tmp - (x+0.5)*np.log(tmp)
    ser=1.000000000190015 # don't you just love these numbers?!
    for j in xrange(6):
        y = y+1
        ser = ser+cof[j]/y
    return -tmp+np.log(2.5066282746310005*ser/x);

def gamm_inc(a,x):
    "Incomple gamma function \gamma; computed from NumRec routine gammp."
    gammap,gln = gammp(a,x)
    return np.exp(gln)*gammap
    
def gammp(a,x):
    "Returns the incomplete gamma function P(a;x). NumRec sect 6.2."
    assert (x > 0 and a >= 0), "Invalid arguments in routine gammp"

    if x < (a+1.0): #Use the series representation
        gamser,gln = _gser(a,x)
        return gamser,gln
    #Use the continued fraction representation
    gammcf,gln = _gcf(a,x)
    return 1.0-gammcf ,gln  #and take its complement.

def _gser(a,x):
    "Series representation of Gamma. NumRec sect 6.1."
    ITMAX=100
    EPS=3.e-7

    gln=gammln(a)
    assert(x>=0),'x < 0 in gser'
    if x == 0 : return 0,gln

    ap = a
    delt = sum = 1./a
    for i in xrange(ITMAX):
        ap=ap+1.
        delt=delt*x/ap
        sum=sum+delt
        if abs(delt) < abs(sum)*EPS: break
    else:
        print 'a too large, ITMAX too small in gser'
    gamser=sum*np.exp(-x+a*np.log(x)-gln)
    return gamser,gln

def _gcf(a,x):
    "Continued fraction representation of Gamma. NumRec sect 6.1"
    ITMAX=100
    EPS=3.e-7
    FPMIN=1.e-30

    gln=gammln(a)
    b=x+1.-a
    c=1./FPMIN
    d=1./b
    h=d
    for i in xrange(1,ITMAX+1):
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if abs(d) < FPMIN: d=FPMIN
        c=b+an/c
        if abs(c) < FPMIN: c=FPMIN
        d=1./d
        delt=d*c
        h=h*delt
        if abs(delt-1.) < EPS: break
    else:
        print 'a too large, ITMAX too small in gcf'
    gammcf=np.exp(-x+a*np.log(x)-gln)*h
    return gammcf,gln

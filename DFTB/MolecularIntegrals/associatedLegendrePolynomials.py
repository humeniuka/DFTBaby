from numpy import *
import matplotlib.pyplot as plt
from scipy.special import lpmn
from mpmath import mpf # needed to evalulate ratios of large factorials
# adapted from Warren Weckesser for use with numpy
# http://projects.scipy.org/scipy/attachment/ticket/1296/assoc_legendre.py

def fact(n):
    """n!"""
    f = 1
    while(n>1):
        f *= n
        n -= 1
    return(f)
    
def fact2(k):
    """n!!"""
    if k <= 0:
        return 1
    f = k
    while k >= 3:
        k -= 2
        f *= k
    return f

def assoc_legendre_p(l, m, x):
    """Compute the associated Legendre polynomial with degree l and order m.

    This function uses the recursion formula in the degree l.
    (Abramowitz & Stegun, Section 8.5.) 
    """
    if m < 0 or m > l:
        raise ValueError('require 0 <= m <= l, but m=%d and l=%d' % (m, l))
    if max(x) > 1.0 or min(x) < -1.0:
        raise ValueError('require -1 <= x <= 1, but x=%f', x)

    # Compute the initial term in the recursion.
    if m == 0:
        one = x/x # array of 1.0's
        p_l = one
    else:
        s = 1
        if m & 1:
            s = -1
        z = sqrt(1.0 - pow(x,2))
        p_l = s * fact2(2*m - 1) * pow(z, m)

    if m == l:
        return p_l
    # Compute the next term in the recursion.
    p_l_prev = p_l
    p_l = x * (2*m + 1) * p_l
    if l == m + 1:
        return p_l
    # Iterate the recursion relation.
    for n in range(m+2, l+1):
        result = (x*(2*n - 1)*p_l - (n + m - 1)*p_l_prev)/(n - m)
        p_l_prev = p_l
        p_l = result

    return result

class CacheFactorial:
    """
    The instances of this class can be used in the same way as a function returning
    the factorial of a number. However, already calculated factorials are cached.
       factorial = CacheFactorial()
       print factorial(10)
       print factorial(1)
       print factorial(10)       # looked up in cache
    The overhead is not really worth the speedup!!!
    """
    def __init__(self, cache_size=100):
        self.cache = [None for n in range(0, cache_size)]
    def __call__(self, n):
        fn = self.cache[n]
        if (fn == None):
            """n! not in cache, so compute it"""
            fn = fact(n)
            self.cache[n] = fn
        return(fn)

class CacheFactorial2:
    def __init__(self, cache_size=100):
        self.cache = [None for n in range(0, cache_size)]
    def __call__(self, n):
        f2n = self.cache[n]
        if (f2n == None):
            f2n = fact2(n)
            self.cache[n] = f2n
        return(f2n)

# implementation with iterators, better use the renormalized version below which
# avoids overflows for high degrees

def increase_degree(P_lminus1_m, P_lminus2_m, x, l, m):
    """ 
    see Abramowitz/Stegun formula 8.5.3 for varying degree 

    P_(m,l)(x) = ((2*l-1) * x * P_(m,l-1)(x) - (l-1+m) * P_(m,l-2)(x))/(l-m)      formula (C)
    
    graphically the effect of formula C is depicted by
        P_(m,l-1) --(C)--> P_(m,l)
    C increases the degree but leaves the order.
    """
    return( ((2*l-1)*x*P_lminus1_m - (l-1+m)*P_lminus2_m)/float(l-m) )

def initial_term(x, l):
    """
    compute associated Legendre polynomial for l=m
    see http://en.wikipedia.org/wiki/Associated_Legendre_polynomials under 'Helpful Identities'

    P_(l,l)(x) = (-1)^l (2l-1)!! (1-x^2)^(l/2)          formula (R)

    Formula R creates a new starting point for an iteration over the degree
         --(R)--> P_(l,l)
    """
    pll = pow(-1,l)*fact2(2*l-1)*pow( sqrt(1.0-pow(x,2)), l)
    return pll
    

def assocLegendre_it(x):
    """create an iterator for the associate Legendre polynomials
    The iteration scheme can be visualized as a matrix where for each call to R a new row is
    inserted which is then filled from left to right by successive calls to C.
    --(R)--> P_(0,0) --(C)--> P_(0,1) --(C)--> P_(0,2) --(C)--> P_(0,3) --(C)--> ....
                     --(R)--> P_(1,1) --(C)--> P_(1,2) --(C)--> P_(1,3) --(C)--> ....
                                      --(R)--> P_(2,2) --(C)--> P_(2,3) --(C)--> ....
                                                          ....   .....        
    The iterator returns the polynomials in the following order:
      P_(0,0),
      P_(1,0),P_(1,1),
      P_(2,0),P_(2,1),P_(2,2),
      ...
      P_(l,0),P_(l,1),...,P_(l,m),...,P_(l,l)
      ...
    """
    if getattr(x, '__iter__', False): # got array
        xf = x.ravel()
        if max(xf) > 1.0 or min(xf) < -1.0:
            raise ValueError('require -1 <= x <= 1, but x=%f', xf)
    else:  # got number
        if x > 1.0 or x < -1.0:
            raise ValueError('require -1 <= x <= 1, but x=%f', x)
    p00 = ones(x.shape)
    """P_(0,0) is just 1.0"""
    yield p00
    lastPcol = [zeros(x.shape)]
    curPcol = [p00]
    l = 1
    while True:
        nextPcol = range(0, l)
        for m in xrange(0, l):
            nextPcol[m] = increase_degree(curPcol[m], lastPcol[m], x, l, m)
            yield nextPcol[m]
        curPcol.append(zeros(x.shape))
        pll = initial_term(x,l)
        yield pll
        nextPcol.append(pll)
        """add P_(l,l) to current column and prepare for the next iteration 
        which will generate P_(l+1,0), P_(l+1,1),...,P_(l+1,l+1)"""
        l += 1
        lastPcol = curPcol
        curPcol = nextPcol

# renormalized associated Legendre polynomials, 
# (2*l-1)!! Pr_(m,l) = P_(m,l)
# The double factorial can be cancelled 

def increase_degree_renorm(P_lminus1_m, P_lminus2_m, x, l, m):
    """ 
    Define renormalized associated Legendre polynmials by (2*l-1)!! * Pr_(m,l) = P_(m,l)
    This avoids overflows because of the double factorial

    Pr_(m,l)(x) = (x * Pr_(m,l-1)(x) - (l-1+m)/((2*l-1)*(2*l-3)) * Pr_(m,l-2)(x))/(l-m)      formula (C)
    
    graphically the effect of formula C is depicted by
        P_(m,l-1) --(C)--> P_(m,l)
    C increases the degree but leaves the order.
    """
    return( (x*P_lminus1_m - (l-1+m)/((2*l-1.0)*(2*l-3.0))*P_lminus2_m)/float(l-m) )

def initial_term_renorm(x, l):
    """
    P_(l,l)(x) = (-1)^l (1-x^2)^(l/2)          formula (R)

    Formula R creates a new starting point for an iteration over the degree
         --(R)--> P_(l,l)
    """
    pll = pow(-1,l)*pow( sqrt(1.0-pow(x,2)), l)
    return pll
    

def assocLegendre_renorm_it(x):
    """create an iterator for the renormalized associate Legendre polynomials
    The iteration scheme can be visualized as a matrix where for each call to R a new row is
    inserted which is then filled from left to right by successive calls to C.
    --(R)--> P_(0,0) --(C)--> P_(0,1) --(C)--> P_(0,2) --(C)--> P_(0,3) --(C)--> ....
                     --(R)--> P_(1,1) --(C)--> P_(1,2) --(C)--> P_(1,3) --(C)--> ....
                                      --(R)--> P_(2,2) --(C)--> P_(2,3) --(C)--> ....
                                                          ....   .....        
    The iterator returns the polynomials in the following order:
      P_(0,0),
      P_(1,0),P_(1,1),
      P_(2,0),P_(2,1),P_(2,2),
      ...
      P_(l,0),P_(l,1),...,P_(l,m),...,P_(l,l)
      ...
    In order to obtain the unrenormalized polynomials the term P_(l,m) has to be 
    multiplied by (2*l-1)!!
    """
    if getattr(x, '__iter__', False): # got array
        xf = x.ravel()
        if max(xf) > 1.0 or min(xf) < -1.0:
            raise ValueError('require -1 <= x <= 1, but x=%f', xf)
    else:  # got number
        if x > 1.0 or x < -1.0:
            raise ValueError('require -1 <= x <= 1, but x=%f', x)
    p00 = ones(x.shape)
    """P_(0,0) is just 1.0"""
    yield p00
    lastPcol = [zeros(x.shape)]
    curPcol = [p00]
    l = 1
    while True:
        nextPcol = range(0, l)
        for m in xrange(0, l):
            nextPcol[m] = increase_degree_renorm(curPcol[m], lastPcol[m], x, l, m)
            yield nextPcol[m]
        curPcol.append(zeros(x.shape))
        pll = initial_term_renorm(x,l)
        yield pll
        nextPcol.append(pll)
        """add P_(l,l) to current column and prepare for the next iteration 
        which will generate P_(l+1,0), P_(l+1,1),...,P_(l+1,l+1)"""
        l += 1
        lastPcol = curPcol
        curPcol = nextPcol


# spherical harmonics

def spherical_harmonics_it(th, phi, outerproducts=True):
    """
    th is the polar coordinate of the points on the grid [0,pi]x[0,2pi]
    phi is the azimuthal coordinate of the points on the grid [0,pi]x[0,2pi]
    
    Returns an iterator to the values of the spherical harmonics
    on the grid [0,pi]x[0,2pi] in the following order
    Y_(0,0),
    Y_(1,0), Y_(1,+1), Y_(1,-1),
    Y_(2,0), Y_(2,+1), Y_(2,-1), Y_(2,+2), Y_(2,-2),
    .....
    i.e. not in the order -m, -(m-1),...,m-1,m  but  0,1,-1,2,-2,...,m-1,-(m-1),m,-m

    The iterator returns a tuple (Ylm, l, m)
    """
    x = cos(th)
    ap_it = assocLegendre_renorm_it(x)
    l = 0
    while True:
        for m in xrange(0, l+1):
            Plm = ap_it.next()
            # factorials = sqrt( fact(l-m)/float(fact(l+m)) ) * fact2(2*l-1)
            """use high precision mpf type to calculate ratios of factorials"""
            factorials = float( sqrt( mpf(fact(l-m))/mpf(fact(l+m)) ) * mpf(fact2(2*l-1)) )
            N = sqrt((2*l+1)/(4.0*pi)) * factorials
            #print "N(%s,%s) = %s" % (l,m,N)
            if outerproducts == True:
                Ylm = N * Plm * exp(1.0j*m*phi) 
            else:
                Ylm = N * outer(Plm, exp(1.0j*m*phi))
            yield (Ylm, l, m)
            if (m > 0):
                Yl_minus_m = pow(-1,m)*Ylm.conjugate()
                yield (Yl_minus_m, l, -m)
        l += 1

def sphvec(th,phi, nchannels):
    """return a vector     array([Y00,Y10,Y11,Y1-1,...])
    with with the spherical harmonics evaluated at the same angle
    the maximum angular momentum is related to the number of channels
    by: nchannels = lmax*(lmax+2)+1
    """
    Ylm_vec = []
    sph_it = spherical_harmonics_it(th,phi)
    for i in range(0, nchannels):
        Ylm,l,m = sph_it.next()
        Ylm_vec.append(Ylm)
    print "nchannels = %s" % nchannels
    assert m == -l 
    assert nchannels == l*(l+2)+1
    return array(Ylm_vec)


# Tests

def test_assocLegendre():
    x = linspace(-1.0, 1.0, 100)
    ap_it = assocLegendre_it(x)
    for l in xrange(0, 100):
        for m in xrange(0, l+1):
            print "l = %s, m = %s" % (l,m)
            Plm = assoc_legendre_p(l, m, x)
#            Plm_scipy = array([lpmn(m,l,xi)[0][-1][-1] for xi in x])
            Plm_it = ap_it.next()
            plt.plot(x, Plm)
#            plt.plot(x, Plm_scipy, ls="-.")
#            plt.plot(x, Plm_it, ls="-.")
    plt.savefig("assoc_legendre_p.png")

def test_assocLegendre_renorm():
    x = linspace(-1.0, 1.0, 100)
    ap_it = assocLegendre_renorm_it(x)
    for l in xrange(0, 100):
        for m in xrange(0, l+1):
            print "l = %s, m = %s" % (l,m)
            Plm = assoc_legendre_p(l, m, x)
            Plm_scipy = array([lpmn(m,l,xi)[0][-1][-1] for xi in x])
            Plm_it = fact2(2*l-1)*ap_it.next()
            dif = sum(abs(Plm_scipy - Plm_it))/sum(abs(Plm_scipy))
            print "dif = %s" % dif
            assert dif < 1e-10
#            plt.plot(x, Plm)
#            plt.plot(x, Plm_scipy, ls="-.")
#            plt.plot(x, Plm_it, ls="-.")
#    plt.savefig("assoc_legendre_renorm.png")


def test_factorialCache():
    factorial = CacheFactorial(cache_size=200)
    for n in xrange(0, 10):
        print factorial(n)
    for n in xrange(0, 10):
        print factorial(n)

def test_spherical_harmonics():
    from scipy.special import sph_harm
    th = linspace(0.000001, pi, 10)
    phi = linspace(0, 2.0*pi, 10, endpoint=False)
    Omega = outer(th, phi)
    TH = outer(th, ones(phi.shape))
    PHI = outer(ones(th.shape), phi)
    Ylm_it = spherical_harmonics_it(TH, PHI)
#    Ylm_it = spherical_harmonics_it(th, phi)
    for l in xrange(0, 50):
        for m in xrange(0, l+1):
            Ylm, l, m = Ylm_it.next()
            Ylm_scipy = sph_harm(m,l,outer(ones(th.shape), phi),outer(th, ones(phi.shape)))
            if (m > 0):
                Yl_minus_m, l, mminus = Ylm_it.next()
                Ylm_scipy_minus_m = sph_harm(-m,l,outer(ones(th.shape), phi),outer(th, ones(phi.shape)))
            difplus = Ylm - Ylm_scipy
            if (m > 0):
                difminus = Yl_minus_m - Ylm_scipy_minus_m
                deltaY = (sqrt(sum(abs(difplus)) + sum(abs(difminus))))/sqrt(sum(abs(Yl_minus_m)) + sum(abs(Ylm)))
            else:
                deltaY = sqrt(sum(abs(difplus)))/sqrt(sum(abs(Ylm)))
            print "|deltaY(%s,%s)| = %s" % (l,m,deltaY)
            assert deltaY < 1e-6
            """
            print Ylm_scipy.real.shape
            plt.imshow(Ylm_scipy.real)
            plt.savefig("/tmp/sph_scipy_%s_%s.png" % (l,m))
            plt.imshow(Ylm.real)
            plt.savefig("/tmp/sph_it_%s_%s.png" % (l,m))
            if (m > 0):
                plt.imshow(Ylm_scipy_minus_m.real)
                plt.savefig("/tmp/sph_scipy_%s_%s.png" % (l,-m))
                plt.imshow(Yl_minus_m.real)
                plt.savefig("/tmp/sph_it_%s_%s.png" % (l,-m))
            assert(deltaY < 10e-5)
            """
            # use the command
            #    animate -delay 100 /tmp/sph_it_*.png & animate -delay 100 /tmp/sph_scipy_*.png &
            # to compare the plots of spherical harmonics generated by scipy with this implementation

if __name__ == "__main__":
    test_assocLegendre_renorm()
    test_factorialCache()
    test_spherical_harmonics()

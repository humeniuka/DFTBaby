#!/usr/bin/env python
"""
Chebyshev interpolation polynomials
"""
# -*- coding: utf-8 -*-
import numpy as np

def chebyshev_it(x):
    """
    generate Chebyshev polynomials of the first kind T_n(x) through the recurrence relation

        T_0(x) = 1
        T_1(x) = x
        T_{n+1}(x) = 2 x T_n(x) - T_{n-1}(x)

    The polynomials are evaluated at the points x

    Parameters
    ----------
    x          : numpy array of shape (n,) with x-values, -1 <= x[i] <= 1

    Returns
    -------
    iterator to numpy arrays with value of Chebyshev polynomials
    T0(x), T1(x), ...
    """
    assert np.all((x >= -1.0) & (x <= 1.0)), "x-value have to be from the interval [-1,1], \n however  x= %s" % x
    # initiate recurrence relation
    T0 = 1.0 + 0*x
    yield T0
    T1 = x
    yield T1
    # T_{n-1}
    Tnm1 = T0
    # T_n
    Tn = T1
    while True:
        # T_{n+1}
        Tnp1 = 2.0*x*Tn - Tnm1
        yield Tnp1
        # shift n+1 -> n
        Tnm1 = Tn
        Tn = Tnp1

class ChebyshevInterpolator(object):
    def __init__(self, n):
        self.n = n
        # sampling points
        k = np.array(range(1, n+1))
        x = np.cos( (k-0.5)/n * np.pi )
        self.points = x
        # evaluate Chebyshev polynomials of degree < n
        # at the sampling points
        self.T = np.zeros((n,n))
        T_it = chebyshev_it(x)
        for i,Ti in enumerate(T_it):
            if i == 0:
                # The first term needs to be halfed
                self.T[i,:] = 0.5 * Ti
            else:
                self.T[i,:] = Ti
            if i == n-1:
                break
        # check that the sampling points are the roots of T_n
        Tn = T_it.next()
        assert np.all(abs(Tn) < 1.0e-10)

    def getSamplingPoints(self):
        """
        Chebyshev-Gauss sampling points

                     k-1/2
            x  = cos(----- pi)       for  k=1,...,n
             k         n

        """
        return self.points
            
    def expansion_coefficients(self, f):
        """
        discrete Chebyshev coefficients for representing the f(x) as a linear
        combination of Chebyshev polynomials
                       n
             f(x) = sum     a  T (x)
                       i=1   i  i

        Parameters
        ----------
        f          :  numpy array with values of function at the sampling points,
                      f[i] = f(x[i])

        Returns
        -------
        coeffs     :  numpy array with n Chebyshev coefficients a_i
        """
        # According to http://mathworld.wolfram.com/Chebyshev-GaussQuadrature.html the expansion
        # coefficients are obtained by Chebyshev Gauss quadrature for the integral
        #
        #            /+1      2  -1/2 
        #  a  = 2/pi | dx (1-x  )     f(x) T (x)
        #   i        /-1                    i
        #               n
        #     = 2/pi sum    w  f(x ) T (x )
        #               k=1  k    k   i  k
        #
        # where w_k = pi/n are the quadrature weights and x_k = cos([k-1/2]/n pi) are the sampling points
        # for k=1,...,n.
        #
        #          n
        #  a  = sum     2/n  T (x )  f(x )
        #   i      k=1        i  k      k
        #
        coeffs = 2.0/self.n * np.dot(self.T, f)
        return coeffs

    def interpolate(self, coeffs, x):
        """
        evaluate the interpolation polynomial given by the Chebyshev coefficients `coeffs`
        at the grid points x

        Parameters
        ----------
        coeffs     :  numpy array of shape (n,) with Chebyshev coefficients as obtained from `expansion_coefficients(...)`
        x          :  numpy array of arbitrary shape with grid points for evaluating the Chebyshev
                      expansion
       
        Returns
        -------
        f          :  numpy array with same shape as x, interpolated values f(x)
        """
        f = 0*x
        # 
        T_it = chebyshev_it(x)
        # evaluate Chebyshev expansion
        #
        #    (interp)            n-1
        #   f         (x)  =  sum     a  T (x)
        #                        i=0   i  i
        #
        # The first term in the expansion (i=0) needs to be halfed. This factor 1/2 has
        # been taken care of by replaing a0 with a0/2 in the computation of the coefficients.
        for i,Ti in enumerate(T_it):
            f += coeffs[i] * Ti
            if i == self.n-1:
                break
        return f
    
    
def plot_chebyshev_polynomials():
    x = np.linspace(-1.0, 1.0, 100)
    
    T_it = chebyshev_it(x)

    import matplotlib.pyplot as plt
    for n,Tn in enumerate(T_it):
        plt.plot(x, Tn, label=r"$T_{%d}$" % n)

        if n > 10:
            break
    plt.xlabel("x")
    plt.ylabel("Chebyshev polynomial")
    plt.legend()
    plt.show()

def abscissas():
    # abscissas in Becke's code
    n = 5
    k = np.array(range(1,n+1))
    # grid points on interval [-1,1]
    zr = k/(n+1.0)
    xr = np.cos(zr * np.pi)

    print xr

    # Chebyshev-Gauss points
    k = np.array(range(0,n+1))
    zr = (2.0*k+1.0)/(2.0*n+2.0)
    xr = np.cos(zr * np.pi)
    print xr

    #
    k = np.array(range(1, n+1))
    zr = (k-0.5)/n
    xr = np.cos(zr * np.pi)
    print xr

def test_chebyshev_interpolation():
    x = np.linspace(-1.0, 1.0, 100)
    def f(x):
        return np.sin(2*np.pi * x)
        #return np.exp(-x**2)

    import matplotlib.pyplot as plt
    plt.plot(x, f(x), lw=3, label=r"f(x)")

    # interpolation with different numbe of sampling points
    for n in range(5, 100, 10):
        cheb = ChebyshevInterpolator(n)

        x_sample = cheb.getSamplingPoints()
        f_sample = f(x_sample)
        coeffs = cheb.expansion_coefficients(f_sample)
        f_interp = cheb.interpolate(coeffs, x)

        # plot original function at the interpolation points,
        # does the interpolation polynomial pass through these points?
        l, = plt.plot(x_sample, f_sample, "o")
        
        # plot interpolated function 
        plt.plot(x, f_interp, ls="--", label=r"f(x) (interp. n = %d)" % n, color=l.get_color())

    plt.legend()
    plt.show()
    
if __name__ == "__main__":
    #plot_chebyshev_polynomials()
    #abscissas()
    test_chebyshev_interpolation()

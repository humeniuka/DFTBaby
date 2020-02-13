#!/usr/bin/env python
from numpy import *
from DFTB.MolecularIntegrals.lebedev_quad_points import LebedevGridPoints, Lebedev_Npts, Lebedev_Lmax, Lebedev_L2max
from DFTB.MolecularIntegrals.associatedLegendrePolynomials import spherical_harmonics_it
from itertools import izip

def outerN(a, b):
    """
    compute the outer product of two arrays a and b, which do not 
    have to be vectors but can have any shape.
    """
    ar = a.ravel()
    br = b.ravel()
    axbr = outer(ar,br)
    axb = axbr.reshape(a.shape + b.shape)
    return axb

def check_lebedev_quadrature(order, tolerance=1e-9):
    th,phi,w = array(LebedevGridPoints[Lebedev_Npts[order]]).transpose()
    sph_it1 = spherical_harmonics_it(th,phi)
    for Ylm1,l1,m1 in sph_it1:
        if l1 > Lebedev_L2max[order]:
            break
        sph_it2 = spherical_harmonics_it(th,phi)
        for Ylm2,l2,m2 in sph_it2:
            if l2 > Lebedev_L2max[order]:
                break
            I = 4.0*pi*sum(w*Ylm1*Ylm2.conjugate())
            print "<%s %s|%s %s> = %s" % (l1,m1,l2,m2,I)
            if l1 == l2 and m1 == m2:
                print "|I-1.0| = %s" % abs(I-1.0)
                assert abs(I-1.0) < tolerance
            else:
                assert abs(I) < tolerance

def test_lebedev_quadrature():
    for order in range(0, len(Lebedev_Npts)):
        print "Test Lebedev grid number %s with %s points" % (order, Lebedev_Npts[order])
        print "============================================"
        check_lebedev_quadrature(order)

def get_lebedev_grid(order):
    """find grid closest to requested order"""
    n = abs(array(Lebedev_L2max) - order).argmin()
    if order != Lebedev_L2max[n]:
        print "No grid for order %s, using grid which integrates up to L2max = %s exactly instead." \
            % (order, Lebedev_L2max[n])
    th,phi,w = array(LebedevGridPoints[Lebedev_Npts[n]]).transpose()
    return th,phi,w

def spherical_wave_expansion_it(f, r, order):
    """find grid closest to requested order"""
    n = abs(array(Lebedev_L2max) - order).argmin()
    if order != Lebedev_L2max[n]:
        print "No grid for order %s, using grid which integrates up to L2max = %s exactly instead." \
            % (order, Lebedev_L2max[n])
    th,phi,w = array(LebedevGridPoints[Lebedev_Npts[n]]).transpose()
    x = outerN(r, sin(th)*cos(phi))
    y = outerN(r, sin(th)*sin(phi))
    z = outerN(r, cos(th))
    fxyz = f(x,y,z)
    sph_it = spherical_harmonics_it(th,phi)
    for Ylm,l,m in sph_it:
        flm = 4.0*pi*sum(w*fxyz*Ylm.conjugate(), axis=-1)
        yield flm, l, m
        if m == -Lebedev_L2max[n]:
            raise StopIteration

def sph_synthesis(flm_it, r, th, phi, lmax=17):
    ylm_it = spherical_harmonics_it(th,phi)
    f = 0.0
    for (flm,l1,m1),(ylm,l,m) in izip(flm_it, ylm_it):
        assert l==l1 and m==m1
        f += flm*ylm
    return f


if __name__ == "__main__":
#    check_lebedev_quadrature(len(Lebedev_Npts)-1)
    test_lebedev_quadrature()

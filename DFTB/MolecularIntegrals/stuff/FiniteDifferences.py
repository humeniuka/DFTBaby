"""
stuff from FiniteDifferences.py that's not needed anymore
"""


def radial_laplacian_boundary_r0(Nr, rm):
    """
    matrix representation of radial part of Laplacian

      __2        1  d^2
      \/ f(r) = --- ---- (r f)
                 r  dr^2

    on a non-equidistant radial grid. 
    It is assumed that f satisfies the boundary condition

        f(r --> 0) = 0             r=0 <-> z=1


    Parameters
    ----------
    Nr            :  int, number of radial grid points
    rm            :  float, scaling factor, 1/2 * Slater radius

    Returns
    -------
    lap_rad       :  nrad x nrad matrix, representation of the
                     Laplacian on the radial grid
    """
    # number of radial grid points
    n = Nr
    # grid is equidistant in z-coordinates
    k = np.array(range(1,n+1))
    # grid points on interval [-1,1]
    zr = k/(n+1.0)
    xr = np.cos(zr * np.pi)
    # radial grid points on interval [0,infinity]
    r = rm * (1+xr)/(1-xr)

    omx = 1-xr
    opx = 1+xr

    # The transformation of partial derivatives from
    # r-coordinates to z-coordinates is accomplished as
    #  d^2 u      d u      d^2 u
    #  ----- = c1 --- + c2 -----
    #  d r^2      d z      d z^2

    # boundaries
    #  z=0 <-> x=1  <-> r=+infinity
    #  z=1 <-> x=-1 <-> r=0

    c1 = 1.0/(4.0*np.pi*rm**2) * omx**(2.5) * (1+opx) / opx**(1.5)
    c2 = 1.0/(4.0*np.pi**2*rm**2) * omx**3 / opx

    # separation between equidistant points
    h = 1.0/(n+1)

    # operators d/dz and d^2/dz^2
    D1 = np.zeros((n,n))
    D2 = np.zeros((n,n))
    # non-centered seven-point formulae for i=0
    # D^1 u_0
    D1[0,0:7] = np.array([ -1764.0, +4320.0,  -5400.0,  +4800.0, -2700.0,  +864.0, -120.0  ])/(720.0*h)
    # D^2 u_0
    D2[0,0:7] = np.array([ +1624.0, -6264.0, +10530.0, -10160.0, +5940.0, -1944.0, +274.0 ])/(360.0*h**2)
    # non-centered seven-point formulae for i=1
    # D^1 u_1
    D1[1,0:7] = np.array([ -120.0, -924.0, +1800.0, -1200.0, +600.0, -180.0, +24.0 ])/(720.0*h)
    # D^2 u_1
    D2[1,0:7] = np.array([ +274.0, -294.0, -510.0, +940.0, -570.0, +186.0, -26.0 ])/(360.0*h**2)
    # non-centered seven-point formulae for i=2
    # D^1 u_2
    D1[2,0:7] = np.array([ +24.0, -288.0, -420.0, +960.0, -360.0, +96.0, -12.0])/(720.0*h)
    # D^2 u_2
    D2[2,0:7] = np.array([ -26.0, +456.0, -840.0, +400.0,  +30.0, -24.0,  +4.0])/(360.0*h**2)
    # centered seven-point formulae for i=3,...,n-4
    A7c1 = np.array([-12.0, +108.0, -540.0,    0.0, +540.0, -108.0, +12.0])/(720.0*h)
    A7c2 = np.array([ +4.0,  -54.0, +540.0, -980.0, +540.0,  -54.0,  +4.0])/(360.0*h**2)
    for i in range(3, n-3):
        D1[i,i-3:i+4] = A7c1
        D2[i,i-3:i+4] = A7c2
    # Continue with centered seven-point formulae for i=n-3,n-2,n-1 assuming that
    # at z=1 (i.e. r=0) the boundary conditions f(n)=0, f(n+1)=0, f(n+2)=0 are fulfilled. 
    for i in range(n-3, n):
        D1[i,i-3:] = A7c1[:n-i+3]
        D2[i,i-3:] = A7c2[:n-i+3]

    # finite difference formula converts differential operators into matrices
    
    # radial part of Laplacian
    #  __2        1  d^2
    #  \/ f(r) = --- ---- (r f)
    #             r  dr^2
    #
    # d^2/dr^2
    d2dr2 = np.dot(np.diag(c1), D1) + np.dot(np.diag(c2), D2)
    #
    lap_rad =  np.dot(np.diag(1.0/r), np.dot(d2dr2, np.diag(r)))

    return lap_rad

def radial_laplacian_boundary_rinf(Nr, rm):
    """
    matrix representation of radial part of Laplacian

      __2        1  d^2
      \/ f(r) = --- ---- (r f)
                 r  dr^2

    on a non-equidistant radial grid. 
    It is assumed that f satisfies the boundary condition

        f(r --> oo) = 0             r=inf <-> z=0


    Parameters
    ----------
    Nr            :  int, number of radial grid points
    rm            :  float, scaling factor, 1/2 * Slater radius

    Returns
    -------
    lap_rad       :  nrad x nrad matrix, representation of the
                     Laplacian on the radial grid
    """
    # number of radial grid points
    n = Nr
    # grid is equidistant in z-coordinates
    k = np.array(range(1,n+1))
    # grid points on interval [-1,1]
    zr = k/(n+1.0)
    xr = np.cos(zr * np.pi)
    # radial grid points on interval [0,infinity]
    r = rm * (1+xr)/(1-xr)

    omx = 1-xr
    opx = 1+xr

    # The transformation of partial derivatives from
    # r-coordinates to z-coordinates is accomplished as
    #  d^2 u      d u      d^2 u
    #  ----- = c1 --- + c2 -----
    #  d r^2      d z      d z^2
    c1 = 1.0/(4.0*np.pi*rm**2) * omx**(2.5) * (1+opx) / opx**(1.5)
    c2 = 1.0/(4.0*np.pi**2*rm**2) * omx**3 / opx

    # separation between equidistant points
    h = 1.0/(n+1)

    # operators d/dz and d^2/dz^2
    D1 = np.zeros((n,n))
    D2 = np.zeros((n,n))

    # centered seven-point formulae for i=3,...,n-4
    A7c1 = np.array([-12.0, +108.0, -540.0,    0.0, +540.0, -108.0, +12.0])/(720.0*h)
    A7c2 = np.array([ +4.0,  -54.0, +540.0, -980.0, +540.0,  -54.0,  +4.0])/(360.0*h**2)

    # Assume boundary conditions f(-3)=0, f(-2)=0, f(-1)=0
    for i in range(0, 3):
        D1[i,0:i+4] = A7c1[3-i:]
        D2[i,0:i+4] = A7c2[3-i:]
    # centered seven-point formulae for i=3,...,n-4
    for i in range(3, n-3):
        D1[i,i-3:i+4] = A7c1
        D2[i,i-3:i+4] = A7c2
    # non-centered seven-point formulae for i=n-3
    # D^1 u_{n-3}   ~   D^1 u_4
    D1[n-3,n-7:] = np.array([+12.0, -96.0, +360.0, -960.0, +420.0, +288.0, -24.0 ])/(720.0*h)
    # D^2 u_{n-3}   ~   D^2 u_4
    D2[n-3,n-7:] = np.array([ +4.0, -24.0,  +30.0, +400.0, -840.0, +456.0, -26.0 ])/(360.0*h**2)
    # non-centered seven-point formulae for i=n-2
    # D^1 u_{n-2}   ~   D^1 u_5
    D1[n-2,n-7:] = np.array([ -24.0, +180.0, -600.0, +1200.0, -1800.0, +924.0, +120.0])/(720.0*h)
    # D^2 u_{n-2}   ~   D^2 u_5
    D2[n-2,n-7:] = np.array([ -26.0, +186.0, -570.0,  +940.0,  -510.0, -294.0, +274.0])/(360.0*h**2)
    # non-centered seven-point formulae for i=n-1
    # D^1 u_{n-1}   ~   D^1 u_6
    D1[n-1,n-7:] = np.array([ +120.0,  -864.0, +2700.0,  -4800.0,  +5400.0, -4320.0, +1764.0])/(720.0*h)
    # D^2 u_{n-1}   ~   D^2 u_6
    D2[n-1,n-7:] = np.array([ +274.0, -1944.0, +5940.0, -10160.0, +10530.0, -6264.0, +1624.0])/(360.0*h**2)

    # finite difference formula converts differential operators into matrices
    
    # radial part of Laplacian
    #  __2        1  d^2
    #  \/ f(r) = --- ---- (r f)
    #             r  dr^2
    #
    # d^2/dr^2
    d2dr2 = np.dot(np.diag(c1), D1) + np.dot(np.diag(c2), D2)
    #
    lap_rad =  np.dot(np.diag(1.0/r), np.dot(d2dr2, np.diag(r)))

    return lap_rad



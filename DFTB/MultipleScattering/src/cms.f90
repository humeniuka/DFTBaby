!
! Continuum Multiple Scattering (CMS) Method 
!
! implementation of Dill and Dehmer's continuum multiple scattering method
!
! Numbers of equations that are mentioned in the source code refer to Ref. [1].
!
! References
! ----------
! [1] D.Dill, J.Dehmer, 
!     "Electron-molecule scattering and molecular photoionization usi!ng the multiple-scattering method",
!     J. Chem. Phys. 61, 692 (1974)
! [2] M.Danos, L.Maximon
!     "Multipole Matrix Elements of the Translation Operator"
!     Journal of Mathematical Physics 6, 766 (1965)
!

subroutine jn_coeffs(l1,m1, l2,m2, vec, jcoeff, ncoeff)
  !
  ! coefficients for expanding spherical Bessel functions of
  ! first and second kind around different center. The new center
  ! differs j from the old center i by the shift vector R_ij.
  !
  ! eqn. 14 in Dill & Dehmer (1974) (Ref. [1])
  !
  ! Parameters
  ! ----------
  ! l1,m1      : angular momentum quantum numbers, 0 <= l1, -l1 <= m1 <= +l1
  ! l2,m2      : angular momentum quantum numbers, 0 <= l2, -l2 <= m1 <= +l2
  ! vec        : cartesian shift vector R_ij, vec = [x,y,z]
  !
  ! Returns
  ! -------
  ! jcoeff,     : complex, coefficients J^{ij}_{L1,L2} and N^{ij}_{L1,L2}
  ! ncoeff        in the notation of Ref. [1]
  !
  implicit none
  ! ... input variables ...
  integer, intent(in) :: l1,m1, l2,m2
  double precision, intent(in) :: vec(3)
  ! ... output variables ...
  complex*16, intent(out) :: jcoeff, ncoeff
  ! ... local variables ...
  ! imaginary unit
  complex*16, parameter :: imag = (0, 1)
  ! pi = 3.14...
  double precision, parameter :: pi = 4.d0*datan(1.d0)
  ! spherical harmonics Y_{l,m}
  double precision :: r, th, ph
  integer :: l, m, lmin, lmax, lrng
  ! offsets where Wigner 3j symbols start
  integer :: jmin0, jmax0, jminm, jmaxm
  complex*16 :: y(l1+l2+1)
  integer :: kfn
  double precision :: eta
  ! spherical Bessel functions j_l(r) and n_l(r)
  double precision :: bessel_j(l1+l2+1), bessel_y(l1+l2+1)
  ! r-derivatives of j_l(r) and n_l(r)
  double precision :: fcp(l1+l2+1), gcp(l1+l2+1)
  integer :: ifail
  double precision :: w3j0(l1+l2+1), w3jm(l1+l2+1)
  complex*16 :: c
  
  ! spherical coordinates of cartesian vector
  call cart2sph(vec(1),vec(2),vec(3), r,th,ph)

  ! The Wigner 3j symbol
  !   (l  l2 l1 )
  !   (m  m2 -m1)
  ! limits the sum over l and m to terms where
  !  (1)  m = m1-m2
  !  (2) max(|m|,|l1-l2|) <= l <= l1+l2
  m = m1-m2
  lmin = max(abs(m), abs(l1-l2))
  lmax = l1+l2

  ! compute spherical harmonics for l=|m|,...,lmin,...,lmax
  ! y(1),y(2),...,y(|m|) are zero, y(|m|+1) = Y_{|m|,|m|} etc.
  call sphharm_const_m(th,ph, m, lmax, y)
  
  ! compute spherical Bessel functions
  lrng = lmax-lmin
  ! select spherical Bessel functions (kfn=1)
  kfn = 1
  call coul90_scalar(r, 0.d0, lmax, kfn, &
       bessel_j, bessel_y, fcp, gcp, &
       ifail)
  
  ! compute Wigner 3j symbols
  !
  !  ( l  l2 l1 )
  !  ( 0  0  0  ) is zero except for jmin0 <= l <= jmax0
  ! Therefore the element w3j0(1) corresponds to l=jmin0
  ! and w3j0(l-jmin0+1) corresponds to l.
  !
  ! Similary
  !  ( l  l2 l1  )
  !  ( m  m2 -m1 ) is zero except for jminm <= l <= jmaxm
  ! Therefore the element w3jm(l-jminm+1) corrsponds to l.
  call wigner3j(w3j0, jmin0, jmax0, l2, l1, 0, 0, 0, ifail)
  call wigner3j(w3jm, jminm, jmaxm, l2, l1, m, m2, -m1, ifail)

  jcoeff = (0.d0, 0.d0)
  ncoeff = (0.d0, 0.d0)
  do l=lmin,lmax
     c = imag**(l+l2-l1) * (-1)**m1 &
          * sqrt(4.d0*pi*(2*l+1)*(2*l1+1)*(2*l2+1)) &
          * w3j0(l-jmin0+1) * w3jm(l-jminm+1) &
          * y(l+1)
     jcoeff = jcoeff + c * bessel_j(l+1)
     ncoeff = ncoeff + c * bessel_y(l+1)
  enddo
end subroutine jn_coeffs

subroutine muffin_tin_potential( &
     vconst, chargeIII, nat, &
     rhoIII, radii, atpos, &
     nspl, kspl, knots, coeffs, &
     npts, x_arr,y_arr,z_arr, &
     v_arr)
  !
  ! evaluate the muffin tin potential at the cartesian coordinates x,y,z
  !
  ! In each region the potential has a form for which the Schroedinger
  ! equation can be solve exactly:
  !
  !   region I     -   atomic potentials, which are spherically symmetric
  !   region II    -   constant potential
  !   region III   -   Coulomb potential -chargeIII/r
  !
  ! At the boundaries between regions the potential has a discontinuity.
  ! The atomic potentials in region I are provided as splines of
  ! degree k with uniformly spaced knots. 
  !
  ! Parameters
  ! ----------
  ! vconst       : float, constant potential in region II (in Hartree)
  ! chargeIII    : float > 0, asymptotic monopole charge of effective potential,
  !                Veff(r) ----> -chargeIII/r  for r --> oo, typically charge=+1 for
  !                a neutral molecule. This is the charge felt by an electron in region III. 
  ! nat          : int > 0, number of atomic centers
  ! rhoIII       : float, radius of molecular sphere separating region II from region III (in bohr)
  ! radii        : float array of shape (nat), radii of atomic spheres separating region I from region II (in bohr),
  !                radii[i] = rho_i
  ! atpos        : float array of shape (3,nat) with cartesian positions of atomic centers (in bohr)
  ! nspl         : integer, giving the total number of knots for
  !                the splines of the atomic radial potentials
  ! kspl         : integer, giving the degree of the splines
  ! knots        : array of shape (nspl,nat), which contains the B-spline knots,
  !                the knots have to be uniformly spaced.
  ! coeffs       : array of shape (nspl,nat), which contains the B-spline coefficients
  ! npts         : integer, number of points
  ! x_arr        : arrays of shape (npts) with cartesian coordinates of point (in bohr)
  ! y_arr
  ! z_arr
  
  ! Returns
  ! -------
  ! v_arr        : array of shape (npts), muffin tin potential V(x,y,z) evaluated at
  !                the grid points
  !
  implicit none
  ! ... input variables ...
  double precision,  intent(in) :: vconst, chargeIII
  integer, intent(in) :: nat
  double precision, intent(in) :: rhoIII, radii(nat), atpos(3,nat)
  integer, intent(in) :: nspl, kspl
  double precision, intent(in) :: knots(nspl,nat), coeffs(nspl,nat)
  integer, intent(in) :: npts
  double precision, intent(in) :: x_arr(npts), y_arr(npts), z_arr(npts)
  ! ... output variables ...
  double precision, intent(out) :: v_arr(npts)
  ! ... local variables ...
  ! enumerate grid points
  integer :: j
  ! enumerate atomic centers
  integer :: i
  ! cartesian coordinates
  double precision :: x,y,z, r
  ! potential V(x,y,z)
  double precision :: v
  ! error code of spline evaluation
  integer :: ier
  
  do j=1,npts
     x = x_arr(j)
     y = y_arr(j)
     z = z_arr(j)
     
     ! by default the value of V in the region II is assigned
     v = vconst
     
     ! region I
     do i=1,nat
        ! distance to center i
        r = sqrt((x-atpos(1,i))**2+(y-atpos(2,i))**2+(z-atpos(3,i))**2)
        if (r < radii(i)) then
           ! evaluate spline of radial potential i
           call splev_uniform(knots(:,i),nspl,coeffs(:,i),kspl, &
                r, v, ier)
           ! The point belongs to region I, so go directly to
           ! assigning this value
           go to 100
        endif
     enddo

     ! region III
     r = sqrt(x*x+y*y+z*z)
     if (rhoIII < r) then
        v = -chargeIII/r
     endif

100  v_arr(j) = v
  enddo
end subroutine muffin_tin_potential

subroutine matching_matrix( &
     energy, vconst, chargeIII, nat, lmax, &
     rhoIII, radii, atpos, &
     nspl, kspl, fIknots, fIcoeffs, &
     debug, &
     matchM, rhs )
  !
  ! construct the matrix for the inhomogeneous system of linear
  ! equations (16) and (17) in Ref. [1]
  !
  ! The coefficients `x` of the CMS wavefunction in the regions I and II
  ! are determined by solving a inhomogeneous system of linear equations
  !
  !      M.x = rhs
  !
  ! The right-hand side contains the coefficients of the wavefunction in region III
  ! and determines its asymptotic form (asymptotic angular momentum). Since there
  ! are `nlm=(lmax+1)^2` angular momentum channels and `nat` atomic centers plus the origin,
  ! the matching matrix has the dimension `ndim=(nat+1)*(lmax+1)^2`, the right hand side
  ! has `ndim` rows and `nlm` columns, one column for each angular momentum channel (l,m).
  !
  ! Parameters
  ! ----------
  ! energy       : float, energy of continuum orbital (in Hartree)
  ! vconst       : float, constant potential in region II (in Hartree)
  ! chargeIII    : float > 0, asymptotic monopole charge of effective potential,
  !                Veff(r) ----> -chargeIII/r  for r --> oo, typically charge=+1 for
  !                a neutral molecule. This is the charge felt by an electron in region III. 
  ! nat          : int > 0, number of atomic centers
  ! lmax         : int > 0, largest angular momentum, l=0,1,...,lmax
  ! rhoIII       : float, radius of molecular sphere separating region II from region III (in bohr)
  ! radii        : float array of shape (nat), radii of atomic spheres separating region I from region II (in bohr),
  !                radii[i] = rho_i
  ! atpos        : float array of shape (3,nat) with cartesian positions of atomic centers (in bohr)
  ! nspl         : integer, giving the total number of knots for
  !                the splines of the atomic radial wavefunctions
  ! kspl         : integer, giving the degree of the splines
  ! fIknots      : array of shape (nspl,lmax+1,nat), which contains the B-spline knots,
  !                the knots have to be uniformly spaced. fIknots(:,l,i) are the knots belonging
  !                to the radial wavefunction with angular momentum `l` on atom `i`
  ! fIcoeffs     : array of shape (nspl,lmax+1,nat), which contains the B-spline coefficients
  !
  ! Optional
  ! --------
  ! debug        : > 0, print additional information
  !
  ! Returns
  ! -------
  ! matchM       : complex array of shape (ndim, ndim), matching matrix
  ! rhs          : complex array of shape (ndim, nlm), right-hand side
  implicit none
  ! ... input variables ...
  double precision, intent(in) :: energy, vconst, chargeIII
  integer, intent(in) :: nat, lmax
  double precision, intent(in) :: rhoIII, radii(nat), atpos(3,nat)
  integer, intent(in) :: nspl, kspl
  double precision, intent(in) :: fIknots(nspl,0:lmax,nat), fIcoeffs(nspl,0:lmax,nat)
  ! ... optional ...
  !f2py integer optional, intent(in) :: debug = 0
  integer, intent(in) :: debug
  ! ... output variables ...
  complex*16, intent(out) :: matchM( (lmax+1)*(lmax+1)*(nat+1), (lmax+1)*(lmax+1)*(nat+1) )
  complex*16, intent(out) :: rhs( (lmax+1)*(lmax+1)*(nat+1), (lmax+1)*(lmax+1) )
  ! ... local variables ...
  !
  ! pi = 3.14...
  double precision, parameter :: pi = 4.d0*datan(1.d0)
  ! fI           : array of shape (lmax+1,nat)
  !                fI(l,i) is the value of the radial wavefunction of atom i with angular
  !                momentum l at the border of region I,
  !                               (i)
  !                    fI[l,i] = f   (r=rho )
  !                               l        i
  ! fIp          : float array of shape (lmax+1,nat)
  !                fIp(l,i) is the r-derivative of the radial wavefunction evaluated at
  !                the border of region I
  !
  !                               d f^(i)_l
  !                    fIp[l,i] = --------- (r = rho )
  !                                 d r             i
  double precision :: fI(0:lmax,nat), fIp(0:lmax,nat)
  ! radius of an atomic sphere
  double precision :: rho_i
  integer :: nlm, ndim
  ! wave vectors in region I and II
  double precision :: k, kappa
  ! Sommerfeld parameter eta = -Z/k
  double precision :: eta
  ! Wronskians
  double precision :: wrQ(0:lmax,nat)
  double precision :: wrR(0:lmax)
  double precision :: wrS(0:lmax)
  double precision :: wronsk_numer, wronsk_denom
  ! enumerates atoms
  integer :: i,j
  ! enumerate rows and columns of matching matrix
  integer :: ii,jj
  ! enumerate angular momentum channels (0,0), ..., (lmax,-lmax)
  integer :: nci, ncj
  ! angular momentum quantum numbers
  integer :: li,mi, lj,mj
  ! enumerates angular momenta
  integer :: l
  ! spherical Bessel functions of 1st and 2nd kind
  double precision :: jII(0:lmax), nII(0:lmax)
  ! and their r-derivatives at a certain radius
  double precision :: jIIp(0:lmax), nIIp(0:lmax)
  ! regular and irregular Coulomb functions
  double precision :: fIII(0:lmax), gIII(0:lmax)
  ! and their r-derivatives
  double precision :: fIIIp(0:lmax), gIIIp(0:lmax)
  ! exit status for calculation Bessel and Coulomb functions
  integer :: ifail
  ! abbreviation for k*rhoIII
  double precision :: x
  ! arrays with l- and m-values
  integer :: l_arr((lmax+1)*(lmax+1)), m_arr((lmax+1)*(lmax+1))
  ! difference vector of atom positions
  double precision :: vecji(3)
  ! coefficients for reexpansion of spherical Bessel functions around
  ! different center
  complex*16 :: jcoeff, ncoeff
  ! sqrt(k/pi)
  double precision :: sqkpi

  ! evaluate atomic radial wavefunctions at the radii of the atomic
  ! spheres. The radial wavefunctions are given as B-splines. 
  do i=1,nat
     rho_i = radii(i)
     do l=0,lmax
        ! evaluate radial wavefunctions
        call splev_uniform(fIknots(1:nspl,l,i),nspl,fIcoeffs(1:nspl,l,i),kspl, &
             rho_i, fI(l,i), ifail)
        ! and the r-derivative
        call splev_deriv_uniform(fIknots(1:nspl,l,i),nspl,fIcoeffs(1:nspl,l,i),kspl, &
             rho_i, fIp(l,i), ifail)
     enddo
  enddo

  ! For each l-value there are 2*l+1 quantum numbers m=-l,...,0,...,l.
  ! For l=0,1,...,lmax there are nlm = sum_l=0^lmax 2*l+1 = (lmax+1)^2
  ! angular momentum channels.
  nlm = (lmax+1)**2
  ! dimension of the matching matrix
  ndim = (nat+1)*nlm
  ! E = 1/2 k^2
  k = sqrt(2.d0*energy)
  ! wavevector in region II
  kappa = sqrt(k*k - 2*vconst)
  if (debug > 0) write(*,*) "length of wavevector in region I & III  k = ", k
  if (debug > 0) write(*,*) "length of wavevector in region II   kappa = ", kappa

  ! compute Coulomb (kfn=0) functions fIII and gIII for l=0,...,lmax
  eta = -chargeIII/k
  call coul90_scalar(k*rhoIII, eta, lmax, 0, &
       fIII, gIII, fIIIp, gIIIp, ifail)
  ! We need the regular Coulomb function
  !  III
  ! f  (k*r) = 1/(k*r) F (eta=-Z/k; k*r)
  !  l                  l
  ! and the irregular Coulomb function
  !  III
  ! g  (k*r) = 1/(k*r) G (eta=-Z/k; k*r)
  !  l                  l
  ! but coul90_scalar returns F_l(k*r), G_l(k*r) and its derivatives
  ! with respect to x=(k*r).
  x = k*rhoIII
  ! Therefore we need to transform the derivatives according to
  !   f(x) = 1/x F(x) ---> f'(x) = -F(x)/x^2 + F'(x)/x
  fIIIp = -fIII/x**2 + fIIIp/x
  gIIIp = -gIII/x**2 + gIIIp/x
  fIII = fIII/x
  gIII = gIII/x
  
  ! intermediate variables needed for constructing M
  if (debug > 0) write(*,*) "computing Wronskians ..."
  ! ratios of Wronskians
  ! --------------------
  ! The definition of the Wronskian of two functions a(r) and b(r)
  ! is
  !   [a,b] = a*b' - a'*b
  ! where the derivative (') is with respect to r.
  ! Because of the chain rule we get
  !   [a(k*r),b(r)] = a(k*r)*b'(r) - k*a'(k*r)*b(r)
  !
  !          [n_l(kappa*rho_i), f_l^i(rho_i)]
  ! wrQ    = --------------------------------
  !    l,i   [j_l(kappa*rho_i), f_l^i(rho_i)]
  do i=1,nat
     ! radius of atomic sphere i
     rho_i = radii(i)
     ! compute spherical Bessel (kfn=1) functions j_l and n_l for l=0,...,lmax
     call coul90_scalar(kappa*rho_i, 0.d0, lmax, 1, &
          jII, nII, jIIp, nIIp, ifail)
     do l=0,lmax
        ! the factor kappa comes from the chain rule
        ! since the derivatives are with respect to kappa*r
        wronsk_numer = nII(l) * fIp(l,i) - kappa*nIIp(l) * fI(l,i)
        wronsk_denom = jII(l) * fIp(l,i) - kappa*jIIp(l) * fI(l,i)
        wrQ(l,i) = wronsk_numer / wronsk_denom
     enddo
  enddo
  !
  !        [j_l(kappa*rhoIII), g_l^III(k*rhoIII)]
  ! wrR  = --------------------------------------
  !    l   [n_l(kappa*rhoIII), g_l^III(k*rhoIII)]
  !
  ! compute spherical Bessel functions at r=rhoIII
  call coul90_scalar(kappa*rhoIII, 0.d0, lmax, 1, &
       jII, nII, jIIp, nIIp, ifail)
  do l=0,lmax
     wronsk_numer = jII(l) * k * gIIIp(l) - kappa * jIIp(l) * gIII(l)
     wronsk_denom = nII(l) * k * gIIIp(l) - kappa * nIIp(l) * gIII(l)
     wrR(l) = wronsk_numer / wronsk_denom
  enddo
  ! needed for constructing right hand side of equation  
  !
  !        [f_l^III(k*rhoIII), g_l^III(k*rhoIII)]
  ! wrS  = --------------------------------------
  !    l   [n_l(kappa*rhoIII), g_l^III(k*rhoIII)]
  do l=0,lmax
     wronsk_numer = fIII(l) * k * gIIIp(l) - k     * fIIIp(l)  * gIII(l)
     wronsk_denom = nII(l)  * k * gIIIp(l) - kappa * nIIp(l)   * gIII(l)
     wrS(l) = wronsk_numer / wronsk_denom
  enddo

  ! ordering of angular momentum channels 
  call enumerate_lm(lmax, l_arr, m_arr)
  
  if (debug > 0) write(*,*) "matching matrix ..."
  ! Fill in elements of matching matrix
  ! Rows and columns are enumerated by multiindices ii and jj,
  ! which stand for (i,li,mi) and (j,lj,mj).
  matchM(:,:) = (0.d0, 0.d0)
  
  ! columns (in Fortran the row index runs faster)
  jj = 1
  ! enumerate atomic centers ji=0 means the origin)
  do j=0,nat
     ! enumerate angular momentum channels (lj,mj) on center j
     do ncj=1,nlm
        lj = l_arr(ncj)
        mj = m_arr(ncj)

        !write(*,*) "filling column ", jj, " (of ", ndim, ")"
        ! rows
        ii = 1
        ! enumerate atomic centers (i=0 means the origin)
        do i=0,nat
           ! enumerate angular momentum channels (li,mi) on center i
           do nci=1,nlm
              li = l_arr(nci)
              mi = m_arr(nci)

              !---------------------------------------------------------------------!
              ! distinction of cases for matrix elements
              if (ii == jj) then
                 ! diagonal elements of matrix
                 if (i == 0) then
                    ! (i=0,j=0) diagonal block, origin
                    matchM(ii,jj) = wrR(li)
                 else
                    ! (i == j, i > 0) diagonal block, atomic centers
                    matchM(ii,jj) = wrQ(li,i)
                 endif
              else
                 ! off-diagonal elements
                 if ((i == 0).and.(j > 0)) then
                    call jn_coeffs(lj,mj, li,mi, -kappa*atpos(:,j), jcoeff, ncoeff)
                    matchM(ii,jj) = jcoeff
                 else if ((j == 0).and.(i > 0)) then
                    call jn_coeffs(lj,mj, li,mi, kappa*atpos(:,i), jcoeff, ncoeff)
                    matchM(ii,jj) = jcoeff
                 else if ((i > 0).and.(j > 0).and.(i /= j)) then
                    ! vector from atom j to atom i
                    vecji = atpos(:,i) - atpos(:,j)
                    call jn_coeffs(lj,mj, li,mi, kappa*vecji, jcoeff, ncoeff)
                    matchM(ii,jj) = ncoeff
                 endif
              endif
              !---------------------------------------------------------------------!
              
              ! increase row counter
              ii = ii+1
           enddo
        enddo
        ! increase column counter
        jj = jj+1
     enddo
  enddo
  
  sqkpi = sqrt(k/pi)
  ! Fill in right-hand side of matching equation
  rhs(:,:) = (0.d0, 0.d0)
  ! column for channel (l,m)
  do ncj=1,nlm
     ! enumerate asymptotic angular momentum channels (l,m) 
     lj = l_arr(ncj)
     mj = m_arr(ncj)
     ! eqn. (23)
     rhs(ncj,ncj) = sqkpi * wrS(lj)
  enddo
  
end subroutine matching_matrix

subroutine wavefunctions( &
     energy, vconst, chargeIII, nat, lmax, &
     rhoIII, radii, atpos, &
     nspl, kspl, fIknots, fIcoeffs, &
     nc, sol, lm0, &
     npts, x_arr, y_arr, z_arr, &
     debug, &
     kmat, wfn_arr )
  !
  ! evaluate the CMS wavefunctions on a grid using the coefficients
  ! obtained by solving the matching equations.
  !
  ! The atomic radial wavefunction for region I have to be
  ! provided as B-splines. 
  !
  ! Parameters
  ! ----------
  ! energy       : float, energy of continuum orbital (in Hartree)
  ! vconst       : float, constant potential in region II (in Hartree)
  ! chargeIII    : float > 0, asymptotic monopole charge of effective potential,
  !                Veff(r) ----> -chargeIII/r  for r --> oo, typically charge=+1 for
  !                a neutral molecule. This is the charge felt by an electron in region III. 
  ! nat          : int > 0, number of atomic centers
  ! lmax         : int > 0, largest angular momentum, l=0,1,...,lmax
  ! rhoIII       : float, radius of molecular sphere separating region II from region III (in bohr)
  ! radii        : float array of shape (nat), radii of atomic spheres separating region I from region II (in bohr),
  !                radii[i] = rho_i
  ! atpos        : float array of shape (3,nat) with cartesian positions of atomic centers (in bohr)
  ! nspl         : integer, giving the total number of knots for
  !                the splines of the atomic radial wavefunctions
  ! kspl         : integer, giving the degree of the splines
  ! fIknots      : array of shape (nspl,lmax+1,nat), which contains the B-spline knots,
  !                the knots have to be uniformly spaced. fIknots(:,l,i) are the knots belonging
  !                to the radial wavefunction with angular momentum `l` on atom `i`
  ! fIcoeffs     : array of shape (nspl,lmax+1,nat), which contains the B-spline coefficients
  ! nc           : integer, number of angular momentum channels,
  !                number of columns in `sol`
  ! sol          : array of shape (ndim,nc) with solution vectors of matching equations, M.sol = rhs
  ! lm0          : integer, index of first solution vector in `sol` (starting at 0)
  ! npts         : integer, number of points
  ! x_arr        : arrays of shape (npts) with cartesian coordinates of point (in bohr)
  ! y_arr
  ! z_arr
  !
  ! Optional
  ! --------
  ! debug        : > 0, print additional information
  !
  ! Returns
  ! -------
  ! kmat         : rows of K-matrix belonging to the solution vectors
  !                in `sol`
  ! wfn_arr      : array of shape (nlm,:) with values of wavefunctions on the grid,
  !                wfn_arr(lm,:) is the wavefunction with asymptotic angular momentum (l,m)
  !
  implicit none
  ! ... input variables ...
  double precision, intent(in) :: energy, vconst, chargeIII
  integer, intent(in) :: nat, lmax
  double precision, intent(in) :: rhoIII, radii(nat), atpos(3,nat)
  integer, intent(in) :: nspl, kspl
  double precision, intent(in) :: fIknots(nspl,0:lmax,nat), fIcoeffs(nspl,0:lmax,nat)
  integer, intent(in) :: nc
  complex*16, intent(in) :: sol((lmax+1)*(lmax+1)*(nat+1),nc)
  integer, intent(in) :: lm0
  integer, intent(in) :: npts
  double precision, intent(in) :: x_arr(npts), y_arr(npts), z_arr(npts)
  ! ... optional ...
  !f2py integer optional, intent(in) :: debug = 0
  integer, intent(in) :: debug
  ! ... output variables ...
  complex*16, intent(out) :: kmat((lmax+1)*(lmax+1),nc)
  complex*16, intent(out) :: wfn_arr(nc,npts)
  ! ... local variables ...
  ! pi = 3.14...
  double precision, parameter :: pi = 4.d0*datan(1.d0)
  integer :: nlm, ndim
  ! wave vectors in region I and II
  double precision :: k, kappa
  ! Sommerfeld parameter eta = -Z/k
  double precision :: eta
  ! enumerate atoms
  integer :: i
  ! enumerate angular momenta
  integer :: l, l1
  ! multiindices (l,m)
  integer :: lm, lm1
  ! enumerate solution vectors
  integer :: ic
  ! arrays with l- and m-values
  integer :: l_arr((lmax+1)*(lmax+1)), m_arr((lmax+1)*(lmax+1))
  ! sqrt(k/pi)
  double precision :: sqkpi
  !
  ! coefficients of wavefunctions
  ! -----------------------------
  !   in region III  (eqn. (8))
  !
  !                       III  III           III  III
  ! Psi   (r) = sum    [ A    f  (k r )  +  B    g (k r ) ] Y  (th ,ph )
  !    III         l,m    l,m  l     0       l,m  l    0     l,m  0   0
  !
  ! where (r0,th0,ph0) are the spherical coordinates relative to the origin
  complex*16 :: AIII((lmax+1)*(lmax+1),nc)
  complex*16 :: BIII((lmax+1)*(lmax+1),nc)
  !  in region I  (eqn. (2))
  !                   Ii   I
  ! Psi (r) = sum    A    f (r ) Y  (th ,ph )
  !    i         l,m  l,m  l  i   l,m  i   i
  !
  ! where (r_i,th_i,ph_i) are the spherical coordinates relative to the atom i
  complex*16 :: AI((lmax+1)*(lmax+1),nat,nc)
  !  in region II (eqn. (7))
  !                    II0                                  nat         IIi
  ! Psi  (r) = sum    A    j (kappa r ) Y  (th ,ph )  +  sum    sum    B    n (kappa r ) Y   (th ,ph )
  !    II         l,m  l,m  l        0   l,m  0   0         i=1    l,m  l,m  l        i   l,m   i   i
  complex*16 :: AII0((lmax+1)*(lmax+1),nc)
  complex*16 :: BII((lmax+1)*(lmax+1),nat,nc)
  ! shape of BII array (nlm,nat)
  integer :: sh(2)
  !
  ! Wronskians
  ! ----------
  ! To find the coefficients A^{I_i}_{l,m} and B^{III}_{l,m} from
  ! the solution of the matching conditions according to eqns.
  ! (18) and (19) we need to store some additional Wronskians
  !
  !                                 i
  ! wronskJF    = [j (kappa*rho ), f (rho )]
  !         l,i     l          i    l    i
  !
  double precision :: wronskJF(0:lmax,nat)
  !
  !                                 III
  ! wronskNG  = [n (kappa*rhoIII), g (kappa*rhoIII)]
  !         l     l                 l
  !        
  !                                 III
  ! wronskNF  = [n (kappa*rhoIII), f (kappa*rhoIII)]
  !         l     l                 l
  !
  double precision :: wronskNG(0:lmax), wronskNF(0:lmax)
  ! radius of atomic sphere i
  double precision :: rho_i
  ! spherical Bessel functions of first kind and second kind
  ! and their r-derivatives
  double precision :: jII(0:lmax), jIIp(0:lmax)
  double precision :: nII(0:lmax), nIIp(0:lmax)
  ! atomic radial functions and their r-derivatives at r=rho_i
  double precision :: fI(0:lmax), fIp(0:lmax)
  ! Coulomb functions and their r-derivatives
  double precision :: fIII(0:lmax), fIIIp(0:lmax)
  double precision :: gIII(0:lmax), gIIIp(0:lmax)
  integer :: ifail
  ! x = k*r
  double precision :: kr
  ! cartesian positions
  double precision :: x,y,z
  ! enumerate grid points
  integer :: j
  ! Which region does the point belong to (1=I,2=II or 3=III)
  integer :: region
  ! distance to origin or atomic center i
  double precision :: r0, ri
  ! cartesian coordinates relative to center i
  double precision :: xi,yi,zi
  ! angles of spherical coordinates relative to origin or center i
  double precision :: th0,ph0, thi,phi
  ! spherical harmonics Ylm
  complex*16 :: ysph((lmax+1)*(lmax+1))
  ! values of wavefunctions at (x,y,z)
  complex*16 :: wfn(nc)
  
  ! For each l-value there are 2*l+1 quantum numbers m=-l,...,0,...,l.
  ! For l=0,1,...,lmax there are nlm = sum_l=0^lmax 2*l+1 = (lmax+1)^2
  ! angular momentum channels.
  nlm = (lmax+1)**2
  ! dimension of the matching matrix
  ndim = (nat+1)*nlm
  ! E = 1/2 k^2
  k = sqrt(2.d0*energy)
  ! wavevector in region II
  kappa = sqrt(k*k - 2*vconst)
  sqkpi = sqrt(k/pi)

  eta = -chargeIII/k
  
  ! ---------- Wronskians --------------------------------------
  ! compute region I,II and III functions at the matching points
  call coul90_scalar(k*rhoIII, eta, lmax, 0, &
       fIII, gIII, fIIIp, gIIIp, ifail)
  ! We need the regular Coulomb function
  !  III
  ! f  (k*r) = 1/(k*r) F (eta=-Z/k; k*r)
  !  l                  l
  ! and the irregular Coulomb function
  !  III
  ! g  (k*r) = 1/(k*r) G (eta=-Z/k; k*r)
  !  l                  l
  ! but coul90_scalar returns F_l(k*r), G_l(k*r) and its derivatives
  ! with respect to x=(k*r).
  kr = k*rhoIII
  ! Therefore we need to transform the derivatives according to
  !   f(x) = 1/x F(x) ---> f'(x) = -F(x)/x^2 + F'(x)/x
  fIIIp = -fIII/kr**2 + fIIIp/kr
  gIIIp = -gIII/kr**2 + gIIIp/kr
  fIII = fIII/kr
  gIII = gIII/kr
  
  ! The jII and jIIp values are not needed.
  call coul90_scalar(kappa*rhoIII, 0.d0, lmax, 1, &
       jII, nII, jIIp, nIIp, ifail)

  do l=0,lmax
     wronskNG(l) = nII(l) * k * gIIIp(l) - kappa * nIIp(l) * gIII(l)
     wronskNF(l) = nII(l) * k * fIIIp(l) - kappa * nIIp(l) * fIII(l)
  enddo

  do i=1,nat
     rho_i = radii(i)
     call coul90_scalar(kappa*rho_i, 0.d0, lmax, 1, &
          jII, nII, jIIp, nIIp, ifail)
     do l=0,lmax
        ! evaluate radial wavefunctions
        call splev_uniform(fIknots(1:nspl,l,i),nspl,fIcoeffs(1:nspl,l,i),kspl, &
             rho_i, fI(l), ifail)
        ! and the r-derivative
        call splev_deriv_uniform(fIknots(1:nspl,l,i),nspl,fIcoeffs(1:nspl,l,i),kspl, &
             rho_i, fIp(l), ifail)
        wronskJF(l,i) = jII(l) * fIp(l) - kappa * jIIp(l) * fI(l)
     enddo
  enddo
  
  ! ----- Wavefunction coefficients ------------------------------
  ! Given the matrix of solution vectors `sol` we have to extract
  ! the coefficients for each region
  !  sol  ->  (AI,AII0,BII_i,AIII,BIII)
  ! ordering of angular momentum channels 
  call enumerate_lm(lmax, l_arr, m_arr)

  ! shape of BII array
  sh(1) = nlm
  sh(2) = nat

  AIII(:,:) = 0.d0
  do lm=lm0+1,nc+lm0
     ic = lm-lm0
     ! AIII coefficients according to eqn. (23)
     AIII(lm,ic) = sqkpi
     ! extract coefficients A^{II_0}_{l,m} and B^{II_i}_{l,m}
     ! from solution vector
     AII0(1:nlm,ic) = sol(1:nlm,ic)
     BII(1:nlm,1:nat,ic) = reshape(sol(nlm+1:ndim,ic), sh)
     ! find coefficients A^{I_i}_{l,m} according to eqn. (18)
     do i=1,nat
        rho_i = radii(i)
        do lm1=1,nlm
           l1 = l_arr(lm1)
           AI(lm1,i,ic) = BII(lm1,i,ic) / (kappa * rho_i**2 * wronskJF(l1,i))
        enddo
     enddo
     ! find coefficients B^{III}_{l,m} according to eqn. (19)
     do lm1=1,nlm
        l1 = l_arr(lm1)
        BIII(lm1,ic) = - ( AII0(lm1,ic) / (kappa * rhoIII**2) &
             + wronskNF(l1) * AIII(lm1,ic) ) / wronskNG(l1)
     enddo
  enddo
  
  ! rows of the K-matrix
  kmat = sqkpi * BIII

  ! --------------- evaluate wavefunctions -----------------
  do j=1,npts
     x = x_arr(j)
     y = y_arr(j)
     z = z_arr(j)

     ! --------------------------------------------------------
     ! First we need to find out which region the point (x,y,z)
     ! belongs to.

     ! distance to origin
     r0 = sqrt(x*x+y*y+z*z)
     if (rhoIII < r0) then
        ! Point lies outside molecular sphere.
        region = 3
     else
        ! By default the point lies in region II unless it lies
        ! in an atomic sphere
        region = 2
        
        atoms: do i=1,nat
           ! distance to atom i
           ri = sqrt((x-atpos(1,i))**2+(y-atpos(2,i))**2+(z-atpos(3,i))**2)
           rho_i = radii(i)
           if (ri < rho_i) then
              ! Point lies in region I and belong to atom i.
              region = 1
              exit atoms
           endif
        enddo atoms
     endif
     !write(*,*) j, "-th point belongs to region ", region
     ! ----------------------------------------------------------

     ! ----------------------------------------------------------
     ! evaluate wavefunction for all channels using the expressions
     ! appropriate for the region where (x,y,z) lies in
     ! 
     wfn(:) = (0.d0, 0.d0)
     
     if (region == 1) then
        !write(*,*) "    and atom ", i
        ! coordinates relative to atom i
        xi = x - atpos(1,i)
        yi = y - atpos(2,i)
        zi = z - atpos(3,i)
        ! convert them to spherical coordinates
        call cart2sph(xi,yi,zi, ri,thi,phi)
        ! spherical harmonics Y_lm
        call sphharm(thi,phi, lmax, ysph)

        ! evaluate radial wavefunctions
        do l=0,lmax
           call splev_uniform(fIknots(1:nspl,l,i),nspl,fIcoeffs(1:nspl,l,i),kspl, &
                ri, fI(l), ifail)
        enddo

        ! compose wavefunction according to eqn. (2)
        do lm=1,nlm
           l = l_arr(lm)
           wfn(:) = wfn(:) + AI(lm,i,:) * fI(l) * ysph(lm)
        enddo
           
     else if (region == 2) then
        ! spherical coordinates relative to origin
        call cart2sph(x,y,z, r0,th0,ph0)
        ! spherical harmonics Y_lm
        call sphharm(th0,ph0, lmax, ysph)
        ! spherical Bessel functions
        call coul90_scalar(kappa*r0, 0.d0, lmax, 1, &
             jII, nII, jIIp, nIIp, ifail)
        ! --- first half of eqn. (7) ---
        do lm=1,nlm
           l = l_arr(lm)
           wfn(:) = wfn(:) + AII0(lm,:) * jII(l) * ysph(lm)
        enddo

        do i=1,nat
           ! coordinates relative to atom i
           xi = x - atpos(1,i)
           yi = y - atpos(2,i)
           zi = z - atpos(3,i)
           ! convert them to spherical coordinates
           call cart2sph(xi,yi,zi, ri,thi,phi)
           ! spherical harmonics Y_lm
           call sphharm(thi,phi, lmax, ysph)
           ! spherical Bessel functions
           call coul90_scalar(kappa*ri, 0.d0, lmax, 1, &
                jII, nII, jIIp, nIIp, ifail)
           ! --- second half of eqn. (7) ----
           do lm=1,nlm
              l = l_arr(lm)
              wfn(:) = wfn(:) + BII(lm,i,:) * nII(l) * ysph(lm)
           enddo
        enddo
     else
        ! region III
        ! spherical coordinates relative to origin
        call cart2sph(x,y,z, r0,th0,ph0)
        ! spherical harmonics Y_lm
        call sphharm(th0,ph0, lmax, ysph)
        ! Coulomb functions
        call coul90_scalar(k*r0, eta, lmax, 0, &
             fIII, gIII, fIIIp, gIIIp, ifail)
        kr = k*r0
        fIII(:) = fIII(:)/kr
        gIII(:) = gIII(:)/kr
        
        ! --- eqn. (8) ---
        do lm=1,nlm
           l = l_arr(lm)
           wfn(:) = wfn(:) + (AIII(lm,:) * fIII(l) + BIII(lm,:) * gIII(l)) * ysph(lm)
        enddo
     endif

     wfn_arr(:,j) = wfn
  enddo
  
end subroutine wavefunctions

subroutine transition_dipoles( &
     energy, vconst, chargeIII, nat, lmax, &
     rhoIII, radii, atpos, &
     nspl, kspl, fIknots, fIcoeffs, &
     nc, sol, lm0, &
     npts, x_arr, y_arr, z_arr, w_arr, &
     norb, orbs_arr, &
     kmat, tdip_arr, norms2, projs2 )
  !
  ! evaluate the CMS wavefunctions on a grid using the coefficients
  ! obtained by solving the matching equations and compute the transition
  ! amplitudes between the bound orbitals (provided in the array `orbs`) and
  ! the CMS continuum orbitals (with standing wave or K-matrix normalization,
  ! see eqn. (20) ).
  !
  !   ( Dx )                 ( x )
  !   ( Dy )   = <Psi (k)  | ( y ) | Psi >
  !   ( Dz )         l,m     ( z )      b
  !
  ! The integration is done by numerical quadrature. The grid points and weights
  ! have to be provided as well as the values of the bound orbitals Psi_b at those
  ! grid points. This functions avoids having to keep the wavefunctions of all
  ! angular momentum channels in memory.
  !
  ! The atomic radial wavefunction for region I have to be
  ! provided as B-splines. 
  !
  ! If the bound orbitals are not eigenfunctions of the same muffin tin potential
  ! as the continuum orbitals, the transition dipoles become dependent on the
  ! choice of the origin. Therefore the projection of the bound orbitals onto
  ! the continuum orbitals is calculated to get an idea of the severity of the
  ! non-orthogonality.
  !
  ! Parameters
  ! ----------
  ! energy       : float, energy of continuum orbital (in Hartree)
  ! vconst       : float, constant potential in region II (in Hartree)
  ! chargeIII    : float > 0, asymptotic monopole charge of effective potential,
  !                Veff(r) ----> -chargeIII/r  for r --> oo, typically charge=+1 for
  !                a neutral molecule. This is the charge felt by an electron in region III. 
  ! nat          : int > 0, number of atomic centers
  ! lmax         : int > 0, largest angular momentum, l=0,1,...,lmax
  ! rhoIII       : float, radius of molecular sphere separating region II from region III (in bohr)
  ! radii        : float array of shape (nat), radii of atomic spheres separating region I from region II (in bohr),
  !                radii[i] = rho_i
  ! atpos        : float array of shape (3,nat) with cartesian positions of atomic centers (in bohr)
  ! nspl         : integer, giving the total number of knots for
  !                the splines of the atomic radial wavefunctions
  ! kspl         : integer, giving the degree of the splines
  ! fIknots      : array of shape (nspl,lmax+1,nat), which contains the B-spline knots,
  !                the knots have to be uniformly spaced. fIknots(:,l,i) are the knots belonging
  !                to the radial wavefunction with angular momentum `l` on atom `i`
  ! fIcoeffs     : array of shape (nspl,lmax+1,nat), which contains the B-spline coefficients
  ! nc           : integer, number of angular momentum channels,
  !                number of columns in `sol`
  ! sol          : array of shape (ndim,nc) with solution vectors of matching equations, M.sol = rhs
  ! lm0          : integer, index of first solution vector in `sol` (starting at 0)
  ! npts         : integer, number of points
  ! x_arr        : arrays of shape (npts) with cartesian coordinates of grid points (in bohr)
  ! y_arr
  ! z_arr
  ! w_arr        : array of shape (npts) with weights w_i for numerical integration on the grid
  !                  /
  !                  | f(x,y,z) dV = sum  w  f(x ,y ,z )
  !                  /                  i  i    i  i  i
  ! norb         : integer, number of bound orbitals
  ! orbs_arr     : array of shape (norb, npts) with values of bound orbitals on the grid
  ! 
  ! Returns
  ! -------
  ! kmat         : rows of K-matrix belonging to the solution vectors
  !                in `sol`
  ! tdip_arr     : array of shape (nc,norb,3) with the transition dipole vectors between bound and
  !                continuum orbitals
  ! norms2       : array of shape (norb) with norms**2 of bound orbitals, If the bound orbitals
  !                are normalized to 1, these values can be used to detect if the resolution of
  !                the grid is sufficient or needs to be increased.
  ! projs2       : array of shape (norb) with the norms**2 of the projection of the bound orbitals
  !                onto the continuum orbitals. If bound and continuum orbitals are eigenfunctions
  !                of the same Hamiltonian, the overlaps should be zero.
  !                                                         2
  !                   projs2(b) = sum     | <Psi   | Psi > |
  !                                   l,m       l,m     b
  !
  implicit none
  ! ... input variables ...
  double precision, intent(in) :: energy, vconst, chargeIII
  integer, intent(in) :: nat, lmax
  double precision, intent(in) :: rhoIII, radii(nat), atpos(3,nat)
  integer, intent(in) :: nspl, kspl
  double precision, intent(in) :: fIknots(nspl,0:lmax,nat), fIcoeffs(nspl,0:lmax,nat)
  integer, intent(in) :: nc
  complex*16, intent(in) :: sol((lmax+1)*(lmax+1)*(nat+1),nc)
  integer, intent(in) :: lm0
  integer, intent(in) :: npts
  double precision, intent(in) :: x_arr(npts), y_arr(npts), z_arr(npts), w_arr(npts)
  integer, intent(in) :: norb
  complex*16, intent(in) :: orbs_arr(norb, npts)
  ! ... output variables ...
  complex*16, intent(out) :: kmat((lmax+1)*(lmax+1),nc)
  complex*16, intent(out) :: tdip_arr(nc,norb,3)
  double precision, intent(out) :: norms2(norb)
  double precision, intent(out) :: projs2(norb)
  ! ... local variables ...
  ! pi = 3.14...
  double precision, parameter :: pi = 4.d0*datan(1.d0)
  integer :: nlm, ndim
  ! wave vectors in region I and II
  double precision :: k, kappa
  ! Sommerfeld parameter eta = -Z/k
  double precision :: eta
  ! enumerate atoms
  integer :: i
  ! enumerate angular momenta
  integer :: l, l1
  ! multiindices (l,m)
  integer :: lm, lm1
  ! enumerate solution vectors
  integer :: ic
  ! arrays with l- and m-values
  integer :: l_arr((lmax+1)*(lmax+1)), m_arr((lmax+1)*(lmax+1))
  ! sqrt(k/pi)
  double precision :: sqkpi
  !
  ! coefficients of wavefunctions
  ! -----------------------------
  !   in region III  (eqn. (8))
  !
  !                       III  III           III  III
  ! Psi   (r) = sum    [ A    f  (k r )  +  B    g (k r ) ] Y  (th ,ph )
  !    III         l,m    l,m  l     0       l,m  l    0     l,m  0   0
  !
  ! where (r0,th0,ph0) are the spherical coordinates relative to the origin
  complex*16 :: AIII((lmax+1)*(lmax+1),nc)
  complex*16 :: BIII((lmax+1)*(lmax+1),nc)
  !  in region I  (eqn. (2))
  !                   I    I
  ! Psi (r) = sum    A    f (r ) Y  (th ,ph )
  !    i         l,m  l,m  l  i   l,m  i   i
  !
  ! where (r_i,th_i,ph_i) are the spherical coordinates relative to the atom i
  complex*16 :: AI((lmax+1)*(lmax+1),nat,nc)
  !  in region II (eqn. (7))
  !                    II0                                  nat         IIi
  ! Psi  (r) = sum    A    j (kappa r ) Y  (th ,ph )  +  sum    sum    B    n (kappa r ) Y   (th ,ph )
  !    II         l,m  l,m  l        0   l,m  0   0         i=1    l,m  l,m  l        i   l,m   i   i
  complex*16 :: AII0((lmax+1)*(lmax+1),nc)
  complex*16 :: BII((lmax+1)*(lmax+1),nat,nc)
  ! shape of BII array (nlm,nat)
  integer :: sh(2)
  !
  ! Wronskians
  ! ----------
  ! To find the coefficients A^{I_i}_{l,m} and B^{III}_{l,m} from
  ! the solution of the matching conditions according to eqns.
  ! (18) and (19) we need to store some additional Wronskians
  !
  !                                 i
  ! wronskJF    = [j (kappa*rho ), f (rho )]
  !         l,i     l          i    l    i
  !
  double precision :: wronskJF(0:lmax,nat)
  !
  !                                 III
  ! wronskNG  = [n (kappa*rhoIII), g (kappa*rhoIII)]
  !         l     l                 l
  !        
  !                                 III
  ! wronskNF  = [n (kappa*rhoIII), f (kappa*rhoIII)]
  !         l     l                 l
  !
  double precision :: wronskNG(0:lmax), wronskNF(0:lmax)
  ! radius of atomic sphere i
  double precision :: rho_i
  ! spherical Bessel functions of first kind and second kind
  ! and their r-derivatives
  double precision :: jII(0:lmax), jIIp(0:lmax)
  double precision :: nII(0:lmax), nIIp(0:lmax)
  ! atomic radial functions and their r-derivatives at r=rho_i
  double precision :: fI(0:lmax), fIp(0:lmax)
  ! Coulomb functions and their r-derivatives
  double precision :: fIII(0:lmax), fIIIp(0:lmax)
  double precision :: gIII(0:lmax), gIIIp(0:lmax)
  integer :: ifail
  ! x = k*r
  double precision :: kr
  ! cartesian positions and weight
  double precision :: x,y,z,w
  ! enumerate grid points
  integer :: j
  ! Which region does the point belong to (1=I,2=II or 3=III)
  integer :: region
  ! distance to origin or atomic center i
  double precision :: r0, ri
  ! cartesian coordinates relative to center i
  double precision :: xi,yi,zi
  ! angles of spherical coordinates relative to origin or center i
  double precision :: th0,ph0, thi,phi
  ! spherical harmonics Ylm
  complex*16 :: ysph((lmax+1)*(lmax+1))
  ! values of continuum wavefunctions at (x,y,z)
  complex*16 :: wfn(nc)
  ! values of bound orbitals at (x,y,z)
  complex*16 :: orbs(norb)
  ! product of bound orbitals and continuum wavefunctions
  complex*16 :: bfprod(nc,norb)
  ! enumerate bound orbitals
  integer :: ib
  ! overlap between bound and continuum wavefunctions
  complex*16 :: olapbf(nc,norb)
  
  ! For each l-value there are 2*l+1 quantum numbers m=-l,...,0,...,l.
  ! For l=0,1,...,lmax there are nlm = sum_l=0^lmax 2*l+1 = (lmax+1)^2
  ! angular momentum channels.
  nlm = (lmax+1)**2
  ! dimension of the matching matrix
  ndim = (nat+1)*nlm
  ! E = 1/2 k^2
  k = sqrt(2.d0*energy)
  ! wavevector in region II
  kappa = sqrt(k*k - 2*vconst)
  sqkpi = sqrt(k/pi)

  eta = -chargeIII/k
  
  ! ---------- Wronskians --------------------------------------
  ! compute region I,II and III functions at the matching points
  call coul90_scalar(k*rhoIII, eta, lmax, 0, &
       fIII, gIII, fIIIp, gIIIp, ifail)
  ! We need the regular Coulomb function
  !  III
  ! f  (k*r) = 1/(k*r) F (eta=-Z/k; k*r)
  !  l                  l
  ! and the irregular Coulomb function
  !  III
  ! g  (k*r) = 1/(k*r) G (eta=-Z/k; k*r)
  !  l                  l
  ! but coul90_scalar returns F_l(k*r), G_l(k*r) and its derivatives
  ! with respect to x=(k*r).
  kr = k*rhoIII
  ! Therefore we need to transform the derivatives according to
  !   f(x) = 1/x F(x) ---> f'(x) = -F(x)/x^2 + F'(x)/x
  fIIIp = -fIII/kr**2 + fIIIp/kr
  gIIIp = -gIII/kr**2 + gIIIp/kr
  fIII = fIII/kr
  gIII = gIII/kr

  ! The jII and jIIp values are not needed.
  call coul90_scalar(kappa*rhoIII, 0.d0, lmax, 1, &
       jII, nII, jIIp, nIIp, ifail)

  do l=0,lmax
     wronskNG(l) = nII(l) * k * gIIIp(l) - kappa * nIIp(l) * gIII(l)
     wronskNF(l) = nII(l) * k * fIIIp(l) - kappa * nIIp(l) * fIII(l)
  enddo

  do i=1,nat
     rho_i = radii(i)
     call coul90_scalar(kappa*rho_i, 0.d0, lmax, 1, &
          jII, nII, jIIp, nIIp, ifail)
     do l=0,lmax
        ! evaluate radial wavefunctions
        call splev_uniform(fIknots(1:nspl,l,i),nspl,fIcoeffs(1:nspl,l,i),kspl, &
             rho_i, fI(l), ifail)
        ! and the r-derivative
        call splev_deriv_uniform(fIknots(1:nspl,l,i),nspl,fIcoeffs(1:nspl,l,i),kspl, &
             rho_i, fIp(l), ifail)
        wronskJF(l,i) = jII(l) * fIp(l) - kappa * jIIp(l) * fI(l)
     enddo
  enddo
  
  ! ----- Wavefunction coefficients ------------------------------
  ! Given the matrix of solution vectors `sol` we have to extract
  ! the coefficients for each region
  !  sol  ->  (AI,AII0,BII_i,AIII,BIII)
  ! ordering of angular momentum channels 
  call enumerate_lm(lmax, l_arr, m_arr)

  ! shape of BII array
  sh(1) = nlm
  sh(2) = nat

  AIII(:,:) = 0.d0
  do lm=lm0+1,nc+lm0
     ic = lm-lm0
     ! AIII coefficients according to eqn. (23)
     AIII(lm,ic) = sqkpi
     ! extract coefficients A^{II_0}_{l,m} and B^{II_i}_{l,m}
     ! from solution vector
     AII0(1:nlm,ic) = sol(1:nlm,ic)
     BII(1:nlm,1:nat,ic) = reshape(sol(nlm+1:ndim,ic), sh)
     ! find coefficients A^{I_i}_{l,m} according to eqn. (18)
     do i=1,nat
        rho_i = radii(i)
        do lm1=1,nlm
           l1 = l_arr(lm1)
           AI(lm1,i,ic) = BII(lm1,i,ic) / (kappa * rho_i**2 * wronskJF(l1,i))
        enddo
     enddo
     ! find coefficients B^{III}_{l,m} according to eqn. (19)
     do lm1=1,nlm
        l1 = l_arr(lm1)
        BIII(lm1,ic) = - ( AII0(lm1,ic) / (kappa * rhoIII**2) &
             + wronskNF(l1) * AIII(lm1,ic) ) / wronskNG(l1)
     enddo
  enddo
  
  ! rows of the K-matrix
  kmat = sqkpi * BIII

  norms2 = 0.d0
  olapbf = 0.d0

  ! initialize transition amplitudes to zero
  tdip_arr(:,:,:) = (0.d0, 0.d0)
  
  ! --------------- evaluate wavefunctions -----------------
  do j=1,npts
     x = x_arr(j)
     y = y_arr(j)
     z = z_arr(j)
     
     ! --------------------------------------------------------
     ! First we need to find out which region the point (x,y,z)
     ! belongs to.

     ! distance to origin
     r0 = sqrt(x*x+y*y+z*z)
     if (rhoIII < r0) then
        ! Point lies outside molecular sphere.
        region = 3
     else
        ! By default the point lies in region II unless it lies
        ! in an atomic sphere
        region = 2
        
        atoms: do i=1,nat
           ! distance to atom i
           ri = sqrt((x-atpos(1,i))**2+(y-atpos(2,i))**2+(z-atpos(3,i))**2)
           rho_i = radii(i)
           if (ri < rho_i) then
              ! Point lies in region I and belong to atom i.
              region = 1
              exit atoms
           endif
        enddo atoms
     endif
     !write(*,*) j, "-th point belongs to region ", region
     ! ----------------------------------------------------------

     ! ----------------------------------------------------------
     ! evaluate wavefunction for all channels using the expressions
     ! appropriate for the region where (x,y,z) lies in
     ! 
     wfn(:) = (0.d0, 0.d0)
     
     if (region == 1) then
        !write(*,*) "    and atom ", i
        ! coordinates relative to atom i
        xi = x - atpos(1,i)
        yi = y - atpos(2,i)
        zi = z - atpos(3,i)
        ! convert them to spherical coordinates
        call cart2sph(xi,yi,zi, ri,thi,phi)
        ! spherical harmonics Y_lm
        call sphharm(thi,phi, lmax, ysph)

        ! evaluate radial wavefunctions
        do l=0,lmax
           call splev_uniform(fIknots(1:nspl,l,i),nspl,fIcoeffs(1:nspl,l,i),kspl, &
                ri, fI(l), ifail)
        enddo

        ! compose wavefunction according to eqn. (2)
        do lm=1,nlm
           l = l_arr(lm)
           wfn(:) = wfn(:) + AI(lm,i,:) * fI(l) * ysph(lm)
        enddo
           
     else if (region == 2) then
        ! spherical coordinates relative to origin
        call cart2sph(x,y,z, r0,th0,ph0)
        ! spherical harmonics Y_lm
        call sphharm(th0,ph0, lmax, ysph)
        ! spherical Bessel functions
        call coul90_scalar(kappa*r0, 0.d0, lmax, 1, &
             jII, nII, jIIp, nIIp, ifail)
        ! --- first half of eqn. (7) ---
        do lm=1,nlm
           l = l_arr(lm)
           wfn(:) = wfn(:) + AII0(lm,:) * jII(l) * ysph(lm)
        enddo

        do i=1,nat
           ! coordinates relative to atom i
           xi = x - atpos(1,i)
           yi = y - atpos(2,i)
           zi = z - atpos(3,i)
           ! convert them to spherical coordinates
           call cart2sph(xi,yi,zi, ri,thi,phi)
           ! spherical harmonics Y_lm
           call sphharm(thi,phi, lmax, ysph)
           ! spherical Bessel functions
           call coul90_scalar(kappa*ri, 0.d0, lmax, 1, &
                jII, nII, jIIp, nIIp, ifail)
           ! --- second half of eqn. (7) ----
           do lm=1,nlm
              l = l_arr(lm)
              wfn(:) = wfn(:) + BII(lm,i,:) * nII(l) * ysph(lm)
           enddo
        enddo
     else
        ! region III
        ! spherical coordinates relative to origin
        call cart2sph(x,y,z, r0,th0,ph0)
        ! spherical harmonics Y_lm
        call sphharm(th0,ph0, lmax, ysph)
        ! Coulomb functions
        call coul90_scalar(k*r0, eta, lmax, 0, &
             fIII, gIII, fIIIp, gIIIp, ifail)
        kr = k*r0
        fIII(:) = fIII(:)/kr
        gIII(:) = gIII(:)/kr
        
        ! --- eqn. (8) ---
        do lm=1,nlm
           l = l_arr(lm)
           wfn(:) = wfn(:) + (AIII(lm,:) * fIII(l) + BIII(lm,:) * gIII(l)) * ysph(lm)
        enddo
     endif

     ! ------- integrate transition dipole moments -----
     ! quadrature weight of point j
     w = w_arr(j)
     ! values of bound orbitals at point j
     orbs = orbs_arr(:,j)
     ! outer product between bound and continuum orbitals at that point
     !  Psi_lm(x,y,z)^* Psi_b(x,y,z)
     do ib=1,norb
        bfprod(:,ib) = conjg(wfn(:)) * orbs(ib)
     enddo
     !   x-component    <Psi_lm|x|Psi_b>
     tdip_arr(:,:,1) = tdip_arr(:,:,1) + w * (x * bfprod)
     !   y-component    <Psi_lm|y|Psi_b>
     tdip_arr(:,:,2) = tdip_arr(:,:,2) + w * (y * bfprod)
     !   z-component    <Psi_lm|z|Psi_b>
     tdip_arr(:,:,3) = tdip_arr(:,:,3) + w * (z * bfprod)

     ! integrate probability densities of bound orbitals
     norms2(:) = norms2(:) + w * (abs(orbs)**2)
     ! overlap of bound onto continuum orbitals
     olapbf(:,:) = olapbf(:,:) + w * bfprod
  enddo
  ! norms of projection onto continuum basis
  projs2 = sum(abs(olapbf)**2, 1)
end subroutine transition_dipoles

subroutine transform_waves_in(energy, chargeIII, &
     kmat, wfn, lmax, npts, &
     wfn_in)
  ! 
  ! transform standing waves into wavefunctions with
  ! incoming wave (-) normalization according to eqns. (28),(29),(30) and (41).
  ! The factor sum_m' (D^l_m',m)^* Y^*_{l,m'}(k') in eqn. (40) is not included. 
  !
  ! nlm=(lmax+1)**2 is the number of angular momentum channels
  !
  ! Parameters
  ! ----------
  ! energy      :  float, kinetic energy of electron (in Hartree)
  ! chargeIII   :  float > 0, asymptotic monopole charge of effective potential
  ! kmat        :  complex array of shape (nlm,nlm) with K-matrix
  ! wfn         :  complex array of shape (nlm,npts) with standing waves
  !                evaluated on `npts` grid points, wfn(lm,:) = Psi_{III,lm}(x,y,z)
  ! lmax        :  integer, highest angular momentum
  ! npts        :  integer, number of grid points
  !
  ! Returns
  ! -------
  ! wfn_in      :  complex array of shape (nlm,npts) with incoming wave (-) normalization
  !                evaluated on the grid points, wfn_in(lm,:) = Psi^(-)_{III,lm}(x,y,z)
  !
  implicit none
  ! ... input variables ...
  double precision, intent(in) :: energy
  double precision, intent(in) :: chargeIII
  complex*16, intent(in) :: kmat((lmax+1)*(lmax+1),(lmax+1)*(lmax+1))
  complex*16, intent(in) :: wfn((lmax+1)*(lmax+1),npts)
  integer, intent(in) :: lmax, npts
  ! ... output variables ...
  complex*16, intent(out) :: wfn_in((lmax+1)*(lmax+1),npts)
  ! ... local variables ...
  ! imaginary unit
  complex*16, parameter :: imag = (0, 1)
  ! number of angular momentum channels
  integer :: nlm
  ! enumerate angular momentum channels
  integer :: lm
  ! enumerate angular momenta
  integer :: l
  ! arrays with l- and m-values for each angular momentum channel lm
  integer :: l_arr((lmax+1)*(lmax+1)), m_arr((lmax+1)*(lmax+1))
  ! inverse of transformation matrix T
  complex*16 :: tinv((lmax+1)*(lmax+1),(lmax+1)*(lmax+1))
  ! electron velocity
  double precision :: k
  ! Sommerfeld parameter eta = -chargeIII/k
  double precision :: eta
  ! real and imaginary part of Gamma(1+l+eta)
  double precision gr, gi
  ! Coulomb phase shifts
  double precision :: sigma(0:lmax)
  ! Lapack info
  integer :: ifail
  ! pivot array
  integer :: ipiv((lmax+1)*(lmax+1))
  external zgesv
  
  ! The transformation from K-matrix normalized standing waves to incoming
  ! wave (-) normalization is given by eqn. (28), (29) and (41)
  !     (-)    l   -i*sigma_l              -1
  !  Psi    = i   e            sum [(1+i*K)  ]      Psi
  !     L                       L'            L,L'     L'
  ! which can be expressed in matrix form as
  !     (-)
  !  Psi    = T . Psi
  !                 
  ! Instead of computing the inverse of (1+i*K) we can the Psi^(-) by solving
  ! the linear system of equations
  !   -1      (-)                            -1                       -l'  +i*sigma_l'
  !  T   . Psi    = Psi            with    [T  ]     = (1 + i*K)     i    e
  !                                             L,L'            L,L'
  ! 

  ! compute Coulomb phase shifts
  k = sqrt(2*energy)
  eta = -chargeIII/k
  do l=0,lmax
     call cgama(l+1.d0, eta, 1, gr, gi)
     sigma(l) = atan2(gi, gr)
  enddo

  nlm = (lmax+1)**2
  ! mapping from angular momentum channels lm to (l,m)
  call enumerate_lm(lmax, l_arr, m_arr)
  
  ! build inverse of transformation matrix T 
  tinv = imag * kmat
  do lm=1,nlm
     ! add identity matrix
     tinv(lm,lm) = tinv(lm,lm) + 1.d0
     
     ! l-value belonging to (l,m)
     l = l_arr(lm)
     ! multiply column by phase from Coulomb shift
     tinv(:,lm) = tinv(:,lm) * imag**(-l) * exp(imag*sigma(l))
  enddo

  ! solve
  !     -1    (-)
  !    T   Psi    = Psi
  !

  ! RHS is destroyed, so make a copy 
  wfn_in = wfn
  call zgesv(nlm, npts, tinv, nlm, ipiv, wfn_in, nlm, ifail)
  
end subroutine transform_waves_in

subroutine transform_tdip_in(energy, chargeIII, &
     kmat, tdip, lmax, norb, &
     tdip_in)
  ! 
  ! transform transition dipoles between between continuum and bound orbitals
  ! from standing wave normalization to incoming wave (-) normalization 
  ! according to eqns. (28),(29),(30) and (41). The factor sum_m' (D^l_m',m)^* Y^*_{l,m'}(k')
  ! in eqn. (40) is not included. 
  !
  ! nlm=(lmax+1)**2 is the number of angular momentum channels
  !
  ! Parameters
  ! ----------
  ! energy      :  float, kinetic energy of electron (in Hartree)
  ! chargeIII   :  float > 0, asymptotic monopole charge of effective potential
  ! kmat        :  complex array of shape (nlm,nlm) with K-matrix
  ! tdip        :  complex array of shape (nlm,norb,3) with cartesian transition dipoles
  !                between K-matrix normalized continuum and bound orbitals:
  !                                          (x)  
  !                  tdip(lm,b,:) = <Psi   | (y) | Psi > 
  !                                     lm   (z)      b
  ! lmax        :  integer, highest angular momentum
  ! norb        :  integer, number of bound orbitals Psi_b
  !
  ! Returns
  ! -------
  ! tdip_in     :  complex array of shape (nlm,norb,3) with spherical transition dipoles
  !                between incoming wave (-) normalized continuum and bound orbitals:
  !                                           (-)
  !                  tdip_in(lm,b,mu+2) = <Psi    | r Y     | Psi >      mu=-1,0,1
  !                                           lm       1,mu      b
  !
  implicit none
  ! ... input variables ...
  double precision, intent(in) :: energy
  double precision, intent(in) :: chargeIII
  complex*16, intent(in) :: kmat((lmax+1)*(lmax+1),(lmax+1)*(lmax+1))
  complex*16, intent(in) :: tdip((lmax+1)*(lmax+1),norb,3)
  integer, intent(in) :: lmax, norb
  ! ... output variables ...
  complex*16, intent(out) :: tdip_in((lmax+1)*(lmax+1),norb,3)
  ! ... local variables ...
  ! imaginary unit
  complex*16, parameter :: imag = (0, 1)
  ! pi = 3.14...
  double precision, parameter :: pi = 4.d0*datan(1.d0)
  ! number of angular momentum channels
  integer :: nlm
  ! enumerate angular momentum channels
  integer :: lm
  ! enumerate angular momenta
  integer :: l
  ! arrays with l- and m-values for each angular momentum channel lm
  integer :: l_arr((lmax+1)*(lmax+1)), m_arr((lmax+1)*(lmax+1))
  ! conjugate of inverse of transformation matrix, [T^(-1)]^*
  complex*16 :: tinvc((lmax+1)*(lmax+1),(lmax+1)*(lmax+1))
  ! enumerate cartesian axes
  integer :: xyz
  ! electron velocity
  double precision :: k
  ! Sommerfeld parameter eta = -chargeIII/k
  double precision :: eta
  ! real and imaginary part of Gamma(1+l+eta)
  double precision gr, gi
  ! Coulomb phase shifts
  double precision :: sigma(0:lmax)
  ! Lapack info
  integer :: ifail
  ! pivot array
  integer :: ipiv((lmax+1)*(lmax+1))
  external zgesv
  ! x,y and z components of transition dipoles
  complex*16 :: tdipx((lmax+1)*(lmax+1),norb)
  complex*16 :: tdipy((lmax+1)*(lmax+1),norb)
  complex*16 :: tdipz((lmax+1)*(lmax+1),norb)
  double precision :: fac
  
  ! The transformation from K-matrix normalized standing waves to incoming
  ! wave (-) normalization is given by eqn. (28), (29) and (41)
  !     (-)    l   -i*sigma_l              -1
  !  Psi    = i   e            sum [(1+i*K)  ]      Psi
  !     L                       L'            L,L'     L'
  ! which can be expressed in matrix form as
  !     (-)
  !  Psi    = T . Psi
  !                 
  ! Instead of computing the inverse of (1+i*K) we can obtain the Psi^(-) by solving
  ! the linear system of equations
  !   -1      (-)                            -1                       -l'  +i*sigma_l'
  !  T   . Psi    = Psi            with    [T  ]     = (1 + i*K)     i    e
  !                                             L,L'            L,L'
  !
  ! The transition dipole moments transform as
  !                               *       l'  -i*sigma_l'      (-)
  ! <Psi |r| Psi >  = sum (1 - i*K )     i   e             <Psi   | r |Psi >
  !     L       b      L'           L,L'                       L'         b
  !
  !                         -1 *        (-)
  !                 = sum [T  ]     <Psi   | r |Psi >
  !                    L'      L,L'     L'         b
  !
  
  ! compute Coulomb phase shifts
  k = sqrt(2*energy)
  eta = -chargeIII/k
  do l=0,lmax
     call cgama(l+1.d0, eta, 1, gr, gi)
     sigma(l) = atan2(gi, gr)
  enddo

  nlm = (lmax+1)**2
  ! mapping from angular momentum channels lm to (l,m)
  call enumerate_lm(lmax, l_arr, m_arr)
  
  ! build conjugate of inverse of transformation matrix T
  !   -1 *              *       l'  -i*sigma_l'
  ! [T  ]     = (1 - i*K )     i   e
  !      L,L'             L,L'
  tinvc = -imag * conjg(kmat)
  do lm=1,nlm
     ! add identity matrix
     tinvc(lm,lm) = tinvc(lm,lm) + 1.d0
     
     ! l-value belonging to (l,m)
     l = l_arr(lm)
     ! multiply column by phase from Coulomb shift
     tinvc(:,lm) = tinvc(:,lm) * imag**(l) * exp(-imag*sigma(l))
  enddo

  ! solve
  !     -1*    
  !    T    tdip    = tdip
  !             in

  ! RHS is destroyed, so make a copy 
  tdip_in = tdip
  do xyz=1,3
     ! transform x-,y- and z-component
     call zgesv(nlm, norb, tinvc, nlm, ipiv, tdip_in(:,:,xyz), nlm, ifail)
  enddo

  ! The cartesian transition dipole moments <x>, <y>, <z> have to be
  ! transformed to <r*Y_{1,-1}>, <r*Y_{1,0}> and <r*Y_{1,+1}> according to
  !       (-)                                          (-)                     (-)
  !    <Psi  | r Y    |Psi > =  1/2 sqrt(3/(2pi)) { <Psi  | x | Psi > - i * <Psi  | y | Psi > }
  !        L      1,-1    b                             L          b            L          b
  !
  !       (-)                                      (-)            
  !    <Psi  | r Y   |Psi >  = 1/2 sqrt(3/pi)  <Psi  | z | Psi > 
  !        L      1,0    b                         L          b  
  !
  !       (-)                                          (-)                     (-)
  !    <Psi  | r Y    |Psi > = -1/2 sqrt(3/(2pi)) { <Psi  | x | Psi > + i * <Psi  | y | Psi > }
  !        L      1,+1    b                             L          b            L          b

  ! cartesian transition dipoles
  tdipx = tdip_in(:,:,1)
  tdipy = tdip_in(:,:,2)
  tdipz = tdip_in(:,:,3)
  fac = 0.5d0 * sqrt(3.d0/(2.d0*pi))
  ! <r*Y_{1,-1}>
  tdip_in(:,:,1) =              fac * (tdipx - imag*tdipy)
  ! <r*Y_{1,0}>
  tdip_in(:,:,2) = sqrt(2.d0) * fac * tdipz
  ! <r*Y_{1,+1}>
  tdip_in(:,:,3) =             -fac * (tdipx + imag*tdipy)
end subroutine transform_tdip_in


subroutine cgama ( x, y, kf, gr, gi )

!*****************************************************************************80
!
!! CGAMA computes the Gamma function for complex argument.
!
!  Discussion:
!
!    This procedcure computes the gamma function g(z) or ln[g(z)]
!    for a complex argument
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    26 July 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, Y, the real and imaginary parts of 
!    the argument Z.
!
!    Input, integer ( kind = 4 ) KF, the function code.
!    0 for ln[g(z)]
!    1 for g(z)
!
!    Output, real ( kind = 8 ) GR, GI, the real and imaginary parts of
!    the selected function.
!
  implicit none

  real ( kind = 8 ), save, dimension ( 10 ) :: a = (/ &
    8.333333333333333D-02, -2.777777777777778D-03, &
    7.936507936507937D-04, -5.952380952380952D-04, &
    8.417508417508418D-04, -1.917526917526918D-03, &
    6.410256410256410D-03, -2.955065359477124D-02, &
    1.796443723688307D-01, -1.39243221690590D+00 /)
  real ( kind = 8 ) g0
  real ( kind = 8 ), intent(out) :: gi
  real ( kind = 8 ) gi1
  real ( kind = 8 ), intent(out) :: gr
  real ( kind = 8 ) gr1
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kf
  integer ( kind = 4 ) na
  real ( kind = 8 ) pi
  real ( kind = 8 ) si
  real ( kind = 8 ) sr
  real ( kind = 8 ) t
  real ( kind = 8 ) th
  real ( kind = 8 ) th1
  real ( kind = 8 ) th2
  real ( kind = 8 ) x
  real ( kind = 8 ) x0
  real ( kind = 8 ) x1
  real ( kind = 8 ) y
  real ( kind = 8 ) y1
  real ( kind = 8 ) z1
  real ( kind = 8 ) z2

  pi = 3.141592653589793D+00

  if ( y == 0.0D+00 .and. x == int ( x ) .and. x <= 0.0D+00 ) then
    gr = 1.0D+300
    gi = 0.0D+00
    return
  else if ( x < 0.0D+00 ) then
    x1 = x
    y1 = y
    x = -x
    y = -y
  end if

  x0 = x

  if ( x <= 7.0D+00 ) then
    na = int ( 7 - x )
    x0 = x + na
  end if

  z1 = sqrt ( x0 * x0 + y * y )
  th = atan ( y / x0 )
  gr = ( x0 - 0.5D+00 ) * log ( z1 ) - th * y - x0 &
    + 0.5D+00 * log ( 2.0D+00 * pi )
  gi = th * ( x0 - 0.5D+00 ) + y * log ( z1 ) - y

  do k = 1, 10
    t = z1 ** ( 1 - 2 * k )
    gr = gr + a(k) * t * cos ( ( 2.0D+00 * k - 1.0D+00 ) * th )
    gi = gi - a(k) * t * sin ( ( 2.0D+00 * k - 1.0D+00 ) * th )
  end do

  if ( x <= 7.0D+00 ) then
    gr1 = 0.0D+00
    gi1 = 0.0D+00
    do j = 0, na - 1
      gr1 = gr1 + 0.5D+00 * log ( ( x + j ) ** 2 + y * y )
      gi1 = gi1 + atan ( y / ( x + j ) )
    end do
    gr = gr - gr1
    gi = gi - gi1
  end if

  if ( x1 < 0.0D+00 ) then
    z1 = sqrt ( x * x + y * y )
    th1 = atan ( y / x )
    sr = - sin ( pi * x ) * cosh ( pi * y )
    si = - cos ( pi * x ) * sinh ( pi * y )
    z2 = sqrt ( sr * sr + si * si )
    th2 = atan ( si / sr )
    if ( sr < 0.0D+00 ) then
      th2 = pi + th2
    end if
    gr = log ( pi / ( z1 * z2 ) ) - gr
    gi = - th1 - th2 - gi
    x = x1
    y = y1
  end if

  if ( kf == 1 ) then
    g0 = exp ( gr )
    gr = g0 * cos ( gi )
    gi = g0 * sin ( gi )
  end if

  return
end subroutine cgama

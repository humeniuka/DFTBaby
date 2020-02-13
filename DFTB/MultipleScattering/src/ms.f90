!
! Multiple-Scattering (MS) Method for bound orbitals
!
! Implementation of Johnson's multiple scattering method according
! to Ref. [1]. Numbers of equations in the source code refer to his 
! publication.
!
! References
! ----------
! [1] K.H. Johnson, "Scattered-Wave Theory of the Chemical Bond" (1973),
!     Advances in quantum chemistry, volume 7, pages 143-185, Elsevier
!

subroutine jn_coeffs(l1,m1, l2,m2, vec, jcoeff, ncoeff)
  !
  ! coefficients for expanding spherical Bessel functions of
  ! first and second kind around different center. The new center
  ! differs j from the old center i by the shift vector R_ij=R_j-R_i.
  !
  ! The expansion theorems for modified spherical Bessel functions of
  ! the first (j_l) and second (n_l) kind are:
  !
  !   j (r-R) Y (r-R) = sum  J (R)  j (r) Y (r)
  !    l       L         L'   L,L'   l'    L'
  !
  !                         
  !   n (r-R) Y (r-R) = sum  N (R)  j (r) Y (r)    for  r < R
  !    l       L         L'   L,L'   l'    L'
  !
  !  and
  !                         
  !   n (r-R) Y (r-R) = sum  J (R)  n (r) Y (r)    for r > R
  !    l       L         L'   L,L'   l'    L'
  !
  ! The coefficients J_{L,L'}(R) and N_{L,L'}(R) are defined in
  ! eqn. 14 in Dill & Dehmer (1974).
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
  ! ncoeff        
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

subroutine ik_coeffs(l1,m1, l2,m2, vec, icoeff, kcoeff)
  !
  ! coefficients for expanding modified spherical Bessel functions of
  ! first (i_l) and second (k_l) kind around a different center. The new center
  ! differs j from the old center i by the shift vector R_ij=R_j-R_i.
  ! The expansion theorems (II.23),(II.24) and (II.25) can be written as:
  !
  !   i (r-R) Y (r-R) = sum  I (R)  i (r) Y (r)
  !    l       L         L'   L,L'   l'    L'
  !
  !                         
  !   k (r-R) Y (r-R) = sum  K (R)  i (r) Y (r)    for  r < R
  !    l       L         L'   L,L'   l'    L'
  !
  !  and
  !                         
  !   k (r-R) Y (r-R) = sum  I (R)  k (r) Y (r)    for r > R
  !    l       L         L'   L,L'   l'    L'
  !
  ! Parameters
  ! ----------
  ! l1,m1      : angular momentum quantum numbers, 0 <= l1, -l1 <= m1 <= +l1
  ! l2,m2      : angular momentum quantum numbers, 0 <= l2, -l2 <= m1 <= +l2
  ! vec        : cartesian shift vector R_ij, vec = [x,y,z]
  !
  ! Returns
  ! -------
  ! icoeff,     : complex, coefficients I(R)    and K(R)
  ! kcoeff                               L1,L2       L1,L2
  !
  implicit none
  ! ... input variables ...
  integer, intent(in) :: l1,m1, l2,m2
  double precision, intent(in) :: vec(3)
  ! ... output variables ...
  complex*16, intent(out) :: icoeff, kcoeff
  ! ... local variables ...
  ! pi = 3.14...
  double precision, parameter :: pi = 4.d0*datan(1.d0)
  ! spherical harmonics Y_{l,m}
  double precision :: r, th, ph
  integer :: l, m, lmin, lmax, lrng
  ! offsets where Wigner 3j symbols start
  integer :: jmin0, jmax0, jminm, jmaxm
  complex*16 :: y(l1+l2+1)
  double precision :: eta
  ! modified spherical Bessel functions i_l(r) and k_l(r)
  double precision :: bessel_i(l1+l2+1), bessel_k(l1+l2+1)
  ! r-derivatives of i_l(r) and k_l(r)
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
  
  ! compute modified spherical Bessel functions of 1st and 2nd kind
  lrng = lmax-lmin
  call mod_bessel(r, lmax, bessel_i, bessel_k, fcp, gcp)
  
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

  icoeff = (0.d0, 0.d0)
  kcoeff = (0.d0, 0.d0)
  do l=lmin,lmax
     c = (-1)**(l+m1) &
          * sqrt(4.d0*pi*(2*l+1)*(2*l1+1)*(2*l2+1)) &
          * w3j0(l-jmin0+1) * w3jm(l-jminm+1) &
          * y(l+1)
     icoeff = icoeff + c * bessel_i(l+1)
     kcoeff = kcoeff + c * bessel_k(l+1)
  enddo
end subroutine ik_coeffs

subroutine jn_or_ik_coeffs(l1,m1, l2,m2, vec, energy, vconst, &
     ji_coeff, nk_coeff)
  !
  ! coefficients for expanding (modified) spherical Bessel functions of
  ! first (j_l or i_l) and second (n_l or k_l) kind around a different center.
  !
  ! Depending on whether the energy is larger than the constant potential
  ! in region II or not the expansion coefficients for j_l and n_l or those
  ! for i_l and k_l are calculated.
  implicit none
  ! ... input variables ...
  integer, intent(in) :: l1,m1, l2,m2
  double precision, intent(in) :: vec(3)
  double precision, intent(in) :: energy, vconst
  ! ... output variables ...
  complex*16, intent(out) :: ji_coeff, nk_coeff

  if (energy >= vconst) then
     ! E >= V_II
     call jn_coeffs(l1,m1, l2,m2, vec, ji_coeff, nk_coeff)
  else
     ! E < V_II
     call ik_coeffs(l1,m1, l2,m2, vec, ji_coeff, nk_coeff)
  endif
end subroutine jn_or_ik_coeffs

subroutine matching_matrix( &
     energy, vconst, chargeIII, nat, lmax, &
     rhoIII, radii, atpos, &
     nsplI, nsplIII, kspl, fIknots, fIcoeffs, fIIIknots, fIIIcoeffs, &
     debug, &
     matchM )
  !
  ! construct the matrix for matching region I,II and III wavefunctions
  ! and their derivatives at the boundaries continuously. The coefficients
  ! {C^III_L} and {C^Ii_L} are eliminated from the system of equations.
  !
  ! The coefficients `x` of the MS wavefunction in the region II
  !
  !  {A^II0_L} and {A^IIi_L}
  !
  ! are determined by solving a homogeneous system of linear equations
  !
  !      M.x = 0
  !
  ! `x` are in fact the eigenvectors belonging to eigenvalue 0. So a necessary
  ! condition for a bound state is det(M) = 0.
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
  ! nsplI        : integer, giving the total number of knots for
  !                the splines of the atomic radial wavefunctions in region I
  ! nsplIII      : integer, giving the total number of knots for
  !                the splines of the radial wavefunctions in region III
  ! kspl         : integer, giving the degree of the splines
  ! fIknots      : array of shape (nsplI,lmax+1,nat), which contains the B-spline knots,
  !                the knots have to be uniformly spaced. fIknots(:,l,i) are the knots belonging
  !                to the radial wavefunction with angular momentum `l` on atom `i`
  ! fIcoeffs     : array of shape (nsplI,lmax+1,nat), which contains the B-spline coefficients
  ! fIIIknots    : array of shape (nsplIII,lmax+1), which contains the B-spline knots,
  !                the knots have to be uniformly spaced. fIIIknots(:,l) are the knots belonging
  !                to the radial wavefunction with angular momentum `l` in region III
  ! fIIIcoeffs   : array of shape (nsplIII,lmax+1), which contains the B-spline coefficients for the
  !                region III wavefunctions
  !
  ! Optional
  ! --------
  ! debug        : > 0, print additional information
  !
  ! Returns
  ! -------
  ! matchM       : complex array of shape (ndim, ndim), matching matrix
  !
  implicit none
  ! ... input variables ...
  double precision, intent(in) :: energy, vconst, chargeIII
  integer, intent(in) :: nat, lmax
  double precision, intent(in) :: rhoIII, radii(nat), atpos(3,nat)
  integer, intent(in) :: nsplI, nsplIII, kspl
  double precision, intent(in) :: fIknots(nsplI,0:lmax,nat), fIcoeffs(nsplI,0:lmax,nat)
  double precision, intent(in) :: fIIIknots(nsplIII,0:lmax), fIIIcoeffs(nsplIII,0:lmax)
  ! ... optional ...
  !f2py integer optional, intent(in) :: debug = 0
  integer, intent(in) :: debug
  ! ... output variables ...
  complex*16, intent(out) :: matchM( (lmax+1)*(lmax+1)*(nat+1), (lmax+1)*(lmax+1)*(nat+1) )
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
  ! fIII         : array of shape (lmax+1)
  !                fIII(l) is the value of the radial wavefunction of region III with angular
  !                momentum l at the border of region III,
  !                               (III)
  !                    fIII[l] = f   (r=rhoIII)
  !                               l        
  ! fIIIp        : float array of shape (lmax+1)
  !                fIIIp(l) is the r-derivative of the radial wavefunction of region III evaluated at
  !                the border of region III
  !
  !                               d f^(II)_l
  !                    fIIIp[l] = --------- (r = rhoIII)
  !                                 d r            
  double precision :: fIII(0:lmax), fIIIp(0:lmax)
  ! radius of an atomic sphere
  double precision :: rho_i
  integer :: nlm, ndim
  ! wave vectors in region I and II
  double precision :: k, kappa
  ! Sommerfeld parameter eta = -Z/k
  double precision :: eta
  ! Wronskians
  double precision :: wr(0:lmax,0:nat)
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
  ! modified spherical Bessel functions of 1st and 2nd kind
  double precision :: iII(0:lmax), kII(0:lmax)
  ! and their r-derivatives at a certain radius
  double precision :: iIIp(0:lmax), kIIp(0:lmax)
  ! abbreviation for k*rhoIII
  double precision :: x
  ! arrays with l- and m-values
  integer :: l_arr((lmax+1)*(lmax+1)), m_arr((lmax+1)*(lmax+1))
  ! difference vector of atom positions
  double precision :: vecij(3)
  ! coefficients for reexpansion of spherical Bessel functions around
  ! different center
  complex*16 :: jcoeff, ncoeff
  ! coefficients for reexpansion of modified spherical Bessel functions around
  ! different center
  complex*16 :: icoeff, kcoeff
  integer :: ifail

  if (energy >= 0.d0) then
     write(*,*) "Energy for bound states has to be < 0! Got ", energy
     stop
  endif
  
  ! evaluate atomic radial wavefunctions at the radii of the atomic
  ! spheres. The radial wavefunctions are given as B-splines. 
  do i=1,nat
     rho_i = radii(i)
     do l=0,lmax
        ! evaluate radial wavefunctions
        call splev_uniform(fIknots(1:nsplI,l,i),nsplI,fIcoeffs(1:nsplI,l,i), &
             kspl, &
             rho_i, fI(l,i), ifail)
        ! and the r-derivative
        call splev_deriv_uniform(fIknots(1:nsplI,l,i),nsplI,fIcoeffs(1:nsplI,l,i), &
             kspl, &
             rho_i, fIp(l,i), ifail)
     enddo
  enddo

  ! evaluate radial wavefunctions of region III at the radius of the molecular sphere (r=rhoIII)
  ! The radial wavefunctions are given as B-splines. 
  do l=0,lmax
     ! evaluate radial wavefunctions
     call splev_uniform(fIIIknots(1:nsplIII,l),nsplIII,fIIIcoeffs(1:nsplIII,l), &
          kspl, &
          rhoIII, fIII(l), ifail)
     ! and the r-derivative
     call splev_deriv_uniform(fIIIknots(1:nsplIII,l),nsplIII,fIIIcoeffs(1:nsplIII,l), &
          kspl, &
          rhoIII, fIIIp(l), ifail)
  enddo

  ! For each l-value there are 2*l+1 quantum numbers m=-l,...,0,...,l.
  ! For l=0,1,...,lmax there are nlm = sum_l=0^lmax 2*l+1 = (lmax+1)^2
  ! angular momentum channels.
  nlm = (lmax+1)**2
  ! dimension of the matching matrix
  ndim = (nat+1)*nlm

  ! The wavefunctions in region II depend on whether the energy E is
  ! larger or smaller than the constant potential V_II.
  ! For E >= V_II the wavefunctions in region II are spherical Bessel functions
  !    j_l(r) and n_l(r),
  ! while for E < V_II they are modified spherical Bessel functions
  !    i_l(r) and k_l(r)
  !
  ! The expressions for E >= V_II can be obtained by those for E < V_II by the
  ! following replacements
  !
  !    E < V_I              E >= V_II
  !   ---------------------------------
  !   I_{L',L}(R)  <---->  J_{L',L}(R)
  !   K_{L',L}(R)  <---->  N_{L',L}(R)
  !   i_l(r)       <---->  j_l(r)
  !   k_l(r)       <---->  n_l(r)
  !
  !write(*,*) "energy= ", energy, " V_II= ", vconst
  if (energy >= vconst) then
     !write(*,*) " E >= V_II "
     kappa = sqrt(2.d0*(energy - vconst))
  else
     !write(*,*) " E < V_II "
     kappa = sqrt(2.d0*(vconst - energy))
  endif
  !write(*,*) "length of wavevector in region II   kappa = ", kappa
  
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
  ! If the energy exceeds the constant potential in region II, i.e. E >= vconst,
  ! 
  !
  !  +     [n_l(kappa*rho_i), f_l^i(rho_i)]      
  ! w    = --------------------------------
  !  l,i   [j_l(kappa*rho_i), f_l^i(rho_i)]
  !
  ! otherwise, i.e. E < vconst,
  !
  !  -     [k_l(kappa*rho_i), f_l^i(rho_i)]      
  ! w    = --------------------------------
  !  l,i   [i_l(kappa*rho_i), f_l^i(rho_i)]
  !
  ! i=0 corresponds to the origin, so f_l^0(r) = f_l^III(r).
  !
  if (energy >= vconst) then
     ! E >= V_II, we have to compute w_{l,i}^+
     
     ! i=0, origin 
     ! compute spherical Bessel (kfn=1) functions j_l and n_l for l=0,...,lmax
     call coul90_scalar(kappa*rhoIII, 0.d0, lmax, 1, &
          jII, nII, jIIp, nIIp, ifail)
     do l=0,lmax
        ! the factor kappa comes from the chain rule
        ! since the derivatives are with respect to kappa*r
        wronsk_numer = nII(l) * fIIIp(l) - kappa*nIIp(l) * fIII(l)
        wronsk_denom = jII(l) * fIIIp(l) - kappa*jIIp(l) * fIII(l)
        ! w_{l,0}^+
        wr(l,0) = wronsk_numer / wronsk_denom
     enddo
     ! i > 0, atomic centers
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
           ! w_{l,i}^+
           wr(l,i) = wronsk_numer / wronsk_denom
        enddo
     enddo
  else
     ! E < V_II, we have to compute w_{l,i}^-

     ! i=0, origin 
     ! compute modified spherical Bessel functions i_l and k_l for l=0,...,lmax
     call mod_bessel(kappa*rhoIII, lmax, &
          iII, kII, iIIp, kIIp, ifail)
     do l=0,lmax
        ! the factor kappa comes from the chain rule
        ! since the derivatives are with respect to kappa*r
        wronsk_numer = kII(l) * fIIIp(l) - kappa*kIIp(l) * fIII(l)
        wronsk_denom = iII(l) * fIIIp(l) - kappa*iIIp(l) * fIII(l)
        ! w_{l,0}^-
        wr(l,0) = wronsk_numer / wronsk_denom
     enddo
     ! i > 0, atomic centers
     do i=1,nat
        ! radius of atomic sphere i
        rho_i = radii(i)
        ! compute modified spherical Bessel functions i_l and k_l for l=0,...,lmax
        call mod_bessel(kappa*rho_i, lmax, &
            iII, kII, iIIp, kIIp, ifail)
        do l=0,lmax
           ! the factor kappa comes from the chain rule
           ! since the derivatives are with respect to kappa*r
           wronsk_numer = kII(l) * fIp(l,i) - kappa*kIIp(l) * fI(l,i)
           wronsk_denom = iII(l) * fIp(l,i) - kappa*iIIp(l) * fI(l,i)
           ! w_{l,i}^-
           wr(l,i) = wronsk_numer / wronsk_denom
        enddo
     enddo
  endif
  
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
                    matchM(ii,jj) = 1.d0/wr(li,0)
                 else
                    ! (i == j, i > 0) diagonal block, atomic centers
                    matchM(ii,jj) = wr(li,i)
                 endif
              else
                 ! off-diagonal elements  i>1 , j>1, i!=j
                 !
                 ! Here we need to distinguish the cases E >= V_II or E < V_II,
                 ! since the expansion theorems for j_l(r), n_l(r) and i_l(r), k_l(r)
                 ! require different coefficients J,N or I,K.
                 !
                 ! NOTE: The additional minus signs in the arguments of I, J, K and N
                 ! for E > V_II were added later and are needed for smooth matching,
                 ! although I don't know why. 
                 ! 
                 if (energy >= vconst) then
                    if ((i == 0).and.(j > 0)) then
                       ! J     (kappa R  )    with  R  = R  - R  = R
                       !  L',L         0j            0j   j    0    j
                       call jn_coeffs(lj,mj, li,mi, -kappa*atpos(:,j), jcoeff, ncoeff)
                       matchM(ii,jj) = jcoeff
                    else if ((j == 0).and.(i > 0)) then
                       ! J     (kappa R  )    with R   = R  - R  = -R
                       !  L',L         i0           i0    0    i     i
                       call jn_coeffs(lj,mj, li,mi, +kappa*atpos(:,i), jcoeff, ncoeff)
                       matchM(ii,jj) = jcoeff
                    else if ((i > 0).and.(j > 0).and.(i /= j)) then
                       ! N    (kappa R  )     with R   = R  - R
                       !  L',L        ij            ij    j    i
                       ! vector from atom i to atom j
                       vecij = atpos(:,j) - atpos(:,i)
                       call jn_coeffs(lj,mj, li,mi, -kappa*vecij, jcoeff, ncoeff)
                       matchM(ii,jj) = ncoeff
                    endif
                 else
                    ! E < V_II
                    if ((i == 0).and.(j > 0)) then
                       ! I     (kappa R  )    with  R  = R  - R  = R
                       !  L',L         0j            0j   j    0    j
                       call ik_coeffs(lj,mj, li,mi, kappa*atpos(:,j), icoeff, kcoeff)
                       matchM(ii,jj) = icoeff
                    else if ((j == 0).and.(i > 0)) then
                       ! I     (kappa R  )    with R   = R  - R  = -R
                       !  L',L         i0           i0    0    i     i
                       call ik_coeffs(lj,mj, li,mi, -kappa*atpos(:,i), icoeff, kcoeff)
                       matchM(ii,jj) = icoeff
                    else if ((i > 0).and.(j > 0).and.(i /= j)) then
                       ! K    (kappa R  )     with R   = R  - R
                       !  L',L        ij            ij    j    i
                       ! vector from atom i to atom j
                       vecij = atpos(:,j) - atpos(:,i)
                       call ik_coeffs(lj,mj, li,mi, kappa*vecij, icoeff, kcoeff)
                       matchM(ii,jj) = kcoeff
                    endif                    
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
  
end subroutine matching_matrix

subroutine normalization( &
     energy, vconst, chargeIII, nat, lmax, &
     rhoIII, radii, atpos, &
     nsplI, nsplIII, kspl, fIknots, fIcoeffs, fIIIknots, fIIIcoeffs, nr, &
     debug, &
     nrm )
  !
  ! construct the coefficients for normalizing the MS wavefunctions.
  !
  ! Given the coefficients `x` of the MS wavefunction in region II
  ! the norm of the wavefunction in regions I and III is given by
  !        I 2        III 2                    2
  !    |Psi |  +  |Psi   |   =  sum  |nrm  x |
  !                                i     i  i
  ! where nrm(:) is a vector of normalization coefficients.
  !
  ! The part of the wavefunction in region II is excluded, because
  ! the integrals are too complicated.
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
  ! nsplI        : integer, giving the total number of knots for
  !                the splines of the atomic radial wavefunctions in region I
  ! nsplIII      : integer, giving the total number of knots for
  !                the splines of the radial wavefunctions in region III
  ! kspl         : integer, giving the degree of the splines
  ! fIknots      : array of shape (nsplI,lmax+1,nat), which contains the B-spline knots,
  !                the knots have to be uniformly spaced. fIknots(:,l,i) are the knots belonging
  !                to the radial wavefunction with angular momentum `l` on atom `i`
  ! fIcoeffs     : array of shape (nsplI,lmax+1,nat), which contains the B-spline coefficients
  ! fIIIknots    : array of shape (nsplIII,lmax+1), which contains the B-spline knots,
  !                the knots have to be uniformly spaced. fIIIknots(:,l) are the knots belonging
  !                to the radial wavefunction with angular momentum `l` in region III
  ! fIIIcoeffs   : array of shape (nsplIII,lmax+1), which contains the B-spline coefficients for the
  !                region III wavefunctions
  ! nr           : int, number of grid points for radial integration using Gauss-Chebyshev quadrature
  !                of the second kind.
  !
  ! Optional
  ! --------
  ! debug        : > 0, print additional information
  !
  ! Returns
  ! -------
  ! nrm          : array of shape (ndim) with normalization coefficients
  !
  implicit none
  ! ... input variables ...
  double precision, intent(in) :: energy, vconst, chargeIII
  integer, intent(in) :: nat, lmax
  double precision, intent(in) :: rhoIII, radii(nat), atpos(3,nat)
  integer, intent(in) :: nsplI, nsplIII, kspl
  double precision, intent(in) :: fIknots(nsplI,0:lmax,nat), fIcoeffs(nsplI,0:lmax,nat)
  double precision, intent(in) :: fIIIknots(nsplIII,0:lmax), fIIIcoeffs(nsplIII,0:lmax)
  integer, intent(in) :: nr
  ! ... optional ...
  !f2py integer optional, intent(in) :: debug = 0
  integer, intent(in) :: debug
  ! ... output variables ...
  double precision, intent(out) :: nrm( (lmax+1)*(lmax+1)*(nat+1) )
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
  ! fIII         : array of shape (lmax+1)
  !                fIII(l) is the value of the radial wavefunction of region III with angular
  !                momentum l at the border of region III,
  !                               (III)
  !                    fIII[l] = f   (r=rhoIII)
  !                               l        
  ! fIIIp        : float array of shape (lmax+1)
  !                fIIIp(l) is the r-derivative of the radial wavefunction of region III evaluated at
  !                the border of region III
  !
  !                               d f^(II)_l
  !                    fIIIp[l] = --------- (r = rhoIII)
  !                                 d r            
  double precision :: fIII(0:lmax), fIIIp(0:lmax)
  ! radius of an atomic sphere
  double precision :: rho_i
  integer :: nlm, ndim
  ! wave vectors in region I and II
  double precision :: kappa
  ! enumerates atoms
  integer :: i,j
  ! enumerate rows and columns
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
  ! modified spherical Bessel functions of 1st and 2nd kind
  double precision :: iII(0:lmax), kII(0:lmax)
  ! and their r-derivatives at a certain radius
  double precision :: iIIp(0:lmax), kIIp(0:lmax)
  ! arrays with l- and m-values
  integer :: l_arr((lmax+1)*(lmax+1)), m_arr((lmax+1)*(lmax+1))
  integer :: ifail
  ! radial integration, nodes and weights
  double precision :: zr(nr), xr(nr), wr(nr)
  ! position of quadrature points and integration limits
  double precision :: x, r, rmin, rmax
  ! volume element, dr/dx, dx and r^2 * dr
  double precision :: drdx, dx, r2dr
  ! enumerate quadrature points
  integer :: n
  
  ! radial integrals of atomic wavefunctions in region I
  double precision :: radint_If(0:lmax,nat)
  !
  !                  /rho_i   2     (i)   2
  ! radint_If(l,i) = |       r   ( f (r) )  dr
  !                  /0             l
  !
  !

  ! radial integrals of wavefunction in region III (rhoIII < r)
  double precision :: radint_IIIf(0:lmax)
  !
  !                  /oo      2     (III)  2
  ! radint_IIIf(l) = |       r   ( f  (r) )  dr
  !                  /rhoIII        l
  !

  ! Wronskians
  double precision :: wronsk_numer, wronsk_denom
  ! Wronskians relate the coefficients in regions I and III to those in region II
  double precision :: wrII(0:lmax,nat)
  !
  !  II     [k_l(kappa*rho_i), i_l(kappa*rho_i)]
  ! w     = ------------------------------------           for E < V_II
  !  l,i      [f_l^i(rho_i), i_l(kappa*rho_i)]
  !
  !  or 
  !
  !  II     [n_l(kappa*rho_i), j_l(kappa*rho_i)]
  ! w     = ------------------------------------           for E >= V_II
  !  l,i      [f_l^i(rho_i), j_l(kappa*rho_i)]
  !
  double precision :: wrIII(0:lmax)
  !
  !  III     [i_l(kappa*rhoIII), k_l(kappa*rhoIII)]
  ! w      = --------------------------------------        for E < V_II
  !   l       [f_l^III(rhoIII), k_l(kappa*rhoIII)]
  !
  !  or
  !
  !  III     [j_l(kappa*rhoIII), n_l(kappa*rhoIII)]
  ! w      = --------------------------------------        for E >= V_II
  !   l       [f_l^III(rhoIII), n_l(kappa*rhoIII)]
  !
  double precision :: nrm2
  
  if (energy >= 0.d0) then
     write(*,*) "Energy for bound states has to be < 0! Got ", energy
     stop
  endif

  ! The wavefunctions in region II depend on whether the energy E is
  ! larger or smaller than the constant potential V_II.
  ! For E >= V_II the wavefunctions in region II are spherical Bessel functions
  !    j_l(r) and n_l(r),
  ! while for E < V_II they are modified spherical Bessel functions
  !    i_l(r) and k_l(r)
  !
  ! The expressions for E >= V_II can be obtained by those for E < V_II by the
  ! following replacements
  !
  !    E < V_I              E >= V_II
  !   ---------------------------------
  !   I_{L',L}(R)  <---->  J_{L',L}(R)
  !   K_{L',L}(R)  <---->  N_{L',L}(R)
  !   i_l(r)       <---->  j_l(r)
  !   k_l(r)       <---->  n_l(r)
  !
  !write(*,*) "energy= ", energy, " V_II= ", vconst
  if (energy >= vconst) then
     !write(*,*) " E >= V_II "
     kappa = sqrt(2.d0*(energy - vconst))
  else
     !write(*,*) " E < V_II "
     kappa = sqrt(2.d0*(vconst - energy))
  endif
  if (debug > 0) write(*,*) "length of wavevector in region II   kappa = ", kappa

  !------------------------------------------------------------------------------!
  if (debug > 0) write(*,*) "computing radial integrals..."
  !
  ! Radial integrals
  ! ----------------
  !
  ! Radial integrals are solved by Gauss-Chebyshev quadrature of the second
  ! kind:
  !
  !   /rmax     2   2          n=nr                    2       2
  !   |     dr r   f (r)  = sum     1/2 (rmax - rmin) r   w   f (r )
  !   /rmin                    n=1                     n   n      n 
  !
  ! The nodes r_n are given by
  !
  !   z_n = n/(N+1)  ,    x_n = cos(pi z_n)   ,   r_n = (1+x_n)/2 * rmax + (1-x_n)/2 * rmin
  !
  ! and the weights are
  !
  !   w_n = pi/(N+1) sqrt((1-x_n)*(1+x_n))
  !

  ! The nodes x_n and weights w_n on the interval [-1,1] are the same independently
  ! of the integration limits and have to be calculated only once
  do n=1,nr
     zr(n) = n/(nr+1.d0)
  enddo
  xr = cos(pi * zr)
  wr = pi/(nr+1.d0) * sqrt((1.d0-xr)*(1.d0+xr))

  ! radial integrals of atomic wavefunctions of region I
  radint_If(0:lmax,1:nat) = 0.d0
  do i=1,nat
     rho_i = radii(i)
     rmin = 0.d0
     rmax = rho_i
     drdx = 0.5d0*(rmax-rmin)
     ! enumerate radial grid points covering the interval [0,rho_i]
     do n=1,nr
        x = xr(n)
        dx = wr(n)
        ! convert node positions from interval [-1,1] to [rmin,rmax]
        r = 0.5d0*((1+x)*rmax + (1-x)*rmin)
        r2dr = r**2 * drdx * dx
        do l=0,lmax
           ! evaluate radial wavefunctions
           call splev_uniform(fIknots(1:nsplI,l,i),nsplI,fIcoeffs(1:nsplI,l,i), &
                kspl, &
                r, fI(l,i), ifail)
           radint_If(l,i) = radint_If(l,i) + fI(l,i)**2 * r2dr
        enddo
     enddo
  enddo

  ! radial integrals of wavefunction outside molecular sphere (rhoIII < r)
  radint_IIIf(0:lmax) = 0.d0
  rmin = rhoIII
  ! Although the integrals extends to r=+oo, the B-splines for f^(III)_l(r) cover
  ! only a finite range. Therefore we choose the outer most knots as rmax. The
  ! interpolation ranges are the same for all l
  rmax = fIIIknots(nsplIII,0)
  drdx = 0.5d0*(rmax-rmin)
  ! enumerate radial grid points covering the interval [rhoIII,rmax]
  do n=1,nr
     x = xr(n)
     dx = wr(n)
     ! convert node positions from interval [-1,1] to [rmin,rmax]
     r = 0.5d0*((1+x)*rmax + (1-x)*rmin)
     r2dr = r**2 * drdx * dx
     do l=0,lmax
        ! evaluate radial wavefunctions
        call splev_uniform(fIIIknots(1:nsplIII,l),nsplIII,fIIIcoeffs(1:nsplIII,l), &
             kspl, &
             r, fIII(l), ifail)
     enddo
     ! sum over quadrature points (at once for l=0,1,...,lmax)
     radint_IIIf(0:lmax) = radint_IIIf(0:lmax) + fIII(0:lmax)**2 * r2dr
  enddo
  !----------------------------------------------------------------------------------!

  ! ---------------------------------------------------------------------------------!
  !
  ! Wronskians
  !

  ! evaluate atomic radial wavefunctions at the radii of the atomic
  ! spheres. The radial wavefunctions are given as B-splines. 
  do i=1,nat
     rho_i = radii(i)
     do l=0,lmax
        ! evaluate radial wavefunctions
        call splev_uniform(fIknots(1:nsplI,l,i),nsplI,fIcoeffs(1:nsplI,l,i), &
             kspl, &
             rho_i, fI(l,i), ifail)
        ! and the r-derivative
        call splev_deriv_uniform(fIknots(1:nsplI,l,i),nsplI,fIcoeffs(1:nsplI,l,i), &
             kspl, &
             rho_i, fIp(l,i), ifail)
     enddo
  enddo

  ! evaluate radial wavefunctions of region III at the radius of the molecular sphere (r=rhoIII)
  ! The radial wavefunctions are given as B-splines. 
  do l=0,lmax
     ! evaluate radial wavefunctions
     call splev_uniform(fIIIknots(1:nsplIII,l),nsplIII,fIIIcoeffs(1:nsplIII,l), &
          kspl, &
          rhoIII, fIII(l), ifail)
     ! and the r-derivative
     call splev_deriv_uniform(fIIIknots(1:nsplIII,l),nsplIII,fIIIcoeffs(1:nsplIII,l), &
          kspl, &
          rhoIII, fIIIp(l), ifail)
  enddo

  ! For each l-value there are 2*l+1 quantum numbers m=-l,...,0,...,l.
  ! For l=0,1,...,lmax there are nlm = sum_l=0^lmax 2*l+1 = (lmax+1)^2
  ! angular momentum channels.
  nlm = (lmax+1)**2
  ! dimension of the matching matrix
  ndim = (nat+1)*nlm
  
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

  ! compute w^{II}
  do i=1,nat
     ! radius of atomic sphere i
     rho_i = radii(i)
     if (energy < vconst) then
        ! E < V_II
        ! evaluate modified spherical Bessel functions i_l and k_l for l=0,...,lmax
        call mod_bessel(kappa*rho_i, lmax, &
             iII, kII, iIIp, kIIp, ifail)
     else
        ! E >= V_II
        ! compute spherical Bessel (kfn=1) functions j_l and n_l for l=0,...,lmax.
        ! In an abuse of notation, j_l and n_l are called i_l and k_l to avoid having
        ! to duplicate code for E < V_II and E >= V_II.
        call coul90_scalar(kappa*rho_i, 0.d0, lmax, 1, &
             iII, kII, iIIp, kIIp, ifail)
     endif
     do l=0,lmax
        ! The factor kappa comes from the chain rule
        ! since the derivatives are with respect to kappa*r
        wronsk_numer = kII(l) * kappa*iIIp(l) - kappa*kIIp(l) * iII(l)
        wronsk_denom = fI(l,i) * kappa*iIIp(l) - fIp(l,i) * iII(l)
        ! w^{II}_{l,i}
        wrII(l,i) = wronsk_numer / wronsk_denom
     enddo
  enddo

  ! and w^{III}
  if (energy < vconst) then
     call mod_bessel(kappa*rhoIII, lmax, &
             iII, kII, iIIp, kIIp, ifail)
  else
     ! j_l and n_l are called i_l and k_l to avoid duplicating code
     call coul90_scalar(kappa*rhoIII, 0.d0, lmax, 1, &
          iII, kII, iIIp, kIIp, ifail)
  endif
  do l=0,lmax
     wronsk_numer = iII(l) * kappa*kIIp(l) - kappa*iIIp(l) * kII(l)
     wronsk_denom = fIII(l) * kappa*kIIp(l) - fIIIp(l) * kII(l)
     wrIII(l)  = wronsk_numer / wronsk_denom
  enddo
  !-----------------------------------------------------------------------!
  
  ! ordering of angular momentum channels 
  call enumerate_lm(lmax, l_arr, m_arr)
  
  if (debug > 0) write(*,*) "computing normalization coefficients ..."
  ! multiindex (i,li,mi)
  ii = 1
  ! enumerate atomic centers (i=0 means the origin)
  do i=0,nat
     ! enumerate angular momentum channels (li,mi) on center i
     do nci=1,nlm
        li = l_arr(nci)
        mi = m_arr(nci)
        
        !---------------------------------------------------------------------!
        if (i == 0) then
           ! norm^2 of wavefunction in region III
           nrm2 = wrIII(li)**2 * radint_IIIf(li)
        else
           ! i > 0
           ! norm of wavefunction in region I
           nrm2 = wrII(li,i)**2 * radint_If(li,i)
        endif
        nrm(ii) = sqrt(nrm2)
        !---------------------------------------------------------------------!
              
        ! increase counter
        ii = ii+1
     enddo
  enddo

end subroutine normalization


subroutine wavefunctions( &
     energy, vconst, chargeIII, nat, lmax, &
     rhoIII, radii, atpos, &
     nsplI, nsplIII, kspl, fIknots, fIcoeffs, fIIIknots, fIIIcoeffs, &
     nc, sol, &
     npts, x_arr, y_arr, z_arr, &
     debug, &
     wfn_arr )
  !
  ! evaluate the MS wavefunctions on a grid using the coefficients
  ! obtained by solving the matching eigenvalue equation.
  !
  ! The atomic radial wavefunction for region I and III have to be
  ! provided as B-splines. 
  !
  ! Parameters
  ! ----------
  ! energy       : float < 0, energy of bound orbital (in Hartree)
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
  ! nsplI        : integer, giving the total number of knots for
  !                the splines of the atomic radial wavefunctions in region I
  ! nsplIII      : integer, giving the total number of knots for
  !                the splines of the radial wavefunctions in region III
  ! kspl         : integer, giving the degree of the splines
  ! fIknots      : array of shape (nsplI,lmax+1,nat), which contains the B-spline knots,
  !                the knots have to be uniformly spaced. fIknots(:,l,i) are the knots belonging
  !                to the radial wavefunction with angular momentum `l` on atom `i`
  ! fIcoeffs     : array of shape (nsplI,lmax+1,nat), which contains the B-spline coefficients
  ! fIIIknots    : array of shape (nsplIII,lmax+1), which contains the B-spline knots,
  !                the knots have to be uniformly spaced. fIIIknots(:,l) are the knots belonging
  !                to the radial wavefunction with angular momentum `l` in region III
  ! fIIIcoeffs   : array of shape (nsplIII,lmax+1), which contains the B-spline coefficients
  ! nc           : int, number of eigenvectors belonging to `energy`, degeneracy of `energy`
  ! sol          : array of shape (ndim,nc) with `nc` eigenvectors belonging to a 0 eigenvalue, M.sol = 0*sol
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
  ! wfn_arr      : array of shape (nc,:) with values of wavefunctions on the grid,
  !                wfn_arr(ic,:) is the ic-th wavefunction
  !
  implicit none
  ! ... input variables ...
  double precision, intent(in) :: energy, vconst, chargeIII
  integer, intent(in) :: nat, lmax
  double precision, intent(in) :: rhoIII, radii(nat), atpos(3,nat)
  integer, intent(in) :: nsplI, nsplIII, kspl
  double precision, intent(in) :: fIknots(nsplI,0:lmax,nat), fIcoeffs(nsplI,0:lmax,nat)
  double precision, intent(in) :: fIIIknots(nsplIII,0:lmax), fIIIcoeffs(nsplIII,0:lmax)
  integer, intent(in) :: nc
  complex*16, intent(in) :: sol((lmax+1)*(lmax+1)*(nat+1),nc)
  integer, intent(in) :: npts
  double precision, intent(in) :: x_arr(npts), y_arr(npts), z_arr(npts)
  ! ... optional ...
  !f2py integer optional, intent(in) :: debug = 0
  integer, intent(in) :: debug
  ! ... output variables ...
  complex*16, intent(out) :: wfn_arr(nc,npts)
  ! ... local variables ...
  integer :: nlm, ndim
  ! wave vectors in region II
  double precision :: kappa
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
  ! for printing coefficients
  integer :: m
  !
  ! coefficients of wavefunctions
  ! -----------------------------
  !   in region III 
  !
  !                       III  III       
  ! Psi   (r) = sum    [ C    f  (k r )  ] Y  (th ,ph )
  !    III         l,m    l,m  l     0     l,m  0   0
  !
  ! where (r0,th0,ph0) are the spherical coordinates relative to the origin
  complex*16 :: CIII((lmax+1)*(lmax+1),nc)
  !  in region I  
  !                   Ii   Ii
  ! Psi (r) = sum    C    f (r ) Y  (th ,ph )
  !    i         l,m  l,m  l  i   l,m  i   i
  !
  ! where (r_i,th_i,ph_i) are the spherical coordinates relative to the atom i
  complex*16 :: CI((lmax+1)*(lmax+1),nat,nc)
  !  in region II  if E >= V_II
  !
  !                    II0                                  nat         IIi
  ! Psi  (r) = sum    A    j (kappa r ) Y  (th ,ph )  +  sum    sum    A    n (kappa r ) Y   (th ,ph )    
  !    II         l,m  l,m  l        0   l,m  0   0         i=1    l,m  l,m  l        i   l,m   i   i
  !
  ! or             if E < V_II
  !
  !                    II0                                  nat         IIi
  ! Psi  (r) = sum    A    i (kappa r ) Y  (th ,ph )  +  sum    sum    A    k (kappa r ) Y   (th ,ph )
  !    II         l,m  l,m  l        0   l,m  0   0         i=1    l,m  l,m  l        i   l,m   i   i
  complex*16 :: AII0((lmax+1)*(lmax+1),nc)
  complex*16 :: AII((lmax+1)*(lmax+1),nat,nc)
  ! shape of AII array (nlm,nat)
  integer :: sh(2)
  !
  ! Wronskians
  ! ----------
  ! To find the coefficients C^{I_i}_{l,m} and C^{III}_{l,m} from
  ! the solution of the matching conditions we need to store some
  ! additional ratios of Wronskians
  !
  ! for E >= V_II:
  ! 
  !  II,+   [n_l(kappa*rho_i), j_l(kappa*rho_i)]
  ! w     = ------------------------------------
  !  l,i      [f_l^i(rho_i), j_l(kappa*rho_i)]
  !
  !  III,+   [j_l(kappa*rhoIII), n_l(kappa*rhoIII)]
  ! w      = --------------------------------------
  !   l       [f_l^III(rhoIII), n_l(kappa*rhoIII)]
  !
  ! and for E < V_II
  ! 
  !  II,-   [k_l(kappa*rho_i), i_l(kappa*rho_i)]
  ! w     = ------------------------------------
  !  l,i      [f_l^i(rho_i), i_l(kappa*rho_i)]
  !
  !  III,-   [i_l(kappa*rhoIII), k_l(kappa*rhoIII)]
  ! w      = --------------------------------------
  !   l       [f_l^III(rhoIII), k_l(kappa*rhoIII)]
  !
  double precision :: wrII(0:lmax,nat)
  double precision :: wrIII(0:lmax)
  double precision :: wronsk_numer, wronsk_denom
  ! radius of atomic sphere i
  double precision :: rho_i
  ! spherical Bessel functions of first kind and second kind
  ! and their r-derivatives
  double precision :: jII(0:lmax), jIIp(0:lmax)
  double precision :: nII(0:lmax), nIIp(0:lmax)
  ! modified spherical Bessel functions of first and second kind
  ! and their r-derivatives
  double precision :: iII(0:lmax), iIIp(0:lmax)
  double precision :: kII(0:lmax), kIIp(0:lmax)
  ! atomic radial functions and their r-derivatives at r=rho_i
  double precision :: fI(0:lmax), fIp(0:lmax)
  ! region III radial wavefunctions and their r-derivatives
  double precision :: fIII(0:lmax), fIIIp(0:lmax)
  integer :: ifail
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
  double precision, parameter :: print_thresh = 1.0e-4
  
  ! For each l-value there are 2*l+1 quantum numbers m=-l,...,0,...,l.
  ! For l=0,1,...,lmax there are nlm = sum_l=0^lmax 2*l+1 = (lmax+1)^2
  ! angular momentum channels.
  nlm = (lmax+1)**2
  ! dimension of the matching matrix
  ndim = (nat+1)*nlm

  ! The wavefunctions in region II depend on whether the energy E is
  ! larger or smaller than the constant potential V_II.
  ! For E >= V_II the wavefunctions in region II are spherical Bessel functions
  !    j_l(r) and n_l(r),
  ! while for E < V_II they are modified spherical Bessel functions
  !    i_l(r) and k_l(r)
  !
  !write(*,*) "energy= ", energy, " V_II= ", vconst
  if (energy >= vconst) then
     !write(*,*) " E >= V_II "
     kappa = sqrt(2.d0*(energy - vconst))
  else
     !write(*,*) " E < V_II "
     kappa = sqrt(2.d0*(vconst - energy))
  endif
  
  ! ---------- Wronskians --------------------------------------
  if (energy >= vconst) then
     ! E >= V_II,
     ! we have to compute w^{II,+}_{l,i}
     do i=1,nat
        ! radius of atomic sphere i
        rho_i = radii(i)
        ! compute spherical Bessel (kfn=1) functions j_l and n_l for l=0,...,lmax
        call coul90_scalar(kappa*rho_i, 0.d0, lmax, 1, &
             jII, nII, jIIp, nIIp, ifail)
        do l=0,lmax
           ! evaluate radial wavefunctions of region I
           call splev_uniform(fIknots(1:nsplI,l,i),nsplI,fIcoeffs(1:nsplI,l,i), &
                kspl, &
                rho_i, fI(l), ifail)
           ! and the r-derivative
           call splev_deriv_uniform(fIknots(1:nsplI,l,i),nsplI,fIcoeffs(1:nsplI,l,i), &
                kspl, &
                rho_i, fIp(l), ifail)
           ! the factor kappa comes from the chain rule
           ! since the derivatives are with respect to kappa*r
           wronsk_numer = nII(l) * kappa*jIIp(l) - kappa*nIIp(l) * jII(l)
           wronsk_denom = fI(l) * kappa*jIIp(l) - fIp(l) * jII(l)
           ! w^{II,+}_{l,i}
           wrII(l,i) = wronsk_numer / wronsk_denom           
        enddo
     enddo
     ! and w^{III,+}_l
     call coul90_scalar(kappa*rhoIII, 0.d0, lmax, 1, &
             jII, nII, jIIp, nIIp, ifail)
     do l=0,lmax
        ! evaluate radial wavefunctions of region III
        call splev_uniform(fIIIknots(1:nsplIII,l),nsplIII,fIIIcoeffs(1:nsplIII,l), &
             kspl, &
             rhoIII, fIII(l), ifail)
        ! and the r-derivative
        call splev_deriv_uniform(fIIIknots(1:nsplIII,l),nsplIII,fIIIcoeffs(1:nsplIII,l), &
             kspl, &
             rhoIII, fIIIp(l), ifail)

        wronsk_numer = jII(l) * kappa*nIIp(l) - kappa*jIIp(l) * nII(l)
        wronsk_denom = fIII(l) * kappa*nIIp(l) - fIIIp(l) * nII(l)
        wrIII(l)  = wronsk_numer / wronsk_denom
     enddo
  else
     ! E < V_II,
     ! we have to compute w^{II,-}_{l,i}
     do i=1,nat
        ! radius of atomic sphere i
        rho_i = radii(i)
        ! compute modified spherical Bessel functions i_l and k_l for l=0,...,lmax
        call mod_bessel(kappa*rho_i, lmax, &
             iII, kII, iIIp, kIIp, ifail)
        do l=0,lmax
           ! evaluate radial wavefunctions of region I
           call splev_uniform(fIknots(1:nsplI,l,i),nsplI,fIcoeffs(1:nsplI,l,i), &
                kspl, &
                rho_i, fI(l), ifail)
           ! and the r-derivative
           call splev_deriv_uniform(fIknots(1:nsplI,l,i),nsplI,fIcoeffs(1:nsplI,l,i), &
                kspl, &
                rho_i, fIp(l), ifail)
           ! the factor kappa comes from the chain rule
           ! since the derivatives are with respect to kappa*r
           wronsk_numer = kII(l) * kappa*iIIp(l) - kappa*kIIp(l) * iII(l)
           wronsk_denom = fI(l) * kappa*iIIp(l) - fIp(l) * iII(l)
           ! w^{II,-}_{l,i}
           wrII(l,i) = wronsk_numer / wronsk_denom           
        enddo
     enddo
     ! and w^{III,-}_l
     call mod_bessel(kappa*rhoIII, lmax, &
             iII, kII, iIIp, kIIp, ifail)
     do l=0,lmax
        ! evaluate radial wavefunctions of region III
        call splev_uniform(fIIIknots(1:nsplIII,l),nsplIII,fIIIcoeffs(1:nsplIII,l), &
             kspl, &
             rhoIII, fIII(l), ifail)
        ! and the r-derivative
        call splev_deriv_uniform(fIIIknots(1:nsplIII,l),nsplIII,fIIIcoeffs(1:nsplIII,l), &
             kspl, &
             rhoIII, fIIIp(l), ifail)

        wronsk_numer = iII(l) * kappa*kIIp(l) - kappa*iIIp(l) * kII(l)
        wronsk_denom = fIII(l) * kappa*kIIp(l) - fIIIp(l) * kII(l)
        wrIII(l)  = wronsk_numer / wronsk_denom
     enddo
  endif
  
  ! ----- Wavefunction coefficients ------------------------------
  ! Given the matrix of solution vectors `sol` we have to extract
  ! the coefficients for each region
  !  sol  ->  (CI,AII0,AII_i,CIII)
  ! ordering of angular momentum channels 
  call enumerate_lm(lmax, l_arr, m_arr)

  ! shape of AII_i array
  sh(1) = nlm
  sh(2) = nat

  do ic=1,nc
     ! extract coefficient in region II from solution
     AII0(1:nlm,ic) = sol(1:nlm,ic)
     AII(1:nlm,1:nat,ic) = reshape(sol(nlm+1:ndim,ic), sh)
     ! find coefficients in region I
     do i=1,nat
        do lm=1,nlm
           l = l_arr(lm)
           CI(lm,i,ic) = wrII(l,i) * AII(lm,i,ic)
        enddo
     enddo
     ! and in region III
     do lm=1,nlm
        l = l_arr(lm)
        CIII(lm,ic) = wrIII(l) * AII0(lm,ic)
     enddo
  enddo

  ! ------- print coefficients CI, CIII, AII0 and AII_i ----
  if (debug > 0) then
     do ic=1,nc
        write(*,*) "solution vector ", ic
        write(*,*) "region I coefficients C^I_{i,l,m}"
        write(*,*) "Atom i             l            m              Re             Im"
        do i=1,nat
           do lm=1,nlm
              l = l_arr(lm)
              m = m_arr(lm)
              if (abs(CI(lm,i,ic)) > print_thresh) then
                 write(*,*) i,l,m,real(CI(lm,i,ic)), aimag(CI(lm,i,ic))
              endif
           enddo
        enddo
        write(*,*) "region II coefficients A^II_{i,l,m}"
        write(*,*) "Atom i      (l,m)      |A^II_{i,l,m}|"
        do i=1,nat
           do lm=1,nlm
              l = l_arr(lm)
              m = m_arr(lm)
              if (abs(AII(lm,i,ic)) > print_thresh) then
                 write(*,*) i,l,m,real(AII(lm,i,ic)),aimag(AII(lm,i,ic))
              endif
           enddo
        enddo
        write(*,*) "region II coefficients A^II0"
        write(*,*) "   (l,m)      |A^II0_{l,m}|"
        do lm=1,nlm
           l = l_arr(lm)
           m = m_arr(lm)
           if (abs(AII0(lm,ic)) > print_thresh) then
              write(*,*) l,m,real(AII0(lm,ic)),aimag(AII0(lm,ic))
           endif
        enddo
        write(*,*) "region III coefficients C^III"
        write(*,*) "   (l,m)      |C^III_{l,m}|"
        do lm=1,nlm
           l = l_arr(lm)
           m = m_arr(lm)
           if (abs(CIII(lm,ic)) > print_thresh) then
              write(*,*) l,m,real(CIII(lm,ic)),aimag(CIII(lm,ic))
           endif
        enddo
        write(*,*) "norm^2 of all A coefficients : ", (sum(abs(AII0(:,ic))**2)+sum(abs(AII(:,:,ic))**2))
        write(*,*) "norm^2 of all C coefficients : ", (sum(abs(CIII(:,ic))**2)+sum(abs(CI(:,:,ic))**2))
        write(*,*) "norm^2 of C^I coefficients   : ", sum(abs(CI(:,:,ic))**2)
        write(*,*) "norm^2 of C^III coefficients : ", sum(abs(CIII(:,ic))**2)
     enddo
  endif
  
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
     ! evaluate wavefunction using the expressions
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

        ! evaluate radial wavefunctions of region I
        do l=0,lmax
           call splev_uniform(fIknots(1:nsplI,l,i),nsplI,fIcoeffs(1:nsplI,l,i), &
                kspl, &
                ri, fI(l), ifail)
        enddo

        ! compose wavefunction 
        do lm=1,nlm
           l = l_arr(lm)
           wfn(:) = wfn(:) + CI(lm,i,:) * fI(l) * ysph(lm)
        enddo
           
     else if (region == 2) then
        ! spherical coordinates relative to origin
        call cart2sph(x,y,z, r0,th0,ph0)
        ! spherical harmonics Y_lm
        call sphharm(th0,ph0, lmax, ysph)
        if (energy >= vconst) then
           ! E > V_II
           ! spherical Bessel functions
           call coul90_scalar(kappa*r0, 0.d0, lmax, 1, &
                jII, nII, jIIp, nIIp, ifail)
           ! --- second half of eqn. (II.16) ---
           do lm=1,nlm
              l = l_arr(lm)
              wfn(:) = wfn(:) + AII0(lm,:) * jII(l) * ysph(lm)
           enddo
        else
           ! E < V_II
           ! modified spherical Bessel functions
           call mod_bessel(kappa*r0, lmax, &
                iII, kII, iIIp, kIIp, ifail)
           ! --- second half of eqn. (II.11) ---
           do lm=1,nlm
              l = l_arr(lm)
              wfn(:) = wfn(:) + AII0(lm,:) * iII(l) * ysph(lm)
           enddo           
        endif
        do i=1,nat
           ! coordinates relative to atom i
           xi = x - atpos(1,i)
           yi = y - atpos(2,i)
           zi = z - atpos(3,i)
           ! convert them to spherical coordinates
           call cart2sph(xi,yi,zi, ri,thi,phi)
           ! spherical harmonics Y_lm
           call sphharm(thi,phi, lmax, ysph)
           if (energy >= vconst) then
              ! E > V_II
              ! spherical Bessel functions
              call coul90_scalar(kappa*ri, 0.d0, lmax, 1, &
                   jII, nII, jIIp, nIIp, ifail)
              ! --- first half of eqn. (II.16) ----
              do lm=1,nlm
                 l = l_arr(lm)
                 wfn(:) = wfn(:) + AII(lm,i,:) * nII(l) * ysph(lm)
              enddo
           else
              ! E < V_II
              ! modified spherical Bessel functions
              call mod_bessel(kappa*ri, lmax, &
                   iII, kII, iIIp, kIIp, ifail)
              ! --- first half of eqn. (II.11) ----
              do lm=1,nlm
                 l = l_arr(lm)
                 wfn(:) = wfn(:) + AII(lm,i,:) * kII(l) * ysph(lm)
              enddo
           endif
        enddo
     else
        ! region III
        ! spherical coordinates relative to origin
        call cart2sph(x,y,z, r0,th0,ph0)
        ! spherical harmonics Y_lm
        call sphharm(th0,ph0, lmax, ysph)
        ! evaluate radial wavefunctions of region III
        do l=0,lmax
           call splev_uniform(fIIIknots(1:nsplIII,l),nsplIII,fIIIcoeffs(1:nsplIII,l), &
                kspl, &
                r0, fIII(l), ifail)
        enddo
        
        ! --- eqn. (II.9) ---
        do lm=1,nlm
           l = l_arr(lm)
           wfn(:) = wfn(:) + CIII(lm,:) * fIII(l) * ysph(lm)
        enddo
     endif

     wfn_arr(:,j) = wfn
  enddo
  
end subroutine wavefunctions

subroutine convert_coefficients_AtoC( &
     energy, vconst, chargeIII, nat, lmax, &
     rhoIII, radii, atpos, &
     nsplI, nsplIII, kspl, fIknots, fIcoeffs, fIIIknots, fIIIcoeffs, &
     nc, solA, &
     debug, &
     solC )
  !
  ! given the A-coefficients (AII0 and AIIi) of a the wavefunctions
  ! compute the C-coefficients (CIII and CIi) using part of the matching equations.
  !
  ! Parameters
  ! ----------
  ! energy       : float < 0, energy of bound orbital (in Hartree)
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
  ! nsplI        : integer, giving the total number of knots for
  !                the splines of the atomic radial wavefunctions in region I
  ! nsplIII      : integer, giving the total number of knots for
  !                the splines of the radial wavefunctions in region III
  ! kspl         : integer, giving the degree of the splines
  ! fIknots      : array of shape (nsplI,lmax+1,nat), which contains the B-spline knots,
  !                the knots have to be uniformly spaced. fIknots(:,l,i) are the knots belonging
  !                to the radial wavefunction with angular momentum `l` on atom `i`
  ! fIcoeffs     : array of shape (nsplI,lmax+1,nat), which contains the B-spline coefficients
  ! fIIIknots    : array of shape (nsplIII,lmax+1), which contains the B-spline knots,
  !                the knots have to be uniformly spaced. fIIIknots(:,l) are the knots belonging
  !                to the radial wavefunction with angular momentum `l` in region III
  ! fIIIcoeffs   : array of shape (nsplIII,lmax+1), which contains the B-spline coefficients
  ! nc           : int, number of eigenvectors belonging to `energy`, degeneracy of `energy`
  ! solA         : array of shape (ndim,nc) with `nc` eigenvectors belonging to a 0 eigenvalue, M.sol = 0*sol
  !                these contain the A-coefficients
  !
  ! Optional
  ! --------
  ! debug        : > 0, print additional information
  !
  ! Returns
  ! -------
  ! solC         : array of shape (ndim,nc) with `nc` eigenvectors,
  !                these contain the C-coefficients
  !
  implicit none
  ! ... input variables ...
  double precision, intent(in) :: energy, vconst, chargeIII
  integer, intent(in) :: nat, lmax
  double precision, intent(in) :: rhoIII, radii(nat), atpos(3,nat)
  integer, intent(in) :: nsplI, nsplIII, kspl
  double precision, intent(in) :: fIknots(nsplI,0:lmax,nat), fIcoeffs(nsplI,0:lmax,nat)
  double precision, intent(in) :: fIIIknots(nsplIII,0:lmax), fIIIcoeffs(nsplIII,0:lmax)
  integer, intent(in) :: nc
  complex*16, intent(in) :: solA((lmax+1)*(lmax+1)*(nat+1),nc)
  ! ... optional ...
  !f2py integer optional, intent(in) :: debug = 0
  integer, intent(in) :: debug
  ! ... output variables ...
  complex*16, intent(out) :: solC((lmax+1)*(lmax+1)*(nat+1),nc)
  ! ... local variables ...
  integer :: nlm, ndim
  ! wave vectors in region II
  double precision :: kappa
  ! enumerate atoms
  integer :: i
  ! enumerate angular momenta
  integer :: l
  ! multiindices (l,m)
  integer :: lm
  ! enumerate solution vectors
  integer :: ic
  ! arrays with l- and m-values
  integer :: l_arr((lmax+1)*(lmax+1)), m_arr((lmax+1)*(lmax+1))
  ! for printing coefficients
  integer :: m
  !
  ! coefficients of wavefunctions
  ! -----------------------------
  !   in region III 
  !
  !                       III  III       
  ! Psi   (r) = sum    [ C    f  (k r )  ] Y  (th ,ph )
  !    III         l,m    l,m  l     0     l,m  0   0
  !
  ! where (r0,th0,ph0) are the spherical coordinates relative to the origin
  complex*16 :: CIII((lmax+1)*(lmax+1),nc)
  !  in region I  
  !                   Ii   Ii
  ! Psi (r) = sum    C    f (r ) Y  (th ,ph )
  !    i         l,m  l,m  l  i   l,m  i   i
  !
  ! where (r_i,th_i,ph_i) are the spherical coordinates relative to the atom i
  complex*16 :: CI((lmax+1)*(lmax+1),nat,nc)
  !  in region II  if E >= V_II
  !
  !                    II0                                  nat         IIi
  ! Psi  (r) = sum    A    j (kappa r ) Y  (th ,ph )  +  sum    sum    A    n (kappa r ) Y   (th ,ph )    
  !    II         l,m  l,m  l        0   l,m  0   0         i=1    l,m  l,m  l        i   l,m   i   i
  !
  ! or             if E < V_II
  !
  !                    II0                                  nat         IIi
  ! Psi  (r) = sum    A    i (kappa r ) Y  (th ,ph )  +  sum    sum    A    k (kappa r ) Y   (th ,ph )
  !    II         l,m  l,m  l        0   l,m  0   0         i=1    l,m  l,m  l        i   l,m   i   i
  complex*16 :: AII0((lmax+1)*(lmax+1),nc)
  complex*16 :: AII((lmax+1)*(lmax+1),nat,nc)
  ! shape of AIIi array (nlm,nat)
  integer :: sh(2)
  !
  ! Wronskians
  ! ----------
  ! To find the coefficients C^{I_i}_{l,m} and C^{III}_{l,m} from
  ! the solution of the matching conditions we need to store some
  ! additional ratios of Wronskians
  !
  ! for E >= V_II:
  ! 
  !  II,+   [n_l(kappa*rho_i), j_l(kappa*rho_i)]
  ! w     = ------------------------------------
  !  l,i      [f_l^i(rho_i), j_l(kappa*rho_i)]
  !
  !  III,+   [j_l(kappa*rhoIII), n_l(kappa*rhoIII)]
  ! w      = --------------------------------------
  !   l       [f_l^III(rhoIII), n_l(kappa*rhoIII)]
  !
  ! and for E < V_II
  ! 
  !  II,-   [k_l(kappa*rho_i), i_l(kappa*rho_i)]
  ! w     = ------------------------------------
  !  l,i      [f_l^i(rho_i), i_l(kappa*rho_i)]
  !
  !  III,-   [i_l(kappa*rhoIII), k_l(kappa*rhoIII)]
  ! w      = --------------------------------------
  !   l       [f_l^III(rhoIII), k_l(kappa*rhoIII)]
  !
  double precision :: wrII(0:lmax,nat)
  double precision :: wrIII(0:lmax)
  double precision :: wronsk_numer, wronsk_denom
  ! radius of atomic sphere i
  double precision :: rho_i
  ! spherical Bessel functions of first kind and second kind
  ! and their r-derivatives
  double precision :: jII(0:lmax), jIIp(0:lmax)
  double precision :: nII(0:lmax), nIIp(0:lmax)
  ! modified spherical Bessel functions of first and second kind
  ! and their r-derivatives
  double precision :: iII(0:lmax), iIIp(0:lmax)
  double precision :: kII(0:lmax), kIIp(0:lmax)
  ! atomic radial functions and their r-derivatives at r=rho_i
  double precision :: fI(0:lmax), fIp(0:lmax)
  ! region III radial wavefunctions and their r-derivatives
  double precision :: fIII(0:lmax), fIIIp(0:lmax)
  integer :: ifail
  double precision, parameter :: print_thresh = 1.0e-4
  
  ! For each l-value there are 2*l+1 quantum numbers m=-l,...,0,...,l.
  ! For l=0,1,...,lmax there are nlm = sum_l=0^lmax 2*l+1 = (lmax+1)^2
  ! angular momentum channels.
  nlm = (lmax+1)**2
  ! dimension of the matching matrix
  ndim = (nat+1)*nlm

  ! The wavefunctions in region II depend on whether the energy E is
  ! larger or smaller than the constant potential V_II.
  ! For E >= V_II the wavefunctions in region II are spherical Bessel functions
  !    j_l(r) and n_l(r),
  ! while for E < V_II they are modified spherical Bessel functions
  !    i_l(r) and k_l(r)
  !
  !write(*,*) "energy= ", energy, " V_II= ", vconst
  if (energy >= vconst) then
     !write(*,*) " E >= V_II "
     kappa = sqrt(2.d0*(energy - vconst))
  else
     !write(*,*) " E < V_II "
     kappa = sqrt(2.d0*(vconst - energy))
  endif
  
  ! ---------- Wronskians --------------------------------------
  if (energy >= vconst) then
     ! E >= V_II,
     ! we have to compute w^{II,+}_{l,i}
     do i=1,nat
        ! radius of atomic sphere i
        rho_i = radii(i)
        ! compute spherical Bessel (kfn=1) functions j_l and n_l for l=0,...,lmax
        call coul90_scalar(kappa*rho_i, 0.d0, lmax, 1, &
             jII, nII, jIIp, nIIp, ifail)
        do l=0,lmax
           ! evaluate radial wavefunctions of region I
           call splev_uniform(fIknots(1:nsplI,l,i),nsplI,fIcoeffs(1:nsplI,l,i), &
                kspl, &
                rho_i, fI(l), ifail)
           ! and the r-derivative
           call splev_deriv_uniform(fIknots(1:nsplI,l,i),nsplI,fIcoeffs(1:nsplI,l,i), &
                kspl, &
                rho_i, fIp(l), ifail)
           ! the factor kappa comes from the chain rule
           ! since the derivatives are with respect to kappa*r
           wronsk_numer = nII(l) * kappa*jIIp(l) - kappa*nIIp(l) * jII(l)
           wronsk_denom = fI(l) * kappa*jIIp(l) - fIp(l) * jII(l)
           ! w^{II,+}_{l,i}
           wrII(l,i) = wronsk_numer / wronsk_denom           
        enddo
     enddo
     ! and w^{III,+}_l
     call coul90_scalar(kappa*rhoIII, 0.d0, lmax, 1, &
             jII, nII, jIIp, nIIp, ifail)
     do l=0,lmax
        ! evaluate radial wavefunctions of region III
        call splev_uniform(fIIIknots(1:nsplIII,l),nsplIII,fIIIcoeffs(1:nsplIII,l), &
             kspl, &
             rhoIII, fIII(l), ifail)
        ! and the r-derivative
        call splev_deriv_uniform(fIIIknots(1:nsplIII,l),nsplIII,fIIIcoeffs(1:nsplIII,l), &
             kspl, &
             rhoIII, fIIIp(l), ifail)

        wronsk_numer = jII(l) * kappa*nIIp(l) - kappa*jIIp(l) * nII(l)
        wronsk_denom = fIII(l) * kappa*nIIp(l) - fIIIp(l) * nII(l)
        wrIII(l)  = wronsk_numer / wronsk_denom
     enddo
  else
     ! E < V_II,
     ! we have to compute w^{II,-}_{l,i}
     do i=1,nat
        ! radius of atomic sphere i
        rho_i = radii(i)
        ! compute modified spherical Bessel functions i_l and k_l for l=0,...,lmax
        call mod_bessel(kappa*rho_i, lmax, &
             iII, kII, iIIp, kIIp, ifail)
        do l=0,lmax
           ! evaluate radial wavefunctions of region I
           call splev_uniform(fIknots(1:nsplI,l,i),nsplI,fIcoeffs(1:nsplI,l,i), &
                kspl, &
                rho_i, fI(l), ifail)
           ! and the r-derivative
           call splev_deriv_uniform(fIknots(1:nsplI,l,i),nsplI,fIcoeffs(1:nsplI,l,i), &
                kspl, &
                rho_i, fIp(l), ifail)
           ! the factor kappa comes from the chain rule
           ! since the derivatives are with respect to kappa*r
           wronsk_numer = kII(l) * kappa*iIIp(l) - kappa*kIIp(l) * iII(l)
           wronsk_denom = fI(l) * kappa*iIIp(l) - fIp(l) * iII(l)
           ! w^{II,-}_{l,i}
           wrII(l,i) = wronsk_numer / wronsk_denom           
        enddo
     enddo
     ! and w^{III,-}_l
     call mod_bessel(kappa*rhoIII, lmax, &
             iII, kII, iIIp, kIIp, ifail)
     do l=0,lmax
        ! evaluate radial wavefunctions of region III
        call splev_uniform(fIIIknots(1:nsplIII,l),nsplIII,fIIIcoeffs(1:nsplIII,l), &
             kspl, &
             rhoIII, fIII(l), ifail)
        ! and the r-derivative
        call splev_deriv_uniform(fIIIknots(1:nsplIII,l),nsplIII,fIIIcoeffs(1:nsplIII,l), &
             kspl, &
             rhoIII, fIIIp(l), ifail)

        wronsk_numer = iII(l) * kappa*kIIp(l) - kappa*iIIp(l) * kII(l)
        wronsk_denom = fIII(l) * kappa*kIIp(l) - fIIIp(l) * kII(l)
        wrIII(l)  = wronsk_numer / wronsk_denom
     enddo
  endif
  
  ! ----- Wavefunction coefficients ------------------------------
  ! Given the matrix of solution vectors `sol` we have to extract
  ! the coefficients for each region
  !  sol  ->  (CI,AII0,AII_i,CIII)
  ! ordering of angular momentum channels 
  call enumerate_lm(lmax, l_arr, m_arr)

  ! shape of AII_i array
  sh(1) = nlm
  sh(2) = nat

  do ic=1,nc
     ! extract coefficient in region II from solution
     AII0(1:nlm,ic) = solA(1:nlm,ic)
     AII(1:nlm,1:nat,ic) = reshape(solA(nlm+1:ndim,ic), sh)
     ! find coefficients in region I
     do i=1,nat
        do lm=1,nlm
           l = l_arr(lm)
           CI(lm,i,ic) = wrII(l,i) * AII(lm,i,ic)
        enddo
     enddo
     ! and in region III
     do lm=1,nlm
        l = l_arr(lm)
        CIII(lm,ic) = wrIII(l) * AII0(lm,ic)
     enddo
  enddo
  
  ! ------- print coefficients CI, CIII, AII0 and AII_i ----
  if (debug > 0) then
     do ic=1,nc
        write(*,*) "solution vector ", ic
        write(*,*) "region I coefficients C^I_{i,l,m}"
        write(*,*) "Atom i             l            m              Re             Im"
        do i=1,nat
           do lm=1,nlm
              l = l_arr(lm)
              m = m_arr(lm)
              if (abs(CI(lm,i,ic)) > print_thresh) then
                 write(*,*) i,l,m,real(CI(lm,i,ic)), aimag(CI(lm,i,ic))
              endif
           enddo
        enddo
        write(*,*) "region II coefficients A^II_{i,l,m}"
        write(*,*) "Atom i      (l,m)      |A^II_{i,l,m}|"
        do i=1,nat
           do lm=1,nlm
              l = l_arr(lm)
              m = m_arr(lm)
              if (abs(AII(lm,i,ic)) > print_thresh) then
                 write(*,*) i,l,m,real(AII(lm,i,ic)),aimag(AII(lm,i,ic))
              endif
           enddo
        enddo
        write(*,*) "region II coefficients A^II0"
        write(*,*) "   (l,m)      |A^II0_{l,m}|"
        do lm=1,nlm
           l = l_arr(lm)
           m = m_arr(lm)
           if (abs(AII0(lm,ic)) > print_thresh) then
              write(*,*) l,m,real(AII0(lm,ic)),aimag(AII0(lm,ic))
           endif
        enddo
        write(*,*) "region III coefficients C^III"
        write(*,*) "   (l,m)      |C^III_{l,m}|"
        do lm=1,nlm
           l = l_arr(lm)
           m = m_arr(lm)
           if (abs(CIII(lm,ic)) > print_thresh) then
              write(*,*) l,m,real(CIII(lm,ic)),aimag(CIII(lm,ic))
           endif
        enddo
        write(*,*) "norm^2 of all A coefficients : ", (sum(abs(AII0(:,ic))**2)+sum(abs(AII(:,:,ic))**2))
        write(*,*) "norm^2 of all C coefficients : ", (sum(abs(CIII(:,ic))**2)+sum(abs(CI(:,:,ic))**2))
        write(*,*) "norm^2 of C^I coefficients   : ", sum(abs(CI(:,:,ic))**2)
        write(*,*) "norm^2 of C^III coefficients : ", sum(abs(CIII(:,ic))**2)
     enddo
  endif

  ! save C-coefficients in a single vector (CIII, CI1, CI2,...,)
  do ic=1,nc
     solC(1:nlm,ic) = CIII(1:nlm,ic)
     solC(nlm+1:ndim,ic) = pack(CI(1:nlm,1:nat,ic), .true.)
  enddo

end subroutine convert_coefficients_AtoC

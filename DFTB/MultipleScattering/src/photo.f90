! 
! Photoelectron angular distribution for isotropic distribution of molecules
! in the gas phase. The total cross section sigma and anisotropy parameters beta1
! and beta2 are obtained according to Ref. [1]
!
! References
! ----------
! [1] B.Ritchie,
! "Theory of the angular distribution of photoelectrons
!  ejected from optically active molecules and molecular ions"
! Phys. Rev. A, vol. 13, num. 4 (1976), p. 1411-1415.
!

subroutine pad_iso(energy, pol, tdip, lmax, norb, &
     pad)
  !
  ! compute the orientation-averaged photoelectron angular distribution (PAD)
  ! for an isotropic ensemble (all orientations are equally probable) of molecules
  ! in the gas phase according to eqn. (11a) in Ref. [1].
  !
  ! The wavefunction of the final state if a function of the position vector r
  ! and the momentum of the ejected electron k. The partial wave resolution of
  ! Psi^(-)(r,k) with respect to k reads:
  !
  !      (-)                                (-)   *            (l) *
  !   Psi   (r,k)  = 4 pi sum sum  sum   Psi     Y (kth,kph)  D      (a,b,g)
  !                          l   m    m'    l,m   l,m'         m',m
  !
  !
  ! The photoangular distribution has the form 
  !                     L=2
  !     PAD(theta) = sum     c  P (cos(theta))       with P : Legendre polynomial
  !                     L=0   L  L                         L
  ! 
  !                  sigma          L=2
  !                = ----- [ 1 + sum     beta  P (cos(theta))
  !                   4 pi          L=1      L  L
  !
  ! theta is the angle between the electron momentum k and the polarization direction
  ! of the electric field.
  !
  ! Parameters
  ! ----------
  ! energy      : float, kinetic energy 1/2 k^2 of photoelectron (in Hartree)
  ! pol         : integer, polarization of ionization radiation
  !               -1 (left), 0 (linear), +1 (right)
  ! tdip        : complex array of shape (norb,3,nlm) with transition dipoles
  !               between bound and (-) continuum orbitals
  !               The transition dipoles between the bound orbital Psi_b(r) and the partial
  !               wave Psi^(-)_{l,m}(r) are defined as
  !                                      (-)
  !                  tdip[b,mu,lm] = <Psi   | r Y    |Psi >    
  !                                      l,m     1,mu    b
  !               where i=1,...,norb and mu=-1,0,1 and lm=1,...,nlm is a multiindex.
  ! lmax        : integer, highest angular momentum in partial wave expansion
  !               of final continuum state, there are nlm = (lmax+1)**2 angular
  !               momentum channels (l,m)
  ! norb        : integer, number of initial bound orbitals
  !
  ! Returns
  ! -------
  ! pad         : float array of shape (norb,3) with total photoionization cross sections
  !               and anisotropy parameters for each bound orbital. For the i-th orbital
  !               the three parameters describing the PAD are 
  !                 sigma = pad(i,1)
  !                 beta1 = pad(i,2)
  !                 beta2 = pad(i,3)
  !               In non-chiral molecules beta1 should be zero.
  !
  implicit none
  ! ... input variables ...
  integer, intent(in) :: norb, lmax
  complex*16, intent(in) :: tdip(norb,-1:1,(lmax+1)*(lmax+1))
  integer, intent(in) :: pol
  double precision, intent(in) :: energy
  ! ... output variables ...
  double precision, intent(out) :: pad(norb,3)
  ! ... local variables ...
  ! pi = 3.14...
  double precision, parameter :: pi = 4.d0*datan(1.d0)
  ! coefficients in the expansion PAD(th) = sum_L c_L P(cos(th))
  complex*16 :: c(norb,0:2)
  ! enumerate P_l
  integer :: l
  ! number of angular momentum channels
  integer :: nlm
  ! multiindices for enumerating angular momentum channels
  integer :: lm1,lm2
  ! angular momentum quantum numbers
  integer :: l1,m1, l2,m2
  ! components of dipole operator r*Y_{1,mu}
  integer :: mu1,mu2
  ! arrays with l- and m-values for mapping multiindex lm -> (l,m)
  integer :: l_arr((lmax+1)*(lmax+1)), m_arr((lmax+1)*(lmax+1))
  ! transition dipoles to specific angular momentum channel
  complex*16 :: tdip1c(norb,-1:1), tdip2(norb,-1:1)
  ! universal function Phi(....,l) for l=0,1,2
  double precision :: phi(0:2)
  ! units of cross section
  double precision :: units
  ! photoelectron velocity 
  double precision :: k
  ! speed of light in atomic units, inverse of fine structure
  ! constant alpha = 1/speed_of_light
  double precision, parameter :: speed_of_light = 137.035999139d0
  
  nlm = (lmax+1)**2
  ! ordering of angular momentum channels 
  call enumerate_lm(lmax, l_arr, m_arr)
  
  c(:,0:2) = (0.d0, 0.d0)
  do lm1=1,nlm
     l1 = l_arr(lm1)
     m1 = m_arr(lm1)
     tdip1c = conjg(tdip(:,:,lm1))
     do lm2=1,nlm
        l2 = l_arr(lm2)
        m2 = m_arr(lm2)
        
        ! Most combinations of (l1,m1,l2,m2) do not contribute
        ! to the cross sections unless the following restrictions
        ! are met:
        !  |l1-l| <= l2 <= l1+l        with l=0,1,2
        if (.not.((max(0,l1-2) <= l2).and.(l2 <= l1+2))) cycle
        !  |m2-m1| <= 2
        if (.not.(abs(m2-m1) <= 2)) cycle

        tdip2 = tdip(:,:,lm2)
        do mu1=-1,1
           do mu2=-1,1
              call phi_iso(pol,l1,m1,mu1, l2,m2,mu2, phi)
              do l=0,2
                 c(:,l) = c(:,l) + tdip1c(:,mu1)*tdip2(:,mu2) * phi(l)
              enddo
           enddo
        enddo        
     enddo
  enddo
  ! photoelectron velocity
  k = sqrt(2*energy)
  ! The cross section is expressed in atomic units (bohr^2)
  units = 8.d0/3.d0 * pi / speed_of_light * k
  c = c * units
  ! NOTE: eqn. (11a) also contains the photon energy E_p as a
  !       factor, which is not included here.

  !
  ! PAD(th) = sigma/(4*pi) [1 + beta1 P1(cos(th)) + beta2 P2(cos(th))]
  !
  ! sigma = c0*(4pi)
  pad(:,1) = dble( c(:,0)*(4.d0*pi) )
  ! beta1 = 4pi/sigma * c1
  pad(:,2) = dble( c(:,1)/c(:,0) )
  ! beta2 = 4pi/sigma * c2
  pad(:,3) = dble( c(:,2)/c(:,0) )

  !
  ! Perform some sanity checks:
  !
  !   The PAD should be a real number.
  if (any( abs(aimag(c)) > 1.0e-10*maxval(abs(c)) )) then
     write(*,*) "BUG: imaginary part of PAD not zero!"
     write(*,*) "  c = ", c
     stop
  endif
  !   PAD(th) should be positive, therefore beta2 has to lie in the range
  !   -1 <= beta2 <= 2
  if (any( (pad(:,3) < -1.d0).or.(2.d0 < pad(:,3)) )) then
     write(*,*) "BUG: beta2 outside valid interval [-1,2]!"
     write(*,*) "  beta2 = ", pad(:,3)
     stop
  endif
  !   Alternatively, the total orientation-averaged cross section can be
  !   computed directly from the transition dipoles as
  !
  !                                      mu=+1     l=oo     m=+l                 2
  !      sigma = (4 pi * units) * 1/3 sum       sum      sum      |tdip[:,mu,lm]| 
  !                                      mu=-1     l=0      m=-l
  !
  !   The factor 1/3 comes from the averaging over orientations.
  if (any( &
       abs( &
       pad(:,1) - 4.d0*pi*units/3.d0 * sum(sum(abs(tdip)**2, 2), 2) &
       )/pad(:,1) > 1.0e-10)) then
     write(*,*) "BUG: sigma from |tdip|^2 and orientation averaging disagree!"
     stop
  endif
end subroutine pad_iso


subroutine phi_iso(p,l1,m1,mu1, l2,m2,mu2, &
     phi)
  ! 
  ! evaluate the universal function in eqn. (11b) of Ref. [1]
  !
  !                                    p+mu2+m2                             1/2
  ! Phi(p;l1,m1,l2,m2;mu1,mu2,l) = (-1)         (2*l+1) * [(2*l1+1)(2*l2+1)]
  !
  !     (l2 l1 l) (1  1  l) (l2  l1      l   ) ( 1     1       l   )
  !  x  (       ) (       ) (                ) (                   )
  !     (0  0  0) (p -p  0) (m2 -m1  -(m2-m1)) (mu2  -mu1  -(m2-m1))
  !
  !        (1)       (2)           (3)                  (4)
  !
  ! The Wigner 3j symbols (1),(2) and (4) impose restrictions onto the input values
  ! for which the products are non-zero. For given integers p,l1,m1,mu1,mu2
  ! only certain values of l,l2 and m2 will lead to a non-zero Phi:
  !
  !   (1)   =>        0 <= l <= 2
  !   (2)   =>   |l1-l| <= l2 <= l1+l
  !   (4)   =>        m2 = mu2-mu1+m1   and  |m2-m1| <= l <= 2
  !
  ! The function Phi(...) is evaluated for l=0,1,2 at once and the result is stored
  ! in the array phi(0:2).
  !
  ! INPUT
  !    p,l1,m1,l2,m2,mu1,mu2 : integers with the ranges
  !                            -1 <= p <= 1,
  !                             0 <= l1, -l1 <= m1 <= l1,
  !                             0 <= l2, -l2 <= m2 <= l2,
  !                             -1 <= mu1 <= 1, -1 <= mu2 <= 1
  ! OUTPUT
  !    phi                   : values of Phi(...) for l=0,1,2
  !
  implicit none
  ! ... input variables ...
  integer, intent(in) :: p,l1,m1,mu1,l2,m2,mu2
  ! ... output variables ...
  double precision, intent(out) :: phi(0:2)
  ! ... local variables ...
  ! Wigner 3j symbols
  double precision :: w3j1(l1+l2+1), w3j2(2+1), w3j3(l1+l2+1), w3j4(2+1)
  ! minimum and maximum values of Wigner 3j output arrays
  integer :: jmin1,jmax1, jmin2,jmax2, jmin3,jmax3, jmin4,jmax4
  ! exit status of Wigner 3j calculation
  integer :: ifail
  double precision :: fac
  ! enumerate Legendre polynomials P_l(cos(th))
  integer :: l

  phi(0:2) = 0.d0
  ! condition (4) determines the value of m2
  if (m2 /= mu2-mu1+m1) then
     ! mu2 does not satisfy condition (4), Phi = 0
     return
  endif
  
  ! compute Wigner 3j symbols (1),(2),(3) and (4)
  !
  !  ( l  l2 l1 )
  !  ( m  m2 m1 ) is zero except for jminm <= l <= jmaxm
  ! Therefore the element w3jm(l-jminm+1) corrsponds to l.
  ! The order of the arguments in calling wigner3j is as follows:
  !   call wigner3j(w3jm, jminm, jmaxm, l2,l1, m, m2, m1, ifail)
  !
  ! (l2 l1 l)   (l l2 l1)
  ! (0  0  0) = (0 0  0 )
  call wigner3j(w3j1, jmin1, jmax1, l2,l1, 0, 0, 0, ifail)
  !
  ! (1  1  l)   (l  1  1)
  ! (p -p  0) = (0  p -p)
  call wigner3j(w3j2, jmin2, jmax2, 1, 1, 0, p, -p, ifail)
  !
  ! (l2  l1      l  )   (   l       l2   l1 )
  ! (m2 -m1  -(m2-m1) = (-(m2-m1)   m2  -m1 )
  call wigner3j(w3j3, jmin3, jmax3, l2,l1, -(m2-m1), m2, -m1, ifail)
  !
  ! (1     1       l  )   (    l      1     1  )
  ! (mu2 -mu1  -(m2-m1) = (-(m2-m1)  mu2  -mu1 )
  call wigner3j(w3j4, jmin4, jmax4, 1, 1, -(m2-m1), mu2, -mu1, ifail)
  
  ! eqn. (11b) without the factor P_l(cos(th))
  fac = (-1.d0)**(p+mu2+m2) * sqrt((2.d0*l1+1.d0)*(2.d0*l2+1.d0))
  ! Because of condition (4) |m2-m1| <= l <= 2
  do l=abs(m2-m1),2
     ! Because of condition (2) |l1-l| <= l2 <= l1+l
     if ((abs(l1-l) <= l2).and.(l2 <= l1+l)) then
        phi(l) = (2.d0*l+1.d0) * fac &
             * w3j1(l-jmin1+1) * w3j2(l-jmin2+1) &
             * w3j3(l-jmin3+1) * w3j4(l-jmin4+1)
     endif
  enddo
end subroutine phi_iso

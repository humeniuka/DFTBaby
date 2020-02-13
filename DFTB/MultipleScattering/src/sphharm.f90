!
! spherical harmonics Y_{l,m}(th,ph)
!

subroutine sphharm_arr(th, ph, lmax, ysph, n)
  !
  ! vectorized version of `sphharm()`, so that spherical harmonics
  ! may be calculated for arrays of angles
  !
  ! ... input variables ...
  ! number of angles
  integer, intent(in) :: n
  double precision, intent(in) :: th(n), ph(n)
  integer, intent(in) :: lmax
  ! ... output variables ...
  complex*16, intent(out) :: ysph((lmax+1)*(lmax+1), n)
  ! ... local variables ...
  integer :: i

  ! compute spherical harmonics Y_{0,0}, ..., Y_{lmax,lmax}
  ! for all angles
  do i=1,n
     call sphharm(th(i), ph(i), lmax, ysph(:,i))
  enddo
  
end subroutine sphharm_arr

subroutine sphharm(th,ph, lmax, ysph)
  !
  ! compute all spherical harmonics Y_{l,m}(th,ph) for
  ! l=0,1,...,lmax and -l <= m <= l recursively.
  !
  ! Parameters
  ! ----------
  ! th         : float, polar coordinate, range [0,pi]
  ! ph         : float, azimuthal coordinate, range [0,2*pi]
  ! lmax       : int, highest l-value
  !
  ! Returns
  ! -------
  ! Ylm        : array of length (lmax+1)^2 with spherical harmonics
  !              in the following order
  !               Y_(0,0),
  !               Y_(1,0), Y_(1,+1), Y_(1,-1),
  !               Y_(2,0), Y_(2,+1), Y_(2,-1), Y_(2,+2), Y_(2,-2),
  !               .....
  !               Y_(lmax,0), ...., Y_{lmax,lmax},Y_{lmax,-lmax}
  !
  !
  ! The following recursion relation allows to increase l for fixed m:
  ! 
  !                 4 l^2-1                           (l-1)^2 - m^2
  !   Y     = sqrt( ------- ) { x * Y        -  sqrt( -------------- ) Y       }          (RR)
  !    l,m          l^2-m^2          l-1,m             4 (l-1)^2 - 1    l-2,m
  !
  ! The recursion is initiated with
  !
  !              m       2 l + 1                          m/2
  !   Y    = (-1)  sqrt( --------- ) (2m-1)!! [(1-x)(1+x)]     exp(i*m*phi)
  !    m,m               4pi (2m)!
  !
  ! and
  !
  !   Y       =  x sqrt(2 m + 3) Y
  !    m+1,m                      m,m
  !
  ! Starting from  Y_{m,m} and Y_{m+1,m} we use the recursion relation (RR) to
  ! find Y_{m+2,m}, Y_{m+3,m}, ..., Y_{lmax,m}.
  !
  ! Reference
  ! ---------
  ! section 6.8 'Spherical Harmonics' in Numerical Recipes in F77
  ! The recursion relation of associated Legendre polynomials was
  ! translated into a relation for spherical harmonics.
  !
  ! input variable
  double precision, intent(in) :: th,ph
  integer, intent(in) :: lmax
  ! output variable
  ! array with spherical harmonics
  complex*16, intent(out) :: ysph((lmax+1)*(lmax+1))
  ! local variables
  complex*16, parameter :: imag = (0, 1)   ! imaginary unit
  ! pi = 3.14...
  double precision, parameter :: pi = 4.d0*datan(1.d0)
  integer :: m,l, m2,l2,lmo2
  ! zero-based array with Y_{l,m} for 0<=m<=l
  complex*16 :: y(0:lmax+1,0:lmax+1)
  double precision :: fac, fac2
  ! x = cos(th)
  double precision :: x
  ! ((1-x)*(1+x))^(1/2)
  double precision :: somx2
  ! ((1-x)*(1+x))^(m/2)
  double precision :: somx2pm
  ! exp(i*ph)
  complex*16 :: expiph
  ! exp(i*m*ph)
  complex*16 :: expimph
  ! coefficients of recursion relation, depend on l and m
  double precision :: a,b
  ! multiindex (l,m)
  integer :: lm

  ! check input
  if (lmax < 0) then
     write(*,*) "ERROR in sphharm: lmax has to be positive, but got lmax = ", lmax
     stop
  endif
  
  x = cos(th)
  somx2 = sqrt((1.d0-x)*(1.d0+x))
  expiph = exp(imag*ph)

  ! initialize factors for m=0
  fac = 1
  fac2 = 1
  somx2pm = 1.d0
  expimph = (1.d0,0.d0)
  do m=0,lmax
     ! Starting from Y_{m,m} and Y_{m+1,m} we generate all
     ! Y_{l,m} for l=m,m+1,...,lmax
     y(m,m) = (-1)**m * sqrt((2.d0*m+1.d0)/(4.d0*pi*fac)) * fac2 &
          * somx2pm * expimph
     if (m < lmax) then
        y(m+1,m) = x * sqrt(2.d0*m+3.d0) * y(m,m)

        m2 = m*m
        do l=m+2,lmax
           l2 = l*l
           lmo2 = (l-1)*(l-1)
           
           a = sqrt(dble(4*l2-1)/dble(l2-m2))
           b = sqrt(dble(lmo2-m2)/dble(4*lmo2-1))
           ! recursion relation
           y(l,m) = a * (x * y(l-1,m) - b * y(l-2,m))
        enddo
     endif

     ! update factors for m --> m+1
     expimph = expimph * expiph
     !  (2*(m+1))! = (2(m+1)) (2m+1) (2m)!
     fac = fac * (2*m+1) * (2*(m+1))
     !  (2(m+1)-1)!! = (2(m+1)-1) (2m-1)!!
     fac2 = fac2 * (2*(m+1)-1)
     somx2pm = somx2 * somx2pm

  enddo

  ! reorder spherical harmonics and compute those with negative m
  ! using Y_{l,-m} = (-1)^m Y_{l,m}^*
  lm = 1
  do l=0,lmax
     do m=0,l
        ysph(lm)   = y(l,m)
        lm = lm+1
        if (m > 0) then
           ysph(lm) = (-1)**m * conjg(y(l,m))
           lm = lm+1
        endif
     enddo
  enddo

end subroutine sphharm
  
subroutine sphharm_const_m(th,ph, mconst, lmax, y)
  !
  ! compute all spherical harmonics Y_{l,m} with constant |m|>0 and
  ! |m| <= l <= lmax.
  !
  ! Parameters
  ! ----------
  ! th         : float, polar coordinate, range [0,pi]
  ! ph         : float, azimuthal coordinate, range [0,2*pi]
  ! mconst     : int, constant m-value
  ! lmax       : int > 0, highest l-value
  !
  ! Returns
  ! -------
  ! Y          : array of length lmax+1, the first |m| entries
  !              are zero, then there are the lmax-|m|+1 spherical
  !              harmonics in the order
  !              Y_{|m|,m}, Y_{|m|+1,m}, Y_{|m|+2,m}, ..., Y_{lmax,m}
  !
  implicit none
  ! input variable
  double precision, intent(in) :: th,ph
  integer, intent(in) :: mconst, lmax
  ! output variable
  ! array with spherical harmonics, zero-based
  complex*16, intent(out) :: y(0:lmax+1)
  ! local variables
  complex*16, parameter :: imag = (0, 1)   ! imaginary unit
  ! pi = 3.14...
  double precision, parameter :: pi = 4.d0*datan(1.d0)
  integer :: m, mm, sign_m, l, m2,l2,lmo2
  double precision :: fac, fac2
  ! x = cos(th)
  double precision :: x
  ! ((1-x)*(1+x))^(m/2)
  double precision :: somx2pm
  ! exp(i*m*ph)
  complex*16 :: expimph
  ! coefficients of recursion relation, depend on l and m
  double precision :: a,b

  ! keep sign of m
  sign_m = sign(1,mconst)
  ! |m|
  m = abs(mconst)

  x = cos(th)
  somx2pm = sqrt((1.d0-x)*(1.d0+x))**m
  expimph = exp(imag*m*ph)

  ! compute (2|m|)! and (2|m|-1)!!
  fac = 1
  fac2 = 1
  do mm=0,m-1
     !  (2*(m+1))! = (2(m+1)) (2m+1) (2m)!
     fac = fac * (2*mm+1) * (2*(mm+1))
     !  (2(m+1)-1)!! = (2(m+1)-1) (2m-1)!!
     fac2 = fac2 * (2*(mm+1)-1)
  enddo

  ! First |m|-1 elements are zero
  y(:lmax) = (0.d0, 0.d0)
  ! Starting from Y_{|m|,|m|} and Y_{|m|+1,|m|} we generate all
  ! Y_{l,|m|} for l=|m|,|m|+1,...,lmax
  y(m) = (-1)**m * sqrt((2.d0*m+1.d0)/(4.d0*pi*fac)) * fac2 &
                 * somx2pm * expimph
  if (m < lmax) then
     y(m+1) = x * sqrt(2.d0*m+3.d0) * y(m)

     m2 = m*m
     do l=m+2,lmax
        l2 = l*l
        lmo2 = (l-1)*(l-1)
        
        a = sqrt(dble(4*l2-1)/dble(l2-m2))
        b = sqrt(dble(lmo2-m2)/dble(4*lmo2-1))
        ! recursion relation
        y(l) = a * (x * y(l-1) - b * y(l-2))
     enddo
  endif

  if (sign_m < 0) then
     ! Y_{l,-|m|} = (-1)^|m| Y_{l,|m|}^*
     y(:lmax) = (-1)**m * conjg(y(:lmax))
  endif
  
end subroutine sphharm_const_m

subroutine enumerate_lm(lmax, l_arr, m_arr)
  !
  ! The angular momentum quantum numbers (l,m) are enumerated in the order
  ! (0,0),
  ! (1,0), (1,+1), (1,-1),
  ! (2,0), (2,+1), (2,-1), (2,+2), (2,-2),
  ! ...., (lmax,-lmax)
  !
  ! Parameters
  ! ----------
  ! lmax      : int >= 0, largest angular momentum
  !
  ! Returns
  ! -------
  ! l_arr  : integer array of length (lmax+1)^2 with l-values
  !          in the above order
  ! m_arr  : integer array of length (lmax+1)^2 with m-values
  !          in the above order
  implicit none
  ! ... input variables ...
  integer, intent(in) :: lmax
  ! ... output variables ...
  integer, intent(out) :: l_arr((lmax+1)*(lmax+1))
  integer, intent(out) :: m_arr((lmax+1)*(lmax+1))
  ! ... local variables ...
  integer :: l,m,lm

  ! multiindex (l,m)
  lm = 1
  do l=0,lmax
     do m=0,l
        l_arr(lm) = l
        m_arr(lm) = m
        ! increase counter
        lm = lm+1
        if (m > 0) then
           l_arr(lm) = l
           m_arr(lm) = -m
           ! increase counter
           lm = lm+1
        endif
     enddo
  enddo
  
  ! consistency check
  if (lm-1 /= (lmax+1)*(lmax+1)) then
     write(*,*) "BUG  lmax= ", lmax
     stop
  endif
  
end subroutine enumerate_lm

! B-Spline interpolation

! These functions are modified versions of the files splev.f and fpbspl.f by Paul Dierckx
!
!   references :
!    de boor c  : on calculating with b-splines, j. approximation theory
!                 6 (1972) 50-62.
!    cox m.g.   : the numerical evaluation of b-splines, j. inst. maths
!                 applics 10 (1972) 134-149.
!    dierckx p. : curve and surface fitting with splines, monographs on
!                 numerical analysis, oxford university press, 1993.
!
!  author of the original subroutines:
!    p.dierckx
!    dept. computer science, k.u.leuven
!    celestijnenlaan 200a, b-3001 heverlee, belgium.
!    e-mail : Paul.Dierckx@cs.kuleuven.ac.be

!     NOTE: I modified these subroutines so that they operate on a single point x
!     instead of an array x(:) of length m, assuming that the knots are spaced
!     uniformly.
!                    A. Humeniuk      


SUBROUTINE splev_uniform(t,n,c,k,x,y,ier)
  ! This function assumes that the knots t are spaced uniformly
  ! so that the interval t_i <= x < t_(i+1) can be found without iterating over all knots
  ! subroutine splev evaluates a spline s(x) of degree k, given in its b-spline representation,
  ! at position x
  !
  !  calling sequence:
  !     call splev_uniform(t,n,c,k,x,y,ier)
  !
  !  input parameters:
  !    t    : array,length n, which contains the position of the knots.
  !    n    : integer, giving the total number of knots of s(x).
  !    c    : array,length n, which contains the b-spline coefficients.
  !    k    : integer, giving the degree of s(x).
  !    x    : contains the point where s(x) must
  !           be evaluated.
  !
  !  output parameter:
  !    y    : gives the value of s(x)
  !
  !    ier  : error flag
  !      ier = 0 : normal return

  IMPLICIT NONE
  ! INPUT
  INTEGER, INTENT(IN) :: n,k
  DOUBLE PRECISION, INTENT(IN) :: t(n),c(n),x
  ! OUTPUT
  DOUBLE PRECISION, INTENT(OUT) :: y
  INTEGER, INTENT(OUT) :: ier
  ! LOCAL VARIABLE
  INTEGER :: k1,nk1,l,ll,j
  DOUBLE PRECISION :: tb,te,arg,dt,sp, h(6)
  ier=0
!  fetch tb and te, the boundaries of the approximation interval.
  k1 = k+1
  nk1 = n-k1
  ! first and last nodes
  tb = t(k1)
  te = t(nk1+1)
  IF (x <= tb) THEN
     arg=tb
     l = k1
  ELSE IF (x >= te) THEN
     arg=te
     l = nk1
  ELSE
     arg = x
     ! find interval such that t(l) <= x < t(l+1)
     dt = t(k1+2)-t(k1+1) ! uniform distance between knots 
     l = INT((x-t(1))/dt)+1
  ENDIF
  ! If l < k, we divide by zero because the interpolating points t(1:k) = 0.0
  IF (l < k) THEN
     l = k1
  ENDIF
  !  evaluate the non-zero b-splines at x.
  call fpbspl(t,n,k,arg,l,h)
  !  find the value of s(x) at x
  sp = 0.d0
  ll = l-k1
  DO j=1,k1
     ll = ll+1
     ! linear combination of b-splines 
     sp = sp+c(ll)*h(j)
  ENDDO
  y = sp  
END SUBROUTINE splev_uniform

INCLUDE "splder.f"

SUBROUTINE splev_deriv_uniform(t,n,c,k,x,dydx,ier)
  ! evaluates the first derivative of a spline
  IMPLICIT NONE
  ! INPUT
  INTEGER, INTENT(IN) :: n,k
  DOUBLE PRECISION, INTENT(IN) :: t(n),c(n),x
  ! OUTPUT
  DOUBLE PRECISION, INTENT(OUT) :: dydx  ! derivative
  INTEGER, INTENT(OUT) :: ier ! error code, 0=OK
  ! LOCAL VARIABLES
  DOUBLE PRECISION :: wrk(n)
  ! splder wants an array
  DOUBLE PRECISION x_arr(1)
  DOUBLE PRECISION dydx_arr(1)
  ! TODO: This is very ugly, I should rewrite splder so that it accepts a scalar x
  x_arr(1)=x
  call splder(t,n,c,k,1,x_arr,dydx_arr,1,wrk,ier)
  dydx=dydx_arr(1)
  
END SUBROUTINE splev_deriv_uniform


PURE SUBROUTINE fpbspl(t,n,k,x,l,h)
  ! subroutine fpbspl evaluates the (k+1) non-zero b-splines of
  ! degree k at t(l) <= x < t(l+1) using the stable recurrence
  ! relation of de boor and cox.
  IMPLICIT NONE
  ! INPUT
  DOUBLE PRECISION, INTENT(IN) :: x
  INTEGER, INTENT(IN) :: n,k,l
  DOUBLE PRECISION, INTENT(IN) :: t(n)
  ! OUTPUT
  DOUBLE PRECISION, INTENT(OUT) :: h(6)
  ! LOCAL VARIABLES
  DOUBLE PRECISION :: f,one
  INTEGER :: i,j,li,lj
  DOUBLE PRECISION :: hh(5)
  one = 0.1e+01
  h(1) = one
  DO j=1,k
     DO i=1,j
        hh(i) = h(i)
     ENDDO
     h(1) = 0.
     DO i=1,j
        li = l+i
        lj = li-j
        f = hh(i)/(t(li)-t(lj))
        h(i) = h(i)+f*(t(li)-x)
        h(i+1) = f*(x-t(lj))
     ENDDO
  ENDDO
END SUBROUTINE fpbspl



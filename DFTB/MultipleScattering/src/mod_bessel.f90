
subroutine mod_bessel(x, nmax, &
     ib, kb, ibp, kbp)
  !
  ! compute the modified spherical Bessel functions of the first
  ! kind,
  !            -n
  !   i (x) = i   j (i*x)
  !    n           n
  ! and second kind,
  !                -n
  !   k (x) = - (i)   ( j (i*x)  + i y (i*x) )
  !    n                 n            n
  !
  ! recursively and their derivatives. j_n and y_n are spherical Bessel
  ! functions of the first and second kind. This definition of k_n(x)
  ! differs from the one in `scipy`:
  !
  !   k (x) = (-1)**n * 2/pi * scipy.special.spherical_kn(n,x)
  !    n
  !
  !
  ! References
  ! ----------
  ! [1] Abramowitz & Stegun
  ! [2] Arfken, Weber & Harris (2013),
  !     "Mathematical Methods for Physicists, Seventh Edition: A Comprehensive Guide"
  !     eqns. (14.195) - (14.198)
  !
  ! Parameters
  ! ----------
  ! x          : float > 0, argument of modified Bessel functions
  ! nmax       : int >= 0, functions are computed for n=0,1,...,nmax
  !
  ! Returns
  ! -------
  ! ib, kb     : vectors of length lmax+1 with modified spherical Bessel
  !              functions of 1st (i_n(x)) and 2nd (k_n(x)) kind.
  ! ibp, kpb   : vectors of length nmax+1 with x-derivative of i_n(x) and k_n(x)
  !
  implicit none
  ! ... input variables ...
  double precision, intent(in) :: x
  integer, intent(in) :: nmax
  ! ... output variables ...
  double precision, intent(out) :: ib(0:nmax), kb(0:nmax)
  double precision, intent(out) :: ibp(0:nmax), kbp(0:nmax)
  ! ... local variables ...
  ! for i_l(x)
  double precision :: sinhx, coshx
  ! scale factor and its derivative
  double precision :: s, sp
  ! for k_l(x)
  double precision :: expmx
  ! xi = 1/x, xi2 = 1/x**2
  double precision :: xi, xi2
  integer :: n
  ! N >> n
  integer :: nn
  integer, parameter :: nn_max = 120
  ! logarithm of machine precision
  double precision, parameter :: logeps = log(1.d-15)
  ! i_n, i_{n+1}, i_{n+2}
  double precision :: i0, i1, i2
  double precision :: tmp

  if (nmax > nn_max/2) then
     write(*,*) "Degree of modified spherical Bessel functions too large, ", nmax, " > ", nn_max/2, "!"
     write(*,*) "Increase parameter nn_max!"
     stop
  endif

  ! Find the initial N for downward iteration. For x < 1,
  !   i (x) ---> 0  for n --> oo
  !    n
  ! We wish to find the smallest N such that i_N = 0 within
  ! machine precision. If we choose a much large N, we will get
  ! overflow errors in the downward iteration. 
  if (x < 1.d0) then
     ! find N such that x^N = machine precision ~ 0
     ! We still need to ensure N >= nmax+2, otherwise the iteration
     ! does not generate all requested orders n=0,1,...,nmax.
     nn = max(nmax+2, int(logeps/log(x)))
     ! N should never exceed the hard coded maximum.
     nn = min(nn, nn_max)
  else
     ! For x > 1, we start at the hard coded maximum.
     nn = nn_max
  endif
  
  xi = 1.d0/x
  xi2 = xi**2

  !
  ! i_n(x) is computed using the three-term-recurrence relation
  !
  !  i (x)  =  i  (x)  +  (2n+3)/x  i (x)
  !   n         n+2                  n+1
  !
  ! in downward direction since the upward direction is unstable.
  ! For small x the modified spherical Bessel function of the first
  ! kind behave like
  !            x^n
  !  i (x) ~ --------
  !   n      (2n+1)!!
  !
  ! We start at some large N >> nmax with the initial values
  !
  !  i  = 1     and      i   = (2n+1)/x i (x)
  !   N                   N-1            N
  !
  ! The initial value for i_N is arbitary only the ratio
  !
  !   i_N(x)/i_{N-1}(x)
  !
  ! matters. We could also choose the initial values such that
  !
  !  i  = x/(2n+1)  and  i   = 1
  !   N                   N-1
  !
  ! At the end the series is scaled such that
  !
  !  i (x) = sinh(x)/x
  !   0
  !
  ! The derivatives i'_n are calculated from
  !
  !  (2n+3) i' (x) = (n+1) i (x)  +  (n+2) i (x)
  !          n+1            n               n+2
  !
  
  ! i_N(x)
  i2 = x/(2*nn+1)
  ! i_{N-1}(x) = (2N+1)/x i_N(x)
  i1 = 1.d0
  ! downward recursion
  do n=nn-2,0,-1
     ! i (x) = i (x)  +  (2*n+3)/x i (x)
     !  n       n+2                 n+1
     i0 = i2 + (2*n+3)*xi*i1
     ! If the downward recursion reaches the desired range, the
     ! values of i_n and i_n' are recorded in the output arrays.
     if (n <= nmax) then
        ib(n) = i0
        if (n <= nmax-1) then
           ! derivatives
           !  i' (x) = [(n+1) i (x) + (n+2) i (x) ] / (2n+3)
           !   n+1             n             n+2
           ibp(n+1) = ((n+1)*i0 + (n+2)*i2)/(2*n+3)
        endif
     endif
     ! shift n -> n-1
     tmp = i1
     i1 = i0
     i2 = tmp
  enddo

  sinhx = sinh(x)
  coshx = cosh(x)
  ! normalize series by matching to the expected value for n=0
  !  i_0(x) = sinh(x)/x
  s = sinhx*xi / i0
  ib = s * ib
  ibp = s * ibp
  
  ! i'_0 has not been calculated in the do-loop.
  ibp(0) = (x*coshx - sinhx)*xi2

  !
  ! k_n(x) is computed using the three term recurrence relation
  !
  !  k (x) = k  (x)  + (2n-1)/x k  (x)
  !   n       n-2                n-1
  !
  ! with the initial values
  !
  !  k (x) = exp(-x)/x   and  k (x) = exp(-x) (1/x + 1/x**2)
  !   0                        1
  !

  expmx = exp(-x)
  ! initiate recursion for k_n(x) and its derivatives k_n'(x)
  kb(0) = expmx * xi
  kbp(0) = -expmx * (xi + xi2)
  if (nmax > 0) then
     kb(1) = -kbp(0)
     kbp(1) = -expmx * (1 + 2*(xi + xi2)) * xi
  endif
  
  ! compute k_n and k'_n by upward recursion for n=2,3,...,lmax
  do n=2,nmax
     kb(n)  = kb(n-2)  + (2*n-1)*xi*kb(n-1)
     kbp(n) = kbp(n-2) + (2*n-1)*(xi*kbp(n-1) - xi2*kb(n-1))
  enddo

  ! The factor (-1)^n differs from the usual definition for k_n(x)
  ! but it makes the expansion theorems (II.23) and (II.24) in
  ! Johnson (1973) "Scattered-Wave Theory of the Chemical Bond"
  ! look more symmetric. 
  do n=0,nmax
     kb(n) = kb(n) * (-1)**n
     kbp(n) = kbp(n) * (-1)**n
  enddo
end subroutine mod_bessel

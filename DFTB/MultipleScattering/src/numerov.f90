!
! The radial Schroedinger equation
!        d^2                        l*(l+1)
!  [ 1/2 ----  +  E  -  V(r)  - 1/2 ------- ] u (r)  =  0
!        dr^2                         r^2      l
!
! is solved for the radial wavefunction u_l(r) = r*R_l(r)
! on an equidistant grid using Numerov's method. The integration
! is performed either from r=0 outward or from r=oo inward.
!
! References
! ----------
! [1] Thijssen, "Computational Physics", chapter 2
!

subroutine numerov_out(r,v, n, E, lmax, &
     u)
  !
  ! integrate radial Schroedinger equation outward
  !
  ! The radial Schroedinger equation has the form
  !                       -2
  !    u''(r) = [ l(l+1) r   +  2 (V(r) - E)] u (r)
  !     l                                      l
  !           = f(r;E,l) u (r)
  !                    l
  !
  ! In terms of the new function
  !
  !    w(r) = [1 - h^2/12 f(r)] u(r)
  !
  ! Numerov's finite difference equation for w_i = w(r(i)) becomes
  !
  !            1 + 5/12 h^2 f_{i-1}
  !    w  = 2 ---------------------- w     -  w       for i=3,4,...,n
  !     i      1 - 1/12 h^2 f_{i-1}   i-1      i-2
  !
  ! For r --> 0 the differential equation is dominated by the l*(l+1) r^-2
  ! term, so that the initial conditions for outward integration are
  !         l+1                l+1
  !   u  = r       and   u  = r
  !    1    1             2    2
  !
  ! Parameters
  ! ----------
  ! r           :  numpy array of shape (n,), equidistant radial grid
  ! v           :  numpy array of shape (n,), potential energy V(r) at the radial grid points,
  !                v(i) = V(r(i))
  ! n           :  integer, number of grid points
  ! E           :  float, energy in Hartree
  ! lmax        :  integer, highest angular momentum
  !
  ! Returns
  ! -------
  ! u           :  numpy array of shape (n,lmax+1), solutions u_l(r) from outward integration,
  !                u(i,l) = u_l(r(i))
  !
  implicit none
  ! ... input variables ...
  double precision, intent(in) :: r(n)
  double precision, intent(in) :: v(n)
  integer, intent(in) :: n
  double precision, intent(in) :: E
  integer, intent(in) :: lmax
  ! ... output variables ...
  double precision, intent(out) :: u(n,0:lmax)
  ! ... local variables ...
  ! enumerate grid points
  integer :: i,ii
  ! enumerate angular momenta l=0,1,...,lmax
  integer :: l
  ! step size h=r(i+1)-r(i)
  double precision :: h
  ! intermediate variables
  double precision :: f(n), fp(n), fm(n), w(n)

  h = r(2)-r(1)
  
  do l=0,lmax
     f = l*(l+1)/r**2 + 2*(v-E)
     fp = 2.d0 + (5.d0/ 6.d0)*h**2 * f
     fm = 1.d0 - (1.d0/12.d0)*h**2 * f
     ! find smallest radius r(ii) such that r(ii)^l is larger
     ! than the machine precision, i.e. r(ii)^l + 1.d0 != 1.d0.
     do ii=1,n-2
        if ((r(ii)**(l+1) + 1.d0) /= 1.d0) then
           exit
        else
           u(ii,l) = r(ii)**(l+1)
        endif
     enddo
     ! initial conditions for r ~ 0
     w(ii)   = fm(ii)   * r(ii)**(l+1)
     w(ii+1) = fm(ii+1) * r(ii+1)**(l+1)
     ! integrate outward
     do i=ii+2,n
        w(i) = fp(i-1)/fm(i-1) * w(i-1) - w(i-2)
     enddo
     u(ii:,l) = w(ii:) / fm(ii:)
  enddo
  
end subroutine numerov_out

subroutine numerov_in(r,v, n, E, lmax, &
     u)
  !
  ! integrate radial Schroedinger equation inward starting from r=r(n).
  !
  ! The radial Schroedinger equation has the form
  !                       -2
  !    u''(r) = [ l(l+1) r   +  2 (V(r) - E)] u (r)
  !     l                                      l
  !           = f(r;E,l) u (r)
  !                    l
  !
  ! In terms of the new function
  !
  !    w(r) = [1 - h^2/12 f(r)] u(r)
  !
  ! Numerov's finite difference equation for w_i = w(r(i)) becomes
  !
  !            1 + 5/12 h^2 f_{i+1}
  !    w  = 2 ---------------------- w     -  w       for i=n-2,n-3,...,1
  !     i      1 - 1/12 h^2 f_{i+1}   i+1      i+2
  !
  ! For r --> oo, V(r) --> 0 and the differential equation is dominated
  ! by the energy term, so that the initial conditions for inward integration
  ! are 
  !
  !   u  = exp(-sqrt(2 |E|) r )   and  u    = exp(-sqrt(2 |E|) r   )
  !    n                     n          n-1                     n-1
  !
  ! assuming E < 0.
  !
  ! Parameters
  ! ----------
  ! r           :  numpy array of shape (n,), equidistant radial grid
  ! v           :  numpy array of shape (n,), potential energy V(r) at the radial grid points,
  !                v(i) = V(r(i))
  ! n           :  integer, number of grid points
  ! E           :  float < 0, energy in Hartree
  ! lmax        :  integer, highest angular momentum
  !
  ! Returns
  ! -------
  ! u           :  numpy array of shape (n,lmax+1), solutions u_l(r) from outward integration,
  !                u(i,l) = u_l(r(i))
  !
  implicit none
  ! ... input variables ...
  double precision, intent(in) :: r(n)
  double precision, intent(in) :: v(n)
  integer, intent(in) :: n
  double precision, intent(in) :: E
  integer, intent(in) :: lmax
  ! ... output variables ...
  double precision, intent(out) :: u(n,0:lmax)
  ! ... local variables ...
  ! enumerate grid points
  integer :: i,ii
  ! enumerate angular momenta l=0,1,...,lmax
  integer :: l
  ! step size h=r(i+1)-r(i)
  double precision :: h
  ! sqrt(2|E|)
  double precision :: sq2e
  ! intermediate variables
  double precision :: f(n), fp(n), fm(n), w(n)

  if (E >= 0) then
     write(*,*) "For inward integration energy has to be positive, but got E= ", E, " !"
     stop
  endif

  if (n < 100) then
     write(*,*) "Too few radial grid points, n=", n, " < 100!"
     stop
  endif
  
  h = r(2)-r(1)
  sq2e = sqrt(2.d0*abs(E))

  ! The outmost 4 points are set to zero, so that a cubic B-spline
  ! fit would extrapolate u(r) = 0 for r > rmax.
  do ii=n,n-4,-1
     u(ii,:) = 0.d0
  enddo
  
  ! Find the largest grid point where exp(-sqrt(2|E|) r) is numerically
  ! different from 0. For large energies and extended grids, the wavefunction
  ! at large radii will be 0 within machine precision.
  do ii=n-5,2,-1
     if (exp(-sq2e*r(ii)) > 0.d0) then
        exit
     else
        u(ii,:) = 0.d0
     endif
  enddo
  
  do l=0,lmax
     f = l*(l+1)/r**2 + 2*(v-E)
     fp = 2.d0 + (5.d0/ 6.d0)*h**2 * f
     fm = 1.d0 - (1.d0/12.d0)*h**2 * f
     ! initial conditions for r --> oo
     w(ii)   = f(ii)   * exp(-sq2e*r(ii))
     w(ii-1) = f(ii-1) * exp(-sq2e*r(ii-1))

     ! integrate inward
     do i=ii-2,1,-1
        w(i) = fp(i+1)/fm(i+1) * w(i+1) - w(i+2)
     enddo
     u(:ii,l) = w(:ii) / fm(:ii)
  enddo
  
end subroutine numerov_in


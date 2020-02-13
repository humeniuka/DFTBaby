!
! conversion between cartesian (x,y,z) and spherical (r,th,ph) coordinates
!

subroutine cart2sph_arr(x,y,z, r,th,ph, n)
  !
  ! vectorized function for transforming arrays of cartesian
  ! coordinates into arrays of spherical coordinates
  !
  implicit none
  ! ... input variables ...
  ! number of elements in arrays
  integer, intent(in) :: n
  ! arrays of cartesian coordinates
  double precision, intent(in) :: x(n), y(n), z(n)
  ! ... output variables ...
  ! arrays of spherical coordinates
  double precision, intent(out) :: r(n), th(n), ph(n)
  ! .. local variables ...
  integer :: i

  do i=1,n
     call cart2sph(x(i),y(i),z(i), r(i),th(i),ph(i))
  enddo
end subroutine cart2sph_arr

subroutine cart2sph(x,y,z, r,th,phi)
  !
  !  convert cartesian to spherical coordinates by inverting 
  !  the equations
  !
  !       x = r*sin(th)*cos(ph)
  !       y = r*sin(th)*sin(ph)
  !       z = r*cos(th)
  !
  !  Parameters
  !  ----------
  !  x,y,z      : 3 floats, cartesian coordinates
  !
  !  Returns
  !  -------
  !  r,th,phi   : 3 floats, spherical coordinates
  !               angles th and ph are in radians
  !
  implicit none
  ! ... input variables ...
  ! cartesian coordinates
  double precision, intent(in) :: x,y,z
  ! ... output variables ...
  ! spherical coordinates
  double precision, intent(out) :: r,th,phi
  ! .. local variables ...
  ! pi = 3.14...
  double precision, parameter :: pi = 4.d0*datan(1.d0)

  r = sqrt(x*x+y*y+z*z)
  ! polar angle theta
  if (r > 0.d0) then
     th = acos(z/r)
  else
     ! for r=0 the angle theta is not well defined, choose th = 0
     th = 0.d0
  endif

  ! azimuthal angle phi
  if (x /= 0.d0) then
     if (y /= 0.d0) then
        phi = atan2(y,x)
        if (phi < 0.d0) then
           ! translate angles from the range [-pi, pi]
           ! to the range [0, 2.0*pi]
           phi = 2.d0*pi + phi
        endif
     else
        ! if y == 0, phi=0 for x positive
        ! and phi=pi for x negative
        phi = (1.d0 - sign(1.d0,x))/2.d0 * pi
     endif
  else
     ! if x == 0, phi = pi/2 if y is positive
     ! and phi=3/2*pi if y is negative
     phi = pi/2.d0 + (1.d0 - sign(1.d0,y))/2.d0 * pi
  endif
end subroutine cart2sph
  

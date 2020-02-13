!
! Fortran 90 wrapper around Ross Barnett's COUL90 module for the
! calculation of Coulomb and Bessel functions
!
! COUL90 was obtained from
!    http://www.fresco.org.uk/programs/barnett/index.htm
!
! References
! ----------
! Comput. Phys. Commun. 8 (1974) 377-395
! Comput. Phys. Commun. 11 (1976) 141-142
! and others (see COUL90.FOR for a complete list)
!

subroutine coul90_scalar(x, eta, lmax, kfn, &
     fc,gc,fcp,gcp)
  !
  ! compute Coulomb functions or spherical Bessel functions
  ! and their derivatives for a range of l-values
  !
  ! Parameters
  ! ----------
  ! x          : float > 0, argument of Coulomb functions
  ! eta        : float, -Z/k
  ! lmax       : int >= 0, functions are computed for l=0,1,...,lmax
  ! kfn        : int, choise of functions to be computed
  !              = 0     real Coulomb functions and derivatives F & G
  !              = 1  spherical Bessel   "       "      "       j & y
  !              = 2  cylindrical Bessel "       "      "       J & Y
  !
  ! Returns
  ! -------
  ! fc, gc     : vectors of length lmax+1 with regular and irregular
  !              Coulomb functions
  ! fcp,gcp    : vectors of length lmax+1 with x-derivative of F and G
  !
  implicit none
  ! ... input variables ...
  double precision, intent(in) :: x
  double precision, intent(in) :: eta
  integer, intent(in) :: lmax
  integer, intent(in) :: kfn
  ! ... output variables ...
  double precision, intent(out) :: fc(lmax+1),gc(lmax+1)
  double precision, intent(out) :: fcp(lmax+1),gcp(lmax+1)
  ! ... local variables ...
  ! error code (see COUL90.FOR)
  !   = 0 : calculation satisfactory
  integer :: ifail

  call coul90(x,1.d0*eta,0.d0,lmax, fc,gc,fcp,gcp, kfn,ifail)
  
  if (ifail /= 0) then
     write(*,*) "ERROR In coul90: ifail = ", ifail
     stop
  endif
  
end subroutine coul90_scalar

subroutine coul90_array(x, n, eta, lmax, kfn, &
     fc,gc,fcp,gcp)
  implicit none
  ! input variables
  double precision, intent(in) :: x(n)
  integer, intent(in) :: n
  double precision, intent(in) :: eta
  integer, intent(in) :: lmax
  integer, intent(in) :: kfn
  ! output variables
  double precision, intent(out) :: fc(lmax+1,n),gc(lmax+1,n)
  double precision, intent(out) :: fcp(lmax+1,n),gcp(lmax+1,n)
  ! local variables
  integer :: i
  ! error code (see COUL90.FOR)
  !   = 0 : calculation satisfactory
  integer :: ifail

  do i=1,n
     call coul90(x(i),1.d0*eta,0.d0,lmax, fc(:,i),gc(:,i),fcp(:,i),gcp(:,i), kfn,ifail)

     if (ifail /= 0) then
        write(*,*) "ERROR In coul90: ifail = ", ifail
        stop
     endif

  enddo
end subroutine coul90_array

      

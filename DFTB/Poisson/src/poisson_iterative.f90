
SUBROUTINE poisson3d(x,y,z, f, u0, nx,ny,nz, u, &
     eps, maxiter)
  !
  ! Purpose:
  ! --------
  ! solve the 3-dimensional Poisson equation
  !  __2
  !  \/ u(x,y,z) = f(x,y,z)
  ! iteratively. The second order partial derivatives are replaced
  ! by finite differences and the resulting linear algebraic equations
  ! are solved using Jacobi iterations.
  !
  ! Input
  ! -----
  ! x,y,z: 1d arrays with grid points along each axis
  ! f: values of source term at the grid points, f[i,j,l] = f(x[i],y[j],z[l])
  ! nx,ny: number of grid points for each axis
  ! u0: initial guess, with correct boundary conditions
  !
  ! Optional
  ! --------
  ! maxiter: maximum number of Jacobi iterations
  ! eps: accuracy for convergence
  !
  ! Output
  ! ------
  ! u: solution of Poisson equation
  !
  ! release global interpreter lock
  !f2py threadsafe
  IMPLICIT NONE
  ! INPUT
  DOUBLE PRECISION, INTENT(IN) :: x(nx), y(ny), z(nz), f(nx,ny,nz), u0(nx,ny,nz)
  INTEGER, INTENT(IN) :: nx,ny,nz
  ! OPTIONAL
  !f2py double precision eps=1.0e-10
  DOUBLE PRECISION, INTENT(IN), OPTIONAL :: eps
  !f2py integer maxiter=1000000
  INTEGER, INTENT(IN), OPTIONAL :: maxiter
  ! OUTPUT
  DOUBLE PRECISION, INTENT(OUT) :: u(nx,ny,nz)

  ! LOCAL
  DOUBLE PRECISION :: dx2(nx), dy2(ny), dz2(nz) ! spacings (x[i]-x[i-1])^2, (y[j]-y[j-1])^2 and (z[l]-z[l-1])^2
  DOUBLE PRECISION :: uk(nx,ny,nz), ukp1(nx,ny,nz)  ! u^(k)_ijl and u^(k+1)_ijl
  DOUBLE PRECISION :: ax(nx,ny,nz), ay(nx,ny,nz), az(nx,ny,nz), af(nx,ny,nz)
  DOUBLE PRECISION :: dr2
  INTEGER :: k ! iteration counter
  INTEGER :: i, j, l   ! enumerate grid points
  DOUBLE PRECISION :: maxvar

  ! spacing for x-axis
  DO i=2,nx
     dx2(i) = (x(i)-x(i-1))**2
  END DO
  dx2(1)=dx2(2)
  ! spacing for y-axis
  DO j=2,ny
     dy2(j) = (y(j)-y(j-1))**2
  END DO
  dy2(1)=dy2(2)
  ! spacing for z-axis
  DO l=2,nz
     dz2(l) = (z(l)-z(l-1))**2
  END DO
  dz2(1)=dz2(2)

  ! factors
  DO l=1,nz
     DO j=1,ny
        DO i=1,nx
           dr2 = dy2(j)*dz2(l)+dx2(i)*dz2(l)+dx2(i)*dy2(j)
           ax(i,j,l) = 0.5d0*dy2(j)*dz2(l)/dr2
           ay(i,j,l) = 0.5d0*dx2(i)*dz2(l)/dr2
           az(i,j,l) = 0.5d0*dx2(i)*dy2(j)/dr2
           af(i,j,l) = 0.5d0*dx2(i)*dy2(j)*dz2(l)/dr2
        END DO
     END DO
  END DO
  ! initial guess and boundary conditions
  k=0
  uk=u0
  ukp1=u0
  DO k=1,maxiter+1
     ! Jacobi iteration
     maxvar = 0.0d0
     DO l=2,nz-1
        DO j=2,ny-1
           DO i=2,nx-1
              ukp1(i,j,l) = &
                   +ax(i,j,l)*(uk(i-1,j  ,l  )+uk(i+1,j  ,l  )) &
                   +ay(i,j,l)*(uk(i  ,j-1,l  )+uk(i  ,j+1,l  )) &
                   +az(i,j,l)*(uk(i  ,j  ,l-1)+uk(i  ,j  ,l+1)) &
                   -af(i,j,l)*f(i,j,l)
              maxvar = max(maxvar, (ukp1(i,j,l)-uk(i,j,l))**2)
           END DO
        END DO
     END DO
     uk = ukp1
     ! show convergence progress
     IF (MOD(k,50) == 0) THEN
        WRITE(*,*) "iter= ", k, "  maxvar= ", maxvar, " (threshold= ", eps**2, ")"
     ENDIF
     ! check for convergence
     IF (maxvar < eps**2) THEN
        WRITE(*,*) "solution of Poisson equation CONVERGED after ", k, " iterations!"
        EXIT
     ENDIF
     ! abort if maximum number of iterations exceeded
     IF (k == maxiter) THEN
        WRITE(*,*) "convergence failure after ", k, " iterations!"
        WRITE(*,*) "last change: max_ijl |u^(k+1)_ijl - u^(k)_ijl|^2 = ", maxvar
        STOP
     ENDIF
  END DO
  ! copy result
  u = ukp1
END SUBROUTINE poisson3d


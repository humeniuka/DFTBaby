!*****************************************************
! analytic gradients for the Thomson problem
!*****************************************************



MODULE thomson
  IMPLICIT NONE
  CONTAINS

!*****************************************************
      SUBROUTINE rij2(x, s, o, R, Z0, Nt, opt, r2, dr2d1, dr2d2)
! Purpose:
! ========
! compute the squares of the distances between Thomson
! points, rij^2, and the gradients of rij^2 w/r/t
! to the spherical or cylindrical coordinates 

      ! release global interpreter lock
      !f2py threadsafe
      IMPLICIT NONE
      ! INPUT:
      ! number of thomson points
      INTEGER, INTENT(IN) :: Nt
      DOUBLE PRECISION, INTENT(IN) :: x( Nt, 2 )! position of each Thomson point
      ! either in cylindrical or spherical coordinates
      CHARACTER, INTENT(IN) :: s( Nt )          ! type of each Thomson point
      ! 'T' => tube, 'C' => cap
      INTEGER, INTENT(IN) :: o( Nt )            ! flag indicating whether this
      ! point should be allowed to move or not
      ! radius of tube
      DOUBLE PRECISION, INTENT(IN) :: R
      ! offset along z-axis
      DOUBLE PRECISION, INTENT(IN) :: Z0
      INTEGER, INTENT(IN) :: opt
      ! OUTPUT:
      ! rij^2, distance squared between points
      DOUBLE PRECISION, INTENT(OUT) :: r2( Nt, Nt ) 
      ! gradients of rij 
      DOUBLE PRECISION, INTENT(OUT) :: dr2d1( Nt, Nt, Nt )
      DOUBLE PRECISION, INTENT(OUT) :: dr2d2( Nt,Nt,Nt ) 

      ! LOCAL VARIABLES:
      INTEGER :: i, j
      DOUBLE PRECISION :: Rd, Rq, Rdq
      DOUBLE PRECISION :: sin1( Nt ), sin2( Nt ), cos1( Nt ), cos2( Nt )
      DOUBLE PRECISION :: si1,sj1,si2,sj2, ci1,cj1,ci2,cj2, cij2, sij2

      Rd = 2*R
      Rq = R**2
      Rdq = 2*Rq

      dr2d1(:,:,:) = 0.0
      dr2d2(:,:,:) = 0.0
      ! precompute sins and cosines
      !$omp parallel do private(j)
      DO j=1, Nt
         sin1(j) = sin(x(j,1))
         sin2(j) = sin(x(j,2))
         cos1(j) = cos(x(j,1))
         cos2(j) = cos(x(j,2))
      ENDDO

      ! rij^2 and gradients
      r2(:,:) = 0.0
      !$omp parallel do private(j,sj1,cj1,sj2,cj2), &
      !$omp& private(i,si1,ci1,si2,ci2,cij2,sij2)
      DO j=1, Nt
         sj1 = sin1(j)
         cj1 = cos1(j)
         sj2 = sin2(j)
         cj2 = cos2(j)
         DO i=1, Nt
            si1 = sin1(i)
            ci1 = cos1(i)
            si2 = sin2(i)
            ci2 = cos2(i)
            cij2 = ci2*cj2 + si2*sj2
            sij2 = si2*cj2 - ci2*sj2

            IF (i == j) CYCLE

            r2(i,j) = 100000000000000.0  ! default value
            IF ((o(i) == 0) .and. (o(j) == 0) .and. (opt == 1)) CYCLE

            IF ((s(i) == 'T') .and. (s(j) == 'T')) THEN
               ! Tube-Tube
               ! distance
               r2(i,j) = Rdq*(1-cij2) + (x(i,1)-x(j,1))**2
               ! derivative with respect to z
               dr2d1(i,j,i) = Rdq*sij2
               dr2d1(i,j,j) = -dr2d1(i,j,i)
               ! derivative with respect to phi
               dr2d2(i,j,i) = 2*(x(i,1)-x(j,1))
               dr2d2(i,j,j) = -dr2d2(i,j,i)
            ELSE IF ((s(i) == 'C') .and. (s(j) == 'C')) THEN
               ! Cap-Cap
               ! distance
               r2(i,j) = Rdq*(1-(si1*sj1*cij2 + ci1*cj1))
               ! derivative with respect to theta               
               dr2d1(i,j,i) = Rdq*(si1*cj1-ci1*sj1*cij2)
               dr2d1(i,j,j) = Rdq*(ci1*sj1-si1*cj1*cij2)
               ! derivative with respect to phi
               dr2d2(i,j,i) = Rdq*si1*sj1*sij2
               dr2d2(i,j,j) = -dr2d2(i,j,i)
            ELSE IF ((s(i) == 'T') .and. (s(j) == 'C')) THEN
               ! Tube-Cap
               ! distance
               r2(i,j) = Rq*(1+sj1**2 - 2*sj1*cij2) + (x(i,1)-Z0-R*cj1)**2
               ! derivative with respect to z
               dr2d1(i,j,i) = 2*(x(i,1)-Z0-R*cj1)
               ! derivative w/r/t theta
               dr2d1(i,j,j) = Rd*(-R*cj1*cij2 + (x(i,1)-Z0)*sj1)
               ! derivative w/r/t phi
               dr2d2(i,j,i) = Rdq*sj1*sij2
               dr2d2(i,j,j) = -dr2d2(i,j,i)
            ELSE IF ((s(i) == 'C') .and. (s(j) == 'T')) THEN
               ! Cap-Tube
               ! distance
               r2(i,j) = Rq*(1+si1**2-2*si1*cij2) + (x(j,1)-Z0-R*ci1)**2
               ! derivative w/r/t z
               dr2d1(i,j,j) = 2*(x(j,1)-Z0-R*ci1)
               ! derivative w/r/t theta
               dr2d1(i,j,i) = Rd*(-R*ci1*cij2+(x(j,1)-Z0)*si1)
               ! derivative w/r/t phi
               dr2d2(i,j,i) = Rdq*si1*sij2
               dr2d2(i,j,j) = -dr2d2(i,j,i)
            ELSE
               WRITE(*,*), "Unknown type for Thomson point, s(",i,")=",s(i)," and s(",j,")=",s(j)
            END IF
         ENDDO
      ENDDO
    END SUBROUTINE rij2
 

!*****************************************************
    SUBROUTINE coulomb_energy(o, r2, dr2d1, dr2d2, Nt, coul, grad)
! Purpose:
! ========
! compute the Coulomb energy and its gradient w/r/t the position of
! the Thomson points
! The distances squared r2 and the partial gradients dr2d1 and dr2d2 should
! have been calculated by the function rij2

      ! release global interpreter lock
      !f2py threadsafe
      IMPLICIT NONE
      ! INPUT:      
      ! number of thomson points
      INTEGER, INTENT(IN) :: Nt
      ! flag indicating whether the i-th point is optimized (1) or not (0)
      INTEGER, INTENT(IN) :: o( Nt )            
      ! rij^2, distance squared between points
      DOUBLE PRECISION, INTENT(IN) :: r2( Nt, Nt ) 
      ! gradients of rij 
      DOUBLE PRECISION, INTENT(IN) :: dr2d1( Nt, Nt, Nt )
      DOUBLE PRECISION, INTENT(IN) :: dr2d2( Nt,Nt,Nt ) 
      ! OUTPUT:
      ! total energy
      DOUBLE PRECISION, INTENT(OUT) :: coul
      DOUBLE PRECISION, INTENT(OUT) :: grad( Nt,2 )
      ! LOCAL VARIABLES:
      INTEGER :: i,j
      DOUBLE PRECISION :: coulj, rij, rinv, rinv3

      coul = 0.0
      grad(:,:) = 0.0

      DO j=1,Nt
         coulj = 0.0
         DO i=1,Nt
            IF (i == j) CYCLE
            rij = SQRT(r2(i,j))
            rinv = 1.0/rij
            coulj = coulj + 0.5*rinv
            ! gradient
            rinv3 = 0.25 * rinv**3
            IF (o(i) == 1) THEN
               grad(i,1) = grad(i,1) - rinv3 * dr2d1(i,j,i)
               grad(i,2) = grad(i,2) - rinv3 * dr2d2(i,j,i)
            ENDIF
            IF (o(j) == 1) THEN
               grad(j,1) = grad(j,1) - rinv3 * dr2d1(i,j,j)
               grad(j,2) = grad(j,2) - rinv3 * dr2d2(i,j,j)
            ENDIF
         ENDDO
         coul = coul + coulj
      ENDDO
    END SUBROUTINE coulomb_energy

!*****************************************************
    SUBROUTINE remove_overlapping(r, Nat, tol, dup)
! Purpose:
! ========
! marks those 3D positions in r that overlap with other positions as duplicates
      ! release global interpreter lock
      !f2py threadsafe
      IMPLICIT NONE
      ! INPUT:      
      INTEGER, INTENT(IN) :: Nat
      DOUBLE PRECISION, INTENT(IN) :: r( Nat, 3 )
      ! maximum separation between positions to be considered distinct
      DOUBLE PRECISION, INTENT(IN) :: tol   
      ! OUTPUT:
      ! dup(i) == 1, means that the i-th position occurs more than once 
      ! and should be removed
      INTEGER, INTENT(OUT) :: dup( Nat )
      ! LOCAL:
      INTEGER :: i,j
      DOUBLE PRECISION :: tol2
      DOUBLE PRECISION :: ri(3), rj(3), rij(3), d2

      tol2 = tol**2

      DO i=1,Nat
         ri = r(i,:)
         DO j=i+1,Nat
            rj = r(j,:)
            rij = ri-rj
            d2 = SUM(rij*rij)
            IF (d2 < tol2) THEN
               dup(j) = 1
            ENDIF
         ENDDO
      ENDDO
    END SUBROUTINE remove_overlapping

END MODULE thomson

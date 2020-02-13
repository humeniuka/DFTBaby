!***************************************************
! implementation of some functions needed for COSMO
! (see DFTB/ImplicitSolvent.py) in Fortran
!***************************************************

SUBROUTINE cosmo_ab_matrices(nats,npts, &
     centers, surface_points, surface_areas, &
     A,B)
  ! Purpose:
  ! ========
  ! compute the matrices A and B needed for constructing
  ! the gamma-matrix due the screening of charges by the
  ! solvent

  ! release global interpreter lock
  !f2py threadsafe
  IMPLICIT NONE
  ! INPUT
  ! number of atoms `npts` and number of surface points `npts`
  ! on the surface of the cavity where the surface charges are induced
  INTEGER, INTENT(IN) :: nats, npts
  ! atomic positions
  DOUBLE PRECISION, INTENT(IN) :: centers(3,nats)
  ! positions on solvent accessible surface
  DOUBLE PRECISION, INTENT(IN) :: surface_points(3,npts)
  ! areas of patches
  DOUBLE PRECISION, INTENT(IN) ::surface_areas(npts)
  ! OUTPUT
  ! matrices
  ! A(k,l) : Coulomb interaction between surface points k and l
  ! B(i,k) : Coulomb interaction between atomic charge i
  !          and surface point k
  DOUBLE PRECISION, INTENT(OUT) :: A(npts,npts), B(nats,npts)
  ! LOCAL VARIABLE
  INTEGER :: k,l  ! enumerates surface points
  INTEGER :: i    ! enumerate atoms
  DOUBLE PRECISION :: tk(3), tl(3)   ! coordinates of surface points
  DOUBLE PRECISION :: ri(3)          ! atomic position
  DOUBLE PRECISION :: dist
  
  ! Coulomb interaction between induced surface charges
  A = 0.d0
  !$omp parallel do private(k,tk,l,tl,dist)
  DO k=1,npts
     ! diagonal corresponds to self-interaction of charges on the
     ! same patch, see eqn. 7b) in Ref. [1]
     A(k,k) = 3.8d0 / surface_areas(k)
     tk = surface_points(:,k)
     DO l=k+1,npts
        ! off-diagonal corresponds to Coulomb interaction of
        ! induced charges at different patches k != l
        tl = surface_points(:,l)
        ! distance between tk and tl
        dist = SQRT(SUM((tk-tl)**2))
        A(k,l) = 1.0/dist
        A(l,k) = A(k,l)
     ENDDO
  ENDDO
  
  ! Coulomb interactions between charges inside cavity (of solute)
  ! and induced charges on the surface of the cavity
  B = 0.d0
  !$omp parallel do private(k,tk,i,ri,dist)
  DO k=1,npts
     tk = surface_points(:,k)
     DO i=1,nats
        ri = centers(:,i)
        ! distance between atomic charge i and induced surface charge k
        dist = SQRT(SUM((ri-tk)**2))
        B(i,k) = 1.0/dist
     ENDDO
  ENDDO

END SUBROUTINE cosmo_ab_matrices

SUBROUTINE cosmo_gamma_gradient(nats,npts, &
     centers, surface_points, parent_atoms, BinvA, f, &
     grad_gamma_solvent)
  ! Purpose
  ! =======
  ! construct the gradient of solvent contribution to the
  ! gamma-matrix
  !
  ! release global interpreter lock
  !f2py threadsafe
  IMPLICIT NONE
  ! INPUT
  ! number of atoms `npts` and number of surface points `npts`
  ! on the surface of the cavity where the surface charges are induced
  INTEGER, INTENT(IN) :: nats, npts
  ! atomic positions
  DOUBLE PRECISION, INTENT(IN) :: centers(3,nats)
  ! positions on solvent accessible surface
  DOUBLE PRECISION, INTENT(IN) :: surface_points(3,npts)
  ! indeces of atoms to which each surface point belongs
  INTEGER, INTENT(IN) :: parent_atoms(npts)
  ! matrix B.A^(-1)
  DOUBLE PRECISION, INTENT(IN) :: BinvA(nats,npts)
  ! scaling factor f=(eps-1)/(eps+x) for scaling the screening energy
  ! from the value for a perfect conductor to a dielectric with
  ! permittivity eps
  DOUBLE PRECISION, INTENT(IN) :: f
  ! OUTPUT
  DOUBLE PRECISION, INTENT(OUT) :: grad_gamma_solvent(3*nats,nats,nats)
  ! LOCAL
  INTEGER :: i,j,k  ! enumerate atoms
  INTEGER :: m,n    ! enumerate surface points
  INTEGER :: xyz    ! enumerate cartesian coordinates xyz=1,2,3
  DOUBLE PRECISION :: tm(3), tn(3) ! coordinates of surface points
  DOUBLE PRECISION :: ri(3)        ! atomic position
  DOUBLE PRECISION :: vec(3)       ! position vector
  DOUBLE PRECISION :: dist
  ! gradient of A matrix w/r/t position of a single atom k
  DOUBLE PRECISION :: gradkA(3,npts,npts)  
  ! gradient of B matrix w/r/t position of a single atom k
  DOUBLE PRECISION :: gradkB(3,nats,npts)

  grad_gamma_solvent = 0.d0
  ! build gradient for each atom position
  !$omp parallel do private(k,gradkA,n,tn,j,m,i,tm,vec,dist,ri,gradkB,xyz)
  DO k=1,nats
     ! construct grad_k A
     ! The surface points tm and tn belong to atoms i and j, respectively.
     !
     !                     tm - tn
     !   grad   A    =  - --------- (delta_{k,i} - delta_{k,j})
     !       rk  m,n      |tm-tn|^3
     !
     gradkA = 0.d0
     DO n=1,npts
        tn = surface_points(:,n)
        j = parent_atoms(n)+1
        DO m=1,npts
           i = parent_atoms(m)+1
           IF (.not.( (i == k).or.(j == k) )) CYCLE
           ! surface points belonging to the same atom exert no force
           ! this atom
           IF (i == j) CYCLE
           tm = surface_points(:,m)
           vec = tm-tn
           dist = SQRT(SUM(vec**2))
           vec = -vec/dist**3
           IF (i == k) THEN
              gradkA(:,m,n) = gradkA(:,m,n) + vec
           ENDIF
           IF (j == k) THEN
              gradkA(:,m,n) = gradkA(:,m,n) - vec
           ENDIF
        ENDDO
     ENDDO
     ! construct grad_k B
     ! Atom i is at position ri, surface point tm belongs to atom j
     !
     !                    tm - ri
     !   grad   B     =  --------- (delta_{k,i} - delta_{k,j})
     !       rk   i,m    |tm-ri|^3
     !
     gradkB = 0.d0
     DO m=1,npts
        tm = surface_points(:,m)
        j = parent_atoms(m)+1
        DO i=1,nats
           IF (.not.( (i == k).or.(j == k)) ) CYCLE
           ! surface point m does not exert force onto its
           ! parent atom
           IF (i == j) CYCLE
           ri = centers(:,i)
           vec = tm-ri
           dist = SQRT(SUM(vec**2))
           vec = vec/dist**3
           IF (i == k) THEN
              gradkB(:,i,m) = gradkB(:,i,m) + vec
           ENDIF
           IF (j == k) THEN
              gradkB(:,i,m) = gradkB(:,i,m) - vec
           ENDIF
        ENDDO
     ENDDO
     ! assemble gradient of gamma_solvent from gradients gradkA and gradkB
     DO xyz=1,3
        grad_gamma_solvent(3*(k-1)+xyz,:,:) = &
              MATMUL(gradkB(xyz,:,:), TRANSPOSE(BinvA)) &
             +MATMUL(BinvA, TRANSPOSE(gradkB(xyz,:,:))) &
             -MATMUL(BinvA, MATMUL(gradkA(xyz,:,:), TRANSPOSE(BinvA)))
     ENDDO
  ENDDO
  ! scale to finite dielectric constant
  grad_gamma_solvent = - f * grad_gamma_solvent
  
END SUBROUTINE cosmo_gamma_gradient


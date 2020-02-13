!***************************************************************
! This file contains code snippets that were too slow in Python
! and therefore have been recoded in Fortran and parallelized
! with OpenMP
!***************************************************************

MODULE tddftb
  IMPLICIT NONE
  CONTAINS

!************************************
    SUBROUTINE trans_charges( &
         orbs_occ, orbs_virt, S, orbsPerAtom, &
         Nat, Norb, nocc, nvirt, &
         q_oo, q_ov, q_vv)
! Purpose:
! ========
! computes the Mulliken transition charges between occupied-occupied
! occupied-virtual and virtual-virtual molecular orbitals.
! Uses point charge approximation of transition densities according to formula (14)
! in Heringer, Niehaus  J Comput Chem 28: 2589-2601 (2007)


      ! release global interpreter lock
      !f2py threadsafe
      IMPLICIT NONE
      ! INPUT:
      INTEGER, INTENT(IN) :: Nat                ! number of atoms
      INTEGER, INTENT(IN) :: Norb               ! total number of atomic orbitals
      INTEGER, INTENT(IN) :: nocc, nvirt        ! number of active occupied and virtual orbitals
      INTEGER, INTENT(IN) :: orbsPerAtom( Nat )  ! number of atomic orbitals per atom
      ! MO coefficients of occupied orbitals, orbs_occ(:,i) contains the coefficients of the i-th occupied MO
      DOUBLE PRECISION, INTENT(IN) :: orbs_occ( Norb, nocc ) 
      ! MO coefficients of virtual orbitals
      DOUBLE PRECISION, INTENT(IN) :: orbs_virt( Norb, nvirt )
      DOUBLE PRECISION, INTENT(IN) :: S(Norb,Norb) ! overlap matrix between atomic orbitals
      ! OUTPUT:
      DOUBLE PRECISION, INTENT(OUT) :: q_oo( Nat, nocc, nocc )   ! transition charges occ-occ
      DOUBLE PRECISION, INTENT(OUT) :: q_ov( Nat, nocc, nvirt )  ! transition charges occ-virt
      DOUBLE PRECISION, INTENT(OUT) :: q_vv( Nat, nvirt, nvirt ) ! transition charges virt-virt
      ! LOCAL VARIABLES:
      INTEGER :: i,j, a,b    ! enumerate occupied and virtual orbitals
      INTEGER :: mu          ! enumerate atomic orbitals
      INTEGER :: muA         ! enumerates atomic orbitals on atom At
      INTEGER :: At          ! enumerates atoms
      DOUBLE PRECISION :: Sc_occ( Norb, nocc )   ! product S.orbs_occ
      DOUBLE PRECISION :: Sc_virt( Norb, nvirt ) ! product S.orbs_virt
      DOUBLE PRECISION :: q_Aij, q_Aia, q_Aab 

      q_oo(:,:,:) = 0.0
      q_ov(:,:,:) = 0.0
      q_vv(:,:,:) = 0.0

      Sc_occ  = MATMUL(S, orbs_occ)
      Sc_virt = MATMUL(S, orbs_virt)

      ! occupied-occupied
      !$omp parallel do private(i,j,mu,At,q_Aij,muA)
      DO i=1,nocc
         DO j=1,nocc
            mu=1
            DO At=1,Nat
               q_Aij = 0.0
               DO muA=1,orbsPerAtom(At)
                  q_Aij =  q_Aij &
                         + orbs_occ(mu,i)*Sc_occ(mu,j) &
                         + orbs_occ(mu,j)*Sc_occ(mu,i)
                  mu = mu+1
               ENDDO
               q_oo(At,i,j) = q_Aij
            ENDDO
         ENDDO
      ENDDO
      q_oo(:,:,:) = 0.5 * q_oo(:,:,:)
      ! occupied-virtual
      !$omp parallel do private(i,a,mu,At,q_Aia,muA)
      DO i=1,nocc
         DO a=1,nvirt
            mu=1
            DO At=1,Nat
               q_Aia = 0.0
               DO muA=1,orbsPerAtom(At)
                  q_Aia =  q_Aia &
                         + orbs_occ(mu,i)*Sc_virt(mu,a) &
                         + orbs_virt(mu,a)*Sc_occ(mu,i)
                  mu = mu+1
               ENDDO
               q_ov(At,i,a) = q_Aia
            ENDDO
         ENDDO
      ENDDO
      q_ov(:,:,:) = 0.5 * q_ov(:,:,:)
      ! virtual-virtual
      !$omp parallel do private(a,b,mu,At,q_Aab,muA)
      DO a=1,nvirt
         DO b=1,nvirt
            mu=1
            DO At=1,Nat
               q_Aab = 0.0
               DO muA=1,orbsPerAtom(At)
                  q_Aab =  q_Aab &
                         + orbs_virt(mu,a)*Sc_virt(mu,b) &
                         + orbs_virt(mu,b)*Sc_virt(mu,a)
                  mu = mu+1
               ENDDO
               q_vv(At,a,b) = q_Aab
            ENDDO
         ENDDO
      ENDDO
      q_vv(:,:,:) = 0.5 * q_vv(:,:,:) 
      
    END SUBROUTINE trans_charges

    SUBROUTINE ApBv(g,g_lr, q_oo, q_ov, q_vv, omega, vs, &
         Nat,nocc,nvirt,nvec, &
         us )
! Purpose:
! ========
! computes the matrix product (A+B).v for all vectors in vs
!
      ! release global interpreter lock
      !f2py threadsafe
      IMPLICIT NONE
      ! INPUT:
      ! gamma matrices with and without long-range correction
      INTEGER, INTENT(IN) :: Nat    ! number of atoms
      INTEGER, INTENT(IN) :: nocc, nvirt  ! number of active occupied and virtuals
      INTEGER, INTENT(IN) :: nvec   ! number of expansion vectors
      DOUBLE PRECISION, INTENT(IN) :: g(Nat,Nat), g_lr(Nat,Nat) 
      DOUBLE PRECISION, INTENT(IN) :: q_oo(Nat,nocc,nocc)   ! transition charges occ-occ
      DOUBLE PRECISION, INTENT(IN) :: q_ov(Nat,nocc,nvirt)  ! transition charges occ-virt
      DOUBLE PRECISION, INTENT(IN) :: q_vv(Nat,nvirt,nvirt) ! transition charges virt-virt
      DOUBLE PRECISION, INTENT(IN) :: omega(nocc,nvirt)     ! orbital energy differences
      DOUBLE PRECISION, INTENT(IN) :: vs(nocc,nvirt,nvec)   ! expansion vectors
      ! OUTPUT:
      DOUBLE PRECISION, INTENT(OUT) :: us(nocc,nvirt,nvec)  ! matrix product (A+B).v
      ! LOCAL VARIABLES:
      INTEGER :: At   ! enumerates atoms
      INTEGER :: l    ! enumerate expansion vectors v
      ! temporary variables, tmpIJ refers to the J-th temporary variable
      ! needed for calculating the I-th term in the product u=(A+B).v
      DOUBLE PRECISION :: ul(nocc,nvirt), vl(nocc,nvirt)
      DOUBLE PRECISION :: tmp21(Nat)   ! sum_(jb) q_B^(jb) v_(jb)
      DOUBLE PRECISION :: tmp22(Nat)
      DOUBLE PRECISION :: tmp31(Nat,nvirt,nocc)
      DOUBLE PRECISION :: tmp32(Nat,nocc,nvirt)
      DOUBLE PRECISION :: tmp33(nocc,nvirt)
      DOUBLE PRECISION :: tmp41(Nat,nvirt,nvirt)
      DOUBLE PRECISION :: tmp42(Nat,nvirt,nvirt)
      DOUBLE PRECISION :: tmp43(nocc,nvirt)

      !$omp parallel do private(l,vl,ul,At), &
      !$omp& private(tmp21,tmp22), &
      !$omp& private(tmp31,tmp32,tmp33), &
      !$omp& private(tmp41,tmp42,tmp43)
      DO l=1,nvec
!         !$omp workshare
         vl = vs(:,:,l)
         ! 1st term - KS orbital energy differences
         ul = omega(:,:)*vl
!         !$omp end workshare
         ! 2nd term - Coulomb
!         !$omp parallel do private(At) reduction(+:tmp21)
         DO At=1,Nat
            tmp21(At) = SUM(q_ov(At,:,:)*vl)
         ENDDO
!         !$omp workshare
         tmp22 = 4*MATMUL(g, tmp21)
!         !$omp end workshare
!         !$omp parallel do private(At) reduction(+:ul)
         DO At=1,Nat
            ul = ul + q_ov(At,:,:)*tmp22(At)
         ENDDO
         ! 3rd term - Exchange
!         !$omp workshare
         tmp31=RESHAPE(&
              MATMUL( &
                RESHAPE(q_vv, (/ Nat*nvirt,nvirt /)), &
                TRANSPOSE(vl)), &
              (/ Nat,nvirt,nocc /))
         tmp32=RESHAPE(RESHAPE(&
              MATMUL(g_lr, RESHAPE(tmp31, (/ Nat,nvirt*nocc /))), &
              (/ Nat,nvirt,nocc /)), (/Nat,nocc,nvirt/), ORDER=(/1,3,2/))
         tmp33=MATMUL(&
              TRANSPOSE(RESHAPE(q_oo, (/ Nat*nocc, nocc /))), &
              RESHAPE(tmp32, (/ Nat*nocc, nvirt /)))
         ul = ul - tmp33         
         ! 4th term - Exchange
         tmp41 = RESHAPE(&
              MATMUL( &
              RESHAPE(RESHAPE(q_ov, (/ Nat,nvirt,nocc /), ORDER=(/ 1,3,2 /)), &
              (/ Nat*nvirt,nocc /)), vl), &
              (/ Nat,nvirt,nvirt /))
         tmp42 = RESHAPE( &
              MATMUL(g_lr, RESHAPE(tmp41, (/ Nat, nvirt*nvirt /))), &
              (/ Nat,nvirt,nvirt /))
         tmp43 = MATMUL( &
              RESHAPE(RESHAPE(q_ov, (/ nocc,Nat,nvirt /), ORDER=(/ 2,1,3 /)), &
               (/ nocc, Nat*nvirt /)), &
              RESHAPE(RESHAPE(tmp42, (/Nat,nvirt,nvirt/), ORDER=(/ 1,3,2 /)), &
              (/ Nat*nvirt, nvirt /)))
         ul = ul - tmp43
         us(:,:,l) = ul
!         !$omp end workshare
      ENDDO
    END SUBROUTINE ApBv
    
    SUBROUTINE AmBv(g,g_lr, q_oo, q_ov, q_vv, omega, vs, &
         Nat,nocc,nvirt,nvec, &
         us )
! Purpose:
! ========
! computes the matrix product (A-B).v for all vectors in vs
!
      ! release global interpreter lock
      !f2py threadsafe
      IMPLICIT NONE
      ! INPUT:
      ! gamma matrices with and without long-range correction
      INTEGER, INTENT(IN) :: Nat    ! number of atoms
      INTEGER, INTENT(IN) :: nocc, nvirt  ! number of active occupied and virtuals
      INTEGER, INTENT(IN) :: nvec   ! number of expansion vectors
      DOUBLE PRECISION, INTENT(IN) :: g(Nat,Nat), g_lr(Nat,Nat) 
      DOUBLE PRECISION, INTENT(IN) :: q_oo(Nat,nocc,nocc)   ! transition charges occ-occ
      DOUBLE PRECISION, INTENT(IN) :: q_ov(Nat,nocc,nvirt)  ! transition charges occ-virt
      DOUBLE PRECISION, INTENT(IN) :: q_vv(Nat,nvirt,nvirt) ! transition charges virt-virt
      DOUBLE PRECISION, INTENT(IN) :: omega(nocc,nvirt)     ! orbital energy differences
      DOUBLE PRECISION, INTENT(IN) :: vs(nocc,nvirt,nvec)   ! expansion vectors
      ! OUTPUT:
      DOUBLE PRECISION, INTENT(OUT) :: us(nocc,nvirt,nvec)  ! matrix product (A-B).v
      ! LOCAL VARIABLES:
      INTEGER :: At   ! enumerates atoms
      INTEGER :: l    ! enumerate expansion vectors v
      ! temporary variables, tmpIJ refers to the J-th temporary variable
      ! needed for calculating the I-th term in the product u=(A+B).v
      DOUBLE PRECISION :: ul(nocc,nvirt), vl(nocc,nvirt)
      DOUBLE PRECISION :: tmp21(Nat,nvirt,nvirt)
      DOUBLE PRECISION :: tmp22(Nat,nvirt,nvirt)
      DOUBLE PRECISION :: tmp23(nocc,nvirt)
      DOUBLE PRECISION :: tmp31(Nat,nvirt,nocc)
      DOUBLE PRECISION :: tmp32(Nat,nvirt,nocc)
      DOUBLE PRECISION :: tmp33(nocc,nvirt)

      !$omp parallel do private(l,vl,ul,tmp21,tmp22,tmp23), &
      !$omp& private(tmp31,tmp32,tmp33) 
      DO l=1,nvec
!         !$omp workshare
         vl = vs(:,:,l)
         ! 1st term, differences in orbital energies
         ul = omega(:,:)*vl
         ! 2nd term, Exchange
         tmp21 = RESHAPE( &
              MATMUL(&
              RESHAPE( &
               RESHAPE(q_ov, (/ Nat, nvirt, nocc/), ORDER=(/ 1,3,2 /)), &
               (/ Nat*nvirt, nocc /)), &
              vl), (/ Nat, nvirt, nvirt /))
         tmp22 = RESHAPE( &
              MATMUL(g_lr, &
              RESHAPE(tmp21, (/ Nat, nvirt*nvirt/))), &
              (/ Nat, nvirt, nvirt /))
         tmp23 = MATMUL( &
              RESHAPE(RESHAPE(q_ov, (/ nocc, Nat, nvirt /), ORDER=(/ 2,1,3 /)), &
              (/ nocc, Nat*nvirt /)), &
              RESHAPE(RESHAPE(tmp22, (/ Nat,nvirt,nvirt /), ORDER=(/ 1,3,2 /)), &
              (/ Nat*nvirt, nvirt /)))
         ul = ul + tmp23
         ! 3rd term, Exchange
         tmp31 = RESHAPE(MATMUL( &
              RESHAPE(q_vv, (/ Nat*nvirt, nvirt /)), &
              TRANSPOSE(vl)), &
              (/ Nat, nvirt, nocc /))
         tmp32 = RESHAPE(MATMUL( &
              g_lr, RESHAPE(tmp31, (/Nat, nvirt*nocc /))), &
              (/ Nat, nvirt,nocc /))
         tmp33 = MATMUL( &
              RESHAPE( &
              RESHAPE(q_oo, (/ nocc, Nat, nocc /), ORDER=(/ 2,1,3 /)), &
                   (/ nocc, Nat*nocc /)), &
              RESHAPE( &
              RESHAPE(tmp32, (/ Nat,nocc,nvirt /), ORDER=(/ 1,3,2 /)), &
                   (/ Nat*nocc, nvirt /)))
         ul = ul - tmp33
         us(:,:,l) = ul
!         !$omp end workshare
      END DO
    END SUBROUTINE AmBv
    
! NOTE: If numpy were parallelized (e.g. when
! compiled with mkl), np.dot and np.tensordot
! should be used instead of the following
! two functions.

    SUBROUTINE dot(v,w, nocc,nvirt, c)
!****************************************
!* compute the dot product between two
!* excitation vectors with shape (nocc,nvirt)
!* Parameters:
!* ===========
!* v,w: numpy arrays (nocc,nvirt)
!* 
!* Returns:
!* ========
!* c: overlap <v|w>
!*
      ! release global interpreter lock
      !f2py threadsafe
      IMPLICIT NONE
      ! INPUT:
      DOUBLE PRECISION, INTENT(IN) :: v(nocc,nvirt)
      DOUBLE PRECISION, INTENT(IN) :: w(nocc,nvirt)
      INTEGER, INTENT(IN) :: nocc, nvirt
      ! OUTPUT:
      DOUBLE PRECISION, INTENT(OUT) :: c
      ! LOCAL VARIABLES:
      DOUBLE PRECISION :: d,ci
      INTEGER :: i,a
      !$omp parallel do private(a,ci,i,d) reduction(+:c)
       DO a=1,nvirt
         ci = 0.0
         DO i=1,nocc
            d = v(i,a)*w(i,a)
            ci = ci+d
         END DO
         c = c+ci
      END DO
    END SUBROUTINE dot

    SUBROUTINE tensordot20(bs, Vb, nocc,nvirt,nn,nm, V)
!******************************************
!* PURPOSE:
!* ========
!* Perform a basis transformation 
!* from Vb (coordinates in basis bs)
!* to V (coordinates in canonical basis)
!*
      ! release global interpreter lock
      !f2py threadsafe
      IMPLICIT NONE
      ! INPUT:
      DOUBLE PRECISION, INTENT(IN) :: bs(nocc,nvirt,nm)
      DOUBLE PRECISION, INTENT(IN) :: Vb(nm,nn)
      INTEGER, INTENT(IN) :: nocc,nvirt
      INTEGER, INTENT(IN) :: nn,nm
      ! OUTPUT:
      DOUBLE PRECISION, INTENT(OUT) :: V(nocc,nvirt,nn)
      ! LOCAL VARIABLES
      INTEGER :: n,m
      DOUBLE PRECISION :: Vn(nocc,nvirt)
      !$omp parallel do private(n,Vn,m)
      DO n=1,nn
         Vn(:,:) = 0.0
         DO m=1,nm
            Vn = Vn + Vb(m,n)*bs(:,:,m)
         END DO
         V(:,:,n) = Vn
      END DO
      
    END SUBROUTINE tensordot20

END MODULE tddftb

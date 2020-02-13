!***************************************************************
! This file contains code snippets that were too slow in Python
! and therefore have been recoded in Fortran
!***************************************************************

MODULE mulliken
  IMPLICIT NONE
  CONTAINS

!************************************
    SUBROUTINE monopoles( &
         P, P0, S, orbsPerAtom, &
         Nat, Norb, &
         q, dq)
! Purpose:
! ========
! compute the (partial) Mulliken charges
!

      ! release global interpreter lock
      !f2py threadsafe
      IMPLICIT NONE
      ! INPUT:
      INTEGER, INTENT(IN) :: Nat                ! number of atoms
      INTEGER, INTENT(IN) :: Norb               ! total number of atomic orbitals
      ! density matrix and reference density matrix
      DOUBLE PRECISION, INTENT(IN) :: P( Norb,Norb ), P0( Norb, Norb ) 
     ! overlap matrix between atomic orbitals
      DOUBLE PRECISION, INTENT(IN) :: S( Norb,Norb ) 
      ! number of atomic orbitals per atom
      INTEGER, INTENT(IN) :: orbsPerAtom( Nat ) 

      ! OUTPUT:
      DOUBLE PRECISION, INTENT(OUT) :: q( Nat )  ! Mulliken charges
      DOUBLE PRECISION, INTENT(OUT) :: dq( Nat ) ! partial Mulliken charges
      ! LOCAL VARIABLES
      INTEGER :: A,B   ! enumerate atoms
      INTEGER :: mu, muA, nu, nuB  ! enumerate orbitals and orbitals on atoms A and B
      DOUBLE PRECISION :: dP( Norb, Norb )   ! P-P0

      dP = P-P0

      ! iterate over atoms A
      mu = 1
      ! WARNING: this loop cannot be parallelized easily because mu is incremented
      ! inside the loop
      DO A=1,Nat
         q(A)  = 0.0
         dq(A) = 0.0
         ! iterate over orbitals on atom A
         DO muA=1,orbsPerAtom(A)
            nu = 1
            ! iterate over atoms B
            DO B=1,Nat
               DO nuB=1,orbsPerAtom(B)
                  !
                   q(A) =  q(A) +  P(mu,nu)*S(mu,nu)
                  dq(A) = dq(A) + dP(mu,nu)*S(mu,nu)
                  !
                  nu = nu+1
               ENDDO
            ENDDO
            mu = mu+1
         ENDDO
      ENDDO
    END SUBROUTINE monopoles

END MODULE mulliken

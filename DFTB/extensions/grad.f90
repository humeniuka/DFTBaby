!********************************************************************
! This module contains some functions needed for assembling
! analytical gradients that were translated from python into Fortran
! and were parallelized with OpenMP.
!********************************************************************

MODULE grad
  IMPLICIT NONE
CONTAINS
  SUBROUTINE Fop(Nat,Norb, &
       S,gradS,g,gradG,v,&
       f)
  ! Purpose:
  ! ========
  ! evaluate the linear operator F_ab[v] 
  !  F_ab[v] = sum_(g,d) d(ab|gd)/dR * v_gd
  ! given in eqn. (90) of excited_state_gradients_Furche.pdf
  
  ! release global interpreter lock
  !f2py threadsafe
  IMPLICIT NONE

  ! INPUT
  ! number of atoms
  INTEGER, INTENT(IN) :: Nat
  ! number of orbitals
  INTEGER, INTENT(IN) :: Norb
  ! overlap matrix
  DOUBLE PRECISION, INTENT(IN) :: S(Norb,Norb)
  ! gradient of overlap matrix 
  DOUBLE PRECISION, INTENT(IN) :: gradS(Norb,Norb,3*Nat)
  ! gamma matrix with AO indeces
  DOUBLE PRECISION, INTENT(IN) :: g(Norb,Norb)
  ! gradient of gamma matrix with AO indeces
  DOUBLE PRECISION, INTENT(IN) :: gradG(Norb,Norb,3*Nat)
  ! vector v
  DOUBLE PRECISION, INTENT(IN) :: v(Norb,Norb)
  ! OUTPUT
  DOUBLE PRECISION, INTENT(OUT) :: f(Norb,Norb,3*Nat)
  ! LOCAL VARIABLES
  INTEGER :: nc  ! enumerates cartesian coordinates
  INTEGER :: a,b ! enumerate orbitals
  DOUBLE PRECISION :: dS(Norb,Norb)  ! one component of the gradient, dS(:,:) = gradS(nc,:,:)
  DOUBLE PRECISION :: dG(Norb,Norb)  ! one component of the gradient, dG(:,:) = gradG(nc,:,:)
  DOUBLE PRECISION :: df(Norb,Norb)  ! one component of results df(:,:) = f(nc,:,:)
  !
  DOUBLE PRECISION :: vp(Norb,Norb), Sv(Norb), gSv(Norb), gdSv(Norb), dgSv(Norb), tmp(Norb)
  
  vp = v+TRANSPOSE(v)
  Sv = SUM(S*vp, 2)
  gSv = MATMUL(g,Sv)
  !$omp parallel do private(nc,dS,dG,gdSv,dgSv,tmp,df,b,a)
  DO nc=1,3*Nat
     ! create temporary variables for component nc 
     dS = gradS(:,:,nc)
     dG = gradG(:,:,nc)
     !
     !!gdSv = MATMUL(g,SUM(dS*vp,2))
     tmp = SUM(dS*vp,2)
     call DGEMV('n',Norb,Norb,1.d0, g, Norb, tmp, 1, 0.d0, gdSv, 1)
     !!dgSv = MATMUL(dG,Sv)
     call DGEMV('n',Norb,Norb,1.d0, dG, Norb, Sv, 1, 0.d0, dgSv, 1)
     df(:,:) = 0.0
     DO b=1,Norb
        DO a=1,Norb
           df(a,b) = dS(a,b)*(gSv(a)+gSv(b)) &
                     +S(a,b)*(dgSv(a)+gdSv(a)+dgSv(b)+gdSv(b))
        ENDDO
     ENDDO
     ! assign result for component nc
     df = df*0.25
     f(:,:,nc) = df
  ENDDO
END SUBROUTINE Fop

SUBROUTINE Flrop(Nat,Norb, &
       S,gradS,g,gradG,v,&
       f)
  ! Purpose:
  ! ========!
  ! evaluate the linear operator F^lr_ab[v] 
  !  F_ab[v] = sum_(g,d) d(ab|gd)_lr/dR * v_gd
  ! given in eqn. (91) of excited_state_gradients_Furche.pdf
  
  ! release global interpreter lock
  !f2py threadsafe
  IMPLICIT NONE

  ! INPUT
  ! number of atoms
  INTEGER, INTENT(IN) :: Nat
  ! number of orbitals
  INTEGER, INTENT(IN) :: Norb
  ! overlap matrix
  DOUBLE PRECISION, INTENT(IN) :: S(Norb,Norb)
  ! gradient of overlap matrix
  DOUBLE PRECISION, INTENT(IN) :: gradS(Norb,Norb,3*Nat)
  ! long-range gamma matrix with AO indeces
  DOUBLE PRECISION, INTENT(IN) :: g(Norb,Norb)
  ! gradient of long-range gamma matrix with AO indeces
  DOUBLE PRECISION, INTENT(IN) :: gradG(Norb,Norb,3*Nat)
  ! vector v
  DOUBLE PRECISION, INTENT(IN) :: v(Norb,Norb)
  ! OUTPUT
  DOUBLE PRECISION, INTENT(OUT) :: f(Norb,Norb,3*Nat)
  ! LOCAL VARIABLES
  INTEGER :: nc  ! enumerates cartesian coordinates
  DOUBLE PRECISION :: dS(Norb,Norb)  ! one component of the gradient, dS(:,:) = gradS(:,:,nc)
  DOUBLE PRECISION :: dG(Norb,Norb)  ! one component of the gradient, dG(:,:) = gradG(:,:,nc)
  DOUBLE PRECISION :: df(Norb,Norb)  ! one component of results df(:,:) = f(:,:,nc)
  !
  DOUBLE PRECISION :: Sv(Norb,Norb), SvT(Norb,Norb), gv(Norb,Norb), vT(Norb,Norb)
  DOUBLE PRECISION :: dSvT(Norb,Norb), dSv(Norb,Norb), dgv(Norb,Norb)
  DOUBLE PRECISION :: tSv(Norb,Norb), tSvg(Norb,Norb), tSgv(Norb,Norb)
  DOUBLE PRECISION :: tmp(Norb,Norb)
  
  Sv = MATMUL(S,v)
  vT = TRANSPOSE(v)
  SvT = MATMUL(S,vT)
  gv = g*v
  !
  tSv = TRANSPOSE(Sv)
  tSvg = TRANSPOSE(Sv*g)
  tSgv = TRANSPOSE(MATMUL(S,gv))
  !$omp parallel do private(nc,dS,dG,dSvT,dSv,dgv,df,tmp)
  DO nc=1,3*Nat
     ! create temporary variables for component nc 
     dS(:,:) = gradS(:,:,nc)
     dG(:,:) = gradG(:,:,nc)
     !
     !!dSvT = MATMUL(dS,vT)
     call DGEMM('n','n',Norb,Norb,Norb,1.d0, dS, Norb, vT, Norb, 0.d0, dSvT, Norb)
     !!dSv = MATMUL(dS,v)
     call DGEMM('n','n',Norb,Norb,Norb,1.d0, dS, Norb, v, Norb, 0.d0, dSv, Norb)
     dgv = dG*v
     df(:,:) = 0.0
     ! 1st term
     !!df = df + g*MATMUL(dS, tSv)
     call DGEMM('n','n',Norb,Norb,Norb,1.d0, dS, Norb, tSv, Norb, 0.d0, tmp, Norb)
     df = df + g*tmp
     ! 2nd term
     !!df = df + MATMUL(dSvT*g, S)
     call DGEMM('n','n',Norb,Norb,Norb,1.d0, dSvT*g, Norb, S, Norb, 1.d0, df, Norb)
     ! 3rd term
     !!df = df + MATMUL(dS, tSvg)
     call DGEMM('n','n',Norb,Norb,Norb,1.d0, dS, Norb, tSvg, Norb, 1.d0, df, Norb)
     ! 4th term
     !!df = df + MATMUL(dS, tSgv)
     call DGEMM('n','n',Norb,Norb,Norb,1.d0, dS, Norb, tSgv, Norb, 1.d0, df, Norb)
     ! 5th term
     !!df = df + g*MATMUL(S, TRANSPOSE(dSv))
     call DGEMM('n','t',Norb,Norb,Norb,1.d0, S, Norb, dSv, Norb, 0.d0, tmp, Norb)
     df = df + g*tmp
     ! 6th term
     !!df = df + MATMUL(SvT*g, TRANSPOSE(dS))
     call DGEMM('n','t',Norb,Norb,Norb,1.d0, SvT*g, Norb, dS, Norb, 1.d0, df, Norb)
     ! 7th term
     !!df = df + MATMUL(S, TRANSPOSE(dSv*g))
     call DGEMM('n','t',Norb,Norb,Norb,1.d0, S, Norb, dSv*g, Norb, 1.d0, df, Norb)     
     ! 8th term
     !!df = df + MATMUL(S, TRANSPOSE(MATMUL(dS,gv)))
     call DGEMM('n','n',Norb,Norb,Norb,1.d0, dS, Norb, gv, Norb, 0.d0, tmp, Norb)
     call DGEMM('n','t',Norb,Norb,Norb,1.d0, S, Norb, tmp, Norb, 1.d0, df, Norb)     
     ! 9th term
     !!df = df + dG*MATMUL(S,tSv)
     call DGEMM('n','n',Norb,Norb,Norb,1.d0, S, Norb, tSv, Norb, 0.d0, tmp, Norb)
     df = df + dG*tmp
     ! 10th term
     !!df = df + MATMUL(SvT*dG , S)
     call DGEMM('n','n',Norb,Norb,Norb,1.d0, SvT*dG, Norb, S, Norb, 1.d0, df, Norb)
     ! 11th term
     !!df = df + MATMUL(S, TRANSPOSE(Sv*dG))
     call DGEMM('n','t',Norb,Norb,Norb,1.d0, S, Norb, Sv*dG, Norb, 1.d0, df, Norb)
     ! 12th term
     !!df = df + MATMUL(S, TRANSPOSE(MATMUL(S, dgv)))
     call DGEMM('n','n',Norb,Norb,Norb,1.d0, S, Norb, dgv, Norb, 0.d0, tmp, Norb)
     call DGEMM('n','t',Norb,Norb,Norb,1.d0, S, Norb, tmp, Norb, 1.d0, df, Norb)
     ! assign result for component nc
     df = df*0.25
     f(:,:,nc) = df
  ENDDO
END SUBROUTINE Flrop

END MODULE grad

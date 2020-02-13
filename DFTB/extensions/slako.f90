!**********************************************************************
! This module contains functions related to Slater-Koster rules
! that were translated into Fortran from python to gain some
! speed and to exploit OpenMP parallelization
!**********************************************************************

MODULE slako
  IMPLICIT NONE
  INTEGER, PARAMETER :: lmax=2  ! dipole rules are available for s,p and d-orbitals
  CONTAINS

  !**********************************************************
  ! The code for the Slater-Koster rules has been generated
  ! automatically using a Mathematica script and is included
  ! here.
  !**********************************************************
  INCLUDE "slako_transformations_lmax2.f90"        ! rules for S and H0
  INCLUDE "slako_transformations_grad_lmax2.f90"   ! rules for gradients of S and H0
  !INCLUDE "slako_transformations_dipole_lmax1.f90"   ! rules for dipole matrix elements
  INCLUDE "slako_transformations_dipole_lmax2.f90"   ! rules for dipole matrix elements
  
  !****************************************************************************************
  ! B-spline interpolation is done using the modified routines from FITPACK
  ! developed by Paul Dierckx.
  ! Reference:
  !   Paul Dierckx, Curve and Surface Fitting with Splines, Oxford University Press, 1993
  !****************************************************************************************
  INCLUDE "splines.f90"
  
  !****************************************************************
SUBROUTINE directional_cosines(pos1,pos2,Nat1,Nat2,  rij,x,y,z)
  ! Purpose:
  ! ========
  ! computes the distances r, and the directional cosines for the atomic
  ! positions in pos1 and pos2
  
  ! release global interpreter lock
  !f2py threadsafe
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: pos1(3,Nat1), pos2(3,Nat2)
  INTEGER, INTENT(IN) :: Nat1, Nat2
  DOUBLE PRECISION, INTENT(OUT) :: rij(Nat1,Nat2), x(Nat1,Nat2), y(Nat1,Nat2), z(Nat1,Nat2)
  ! local variables
  INTEGER :: i,j
  DOUBLE PRECISION :: r,xc,yc,zc
  
  !$omp parallel do private(i,j,xc,yc,zc,r)
  DO j=1,Nat2
     DO i=1,Nat1
        xc=pos2(1,i)-pos1(1,j)
        yc=pos2(2,i)-pos1(2,j)
        zc=pos2(3,i)-pos1(3,j)
        r=SQRT(xc*xc+yc*yc+zc*zc)
        ! distance
        rij(i,j)=r
        ! directional cosine
        IF (r > 0.0) THEN
           x(i,j)=xc/r
           y(i,j)=yc/r
           z(i,j)=zc/r
        ELSE
           x(i,j)=0.0
           y(i,j)=0.0
           z(i,j)=1.0
        ENDIF
        
     ENDDO
  ENDDO
END SUBROUTINE directional_cosines

!*****************************************************************
! computes the matrix elements of the reference Hamiltonian H0
! and the overlap matrix between valence orbitals using Slater Koster rules
!*****************************************************************

SUBROUTINE H0andS(&
     Norb,Nat,&
     ind1,typ1,ls1,ms1,&
     ind2,typ2,ls2,ms2,&
     r,x,y,z,&
     Mprox,&
     orben,&
     deg,tab_filled_SH0,&
     S_knots,S_coefs,H_knots,H_coefs,&
     Npair,Nsk,Npts,&
     S,H0)
  ! release global interpreter lock
  !f2py threadsafe
  IMPLICIT NONE

  ! INPUT
  ! number of atoms in geometries 1 and 2
  INTEGER, INTENT(IN) :: Nat
  ! number of valence orbitals in geometries 1 and 2
  INTEGER, INTENT(IN) :: Norb
  ! indeces of the atoms belonging to the orbitals
  INTEGER, INTENT(IN) :: ind1(Norb),ind2(Norb)
  ! atom types of the atoms belonging to the orbitals
  INTEGER, INTENT(IN) :: typ1(Norb),typ2(Norb)
  ! quantum numbers l,m of the orbitals
  INTEGER, INTENT(IN) :: ls1(Norb),ms1(Norb), ls2(Norb),ms2(Norb)
  ! distances and directional cosines between atoms in geometries 1 and 2
  DOUBLE PRECISION, INTENT(IN) :: r(Nat,Nat), &
       x(Nat,Nat), y(Nat,Nat), z(Nat,Nat)
  ! proximity matrix M, if M[i,j] == 1, atoms i and j are close enough
  ! for the orbitals on them to interact
  INTEGER, INTENT(IN) :: Mprox(Nat,Nat)
  ! energies of orbitals, taken from a full DFT calculation
  DOUBLE PRECISION, INTENT(IN) :: orben(Norb)
  ! degree of the B-splines
  INTEGER, INTENT(IN) :: deg
  ! If the SK table t of a atom pair with index ipair contains data, the flag
  ! tab_filled_SH0[ipair,t] is set to 1
  INTEGER, INTENT(IN) :: tab_filled_SH0(Npair,Nsk)
  ! spline data with knots and coefficients for S and H0
  INTEGER, INTENT(IN) :: Npair ! number of unique atom pairs
  INTEGER, INTENT(IN) :: Nsk   ! maximum number of SK tables (for some pairs, like h-h, most tables are empty)
  INTEGER, INTENT(IN) :: Npts  ! number of knots in a B-spline
  !
  DOUBLE PRECISION, INTENT(IN) :: &
       S_knots(Npair,Nsk,Npts), S_coefs(Npair,Nsk,Npts),&
       H_knots(Npair,Nsk,Npts), H_coefs(Npair,Nsk,Npts)

  ! OUTPUT
  ! matrix elements of overlap S and reference Hamiltonian H0
  DOUBLE PRECISION, INTENT(OUT) :: S(Norb,Norb), H0(Norb,Norb)
       
  ! LOCAL VARIABLES
  INTEGER :: mu,nu  ! enumerate orbitals
  INTEGER :: l1,m1,l2,m2 ! quantum numbers of orbitals mu and nu
  INTEGER :: i,j,ip1,jp1    ! enumerate atom indeces
  INTEGER :: typi,typj ! atom types
  INTEGER :: typMin,typMax ! ordered atom types
  INTEGER :: tab    ! enumerate SK tables
  INTEGER :: ipair  ! encodes the atom pair 
  DOUBLE PRECISION :: Sinterp(Nsk), Hinterp(Nsk)  ! values of interpolated SK tables
  DOUBLE PRECISION :: rij,xij,yij,zij  ! r(i,j), x(i,j), etc.
  INTEGER :: ier ! return code of FITPACK routines

  !$omp parallel do private(nu,j,typj,l2,m2,mu,i,typi,l1,m1), &
  !$omp& private(ip1,jp1,rij,xij,yij,zij,typMin,typMax,ipair), &
  !$omp& private(tab,Sinterp,Hinterp,ier) schedule(dynamic,1)
  DO nu=1,Norb
     ! Which atom does orbital mu sit on?
     j=ind2(nu)
     ! What's the type of this atom?
     typj=typ2(nu)
     ! What are the quantum numbers of orbital mu
     l2=ls2(nu)
     m2=ms2(nu)
     DO mu=nu,Norb
        ! Which atom does orbital mu sit on?
        i=ind1(mu)
        ! What's the type of this atom?
        typi=typ1(mu)
        ! What are the quantum numbers of orbital mu
        l1=ls1(mu)
        m1=ms1(mu)
        ! If the atoms are too far apart, skip the calculation
        ip1=i+1 ! Fortran indeces start at 1
        jp1=j+1
        IF (Mprox(ip1,jp1) == 0) CYCLE
        ! One triangle of the symmetric matrices S and H0 is filled
        IF (mu == nu) THEN
           ! orbitals are normalized to 1
           S(mu,nu) = 1.0
           ! use the true single particle orbital energies for the diagonal elements
           ! This assignment is questionable if the two molecular geometries are
           ! different.
           H0(mu,nu) = orben(mu)
        ELSE
           ! mu > nu
           rij=r(jp1,ip1)
           xij=x(jp1,ip1)
           yij=y(jp1,ip1)
           zij=z(jp1,ip1)

           IF (i /= j) THEN
              ! interpolate Slater-Koster tables
              ! find the chunk of the tables belonging to atom pair (i,j)
              ! When ipair is calculated from (typi,typj) the tuple has to be ordered
              typMin=min(typi,typj)
              typMax=max(typi,typj)
              
              ipair=(typMax*(typMax+1))/2+typMin + 1 ! Fortran indeces start at 1
              DO tab=1,Nsk
                 IF (tab_filled_SH0(ipair,tab) == 1) THEN
                    ! this table has useful data, interpolate it
                    !
                    ! assuming the knots of the B-spline are uniformly spaced
                    call splev_uniform(S_knots(ipair,tab,:),Npts,S_coefs(ipair,tab,:),deg,&
                         rij, Sinterp(tab), ier)
                    call splev_uniform(H_knots(ipair,tab,:),Npts,H_coefs(ipair,tab,:),deg,&
                         rij, Hinterp(tab), ier)
                 ENDIF
              ENDDO
    
           ENDIF
           IF (typi <= typj) THEN
              !
              IF (i == j) THEN
                 ! different orbitals on the same atom should be orthogonal
                 ! Here mu != nu
                 S(mu,nu) = 0.0
                 H0(mu,nu) = 0.0
              ELSE
                 S(mu,nu) = slako_transformation(rij,xij,yij,zij,Sinterp,Nsk,&
                      l1,m1,l2,m2)
                 H0(mu,nu) = slako_transformation(rij,xij,yij,zij,Hinterp,Nsk,&
                      l1,m1,l2,m2)
              ENDIF
           ELSE
              ! swap atoms if typj > typi
              S(mu,nu) = slako_transformation(rij,-xij,-yij,-zij,Sinterp,Nsk,&
                   l2,m2,l1,m1)
              H0(mu,nu) = slako_transformation(rij,-xij,-yij,-zij,Hinterp,Nsk,&
                   l2,m2,l1,m1)
           ENDIF
        ENDIF
     ENDDO
  ENDDO
  ! make matrix symmetric by filling in the other triangle
  !$omp parallel do private(nu,mu) schedule(dynamic,1)
  DO nu=1,Norb
     DO mu=1,nu-1
        S(mu,nu) = S(nu,mu)
        H0(mu,nu) = H0(nu,mu)
     ENDDO
  ENDDO
END SUBROUTINE H0andS

!*****************************************************************
! computes the matrix elements of the dipole operator between
! valence orbitals using Slater Koster rules
!*****************************************************************

SUBROUTINE DipoleMatrix(&
     Norb,Nat,&
     ind1,typ1,ls1,ms1,&
     ind2,typ2,ls2,ms2,&
     r,x,y,z,pos,&
     Mprox,&
     S,&
     deg,tab_filled_D,&
     D_knots,D_coefs,&
     Npair,NskD,Npts,&
     Dipole)
  ! release global interpreter lock
  !f2py threadsafe
  IMPLICIT NONE

  ! INPUT
  ! number of atoms in geometries 1 and 2
  INTEGER, INTENT(IN) :: Nat
  ! number of valence orbitals in geometries 1 and 2
  INTEGER, INTENT(IN) :: Norb
  ! indeces of the atoms belonging to the orbitals
  INTEGER, INTENT(IN) :: ind1(Norb),ind2(Norb)
  ! atom types of the atoms belonging to the orbitals
  INTEGER, INTENT(IN) :: typ1(Norb),typ2(Norb)
  ! quantum numbers l,m of the orbitals
  INTEGER, INTENT(IN) :: ls1(Norb),ms1(Norb), ls2(Norb),ms2(Norb)
  ! distances and directional cosines between atoms in geometries 1 and 2
  DOUBLE PRECISION, INTENT(IN) :: r(Nat,Nat), &
       x(Nat,Nat), y(Nat,Nat), z(Nat,Nat)
  ! position vectors of atoms
  DOUBLE PRECISION, INTENT(IN) :: pos(3,Nat)
  ! proximity matrix M, if M[i,j] == 1, atoms i and j are close enough
  ! for the orbitals on them to interact
  INTEGER, INTENT(IN) :: Mprox(Nat,Nat)
  ! overlap matrix
  DOUBLE PRECISION, INTENT(IN) :: S(Norb,Norb)
  ! degree of the B-splines
  INTEGER, INTENT(IN) :: deg
  ! If the SK table t of a atom pair with index ipair contains data, the flag
  ! tab_filled_D[ipair,t] is set to 1
  INTEGER, INTENT(IN) :: tab_filled_D(Npair,NskD)
  ! spline data with knots and coefficients for S and H0
  INTEGER, INTENT(IN) :: Npair ! number of unique atom pairs
  INTEGER, INTENT(IN) :: NskD  ! maximum number of SK tables for D (for some pairs, like h-h, most tables are empty)
  INTEGER, INTENT(IN) :: Npts  ! number of knots in a B-spline
  !
  DOUBLE PRECISION, INTENT(IN) :: &
       D_knots(Npair,NskD,Npts), D_coefs(Npair,NskD,Npts)

  ! OUTPUT
  ! matrix elements of dipole operator
  DOUBLE PRECISION, INTENT(OUT) :: Dipole(3,Norb,Norb)
       
  ! LOCAL VARIABLES
  INTEGER :: mu,nu  ! enumerate orbitals
  INTEGER :: l1,m1,l2,m2 ! quantum numbers of orbitals mu and nu
  INTEGER :: i,j,ip1,jp1    ! enumerate atom indeces
  INTEGER :: typi,typj ! atom types
  INTEGER :: typMin,typMax ! ordered atom types
  INTEGER :: tab    ! enumerate SK tables
  INTEGER :: ipair  ! encodes the atom pair 
  DOUBLE PRECISION :: Dinterp(NskD)  ! values of interpolated SK tables
  DOUBLE PRECISION :: dip(3) ! dipole vector
  DOUBLE PRECISION :: rij,xij,yij,zij  ! r(i,j), x(i,j), etc.
  INTEGER :: ier ! return code of FITPACK routines, not used

  !$omp parallel do private(nu,j,typj,l2,m2,mu,i,typi,l1,m1), &
  !$omp& private(ip1,jp1,rij,xij,yij,zij,typMin,typMax,ipair), &
  !$omp& private(tab,Dinterp,ier,dip) schedule(dynamic,1)
  DO nu=1,Norb
     ! Which atom does orbital mu sit on?
     j=ind2(nu)
     ! What's the type of this atom?
     typj=typ2(nu)
     ! What are the quantum numbers of orbital mu
     l2=ls2(nu)
     !
     IF (l2 > lmax) CYCLE  ! dipole matrix elements between orbitals with l > 1 are not implemented
     m2=ms2(nu)
     DO mu=nu,Norb
        ! Which atom does orbital mu sit on?
        i=ind1(mu)
        ! What's the type of this atom?
        typi=typ1(mu)
        ! What are the quantum numbers of orbital mu
        l1=ls1(mu)
        IF (l1 > lmax) CYCLE  ! dipole matrix elements between orbitals with l > 1 are not implemented
        m1=ms1(mu)
        ! If the atoms are too far apart, skip the calculation
        ip1=i+1 ! Fortran indeces start at 1
        jp1=j+1
        IF (Mprox(ip1,jp1) == 0) CYCLE
        ! One triangle of the symmetric matrix Dipole(3,:,:) is filled
        ! mu >= nu
        rij=r(jp1,ip1)
        xij=x(jp1,ip1)
        yij=y(jp1,ip1)
        zij=z(jp1,ip1)

        ! HACK:
        ! The Slater-Koster rules for dipole matrix elements between d-orbitals
        ! contain divisions by (x**2+y**2), which may raise a division by zero error.
        ! To avoid this, x and y are shifted slightly away from 0 
        IF ((l1 > 1) .OR. (l2 > 1)) THEN
           IF (xij**2+yij**2 < 1.0e-20) THEN
              xij = xij + 1.0e-20
              yij = yij + 1.0e-20
           ENDIF
        ENDIF
          
        ! interpolate Slater-Koster tables
        ! find the chunk of the tables belonging to atom pair (i,j)
        ! When ipair is calculated from (typi,typj) the tuple has to be ordered
        typMin=min(typi,typj)
        typMax=max(typi,typj)
        
        ipair=(typMax*(typMax+1))/2+typMin + 1 ! Fortran indeces start at 1
        DO tab=1,NskD
           IF (tab_filled_D(ipair,tab) == 1) THEN
              ! this table has useful data, interpolate it
              !
              ! assuming the knots of the B-spline are uniformly spaced
              call splev_uniform(D_knots(ipair,tab,:),Npts,D_coefs(ipair,tab,:),deg,&
                   rij, Dinterp(tab), ier)
           ENDIF
        ENDDO
        IF (typi <= typj) THEN
           ! dipole
           dip(1) = slako_transformation_dipole(rij,xij,yij,zij,Dinterp,NskD,&
                l1,m1,+1,l2,m2)
           dip(2) = slako_transformation_dipole(rij,xij,yij,zij,Dinterp,NskD,&
                l1,m1,-1,l2,m2)
           dip(3) = slako_transformation_dipole(rij,xij,yij,zij,Dinterp,NskD,&
                l1,m1,0,l2,m2)
           ! <nu(r-R1)|r|mu(r-R2)> = <nu(r)|r|mu(r-(R2-R1))> + R1*<nu(r)|mu(r-(R2-R1))> = D + R1*S
           Dipole(:,mu,nu) = dip(:) + pos(:,jp1)*S(mu,nu)
        ELSE
           ! swap atoms if typj > typi
           dip(1) = slako_transformation_dipole(rij,-xij,-yij,-zij,Dinterp,NskD,&
                l2,m2,+1,l1,m1)
           dip(2) = slako_transformation_dipole(rij,-xij,-yij,-zij,Dinterp,NskD,&
                l2,m2,-1,l1,m1)
           dip(3) = slako_transformation_dipole(rij,-xij,-yij,-zij,Dinterp,NskD,&
                l2,m2,0,l1,m1)
           ! <nu(r-R1)|r|mu(r-R2)> = <nu(r)|r|mu(r-(R2-R1))> + R1*<nu(r)|mu(r-(R2-R1))> = D + R1*S
           Dipole(:,mu,nu) = dip(:) + pos(:,ip1)*S(mu,nu)           
        ENDIF
     ENDDO
  ENDDO
  ! make matrix symmetric by filling in the other triangle
  !$omp parallel do private(nu,mu) schedule(dynamic,1)
  DO nu=1,Norb
     DO mu=1,nu-1
        Dipole(:,mu,nu) = Dipole(:,nu,mu)
     ENDDO
  ENDDO
END SUBROUTINE DipoleMatrix

!***** GRADIENTS OF MATRIX ELEMENTS **********
!*****************************************************************
! computes the GRADIENTS of the matrix elements of the reference Hamiltonian H0
! and the overlap matrix between valence orbitals using Slater Koster rules
!*****************************************************************

SUBROUTINE gradients_H0andS(&
     Norb,Nat,&
     ind,typ,ls,ms,&
     r,x,y,z,&
     Mprox,&
     deg,tab_filled_SH0,&
     S_knots,S_coefs,H_knots,H_coefs,&
     Npair,Nsk,Npts,&
     gradS,gradH0)
  ! release global interpreter lock
  !f2py threadsafe
  IMPLICIT NONE

  ! INPUT
  ! number of atoms in geometry
  INTEGER, INTENT(IN) :: Nat
  ! number of valence orbitals
  INTEGER, INTENT(IN) :: Norb
  ! indeces of the atoms belonging to the orbitals
  INTEGER, INTENT(IN) :: ind(Norb)
  ! atom types of the atoms belonging to the orbitals
  INTEGER, INTENT(IN) :: typ(Norb)
  ! quantum numbers l,m of the orbitals
  INTEGER, INTENT(IN) :: ls(Norb),ms(Norb)
  ! distances and directional cosines between atoms
  DOUBLE PRECISION, INTENT(IN) :: r(Nat,Nat), &
       x(Nat,Nat), y(Nat,Nat), z(Nat,Nat)
  ! proximity matrix M, if M[i,j] == 1, atoms i and j are close enough
  ! for the orbitals on them to interact
  INTEGER, INTENT(IN) :: Mprox(Nat,Nat)
  ! degree of the B-splines
  INTEGER, INTENT(IN) :: deg
  ! If the SK table t of a atom pair with index ipair contains data, the flag
  ! tab_filled_SH0[ipair,t] is set to 1
  INTEGER, INTENT(IN) :: tab_filled_SH0(Npair,Nsk)
  ! spline data with knots and coefficients for S and H0
  INTEGER, INTENT(IN) :: Npair ! number of unique atom pairs
  INTEGER, INTENT(IN) :: Nsk   ! maximum number of SK tables (for some pairs, like h-h, most tables are empty)
  INTEGER, INTENT(IN) :: Npts  ! number of knots in a B-spline
  !
  DOUBLE PRECISION, INTENT(IN) :: &
       S_knots(Npair,Nsk,Npts), S_coefs(Npair,Nsk,Npts),&
       H_knots(Npair,Nsk,Npts), H_coefs(Npair,Nsk,Npts)

  ! OUTPUT
  ! gradients matrix elements of overlap S and reference Hamiltonian H0
  DOUBLE PRECISION, INTENT(OUT) :: gradS(3*Nat,Norb,Norb), gradH0(3*Nat,Norb,Norb)
       
  ! LOCAL VARIABLES
  INTEGER :: mu,nu  ! enumerate orbitals
  INTEGER :: l1,m1,l2,m2 ! quantum numbers of orbitals mu and nu
  INTEGER :: i,j,ip1,jp1    ! enumerate atom indeces
  INTEGER :: c3b,c3e  ! indeces where the 3D coordinate vector of atom j starts and ends
  INTEGER :: typi,typj ! atom types
  INTEGER :: typMin,typMax ! ordered atom types
  INTEGER :: tab    ! enumerate SK tables
  INTEGER :: ipair  ! encodes the atom pair 
  DOUBLE PRECISION :: Sinterp(Nsk), Hinterp(Nsk)  ! interpolated values of SK tables
  DOUBLE PRECISION :: Sinterp_deriv(Nsk), Hinterp_deriv(Nsk)  ! interpolated derivatives of SK tables
  DOUBLE PRECISION :: s_deriv(3), h0_deriv(3)  ! gradient contribution from atom pair
  DOUBLE PRECISION :: rij,xij,yij,zij  ! r(i,j), x(i,j), etc.
  INTEGER :: ier ! return code of FITPACK routines

  gradS(:,:,:) = 0.0
  gradH0(:,:,:) = 0.0
  !$omp parallel do private(nu,j,typj,l2,m2,c3b,c3e,mu,i,typi,l1,m1), &
  !$omp& private(ip1,jp1,rij,xij,yij,zij,typMin,typMax,ipair), &
  !$omp& private(tab,Sinterp,Hinterp,Sinterp_deriv,Hinterp_deriv), &
  !$omp& private(s_deriv,h0_deriv,ier) schedule(dynamic,1)
  DO nu=1,Norb
     ! Which atom does orbital mu sit on?
     j=ind(nu)
     ! What's the type of this atom?
     typj=typ(nu)
     ! What are the quantum numbers of orbital mu
     l2=ls(nu)
     m2=ms(nu)
     ! coordinates x,y,z of atom j start at c3b and end at c3e
     c3b=3*j+1
     c3e=3*(j+1)
     DO mu=1,Norb
        ! Which atom does orbital mu sit on?
        i=ind(mu)
        ! What's the type of this atom?
        typi=typ(mu)
        ! What are the quantum numbers of orbital mu
        l1=ls(mu)
        m1=ms(mu)
        ! If the atoms are too far apart, skip the calculation
        ip1=i+1 ! Fortran indeces start at 1
        jp1=j+1
        IF (Mprox(ip1,jp1) == 0) CYCLE
        ! One triangle of the symmetric matrices S and H0 is filled
        IF (mu == nu) THEN
           ! gradient of constant S(mu,mu) = 1 and H0(mu,mu) = orbe(mu) is 0
           ! NO CONTRIBUTION
        ELSE
           ! mu > nu
           rij=r(jp1,ip1)
           xij=x(jp1,ip1)
           yij=y(jp1,ip1)
           zij=z(jp1,ip1)

           IF (i /= j) THEN
              ! interpolate Slater-Koster tables
              ! find the chunk of the tables belonging to atom pair (i,j)
              ! When ipair is calculated from (typi,typj) the tuple has to be ordered
              typMin=min(typi,typj)
              typMax=max(typi,typj)
              
              ipair=(typMax*(typMax+1))/2+typMin + 1 ! Fortran indeces start at 1
              DO tab=1,Nsk
                 IF (tab_filled_SH0(ipair,tab) == 1) THEN
                    ! this table has useful data, interpolate it
                    !
                    ! assuming the knots of the B-spline are uniformly spaced
                    call splev_uniform(S_knots(ipair,tab,:),Npts,S_coefs(ipair,tab,:),deg,&
                         rij, Sinterp(tab), ier)
                    call splev_uniform(H_knots(ipair,tab,:),Npts,H_coefs(ipair,tab,:),deg,&
                         rij, Hinterp(tab), ier)
                    ! find derivatives of B-splines
                    call splev_deriv_uniform(S_knots(ipair,tab,:),Npts,S_coefs(ipair,tab,:),deg,&
                         rij, Sinterp_deriv(tab), ier)
                    call splev_deriv_uniform(H_knots(ipair,tab,:),Npts,H_coefs(ipair,tab,:),deg,&
                         rij, Hinterp_deriv(tab), ier)                    
                 ENDIF
              ENDDO
           ENDIF
           IF (typi <= typj) THEN
              !
              IF (i == j) THEN
                 ! different orbitals on the same atom
                 ! Here mu != nu
                 ! NO CONTRIBUTION
                 s_deriv(:) = 0.0
                 h0_deriv(:) = 0.0
              ELSE
                 s_deriv(:) = -slako_transformation_gradient(&
                      rij,xij,yij,zij,&
                      Sinterp,Sinterp_deriv,Nsk,&
                      l1,m1,l2,m2)
                 h0_deriv(:) = -slako_transformation_gradient(&
                      rij,xij,yij,zij,&
                      Hinterp,Hinterp_deriv,Nsk,&
                      l1,m1,l2,m2)
              ENDIF
           ELSE
              ! swap atoms if typj > typi
              s_deriv(:) = slako_transformation_gradient(&
                   rij,-xij,-yij,-zij,&
                   Sinterp,Sinterp_deriv,Nsk,&
                   l2,m2,l1,m1)
              h0_deriv(:) = slako_transformation_gradient(&
                   rij,-xij,-yij,-zij,&
                   Hinterp,Hinterp_deriv,Nsk,&
                   l2,m2,l1,m1)
           ENDIF
           ! contributions are SUBTRACTED because we loop over j, while
           ! in the python code we loop over i and ADD them
           gradS(c3b:c3e,mu,nu) = gradS(c3b:c3e,mu,nu) - s_deriv(:)
           gradH0(c3b:c3e,mu,nu) = gradH0(c3b:c3e,mu,nu) - h0_deriv(:)
           ! S and H0 are hermitian/symmetric
           gradS(c3b:c3e,nu,mu) = gradS(c3b:c3e,nu,mu) - s_deriv(:)
           gradH0(c3b:c3e,nu,mu) = gradH0(c3b:c3e,nu,mu) - h0_deriv(:)           
        ENDIF
     ENDDO
  ENDDO
END SUBROUTINE gradients_H0andS

!*****************************************************************
! computes the overlap matrix between two sets of atoms using
! Slater-Koster rules
!*****************************************************************

SUBROUTINE overlap12(&
     Norb1,Norb2,Nat1,Nat2,&
     ind1,typ1,ls1,ms1,&
     ind2,typ2,ls2,ms2,&
     r,x,y,z,&
     deg,tab_filled_SH0,&
     S_knots,S_coefs,&
     Npair,Nsk,Npts,&
     S)
  ! release global interpreter lock
  !f2py threadsafe
  IMPLICIT NONE

  ! INPUT
  ! number of atoms in geometries 1 and 2
  INTEGER, INTENT(IN) :: Nat1, Nat2
  ! number of valence orbitals in geometries 1 and 2
  INTEGER, INTENT(IN) :: Norb1, Norb2
  ! indeces of the atoms belonging to the orbitals
  INTEGER, INTENT(IN) :: ind1(Norb1),ind2(Norb2)
  ! atom types of the atoms belonging to the orbitals
  INTEGER, INTENT(IN) :: typ1(Norb1),typ2(Norb2)
  ! quantum numbers l,m of the orbitals
  INTEGER, INTENT(IN) :: ls1(Norb1),ms1(Norb1), ls2(Norb2),ms2(Norb2)
  ! distances and directional cosines between atoms in geometries 1 and 2
  DOUBLE PRECISION, INTENT(IN) :: r(Nat1,Nat2), &
       x(Nat1,Nat2), y(Nat1,Nat2), z(Nat1,Nat2)
  ! degree of the B-splines
  INTEGER, INTENT(IN) :: deg
  ! If the SK table t of a atom pair with index ipair contains data, the flag
  ! tab_filled_SH0[ipair,t] is set to 1
  INTEGER, INTENT(IN) :: tab_filled_SH0(Npair,Nsk)
  ! spline data with knots and coefficients for S and H0
  INTEGER, INTENT(IN) :: Npair ! number of unique atom pairs
  INTEGER, INTENT(IN) :: Nsk   ! maximum number of SK tables (for some pairs, like h-h, most tables are empty)
  INTEGER, INTENT(IN) :: Npts  ! number of knots in a B-spline
  !
  DOUBLE PRECISION, INTENT(IN) :: &
       S_knots(Npair,Nsk,Npts), S_coefs(Npair,Nsk,Npts)

  ! OUTPUT
  ! matrix elements of overlap S 
  DOUBLE PRECISION, INTENT(OUT) :: S(Norb1,Norb2)
       
  ! LOCAL VARIABLES
  INTEGER :: mu,nu  ! enumerate orbitals
  INTEGER :: l1,m1,l2,m2 ! quantum numbers of orbitals mu and nu
  INTEGER :: i,j,ip1,jp1    ! enumerate atom indeces
  INTEGER :: typi,typj ! atom types
  INTEGER :: typMin,typMax ! ordered atom types
  INTEGER :: tab    ! enumerate SK tables
  INTEGER :: ipair  ! encodes the atom pair 
  DOUBLE PRECISION :: Sinterp(Nsk), Hinterp(Nsk)  ! values of interpolated SK tables
  DOUBLE PRECISION :: rij,xij,yij,zij  ! r(i,j), x(i,j), etc.
  INTEGER :: ier ! return code of FITPACK routines

  !$omp parallel do private(nu,j,typj,l2,m2,mu,i,typi,l1,m1), &
  !$omp& private(ip1,jp1,rij,xij,yij,zij,typMin,typMax,ipair), &
  !$omp& private(tab,Sinterp,ier) schedule(dynamic,1)
  DO nu=1,Norb2
     ! Which atom does orbital mu sit on?
     j=ind2(nu)
     ! What's the type of this atom?
     typj=typ2(nu)
     ! What are the quantum numbers of orbital mu
     l2=ls2(nu)
     m2=ms2(nu)
     DO mu=1,Norb1
        ! Which atom does orbital mu sit on?
        i=ind1(mu)
        ! What's the type of this atom?
        typi=typ1(mu)
        ! What are the quantum numbers of orbital mu
        l1=ls1(mu)
        m1=ms1(mu)
        ! If the atoms are too far apart, skip the calculation
        ip1=i+1 ! Fortran indeces start at 1
        jp1=j+1

        ! directional cosines and distance between atom i in geometry 1 and atom j in geometry 2
        rij=r(jp1,ip1)
        xij=x(jp1,ip1)
        yij=y(jp1,ip1)
        zij=z(jp1,ip1)

        ! interpolate Slater-Koster tables
        ! find the chunk of the tables belonging to atom pair (i,j)
        ! When ipair is calculated from (typi,typj) the tuple has to be ordered
        typMin=min(typi,typj)
        typMax=max(typi,typj)
        
        ipair=(typMax*(typMax+1))/2+typMin + 1 ! Fortran indeces start at 1
        DO tab=1,Nsk
           IF (tab_filled_SH0(ipair,tab) == 1) THEN
              ! this table has useful data, interpolate it
              ! assuming the knots of the B-spline are uniformly spaced
              call splev_uniform(S_knots(ipair,tab,:),Npts,S_coefs(ipair,tab,:),deg,&
                   rij, Sinterp(tab), ier)
           ENDIF
        ENDDO

        IF (typi <= typj) THEN
           S(mu,nu) = slako_transformation(rij,xij,yij,zij,Sinterp,Nsk,&
                l1,m1,l2,m2)
        ELSE
              ! swap atoms if typj > typi
           S(mu,nu) = slako_transformation(rij,-xij,-yij,-zij,Sinterp,Nsk,&
                l2,m2,l1,m1)
        ENDIF
     ENDDO
  ENDDO
END SUBROUTINE overlap12


END MODULE slako

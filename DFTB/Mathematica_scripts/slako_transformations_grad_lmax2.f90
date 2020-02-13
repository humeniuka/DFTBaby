! This file has been generated automatically


! transformation rules for GRADIENTS of matrix elements
FUNCTION slako_transformation_gradient(r,x,y,z, SorH, SorHd, N, l1,m1,l2,m2)
   IMPLICIT NONE
   DOUBLE PRECISION slako_transformation_gradient(3) ! function returns 3d vector
   ! x,y,z are directional cosines, r is the distance between the two centers
   DOUBLE PRECISION, INTENT(IN) :: r,x,y,z
   INTEGER, INTENT(IN) :: N  ! length of arrays SorH and SorHg
   ! values of the N Slater-Koster tables for S or H0 evaluated at distance r
   DOUBLE PRECISION, INTENT(IN) :: SorH(N)
   ! values of the derivatives of S or H0 evaluated at distance r
   DOUBLE PRECISION, INTENT(IN) :: SorHd(N)
   ! orbital qm numbers for center 1 and center 2
   INTEGER, INTENT(IN) :: l1,m1,l2,m2
   ! Local Variables
   INTEGER, PARAMETER :: lmax=2
   DOUBLE PRECISION, PARAMETER :: sqrt3=1.7320508075688772
   ! Result grad[S(x,y,z)] or grad[H(x,y,z)] after applying SK rules
   DOUBLE PRECISION :: grad(3)
   INTEGER :: i  ! index that encodes the tuple (l1,m1,l2,m2)

   ! First we need to transform the tuple (l1,m1,l2,m2) into a unique integer
   ! so that the compiler can build a branching table for each case.
   ! Valid ranges for qm numbers: 0 <= l1,l2 <= lmax, -lmax <= m1,m2 <= lmax
   i = l1*1000+(lmax+m1)*100+l2*10+(lmax+m2)
!$OMP PARALLEL NUM_THREADS(3)
!$OMP SECTIONS
!$OMP SECTION
   SELECT CASE(i)
	CASE (202)
	! (0,0,0,0)
	 grad(1) = x*SorHd(1)
	CASE (211)
	! (0,0,1,-1)
	 grad(1) = x*y*(-(SorH(3)/r) + SorHd(3))
	CASE (212)
	! (0,0,1,0)
	 grad(1) = x*z*(-(SorH(3)/r) + SorHd(3))
	CASE (213)
	! (0,0,1,1)
	 grad(1) = -(((-1 + x**2)*SorH(3))/r) + x**2*SorHd(3)
	CASE (220)
	! (0,0,2,-2)
	 grad(1) = (sqrt3*y*((1 - 2*x**2)*SorH(4) + r*x**2*SorHd(4)))/r
	CASE (221)
	! (0,0,2,-1)
	 grad(1) = (sqrt3*x*y*z*(-2*SorH(4) + r*SorHd(4)))/r
	CASE (222)
	! (0,0,2,0)
	 grad(1) = (x*(-1 + x**2 + y**2 - 2*z**2)*SorH(4))/r - (x*(x**2 + y**2 - &
2*z**2)*SorHd(4))/2.
	CASE (223)
	! (0,0,2,1)
	 grad(1) = (sqrt3*z*((1 - 2*x**2)*SorH(4) + r*x**2*SorHd(4)))/r
	CASE (224)
	! (0,0,2,2)
	 grad(1) = (sqrt3*x*(2*(1 - x**2 + y**2)*SorH(4) + r*(x - y)*(x + y)*SorHd(4)))/(2.*r)
	CASE (1102)
	! (1,-1,0,0)
	 grad(1) = x*y*(-(SorH(5)/r) + SorHd(5))
	CASE (1111)
	! (1,-1,1,-1)
	 grad(1) = (x*(-2*(-1 + x**2 + z**2)*SorH(6) + r*(x**2 + z**2)*SorHd(6) + &
y**2*(-2*SorH(7) + r*SorHd(7))))/r
	CASE (1112)
	! (1,-1,1,0)
	 grad(1) = (x*y*z*(2*SorH(6) - 2*SorH(7) + r*(-SorHd(6) + SorHd(7))))/r
	CASE (1113)
	! (1,-1,1,1)
	 grad(1) = (y*((-1 + 2*x**2)*SorH(6) + SorH(7) + x**2*(-2*SorH(7) + r*(-SorHd(6) + &
SorHd(7)))))/r
	CASE (1120)
	! (1,-1,2,-2)
	 grad(1) = ((-y**2 + z**2 - 3*x**2*(-1 + x**2 - y**2 + z**2))*SorH(8) + r*x**2*(x**2 - &
y**2 + z**2)*SorHd(8) + sqrt3*y**2*((1 - 3*x**2)*SorH(9) + r*x**2*SorHd(9)))/r
	CASE (1121)
	! (1,-1,2,-1)
	 grad(1) = (x*z*((2 - 3*x**2 + 3*y**2 - 3*z**2)*SorH(8) + r*(x**2 - y**2 + z**2)*SorHd(8) &
+ sqrt3*y**2*(-3*SorH(9) + r*SorHd(9))))/r
	CASE (1122)
	! (1,-1,2,0)
	 grad(1) = -(x*y*((2 - 3*x**2 - 3*y**2)*SorH(9) + 2*z**2*(-3*sqrt3*SorH(8) + 3*SorH(9) + &
r*sqrt3*SorHd(8)) + r*(x**2 + y**2 - 2*z**2)*SorHd(9)))/(2.*r)
	CASE (1123)
	! (1,-1,2,1)
	 grad(1) = (y*z*((-2 + 6*x**2)*SorH(8) + sqrt3*(1 - 3*x**2)*SorH(9) + r*x**2*(-2*SorHd(8) &
+ sqrt3*SorHd(9))))/r
	CASE (1124)
	! (1,-1,2,2)
	 grad(1) = (x*y*(2*(-4 + 6*x**2 + 3*z**2)*SorH(8) + sqrt3*(2 - 3*x**2 + 3*y**2)*SorH(9) - &
2*r*(2*x**2 + z**2)*SorHd(8) + r*sqrt3*(x - y)*(x + y)*SorHd(9)))/(2.*r)
	CASE (1202)
	! (1,0,0,0)
	 grad(1) = x*z*(-(SorH(5)/r) + SorHd(5))
	CASE (1211)
	! (1,0,1,-1)
	 grad(1) = (x*y*z*(2*SorH(6) - 2*SorH(7) + r*(-SorHd(6) + SorHd(7))))/r
	CASE (1212)
	! (1,0,1,0)
	 grad(1) = (x*(-2*(-1 + x**2 + y**2)*SorH(6) + r*(x**2 + y**2)*SorHd(6) + &
z**2*(-2*SorH(7) + r*SorHd(7))))/r
	CASE (1213)
	! (1,0,1,1)
	 grad(1) = (z*((-1 + 2*x**2)*SorH(6) + SorH(7) + x**2*(-2*SorH(7) + r*(-SorHd(6) + &
SorHd(7)))))/r
	CASE (1220)
	! (1,0,2,-2)
	 grad(1) = (y*z*((-2 + 6*x**2)*SorH(8) + sqrt3*(1 - 3*x**2)*SorH(9) + r*x**2*(-2*SorHd(8) &
+ sqrt3*SorHd(9))))/r
	CASE (1221)
	! (1,0,2,-1)
	 grad(1) = (x*y*((2 - 3*x**2 - 3*y**2 + 3*z**2)*SorH(8) + r*(x**2 + y**2 - z**2)*SorHd(8) &
+ sqrt3*z**2*(-3*SorH(9) + r*SorHd(9))))/r
	CASE (1222)
	! (1,0,2,0)
	 grad(1) = -(x*z*(2*sqrt3*(-2 + 3*x**2 + 3*y**2)*SorH(8) + (2 - 3*x**2 - 3*y**2 + &
6*z**2)*SorH(9) + r*(-2*sqrt3*(x**2 + y**2)*SorHd(8) + (x**2 + y**2 - &
2*z**2)*SorHd(9))))/(2.*r)
	CASE (1223)
	! (1,0,2,1)
	 grad(1) = (((y - z)*(y + z) - 3*x**2*(-1 + x**2 + y**2 - z**2))*SorH(8) + r*x**2*(x**2 + &
y**2 - z**2)*SorHd(8) + sqrt3*z**2*((1 - 3*x**2)*SorH(9) + r*x**2*SorHd(9)))/r
	CASE (1224)
	! (1,0,2,2)
	 grad(1) = (x*z*((-4 + 6*x**2 - 6*y**2)*SorH(8) + 2*sqrt3*SorH(9) - (x - y)*(x + &
y)*(3*sqrt3*SorH(9) + 2*r*SorHd(8) - r*sqrt3*SorHd(9))))/(2.*r)
	CASE (1302)
	! (1,1,0,0)
	 grad(1) = -(((-1 + x**2)*SorH(5))/r) + x**2*SorHd(5)
	CASE (1311)
	! (1,1,1,-1)
	 grad(1) = (y*((-1 + 2*x**2)*SorH(6) + SorH(7) + x**2*(-2*SorH(7) + r*(-SorHd(6) + &
SorHd(7)))))/r
	CASE (1312)
	! (1,1,1,0)
	 grad(1) = (z*((-1 + 2*x**2)*SorH(6) + SorH(7) + x**2*(-2*SorH(7) + r*(-SorHd(6) + &
SorHd(7)))))/r
	CASE (1313)
	! (1,1,1,1)
	 grad(1) = (x*(-2*(-1 + x**2)*SorH(7) - (y**2 + z**2)*(2*SorH(6) - r*SorHd(6)) + &
r*x**2*SorHd(7)))/r
	CASE (1320)
	! (1,1,2,-2)
	 grad(1) = (x*y*((-2 + 3*x**2 - 3*y**2 - 3*z**2)*SorH(8) + sqrt3*(2 - 3*x**2)*SorH(9) + &
r*(-x**2 + y**2 + z**2)*SorHd(8) + r*sqrt3*x**2*SorHd(9)))/r
	CASE (1321)
	! (1,1,2,-1)
	 grad(1) = (y*z*((-2 + 6*x**2)*SorH(8) + sqrt3*(1 - 3*x**2)*SorH(9) + r*x**2*(-2*SorHd(8) &
+ sqrt3*SorHd(9))))/r
	CASE (1322)
	! (1,1,2,0)
	 grad(1) = -(2*sqrt3*(1 - 3*x**2)*z**2*SorH(8) + (y**2 - 2*z**2 - 3*x**2*(-1 + x**2 + &
y**2 - 2*z**2))*SorH(9) + r*x**2*(2*sqrt3*z**2*SorHd(8) + (x**2 + y**2 - &
2*z**2)*SorHd(9)))/(2.*r)
	CASE (1323)
	! (1,1,2,1)
	 grad(1) = (x*z*((-2 + 3*x**2 - 3*y**2 - 3*z**2)*SorH(8) + sqrt3*(2 - 3*x**2)*SorH(9) + &
r*(-x**2 + y**2 + z**2)*SorHd(8) + r*sqrt3*x**2*SorHd(9)))/r
	CASE (1324)
	! (1,1,2,2)
	 grad(1) = (-2*(-1 + 3*x**2)*(2*y**2 + z**2)*SorH(8) - sqrt3*(3*x**4 + y**2 - 3*x**2*(1 + &
y**2))*SorH(9) + r*x**2*(2*(2*y**2 + z**2)*SorHd(8) + sqrt3*(x - y)*(x + &
y)*SorHd(9)))/(2.*r)
	CASE (2002)
	! (2,-2,0,0)
	 grad(1) = (sqrt3*y*((1 - 2*x**2)*SorH(10) + r*x**2*SorHd(10)))/r
	CASE (2011)
	! (2,-2,1,-1)
	 grad(1) = ((-y**2 + z**2 - 3*x**2*(-1 + x**2 - y**2 + z**2))*SorH(11) + r*x**2*(x**2 - &
y**2 + z**2)*SorHd(11) + sqrt3*y**2*((1 - 3*x**2)*SorH(12) + &
r*x**2*SorHd(12)))/r
	CASE (2012)
	! (2,-2,1,0)
	 grad(1) = (y*z*((-2 + 6*x**2)*SorH(11) + sqrt3*(1 - 3*x**2)*SorH(12) + &
r*x**2*(-2*SorHd(11) + sqrt3*SorHd(12))))/r
	CASE (2013)
	! (2,-2,1,1)
	 grad(1) = (x*y*((-2 + 3*x**2 - 3*y**2 - 3*z**2)*SorH(11) + sqrt3*(2 - 3*x**2)*SorH(12) + &
r*(-x**2 + y**2 + z**2)*SorHd(11) + r*sqrt3*x**2*SorHd(12)))/r
	CASE (2020)
	! (2,-2,2,-2)
	 grad(1) = (x*(-2*(y**2 + z**2)*(-1 + 2*x**2 + 2*z**2)*SorH(13) + 2*(-2*x**4 + z**2 + &
x**2*(2 + 4*y**2 - 2*z**2) - 2*y**2*(1 + y**2 + z**2))*SorH(14) + &
6*y**2*SorH(15) - 12*x**2*y**2*SorH(15) + r*x**2*y**2*SorHd(13) + &
r*x**2*z**2*SorHd(13) + r*y**2*z**2*SorHd(13) + r*z**4*SorHd(13) + &
r*x**4*SorHd(14) - 2*r*x**2*y**2*SorHd(14) + r*y**4*SorHd(14) + &
r*x**2*z**2*SorHd(14) + r*y**2*z**2*SorHd(14) + 3*r*x**2*y**2*SorHd(15)))/r
	CASE (2021)
	! (2,-2,2,-1)
	 grad(1) = (z*((4*x**4 - z**2 + x**2*(-3 + 4*z**2))*SorH(13) + (-3*y**2 + z**2 + x**2*(3 &
- 4*x**2 + 12*y**2 - 4*z**2))*SorH(14) + 3*y**2*SorH(15) - &
12*x**2*y**2*SorH(15) - r*x**4*SorHd(13) - r*x**2*z**2*SorHd(13) + &
r*x**4*SorHd(14) - 3*r*x**2*y**2*SorHd(14) + r*x**2*z**2*SorHd(14) + &
3*r*x**2*y**2*SorHd(15)))/r
	CASE (2022)
	! (2,-2,2,0)
	 grad(1) = (sqrt3*y*((-4*x**4 + y**2 + 2*z**2 + x**2*(3 - 4*y**2 - 8*z**2))*SorH(13) + &
4*(-1 + 4*x**2)*z**2*SorH(14) - 3*x**2*SorH(15) + 4*x**4*SorH(15) - &
y**2*SorH(15) + 4*x**2*y**2*SorH(15) + 2*z**2*SorH(15) - 8*x**2*z**2*SorH(15) &
+ r*x**4*SorHd(13) + r*x**2*y**2*SorHd(13) + 2*r*x**2*z**2*SorHd(13) - &
4*r*x**2*z**2*SorHd(14) - r*x**2*(x**2 + y**2 - 2*z**2)*SorHd(15)))/(2.*r)
	CASE (2023)
	! (2,-2,2,1)
	 grad(1) = (x*y*z*(4*(y**2 + z**2)*SorH(13) - 4*(y**2 + z**2)*SorH(14) + 6*(-1 + &
2*x**2)*(SorH(14) - SorH(15)) + r*(-((y**2 + z**2)*SorHd(13)) + (-3*x**2 + &
y**2 + z**2)*SorHd(14) + 3*x**2*SorHd(15))))/r
	CASE (2024)
	! (2,-2,2,2)
	 grad(1) = (y*(-((4*x**4 + y**2 - x**2*(3 + 4*y**2))*SorH(13)) + 4*(y**2 + x**2*(-3 + &
4*x**2 - 4*y**2))*SorH(14) + 9*x**2*SorH(15) - 12*x**4*SorH(15) - &
3*y**2*SorH(15) + 12*x**2*y**2*SorH(15) + r*x**4*SorHd(13) - &
r*x**2*y**2*SorHd(13) - 4*r*x**4*SorHd(14) + 4*r*x**2*y**2*SorHd(14) + &
3*r*x**2*(x - y)*(x + y)*SorHd(15)))/(2.*r)
	CASE (2102)
	! (2,-1,0,0)
	 grad(1) = (sqrt3*x*y*z*(-2*SorH(10) + r*SorHd(10)))/r
	CASE (2111)
	! (2,-1,1,-1)
	 grad(1) = (x*z*((2 - 3*x**2 + 3*y**2 - 3*z**2)*SorH(11) + r*(x**2 - y**2 + &
z**2)*SorHd(11) + sqrt3*y**2*(-3*SorH(12) + r*SorHd(12))))/r
	CASE (2112)
	! (2,-1,1,0)
	 grad(1) = (x*y*((2 - 3*x**2 - 3*y**2 + 3*z**2)*SorH(11) + r*(x**2 + y**2 - &
z**2)*SorHd(11) + sqrt3*z**2*(-3*SorH(12) + r*SorHd(12))))/r
	CASE (2113)
	! (2,-1,1,1)
	 grad(1) = (y*z*((-2 + 6*x**2)*SorH(11) + sqrt3*(1 - 3*x**2)*SorH(12) + &
r*x**2*(-2*SorHd(11) + sqrt3*SorHd(12))))/r
	CASE (2120)
	! (2,-1,2,-2)
	 grad(1) = (z*((4*x**4 - z**2 + x**2*(-3 + 4*z**2))*SorH(13) + (-3*y**2 + z**2 + x**2*(3 &
- 4*x**2 + 12*y**2 - 4*z**2))*SorH(14) + 3*y**2*SorH(15) - &
12*x**2*y**2*SorH(15) - r*x**4*SorHd(13) - r*x**2*z**2*SorHd(13) + &
r*x**4*SorHd(14) - 3*r*x**2*y**2*SorHd(14) + r*x**2*z**2*SorHd(14) + &
3*r*x**2*y**2*SorHd(15)))/r
	CASE (2121)
	! (2,-1,2,-1)
	 grad(1) = (x*(2*(y**2 - 2*x**2*(-1 + x**2 + y**2) + z**2 - 2*(x**2 + &
y**2)*z**2)*SorH(13) + 2*(-2*y**4 + z**2 - 2*z**2*(x**2 + z**2) + y**2*(1 - &
2*x**2 + 4*z**2))*SorH(14) - 12*y**2*z**2*SorH(15) + r*x**4*SorHd(13) + &
r*x**2*y**2*SorHd(13) + r*x**2*z**2*SorHd(13) + r*y**2*z**2*SorHd(13) + &
r*x**2*y**2*SorHd(14) + r*y**4*SorHd(14) + r*x**2*z**2*SorHd(14) - &
2*r*y**2*z**2*SorHd(14) + r*z**4*SorHd(14) + 3*r*y**2*z**2*SorHd(15)))/r
	CASE (2122)
	! (2,-1,2,0)
	 grad(1) = (sqrt3*x*y*z*((-2 + 4*x**2 + 4*y**2)*SorH(13) + (4 - 8*x**2 - 8*y**2 + &
8*z**2)*SorH(14) - 2*SorH(15) + 4*x**2*SorH(15) + 4*y**2*SorH(15) - &
8*z**2*SorH(15) - r*x**2*SorHd(13) - r*y**2*SorHd(13) + 2*r*x**2*SorHd(14) + &
2*r*y**2*SorHd(14) - 2*r*z**2*SorHd(14) - r*(x**2 + y**2 - &
2*z**2)*SorHd(15)))/(2.*r)
	CASE (2123)
	! (2,-1,2,1)
	 grad(1) = (y*((4*x**4 - y**2 + x**2*(-3 + 4*y**2))*SorH(13) + (y**2 - 3*z**2 + x**2*(3 - &
4*x**2 - 4*y**2 + 12*z**2))*SorH(14) + 3*z**2*SorH(15) - 12*x**2*z**2*SorH(15) &
- r*x**4*SorHd(13) - r*x**2*y**2*SorHd(13) + r*x**4*SorHd(14) + &
r*x**2*y**2*SorHd(14) - 3*r*x**2*z**2*SorHd(14) + 3*r*x**2*z**2*SorHd(15)))/r
	CASE (2124)
	! (2,-1,2,2)
	 grad(1) = (x*y*z*(-2*(-3 + 6*x**2 + 2*y**2 + 4*z**2)*SorH(13) + 4*(-3 + 6*x**2 - 2*y**2 &
+ 2*z**2)*SorH(14) + 6*SorH(15) - 12*x**2*SorH(15) + 12*y**2*SorH(15) + &
3*r*x**2*SorHd(13) + r*y**2*SorHd(13) + 2*r*z**2*SorHd(13) - &
6*r*x**2*SorHd(14) + 2*r*y**2*SorHd(14) - 2*r*z**2*SorHd(14) + 3*r*(x - y)*(x &
+ y)*SorHd(15)))/(2.*r)
	CASE (2202)
	! (2,0,0,0)
	 grad(1) = (x*(-1 + x**2 + y**2 - 2*z**2)*SorH(10))/r - (x*(x**2 + y**2 - &
2*z**2)*SorHd(10))/2.
	CASE (2211)
	! (2,0,1,-1)
	 grad(1) = -(x*y*((2 - 3*x**2 - 3*y**2)*SorH(12) + 2*z**2*(-3*sqrt3*SorH(11) + 3*SorH(12) &
+ r*sqrt3*SorHd(11)) + r*(x**2 + y**2 - 2*z**2)*SorHd(12)))/(2.*r)
	CASE (2212)
	! (2,0,1,0)
	 grad(1) = -(x*z*(2*sqrt3*(-2 + 3*x**2 + 3*y**2)*SorH(11) + (2 - 3*x**2 - 3*y**2 + &
6*z**2)*SorH(12) + r*(-2*sqrt3*(x**2 + y**2)*SorHd(11) + (x**2 + y**2 - &
2*z**2)*SorHd(12))))/(2.*r)
	CASE (2213)
	! (2,0,1,1)
	 grad(1) = -(2*sqrt3*(1 - 3*x**2)*z**2*SorH(11) + (y**2 - 2*z**2 - 3*x**2*(-1 + x**2 + &
y**2 - 2*z**2))*SorH(12) + r*x**2*(2*sqrt3*z**2*SorHd(11) + (x**2 + y**2 - &
2*z**2)*SorHd(12)))/(2.*r)
	CASE (2220)
	! (2,0,2,-2)
	 grad(1) = (sqrt3*y*((-4*x**4 + y**2 + 2*z**2 + x**2*(3 - 4*y**2 - 8*z**2))*SorH(13) + &
4*(-1 + 4*x**2)*z**2*SorH(14) - 3*x**2*SorH(15) + 4*x**4*SorH(15) - &
y**2*SorH(15) + 4*x**2*y**2*SorH(15) + 2*z**2*SorH(15) - 8*x**2*z**2*SorH(15) &
+ r*x**4*SorHd(13) + r*x**2*y**2*SorHd(13) + 2*r*x**2*z**2*SorHd(13) - &
4*r*x**2*z**2*SorHd(14) - r*x**2*(x**2 + y**2 - 2*z**2)*SorHd(15)))/(2.*r)
	CASE (2221)
	! (2,0,2,-1)
	 grad(1) = (sqrt3*x*y*z*((-2 + 4*x**2 + 4*y**2)*SorH(13) + (4 - 8*x**2 - 8*y**2 + &
8*z**2)*SorH(14) - 2*SorH(15) + 4*x**2*SorH(15) + 4*y**2*SorH(15) - &
8*z**2*SorH(15) - r*x**2*SorHd(13) - r*y**2*SorHd(13) + 2*r*x**2*SorHd(14) + &
2*r*y**2*SorHd(14) - 2*r*z**2*SorHd(14) - r*(x**2 + y**2 - &
2*z**2)*SorHd(15)))/(2.*r)
	CASE (2222)
	! (2,0,2,0)
	 grad(1) = (x*(-12*(-1 + x**2 + y**2)*(x**2 + y**2)*SorH(13) - 24*(-1 + 2*x**2 + &
2*y**2)*z**2*SorH(14) + 4*x**2*SorH(15) - 4*x**4*SorH(15) + 4*y**2*SorH(15) - &
8*x**2*y**2*SorH(15) - 4*y**4*SorH(15) - 8*z**2*SorH(15) + &
16*x**2*z**2*SorH(15) + 16*y**2*z**2*SorH(15) - 16*z**4*SorH(15) + &
3*r*x**4*SorHd(13) + 6*r*x**2*y**2*SorHd(13) + 3*r*y**4*SorHd(13) + &
12*r*x**2*z**2*SorHd(14) + 12*r*y**2*z**2*SorHd(14) + r*(x**2 + y**2 - &
2*z**2)**2*SorHd(15)))/(4.*r)
	CASE (2223)
	! (2,0,2,1)
	 grad(1) = (sqrt3*z*((4*x**4 - y**2 + x**2*(-3 + 4*y**2))*SorH(13) + 2*((y - z)*(y + z) + &
x**2*(3 - 4*x**2 - 4*y**2 + 4*z**2))*SorH(14) - 3*x**2*SorH(15) + &
4*x**4*SorH(15) - y**2*SorH(15) + 4*x**2*y**2*SorH(15) + 2*z**2*SorH(15) - &
8*x**2*z**2*SorH(15) - r*x**4*SorHd(13) - r*x**2*y**2*SorHd(13) + &
2*r*x**4*SorHd(14) + 2*r*x**2*y**2*SorHd(14) - 2*r*x**2*z**2*SorHd(14) - &
r*x**2*(x**2 + y**2 - 2*z**2)*SorHd(15)))/(2.*r)
	CASE (2224)
	! (2,0,2,2)
	 grad(1) = (sqrt3*x*(4*(-x**4 + y**4 + z**2 + 2*y**2*z**2 + x**2*(1 - 2*z**2))*SorH(13) + &
8*(-1 + 2*x**2 - 2*y**2)*z**2*SorH(14) - 4*x**2*SorH(15) + 4*x**4*SorH(15) - &
4*y**4*SorH(15) + 4*z**2*SorH(15) - 8*x**2*z**2*SorH(15) + &
8*y**2*z**2*SorH(15) + r*x**4*SorHd(13) - r*y**4*SorHd(13) + &
2*r*x**2*z**2*SorHd(13) - 2*r*y**2*z**2*SorHd(13) - 4*r*x**2*z**2*SorHd(14) + &
4*r*y**2*z**2*SorHd(14) - r*(x - y)*(x + y)*(x**2 + y**2 - &
2*z**2)*SorHd(15)))/(4.*r)
	CASE (2302)
	! (2,1,0,0)
	 grad(1) = (sqrt3*z*((1 - 2*x**2)*SorH(10) + r*x**2*SorHd(10)))/r
	CASE (2311)
	! (2,1,1,-1)
	 grad(1) = (y*z*((-2 + 6*x**2)*SorH(11) + sqrt3*(1 - 3*x**2)*SorH(12) + &
r*x**2*(-2*SorHd(11) + sqrt3*SorHd(12))))/r
	CASE (2312)
	! (2,1,1,0)
	 grad(1) = (((y - z)*(y + z) - 3*x**2*(-1 + x**2 + y**2 - z**2))*SorH(11) + r*x**2*(x**2 &
+ y**2 - z**2)*SorHd(11) + sqrt3*z**2*((1 - 3*x**2)*SorH(12) + &
r*x**2*SorHd(12)))/r
	CASE (2313)
	! (2,1,1,1)
	 grad(1) = (x*z*((-2 + 3*x**2 - 3*y**2 - 3*z**2)*SorH(11) + sqrt3*(2 - 3*x**2)*SorH(12) + &
r*(-x**2 + y**2 + z**2)*SorHd(11) + r*sqrt3*x**2*SorHd(12)))/r
	CASE (2320)
	! (2,1,2,-2)
	 grad(1) = (x*y*z*(4*(y**2 + z**2)*SorH(13) - 4*(y**2 + z**2)*SorH(14) + 6*(-1 + &
2*x**2)*(SorH(14) - SorH(15)) + r*(-((y**2 + z**2)*SorHd(13)) + (-3*x**2 + &
y**2 + z**2)*SorHd(14) + 3*x**2*SorHd(15))))/r
	CASE (2321)
	! (2,1,2,-1)
	 grad(1) = (y*((4*x**4 - y**2 + x**2*(-3 + 4*y**2))*SorH(13) + (y**2 - 3*z**2 + x**2*(3 - &
4*x**2 - 4*y**2 + 12*z**2))*SorH(14) + 3*z**2*SorH(15) - 12*x**2*z**2*SorH(15) &
- r*x**4*SorHd(13) - r*x**2*y**2*SorHd(13) + r*x**4*SorHd(14) + &
r*x**2*y**2*SorHd(14) - 3*r*x**2*z**2*SorHd(14) + 3*r*x**2*z**2*SorHd(15)))/r
	CASE (2322)
	! (2,1,2,0)
	 grad(1) = (sqrt3*z*((4*x**4 - y**2 + x**2*(-3 + 4*y**2))*SorH(13) + 2*((y - z)*(y + z) + &
x**2*(3 - 4*x**2 - 4*y**2 + 4*z**2))*SorH(14) - 3*x**2*SorH(15) + &
4*x**4*SorH(15) - y**2*SorH(15) + 4*x**2*y**2*SorH(15) + 2*z**2*SorH(15) - &
8*x**2*z**2*SorH(15) - r*x**4*SorHd(13) - r*x**2*y**2*SorHd(13) + &
2*r*x**4*SorHd(14) + 2*r*x**2*y**2*SorHd(14) - 2*r*x**2*z**2*SorHd(14) - &
r*x**2*(x**2 + y**2 - 2*z**2)*SorHd(15)))/(2.*r)
	CASE (2323)
	! (2,1,2,1)
	 grad(1) = (x*(-2*(-1 + 2*x**2 + 2*y**2)*(y**2 + z**2)*SorH(13) + 2*(y**2 - 2*z**2 - &
2*(x**2*(-1 + x**2 + y**2) + (-2*x**2 + y**2)*z**2 + z**4))*SorH(14) + &
6*z**2*SorH(15) - 12*x**2*z**2*SorH(15) + r*x**2*y**2*SorHd(13) + &
r*y**4*SorHd(13) + r*x**2*z**2*SorHd(13) + r*y**2*z**2*SorHd(13) + &
r*x**4*SorHd(14) + r*x**2*y**2*SorHd(14) - 2*r*x**2*z**2*SorHd(14) + &
r*y**2*z**2*SorHd(14) + r*z**4*SorHd(14) + 3*r*x**2*z**2*SorHd(15)))/r
	CASE (2324)
	! (2,1,2,2)
	 grad(1) = -(z*((-4*x**4 + 3*y**2 + 2*z**2 + x**2*(3 - 12*y**2 - 8*z**2))*SorH(13) - &
2*(4*x**4 + 3*y**2 + z**2 - x**2*(3 + 12*y**2 + 4*z**2))*SorH(14) - &
9*x**2*SorH(15) + 12*x**4*SorH(15) + 3*y**2*SorH(15) - 12*x**2*y**2*SorH(15) + &
r*x**4*SorHd(13) + 3*r*x**2*y**2*SorHd(13) + 2*r*x**2*z**2*SorHd(13) + &
2*r*x**4*SorHd(14) - 6*r*x**2*y**2*SorHd(14) - 2*r*x**2*z**2*SorHd(14) - &
3*r*x**2*(x - y)*(x + y)*SorHd(15)))/(2.*r)
	CASE (2402)
	! (2,2,0,0)
	 grad(1) = (sqrt3*x*(2*(1 - x**2 + y**2)*SorH(10) + r*(x - y)*(x + y)*SorHd(10)))/(2.*r)
	CASE (2411)
	! (2,2,1,-1)
	 grad(1) = (x*y*(2*(-4 + 6*x**2 + 3*z**2)*SorH(11) + sqrt3*(2 - 3*x**2 + 3*y**2)*SorH(12) &
- 2*r*(2*x**2 + z**2)*SorHd(11) + r*sqrt3*(x - y)*(x + y)*SorHd(12)))/(2.*r)
	CASE (2412)
	! (2,2,1,0)
	 grad(1) = (x*z*((-4 + 6*x**2 - 6*y**2)*SorH(11) + sqrt3*(2 - 3*x**2 + 3*y**2)*SorH(12) - &
r*(x - y)*(x + y)*(2*SorHd(11) - sqrt3*SorHd(12))))/(2.*r)
	CASE (2413)
	! (2,2,1,1)
	 grad(1) = (-2*(-1 + 3*x**2)*(2*y**2 + z**2)*SorH(11) - sqrt3*(3*x**4 + y**2 - 3*x**2*(1 &
+ y**2))*SorH(12) + r*x**2*(2*(2*y**2 + z**2)*SorHd(11) + sqrt3*(x - y)*(x + &
y)*SorHd(12)))/(2.*r)
	CASE (2420)
	! (2,2,2,-2)
	 grad(1) = (y*(-((4*x**4 + y**2 - x**2*(3 + 4*y**2))*SorH(13)) + 4*(y**2 + x**2*(-3 + &
4*x**2 - 4*y**2))*SorH(14) + 9*x**2*SorH(15) - 12*x**4*SorH(15) - &
3*y**2*SorH(15) + 12*x**2*y**2*SorH(15) + r*x**4*SorHd(13) - &
r*x**2*y**2*SorHd(13) - 4*r*x**4*SorHd(14) + 4*r*x**2*y**2*SorHd(14) + &
3*r*x**2*(x - y)*(x + y)*SorHd(15)))/(2.*r)
	CASE (2421)
	! (2,2,2,-1)
	 grad(1) = (x*y*z*(-2*(-3 + 6*x**2 + 2*y**2 + 4*z**2)*SorH(13) + 4*(-3 + 6*x**2 - 2*y**2 &
+ 2*z**2)*SorH(14) + 6*SorH(15) - 12*x**2*SorH(15) + 12*y**2*SorH(15) + &
3*r*x**2*SorHd(13) + r*y**2*SorHd(13) + 2*r*z**2*SorHd(13) - &
6*r*x**2*SorHd(14) + 2*r*y**2*SorHd(14) - 2*r*z**2*SorHd(14) + 3*r*(x - y)*(x &
+ y)*SorHd(15)))/(2.*r)
	CASE (2422)
	! (2,2,2,0)
	 grad(1) = (sqrt3*x*(4*(-x**4 + y**4 + z**2 + 2*y**2*z**2 + x**2*(1 - 2*z**2))*SorH(13) + &
8*(-1 + 2*x**2 - 2*y**2)*z**2*SorH(14) - 4*x**2*SorH(15) + 4*x**4*SorH(15) - &
4*y**4*SorH(15) + 4*z**2*SorH(15) - 8*x**2*z**2*SorH(15) + &
8*y**2*z**2*SorH(15) + r*x**4*SorHd(13) - r*y**4*SorHd(13) + &
2*r*x**2*z**2*SorHd(13) - 2*r*y**2*z**2*SorHd(13) - 4*r*x**2*z**2*SorHd(14) + &
4*r*y**2*z**2*SorHd(14) - r*(x - y)*(x + y)*(x**2 + y**2 - &
2*z**2)*SorHd(15)))/(4.*r)
	CASE (2423)
	! (2,2,2,1)
	 grad(1) = -(z*((-4*x**4 + 3*y**2 + 2*z**2 + x**2*(3 - 12*y**2 - 8*z**2))*SorH(13) - &
2*(4*x**4 + 3*y**2 + z**2 - x**2*(3 + 12*y**2 + 4*z**2))*SorH(14) - &
9*x**2*SorH(15) + 12*x**4*SorH(15) + 3*y**2*SorH(15) - 12*x**2*y**2*SorH(15) + &
r*x**4*SorHd(13) + 3*r*x**2*y**2*SorHd(13) + 2*r*x**2*z**2*SorHd(13) + &
2*r*x**4*SorHd(14) - 6*r*x**2*y**2*SorHd(14) - 2*r*x**2*z**2*SorHd(14) - &
3*r*x**2*(x - y)*(x + y)*SorHd(15)))/(2.*r)
	CASE (2424)
	! (2,2,2,2)
	 grad(1) = (x*(-4*(x**4 + y**2 - 2*z**2 + (y**2 + 2*z**2)**2 + x**2*(-1 - 2*y**2 + &
4*z**2))*SorH(13) - 8*((-1 + 2*x**2)*z**2 + 2*y**2*(-2 + 4*x**2 + &
z**2))*SorH(14) + 12*x**2*SorH(15) - 12*x**4*SorH(15) - 12*y**2*SorH(15) + &
24*x**2*y**2*SorH(15) - 12*y**4*SorH(15) + r*x**4*SorHd(13) - &
2*r*x**2*y**2*SorHd(13) + r*y**4*SorHd(13) + 4*r*x**2*z**2*SorHd(13) + &
4*r*y**2*z**2*SorHd(13) + 4*r*z**4*SorHd(13) + 16*r*x**2*y**2*SorHd(14) + &
4*r*x**2*z**2*SorHd(14) + 4*r*y**2*z**2*SorHd(14) + 3*r*(x**2 - &
y**2)**2*SorHd(15)))/(4.*r)
  CASE DEFAULT
	 WRITE(*,*) "BUG: No case for i=",i," which encodes (l1,m1,l2,m2)=(",l1,",",m1,",",l2,",",m2,") !"
	 STOP
   END SELECT
!$OMP SECTION
   SELECT CASE(i)
	CASE (202)
	! (0,0,0,0)
	 grad(2) = y*SorHd(1)
	CASE (211)
	! (0,0,1,-1)
	 grad(2) = -(((-1 + y**2)*SorH(3))/r) + y**2*SorHd(3)
	CASE (212)
	! (0,0,1,0)
	 grad(2) = y*z*(-(SorH(3)/r) + SorHd(3))
	CASE (213)
	! (0,0,1,1)
	 grad(2) = x*y*(-(SorH(3)/r) + SorHd(3))
	CASE (220)
	! (0,0,2,-2)
	 grad(2) = (sqrt3*x*((1 - 2*y**2)*SorH(4) + r*y**2*SorHd(4)))/r
	CASE (221)
	! (0,0,2,-1)
	 grad(2) = (sqrt3*z*((1 - 2*y**2)*SorH(4) + r*y**2*SorHd(4)))/r
	CASE (222)
	! (0,0,2,0)
	 grad(2) = (y*(-1 + x**2 + y**2 - 2*z**2)*SorH(4))/r - (y*(x**2 + y**2 - &
2*z**2)*SorHd(4))/2.
	CASE (223)
	! (0,0,2,1)
	 grad(2) = (sqrt3*x*y*z*(-2*SorH(4) + r*SorHd(4)))/r
	CASE (224)
	! (0,0,2,2)
	 grad(2) = (sqrt3*y*(-2*(1 + x**2 - y**2)*SorH(4) + r*(x - y)*(x + y)*SorHd(4)))/(2.*r)
	CASE (1102)
	! (1,-1,0,0)
	 grad(2) = -(((-1 + y**2)*SorH(5))/r) + y**2*SorHd(5)
	CASE (1111)
	! (1,-1,1,-1)
	 grad(2) = (y*(-2*(-1 + y**2)*SorH(7) - (x**2 + z**2)*(2*SorH(6) - r*SorHd(6)) + &
r*y**2*SorHd(7)))/r
	CASE (1112)
	! (1,-1,1,0)
	 grad(2) = (z*((-1 + 2*y**2)*SorH(6) + SorH(7) + y**2*(-2*SorH(7) + r*(-SorHd(6) + &
SorHd(7)))))/r
	CASE (1113)
	! (1,-1,1,1)
	 grad(2) = (x*((-1 + 2*y**2)*SorH(6) + SorH(7) + y**2*(-2*SorH(7) + r*(-SorHd(6) + &
SorHd(7)))))/r
	CASE (1120)
	! (1,-1,2,-2)
	 grad(2) = (x*y*(-((2 + 3*x**2 - 3*y**2 + 3*z**2)*SorH(8)) + sqrt3*(2 - 3*y**2)*SorH(9) + &
r*(x**2 - y**2 + z**2)*SorHd(8) + r*sqrt3*y**2*SorHd(9)))/r
	CASE (1121)
	! (1,-1,2,-1)
	 grad(2) = (y*z*(-((2 + 3*x**2 - 3*y**2 + 3*z**2)*SorH(8)) + sqrt3*(2 - 3*y**2)*SorH(9) + &
r*(x**2 - y**2 + z**2)*SorHd(8) + r*sqrt3*y**2*SorHd(9)))/r
	CASE (1122)
	! (1,-1,2,0)
	 grad(2) = -(2*sqrt3*(1 - 3*y**2)*z**2*SorH(8) + x**2*SorH(9) + (-2*z**2 - 3*y**2*(-1 + &
x**2 + y**2 - 2*z**2))*SorH(9) + r*y**2*(2*sqrt3*z**2*SorHd(8) + (x**2 + y**2 &
- 2*z**2)*SorHd(9)))/(2.*r)
	CASE (1123)
	! (1,-1,2,1)
	 grad(2) = (x*z*((-2 + 6*y**2)*SorH(8) + sqrt3*(1 - 3*y**2)*SorH(9) + r*y**2*(-2*SorHd(8) &
+ sqrt3*SorHd(9))))/r
	CASE (1124)
	! (1,-1,2,2)
	 grad(2) = (2*(-1 + 3*y**2)*(2*x**2 + z**2)*SorH(8) + sqrt3*(x**2 - 3*(1 + x**2)*y**2 + &
3*y**4)*SorH(9) + r*y**2*(-2*(2*x**2 + z**2)*SorHd(8) + sqrt3*(x - y)*(x + &
y)*SorHd(9)))/(2.*r)
	CASE (1202)
	! (1,0,0,0)
	 grad(2) = y*z*(-(SorH(5)/r) + SorHd(5))
	CASE (1211)
	! (1,0,1,-1)
	 grad(2) = (z*((-1 + 2*y**2)*SorH(6) + SorH(7) + y**2*(-2*SorH(7) + r*(-SorHd(6) + &
SorHd(7)))))/r
	CASE (1212)
	! (1,0,1,0)
	 grad(2) = (y*(-2*(-1 + x**2 + y**2)*SorH(6) + r*(x**2 + y**2)*SorHd(6) + &
z**2*(-2*SorH(7) + r*SorHd(7))))/r
	CASE (1213)
	! (1,0,1,1)
	 grad(2) = (x*y*z*(2*SorH(6) - 2*SorH(7) + r*(-SorHd(6) + SorHd(7))))/r
	CASE (1220)
	! (1,0,2,-2)
	 grad(2) = (x*z*((-2 + 6*y**2)*SorH(8) + sqrt3*(1 - 3*y**2)*SorH(9) + r*y**2*(-2*SorHd(8) &
+ sqrt3*SorHd(9))))/r
	CASE (1221)
	! (1,0,2,-1)
	 grad(2) = ((x**2*(1 - 3*y**2) - z**2 + 3*y**2*(1 - y**2 + z**2))*SorH(8) + r*y**2*(x**2 &
+ y**2 - z**2)*SorHd(8) + sqrt3*z**2*((1 - 3*y**2)*SorH(9) + &
r*y**2*SorHd(9)))/r
	CASE (1222)
	! (1,0,2,0)
	 grad(2) = -(y*z*(2*sqrt3*(-2 + 3*x**2 + 3*y**2)*SorH(8) + (2 - 3*x**2 - 3*y**2 + &
6*z**2)*SorH(9) + r*(-2*sqrt3*(x**2 + y**2)*SorHd(8) + (x**2 + y**2 - &
2*z**2)*SorHd(9))))/(2.*r)
	CASE (1223)
	! (1,0,2,1)
	 grad(2) = (x*y*((2 - 3*x**2 - 3*y**2 + 3*z**2)*SorH(8) + r*(x**2 + y**2 - z**2)*SorHd(8) &
+ sqrt3*z**2*(-3*SorH(9) + r*SorHd(9))))/r
	CASE (1224)
	! (1,0,2,2)
	 grad(2) = (y*z*((4 + 6*x**2 - 6*y**2)*SorH(8) - 2*sqrt3*SorH(9) - (x - y)*(x + &
y)*(3*sqrt3*SorH(9) + 2*r*SorHd(8) - r*sqrt3*SorHd(9))))/(2.*r)
	CASE (1302)
	! (1,1,0,0)
	 grad(2) = x*y*(-(SorH(5)/r) + SorHd(5))
	CASE (1311)
	! (1,1,1,-1)
	 grad(2) = (x*((-1 + 2*y**2)*SorH(6) + SorH(7) + y**2*(-2*SorH(7) + r*(-SorHd(6) + &
SorHd(7)))))/r
	CASE (1312)
	! (1,1,1,0)
	 grad(2) = (x*y*z*(2*SorH(6) - 2*SorH(7) + r*(-SorHd(6) + SorHd(7))))/r
	CASE (1313)
	! (1,1,1,1)
	 grad(2) = (y*(-2*(-1 + y**2 + z**2)*SorH(6) + r*(y**2 + z**2)*SorHd(6) + &
x**2*(-2*SorH(7) + r*SorHd(7))))/r
	CASE (1320)
	! (1,1,2,-2)
	 grad(2) = ((x**2*(-1 + 3*y**2) + z**2 - 3*y**2*(-1 + y**2 + z**2))*SorH(8) + &
r*y**2*(-x**2 + y**2 + z**2)*SorHd(8) + sqrt3*x**2*((1 - 3*y**2)*SorH(9) + &
r*y**2*SorHd(9)))/r
	CASE (1321)
	! (1,1,2,-1)
	 grad(2) = (x*z*((-2 + 6*y**2)*SorH(8) + sqrt3*(1 - 3*y**2)*SorH(9) + r*y**2*(-2*SorHd(8) &
+ sqrt3*SorHd(9))))/r
	CASE (1322)
	! (1,1,2,0)
	 grad(2) = -(x*y*((2 - 3*x**2 - 3*y**2)*SorH(9) + 2*z**2*(-3*sqrt3*SorH(8) + 3*SorH(9) + &
r*sqrt3*SorHd(8)) + r*(x**2 + y**2 - 2*z**2)*SorHd(9)))/(2.*r)
	CASE (1323)
	! (1,1,2,1)
	 grad(2) = (y*z*((2 + 3*x**2 - 3*y**2 - 3*z**2)*SorH(8) + r*(-x**2 + y**2 + &
z**2)*SorHd(8) + sqrt3*x**2*(-3*SorH(9) + r*SorHd(9))))/r
	CASE (1324)
	! (1,1,2,2)
	 grad(2) = (x*y*(-2*(-4 + 6*y**2 + 3*z**2)*SorH(8) + sqrt3*(-2 - 3*x**2 + 3*y**2)*SorH(9) &
+ 2*r*(2*y**2 + z**2)*SorHd(8) + r*sqrt3*(x - y)*(x + y)*SorHd(9)))/(2.*r)
	CASE (2002)
	! (2,-2,0,0)
	 grad(2) = (sqrt3*x*((1 - 2*y**2)*SorH(10) + r*y**2*SorHd(10)))/r
	CASE (2011)
	! (2,-2,1,-1)
	 grad(2) = (x*y*(-((2 + 3*x**2 - 3*y**2 + 3*z**2)*SorH(11)) + sqrt3*(2 - 3*y**2)*SorH(12) &
+ r*(x**2 - y**2 + z**2)*SorHd(11) + r*sqrt3*y**2*SorHd(12)))/r
	CASE (2012)
	! (2,-2,1,0)
	 grad(2) = (x*z*((-2 + 6*y**2)*SorH(11) + sqrt3*(1 - 3*y**2)*SorH(12) + &
r*y**2*(-2*SorHd(11) + sqrt3*SorHd(12))))/r
	CASE (2013)
	! (2,-2,1,1)
	 grad(2) = ((x**2*(-1 + 3*y**2) + z**2 - 3*y**2*(-1 + y**2 + z**2))*SorH(11) + &
r*y**2*(-x**2 + y**2 + z**2)*SorHd(11) + sqrt3*x**2*((1 - 3*y**2)*SorH(12) + &
r*y**2*SorHd(12)))/r
	CASE (2020)
	! (2,-2,2,-2)
	 grad(2) = (y*(-2*(x**2 + z**2)*(-1 + 2*y**2 + 2*z**2)*SorH(13) + 2*(-2*(x - y)*(x + &
y)*(1 + x**2 - y**2) + (1 - 2*x**2 - 2*y**2)*z**2)*SorH(14) + 6*x**2*SorH(15) &
- 12*x**2*y**2*SorH(15) + r*x**2*y**2*SorHd(13) + r*x**2*z**2*SorHd(13) + &
r*y**2*z**2*SorHd(13) + r*z**4*SorHd(13) + r*x**4*SorHd(14) - &
2*r*x**2*y**2*SorHd(14) + r*y**4*SorHd(14) + r*x**2*z**2*SorHd(14) + &
r*y**2*z**2*SorHd(14) + 3*r*x**2*y**2*SorHd(15)))/r
	CASE (2021)
	! (2,-2,2,-1)
	 grad(2) = (x*y*z*(4*(x**2 + z**2)*SorH(13) - 2*(3 + 2*x**2 - 6*y**2 + 2*z**2)*SorH(14) + &
6*SorH(15) - 12*y**2*SorH(15) - r*(x**2 + z**2)*SorHd(13) + r*x**2*SorHd(14) - &
3*r*y**2*SorHd(14) + r*z**2*SorHd(14) + 3*r*y**2*SorHd(15)))/r
	CASE (2022)
	! (2,-2,2,0)
	 grad(2) = (sqrt3*x*((x**2*(1 - 4*y**2) + 2*z**2 + y**2*(3 - 4*y**2 - 8*z**2))*SorH(13) + &
4*(-1 + 4*y**2)*z**2*SorH(14) - x**2*SorH(15) - 3*y**2*SorH(15) + &
4*x**2*y**2*SorH(15) + 4*y**4*SorH(15) + 2*z**2*SorH(15) - &
8*y**2*z**2*SorH(15) + r*x**2*y**2*SorHd(13) + r*y**4*SorHd(13) + &
2*r*y**2*z**2*SorHd(13) - 4*r*y**2*z**2*SorHd(14) - r*y**2*(x**2 + y**2 - &
2*z**2)*SorHd(15)))/(2.*r)
	CASE (2023)
	! (2,-2,2,1)
	 grad(2) = (z*((4*y**4 - z**2 + y**2*(-3 + 4*z**2))*SorH(13) + (3*x**2*(-1 + 4*y**2) + &
z**2 + y**2*(3 - 4*y**2 - 4*z**2))*SorH(14) + 3*x**2*SorH(15) - &
12*x**2*y**2*SorH(15) - r*y**4*SorHd(13) - r*y**2*z**2*SorHd(13) - &
3*r*x**2*y**2*SorHd(14) + r*y**4*SorHd(14) + r*y**2*z**2*SorHd(14) + &
3*r*x**2*y**2*SorHd(15)))/r
	CASE (2024)
	! (2,-2,2,2)
	 grad(2) = (x*((x**2 - (3 + 4*x**2)*y**2 + 4*y**4)*SorH(13) + 4*(-x**2 + (3 + &
4*x**2)*y**2 - 4*y**4)*SorH(14) + 3*x**2*SorH(15) - 9*y**2*SorH(15) - &
12*x**2*y**2*SorH(15) + 12*y**4*SorH(15) + r*x**2*y**2*SorHd(13) - &
r*y**4*SorHd(13) - 4*r*x**2*y**2*SorHd(14) + 4*r*y**4*SorHd(14) + 3*r*(x - &
y)*y**2*(x + y)*SorHd(15)))/(2.*r)
	CASE (2102)
	! (2,-1,0,0)
	 grad(2) = (sqrt3*z*((1 - 2*y**2)*SorH(10) + r*y**2*SorHd(10)))/r
	CASE (2111)
	! (2,-1,1,-1)
	 grad(2) = (y*z*(-((2 + 3*x**2 - 3*y**2 + 3*z**2)*SorH(11)) + sqrt3*(2 - 3*y**2)*SorH(12) &
+ r*(x**2 - y**2 + z**2)*SorHd(11) + r*sqrt3*y**2*SorHd(12)))/r
	CASE (2112)
	! (2,-1,1,0)
	 grad(2) = ((x**2*(1 - 3*y**2) - z**2 + 3*y**2*(1 - y**2 + z**2))*SorH(11) + r*y**2*(x**2 &
+ y**2 - z**2)*SorHd(11) + sqrt3*z**2*((1 - 3*y**2)*SorH(12) + &
r*y**2*SorHd(12)))/r
	CASE (2113)
	! (2,-1,1,1)
	 grad(2) = (x*z*((-2 + 6*y**2)*SorH(11) + sqrt3*(1 - 3*y**2)*SorH(12) + &
r*y**2*(-2*SorHd(11) + sqrt3*SorHd(12))))/r
	CASE (2120)
	! (2,-1,2,-2)
	 grad(2) = (x*y*z*(4*(x**2 + z**2)*SorH(13) - 2*(3 + 2*x**2 - 6*y**2 + 2*z**2)*SorH(14) + &
6*SorH(15) - 12*y**2*SorH(15) - r*(x**2 + z**2)*SorHd(13) + r*x**2*SorHd(14) - &
3*r*y**2*SorHd(14) + r*z**2*SorHd(14) + 3*r*y**2*SorHd(15)))/r
	CASE (2121)
	! (2,-1,2,-1)
	 grad(2) = (y*(-2*(-1 + 2*x**2 + 2*y**2)*(x**2 + z**2)*SorH(13) + 2*(x**2*(1 - 2*y**2 - &
2*z**2) - 2*(y - z)*(y + z)*(-1 + y**2 - z**2))*SorH(14) + 6*z**2*SorH(15) - &
12*y**2*z**2*SorH(15) + r*x**4*SorHd(13) + r*x**2*y**2*SorHd(13) + &
r*x**2*z**2*SorHd(13) + r*y**2*z**2*SorHd(13) + r*x**2*y**2*SorHd(14) + &
r*y**4*SorHd(14) + r*x**2*z**2*SorHd(14) - 2*r*y**2*z**2*SorHd(14) + &
r*z**4*SorHd(14) + 3*r*y**2*z**2*SorHd(15)))/r
	CASE (2122)
	! (2,-1,2,0)
	 grad(2) = (sqrt3*z*((-x**2 + (-3 + 4*x**2)*y**2 + 4*y**4)*SorH(13) + 2*(x**2 + 3*y**2 - &
4*x**2*y**2 - 4*y**4 + (-1 + 4*y**2)*z**2)*SorH(14) - x**2*SorH(15) - &
3*y**2*SorH(15) + 4*x**2*y**2*SorH(15) + 4*y**4*SorH(15) + 2*z**2*SorH(15) - &
8*y**2*z**2*SorH(15) - r*x**2*y**2*SorHd(13) - r*y**4*SorHd(13) + &
2*r*x**2*y**2*SorHd(14) + 2*r*y**4*SorHd(14) - 2*r*y**2*z**2*SorHd(14) - &
r*y**2*(x**2 + y**2 - 2*z**2)*SorHd(15)))/(2.*r)
	CASE (2123)
	! (2,-1,2,1)
	 grad(2) = (x*((-x**2 + (-3 + 4*x**2)*y**2 + 4*y**4)*SorH(13) + (x**2 + 3*y**2 - &
4*x**2*y**2 - 4*y**4 + 3*(-1 + 4*y**2)*z**2)*SorH(14) + 3*z**2*SorH(15) - &
12*y**2*z**2*SorH(15) - r*x**2*y**2*SorHd(13) - r*y**4*SorHd(13) + &
r*x**2*y**2*SorHd(14) + r*y**4*SorHd(14) - 3*r*y**2*z**2*SorHd(14) + &
3*r*y**2*z**2*SorHd(15)))/r
	CASE (2124)
	! (2,-1,2,2)
	 grad(2) = (z*((-4*y**4 + x**2*(3 - 12*y**2) + 2*z**2 + y**2*(3 - 8*z**2))*SorH(13) + &
(-8*y**4 + 6*x**2*(-1 + 4*y**2) - 2*z**2 + y**2*(6 + 8*z**2))*SorH(14) + &
3*x**2*SorH(15) - 9*y**2*SorH(15) - 12*x**2*y**2*SorH(15) + 12*y**4*SorH(15) + &
3*r*x**2*y**2*SorHd(13) + r*y**4*SorHd(13) + 2*r*y**2*z**2*SorHd(13) - &
6*r*x**2*y**2*SorHd(14) + 2*r*y**4*SorHd(14) - 2*r*y**2*z**2*SorHd(14) + &
3*r*(x - y)*y**2*(x + y)*SorHd(15)))/(2.*r)
	CASE (2202)
	! (2,0,0,0)
	 grad(2) = (y*(-1 + x**2 + y**2 - 2*z**2)*SorH(10))/r - (y*(x**2 + y**2 - &
2*z**2)*SorHd(10))/2.
	CASE (2211)
	! (2,0,1,-1)
	 grad(2) = -(2*sqrt3*(1 - 3*y**2)*z**2*SorH(11) + (-3*y**4 + x**2*(1 - 3*y**2) - 2*z**2 + &
y**2*(3 + 6*z**2))*SorH(12) + r*y**2*(2*sqrt3*z**2*SorHd(11) + (x**2 + y**2 - &
2*z**2)*SorHd(12)))/(2.*r)
	CASE (2212)
	! (2,0,1,0)
	 grad(2) = -(y*z*(2*sqrt3*(-2 + 3*x**2 + 3*y**2)*SorH(11) + (2 - 3*x**2 - 3*y**2 + &
6*z**2)*SorH(12) + r*(-2*sqrt3*(x**2 + y**2)*SorHd(11) + (x**2 + y**2 - &
2*z**2)*SorHd(12))))/(2.*r)
	CASE (2213)
	! (2,0,1,1)
	 grad(2) = -(x*y*((2 - 3*x**2 - 3*y**2)*SorH(12) + 2*z**2*(-3*sqrt3*SorH(11) + 3*SorH(12) &
+ r*sqrt3*SorHd(11)) + r*(x**2 + y**2 - 2*z**2)*SorHd(12)))/(2.*r)
	CASE (2220)
	! (2,0,2,-2)
	 grad(2) = (sqrt3*x*((x**2*(1 - 4*y**2) + 2*z**2 + y**2*(3 - 4*y**2 - 8*z**2))*SorH(13) + &
4*(-1 + 4*y**2)*z**2*SorH(14) - x**2*SorH(15) - 3*y**2*SorH(15) + &
4*x**2*y**2*SorH(15) + 4*y**4*SorH(15) + 2*z**2*SorH(15) - &
8*y**2*z**2*SorH(15) + r*x**2*y**2*SorHd(13) + r*y**4*SorHd(13) + &
2*r*y**2*z**2*SorHd(13) - 4*r*y**2*z**2*SorHd(14) - r*y**2*(x**2 + y**2 - &
2*z**2)*SorHd(15)))/(2.*r)
	CASE (2221)
	! (2,0,2,-1)
	 grad(2) = (sqrt3*z*((-x**2 + (-3 + 4*x**2)*y**2 + 4*y**4)*SorH(13) + 2*(x**2 + 3*y**2 - &
4*x**2*y**2 - 4*y**4 + (-1 + 4*y**2)*z**2)*SorH(14) - x**2*SorH(15) - &
3*y**2*SorH(15) + 4*x**2*y**2*SorH(15) + 4*y**4*SorH(15) + 2*z**2*SorH(15) - &
8*y**2*z**2*SorH(15) - r*x**2*y**2*SorHd(13) - r*y**4*SorHd(13) + &
2*r*x**2*y**2*SorHd(14) + 2*r*y**4*SorHd(14) - 2*r*y**2*z**2*SorHd(14) - &
r*y**2*(x**2 + y**2 - 2*z**2)*SorHd(15)))/(2.*r)
	CASE (2222)
	! (2,0,2,0)
	 grad(2) = (y*(-12*(-1 + x**2 + y**2)*(x**2 + y**2)*SorH(13) - 24*(-1 + 2*x**2 + &
2*y**2)*z**2*SorH(14) + 4*x**2*SorH(15) - 4*x**4*SorH(15) + 4*y**2*SorH(15) - &
8*x**2*y**2*SorH(15) - 4*y**4*SorH(15) - 8*z**2*SorH(15) + &
16*x**2*z**2*SorH(15) + 16*y**2*z**2*SorH(15) - 16*z**4*SorH(15) + &
3*r*x**4*SorHd(13) + 6*r*x**2*y**2*SorHd(13) + 3*r*y**4*SorHd(13) + &
12*r*x**2*z**2*SorHd(14) + 12*r*y**2*z**2*SorHd(14) + r*(x**2 + y**2 - &
2*z**2)**2*SorHd(15)))/(4.*r)
	CASE (2223)
	! (2,0,2,1)
	 grad(2) = (sqrt3*x*y*z*((-2 + 4*x**2 + 4*y**2)*SorH(13) + (4 - 8*x**2 - 8*y**2 + &
8*z**2)*SorH(14) - 2*SorH(15) + 4*x**2*SorH(15) + 4*y**2*SorH(15) - &
8*z**2*SorH(15) - r*x**2*SorHd(13) - r*y**2*SorHd(13) + 2*r*x**2*SorHd(14) + &
2*r*y**2*SorHd(14) - 2*r*z**2*SorHd(14) - r*(x**2 + y**2 - &
2*z**2)*SorHd(15)))/(2.*r)
	CASE (2224)
	! (2,0,2,2)
	 grad(2) = (sqrt3*y*(-4*(x**4 + y**2 - y**4 + (1 + 2*x**2 - 2*y**2)*z**2)*SorH(13) + 8*(1 &
+ 2*x**2 - 2*y**2)*z**2*SorH(14) + 4*x**4*SorH(15) + 4*y**2*SorH(15) - &
4*y**4*SorH(15) - 4*z**2*SorH(15) - 8*x**2*z**2*SorH(15) + &
8*y**2*z**2*SorH(15) + r*x**4*SorHd(13) - r*y**4*SorHd(13) + &
2*r*x**2*z**2*SorHd(13) - 2*r*y**2*z**2*SorHd(13) - 4*r*x**2*z**2*SorHd(14) + &
4*r*y**2*z**2*SorHd(14) - r*(x - y)*(x + y)*(x**2 + y**2 - &
2*z**2)*SorHd(15)))/(4.*r)
	CASE (2302)
	! (2,1,0,0)
	 grad(2) = (sqrt3*x*y*z*(-2*SorH(10) + r*SorHd(10)))/r
	CASE (2311)
	! (2,1,1,-1)
	 grad(2) = (x*z*((-2 + 6*y**2)*SorH(11) + sqrt3*(1 - 3*y**2)*SorH(12) + &
r*y**2*(-2*SorHd(11) + sqrt3*SorHd(12))))/r
	CASE (2312)
	! (2,1,1,0)
	 grad(2) = (x*y*((2 - 3*x**2 - 3*y**2 + 3*z**2)*SorH(11) + r*(x**2 + y**2 - &
z**2)*SorHd(11) + sqrt3*z**2*(-3*SorH(12) + r*SorHd(12))))/r
	CASE (2313)
	! (2,1,1,1)
	 grad(2) = (y*z*((2 + 3*x**2 - 3*y**2 - 3*z**2)*SorH(11) + r*(-x**2 + y**2 + &
z**2)*SorHd(11) + sqrt3*x**2*(-3*SorH(12) + r*SorHd(12))))/r
	CASE (2320)
	! (2,1,2,-2)
	 grad(2) = (z*((4*y**4 - z**2 + y**2*(-3 + 4*z**2))*SorH(13) + (3*x**2*(-1 + 4*y**2) + &
z**2 + y**2*(3 - 4*y**2 - 4*z**2))*SorH(14) + 3*x**2*SorH(15) - &
12*x**2*y**2*SorH(15) - r*y**4*SorHd(13) - r*y**2*z**2*SorHd(13) - &
3*r*x**2*y**2*SorHd(14) + r*y**4*SorHd(14) + r*y**2*z**2*SorHd(14) + &
3*r*x**2*y**2*SorHd(15)))/r
	CASE (2321)
	! (2,1,2,-1)
	 grad(2) = (x*((-x**2 + (-3 + 4*x**2)*y**2 + 4*y**4)*SorH(13) + (x**2 + 3*y**2 - &
4*x**2*y**2 - 4*y**4 + 3*(-1 + 4*y**2)*z**2)*SorH(14) + 3*z**2*SorH(15) - &
12*y**2*z**2*SorH(15) - r*x**2*y**2*SorHd(13) - r*y**4*SorHd(13) + &
r*x**2*y**2*SorHd(14) + r*y**4*SorHd(14) - 3*r*y**2*z**2*SorHd(14) + &
3*r*y**2*z**2*SorHd(15)))/r
	CASE (2322)
	! (2,1,2,0)
	 grad(2) = (sqrt3*x*y*z*((-2 + 4*x**2 + 4*y**2)*SorH(13) + (4 - 8*x**2 - 8*y**2 + &
8*z**2)*SorH(14) - 2*SorH(15) + 4*x**2*SorH(15) + 4*y**2*SorH(15) - &
8*z**2*SorH(15) - r*x**2*SorHd(13) - r*y**2*SorHd(13) + 2*r*x**2*SorHd(14) + &
2*r*y**2*SorHd(14) - 2*r*z**2*SorHd(14) - r*(x**2 + y**2 - &
2*z**2)*SorHd(15)))/(2.*r)
	CASE (2323)
	! (2,1,2,1)
	 grad(2) = (y*(2*(z**2 + x**2*(1 - 2*y**2 - 2*z**2) - 2*y**2*(-1 + y**2 + z**2))*SorH(13) &
+ 2*(-2*x**4 + z**2 - 2*z**2*(y**2 + z**2) + x**2*(1 - 2*y**2 + &
4*z**2))*SorH(14) - 12*x**2*z**2*SorH(15) + r*x**2*y**2*SorHd(13) + &
r*y**4*SorHd(13) + r*x**2*z**2*SorHd(13) + r*y**2*z**2*SorHd(13) + &
r*x**4*SorHd(14) + r*x**2*y**2*SorHd(14) - 2*r*x**2*z**2*SorHd(14) + &
r*y**2*z**2*SorHd(14) + r*z**4*SorHd(14) + 3*r*x**2*z**2*SorHd(15)))/r
	CASE (2324)
	! (2,1,2,2)
	 grad(2) = (x*y*z*(2*(-3 + 2*x**2 + 6*y**2 + 4*z**2)*SorH(13) + 4*(3 + 2*x**2 - 6*y**2 - &
2*z**2)*SorH(14) - 6*SorH(15) - 12*x**2*SorH(15) + 12*y**2*SorH(15) - &
r*x**2*SorHd(13) - 3*r*y**2*SorHd(13) - 2*r*z**2*SorHd(13) - &
2*r*x**2*SorHd(14) + 6*r*y**2*SorHd(14) + 2*r*z**2*SorHd(14) + 3*r*(x - y)*(x &
+ y)*SorHd(15)))/(2.*r)
	CASE (2402)
	! (2,2,0,0)
	 grad(2) = (sqrt3*y*(-2*(1 + x**2 - y**2)*SorH(10) + r*(x - y)*(x + y)*SorHd(10)))/(2.*r)
	CASE (2411)
	! (2,2,1,-1)
	 grad(2) = (2*(-1 + 3*y**2)*(2*x**2 + z**2)*SorH(11) + sqrt3*(x**2 - 3*(1 + x**2)*y**2 + &
3*y**4)*SorH(12) + r*y**2*(-2*(2*x**2 + z**2)*SorHd(11) + sqrt3*(x - y)*(x + &
y)*SorHd(12)))/(2.*r)
	CASE (2412)
	! (2,2,1,0)
	 grad(2) = -(y*z*((-4 - 6*x**2 + 6*y**2)*SorH(11) + sqrt3*(2 + 3*x**2 - 3*y**2)*SorH(12) &
+ r*(x - y)*(x + y)*(2*SorHd(11) - sqrt3*SorHd(12))))/(2.*r)
	CASE (2413)
	! (2,2,1,1)
	 grad(2) = (x*y*(-2*(-4 + 6*y**2 + 3*z**2)*SorH(11) + sqrt3*(-2 - 3*x**2 + &
3*y**2)*SorH(12) + 2*r*(2*y**2 + z**2)*SorHd(11) + r*sqrt3*(x - y)*(x + &
y)*SorHd(12)))/(2.*r)
	CASE (2420)
	! (2,2,2,-2)
	 grad(2) = (x*((x**2 - (3 + 4*x**2)*y**2 + 4*y**4)*SorH(13) + 4*(-x**2 + (3 + &
4*x**2)*y**2 - 4*y**4)*SorH(14) + 3*x**2*SorH(15) - 9*y**2*SorH(15) - &
12*x**2*y**2*SorH(15) + 12*y**4*SorH(15) + r*x**2*y**2*SorHd(13) - &
r*y**4*SorHd(13) - 4*r*x**2*y**2*SorHd(14) + 4*r*y**4*SorHd(14) + 3*r*(x - &
y)*y**2*(x + y)*SorHd(15)))/(2.*r)
	CASE (2421)
	! (2,2,2,-1)
	 grad(2) = (z*((-4*y**4 + x**2*(3 - 12*y**2) + 2*z**2 + y**2*(3 - 8*z**2))*SorH(13) + &
(-8*y**4 + 6*x**2*(-1 + 4*y**2) - 2*z**2 + y**2*(6 + 8*z**2))*SorH(14) + &
3*x**2*SorH(15) - 9*y**2*SorH(15) - 12*x**2*y**2*SorH(15) + 12*y**4*SorH(15) + &
3*r*x**2*y**2*SorHd(13) + r*y**4*SorHd(13) + 2*r*y**2*z**2*SorHd(13) - &
6*r*x**2*y**2*SorHd(14) + 2*r*y**4*SorHd(14) - 2*r*y**2*z**2*SorHd(14) + &
3*r*(x - y)*y**2*(x + y)*SorHd(15)))/(2.*r)
	CASE (2422)
	! (2,2,2,0)
	 grad(2) = (sqrt3*y*(-4*(x**4 + y**2 - y**4 + (1 + 2*x**2 - 2*y**2)*z**2)*SorH(13) + 8*(1 &
+ 2*x**2 - 2*y**2)*z**2*SorH(14) + 4*x**4*SorH(15) + 4*y**2*SorH(15) - &
4*y**4*SorH(15) - 4*z**2*SorH(15) - 8*x**2*z**2*SorH(15) + &
8*y**2*z**2*SorH(15) + r*x**4*SorHd(13) - r*y**4*SorHd(13) + &
2*r*x**2*z**2*SorHd(13) - 2*r*y**2*z**2*SorHd(13) - 4*r*x**2*z**2*SorHd(14) + &
4*r*y**2*z**2*SorHd(14) - r*(x - y)*(x + y)*(x**2 + y**2 - &
2*z**2)*SorHd(15)))/(4.*r)
	CASE (2423)
	! (2,2,2,1)
	 grad(2) = (x*y*z*(2*(-3 + 2*x**2 + 6*y**2 + 4*z**2)*SorH(13) + 4*(3 + 2*x**2 - 6*y**2 - &
2*z**2)*SorH(14) - 6*SorH(15) - 12*x**2*SorH(15) + 12*y**2*SorH(15) - &
r*x**2*SorHd(13) - 3*r*y**2*SorHd(13) - 2*r*z**2*SorHd(13) - &
2*r*x**2*SorHd(14) + 6*r*y**2*SorHd(14) + 2*r*z**2*SorHd(14) + 3*r*(x - y)*(x &
+ y)*SorHd(15)))/(2.*r)
	CASE (2424)
	! (2,2,2,2)
	 grad(2) = (y*(-4*(x**4 + (-1 + y**2 + 2*z**2)*(y**2 + 2*z**2) + x**2*(1 - 2*y**2 + &
4*z**2))*SorH(13) - 8*((-1 + 2*y**2)*z**2 + 2*x**2*(-2 + 4*y**2 + &
z**2))*SorH(14) - 12*x**2*SorH(15) - 12*x**4*SorH(15) + 12*y**2*SorH(15) + &
24*x**2*y**2*SorH(15) - 12*y**4*SorH(15) + r*x**4*SorHd(13) - &
2*r*x**2*y**2*SorHd(13) + r*y**4*SorHd(13) + 4*r*x**2*z**2*SorHd(13) + &
4*r*y**2*z**2*SorHd(13) + 4*r*z**4*SorHd(13) + 16*r*x**2*y**2*SorHd(14) + &
4*r*x**2*z**2*SorHd(14) + 4*r*y**2*z**2*SorHd(14) + 3*r*(x**2 - &
y**2)**2*SorHd(15)))/(4.*r)
  CASE DEFAULT
	 WRITE(*,*) "BUG: No case for i=",i," which encodes (l1,m1,l2,m2)=(",l1,",",m1,",",l2,",",m2,") !"
	 STOP
   END SELECT
!$OMP SECTION
   SELECT CASE(i)
	CASE (202)
	! (0,0,0,0)
	 grad(3) = z*SorHd(1)
	CASE (211)
	! (0,0,1,-1)
	 grad(3) = y*z*(-(SorH(3)/r) + SorHd(3))
	CASE (212)
	! (0,0,1,0)
	 grad(3) = -(((-1 + z**2)*SorH(3))/r) + z**2*SorHd(3)
	CASE (213)
	! (0,0,1,1)
	 grad(3) = x*z*(-(SorH(3)/r) + SorHd(3))
	CASE (220)
	! (0,0,2,-2)
	 grad(3) = (sqrt3*x*y*z*(-2*SorH(4) + r*SorHd(4)))/r
	CASE (221)
	! (0,0,2,-1)
	 grad(3) = (sqrt3*y*((1 - 2*z**2)*SorH(4) + r*z**2*SorHd(4)))/r
	CASE (222)
	! (0,0,2,0)
	 grad(3) = (z*(2*(2 + x**2 + y**2 - 2*z**2)*SorH(4) - r*(x**2 + y**2 - &
2*z**2)*SorHd(4)))/(2.*r)
	CASE (223)
	! (0,0,2,1)
	 grad(3) = (sqrt3*x*((1 - 2*z**2)*SorH(4) + r*z**2*SorHd(4)))/r
	CASE (224)
	! (0,0,2,2)
	 grad(3) = (sqrt3*(x - y)*(x + y)*z*(-2*SorH(4) + r*SorHd(4)))/(2.*r)
	CASE (1102)
	! (1,-1,0,0)
	 grad(3) = y*z*(-(SorH(5)/r) + SorHd(5))
	CASE (1111)
	! (1,-1,1,-1)
	 grad(3) = (z*(-2*(-1 + x**2 + z**2)*SorH(6) + r*(x**2 + z**2)*SorHd(6) + &
y**2*(-2*SorH(7) + r*SorHd(7))))/r
	CASE (1112)
	! (1,-1,1,0)
	 grad(3) = (y*((-1 + 2*z**2)*SorH(6) + SorH(7) + z**2*(-2*SorH(7) + r*(-SorHd(6) + &
SorHd(7)))))/r
	CASE (1113)
	! (1,-1,1,1)
	 grad(3) = (x*y*z*(2*SorH(6) - 2*SorH(7) + r*(-SorHd(6) + SorHd(7))))/r
	CASE (1120)
	! (1,-1,2,-2)
	 grad(3) = (x*z*((2 - 3*x**2 + 3*y**2 - 3*z**2)*SorH(8) + r*(x**2 - y**2 + z**2)*SorHd(8) &
+ sqrt3*y**2*(-3*SorH(9) + r*SorHd(9))))/r
	CASE (1121)
	! (1,-1,2,-1)
	 grad(3) = ((x**2 - y**2 + 3*(1 - x**2 + y**2)*z**2 - 3*z**4)*SorH(8) + r*z**2*(x**2 - &
y**2 + z**2)*SorHd(8) + sqrt3*y**2*((1 - 3*z**2)*SorH(9) + r*z**2*SorHd(9)))/r
	CASE (1122)
	! (1,-1,2,0)
	 grad(3) = (y*z*(2*sqrt3*(-2 + 3*z**2)*SorH(8) + (4 + 3*x**2 + 3*y**2 - 6*z**2)*SorH(9) - &
r*(2*sqrt3*z**2*SorHd(8) + (x**2 + y**2 - 2*z**2)*SorHd(9))))/(2.*r)
	CASE (1123)
	! (1,-1,2,1)
	 grad(3) = (x*y*((-2 + 6*z**2)*SorH(8) + sqrt3*(1 - 3*z**2)*SorH(9) + r*z**2*(-2*SorHd(8) &
+ sqrt3*SorHd(9))))/r
	CASE (1124)
	! (1,-1,2,2)
	 grad(3) = (y*z*(2*(-2 + 6*x**2 + 3*z**2)*SorH(8) + 3*sqrt3*(-x**2 + y**2)*SorH(9) - &
2*r*(2*x**2 + z**2)*SorHd(8) + r*sqrt3*(x - y)*(x + y)*SorHd(9)))/(2.*r)
	CASE (1202)
	! (1,0,0,0)
	 grad(3) = -(((-1 + z**2)*SorH(5))/r) + z**2*SorHd(5)
	CASE (1211)
	! (1,0,1,-1)
	 grad(3) = (y*((-1 + 2*z**2)*SorH(6) + SorH(7) + z**2*(-2*SorH(7) + r*(-SorHd(6) + &
SorHd(7)))))/r
	CASE (1212)
	! (1,0,1,0)
	 grad(3) = (z*(-2*(-1 + z**2)*SorH(7) - (x**2 + y**2)*(2*SorH(6) - r*SorHd(6)) + &
r*z**2*SorHd(7)))/r
	CASE (1213)
	! (1,0,1,1)
	 grad(3) = (x*((-1 + 2*z**2)*SorH(6) + SorH(7) + z**2*(-2*SorH(7) + r*(-SorHd(6) + &
SorHd(7)))))/r
	CASE (1220)
	! (1,0,2,-2)
	 grad(3) = (x*y*((-2 + 6*z**2)*SorH(8) + sqrt3*(1 - 3*z**2)*SorH(9) + r*z**2*(-2*SorHd(8) &
+ sqrt3*SorHd(9))))/r
	CASE (1221)
	! (1,0,2,-1)
	 grad(3) = (y*z*(-((2 + 3*x**2 + 3*y**2 - 3*z**2)*SorH(8)) + sqrt3*(2 - 3*z**2)*SorH(9) + &
r*(x**2 + y**2 - z**2)*SorHd(8) + r*sqrt3*z**2*SorHd(9)))/r
	CASE (1222)
	! (1,0,2,0)
	 grad(3) = (-2*sqrt3*(x**2 + y**2)*(-1 + 3*z**2)*SorH(8) + (-x**2 - y**2 + 3*(2 + x**2 + &
y**2)*z**2 - 6*z**4)*SorH(9) + r*(x**2 + y**2)*z**2*(2*sqrt3*SorHd(8) - &
SorHd(9)) + 2*r*z**4*SorHd(9))/(2.*r)
	CASE (1223)
	! (1,0,2,1)
	 grad(3) = (x*z*(-((2 + 3*x**2 + 3*y**2 - 3*z**2)*SorH(8)) + sqrt3*(2 - 3*z**2)*SorH(9) + &
r*(x**2 + y**2 - z**2)*SorHd(8) + r*sqrt3*z**2*SorHd(9)))/r
	CASE (1224)
	! (1,0,2,2)
	 grad(3) = ((x - y)*(x + y)*((-2 + 6*z**2)*SorH(8) + sqrt3*(1 - 3*z**2)*SorH(9) + &
r*z**2*(-2*SorHd(8) + sqrt3*SorHd(9))))/(2.*r)
	CASE (1302)
	! (1,1,0,0)
	 grad(3) = x*z*(-(SorH(5)/r) + SorHd(5))
	CASE (1311)
	! (1,1,1,-1)
	 grad(3) = (x*y*z*(2*SorH(6) - 2*SorH(7) + r*(-SorHd(6) + SorHd(7))))/r
	CASE (1312)
	! (1,1,1,0)
	 grad(3) = (x*((-1 + 2*z**2)*SorH(6) + SorH(7) + z**2*(-2*SorH(7) + r*(-SorHd(6) + &
SorHd(7)))))/r
	CASE (1313)
	! (1,1,1,1)
	 grad(3) = (z*(-2*(-1 + y**2 + z**2)*SorH(6) + r*(y**2 + z**2)*SorHd(6) + &
x**2*(-2*SorH(7) + r*SorHd(7))))/r
	CASE (1320)
	! (1,1,2,-2)
	 grad(3) = (y*z*((2 + 3*x**2 - 3*y**2 - 3*z**2)*SorH(8) + r*(-x**2 + y**2 + &
z**2)*SorHd(8) + sqrt3*x**2*(-3*SorH(9) + r*SorHd(9))))/r
	CASE (1321)
	! (1,1,2,-1)
	 grad(3) = (x*y*((-2 + 6*z**2)*SorH(8) + sqrt3*(1 - 3*z**2)*SorH(9) + r*z**2*(-2*SorHd(8) &
+ sqrt3*SorHd(9))))/r
	CASE (1322)
	! (1,1,2,0)
	 grad(3) = (x*z*(2*sqrt3*(-2 + 3*z**2)*SorH(8) + (4 + 3*x**2 + 3*y**2 - 6*z**2)*SorH(9) - &
r*(2*sqrt3*z**2*SorHd(8) + (x**2 + y**2 - 2*z**2)*SorHd(9))))/(2.*r)
	CASE (1323)
	! (1,1,2,1)
	 grad(3) = ((-x**2 + y**2 + 3*(1 + x**2 - y**2)*z**2 - 3*z**4)*SorH(8) + r*z**2*(-x**2 + &
y**2 + z**2)*SorHd(8) + sqrt3*x**2*((1 - 3*z**2)*SorH(9) + r*z**2*SorHd(9)))/r
	CASE (1324)
	! (1,1,2,2)
	 grad(3) = (x*z*(-2*(-2 + 6*y**2 + 3*z**2)*SorH(8) + 3*sqrt3*(-x**2 + y**2)*SorH(9) + &
2*r*(2*y**2 + z**2)*SorHd(8) + r*sqrt3*(x - y)*(x + y)*SorHd(9)))/(2.*r)
	CASE (2002)
	! (2,-2,0,0)
	 grad(3) = (sqrt3*x*y*z*(-2*SorH(10) + r*SorHd(10)))/r
	CASE (2011)
	! (2,-2,1,-1)
	 grad(3) = (x*z*((2 - 3*x**2 + 3*y**2 - 3*z**2)*SorH(11) + r*(x**2 - y**2 + &
z**2)*SorHd(11) + sqrt3*y**2*(-3*SorH(12) + r*SorHd(12))))/r
	CASE (2012)
	! (2,-2,1,0)
	 grad(3) = (x*y*((-2 + 6*z**2)*SorH(11) + sqrt3*(1 - 3*z**2)*SorH(12) + &
r*z**2*(-2*SorHd(11) + sqrt3*SorHd(12))))/r
	CASE (2013)
	! (2,-2,1,1)
	 grad(3) = (y*z*((2 + 3*x**2 - 3*y**2 - 3*z**2)*SorH(11) + r*(-x**2 + y**2 + &
z**2)*SorHd(11) + sqrt3*x**2*(-3*SorH(12) + r*SorHd(12))))/r
	CASE (2020)
	! (2,-2,2,-2)
	 grad(3) = (z*(2*(x**2 + y**2 - 2*x**2*y**2 - 2*(-1 + x**2 + y**2)*z**2 - &
2*z**4)*SorH(13) + 2*(x**2 - 2*x**4 + y**2 + 4*x**2*y**2 - 2*y**4 - 2*(x**2 + &
y**2)*z**2)*SorH(14) - 12*x**2*y**2*SorH(15) + r*x**2*y**2*SorHd(13) + &
r*x**2*z**2*SorHd(13) + r*y**2*z**2*SorHd(13) + r*z**4*SorHd(13) + &
r*x**4*SorHd(14) - 2*r*x**2*y**2*SorHd(14) + r*y**4*SorHd(14) + &
r*x**2*z**2*SorHd(14) + r*y**2*z**2*SorHd(14) + 3*r*x**2*y**2*SorHd(15)))/r
	CASE (2021)
	! (2,-2,2,-1)
	 grad(3) = (x*((-x**2 + (-3 + 4*x**2)*z**2 + 4*z**4)*SorH(13) + (x**2 - 3*y**2 + (3 - &
4*x**2 + 12*y**2)*z**2 - 4*z**4)*SorH(14) + 3*y**2*SorH(15) - &
12*y**2*z**2*SorH(15) - r*x**2*z**2*SorHd(13) - r*z**4*SorHd(13) + &
r*x**2*z**2*SorHd(14) - 3*r*y**2*z**2*SorHd(14) + r*z**4*SorHd(14) + &
3*r*y**2*z**2*SorHd(15)))/r
	CASE (2022)
	! (2,-2,2,0)
	 grad(3) = (sqrt3*x*y*z*(-4*(-1 + x**2 + y**2 + 2*z**2)*SorH(13) - 8*SorH(14) + &
4*SorH(15) + (x**2 + y**2)*(4*SorH(15) + r*SorHd(13)) + 2*z**2*(8*SorH(14) - &
4*SorH(15) + r*(SorHd(13) - 2*SorHd(14))) - r*(x**2 + y**2 - &
2*z**2)*SorHd(15)))/(2.*r)
	CASE (2023)
	! (2,-2,2,1)
	 grad(3) = (y*((-y**2 + (-3 + 4*y**2)*z**2 + 4*z**4)*SorH(13) + (-3*x**2 + y**2 + (3 + &
12*x**2 - 4*y**2)*z**2 - 4*z**4)*SorH(14) + 3*x**2*SorH(15) - &
12*x**2*z**2*SorH(15) - r*y**2*z**2*SorHd(13) - r*z**4*SorHd(13) - &
3*r*x**2*z**2*SorHd(14) + r*y**2*z**2*SorHd(14) + r*z**4*SorHd(14) + &
3*r*x**2*z**2*SorHd(15)))/r
	CASE (2024)
	! (2,-2,2,2)
	 grad(3) = (x*(x - y)*y*(x + y)*z*(-4*SorH(13) + 16*SorH(14) - 12*SorH(15) + r*(SorHd(13) &
- 4*SorHd(14) + 3*SorHd(15))))/(2.*r)
	CASE (2102)
	! (2,-1,0,0)
	 grad(3) = (sqrt3*y*((1 - 2*z**2)*SorH(10) + r*z**2*SorHd(10)))/r
	CASE (2111)
	! (2,-1,1,-1)
	 grad(3) = ((x**2 - y**2 + 3*(1 - x**2 + y**2)*z**2 - 3*z**4)*SorH(11) + r*z**2*(x**2 - &
y**2 + z**2)*SorHd(11) + sqrt3*y**2*((1 - 3*z**2)*SorH(12) + &
r*z**2*SorHd(12)))/r
	CASE (2112)
	! (2,-1,1,0)
	 grad(3) = (y*z*(-((2 + 3*x**2 + 3*y**2 - 3*z**2)*SorH(11)) + sqrt3*(2 - 3*z**2)*SorH(12) &
+ r*(x**2 + y**2 - z**2)*SorHd(11) + r*sqrt3*z**2*SorHd(12)))/r
	CASE (2113)
	! (2,-1,1,1)
	 grad(3) = (x*y*((-2 + 6*z**2)*SorH(11) + sqrt3*(1 - 3*z**2)*SorH(12) + &
r*z**2*(-2*SorHd(11) + sqrt3*SorHd(12))))/r
	CASE (2120)
	! (2,-1,2,-2)
	 grad(3) = (x*((-x**2 + (-3 + 4*x**2)*z**2 + 4*z**4)*SorH(13) + (x**2 - 3*y**2 + (3 - &
4*x**2 + 12*y**2)*z**2 - 4*z**4)*SorH(14) + 3*y**2*SorH(15) - &
12*y**2*z**2*SorH(15) - r*x**2*z**2*SorHd(13) - r*z**4*SorHd(13) + &
r*x**2*z**2*SorHd(14) - 3*r*y**2*z**2*SorHd(14) + r*z**4*SorHd(14) + &
3*r*y**2*z**2*SorHd(15)))/r
	CASE (2121)
	! (2,-1,2,-1)
	 grad(3) = (z*(-2*(x**2 + y**2)*(-1 + 2*x**2 + 2*z**2)*SorH(13) + 2*(x**2*(1 - 2*y**2 - &
2*z**2) - 2*(y - z)*(y + z)*(1 + y**2 - z**2))*SorH(14) + 6*y**2*SorH(15) - &
12*y**2*z**2*SorH(15) + r*x**4*SorHd(13) + r*x**2*y**2*SorHd(13) + &
r*x**2*z**2*SorHd(13) + r*y**2*z**2*SorHd(13) + r*x**2*y**2*SorHd(14) + &
r*y**4*SorHd(14) + r*x**2*z**2*SorHd(14) - 2*r*y**2*z**2*SorHd(14) + &
r*z**4*SorHd(14) + 3*r*y**2*z**2*SorHd(15)))/r
	CASE (2122)
	! (2,-1,2,0)
	 grad(3) = (sqrt3*y*((x**2 + y**2)*(-1 + 4*z**2)*SorH(13) + 2*(x**2 + y**2 - (3 + 4*x**2 &
+ 4*y**2)*z**2 + 4*z**4)*SorH(14) - x**2*SorH(15) - y**2*SorH(15) + &
6*z**2*SorH(15) + 4*x**2*z**2*SorH(15) + 4*y**2*z**2*SorH(15) - &
8*z**4*SorH(15) - r*x**2*z**2*SorHd(13) - r*y**2*z**2*SorHd(13) + &
2*r*x**2*z**2*SorHd(14) + 2*r*y**2*z**2*SorHd(14) - 2*r*z**4*SorHd(14) - &
r*z**2*(x**2 + y**2 - 2*z**2)*SorHd(15)))/(2.*r)
	CASE (2123)
	! (2,-1,2,1)
	 grad(3) = (x*y*z*(4*(x**2 + y**2)*SorH(13) - 2*(3 + 2*x**2 + 2*y**2 - 6*z**2)*SorH(14) + &
6*SorH(15) + r*(x**2 + y**2)*(-SorHd(13) + SorHd(14)) - 3*z**2*(4*SorH(15) + &
r*(SorHd(14) - SorHd(15)))))/r
	CASE (2124)
	! (2,-1,2,2)
	 grad(3) = (y*((3*x**2 + y**2 - 2*(-3 + 6*x**2 + 2*y**2)*z**2 - 8*z**4)*SorH(13) + &
2*(-3*x**2 + y**2 + (-3 + 12*x**2 - 4*y**2)*z**2 + 4*z**4)*SorH(14) + &
3*x**2*SorH(15) - 3*y**2*SorH(15) - 12*x**2*z**2*SorH(15) + &
12*y**2*z**2*SorH(15) + 3*r*x**2*z**2*SorHd(13) + r*y**2*z**2*SorHd(13) + &
2*r*z**4*SorHd(13) - 6*r*x**2*z**2*SorHd(14) + 2*r*y**2*z**2*SorHd(14) - &
2*r*z**4*SorHd(14) + 3*r*(x - y)*(x + y)*z**2*SorHd(15)))/(2.*r)
	CASE (2202)
	! (2,0,0,0)
	 grad(3) = (z*(2*(2 + x**2 + y**2 - 2*z**2)*SorH(10) - r*(x**2 + y**2 - &
2*z**2)*SorHd(10)))/(2.*r)
	CASE (2211)
	! (2,0,1,-1)
	 grad(3) = (y*z*(2*sqrt3*(-2 + 3*z**2)*SorH(11) + (4 + 3*x**2 + 3*y**2 - 6*z**2)*SorH(12) &
- r*(2*sqrt3*z**2*SorHd(11) + (x**2 + y**2 - 2*z**2)*SorHd(12))))/(2.*r)
	CASE (2212)
	! (2,0,1,0)
	 grad(3) = (-2*sqrt3*(x**2 + y**2)*(-1 + 3*z**2)*SorH(11) + (-x**2 - y**2 + 3*(2 + x**2 + &
y**2)*z**2 - 6*z**4)*SorH(12) + r*(x**2 + y**2)*z**2*(2*sqrt3*SorHd(11) - &
SorHd(12)) + 2*r*z**4*SorHd(12))/(2.*r)
	CASE (2213)
	! (2,0,1,1)
	 grad(3) = (x*z*(2*sqrt3*(-2 + 3*z**2)*SorH(11) + (4 + 3*x**2 + 3*y**2 - 6*z**2)*SorH(12) &
- r*(2*sqrt3*z**2*SorHd(11) + (x**2 + y**2 - 2*z**2)*SorHd(12))))/(2.*r)
	CASE (2220)
	! (2,0,2,-2)
	 grad(3) = (sqrt3*x*y*z*(-4*(-1 + x**2 + y**2 + 2*z**2)*SorH(13) - 8*SorH(14) + &
4*SorH(15) + (x**2 + y**2)*(4*SorH(15) + r*SorHd(13)) + 2*z**2*(8*SorH(14) - &
4*SorH(15) + r*(SorHd(13) - 2*SorHd(14))) - r*(x**2 + y**2 - &
2*z**2)*SorHd(15)))/(2.*r)
	CASE (2221)
	! (2,0,2,-1)
	 grad(3) = (sqrt3*y*((x**2 + y**2)*(-1 + 4*z**2)*SorH(13) + 2*(x**2 + y**2 - (3 + 4*x**2 &
+ 4*y**2)*z**2 + 4*z**4)*SorH(14) - x**2*SorH(15) - y**2*SorH(15) + &
6*z**2*SorH(15) + 4*x**2*z**2*SorH(15) + 4*y**2*z**2*SorH(15) - &
8*z**4*SorH(15) - r*x**2*z**2*SorHd(13) - r*y**2*z**2*SorHd(13) + &
2*r*x**2*z**2*SorHd(14) + 2*r*y**2*z**2*SorHd(14) - 2*r*z**4*SorHd(14) - &
r*z**2*(x**2 + y**2 - 2*z**2)*SorHd(15)))/(2.*r)
	CASE (2222)
	! (2,0,2,0)
	 grad(3) = (z*(-8*(x**2 + y**2 - 2*z**2)*SorH(15) - 4*(3*(x**2 + y**2)*((x**2 + &
y**2)*SorH(13) + 4*z**2*SorH(14)) + (x**2 + y**2 - 2*z**2)**2*SorH(15)) + &
3*(x**2 + y**2)*(8*SorH(14) + r*(x**2 + y**2)*SorHd(13) + 4*r*z**2*SorHd(14)) &
+ r*(x**2 + y**2 - 2*z**2)**2*SorHd(15)))/(4.*r)
	CASE (2223)
	! (2,0,2,1)
	 grad(3) = (sqrt3*x*((x**2 + y**2)*(-1 + 4*z**2)*SorH(13) + 2*(x**2 + y**2 - (3 + 4*x**2 &
+ 4*y**2)*z**2 + 4*z**4)*SorH(14) - x**2*SorH(15) - y**2*SorH(15) + &
6*z**2*SorH(15) + 4*x**2*z**2*SorH(15) + 4*y**2*z**2*SorH(15) - &
8*z**4*SorH(15) - r*x**2*z**2*SorHd(13) - r*y**2*z**2*SorHd(13) + &
2*r*x**2*z**2*SorHd(14) + 2*r*y**2*z**2*SorHd(14) - 2*r*z**4*SorHd(14) - &
r*z**2*(x**2 + y**2 - 2*z**2)*SorHd(15)))/(2.*r)
	CASE (2224)
	! (2,0,2,2)
	 grad(3) = (sqrt3*(x - y)*(x + y)*z*(-4*(-1 + x**2 + y**2 + 2*z**2)*SorH(13) - 8*SorH(14) &
+ 4*SorH(15) + (x**2 + y**2)*(4*SorH(15) + r*SorHd(13)) + 2*z**2*(8*SorH(14) - &
4*SorH(15) + r*(SorHd(13) - 2*SorHd(14))) - r*(x**2 + y**2 - &
2*z**2)*SorHd(15)))/(4.*r)
	CASE (2302)
	! (2,1,0,0)
	 grad(3) = (sqrt3*x*((1 - 2*z**2)*SorH(10) + r*z**2*SorHd(10)))/r
	CASE (2311)
	! (2,1,1,-1)
	 grad(3) = (x*y*((-2 + 6*z**2)*SorH(11) + sqrt3*(1 - 3*z**2)*SorH(12) + &
r*z**2*(-2*SorHd(11) + sqrt3*SorHd(12))))/r
	CASE (2312)
	! (2,1,1,0)
	 grad(3) = (x*z*(-((2 + 3*x**2 + 3*y**2 - 3*z**2)*SorH(11)) + sqrt3*(2 - 3*z**2)*SorH(12) &
+ r*(x**2 + y**2 - z**2)*SorHd(11) + r*sqrt3*z**2*SorHd(12)))/r
	CASE (2313)
	! (2,1,1,1)
	 grad(3) = ((-x**2 + y**2 + 3*(1 + x**2 - y**2)*z**2 - 3*z**4)*SorH(11) + r*z**2*(-x**2 + &
y**2 + z**2)*SorHd(11) + sqrt3*x**2*((1 - 3*z**2)*SorH(12) + &
r*z**2*SorHd(12)))/r
	CASE (2320)
	! (2,1,2,-2)
	 grad(3) = (y*((-y**2 + (-3 + 4*y**2)*z**2 + 4*z**4)*SorH(13) + (-3*x**2 + y**2 + (3 + &
12*x**2 - 4*y**2)*z**2 - 4*z**4)*SorH(14) + 3*x**2*SorH(15) - &
12*x**2*z**2*SorH(15) - r*y**2*z**2*SorHd(13) - r*z**4*SorHd(13) - &
3*r*x**2*z**2*SorHd(14) + r*y**2*z**2*SorHd(14) + r*z**4*SorHd(14) + &
3*r*x**2*z**2*SorHd(15)))/r
	CASE (2321)
	! (2,1,2,-1)
	 grad(3) = (x*y*z*(4*(x**2 + y**2)*SorH(13) - 2*(3 + 2*x**2 + 2*y**2 - 6*z**2)*SorH(14) + &
6*SorH(15) + r*(x**2 + y**2)*(-SorHd(13) + SorHd(14)) - 3*z**2*(4*SorH(15) + &
r*(SorHd(14) - SorHd(15)))))/r
	CASE (2322)
	! (2,1,2,0)
	 grad(3) = (sqrt3*x*((x**2 + y**2)*(-1 + 4*z**2)*SorH(13) + 2*(x**2 + y**2 - (3 + 4*x**2 &
+ 4*y**2)*z**2 + 4*z**4)*SorH(14) - x**2*SorH(15) - y**2*SorH(15) + &
6*z**2*SorH(15) + 4*x**2*z**2*SorH(15) + 4*y**2*z**2*SorH(15) - &
8*z**4*SorH(15) - r*x**2*z**2*SorHd(13) - r*y**2*z**2*SorHd(13) + &
2*r*x**2*z**2*SorHd(14) + 2*r*y**2*z**2*SorHd(14) - 2*r*z**4*SorHd(14) - &
r*z**2*(x**2 + y**2 - 2*z**2)*SorHd(15)))/(2.*r)
	CASE (2323)
	! (2,1,2,1)
	 grad(3) = (z*(-2*(x**2 + y**2)*(-1 + 2*y**2 + 2*z**2)*SorH(13) + 2*(y**2 + 2*z**2 - &
2*(x**4 + x**2*(1 + y**2 - 2*z**2) + z**2*(y**2 + z**2)))*SorH(14) + &
6*x**2*SorH(15) - 12*x**2*z**2*SorH(15) + r*x**2*y**2*SorHd(13) + &
r*y**4*SorHd(13) + r*x**2*z**2*SorHd(13) + r*y**2*z**2*SorHd(13) + &
r*x**4*SorHd(14) + r*x**2*y**2*SorHd(14) - 2*r*x**2*z**2*SorHd(14) + &
r*y**2*z**2*SorHd(14) + r*z**4*SorHd(14) + 3*r*x**2*z**2*SorHd(15)))/r
	CASE (2324)
	! (2,1,2,2)
	 grad(3) = -(x*((x**2 + 3*y**2 - 2*(-3 + 2*x**2 + 6*y**2)*z**2 - 8*z**4)*SorH(13) + &
2*(x**2 - 3*y**2 + (-3 - 4*x**2 + 12*y**2)*z**2 + 4*z**4)*SorH(14) - &
3*x**2*SorH(15) + 3*y**2*SorH(15) + 12*x**2*z**2*SorH(15) - &
12*y**2*z**2*SorH(15) + r*x**2*z**2*SorHd(13) + 3*r*y**2*z**2*SorHd(13) + &
2*r*z**4*SorHd(13) + 2*r*x**2*z**2*SorHd(14) - 6*r*y**2*z**2*SorHd(14) - &
2*r*z**4*SorHd(14) - 3*r*(x - y)*(x + y)*z**2*SorHd(15)))/(2.*r)
	CASE (2402)
	! (2,2,0,0)
	 grad(3) = (sqrt3*(x - y)*(x + y)*z*(-2*SorH(10) + r*SorHd(10)))/(2.*r)
	CASE (2411)
	! (2,2,1,-1)
	 grad(3) = (y*z*(2*(-2 + 6*x**2 + 3*z**2)*SorH(11) + 3*sqrt3*(-x**2 + y**2)*SorH(12) - &
2*r*(2*x**2 + z**2)*SorHd(11) + r*sqrt3*(x - y)*(x + y)*SorHd(12)))/(2.*r)
	CASE (2412)
	! (2,2,1,0)
	 grad(3) = ((x - y)*(x + y)*((-2 + 6*z**2)*SorH(11) + sqrt3*(1 - 3*z**2)*SorH(12) + &
r*z**2*(-2*SorHd(11) + sqrt3*SorHd(12))))/(2.*r)
	CASE (2413)
	! (2,2,1,1)
	 grad(3) = (x*z*(-2*(-2 + 6*y**2 + 3*z**2)*SorH(11) + 3*sqrt3*(-x**2 + y**2)*SorH(12) + &
2*r*(2*y**2 + z**2)*SorHd(11) + r*sqrt3*(x - y)*(x + y)*SorHd(12)))/(2.*r)
	CASE (2420)
	! (2,2,2,-2)
	 grad(3) = (x*(x - y)*y*(x + y)*z*(-4*SorH(13) + 16*SorH(14) - 12*SorH(15) + r*(SorHd(13) &
- 4*SorHd(14) + 3*SorHd(15))))/(2.*r)
	CASE (2421)
	! (2,2,2,-1)
	 grad(3) = (y*((3*x**2 + y**2 - 2*(-3 + 6*x**2 + 2*y**2)*z**2 - 8*z**4)*SorH(13) + &
2*(-3*x**2 + y**2 + (-3 + 12*x**2 - 4*y**2)*z**2 + 4*z**4)*SorH(14) + &
3*x**2*SorH(15) - 3*y**2*SorH(15) - 12*x**2*z**2*SorH(15) + &
12*y**2*z**2*SorH(15) + 3*r*x**2*z**2*SorHd(13) + r*y**2*z**2*SorHd(13) + &
2*r*z**4*SorHd(13) - 6*r*x**2*z**2*SorHd(14) + 2*r*y**2*z**2*SorHd(14) - &
2*r*z**4*SorHd(14) + 3*r*(x - y)*(x + y)*z**2*SorHd(15)))/(2.*r)
	CASE (2422)
	! (2,2,2,0)
	 grad(3) = (sqrt3*(x - y)*(x + y)*z*(-4*(-1 + x**2 + y**2 + 2*z**2)*SorH(13) - 8*SorH(14) &
+ 4*SorH(15) + (x**2 + y**2)*(4*SorH(15) + r*SorHd(13)) + 2*z**2*(8*SorH(14) - &
4*SorH(15) + r*(SorHd(13) - 2*SorHd(14))) - r*(x**2 + y**2 - &
2*z**2)*SorHd(15)))/(4.*r)
	CASE (2423)
	! (2,2,2,1)
	 grad(3) = -(x*((x**2 + 3*y**2 - 2*(-3 + 2*x**2 + 6*y**2)*z**2 - 8*z**4)*SorH(13) + &
2*(x**2 - 3*y**2 + (-3 - 4*x**2 + 12*y**2)*z**2 + 4*z**4)*SorH(14) - &
3*x**2*SorH(15) + 3*y**2*SorH(15) + 12*x**2*z**2*SorH(15) - &
12*y**2*z**2*SorH(15) + r*x**2*z**2*SorHd(13) + 3*r*y**2*z**2*SorHd(13) + &
2*r*z**4*SorHd(13) + 2*r*x**2*z**2*SorHd(14) - 6*r*y**2*z**2*SorHd(14) - &
2*r*z**4*SorHd(14) - 3*r*(x - y)*(x + y)*z**2*SorHd(15)))/(2.*r)
	CASE (2424)
	! (2,2,2,2)
	 grad(3) = (z*(-4*(x**4 - 2*x**2*(1 + y**2 - 2*z**2) + (-2 + y**2 + 2*z**2)*(y**2 + &
2*z**2))*SorH(13) - 8*(y**2*(-1 + 2*z**2) + x**2*(-1 + 8*y**2 + &
2*z**2))*SorH(14) - 12*x**4*SorH(15) + 24*x**2*y**2*SorH(15) - &
12*y**4*SorH(15) + r*x**4*SorHd(13) - 2*r*x**2*y**2*SorHd(13) + &
r*y**4*SorHd(13) + 4*r*x**2*z**2*SorHd(13) + 4*r*y**2*z**2*SorHd(13) + &
4*r*z**4*SorHd(13) + 16*r*x**2*y**2*SorHd(14) + 4*r*x**2*z**2*SorHd(14) + &
4*r*y**2*z**2*SorHd(14) + 3*r*(x**2 - y**2)**2*SorHd(15)))/(4.*r)
  CASE DEFAULT
	 WRITE(*,*) "BUG: No case for i=",i," which encodes (l1,m1,l2,m2)=(",l1,",",m1,",",l2,",",m2,") !"
	 STOP
   END SELECT
!$OMP END SECTIONS
!$OMP END PARALLEL
   slako_transformation_gradient(:) = grad(:)
END FUNCTION

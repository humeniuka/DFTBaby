! This file has been generated automatically


! transformation rules for matrix elements
DOUBLE PRECISION FUNCTION slako_transformation(r,x,y,z, SorH, N, l1,m1,l2,m2)
   IMPLICIT NONE
   ! x,y,z are directional cosines, r is the distance between the two centers
   DOUBLE PRECISION, INTENT(IN) :: r,x,y,z
   INTEGER, INTENT(IN) :: N  ! length of array SorH
   ! values of the N Slater-Koster tables for S or H0 evaluated at distance r
   DOUBLE PRECISION, INTENT(IN) :: SorH(N)
   ! orbital qm numbers for center 1 and center 2
   INTEGER, INTENT(IN) :: l1,m1,l2,m2
   ! Local Variables
   INTEGER, PARAMETER :: lmax=2
   DOUBLE PRECISION, PARAMETER :: sqrt3=1.7320508075688772
   ! Result S(x,y,z) or H(x,y,z) after applying SK rules
   DOUBLE PRECISION :: res
   INTEGER :: i  ! index that encodes the tuple (l1,m1,l2,m2)

   ! First we need to transform the tuple (l1,m1,l2,m2) into a unique integer
   ! so that the compiler can build a branching table for each case.
   ! Valid ranges for qm numbers: 0 <= l1,l2 <= lmax, -lmax <= m1,m2 <= lmax
   i = l1*1000+(lmax+m1)*100+l2*10+(lmax+m2)
   SELECT CASE(i)
	CASE (202)
	! (0,0,0,0)
	 res = SorH(1)
	CASE (211)
	! (0,0,1,-1)
	 res = y*SorH(3)
	CASE (212)
	! (0,0,1,0)
	 res = z*SorH(3)
	CASE (213)
	! (0,0,1,1)
	 res = x*SorH(3)
	CASE (220)
	! (0,0,2,-2)
	 res = x*y*SorH(4)*sqrt3
	CASE (221)
	! (0,0,2,-1)
	 res = y*z*SorH(4)*sqrt3
	CASE (222)
	! (0,0,2,0)
	 res = -((x**2 + y**2 - 2*z**2)*SorH(4))/2.
	CASE (223)
	! (0,0,2,1)
	 res = x*z*SorH(4)*sqrt3
	CASE (224)
	! (0,0,2,2)
	 res = ((x - y)*(x + y)*SorH(4)*sqrt3)/2.
	CASE (1102)
	! (1,-1,0,0)
	 res = y*SorH(5)
	CASE (1111)
	! (1,-1,1,-1)
	 res = (x**2 + z**2)*SorH(6) + y**2*SorH(7)
	CASE (1112)
	! (1,-1,1,0)
	 res = y*z*(-SorH(6) + SorH(7))
	CASE (1113)
	! (1,-1,1,1)
	 res = x*y*(-SorH(6) + SorH(7))
	CASE (1120)
	! (1,-1,2,-2)
	 res = x*((x**2 - y**2 + z**2)*SorH(8) + y**2*SorH(9)*sqrt3)
	CASE (1121)
	! (1,-1,2,-1)
	 res = z*((x**2 - y**2 + z**2)*SorH(8) + y**2*SorH(9)*sqrt3)
	CASE (1122)
	! (1,-1,2,0)
	 res = -(y*((x**2 + y**2 - 2*z**2)*SorH(9) + 2*z**2*SorH(8)*sqrt3))/2.
	CASE (1123)
	! (1,-1,2,1)
	 res = x*y*z*(-2*SorH(8) + SorH(9)*sqrt3)
	CASE (1124)
	! (1,-1,2,2)
	 res = -(y*(2*(2*x**2 + z**2)*SorH(8) + (-x**2 + y**2)*SorH(9)*sqrt3))/2.
	CASE (1202)
	! (1,0,0,0)
	 res = z*SorH(5)
	CASE (1211)
	! (1,0,1,-1)
	 res = y*z*(-SorH(6) + SorH(7))
	CASE (1212)
	! (1,0,1,0)
	 res = (x**2 + y**2)*SorH(6) + z**2*SorH(7)
	CASE (1213)
	! (1,0,1,1)
	 res = x*z*(-SorH(6) + SorH(7))
	CASE (1220)
	! (1,0,2,-2)
	 res = x*y*z*(-2*SorH(8) + SorH(9)*sqrt3)
	CASE (1221)
	! (1,0,2,-1)
	 res = y*((x**2 + y**2 - z**2)*SorH(8) + z**2*SorH(9)*sqrt3)
	CASE (1222)
	! (1,0,2,0)
	 res = z**3*SorH(9) - ((x**2 + y**2)*z*(SorH(9) - 2*SorH(8)*sqrt3))/2.
	CASE (1223)
	! (1,0,2,1)
	 res = x*((x**2 + y**2 - z**2)*SorH(8) + z**2*SorH(9)*sqrt3)
	CASE (1224)
	! (1,0,2,2)
	 res = -((x - y)*(x + y)*z*(2*SorH(8) - SorH(9)*sqrt3))/2.
	CASE (1302)
	! (1,1,0,0)
	 res = x*SorH(5)
	CASE (1311)
	! (1,1,1,-1)
	 res = x*y*(-SorH(6) + SorH(7))
	CASE (1312)
	! (1,1,1,0)
	 res = x*z*(-SorH(6) + SorH(7))
	CASE (1313)
	! (1,1,1,1)
	 res = (y**2 + z**2)*SorH(6) + x**2*SorH(7)
	CASE (1320)
	! (1,1,2,-2)
	 res = y*((-x**2 + y**2 + z**2)*SorH(8) + x**2*SorH(9)*sqrt3)
	CASE (1321)
	! (1,1,2,-1)
	 res = x*y*z*(-2*SorH(8) + SorH(9)*sqrt3)
	CASE (1322)
	! (1,1,2,0)
	 res = -(x*((x**2 + y**2 - 2*z**2)*SorH(9) + 2*z**2*SorH(8)*sqrt3))/2.
	CASE (1323)
	! (1,1,2,1)
	 res = z*((-x**2 + y**2 + z**2)*SorH(8) + x**2*SorH(9)*sqrt3)
	CASE (1324)
	! (1,1,2,2)
	 res = x*(2*y**2 + z**2)*SorH(8) + (x*(x - y)*(x + y)*SorH(9)*sqrt3)/2.
	CASE (2002)
	! (2,-2,0,0)
	 res = x*y*SorH(10)*sqrt3
	CASE (2011)
	! (2,-2,1,-1)
	 res = x*((x**2 - y**2 + z**2)*SorH(11) + y**2*SorH(12)*sqrt3)
	CASE (2012)
	! (2,-2,1,0)
	 res = x*y*z*(-2*SorH(11) + SorH(12)*sqrt3)
	CASE (2013)
	! (2,-2,1,1)
	 res = y*((-x**2 + y**2 + z**2)*SorH(11) + x**2*SorH(12)*sqrt3)
	CASE (2020)
	! (2,-2,2,-2)
	 res = (x**2 + z**2)*(y**2 + z**2)*SorH(13) + ((x**2 - y**2)**2 + (x**2 + y**2)*z**2)*SorH(14) + 3*x**2*y**2*SorH(15)
	CASE (2021)
	! (2,-2,2,-1)
	 res = x*z*(-((x**2 + z**2)*SorH(13)) + (x**2 - 3*y**2 + z**2)*SorH(14) + 3*y**2*SorH(15))
	CASE (2022)
	! (2,-2,2,0)
	 res = (x*y*((x**2 + y**2 + 2*z**2)*SorH(13) - 4*z**2*SorH(14) - (x**2 + y**2 - 2*z**2)*SorH(15))*sqrt3)/2.
	CASE (2023)
	! (2,-2,2,1)
	 res = y*z*(-((y**2 + z**2)*(SorH(13) - SorH(14))) + 3*x**2*(-SorH(14) + SorH(15)))
	CASE (2024)
	! (2,-2,2,2)
	 res = (x*(x - y)*y*(x + y)*(SorH(13) - 4*SorH(14) + 3*SorH(15)))/2.
	CASE (2102)
	! (2,-1,0,0)
	 res = y*z*SorH(10)*sqrt3
	CASE (2111)
	! (2,-1,1,-1)
	 res = z*((x**2 - y**2 + z**2)*SorH(11) + y**2*SorH(12)*sqrt3)
	CASE (2112)
	! (2,-1,1,0)
	 res = y*((x**2 + y**2 - z**2)*SorH(11) + z**2*SorH(12)*sqrt3)
	CASE (2113)
	! (2,-1,1,1)
	 res = x*y*z*(-2*SorH(11) + SorH(12)*sqrt3)
	CASE (2120)
	! (2,-1,2,-2)
	 res = x*z*(-((x**2 + z**2)*SorH(13)) + (x**2 - 3*y**2 + z**2)*SorH(14) + 3*y**2*SorH(15))
	CASE (2121)
	! (2,-1,2,-1)
	 res = (x**2 + y**2)*(x**2 + z**2)*SorH(13) + ((y**2 - z**2)**2 + x**2*(y**2 + z**2))*SorH(14) + 3*y**2*z**2*SorH(15)
	CASE (2122)
	! (2,-1,2,0)
	 res = -(y*z*((x**2 + y**2)*SorH(13) - 2*(x**2 + y**2 - z**2)*SorH(14) + (x**2 + y**2 - 2*z**2)*SorH(15))*sqrt3)/2.
	CASE (2123)
	! (2,-1,2,1)
	 res = x*y*(-((x**2 + y**2)*SorH(13)) + (x**2 + y**2 - 3*z**2)*SorH(14) + 3*z**2*SorH(15))
	CASE (2124)
	! (2,-1,2,2)
	 res = (y*z*((3*x**2 + y**2 + 2*z**2)*SorH(13) - 2*(3*x**2 - y**2 + z**2)*SorH(14) + 3*(x - y)*(x + y)*SorH(15)))/2.
	CASE (2202)
	! (2,0,0,0)
	 res = -((x**2 + y**2 - 2*z**2)*SorH(10))/2.
	CASE (2211)
	! (2,0,1,-1)
	 res = -(y*((x**2 + y**2 - 2*z**2)*SorH(12) + 2*z**2*SorH(11)*sqrt3))/2.
	CASE (2212)
	! (2,0,1,0)
	 res = z**3*SorH(12) - ((x**2 + y**2)*z*(SorH(12) - 2*SorH(11)*sqrt3))/2.
	CASE (2213)
	! (2,0,1,1)
	 res = -(x*((x**2 + y**2 - 2*z**2)*SorH(12) + 2*z**2*SorH(11)*sqrt3))/2.
	CASE (2220)
	! (2,0,2,-2)
	 res = (x*y*((x**2 + y**2 + 2*z**2)*SorH(13) - 4*z**2*SorH(14) - (x**2 + y**2 - 2*z**2)*SorH(15))*sqrt3)/2.
	CASE (2221)
	! (2,0,2,-1)
	 res = -(y*z*((x**2 + y**2)*SorH(13) - 2*(x**2 + y**2 - z**2)*SorH(14) + (x**2 + y**2 - 2*z**2)*SorH(15))*sqrt3)/2.
	CASE (2222)
	! (2,0,2,0)
	 res = (3*(x**2 + y**2)**2*SorH(13) + 12*(x**2 + y**2)*z**2*SorH(14) + (x**2 + y**2 - 2*z**2)**2*SorH(15))/4.
	CASE (2223)
	! (2,0,2,1)
	 res = -(x*z*((x**2 + y**2)*SorH(13) - 2*(x**2 + y**2 - z**2)*SorH(14) + (x**2 + y**2 - 2*z**2)*SorH(15))*sqrt3)/2.
	CASE (2224)
	! (2,0,2,2)
	 res = ((x - y)*(x + y)*((x**2 + y**2 + 2*z**2)*SorH(13) - 4*z**2*SorH(14) - (x**2 + y**2 - 2*z**2)*SorH(15))*sqrt3)/4.
	CASE (2302)
	! (2,1,0,0)
	 res = x*z*SorH(10)*sqrt3
	CASE (2311)
	! (2,1,1,-1)
	 res = x*y*z*(-2*SorH(11) + SorH(12)*sqrt3)
	CASE (2312)
	! (2,1,1,0)
	 res = x*((x**2 + y**2 - z**2)*SorH(11) + z**2*SorH(12)*sqrt3)
	CASE (2313)
	! (2,1,1,1)
	 res = z*((-x**2 + y**2 + z**2)*SorH(11) + x**2*SorH(12)*sqrt3)
	CASE (2320)
	! (2,1,2,-2)
	 res = y*z*(-((y**2 + z**2)*(SorH(13) - SorH(14))) + 3*x**2*(-SorH(14) + SorH(15)))
	CASE (2321)
	! (2,1,2,-1)
	 res = x*y*(-((x**2 + y**2)*SorH(13)) + (x**2 + y**2 - 3*z**2)*SorH(14) + 3*z**2*SorH(15))
	CASE (2322)
	! (2,1,2,0)
	 res = -(x*z*((x**2 + y**2)*SorH(13) - 2*(x**2 + y**2 - z**2)*SorH(14) + (x**2 + y**2 - 2*z**2)*SorH(15))*sqrt3)/2.
	CASE (2323)
	! (2,1,2,1)
         res = (x**2 + y**2)*(y**2 + z**2)*SorH(13) + (x**4 + x**2*(y**2 - 2*z**2) &
                 + z**2*(y**2 + z**2))*SorH(14) + 3*x**2*z**2*SorH(15)
	CASE (2324)
	! (2,1,2,2)
	 res = -(x*z*((x**2 + 3*y**2 + 2*z**2)*SorH(13) + 2*(x**2 - 3*y**2 - z**2)*SorH(14) + 3*(-x**2 + y**2)*SorH(15)))/2.
	CASE (2402)
	! (2,2,0,0)
	 res = ((x - y)*(x + y)*SorH(10)*sqrt3)/2.
	CASE (2411)
	! (2,2,1,-1)
	 res = -(y*(2*(2*x**2 + z**2)*SorH(11) + (-x**2 + y**2)*SorH(12)*sqrt3))/2.
	CASE (2412)
	! (2,2,1,0)
	 res = -((x - y)*(x + y)*z*(2*SorH(11) - SorH(12)*sqrt3))/2.
	CASE (2413)
	! (2,2,1,1)
	 res = x*(2*y**2 + z**2)*SorH(11) + (x*(x - y)*(x + y)*SorH(12)*sqrt3)/2.
	CASE (2420)
	! (2,2,2,-2)
	 res = (x*(x - y)*y*(x + y)*(SorH(13) - 4*SorH(14) + 3*SorH(15)))/2.
	CASE (2421)
	! (2,2,2,-1)
	 res = (y*z*((3*x**2 + y**2 + 2*z**2)*SorH(13) - 2*(3*x**2 - y**2 + z**2)*SorH(14) + 3*(x - y)*(x + y)*SorH(15)))/2.
	CASE (2422)
	! (2,2,2,0)
	 res = ((x - y)*(x + y)*((x**2 + y**2 + 2*z**2)*SorH(13) - 4*z**2*SorH(14) - (x**2 + y**2 - 2*z**2)*SorH(15))*sqrt3)/4.
	CASE (2423)
	! (2,2,2,1)
	 res = -(x*z*((x**2 + 3*y**2 + 2*z**2)*SorH(13) + 2*(x**2 - 3*y**2 - z**2)*SorH(14) + 3*(-x**2 + y**2)*SorH(15)))/2.
	CASE (2424)
	! (2,2,2,2)
         res = (((x - y)**2 + 2*z**2)*((x + y)**2 + 2*z**2)*SorH(13) + 4*(4*x**2*y**2 + &
                 (x**2 + y**2)*z**2)*SorH(14) + 3*(x**2 - y**2)**2*SorH(15))/4.
      CASE DEFAULT
	 WRITE(*,*) "BUG: No case for i=",i," which encodes (l1,m1,l2,m2)=(",l1,",",m1,",",l2,",",m2,") !"
	 STOP
   END SELECT
   slako_transformation = res
END FUNCTION

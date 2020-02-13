lmax = 1
# This file has been generated automatically
from numpy import sin,cos,sqrt,pi
phi = [ \
	lambda th1,th2: 0, \
	lambda th1,th2: (sin(th1)*sin(th2)*sqrt(3))/4., \
	lambda th1,th2: cos(th1)/2., \
	lambda th1,th2: (cos(th1)*cos(th2)*sqrt(3))/2., \
	lambda th1,th2: (sin(th1)**2*sqrt(3))/4., \
	lambda th1,th2: (3*cos(th2)*sin(th1)**2)/4., \
	lambda th1,th2: (3*sin(2*th1)*sin(th2))/8., \
	lambda th1,th2: (cos(th1)**2*sqrt(3))/2., \
	lambda th1,th2: (3*cos(th1)**2*cos(th2))/2., \
]

angular_phi = {\
	(0,0,1,-1,1,-1) : phi[1],	\
	(0,0,1,0,0,0) : phi[2],	\
	(0,0,1,0,1,0) : phi[3],	\
	(0,0,1,1,1,1) : phi[1],	\
	(1,-1,1,-1,0,0) : phi[4],	\
	(1,-1,1,-1,1,0) : phi[5],	\
	(1,0,1,-1,1,-1) : phi[6],	\
	(1,0,1,0,0,0) : phi[7],	\
	(1,-1,1,0,1,-1) : phi[6],	\
	(1,0,1,0,1,0) : phi[8],	\
	(1,1,1,0,1,1) : phi[6],	\
	(1,1,1,1,0,0) : phi[4],	\
	(1,0,1,1,1,1) : phi[6],	\
	(1,1,1,1,1,0) : phi[5],	\
}

tau2index = {\
	(0,0,1,-1,1,-1) : 1,\
	(0,0,1,0,0,0) : 2,\
	(0,0,1,0,1,0) : 3,\
	(0,0,1,1,1,1) : 1,\
	(1,-1,1,-1,0,0) : 4,\
	(1,-1,1,-1,1,0) : 5,\
	(1,0,1,-1,1,-1) : 6,\
	(1,0,1,0,0,0) : 7,\
	(1,-1,1,0,1,-1) : 6,\
	(1,0,1,0,1,0) : 8,\
	(1,1,1,0,1,1) : 6,\
	(1,1,1,1,0,0) : 4,\
	(1,0,1,1,1,1) : 6,\
	(1,1,1,1,1,0) : 5,\
}

index2tau = {\
	1 : (0,0,1,-1,1,-1),\
	2 : (0,0,1,0,0,0),\
	3 : (0,0,1,0,1,0),\
	4 : (1,-1,1,-1,0,0),\
	5 : (1,-1,1,-1,1,0),\
	6 : (1,-1,1,0,1,-1),\
	7 : (1,0,1,0,0,0),\
	8 : (1,0,1,0,1,0),\
}

# transformation rules for matrix elements
# x,y,z are directional cosines, r is the distance between the two centers
slako_transformations = {\
	(0,0,1,-1,0,0) : lambda r,x,y,z,Dipole: y*Dipole(2,r), \
	(0,0,1,-1,1,-1) : lambda r,x,y,z,Dipole: (x**2 + z**2)*Dipole(1,r) + y**2*Dipole(3,r), \
	(0,0,1,-1,1,0) : lambda r,x,y,z,Dipole: y*z*(-Dipole(1,r) + Dipole(3,r)), \
	(0,0,1,-1,1,1) : lambda r,x,y,z,Dipole: x*y*(-Dipole(1,r) + Dipole(3,r)), \
	(0,0,1,0,0,0) : lambda r,x,y,z,Dipole: z*Dipole(2,r), \
	(0,0,1,0,1,-1) : lambda r,x,y,z,Dipole: y*z*(-Dipole(1,r) + Dipole(3,r)), \
	(0,0,1,0,1,0) : lambda r,x,y,z,Dipole: (x**2 + y**2)*Dipole(1,r) + z**2*Dipole(3,r), \
	(0,0,1,0,1,1) : lambda r,x,y,z,Dipole: x*z*(-Dipole(1,r) + Dipole(3,r)), \
	(0,0,1,1,0,0) : lambda r,x,y,z,Dipole: x*Dipole(2,r), \
	(0,0,1,1,1,-1) : lambda r,x,y,z,Dipole: x*y*(-Dipole(1,r) + Dipole(3,r)), \
	(0,0,1,1,1,0) : lambda r,x,y,z,Dipole: x*z*(-Dipole(1,r) + Dipole(3,r)), \
	(0,0,1,1,1,1) : lambda r,x,y,z,Dipole: (y**2 + z**2)*Dipole(1,r) + x**2*Dipole(3,r), \
	(1,-1,1,-1,0,0) : lambda r,x,y,z,Dipole: (x**2 + z**2)*Dipole(4,r) + y**2*Dipole(7,r), \
	(1,-1,1,-1,1,-1) : lambda r,x,y,z,Dipole: y*(x**2 + z**2)*(Dipole(5,r) + 2*Dipole(6,r)) + y**3*Dipole(8,r), \
	(1,-1,1,-1,1,0) : lambda r,x,y,z,Dipole: z*((x**2 + z**2)*Dipole(5,r) + y**2*(-2*Dipole(6,r) + Dipole(8,r))), \
	(1,-1,1,-1,1,1) : lambda r,x,y,z,Dipole: x*((x**2 + z**2)*Dipole(5,r) + y**2*(-2*Dipole(6,r) + Dipole(8,r))), \
	(1,-1,1,0,0,0) : lambda r,x,y,z,Dipole: y*z*(-Dipole(4,r) + Dipole(7,r)), \
	(1,-1,1,0,1,-1) : lambda r,x,y,z,Dipole: z*(-(y**2*Dipole(5,r)) + (x**2 - y**2 + z**2)*Dipole(6,r) + y**2*Dipole(8,r)), \
	(1,-1,1,0,1,0) : lambda r,x,y,z,Dipole: y*((x**2 + y**2)*Dipole(6,r) - z**2*(Dipole(5,r) + Dipole(6,r) - Dipole(8,r))), \
	(1,-1,1,0,1,1) : lambda r,x,y,z,Dipole: x*y*z*(-Dipole(5,r) - 2*Dipole(6,r) + Dipole(8,r)), \
	(1,-1,1,1,0,0) : lambda r,x,y,z,Dipole: x*y*(-Dipole(4,r) + Dipole(7,r)), \
	(1,-1,1,1,1,-1) : lambda r,x,y,z,Dipole: x*(-(y**2*Dipole(5,r)) + (x**2 - y**2 + z**2)*Dipole(6,r) + y**2*Dipole(8,r)), \
	(1,-1,1,1,1,0) : lambda r,x,y,z,Dipole: x*y*z*(-Dipole(5,r) - 2*Dipole(6,r) + Dipole(8,r)), \
	(1,-1,1,1,1,1) : lambda r,x,y,z,Dipole: y*(-(x**2*Dipole(5,r)) + (-x**2 + y**2 + z**2)*Dipole(6,r) + x**2*Dipole(8,r)), \
	(1,0,1,-1,0,0) : lambda r,x,y,z,Dipole: y*z*(-Dipole(4,r) + Dipole(7,r)), \
	(1,0,1,-1,1,-1) : lambda r,x,y,z,Dipole: z*(-(y**2*Dipole(5,r)) + (x**2 - y**2 + z**2)*Dipole(6,r) + y**2*Dipole(8,r)), \
	(1,0,1,-1,1,0) : lambda r,x,y,z,Dipole: y*((x**2 + y**2)*Dipole(6,r) - z**2*(Dipole(5,r) + Dipole(6,r) - Dipole(8,r))), \
	(1,0,1,-1,1,1) : lambda r,x,y,z,Dipole: x*y*z*(-Dipole(5,r) - 2*Dipole(6,r) + Dipole(8,r)), \
	(1,0,1,0,0,0) : lambda r,x,y,z,Dipole: (x**2 + y**2)*Dipole(4,r) + z**2*Dipole(7,r), \
	(1,0,1,0,1,-1) : lambda r,x,y,z,Dipole: y*((x**2 + y**2)*Dipole(5,r) + z**2*(-2*Dipole(6,r) + Dipole(8,r))), \
	(1,0,1,0,1,0) : lambda r,x,y,z,Dipole: (x**2 + y**2)*z*(Dipole(5,r) + 2*Dipole(6,r)) + z**3*Dipole(8,r), \
	(1,0,1,0,1,1) : lambda r,x,y,z,Dipole: x*((x**2 + y**2)*Dipole(5,r) + z**2*(-2*Dipole(6,r) + Dipole(8,r))), \
	(1,0,1,1,0,0) : lambda r,x,y,z,Dipole: x*z*(-Dipole(4,r) + Dipole(7,r)), \
	(1,0,1,1,1,-1) : lambda r,x,y,z,Dipole: x*y*z*(-Dipole(5,r) - 2*Dipole(6,r) + Dipole(8,r)), \
	(1,0,1,1,1,0) : lambda r,x,y,z,Dipole: x*((x**2 + y**2)*Dipole(6,r) - z**2*(Dipole(5,r) + Dipole(6,r) - Dipole(8,r))), \
	(1,0,1,1,1,1) : lambda r,x,y,z,Dipole: z*(-(x**2*Dipole(5,r)) + (-x**2 + y**2 + z**2)*Dipole(6,r) + x**2*Dipole(8,r)), \
	(1,1,1,-1,0,0) : lambda r,x,y,z,Dipole: x*y*(-Dipole(4,r) + Dipole(7,r)), \
	(1,1,1,-1,1,-1) : lambda r,x,y,z,Dipole: x*(-(y**2*Dipole(5,r)) + (x**2 - y**2 + z**2)*Dipole(6,r) + y**2*Dipole(8,r)), \
	(1,1,1,-1,1,0) : lambda r,x,y,z,Dipole: x*y*z*(-Dipole(5,r) - 2*Dipole(6,r) + Dipole(8,r)), \
	(1,1,1,-1,1,1) : lambda r,x,y,z,Dipole: y*(-(x**2*Dipole(5,r)) + (-x**2 + y**2 + z**2)*Dipole(6,r) + x**2*Dipole(8,r)), \
	(1,1,1,0,0,0) : lambda r,x,y,z,Dipole: x*z*(-Dipole(4,r) + Dipole(7,r)), \
	(1,1,1,0,1,-1) : lambda r,x,y,z,Dipole: x*y*z*(-Dipole(5,r) - 2*Dipole(6,r) + Dipole(8,r)), \
	(1,1,1,0,1,0) : lambda r,x,y,z,Dipole: x*((x**2 + y**2)*Dipole(6,r) - z**2*(Dipole(5,r) + Dipole(6,r) - Dipole(8,r))), \
	(1,1,1,0,1,1) : lambda r,x,y,z,Dipole: z*(-(x**2*Dipole(5,r)) + (-x**2 + y**2 + z**2)*Dipole(6,r) + x**2*Dipole(8,r)), \
	(1,1,1,1,0,0) : lambda r,x,y,z,Dipole: (y**2 + z**2)*Dipole(4,r) + x**2*Dipole(7,r), \
	(1,1,1,1,1,-1) : lambda r,x,y,z,Dipole: y*((y**2 + z**2)*Dipole(5,r) + x**2*(-2*Dipole(6,r) + Dipole(8,r))), \
	(1,1,1,1,1,0) : lambda r,x,y,z,Dipole: z*((y**2 + z**2)*Dipole(5,r) + x**2*(-2*Dipole(6,r) + Dipole(8,r))), \
	(1,1,1,1,1,1) : lambda r,x,y,z,Dipole: x*(y**2 + z**2)*(Dipole(5,r) + 2*Dipole(6,r)) + x**3*Dipole(8,r), \
}

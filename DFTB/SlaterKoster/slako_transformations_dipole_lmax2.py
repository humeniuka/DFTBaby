# This file has been generated automatically
lmax = 2
from numpy import sin,cos,sqrt,pi
phi = [ \
	lambda th1,th2: 0, \
	lambda th1,th2: (sin(th1)*sin(th2)*sqrt(3))/4., \
	lambda th1,th2: (sin(th1)*sin(2*th2)*sqrt(15))/8., \
	lambda th1,th2: cos(th1)/2., \
	lambda th1,th2: (cos(th1)*cos(th2)*sqrt(3))/2., \
	lambda th1,th2: (cos(th1)*(-1 + 3*cos(th2)**2)*sqrt(5))/4., \
	lambda th1,th2: (sin(th1)**2*sqrt(3))/4., \
	lambda th1,th2: (3*cos(th2)*sin(th1)**2)/4., \
	lambda th1,th2: (3*sin(2*th1)*sin(th2))/8., \
	lambda th1,th2: ((1 + 3*cos(2*th2))*sin(th1)**2*sqrt(15))/16., \
	lambda th1,th2: (-3*sin(th1)**2*sin(th2)**2*sqrt(5))/16., \
	lambda th1,th2: (3*cos(th1)*cos(th2)*sin(th1)*sin(th2)*sqrt(5))/4., \
	lambda th1,th2: (3*sin(th1)**2*sin(th2)**2*sqrt(5))/16., \
	lambda th1,th2: (cos(th1)**2*sqrt(3))/2., \
	lambda th1,th2: (3*cos(th1)**2*cos(th2))/2., \
	lambda th1,th2: (cos(th1)**2*(-1 + 3*cos(th2)**2)*sqrt(15))/4., \
	lambda th1,th2: (cos(th1)*sin(th1)**2*sqrt(15))/4., \
	lambda th1,th2: (3*sin(th1)**3*sin(th2)*sqrt(5))/16., \
	lambda th1,th2: (3*cos(th1)*cos(th2)*sin(th1)**2*sqrt(5))/4., \
	lambda th1,th2: ((1 + 3*cos(2*th1))*sin(th1)*sin(th2)*sqrt(15))/16., \
	lambda th1,th2: (-3*sin(th1)**3*sin(th2)*sqrt(5))/16., \
	lambda th1,th2: (15*cos(th2)*sin(th1)**3*sin(th2))/16., \
	lambda th1,th2: (5*cos(th1)*(1 + 3*cos(2*th2))*sin(th1)**2*sqrt(3))/16., \
	lambda th1,th2: (-15*cos(th1)*sin(th1)**2*sin(th2)**2)/16., \
	lambda th1,th2: (5*(1 + 3*cos(2*th1))*sin(th1)*sin(2*th2)*sqrt(3))/32., \
	lambda th1,th2: (15*cos(th1)*sin(th1)**2*sin(th2)**2)/16., \
	lambda th1,th2: (-15*cos(th2)*sin(th1)**3*sin(th2))/16., \
	lambda th1,th2: (cos(th1)*(-1 + 3*cos(th1)**2)*sqrt(5))/4., \
	lambda th1,th2: (3*cos(th1)**2*sin(th1)*sin(th2)*sqrt(5))/4., \
	lambda th1,th2: (cos(th1)*(-1 + 3*cos(th1)**2)*cos(th2)*sqrt(15))/4., \
	lambda th1,th2: (15*cos(th1)**2*sin(th1)*sin(2*th2))/8., \
	lambda th1,th2: (5*cos(th1)*(-1 + 3*cos(th1)**2)*(-1 + 3*cos(th2)**2))/8., \
]

angular_phi = {\
	(0,0,1,-1,1,-1) : phi[1],	\
	(0,0,1,-1,2,-1) : phi[2],	\
	(0,0,1,0,0,0) : phi[3],	\
	(0,0,1,0,1,0) : phi[4],	\
	(0,0,1,0,2,0) : phi[5],	\
	(0,0,1,1,1,1) : phi[1],	\
	(0,0,1,1,2,1) : phi[2],	\
	(1,-1,1,-1,0,0) : phi[6],	\
	(1,-1,1,-1,1,0) : phi[7],	\
	(1,0,1,-1,1,-1) : phi[8],	\
	(1,-1,1,-1,2,0) : phi[9],	\
	(1,-1,1,-1,2,2) : phi[10],	\
	(1,0,1,-1,2,-1) : phi[11],	\
	(1,1,1,-1,2,-2) : phi[12],	\
	(1,0,1,0,0,0) : phi[13],	\
	(1,-1,1,0,1,-1) : phi[8],	\
	(1,0,1,0,1,0) : phi[14],	\
	(1,1,1,0,1,1) : phi[8],	\
	(1,-1,1,0,2,-1) : phi[11],	\
	(1,0,1,0,2,0) : phi[15],	\
	(1,1,1,0,2,1) : phi[11],	\
	(1,1,1,1,0,0) : phi[6],	\
	(1,0,1,1,1,1) : phi[8],	\
	(1,1,1,1,1,0) : phi[7],	\
	(1,-1,1,1,2,-2) : phi[12],	\
	(1,0,1,1,2,1) : phi[11],	\
	(1,1,1,1,2,0) : phi[9],	\
	(1,1,1,1,2,2) : phi[12],	\
	(2,-1,1,-1,0,0) : phi[16],	\
	(2,-2,1,-1,1,1) : phi[17],	\
	(2,-1,1,-1,1,0) : phi[18],	\
	(2,0,1,-1,1,-1) : phi[19],	\
	(2,2,1,-1,1,-1) : phi[20],	\
	(2,-2,1,-1,2,1) : phi[21],	\
	(2,-1,1,-1,2,0) : phi[22],	\
	(2,-1,1,-1,2,2) : phi[23],	\
	(2,0,1,-1,2,-1) : phi[24],	\
	(2,1,1,-1,2,-2) : phi[25],	\
	(2,2,1,-1,2,-1) : phi[26],	\
	(2,0,1,0,0,0) : phi[27],	\
	(2,-1,1,0,1,-1) : phi[28],	\
	(2,0,1,0,1,0) : phi[29],	\
	(2,1,1,0,1,1) : phi[28],	\
	(2,-2,1,0,2,-2) : phi[25],	\
	(2,-1,1,0,2,-1) : phi[30],	\
	(2,0,1,0,2,0) : phi[31],	\
	(2,1,1,0,2,1) : phi[30],	\
	(2,2,1,0,2,2) : phi[25],	\
	(2,1,1,1,0,0) : phi[16],	\
	(2,-2,1,1,1,-1) : phi[17],	\
	(2,0,1,1,1,1) : phi[19],	\
	(2,1,1,1,1,0) : phi[18],	\
	(2,2,1,1,1,1) : phi[17],	\
	(2,-2,1,1,2,-1) : phi[21],	\
	(2,-1,1,1,2,-2) : phi[25],	\
	(2,0,1,1,2,1) : phi[24],	\
	(2,1,1,1,2,0) : phi[22],	\
	(2,1,1,1,2,2) : phi[25],	\
	(2,2,1,1,2,1) : phi[21],	\
}

tau2index = {\
	(0,0,1,-1,1,-1) : 1,\
	(0,0,1,-1,2,-1) : 2,\
	(0,0,1,0,0,0) : 3,\
	(0,0,1,0,1,0) : 4,\
	(0,0,1,0,2,0) : 5,\
	(0,0,1,1,1,1) : 1,\
	(0,0,1,1,2,1) : 2,\
	(1,-1,1,-1,0,0) : 6,\
	(1,-1,1,-1,1,0) : 7,\
	(1,0,1,-1,1,-1) : 8,\
	(1,-1,1,-1,2,0) : 9,\
	(1,-1,1,-1,2,2) : 10,\
	(1,0,1,-1,2,-1) : 11,\
	(1,1,1,-1,2,-2) : 12,\
	(1,0,1,0,0,0) : 13,\
	(1,-1,1,0,1,-1) : 8,\
	(1,0,1,0,1,0) : 14,\
	(1,1,1,0,1,1) : 8,\
	(1,-1,1,0,2,-1) : 11,\
	(1,0,1,0,2,0) : 15,\
	(1,1,1,0,2,1) : 11,\
	(1,1,1,1,0,0) : 6,\
	(1,0,1,1,1,1) : 8,\
	(1,1,1,1,1,0) : 7,\
	(1,-1,1,1,2,-2) : 12,\
	(1,0,1,1,2,1) : 11,\
	(1,1,1,1,2,0) : 9,\
	(1,1,1,1,2,2) : 12,\
	(2,-1,1,-1,0,0) : 16,\
	(2,-2,1,-1,1,1) : 17,\
	(2,-1,1,-1,1,0) : 18,\
	(2,0,1,-1,1,-1) : 19,\
	(2,2,1,-1,1,-1) : 20,\
	(2,-2,1,-1,2,1) : 21,\
	(2,-1,1,-1,2,0) : 22,\
	(2,-1,1,-1,2,2) : 23,\
	(2,0,1,-1,2,-1) : 24,\
	(2,1,1,-1,2,-2) : 25,\
	(2,2,1,-1,2,-1) : 26,\
	(2,0,1,0,0,0) : 27,\
	(2,-1,1,0,1,-1) : 28,\
	(2,0,1,0,1,0) : 29,\
	(2,1,1,0,1,1) : 28,\
	(2,-2,1,0,2,-2) : 25,\
	(2,-1,1,0,2,-1) : 30,\
	(2,0,1,0,2,0) : 31,\
	(2,1,1,0,2,1) : 30,\
	(2,2,1,0,2,2) : 25,\
	(2,1,1,1,0,0) : 16,\
	(2,-2,1,1,1,-1) : 17,\
	(2,0,1,1,1,1) : 19,\
	(2,1,1,1,1,0) : 18,\
	(2,2,1,1,1,1) : 17,\
	(2,-2,1,1,2,-1) : 21,\
	(2,-1,1,1,2,-2) : 25,\
	(2,0,1,1,2,1) : 24,\
	(2,1,1,1,2,0) : 22,\
	(2,1,1,1,2,2) : 25,\
	(2,2,1,1,2,1) : 21,\
}

index2tau = {\
	1 : (0,0,1,-1,1,-1),\
	2 : (0,0,1,-1,2,-1),\
	3 : (0,0,1,0,0,0),\
	4 : (0,0,1,0,1,0),\
	5 : (0,0,1,0,2,0),\
	6 : (1,-1,1,-1,0,0),\
	7 : (1,-1,1,-1,1,0),\
	9 : (1,-1,1,-1,2,0),\
	10 : (1,-1,1,-1,2,2),\
	8 : (1,-1,1,0,1,-1),\
	11 : (1,-1,1,0,2,-1),\
	12 : (1,-1,1,1,2,-2),\
	13 : (1,0,1,0,0,0),\
	14 : (1,0,1,0,1,0),\
	15 : (1,0,1,0,2,0),\
	17 : (2,-2,1,-1,1,1),\
	21 : (2,-2,1,-1,2,1),\
	25 : (2,-2,1,0,2,-2),\
	16 : (2,-1,1,-1,0,0),\
	18 : (2,-1,1,-1,1,0),\
	22 : (2,-1,1,-1,2,0),\
	23 : (2,-1,1,-1,2,2),\
	28 : (2,-1,1,0,1,-1),\
	30 : (2,-1,1,0,2,-1),\
	19 : (2,0,1,-1,1,-1),\
	24 : (2,0,1,-1,2,-1),\
	27 : (2,0,1,0,0,0),\
	29 : (2,0,1,0,1,0),\
	31 : (2,0,1,0,2,0),\
	20 : (2,2,1,-1,1,-1),\
	26 : (2,2,1,-1,2,-1),\
}

# transformation rules for matrix elements
# x,y,z are directional cosines, r is the distance between the two centers
slako_transformations = {\
	(0,0,1,-1,0,0) : lambda r,x,y,z,Dipole: y*Dipole(3,r), \
	(0,0,1,-1,1,-1) : lambda r,x,y,z,Dipole: (x**2 + z**2)*Dipole(1,r) + y**2*Dipole(4,r), \
	(0,0,1,-1,1,0) : lambda r,x,y,z,Dipole: y*z*(-Dipole(1,r) + Dipole(4,r)), \
	(0,0,1,-1,1,1) : lambda r,x,y,z,Dipole: x*y*(-Dipole(1,r) + Dipole(4,r)), \
	(0,0,1,-1,2,-2) : lambda r,x,y,z,Dipole: x*((x**2 - y**2 + z**2)*Dipole(2,r) + y**2*Dipole(5,r)*sqrt(3)), \
	(0,0,1,-1,2,-1) : lambda r,x,y,z,Dipole: z*((x**2 - y**2 + z**2)*Dipole(2,r) + y**2*Dipole(5,r)*sqrt(3)), \
	(0,0,1,-1,2,0) : lambda r,x,y,z,Dipole: -(y*((x**2 + y**2 - 2*z**2)*Dipole(5,r) + 2*z**2*Dipole(2,r)*sqrt(3)))/2., \
	(0,0,1,-1,2,1) : lambda r,x,y,z,Dipole: x*y*z*(-2*Dipole(2,r) + Dipole(5,r)*sqrt(3)), \
	(0,0,1,-1,2,2) : lambda r,x,y,z,Dipole: -(y*(2*(2*x**2 + z**2)*Dipole(2,r) + (-x**2 + y**2)*Dipole(5,r)*sqrt(3)))/2., \
	(0,0,1,0,0,0) : lambda r,x,y,z,Dipole: z*Dipole(3,r), \
	(0,0,1,0,1,-1) : lambda r,x,y,z,Dipole: y*z*(-Dipole(1,r) + Dipole(4,r)), \
	(0,0,1,0,1,0) : lambda r,x,y,z,Dipole: (x**2 + y**2)*Dipole(1,r) + z**2*Dipole(4,r), \
	(0,0,1,0,1,1) : lambda r,x,y,z,Dipole: x*z*(-Dipole(1,r) + Dipole(4,r)), \
	(0,0,1,0,2,-2) : lambda r,x,y,z,Dipole: x*y*z*(-2*Dipole(2,r) + Dipole(5,r)*sqrt(3)), \
	(0,0,1,0,2,-1) : lambda r,x,y,z,Dipole: y*((x**2 + y**2 - z**2)*Dipole(2,r) + z**2*Dipole(5,r)*sqrt(3)), \
	(0,0,1,0,2,0) : lambda r,x,y,z,Dipole: z**3*Dipole(5,r) - ((x**2 + y**2)*z*(Dipole(5,r) - 2*Dipole(2,r)*sqrt(3)))/2., \
	(0,0,1,0,2,1) : lambda r,x,y,z,Dipole: x*((x**2 + y**2 - z**2)*Dipole(2,r) + z**2*Dipole(5,r)*sqrt(3)), \
	(0,0,1,0,2,2) : lambda r,x,y,z,Dipole: -((x - y)*(x + y)*z*(2*Dipole(2,r) - Dipole(5,r)*sqrt(3)))/2., \
	(0,0,1,1,0,0) : lambda r,x,y,z,Dipole: x*Dipole(3,r), \
	(0,0,1,1,1,-1) : lambda r,x,y,z,Dipole: x*y*(-Dipole(1,r) + Dipole(4,r)), \
	(0,0,1,1,1,0) : lambda r,x,y,z,Dipole: x*z*(-Dipole(1,r) + Dipole(4,r)), \
	(0,0,1,1,1,1) : lambda r,x,y,z,Dipole: (y**2 + z**2)*Dipole(1,r) + x**2*Dipole(4,r), \
	(0,0,1,1,2,-2) : lambda r,x,y,z,Dipole: y*((-x**2 + y**2 + z**2)*Dipole(2,r) + x**2*Dipole(5,r)*sqrt(3)), \
	(0,0,1,1,2,-1) : lambda r,x,y,z,Dipole: x*y*z*(-2*Dipole(2,r) + Dipole(5,r)*sqrt(3)), \
	(0,0,1,1,2,0) : lambda r,x,y,z,Dipole: -(x*((x**2 + y**2 - 2*z**2)*Dipole(5,r) + 2*z**2*Dipole(2,r)*sqrt(3)))/2., \
	(0,0,1,1,2,1) : lambda r,x,y,z,Dipole: z*((-x**2 + y**2 + z**2)*Dipole(2,r) + x**2*Dipole(5,r)*sqrt(3)), \
	(0,0,1,1,2,2) : lambda r,x,y,z,Dipole: x*(2*y**2 + z**2)*Dipole(2,r) + (x*(x - y)*(x + y)*Dipole(5,r)*sqrt(3))/2., \
	(1,-1,1,-1,0,0) : lambda r,x,y,z,Dipole: (x**2 + z**2)*Dipole(6,r) + y**2*Dipole(13,r), \
	(1,-1,1,-1,1,-1) : lambda r,x,y,z,Dipole: y*(x**2 + z**2)*(Dipole(7,r) + 2*Dipole(8,r)) + y**3*Dipole(14,r), \
	(1,-1,1,-1,1,0) : lambda r,x,y,z,Dipole: z*((x**2 + z**2)*Dipole(7,r) + y**2*(-2*Dipole(8,r) + Dipole(14,r))), \
	(1,-1,1,-1,1,1) : lambda r,x,y,z,Dipole: x*((x**2 + z**2)*Dipole(7,r) + y**2*(-2*Dipole(8,r) + Dipole(14,r))), \
	(1,-1,1,-1,2,-2) : lambda r,x,y,z,Dipole: (x*y*(-(y**2*z**2*(x**2 + y**2 + 2*z**2)*Dipole(10,r)) + 2*(x**2 + y**2)**2*(x**2 - y**2 + z**2)*Dipole(11,r) - (x**4 + x**2*y**2 + 2*y**2*z**2)*Dipole(12,r) + (x**2 + y**2)**2*((x**2 + z**2)*Dipole(9,r) + y**2*Dipole(15,r))*sqrt(3)))/(x**2 + y**2)**2, \
	(1,-1,1,-1,2,-1) : lambda r,x,y,z,Dipole: (y*z*(y**2*z**2*Dipole(10,r) + 2*(x**2 + y**2)*(x**2 - y**2 + z**2)*Dipole(11,r) - x**2*Dipole(12,r) + (x**2 + y**2)*((x**2 + z**2)*Dipole(9,r) + y**2*Dipole(15,r))*sqrt(3)))/(x**2 + y**2), \
	(1,-1,1,-1,2,0) : lambda r,x,y,z,Dipole: (-((x**2 + y**2 - 2*z**2)*((x**2 + z**2)*Dipole(9,r) + y**2*Dipole(15,r))) - (y**2*z**2*(Dipole(10,r) + 4*Dipole(11,r)) + x**2*(x**2 + y**2 + z**2)*Dipole(12,r))*sqrt(3))/2., \
	(1,-1,1,-1,2,1) : lambda r,x,y,z,Dipole: (x*z*(y**2*z**2*Dipole(10,r) - 4*y**2*(x**2 + y**2)*Dipole(11,r) + (x**2 + 2*y**2)*Dipole(12,r) + (x**2 + y**2)*((x**2 + z**2)*Dipole(9,r) + y**2*Dipole(15,r))*sqrt(3)))/(x**2 + y**2), \
	(1,-1,1,-1,2,2) : lambda r,x,y,z,Dipole: (y**2*(-x + y)*(x + y)*z**2*(x**2 + y**2 + 2*z**2)*Dipole(10,r) - 4*y**2*(x**2 + y**2)**2*(2*x**2 + z**2)*Dipole(11,r) - x**2*(x**4 - y**4 + 2*(x**2 + 3*y**2)*z**2)*Dipole(12,r) + (x - y)*(x + y)*(x**2 + y**2)**2*((x**2 + z**2)*Dipole(9,r) + y**2*Dipole(15,r))*sqrt(3))/(2.*(x**2 + y**2)**2), \
	(1,-1,1,0,0,0) : lambda r,x,y,z,Dipole: y*z*(-Dipole(6,r) + Dipole(13,r)), \
	(1,-1,1,0,1,-1) : lambda r,x,y,z,Dipole: z*(-(y**2*Dipole(7,r)) + (x**2 - y**2 + z**2)*Dipole(8,r) + y**2*Dipole(14,r)), \
	(1,-1,1,0,1,0) : lambda r,x,y,z,Dipole: y*((x**2 + y**2)*Dipole(8,r) - z**2*(Dipole(7,r) + Dipole(8,r) - Dipole(14,r))), \
	(1,-1,1,0,1,1) : lambda r,x,y,z,Dipole: x*y*z*(-Dipole(7,r) - 2*Dipole(8,r) + Dipole(14,r)), \
	(1,-1,1,0,2,-2) : lambda r,x,y,z,Dipole: (x*z*(y**2*(x**2 + y**2 + 2*z**2)*Dipole(10,r) + (x**2 + y**2)*(x**2 - 3*y**2 + z**2)*Dipole(11,r) + (-x**2 + y**2)*Dipole(12,r) - y**2*(x**2 + y**2)*(Dipole(9,r) - Dipole(15,r))*sqrt(3)))/(x**2 + y**2), \
	(1,-1,1,0,2,-1) : lambda r,x,y,z,Dipole: -(y**2*z**2*Dipole(10,r)) + ((y**2 - z**2)**2 + x**2*(y**2 + z**2))*Dipole(11,r) + x**2*(x**2 + y**2 + z**2)*Dipole(12,r) + y**2*z**2*(-Dipole(9,r) + Dipole(15,r))*sqrt(3), \
	(1,-1,1,0,2,0) : lambda r,x,y,z,Dipole: (y*z*((x**2 + y**2 - 2*z**2)*(Dipole(9,r) - Dipole(15,r)) + ((x**2 + y**2)*Dipole(10,r) + 2*(x**2 + y**2 - z**2)*Dipole(11,r))*sqrt(3)))/2., \
	(1,-1,1,0,2,1) : lambda r,x,y,z,Dipole: x*y*((x**2 + y**2)*(Dipole(11,r) - Dipole(12,r)) - z**2*(Dipole(10,r) + 3*Dipole(11,r) + Dipole(12,r) + (Dipole(9,r) - Dipole(15,r))*sqrt(3))), \
	(1,-1,1,0,2,2) : lambda r,x,y,z,Dipole: (y*z*((x - y)*(x + y)*(x**2 + y**2 + 2*z**2)*Dipole(10,r) - 2*(x**2 + y**2)*(3*x**2 - y**2 + z**2)*Dipole(11,r) + 4*x**2*Dipole(12,r) - (x**4 - y**4)*(Dipole(9,r) - Dipole(15,r))*sqrt(3)))/(2.*(x**2 + y**2)), \
	(1,-1,1,1,0,0) : lambda r,x,y,z,Dipole: x*y*(-Dipole(6,r) + Dipole(13,r)), \
	(1,-1,1,1,1,-1) : lambda r,x,y,z,Dipole: x*(-(y**2*Dipole(7,r)) + (x**2 - y**2 + z**2)*Dipole(8,r) + y**2*Dipole(14,r)), \
	(1,-1,1,1,1,0) : lambda r,x,y,z,Dipole: x*y*z*(-Dipole(7,r) - 2*Dipole(8,r) + Dipole(14,r)), \
	(1,-1,1,1,1,1) : lambda r,x,y,z,Dipole: y*(-(x**2*Dipole(7,r)) + (-x**2 + y**2 + z**2)*Dipole(8,r) + x**2*Dipole(14,r)), \
	(1,-1,1,1,2,-2) : lambda r,x,y,z,Dipole: (-(x**2*y**2*z**2*(x**2 + y**2 + 2*z**2)*Dipole(10,r)) + (x**2 + y**2)**2*((x**2 - y**2)**2 + (x**2 + y**2)*z**2)*Dipole(11,r) + (x**2*y**2*(x**2 + y**2) + (x**4 + y**4)*z**2)*Dipole(12,r) - y**2*(x**3 + x*y**2)**2*(Dipole(9,r) - Dipole(15,r))*sqrt(3))/(x**2 + y**2)**2, \
	(1,-1,1,1,2,-1) : lambda r,x,y,z,Dipole: (x*z*(y**2*z**2*Dipole(10,r) + (x**2 + y**2)*(x**2 - 3*y**2 + z**2)*Dipole(11,r) - x**2*Dipole(12,r) - y**2*(x**2 + y**2)*(Dipole(9,r) - Dipole(15,r))*sqrt(3)))/(x**2 + y**2), \
	(1,-1,1,1,2,0) : lambda r,x,y,z,Dipole: (x*y*((x**2 + y**2 - 2*z**2)*(Dipole(9,r) - Dipole(15,r)) + (-(z**2*(Dipole(10,r) + 4*Dipole(11,r))) + (x**2 + y**2 + z**2)*Dipole(12,r))*sqrt(3)))/2., \
	(1,-1,1,1,2,1) : lambda r,x,y,z,Dipole: (y*z*(x**2*z**2*Dipole(10,r) - y**2*Dipole(12,r) - x**2*(x**2 + y**2)*Dipole(9,r)*sqrt(3) + (x**2 + y**2)*((-3*x**2 + y**2 + z**2)*Dipole(11,r) + x**2*Dipole(15,r)*sqrt(3))))/(x**2 + y**2), \
	(1,-1,1,1,2,2) : lambda r,x,y,z,Dipole: (x*(x - y)*y*(x + y)*(-(z**2*(x**2 + y**2 + 2*z**2)*Dipole(10,r)) - 4*(x**2 + y**2)**2*Dipole(11,r) + (x**2 + y**2 - 2*z**2)*Dipole(12,r) - (x**2 + y**2)**2*(Dipole(9,r) - Dipole(15,r))*sqrt(3)))/(2.*(x**2 + y**2)**2), \
	(1,0,1,-1,0,0) : lambda r,x,y,z,Dipole: y*z*(-Dipole(6,r) + Dipole(13,r)), \
	(1,0,1,-1,1,-1) : lambda r,x,y,z,Dipole: z*(-(y**2*Dipole(7,r)) + (x**2 - y**2 + z**2)*Dipole(8,r) + y**2*Dipole(14,r)), \
	(1,0,1,-1,1,0) : lambda r,x,y,z,Dipole: y*((x**2 + y**2)*Dipole(8,r) - z**2*(Dipole(7,r) + Dipole(8,r) - Dipole(14,r))), \
	(1,0,1,-1,1,1) : lambda r,x,y,z,Dipole: x*y*z*(-Dipole(7,r) - 2*Dipole(8,r) + Dipole(14,r)), \
	(1,0,1,-1,2,-2) : lambda r,x,y,z,Dipole: (x*z*(y**2*(x**2 + y**2 + 2*z**2)*Dipole(10,r) + (x**2 + y**2)*(x**2 - 3*y**2 + z**2)*Dipole(11,r) + (-x**2 + y**2)*Dipole(12,r) - y**2*(x**2 + y**2)*(Dipole(9,r) - Dipole(15,r))*sqrt(3)))/(x**2 + y**2), \
	(1,0,1,-1,2,-1) : lambda r,x,y,z,Dipole: -(y**2*z**2*Dipole(10,r)) + ((y**2 - z**2)**2 + x**2*(y**2 + z**2))*Dipole(11,r) + x**2*(x**2 + y**2 + z**2)*Dipole(12,r) + y**2*z**2*(-Dipole(9,r) + Dipole(15,r))*sqrt(3), \
	(1,0,1,-1,2,0) : lambda r,x,y,z,Dipole: (y*z*((x**2 + y**2 - 2*z**2)*(Dipole(9,r) - Dipole(15,r)) + ((x**2 + y**2)*Dipole(10,r) + 2*(x**2 + y**2 - z**2)*Dipole(11,r))*sqrt(3)))/2., \
	(1,0,1,-1,2,1) : lambda r,x,y,z,Dipole: x*y*((x**2 + y**2)*(Dipole(11,r) - Dipole(12,r)) - z**2*(Dipole(10,r) + 3*Dipole(11,r) + Dipole(12,r) + (Dipole(9,r) - Dipole(15,r))*sqrt(3))), \
	(1,0,1,-1,2,2) : lambda r,x,y,z,Dipole: (y*z*((x - y)*(x + y)*(x**2 + y**2 + 2*z**2)*Dipole(10,r) - 2*(x**2 + y**2)*(3*x**2 - y**2 + z**2)*Dipole(11,r) + 4*x**2*Dipole(12,r) - (x**4 - y**4)*(Dipole(9,r) - Dipole(15,r))*sqrt(3)))/(2.*(x**2 + y**2)), \
	(1,0,1,0,0,0) : lambda r,x,y,z,Dipole: (x**2 + y**2)*Dipole(6,r) + z**2*Dipole(13,r), \
	(1,0,1,0,1,-1) : lambda r,x,y,z,Dipole: y*((x**2 + y**2)*Dipole(7,r) + z**2*(-2*Dipole(8,r) + Dipole(14,r))), \
	(1,0,1,0,1,0) : lambda r,x,y,z,Dipole: (x**2 + y**2)*z*(Dipole(7,r) + 2*Dipole(8,r)) + z**3*Dipole(14,r), \
	(1,0,1,0,1,1) : lambda r,x,y,z,Dipole: x*((x**2 + y**2)*Dipole(7,r) + z**2*(-2*Dipole(8,r) + Dipole(14,r))), \
	(1,0,1,0,2,-2) : lambda r,x,y,z,Dipole: x*y*(-((x**2 + y**2 + 2*z**2)*Dipole(10,r)) + (x**2 + y**2)*Dipole(9,r)*sqrt(3) + z**2*(-4*Dipole(11,r) + Dipole(15,r)*sqrt(3))), \
	(1,0,1,0,2,-1) : lambda r,x,y,z,Dipole: y*(x**2 + y**2)*z*(Dipole(10,r) + 2*Dipole(11,r) + Dipole(9,r)*sqrt(3)) + y*z**3*(-2*Dipole(11,r) + Dipole(15,r)*sqrt(3)), \
	(1,0,1,0,2,0) : lambda r,x,y,z,Dipole: (-((x**2 + y**2 - 2*z**2)*((x**2 + y**2)*Dipole(9,r) + z**2*Dipole(15,r))) - (x**2 + y**2)*((x**2 + y**2)*Dipole(10,r) - 4*z**2*Dipole(11,r))*sqrt(3))/2., \
	(1,0,1,0,2,1) : lambda r,x,y,z,Dipole: x*(x**2 + y**2)*z*(Dipole(10,r) + 2*Dipole(11,r) + Dipole(9,r)*sqrt(3)) + x*z**3*(-2*Dipole(11,r) + Dipole(15,r)*sqrt(3)), \
	(1,0,1,0,2,2) : lambda r,x,y,z,Dipole: ((x - y)*(x + y)*(-((x**2 + y**2 + 2*z**2)*Dipole(10,r)) + (x**2 + y**2)*Dipole(9,r)*sqrt(3) + z**2*(-4*Dipole(11,r) + Dipole(15,r)*sqrt(3))))/2., \
	(1,0,1,1,0,0) : lambda r,x,y,z,Dipole: x*z*(-Dipole(6,r) + Dipole(13,r)), \
	(1,0,1,1,1,-1) : lambda r,x,y,z,Dipole: x*y*z*(-Dipole(7,r) - 2*Dipole(8,r) + Dipole(14,r)), \
	(1,0,1,1,1,0) : lambda r,x,y,z,Dipole: x*((x**2 + y**2)*Dipole(8,r) - z**2*(Dipole(7,r) + Dipole(8,r) - Dipole(14,r))), \
	(1,0,1,1,1,1) : lambda r,x,y,z,Dipole: z*(-(x**2*Dipole(7,r)) + (-x**2 + y**2 + z**2)*Dipole(8,r) + x**2*Dipole(14,r)), \
	(1,0,1,1,2,-2) : lambda r,x,y,z,Dipole: (y*z*(x**2*(x**2 + y**2 + 2*z**2)*Dipole(10,r) + (x - y)*(x + y)*Dipole(12,r) - x**2*(x**2 + y**2)*Dipole(9,r)*sqrt(3) + (x**2 + y**2)*((-3*x**2 + y**2 + z**2)*Dipole(11,r) + x**2*Dipole(15,r)*sqrt(3))))/(x**2 + y**2), \
	(1,0,1,1,2,-1) : lambda r,x,y,z,Dipole: x*y*((x**2 + y**2)*(Dipole(11,r) - Dipole(12,r)) - z**2*(Dipole(10,r) + 3*Dipole(11,r) + Dipole(12,r) + (Dipole(9,r) - Dipole(15,r))*sqrt(3))), \
	(1,0,1,1,2,0) : lambda r,x,y,z,Dipole: (x*z*((x**2 + y**2 - 2*z**2)*(Dipole(9,r) - Dipole(15,r)) + ((x**2 + y**2)*Dipole(10,r) + 2*(x**2 + y**2 - z**2)*Dipole(11,r))*sqrt(3)))/2., \
	(1,0,1,1,2,1) : lambda r,x,y,z,Dipole: x**4*Dipole(11,r) + (y**2 + z**2)*(z**2*Dipole(11,r) + y**2*Dipole(12,r)) + x**2*(y**2*(Dipole(11,r) + Dipole(12,r)) - z**2*(Dipole(10,r) + 2*Dipole(11,r) + (Dipole(9,r) - Dipole(15,r))*sqrt(3))), \
	(1,0,1,1,2,2) : lambda r,x,y,z,Dipole: (x*z*((x - y)*(x + y)*(x**2 + y**2 + 2*z**2)*Dipole(10,r) - 2*(x**2 + y**2)*(x**2 - 3*y**2 - z**2)*Dipole(11,r) - 4*y**2*Dipole(12,r) + (-x**4 + y**4)*(Dipole(9,r) - Dipole(15,r))*sqrt(3)))/(2.*(x**2 + y**2)), \
	(1,1,1,-1,0,0) : lambda r,x,y,z,Dipole: x*y*(-Dipole(6,r) + Dipole(13,r)), \
	(1,1,1,-1,1,-1) : lambda r,x,y,z,Dipole: x*(-(y**2*Dipole(7,r)) + (x**2 - y**2 + z**2)*Dipole(8,r) + y**2*Dipole(14,r)), \
	(1,1,1,-1,1,0) : lambda r,x,y,z,Dipole: x*y*z*(-Dipole(7,r) - 2*Dipole(8,r) + Dipole(14,r)), \
	(1,1,1,-1,1,1) : lambda r,x,y,z,Dipole: y*(-(x**2*Dipole(7,r)) + (-x**2 + y**2 + z**2)*Dipole(8,r) + x**2*Dipole(14,r)), \
	(1,1,1,-1,2,-2) : lambda r,x,y,z,Dipole: (-(x**2*y**2*z**2*(x**2 + y**2 + 2*z**2)*Dipole(10,r)) + (x**2 + y**2)**2*((x**2 - y**2)**2 + (x**2 + y**2)*z**2)*Dipole(11,r) + (x**2*y**2*(x**2 + y**2) + (x**4 + y**4)*z**2)*Dipole(12,r) - y**2*(x**3 + x*y**2)**2*(Dipole(9,r) - Dipole(15,r))*sqrt(3))/(x**2 + y**2)**2, \
	(1,1,1,-1,2,-1) : lambda r,x,y,z,Dipole: (x*z*(y**2*z**2*Dipole(10,r) + (x**2 + y**2)*(x**2 - 3*y**2 + z**2)*Dipole(11,r) - x**2*Dipole(12,r) - y**2*(x**2 + y**2)*(Dipole(9,r) - Dipole(15,r))*sqrt(3)))/(x**2 + y**2), \
	(1,1,1,-1,2,0) : lambda r,x,y,z,Dipole: (x*y*((x**2 + y**2 - 2*z**2)*(Dipole(9,r) - Dipole(15,r)) + (-(z**2*(Dipole(10,r) + 4*Dipole(11,r))) + (x**2 + y**2 + z**2)*Dipole(12,r))*sqrt(3)))/2., \
	(1,1,1,-1,2,1) : lambda r,x,y,z,Dipole: (y*z*(x**2*z**2*Dipole(10,r) - y**2*Dipole(12,r) - x**2*(x**2 + y**2)*Dipole(9,r)*sqrt(3) + (x**2 + y**2)*((-3*x**2 + y**2 + z**2)*Dipole(11,r) + x**2*Dipole(15,r)*sqrt(3))))/(x**2 + y**2), \
	(1,1,1,-1,2,2) : lambda r,x,y,z,Dipole: (x*(x - y)*y*(x + y)*(-(z**2*(x**2 + y**2 + 2*z**2)*Dipole(10,r)) - 4*(x**2 + y**2)**2*Dipole(11,r) + (x**2 + y**2 - 2*z**2)*Dipole(12,r) - (x**2 + y**2)**2*(Dipole(9,r) - Dipole(15,r))*sqrt(3)))/(2.*(x**2 + y**2)**2), \
	(1,1,1,0,0,0) : lambda r,x,y,z,Dipole: x*z*(-Dipole(6,r) + Dipole(13,r)), \
	(1,1,1,0,1,-1) : lambda r,x,y,z,Dipole: x*y*z*(-Dipole(7,r) - 2*Dipole(8,r) + Dipole(14,r)), \
	(1,1,1,0,1,0) : lambda r,x,y,z,Dipole: x*((x**2 + y**2)*Dipole(8,r) - z**2*(Dipole(7,r) + Dipole(8,r) - Dipole(14,r))), \
	(1,1,1,0,1,1) : lambda r,x,y,z,Dipole: z*(-(x**2*Dipole(7,r)) + (-x**2 + y**2 + z**2)*Dipole(8,r) + x**2*Dipole(14,r)), \
	(1,1,1,0,2,-2) : lambda r,x,y,z,Dipole: (y*z*(x**2*(x**2 + y**2 + 2*z**2)*Dipole(10,r) + (x - y)*(x + y)*Dipole(12,r) - x**2*(x**2 + y**2)*Dipole(9,r)*sqrt(3) + (x**2 + y**2)*((-3*x**2 + y**2 + z**2)*Dipole(11,r) + x**2*Dipole(15,r)*sqrt(3))))/(x**2 + y**2), \
	(1,1,1,0,2,-1) : lambda r,x,y,z,Dipole: x*y*((x**2 + y**2)*(Dipole(11,r) - Dipole(12,r)) - z**2*(Dipole(10,r) + 3*Dipole(11,r) + Dipole(12,r) + (Dipole(9,r) - Dipole(15,r))*sqrt(3))), \
	(1,1,1,0,2,0) : lambda r,x,y,z,Dipole: (x*z*((x**2 + y**2 - 2*z**2)*(Dipole(9,r) - Dipole(15,r)) + ((x**2 + y**2)*Dipole(10,r) + 2*(x**2 + y**2 - z**2)*Dipole(11,r))*sqrt(3)))/2., \
	(1,1,1,0,2,1) : lambda r,x,y,z,Dipole: x**4*Dipole(11,r) + (y**2 + z**2)*(z**2*Dipole(11,r) + y**2*Dipole(12,r)) + x**2*(y**2*(Dipole(11,r) + Dipole(12,r)) - z**2*(Dipole(10,r) + 2*Dipole(11,r) + (Dipole(9,r) - Dipole(15,r))*sqrt(3))), \
	(1,1,1,0,2,2) : lambda r,x,y,z,Dipole: (x*z*((x - y)*(x + y)*(x**2 + y**2 + 2*z**2)*Dipole(10,r) - 2*(x**2 + y**2)*(x**2 - 3*y**2 - z**2)*Dipole(11,r) - 4*y**2*Dipole(12,r) + (-x**4 + y**4)*(Dipole(9,r) - Dipole(15,r))*sqrt(3)))/(2.*(x**2 + y**2)), \
	(1,1,1,1,0,0) : lambda r,x,y,z,Dipole: (y**2 + z**2)*Dipole(6,r) + x**2*Dipole(13,r), \
	(1,1,1,1,1,-1) : lambda r,x,y,z,Dipole: y*((y**2 + z**2)*Dipole(7,r) + x**2*(-2*Dipole(8,r) + Dipole(14,r))), \
	(1,1,1,1,1,0) : lambda r,x,y,z,Dipole: z*((y**2 + z**2)*Dipole(7,r) + x**2*(-2*Dipole(8,r) + Dipole(14,r))), \
	(1,1,1,1,1,1) : lambda r,x,y,z,Dipole: x*(y**2 + z**2)*(Dipole(7,r) + 2*Dipole(8,r)) + x**3*Dipole(14,r), \
	(1,1,1,1,2,-2) : lambda r,x,y,z,Dipole: (x*y*(-(x**2*z**2*(x**2 + y**2 + 2*z**2)*Dipole(10,r)) - (y**4 + x**2*(y**2 + 2*z**2))*Dipole(12,r) + (x**2 + y**2)**2*(y**2 + z**2)*Dipole(9,r)*sqrt(3) + (x**2 + y**2)**2*(2*(-x**2 + y**2 + z**2)*Dipole(11,r) + x**2*Dipole(15,r)*sqrt(3))))/(x**2 + y**2)**2, \
	(1,1,1,1,2,-1) : lambda r,x,y,z,Dipole: (y*z*(x**2*z**2*Dipole(10,r) - 4*x**2*(x**2 + y**2)*Dipole(11,r) + (2*x**2 + y**2)*Dipole(12,r) + (x**2 + y**2)*((y**2 + z**2)*Dipole(9,r) + x**2*Dipole(15,r))*sqrt(3)))/(x**2 + y**2), \
	(1,1,1,1,2,0) : lambda r,x,y,z,Dipole: (-((x**2 + y**2 - 2*z**2)*((y**2 + z**2)*Dipole(9,r) + x**2*Dipole(15,r))) - (x**2*z**2*(Dipole(10,r) + 4*Dipole(11,r)) + y**2*(x**2 + y**2 + z**2)*Dipole(12,r))*sqrt(3))/2., \
	(1,1,1,1,2,1) : lambda r,x,y,z,Dipole: (x*z*(x**2*z**2*Dipole(10,r) + 2*(-x**4 + y**4 + (x**2 + y**2)*z**2)*Dipole(11,r) - y**2*Dipole(12,r) + (x**2 + y**2)*((y**2 + z**2)*Dipole(9,r) + x**2*Dipole(15,r))*sqrt(3)))/(x**2 + y**2), \
	(1,1,1,1,2,2) : lambda r,x,y,z,Dipole: (-(x**2*(x - y)*(x + y)*z**2*(x**2 + y**2 + 2*z**2)*Dipole(10,r)) + 4*x**2*(x**2 + y**2)**2*(2*y**2 + z**2)*Dipole(11,r) + y**2*(-x**4 + y**4 + 2*(3*x**2 + y**2)*z**2)*Dipole(12,r) + (x - y)*(x + y)*(x**2 + y**2)**2*((y**2 + z**2)*Dipole(9,r) + x**2*Dipole(15,r))*sqrt(3))/(2.*(x**2 + y**2)**2), \
	(2,-2,1,-1,0,0) : lambda r,x,y,z,Dipole: x*((x**2 - y**2 + z**2)*Dipole(16,r) + y**2*Dipole(27,r)*sqrt(3)), \
	(2,-2,1,-1,1,-1) : lambda r,x,y,z,Dipole: (x*y*(-((x**4 + x**2*y**2 + 2*y**2*z**2)*Dipole(17,r)) + (x**2 + y**2)**2*(x**2 - y**2 + z**2)*Dipole(18,r) - y**2*z**2*(x**2 + y**2 + 2*z**2)*Dipole(20,r) + (x**2 + y**2)**2*((x**2 - y**2 + z**2)*Dipole(28,r) + ((x**2 + z**2)*Dipole(19,r) + y**2*Dipole(29,r))*sqrt(3))))/(x**2 + y**2)**2, \
	(2,-2,1,-1,1,0) : lambda r,x,y,z,Dipole: (x*z*((-x**2 + y**2)*Dipole(17,r) + (x**2 + y**2)*(x**2 - y**2 + z**2)*Dipole(18,r) + y**2*((x**2 + y**2 + 2*z**2)*Dipole(20,r) - (x**2 + y**2)*(2*Dipole(28,r) + (Dipole(19,r) - Dipole(29,r))*sqrt(3)))))/(x**2 + y**2), \
	(2,-2,1,-1,1,1) : lambda r,x,y,z,Dipole: ((x**2*y**4 + y**4*z**2 + x**4*(y**2 + z**2))*Dipole(17,r) + x**2*(x**2 + y**2)**2*(x**2 - y**2 + z**2)*Dipole(18,r) + y**2*(-(x**2*z**2*(x**2 + y**2 + 2*z**2)*Dipole(20,r)) + (x**2 + y**2)**2*((-x**2 + y**2 + z**2)*Dipole(28,r) + x**2*(-Dipole(19,r) + Dipole(29,r))*sqrt(3))))/(x**2 + y**2)**2, \
	(2,-2,1,-1,2,-2) : lambda r,x,y,z,Dipole: (y*(-((x - y)*(x + y)*(x**4 + y**2*z**2 + x**2*(y - z)*(y + z))*Dipole(21,r)) + y**4*(y**2 + z**2)*(2*z**2*Dipole(25,r) + y**2*Dipole(30,r)) + x**8*(-Dipole(25,r) + Dipole(30,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3)) + x**6*(z**2*Dipole(25,r) + z**2*Dipole(30,r) + 3*y**2*Dipole(31,r) + (y**2 + z**2)*(Dipole(22,r) + Dipole(24,r))*sqrt(3)) + x**2*y**2*(-4*z**4*(Dipole(23,r) + Dipole(26,r)) + y**4*(2*Dipole(25,r) + 3*Dipole(31,r) - (Dipole(22,r) + Dipole(24,r))*sqrt(3)) + y**2*z**2*(-2*Dipole(23,r) + 3*Dipole(25,r) - 2*Dipole(26,r) + 3*Dipole(30,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3))) + x**4*(2*z**4*Dipole(25,r) - y**4*(-3*Dipole(25,r) + 2*Dipole(30,r) - 6*Dipole(31,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3)) + y**2*z**2*(-2*Dipole(23,r) + 2*Dipole(25,r) - 2*Dipole(26,r) + 3*Dipole(30,r) + 2*(Dipole(22,r) + Dipole(24,r))*sqrt(3)))))/(x**2 + y**2)**2, \
	(2,-2,1,-1,2,-1) : lambda r,x,y,z,Dipole: (x*y*z*((-2*x**4 - x**2*y**2 + y**4 - 2*y**2*z**2)*Dipole(21,r) - 2*y**2*z**4*Dipole(26,r) + y**6*(Dipole(26,r) - 3*Dipole(30,r) + 3*Dipole(31,r) - (Dipole(22,r) + Dipole(24,r))*sqrt(3)) + x**6*(-3*Dipole(25,r) + Dipole(30,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3)) + y**4*z**2*(2*Dipole(23,r) - Dipole(25,r) + Dipole(26,r) + Dipole(30,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3)) + x**4*(-3*(2*y**2 + z**2)*Dipole(25,r) + y**2*Dipole(26,r) - y**2*Dipole(30,r) + z**2*Dipole(30,r) + 3*y**2*Dipole(31,r) + (y**2 + z**2)*(Dipole(22,r) + Dipole(24,r))*sqrt(3)) + x**2*y**2*(-(y**2*(3*Dipole(25,r) - 2*Dipole(26,r) + 5*Dipole(30,r) - 6*Dipole(31,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3))) + z**2*(2*Dipole(23,r) - 4*Dipole(25,r) + Dipole(26,r) + 2*Dipole(30,r) + 2*(Dipole(22,r) + Dipole(24,r))*sqrt(3)))))/(x**2 + y**2)**2, \
	(2,-2,1,-1,2,0) : lambda r,x,y,z,Dipole: (x*(-((x**2 + y**2)*(x**2 + y**2 - 2*z**2)*(x**2 - y**2 + z**2)*Dipole(22,r)) + 2*(-x**2 + y**2)*z**2*Dipole(21,r)*sqrt(3) - x**6*Dipole(25,r)*sqrt(3) + 4*y**2*z**4*Dipole(26,r)*sqrt(3) + y**6*(2*Dipole(25,r) - Dipole(31,r))*sqrt(3) - x**4*(z**2*Dipole(25,r) + y**2*Dipole(31,r))*sqrt(3) + y**4*z**2*(-6*Dipole(24,r) + (-2*Dipole(23,r) + 3*Dipole(25,r) + 2*(Dipole(26,r) - 2*Dipole(30,r) + Dipole(31,r)))*sqrt(3)) + x**2*y**2*(y**2*(3*Dipole(25,r) - 2*Dipole(31,r))*sqrt(3) + 2*z**2*(-3*Dipole(24,r) + (-Dipole(23,r) + Dipole(25,r) + Dipole(26,r) - 2*Dipole(30,r) + Dipole(31,r))*sqrt(3)))))/(2.*(x**2 + y**2)), \
	(2,-2,1,-1,2,1) : lambda r,x,y,z,Dipole: (z*((-x**6 + 2*x**2*y**4 + y**4*z**2 + x**4*(y**2 + z**2))*Dipole(21,r) - y**6*(y**2 + z**2)*(2*Dipole(25,r) - Dipole(30,r)) + x**8*(Dipole(25,r) + Dipole(22,r)*sqrt(3)) + x**4*y**2*(z**2*(2*Dipole(23,r) + 2*Dipole(25,r) + Dipole(26,r) + Dipole(30,r)) - (y**2 - 2*z**2)*Dipole(22,r)*sqrt(3) + y**2*(3*Dipole(25,r) + 2*Dipole(26,r) - 5*Dipole(30,r) + 6*Dipole(31,r) - 4*Dipole(24,r)*sqrt(3))) + x**6*((4*y**2 + z**2)*Dipole(25,r) + (y**2 + z**2)*Dipole(22,r)*sqrt(3) + y**2*(Dipole(26,r) - 3*Dipole(30,r) + 3*Dipole(31,r) - 2*Dipole(24,r)*sqrt(3))) + x**2*y**2*(-2*z**4*Dipole(26,r) + y**2*z**2*(2*Dipole(23,r) - Dipole(25,r) + Dipole(26,r) + 2*Dipole(30,r) + Dipole(22,r)*sqrt(3)) - y**4*(2*Dipole(25,r) - Dipole(26,r) + Dipole(30,r) - 3*Dipole(31,r) + (Dipole(22,r) + 2*Dipole(24,r))*sqrt(3)))))/(x**2 + y**2)**2, \
	(2,-2,1,-1,2,2) : lambda r,x,y,z,Dipole: (x*(2*(2*x**2*y**2*(x**2 + y**2) + (x**4 + 3*y**4)*z**2)*Dipole(21,r) + 2*y**4*z**4*(2*Dipole(23,r) + Dipole(25,r) + 2*Dipole(26,r)) + x**8*(-Dipole(25,r) + Dipole(22,r)*sqrt(3)) + y**8*(-2*Dipole(25,r) + 4*Dipole(30,r) - 3*Dipole(31,r) + Dipole(22,r)*sqrt(3)) + y**6*z**2*(2*Dipole(23,r) + Dipole(25,r) + 2*Dipole(26,r) - (Dipole(22,r) + 2*Dipole(24,r))*sqrt(3)) + x**6*(z**2*(-3*Dipole(25,r) + Dipole(22,r)*sqrt(3)) + y**2*(Dipole(25,r) - 4*Dipole(30,r) + 3*Dipole(31,r) - 4*Dipole(24,r)*sqrt(3))) + x**2*y**2*(-4*z**4*(Dipole(23,r) + 2*Dipole(25,r) + Dipole(26,r)) - y**4*(Dipole(25,r) - 4*Dipole(30,r) + 3*Dipole(31,r) + 4*Dipole(24,r)*sqrt(3)) - y**2*z**2*(5*Dipole(25,r) + (Dipole(22,r) + 4*Dipole(24,r))*sqrt(3))) - x**4*(2*z**4*Dipole(25,r) + y**2*z**2*(2*Dipole(23,r) + 9*Dipole(25,r) + 2*Dipole(26,r) - Dipole(22,r)*sqrt(3) + 2*Dipole(24,r)*sqrt(3)) + y**4*(-3*Dipole(25,r) + 4*Dipole(30,r) - 3*Dipole(31,r) + 2*(Dipole(22,r) + 4*Dipole(24,r))*sqrt(3)))))/(2.*(x**2 + y**2)**2), \
	(2,-2,1,0,0,0) : lambda r,x,y,z,Dipole: x*y*z*(-2*Dipole(16,r) + Dipole(27,r)*sqrt(3)), \
	(2,-2,1,0,1,-1) : lambda r,x,y,z,Dipole: (x*z*((-x**2 + y**2)*Dipole(17,r) - 2*y**2*(x**2 + y**2)*Dipole(18,r) + (x**2 + y**2)*(x**2 - y**2 + z**2)*Dipole(28,r) - y**2*(x**2 + y**2)*Dipole(19,r)*sqrt(3) + y**2*((x**2 + y**2 + 2*z**2)*Dipole(20,r) + (x**2 + y**2)*Dipole(29,r)*sqrt(3))))/(x**2 + y**2), \
	(2,-2,1,0,1,0) : lambda r,x,y,z,Dipole: x*y*((x**2 + y**2)*(-Dipole(20,r) + Dipole(19,r)*sqrt(3)) + z**2*(-2*(Dipole(18,r) + Dipole(20,r) + Dipole(28,r)) + Dipole(29,r)*sqrt(3))), \
	(2,-2,1,0,1,1) : lambda r,x,y,z,Dipole: (y*z*((x - y)*(x + y)*Dipole(17,r) - 2*x**2*(x**2 + y**2)*Dipole(18,r) + (x**2 + y**2)*(-x**2 + y**2 + z**2)*Dipole(28,r) + x**2*((x**2 + y**2 + 2*z**2)*Dipole(20,r) - (x**2 + y**2)*(Dipole(19,r) - Dipole(29,r))*sqrt(3))))/(x**2 + y**2), \
	(2,-2,1,0,2,-2) : lambda r,x,y,z,Dipole: (z*(-((x**2 - y**2)**2*Dipole(21,r)) + 2*x**2*y**2*(x**2 + y**2 + 2*z**2)*Dipole(23,r) - x**6*Dipole(25,r) + 2*x**4*y**2*Dipole(25,r) + 2*x**2*y**4*Dipole(25,r) - y**6*Dipole(25,r) + 4*x**2*y**2*z**2*Dipole(25,r) + x**2*z**4*Dipole(25,r) + y**2*z**4*Dipole(25,r) + 2*x**4*y**2*Dipole(26,r) + 2*x**2*y**4*Dipole(26,r) + 4*x**2*y**2*z**2*Dipole(26,r) + x**6*Dipole(30,r) - x**4*y**2*Dipole(30,r) - x**2*y**4*Dipole(30,r) + y**6*Dipole(30,r) + x**4*z**2*Dipole(30,r) + 2*x**2*y**2*z**2*Dipole(30,r) + y**4*z**2*Dipole(30,r) + 3*x**4*y**2*Dipole(31,r) + 3*x**2*y**4*Dipole(31,r) - 2*x**2*y**2*(x**2 + y**2)*(Dipole(22,r) + Dipole(24,r))*sqrt(3)))/(x**2 + y**2), \
	(2,-2,1,0,2,-1) : lambda r,x,y,z,Dipole: (x*((-x**2 + y**2)*z**2*Dipole(21,r) - 2*y**2*(x**2 + y**2)*z**2*Dipole(23,r) + x**6*Dipole(25,r) + x**4*y**2*Dipole(25,r) - x**2*y**4*Dipole(25,r) - y**6*Dipole(25,r) - x**2*y**2*z**2*Dipole(25,r) - y**4*z**2*Dipole(25,r) - x**2*z**4*Dipole(25,r) - y**2*z**4*Dipole(25,r) - x**4*y**2*Dipole(26,r) - 2*x**2*y**4*Dipole(26,r) - y**6*Dipole(26,r) - x**2*y**2*z**2*Dipole(26,r) - y**4*z**2*Dipole(26,r) + 2*y**2*z**4*Dipole(26,r) + x**4*z**2*Dipole(30,r) - 2*x**2*y**2*z**2*Dipole(30,r) - 3*y**4*z**2*Dipole(30,r) + x**2*z**4*Dipole(30,r) + y**2*z**4*Dipole(30,r) + 3*x**2*y**2*z**2*Dipole(31,r) + 3*y**4*z**2*Dipole(31,r) + y**2*(x**2 + y**2)*(-2*z**2*Dipole(22,r) + (x**2 + y**2 - z**2)*Dipole(24,r))*sqrt(3)))/(x**2 + y**2), \
	(2,-2,1,0,2,0) : lambda r,x,y,z,Dipole: (x*y*z*(2*(x**2 + y**2 - 2*z**2)*Dipole(22,r) + 6*(x**2 + y**2)*Dipole(24,r) + (2*(x**2 + y**2)*Dipole(23,r) + (x**2 + y**2 + 2*z**2)*Dipole(25,r) - 2*(x**2 + y**2 + 2*z**2)*Dipole(26,r) - 4*z**2*Dipole(30,r) - (x**2 + y**2 - 2*z**2)*Dipole(31,r))*sqrt(3)))/2., \
	(2,-2,1,0,2,1) : lambda r,x,y,z,Dipole: (y*((x - y)*(x + y)*z**2*Dipole(21,r) + y**2*((y**4 - z**4)*Dipole(25,r) + z**2*(y**2 + z**2)*Dipole(30,r)) - x**6*(Dipole(25,r) + Dipole(26,r) - Dipole(24,r)*sqrt(3)) - x**2*(z**4*(Dipole(25,r) - 2*Dipole(26,r) - Dipole(30,r)) - y**4*(Dipole(25,r) - Dipole(26,r) + Dipole(24,r)*sqrt(3)) + y**2*z**2*(2*Dipole(23,r) + Dipole(25,r) + Dipole(26,r) + 2*Dipole(30,r) - 3*Dipole(31,r) + (2*Dipole(22,r) + Dipole(24,r))*sqrt(3))) - x**4*(y**2*(Dipole(25,r) + 2*Dipole(26,r) - 2*Dipole(24,r)*sqrt(3)) + z**2*(2*Dipole(23,r) + Dipole(25,r) + Dipole(26,r) + 3*Dipole(30,r) - 3*Dipole(31,r) + (2*Dipole(22,r) + Dipole(24,r))*sqrt(3)))))/(x**2 + y**2), \
	(2,-2,1,0,2,2) : lambda r,x,y,z,Dipole: (x*(x - y)*y*(x + y)*z*(4*Dipole(21,r) + 2*(x**2 + y**2 + 2*z**2)*Dipole(23,r) + (5*(x**2 + y**2) + 4*z**2)*Dipole(25,r) + 2*x**2*Dipole(26,r) + 2*y**2*Dipole(26,r) + 4*z**2*Dipole(26,r) - 4*x**2*Dipole(30,r) - 4*y**2*Dipole(30,r) + 3*x**2*Dipole(31,r) + 3*y**2*Dipole(31,r) - 2*(x**2 + y**2)*(Dipole(22,r) + Dipole(24,r))*sqrt(3)))/(2.*(x**2 + y**2)), \
	(2,-2,1,1,0,0) : lambda r,x,y,z,Dipole: y*((-x**2 + y**2 + z**2)*Dipole(16,r) + x**2*Dipole(27,r)*sqrt(3)), \
	(2,-2,1,1,1,-1) : lambda r,x,y,z,Dipole: ((x**2*y**4 + y**4*z**2 + x**4*(y**2 + z**2))*Dipole(17,r) + y**2*(x**2 + y**2)**2*(-x**2 + y**2 + z**2)*Dipole(18,r) + x**2*(-(y**2*z**2*(x**2 + y**2 + 2*z**2)*Dipole(20,r)) + (x**2 + y**2)**2*((x**2 - y**2 + z**2)*Dipole(28,r) + y**2*(-Dipole(19,r) + Dipole(29,r))*sqrt(3))))/(x**2 + y**2)**2, \
	(2,-2,1,1,1,0) : lambda r,x,y,z,Dipole: (y*z*((x - y)*(x + y)*Dipole(17,r) + (-x**4 + y**4 + (x**2 + y**2)*z**2)*Dipole(18,r) + x**2*((x**2 + y**2 + 2*z**2)*Dipole(20,r) - (x**2 + y**2)*(2*Dipole(28,r) + (Dipole(19,r) - Dipole(29,r))*sqrt(3)))))/(x**2 + y**2), \
	(2,-2,1,1,1,1) : lambda r,x,y,z,Dipole: (x*y*(-((y**4 + x**2*(y**2 + 2*z**2))*Dipole(17,r)) - (x**2 + y**2)**2*(x**2 - y**2 - z**2)*Dipole(18,r) - x**2*z**2*(x**2 + y**2 + 2*z**2)*Dipole(20,r) + (x**2 + y**2)**2*((y**2 + z**2)*(Dipole(28,r) + Dipole(19,r)*sqrt(3)) + x**2*(-Dipole(28,r) + Dipole(29,r)*sqrt(3)))))/(x**2 + y**2)**2, \
	(2,-2,1,1,2,-2) : lambda r,x,y,z,Dipole: (x*((x - y)*(x + y)*(y**4 - y**2*z**2 + x**2*(y**2 + z**2))*Dipole(21,r) + x**8*Dipole(30,r) + x**6*(2*(y**2 + z**2)*Dipole(25,r) + z**2*Dipole(30,r) + y**2*(3*Dipole(31,r) - (Dipole(22,r) + Dipole(24,r))*sqrt(3))) + y**4*(y**2 + z**2)*(2*z**2*Dipole(25,r) + y**2*(-Dipole(25,r) + Dipole(30,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3))) + x**4*(2*z**4*Dipole(25,r) + y**2*z**2*(-2*Dipole(23,r) + 3*Dipole(25,r) - 2*Dipole(26,r) + 3*Dipole(30,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3)) - y**4*(-3*Dipole(25,r) + 2*Dipole(30,r) - 6*Dipole(31,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3))) + x**2*y**2*(-4*z**4*(Dipole(23,r) + Dipole(26,r)) + y**4*(3*Dipole(31,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3)) + y**2*z**2*(-2*Dipole(23,r) + 2*Dipole(25,r) - 2*Dipole(26,r) + 3*Dipole(30,r) + 2*(Dipole(22,r) + Dipole(24,r))*sqrt(3)))))/(x**2 + y**2)**2, \
	(2,-2,1,1,2,-1) : lambda r,x,y,z,Dipole: (z*((2*x**4*y**2 + x**2*y**4 - y**6 + (x**4 + y**4)*z**2)*Dipole(21,r) + x**8*(-2*Dipole(25,r) + Dipole(30,r)) + y**6*(y**2 + z**2)*(Dipole(25,r) + Dipole(22,r)*sqrt(3)) + x**2*y**2*(-2*z**4*Dipole(26,r) + y**2*z**2*(2*Dipole(23,r) + 2*Dipole(25,r) + Dipole(26,r) + Dipole(30,r) + 2*Dipole(22,r)*sqrt(3)) + y**4*(4*Dipole(25,r) + Dipole(26,r) - 3*Dipole(30,r) + 3*Dipole(31,r) + (Dipole(22,r) - 2*Dipole(24,r))*sqrt(3))) + x**4*y**2*(z**2*(2*Dipole(23,r) - Dipole(25,r) + Dipole(26,r) + 2*Dipole(30,r)) + (-y**2 + z**2)*Dipole(22,r)*sqrt(3) + y**2*(3*Dipole(25,r) + 2*Dipole(26,r) - 5*Dipole(30,r) + 6*Dipole(31,r) - 4*Dipole(24,r)*sqrt(3))) - x**6*(z**2*(2*Dipole(25,r) - Dipole(30,r)) + y**2*(2*Dipole(25,r) - Dipole(26,r) + Dipole(30,r) - 3*Dipole(31,r) + (Dipole(22,r) + 2*Dipole(24,r))*sqrt(3)))))/(x**2 + y**2)**2, \
	(2,-2,1,1,2,0) : lambda r,x,y,z,Dipole: (y*((x**2 + y**2)*(x**2 + y**2 - 2*z**2)*(x**2 - y**2 - z**2)*Dipole(22,r) + 2*(x - y)*(x + y)*z**2*Dipole(21,r)*sqrt(3) - y**4*(y**2 + z**2)*Dipole(25,r)*sqrt(3) + x**6*(2*Dipole(25,r) - Dipole(31,r))*sqrt(3) + x**4*(-6*z**2*Dipole(24,r) + (y**2*(3*Dipole(25,r) - 2*Dipole(31,r)) + z**2*(-2*Dipole(23,r) + 3*Dipole(25,r) + 2*(Dipole(26,r) - 2*Dipole(30,r) + Dipole(31,r))))*sqrt(3)) + x**2*(4*z**4*Dipole(26,r)*sqrt(3) - y**4*Dipole(31,r)*sqrt(3) + 2*y**2*z**2*(-3*Dipole(24,r) + (-Dipole(23,r) + Dipole(25,r) + Dipole(26,r) - 2*Dipole(30,r) + Dipole(31,r))*sqrt(3)))))/(2.*(x**2 + y**2)), \
	(2,-2,1,1,2,1) : lambda r,x,y,z,Dipole: (x*y*z*((x**4 - 2*y**4 - x**2*(y**2 + 2*z**2))*Dipole(21,r) - 2*x**2*z**4*Dipole(26,r) + (x**2 + y**2)**2*(-3*y**2*Dipole(25,r) + x**2*Dipole(26,r) - 3*x**2*Dipole(30,r) + y**2*Dipole(30,r) + 3*x**2*Dipole(31,r) + (-x + y)*(x + y)*(Dipole(22,r) + Dipole(24,r))*sqrt(3)) + (x**2 + y**2)*z**2*(2*x**2*Dipole(23,r) - (x**2 + 3*y**2)*Dipole(25,r) + x**2*Dipole(26,r) + x**2*Dipole(30,r) + y**2*Dipole(30,r) + (x**2 + y**2)*(Dipole(22,r) + Dipole(24,r))*sqrt(3))))/(x**2 + y**2)**2, \
	(2,-2,1,1,2,2) : lambda r,x,y,z,Dipole: (y*(-2*(2*x**2*y**2*(x**2 + y**2) + (3*x**4 + y**4)*z**2)*Dipole(21,r) + 2*z**4*(2*x**2*(-x**2 + y**2)*Dipole(23,r) + (-x**4 + 4*x**2*y**2 + y**4)*Dipole(25,r) + 2*x**2*(-x**2 + y**2)*Dipole(26,r)) - (x**2 + y**2)**2*(-((x - y)*(x + y)*((2*x**2 - y**2)*Dipole(25,r) + x**2*(-4*Dipole(30,r) + 3*Dipole(31,r)))) + (x**2 - y**2)**2*Dipole(22,r)*sqrt(3) - 4*x**2*y**2*Dipole(24,r)*sqrt(3)) + (x**2 + y**2)*z**2*(y**4*(3*Dipole(25,r) - Dipole(22,r)*sqrt(3)) + 2*x**2*y**2*(Dipole(23,r) + 3*Dipole(25,r) + Dipole(26,r) + Dipole(24,r)*sqrt(3)) - x**4*(2*Dipole(23,r) + Dipole(25,r) + 2*Dipole(26,r) - (Dipole(22,r) + 2*Dipole(24,r))*sqrt(3)))))/(2.*(x**2 + y**2)**2), \
	(2,-1,1,-1,0,0) : lambda r,x,y,z,Dipole: z*((x**2 - y**2 + z**2)*Dipole(16,r) + y**2*Dipole(27,r)*sqrt(3)), \
	(2,-1,1,-1,1,-1) : lambda r,x,y,z,Dipole: (y*z*(-(x**2*Dipole(17,r)) + (x**2 + y**2)*(x**2 - y**2 + z**2)*Dipole(18,r) + y**2*z**2*Dipole(20,r) + x**4*Dipole(28,r) - y**4*Dipole(28,r) + x**2*z**2*Dipole(28,r) + y**2*z**2*Dipole(28,r) + (x**2 + y**2)*((x**2 + z**2)*Dipole(19,r) + y**2*Dipole(29,r))*sqrt(3)))/(x**2 + y**2), \
	(2,-1,1,-1,1,0) : lambda r,x,y,z,Dipole: x**2*Dipole(17,r) + z**2*(x**2 - y**2 + z**2)*Dipole(18,r) + y**2*((x**2 + y**2)*Dipole(28,r) - z**2*(Dipole(20,r) + Dipole(28,r) + (Dipole(19,r) - Dipole(29,r))*sqrt(3))), \
	(2,-1,1,-1,1,1) : lambda r,x,y,z,Dipole: (x*z*(-(x**2*Dipole(17,r)) + (x**2 + y**2)*(x**2 - y**2 + z**2)*Dipole(18,r) + y**2*(z**2*Dipole(20,r) - (x**2 + y**2)*(2*Dipole(28,r) + (Dipole(19,r) - Dipole(29,r))*sqrt(3)))))/(x**2 + y**2), \
	(2,-1,1,-1,2,-2) : lambda r,x,y,z,Dipole: (x*y*z*(-2*x**2*(x**2 + y**2)*Dipole(21,r) - 2*y**2*z**4*(Dipole(23,r) + Dipole(25,r)) + y**6*(Dipole(23,r) + Dipole(25,r) - 3*Dipole(30,r) + 3*Dipole(31,r) - (Dipole(22,r) + Dipole(24,r))*sqrt(3)) + x**6*(-3*Dipole(25,r) + Dipole(30,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3)) + y**4*z**2*(Dipole(23,r) - 2*Dipole(25,r) + 2*Dipole(26,r) + Dipole(30,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3)) + x**4*(z**2*(-3*Dipole(25,r) + Dipole(30,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3)) + y**2*(Dipole(23,r) - 5*Dipole(25,r) - Dipole(30,r) + 3*Dipole(31,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3))) + x**2*y**2*(-(y**2*(-2*Dipole(23,r) + Dipole(25,r) + 5*Dipole(30,r) - 6*Dipole(31,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3))) + z**2*(Dipole(23,r) - 5*Dipole(25,r) + 2*(Dipole(26,r) + Dipole(30,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3))))))/(x**2 + y**2)**2, \
	(2,-1,1,-1,2,-1) : lambda r,x,y,z,Dipole: (y*(x**2*((x**2 + y**2)**2 - z**4)*Dipole(21,r) + 2*x**6*Dipole(25,r) + y**6*Dipole(30,r) + y**2*z**4*(Dipole(23,r) + Dipole(26,r) + Dipole(30,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3)) - y**4*z**2*(Dipole(23,r) - Dipole(25,r) + Dipole(26,r) + 2*Dipole(30,r) - 3*Dipole(31,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3)) + x**2*(2*y**4*(Dipole(25,r) + Dipole(30,r)) - y**2*z**2*(Dipole(23,r) - 2*Dipole(25,r) + Dipole(26,r) + Dipole(30,r) - 3*Dipole(31,r)) + z**4*(-Dipole(25,r) + Dipole(30,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3))) + x**4*(y**2*(4*Dipole(25,r) + Dipole(30,r)) + z**2*(Dipole(25,r) + Dipole(30,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3)))))/(x**2 + y**2), \
	(2,-1,1,-1,2,0) : lambda r,x,y,z,Dipole: (z*(-((x**2 + y**2 - 2*z**2)*(x**2 - y**2 + z**2)*Dipole(22,r)) - 6*y**2*z**2*Dipole(24,r) + 2*x**2*Dipole(21,r)*sqrt(3) - (-(y**2*(x**2 + y**2 - z**2)*Dipole(23,r)) + ((x**2 + y**2)**2 + x**2*z**2)*Dipole(25,r) + y**2*(2*z**2*(Dipole(26,r) + Dipole(30,r) - Dipole(31,r)) - (x**2 + y**2)*(2*Dipole(30,r) - Dipole(31,r))))*sqrt(3)))/2., \
	(2,-1,1,-1,2,1) : lambda r,x,y,z,Dipole: (x*(x**2*((x**2 + y**2)**2 - z**4)*Dipole(21,r) - y**2*(x**2 + y**2)**2*(2*Dipole(25,r) - Dipole(30,r)) + z**4*(x**2*Dipole(25,r) + y**2*(Dipole(23,r) + 2*Dipole(25,r) + Dipole(26,r)) + (x**2 + y**2)*Dipole(22,r)*sqrt(3)) + (x**2 + y**2)*z**2*((x**2 + y**2)*Dipole(25,r) + (x - y)*(x + y)*Dipole(22,r)*sqrt(3) - y**2*(Dipole(23,r) + Dipole(26,r) + 3*Dipole(30,r) - 3*Dipole(31,r) + 2*Dipole(24,r)*sqrt(3)))))/(x**2 + y**2), \
	(2,-1,1,-1,2,2) : lambda r,x,y,z,Dipole: -(z*(2*(x**6 - x**2*y**4)*Dipole(21,r) - 2*y**4*z**4*Dipole(23,r) + x**8*(Dipole(25,r) - Dipole(22,r)*sqrt(3)) + y**8*(Dipole(23,r) - Dipole(25,r) - 2*Dipole(30,r) + 3*Dipole(31,r) - Dipole(22,r)*sqrt(3)) + y**6*z**2*(Dipole(23,r) - 2*Dipole(25,r) + 2*(Dipole(26,r) + Dipole(30,r)) + (Dipole(22,r) + 2*Dipole(24,r))*sqrt(3)) - x**6*(z**2*(-3*Dipole(25,r) + Dipole(22,r)*sqrt(3)) + y**2*(Dipole(23,r) + 6*Dipole(25,r) - 6*Dipole(30,r) + 3*Dipole(31,r) - 4*Dipole(24,r)*sqrt(3))) + x**2*y**2*(2*z**4*(Dipole(23,r) + 3*Dipole(25,r)) + y**4*(Dipole(23,r) - 10*Dipole(25,r) + 2*Dipole(30,r) + 3*Dipole(31,r) + 4*Dipole(24,r)*sqrt(3)) + y**2*z**2*(-3*Dipole(25,r) + 4*Dipole(30,r) + (Dipole(22,r) + 4*Dipole(24,r))*sqrt(3))) + x**4*(2*z**4*Dipole(25,r) - y**4*(Dipole(23,r) + 16*Dipole(25,r) - 10*Dipole(30,r) + 3*Dipole(31,r) - 2*(Dipole(22,r) + 4*Dipole(24,r))*sqrt(3)) - y**2*z**2*(Dipole(23,r) + Dipole(22,r)*sqrt(3) - 2*(Dipole(25,r) - Dipole(26,r) + Dipole(30,r) + Dipole(24,r)*sqrt(3))))))/(2.*(x**2 + y**2)**2), \
	(2,-1,1,0,0,0) : lambda r,x,y,z,Dipole: y*((x**2 + y**2 - z**2)*Dipole(16,r) + z**2*Dipole(27,r)*sqrt(3)), \
	(2,-1,1,0,1,-1) : lambda r,x,y,z,Dipole: x**2*Dipole(17,r) + y**2*(x**2 + y**2 - z**2)*Dipole(18,r) + z**2*((x**2 + z**2)*Dipole(28,r) - y**2*(Dipole(20,r) + Dipole(28,r) + (Dipole(19,r) - Dipole(29,r))*sqrt(3))), \
	(2,-1,1,0,1,0) : lambda r,x,y,z,Dipole: y*z*((x**2 + y**2 - z**2)*Dipole(18,r) + (x**2 + y**2)*Dipole(20,r) + (x**2 + y**2 - z**2)*Dipole(28,r) + ((x**2 + y**2)*Dipole(19,r) + z**2*Dipole(29,r))*sqrt(3)), \
	(2,-1,1,0,1,1) : lambda r,x,y,z,Dipole: x*y*(-Dipole(17,r) + (x**2 + y**2 - z**2)*Dipole(18,r) - z**2*(Dipole(20,r) + 2*Dipole(28,r) + (Dipole(19,r) - Dipole(29,r))*sqrt(3))), \
	(2,-1,1,0,2,-2) : lambda r,x,y,z,Dipole: (x*((x**4 - y**4)*Dipole(21,r) - y**2*(x**2 + y**2 - z**2)*(x**2 + y**2 + 2*z**2)*Dipole(23,r) + z**2*((-2*x**4 + y**4 - x**2*(y**2 + 2*z**2))*Dipole(25,r) + (x**2 + y**2)*(-2*y**2*Dipole(26,r) + (x**2 - 3*y**2 + z**2)*Dipole(30,r) + 3*y**2*Dipole(31,r))) + y**2*(x**2 + y**2)*((x**2 + y**2 - z**2)*Dipole(22,r) - 2*z**2*Dipole(24,r))*sqrt(3)))/(x**2 + y**2), \
	(2,-1,1,0,2,-1) : lambda r,x,y,z,Dipole: z*(x**2*Dipole(21,r) + (x**2 + z**2)*(2*x**2*Dipole(25,r) + z**2*Dipole(30,r)) + y**4*(Dipole(23,r) + Dipole(26,r) + Dipole(30,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3)) + y**2*(x**2*(Dipole(23,r) + 2*Dipole(25,r) + Dipole(26,r) + Dipole(30,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3)) - z**2*(Dipole(23,r) - Dipole(25,r) + Dipole(26,r) + 2*Dipole(30,r) - 3*Dipole(31,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3)))), \
	(2,-1,1,0,2,0) : lambda r,x,y,z,Dipole: (y*(-((x**2 + y**2 - 2*z**2)*(x**2 + y**2 - z**2)*Dipole(22,r)) + 6*(x**2 + y**2)*z**2*Dipole(24,r) - ((x**2 + y**2)*(x**2 + y**2 - z**2)*Dipole(23,r) + z**2*((x**2 + y**2)*Dipole(25,r) - 2*(x**2 + y**2)*Dipole(26,r) - 2*(x**2 + y**2 - z**2)*Dipole(30,r) + (x**2 + y**2 - 2*z**2)*Dipole(31,r)))*sqrt(3)))/2., \
	(2,-1,1,0,2,1) : lambda r,x,y,z,Dipole: x*y*z*(-Dipole(21,r) + (x**2 + y**2 - z**2)*Dipole(23,r) - (2*(x**2 + y**2) + z**2)*Dipole(25,r) + x**2*Dipole(26,r) + y**2*Dipole(26,r) - z**2*Dipole(26,r) + x**2*Dipole(30,r) + y**2*Dipole(30,r) - 3*z**2*Dipole(30,r) + 3*z**2*Dipole(31,r) + (x**2 + y**2 - z**2)*(Dipole(22,r) + Dipole(24,r))*sqrt(3)), \
	(2,-1,1,0,2,2) : lambda r,x,y,z,Dipole: (y*(-4*x**2*(x**2 + y**2)*Dipole(21,r) - (x - y)*(x + y)*(x**2 + y**2 - z**2)*(x**2 + y**2 + 2*z**2)*Dipole(23,r) + z**2*((7*x**4 + 8*x**2*y**2 + y**4 + 2*(3*x**2 + y**2)*z**2)*Dipole(25,r) + 2*(-x**4 + y**4)*Dipole(26,r) - 2*(x**2 + y**2)*(3*x**2 - y**2 + z**2)*Dipole(30,r) + 3*(x**4 - y**4)*Dipole(31,r)) + (x**4 - y**4)*((x**2 + y**2 - z**2)*Dipole(22,r) - 2*z**2*Dipole(24,r))*sqrt(3)))/(2.*(x**2 + y**2)), \
	(2,-1,1,1,0,0) : lambda r,x,y,z,Dipole: x*y*z*(-2*Dipole(16,r) + Dipole(27,r)*sqrt(3)), \
	(2,-1,1,1,1,-1) : lambda r,x,y,z,Dipole: (x*z*(-(x**2*Dipole(17,r)) - 2*y**2*(x**2 + y**2)*Dipole(18,r) + (x**2 + y**2)*(x**2 - y**2 + z**2)*Dipole(28,r) - y**2*(x**2 + y**2)*Dipole(19,r)*sqrt(3) + y**2*(z**2*Dipole(20,r) + (x**2 + y**2)*Dipole(29,r)*sqrt(3))))/(x**2 + y**2), \
	(2,-1,1,1,1,0) : lambda r,x,y,z,Dipole: -(x*y*(Dipole(17,r) - (x**2 + y**2)*Dipole(28,r) + z**2*(2*Dipole(18,r) + Dipole(20,r) + Dipole(28,r) + (Dipole(19,r) - Dipole(29,r))*sqrt(3)))), \
	(2,-1,1,1,1,1) : lambda r,x,y,z,Dipole: (y*z*((2*x**2 + y**2)*Dipole(17,r) - x**2*(-(z**2*Dipole(20,r)) + 2*(x**2 + y**2)*(Dipole(18,r) + Dipole(28,r))) + (x**2 + y**2)*((y**2 + z**2)*Dipole(19,r) + x**2*Dipole(29,r))*sqrt(3)))/(x**2 + y**2), \
	(2,-1,1,1,2,-2) : lambda r,x,y,z,Dipole: (z*((-x**6 + x**4*y**2 + 3*x**2*y**4 + y**6)*Dipole(21,r) + z**4*(-2*x**2*y**2*Dipole(23,r) + (x**4 + y**4)*Dipole(25,r)) + (x**2 + y**2)*z**2*(x**2*y**2*Dipole(23,r) + x**2*(y**2*(Dipole(25,r) + 2*Dipole(26,r)) + (x**2 + y**2)*Dipole(30,r)) + y**2*(x**2 + y**2)*Dipole(24,r)*sqrt(3)) + (x**2 + y**2)**2*(x**4*(-Dipole(25,r) + Dipole(30,r)) + y**4*(-Dipole(25,r) + Dipole(24,r)*sqrt(3)) + x**2*y**2*(Dipole(23,r) + 2*Dipole(25,r) - 3*Dipole(30,r) + 3*Dipole(31,r) - (2*Dipole(22,r) + Dipole(24,r))*sqrt(3)))))/(x**2 + y**2)**2, \
	(2,-1,1,1,2,-1) : lambda r,x,y,z,Dipole: (x*(-((y**4 + x**2*(y**2 + z**2))*Dipole(21,r)) + x**6*Dipole(25,r) + y**6*(-Dipole(25,r) + Dipole(30,r)) + y**2*z**4*(Dipole(23,r) + Dipole(26,r) + Dipole(30,r)) + x**4*(y**2*Dipole(25,r) + (y**2 + z**2)*Dipole(30,r)) - y**4*z**2*(Dipole(23,r) + Dipole(26,r) + 2*Dipole(30,r) - 3*Dipole(31,r) + 2*(Dipole(22,r) + Dipole(24,r))*sqrt(3)) - x**2*(y**4*(Dipole(25,r) - 2*Dipole(30,r)) + z**4*(Dipole(25,r) - Dipole(30,r)) + y**2*z**2*(Dipole(23,r) + Dipole(26,r) + Dipole(30,r) - 3*Dipole(31,r) + 2*(Dipole(22,r) + Dipole(24,r))*sqrt(3)))))/(x**2 + y**2), \
	(2,-1,1,1,2,0) : lambda r,x,y,z,Dipole: (x*y*z*(-(z**2*(4*Dipole(22,r) + 6*Dipole(24,r) + (2*Dipole(21,r) + Dipole(23,r) - Dipole(25,r) + 2*(Dipole(26,r) + Dipole(30,r) - Dipole(31,r)))*sqrt(3))) + (x**2 + y**2)*(2*Dipole(22,r) + (-2*Dipole(21,r) + Dipole(23,r) + 2*Dipole(30,r) - Dipole(31,r))*sqrt(3))))/2., \
	(2,-1,1,1,2,1) : lambda r,x,y,z,Dipole: (y*(-((x**4 - y**2*z**2 + x**2*(y**2 - 2*z**2))*Dipole(21,r)) + (x**2 + y**2)**2*((-x**2 + y**2)*Dipole(25,r) + x**2*Dipole(30,r)) + z**4*(-(y**2*Dipole(25,r)) + x**2*(Dipole(23,r) + Dipole(26,r)) + (x**2 + y**2)*Dipole(24,r)*sqrt(3)) - (x**2 + y**2)*z**2*(-(y**2*Dipole(24,r)*sqrt(3)) + x**2*(Dipole(23,r) + Dipole(26,r) + 3*Dipole(30,r) - 3*Dipole(31,r) + (2*Dipole(22,r) + Dipole(24,r))*sqrt(3)))))/(x**2 + y**2), \
	(2,-1,1,1,2,2) : lambda r,x,y,z,Dipole: (x*y*z*(2*(x**2 + y**2)*(3*x**2 + y**2)*Dipole(21,r) + 2*y**2*z**4*(Dipole(23,r) + Dipole(25,r)) + x**6*(Dipole(23,r) + 4*Dipole(25,r) - 6*Dipole(30,r) + 3*Dipole(31,r) - 2*Dipole(22,r)*sqrt(3)) - y**6*(Dipole(23,r) + 4*Dipole(25,r) - 2*Dipole(30,r) + 3*Dipole(31,r) - 2*(Dipole(22,r) + 2*Dipole(24,r))*sqrt(3)) - y**4*z**2*(Dipole(23,r) + Dipole(25,r) + 2*(Dipole(26,r) + Dipole(30,r) - Dipole(24,r)*sqrt(3))) + x**2*(-2*z**4*(Dipole(23,r) + Dipole(25,r)) - 4*y**2*z**2*(Dipole(30,r) - Dipole(24,r)*sqrt(3)) - y**4*(Dipole(23,r) + 4*Dipole(25,r) + 2*Dipole(30,r) + 3*Dipole(31,r) - 2*(Dipole(22,r) + 4*Dipole(24,r))*sqrt(3))) + x**4*(y**2*(Dipole(23,r) + 4*Dipole(25,r) - 10*Dipole(30,r) + 3*Dipole(31,r) - 2*Dipole(22,r)*sqrt(3) + 4*Dipole(24,r)*sqrt(3)) + z**2*(Dipole(23,r) + Dipole(25,r) + 2*(Dipole(26,r) - Dipole(30,r) + Dipole(24,r)*sqrt(3))))))/(2.*(x**2 + y**2)**2), \
	(2,0,1,-1,0,0) : lambda r,x,y,z,Dipole: -(y*((x**2 + y**2 - 2*z**2)*Dipole(27,r) + 2*z**2*Dipole(16,r)*sqrt(3)))/2., \
	(2,0,1,-1,1,-1) : lambda r,x,y,z,Dipole: (-((x**2 + y**2 - 2*z**2)*((x**2 + z**2)*Dipole(19,r) + y**2*Dipole(29,r))) - (x**2*Dipole(17,r) + y**2*z**2*(2*Dipole(18,r) + Dipole(20,r) + 2*Dipole(28,r)))*sqrt(3))/2., \
	(2,0,1,-1,1,0) : lambda r,x,y,z,Dipole: -(y*z**3*(Dipole(19,r) - Dipole(29,r) + Dipole(18,r)*sqrt(3))) + (y*(x**2 + y**2)*z*(Dipole(19,r) - Dipole(29,r) + (Dipole(20,r) + 2*Dipole(28,r))*sqrt(3)))/2., \
	(2,0,1,-1,1,1) : lambda r,x,y,z,Dipole: (x*y*((x**2 + y**2 - 2*z**2)*(Dipole(19,r) - Dipole(29,r)) + (Dipole(17,r) - z**2*(2*Dipole(18,r) + Dipole(20,r) + 2*Dipole(28,r)))*sqrt(3)))/2., \
	(2,0,1,-1,2,-2) : lambda r,x,y,z,Dipole: (x*((-x**4 + x**2*z**2 + 2*z**4)*Dipole(24,r) - x**2*(x**2 + z**2)*Dipole(21,r)*sqrt(3) + 4*z**4*Dipole(23,r)*sqrt(3) + 2*z**2*(-x + z)*(x + z)*Dipole(25,r)*sqrt(3) - (4*x**2*z**4*(Dipole(23,r) + Dipole(25,r))*sqrt(3))/(x**2 + y**2) + y**4*(Dipole(24,r) + (Dipole(21,r) + Dipole(25,r) - Dipole(31,r))*sqrt(3)) + y**2*(x**2*(Dipole(25,r) - Dipole(31,r))*sqrt(3) + z**2*(-6*Dipole(22,r) - 3*Dipole(24,r) + (Dipole(21,r) + 2*(Dipole(23,r) + 2*Dipole(25,r) - Dipole(26,r) - 2*Dipole(30,r) + Dipole(31,r)))*sqrt(3)))))/2., \
	(2,0,1,-1,2,-1) : lambda r,x,y,z,Dipole: (2*z**5*Dipole(24,r) - z**3*(6*y**2*Dipole(22,r) - (x**2 - 3*y**2)*Dipole(24,r) + (x**2*(Dipole(21,r) - 2*Dipole(25,r)) + y**2*(2*Dipole(23,r) + Dipole(26,r) + 2*Dipole(30,r) - 2*Dipole(31,r)))*sqrt(3)) - (x**2 + y**2)*z*((x - y)*(x + y)*Dipole(24,r) + x**2*Dipole(21,r)*sqrt(3) + (-2*x**2*Dipole(25,r) + y**2*(Dipole(25,r) - Dipole(26,r) - 2*Dipole(30,r) + Dipole(31,r)))*sqrt(3)))/2., \
	(2,0,1,-1,2,0) : lambda r,x,y,z,Dipole: (y*((x**2 + y**2)**2*(3*Dipole(25,r) + Dipole(31,r)) + 4*z**4*(Dipole(31,r) - (Dipole(22,r) + Dipole(24,r))*sqrt(3)) + 2*(x**2 + y**2)*z**2*(3*Dipole(23,r) + 3*Dipole(26,r) + 6*Dipole(30,r) - 2*Dipole(31,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3))))/4., \
	(2,0,1,-1,2,1) : lambda r,x,y,z,Dipole: (x*y*z*(Dipole(21,r)*sqrt(3) - z**2*(6*Dipole(22,r) + 4*Dipole(24,r) + (2*Dipole(23,r) + 2*Dipole(25,r) + Dipole(26,r) + 2*Dipole(30,r) - 2*Dipole(31,r))*sqrt(3)) + (x**2 + y**2)*(2*Dipole(24,r) + (-3*Dipole(25,r) + Dipole(26,r) + 2*Dipole(30,r) - Dipole(31,r))*sqrt(3))))/2., \
	(2,0,1,-1,2,2) : lambda r,x,y,z,Dipole: (2*y*(x**2 + y**2)*(3*(-x + y)*(x + y)*z**2*Dipole(22,r) + (x**2 + y**2 - 2*z**2)*(2*x**2 + z**2)*Dipole(24,r)) + y*(4*x**2*(x**2 + y**2)*(x**2 + y**2 + z**2)*Dipole(21,r) + 2*(x - y)*(x + y)*z**2*(x**2 + y**2 + 2*z**2)*Dipole(23,r) + x**6*Dipole(25,r) + x**4*y**2*Dipole(25,r) - x**2*y**4*Dipole(25,r) - y**6*Dipole(25,r) + 10*x**4*z**2*Dipole(25,r) + 8*x**2*y**2*z**2*Dipole(25,r) - 2*y**4*z**2*Dipole(25,r) + 8*x**2*z**4*Dipole(25,r) - 2*x**4*z**2*Dipole(26,r) + 2*y**4*z**2*Dipole(26,r) - 4*x**4*z**2*Dipole(30,r) + 4*y**4*z**2*Dipole(30,r) - (x**4 - y**4)*(x**2 + y**2 - 2*z**2)*Dipole(31,r))*sqrt(3))/(4.*(x**2 + y**2)), \
	(2,0,1,0,0,0) : lambda r,x,y,z,Dipole: z**3*Dipole(27,r) - ((x**2 + y**2)*z*(Dipole(27,r) - 2*Dipole(16,r)*sqrt(3)))/2., \
	(2,0,1,0,1,-1) : lambda r,x,y,z,Dipole: (y*(x**2 + y**2)*z*(Dipole(19,r) - Dipole(29,r) + (2*Dipole(18,r) + Dipole(20,r))*sqrt(3)))/2. - y*z**3*(Dipole(19,r) - Dipole(29,r) + Dipole(28,r)*sqrt(3)), \
	(2,0,1,0,1,0) : lambda r,x,y,z,Dipole: (-((x**2 + y**2 - 2*z**2)*((x**2 + y**2)*Dipole(19,r) + z**2*Dipole(29,r))) - (x**2 + y**2)*((x**2 + y**2)*Dipole(20,r) - 2*z**2*(Dipole(18,r) + Dipole(28,r)))*sqrt(3))/2., \
	(2,0,1,0,1,1) : lambda r,x,y,z,Dipole: (x*(x**2 + y**2)*z*(Dipole(19,r) - Dipole(29,r) + (2*Dipole(18,r) + Dipole(20,r))*sqrt(3)))/2. - x*z**3*(Dipole(19,r) - Dipole(29,r) + Dipole(28,r)*sqrt(3)), \
	(2,0,1,0,2,-2) : lambda r,x,y,z,Dipole: (x*y*z*(6*(x**2 + y**2)*Dipole(22,r) + 2*(x**2 + y**2 - 2*z**2)*Dipole(24,r) + (-2*(x**2 + y**2 + 2*z**2)*Dipole(23,r) + (x**2 + y**2 + 2*z**2)*Dipole(25,r) + 2*(x**2 + y**2)*Dipole(26,r) - 4*z**2*Dipole(30,r) - (x**2 + y**2 - 2*z**2)*Dipole(31,r))*sqrt(3)))/2., \
	(2,0,1,0,2,-1) : lambda r,x,y,z,Dipole: -(y*(-6*(x**2 + y**2)*z**2*Dipole(22,r) + (x**2 + y**2 - 2*z**2)*(x**2 + y**2 - z**2)*Dipole(24,r) - 2*(x**2 + y**2)*z**2*Dipole(23,r)*sqrt(3) + ((x**2 + y**2)*z**2*Dipole(25,r) + (x**2 + y**2 - z**2)*((x**2 + y**2)*Dipole(26,r) - 2*z**2*Dipole(30,r)) + z**2*(x**2 + y**2 - 2*z**2)*Dipole(31,r))*sqrt(3)))/2., \
	(2,0,1,0,2,0) : lambda r,x,y,z,Dipole: z**5*Dipole(31,r) + (x**2 + y**2)*z**3*(3*Dipole(30,r) - Dipole(31,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3)) - ((x**2 + y**2)**2*z*(6*Dipole(23,r) - 3*Dipole(25,r) + 6*Dipole(26,r) - Dipole(31,r) + 2*(Dipole(22,r) + Dipole(24,r))*sqrt(3)))/4., \
	(2,0,1,0,2,1) : lambda r,x,y,z,Dipole: -(x*(-6*(x**2 + y**2)*z**2*Dipole(22,r) + (x**2 + y**2 - 2*z**2)*(x**2 + y**2 - z**2)*Dipole(24,r) - 2*(x**2 + y**2)*z**2*Dipole(23,r)*sqrt(3) + ((x**2 + y**2)*z**2*Dipole(25,r) + (x**2 + y**2 - z**2)*((x**2 + y**2)*Dipole(26,r) - 2*z**2*Dipole(30,r)) + z**2*(x**2 + y**2 - 2*z**2)*Dipole(31,r))*sqrt(3)))/2., \
	(2,0,1,0,2,2) : lambda r,x,y,z,Dipole: ((x - y)*(x + y)*z*(6*(x**2 + y**2)*Dipole(22,r) + 2*(x**2 + y**2 - 2*z**2)*Dipole(24,r) + (-2*(x**2 + y**2 + 2*z**2)*Dipole(23,r) + (x**2 + y**2 + 2*z**2)*Dipole(25,r) + 2*(x**2 + y**2)*Dipole(26,r) - 4*z**2*Dipole(30,r) - (x**2 + y**2 - 2*z**2)*Dipole(31,r))*sqrt(3)))/4., \
	(2,0,1,1,0,0) : lambda r,x,y,z,Dipole: -(x*((x**2 + y**2 - 2*z**2)*Dipole(27,r) + 2*z**2*Dipole(16,r)*sqrt(3)))/2., \
	(2,0,1,1,1,-1) : lambda r,x,y,z,Dipole: (x*y*((x**2 + y**2 - 2*z**2)*(Dipole(19,r) - Dipole(29,r)) + (Dipole(17,r) - z**2*(2*Dipole(18,r) + Dipole(20,r) + 2*Dipole(28,r)))*sqrt(3)))/2., \
	(2,0,1,1,1,0) : lambda r,x,y,z,Dipole: -(x*z**3*(Dipole(19,r) - Dipole(29,r) + Dipole(18,r)*sqrt(3))) + (x*(x**2 + y**2)*z*(Dipole(19,r) - Dipole(29,r) + (Dipole(20,r) + 2*Dipole(28,r))*sqrt(3)))/2., \
	(2,0,1,1,1,1) : lambda r,x,y,z,Dipole: (-((x**2 + y**2 - 2*z**2)*((y**2 + z**2)*Dipole(19,r) + x**2*Dipole(29,r))) - (y**2*Dipole(17,r) + x**2*z**2*(2*Dipole(18,r) + Dipole(20,r) + 2*Dipole(28,r)))*sqrt(3))/2., \
	(2,0,1,1,2,-2) : lambda r,x,y,z,Dipole: (y*(2*z**4*((x**2 + y**2)*Dipole(24,r) + (2*x**2*Dipole(23,r) + (x - y)*(x + y)*Dipole(25,r))*sqrt(3)) + (x**2 + y**2)**2*((x - y)*(x + y)*Dipole(24,r) + ((x - y)*(x + y)*Dipole(21,r) + x**2*(Dipole(25,r) - Dipole(31,r)))*sqrt(3)) + (x**2 + y**2)*z**2*(-6*x**2*Dipole(22,r) + (-3*x**2 + y**2)*Dipole(24,r) + (x - y)*(x + y)*Dipole(21,r)*sqrt(3) + 2*(-(y**2*Dipole(25,r)) + x**2*(Dipole(23,r) + 2*Dipole(25,r) - Dipole(26,r) - 2*Dipole(30,r) + Dipole(31,r)))*sqrt(3))))/(2.*(x**2 + y**2)), \
	(2,0,1,1,2,-1) : lambda r,x,y,z,Dipole: (x*y*z*(Dipole(21,r)*sqrt(3) - z**2*(6*Dipole(22,r) + 4*Dipole(24,r) + (2*Dipole(23,r) + 2*Dipole(25,r) + Dipole(26,r) + 2*Dipole(30,r) - 2*Dipole(31,r))*sqrt(3)) + (x**2 + y**2)*(2*Dipole(24,r) + (-3*Dipole(25,r) + Dipole(26,r) + 2*Dipole(30,r) - Dipole(31,r))*sqrt(3))))/2., \
	(2,0,1,1,2,0) : lambda r,x,y,z,Dipole: (x*((x**2 + y**2)**2*(3*Dipole(25,r) + Dipole(31,r)) + 4*z**4*(Dipole(31,r) - (Dipole(22,r) + Dipole(24,r))*sqrt(3)) + 2*(x**2 + y**2)*z**2*(3*Dipole(23,r) + 3*Dipole(26,r) + 6*Dipole(30,r) - 2*Dipole(31,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3))))/4., \
	(2,0,1,1,2,1) : lambda r,x,y,z,Dipole: (2*z**5*Dipole(24,r) - z**3*(6*x**2*Dipole(22,r) + (3*x**2 - y**2)*Dipole(24,r) + (y**2*(Dipole(21,r) - 2*Dipole(25,r)) + x**2*(2*Dipole(23,r) + Dipole(26,r) + 2*Dipole(30,r) - 2*Dipole(31,r)))*sqrt(3)) - (x**2 + y**2)*z*((-x**2 + y**2)*Dipole(24,r) + (y**2*(Dipole(21,r) - 2*Dipole(25,r)) + x**2*(Dipole(25,r) - Dipole(26,r) - 2*Dipole(30,r) + Dipole(31,r)))*sqrt(3)))/2., \
	(2,0,1,1,2,2) : lambda r,x,y,z,Dipole: (x*(6*(-x**4 + y**4)*z**2*Dipole(22,r) - 4*y**2*(x**2 + y**2)*Dipole(21,r)*sqrt(3) + 4*z**4*((x**2 + y**2)*Dipole(24,r) + ((x - y)*(x + y)*Dipole(23,r) - 2*y**2*Dipole(25,r))*sqrt(3)) + 2*(x**2 + y**2)*z**2*(-((x**2 - 3*y**2)*Dipole(24,r)) + (x - y)*(x + y)*Dipole(23,r)*sqrt(3) + ((x**2 - 5*y**2)*Dipole(25,r) - (x - y)*(x + y)*(Dipole(26,r) + 2*Dipole(30,r) - Dipole(31,r)))*sqrt(3)) + (x**2 + y**2)**2*(-4*y**2*Dipole(24,r) + (x - y)*(x + y)*(Dipole(25,r) - Dipole(31,r))*sqrt(3))))/(4.*(x**2 + y**2)), \
	(2,1,1,-1,0,0) : lambda r,x,y,z,Dipole: x*y*z*(-2*Dipole(16,r) + Dipole(27,r)*sqrt(3)), \
	(2,1,1,-1,1,-1) : lambda r,x,y,z,Dipole: (x*z*((x**2 + 2*y**2)*Dipole(17,r) - y**2*(-(z**2*Dipole(20,r)) + 2*(x**2 + y**2)*(Dipole(18,r) + Dipole(28,r))) + (x**2 + y**2)*((x**2 + z**2)*Dipole(19,r) + y**2*Dipole(29,r))*sqrt(3)))/(x**2 + y**2), \
	(2,1,1,-1,1,0) : lambda r,x,y,z,Dipole: -(x*y*(Dipole(17,r) - (x**2 + y**2)*Dipole(28,r) + z**2*(2*Dipole(18,r) + Dipole(20,r) + Dipole(28,r) + (Dipole(19,r) - Dipole(29,r))*sqrt(3)))), \
	(2,1,1,-1,1,1) : lambda r,x,y,z,Dipole: (y*z*(-(y**2*Dipole(17,r)) - 2*x**2*(x**2 + y**2)*Dipole(18,r) + (x**2 + y**2)*(-x**2 + y**2 + z**2)*Dipole(28,r) + x**2*(z**2*Dipole(20,r) - (x**2 + y**2)*(Dipole(19,r) - Dipole(29,r))*sqrt(3))))/(x**2 + y**2), \
	(2,1,1,-1,2,-2) : lambda r,x,y,z,Dipole: (z*((x**6 + 3*x**4*y**2 + x**2*y**4 - y**6)*Dipole(21,r) + y**4*(-y**4 + z**4)*Dipole(25,r) + y**6*(y**2 + z**2)*Dipole(30,r) + x**8*(-Dipole(25,r) + Dipole(24,r)*sqrt(3)) + x**6*(z**2*Dipole(24,r)*sqrt(3) + y**2*(Dipole(23,r) - 3*Dipole(30,r) + 3*Dipole(31,r) + (-2*Dipole(22,r) + Dipole(24,r))*sqrt(3))) + x**2*y**2*(-2*z**4*Dipole(23,r) + y**2*z**2*(Dipole(23,r) + Dipole(25,r) + 2*(Dipole(26,r) + Dipole(30,r)) + Dipole(24,r)*sqrt(3)) + y**4*(Dipole(23,r) - Dipole(30,r) + 3*Dipole(31,r) - (2*Dipole(22,r) + Dipole(24,r))*sqrt(3))) + x**4*(z**4*Dipole(25,r) + y**2*z**2*(Dipole(23,r) + Dipole(25,r) + 2*Dipole(26,r) + Dipole(30,r) + 2*Dipole(24,r)*sqrt(3)) + y**4*(2*Dipole(23,r) + 2*Dipole(25,r) - 5*Dipole(30,r) + 6*Dipole(31,r) - (4*Dipole(22,r) + Dipole(24,r))*sqrt(3)))))/(x**2 + y**2)**2, \
	(2,1,1,-1,2,-1) : lambda r,x,y,z,Dipole: (x*(-((y**4 - 2*y**2*z**2 + x**2*(y - z)*(y + z))*Dipole(21,r)) + x**6*Dipole(25,r) + y**6*(-Dipole(25,r) + Dipole(30,r)) + y**2*z**4*(Dipole(23,r) + Dipole(26,r) + Dipole(24,r)*sqrt(3)) + x**4*(y**2*(Dipole(25,r) + Dipole(30,r)) + z**2*Dipole(24,r)*sqrt(3)) - y**4*z**2*(Dipole(23,r) + Dipole(26,r) + 3*Dipole(30,r) - 3*Dipole(31,r) + (2*Dipole(22,r) + Dipole(24,r))*sqrt(3)) - x**2*(y**4*(Dipole(25,r) - 2*Dipole(30,r)) + y**2*z**2*(Dipole(23,r) + Dipole(26,r) + 3*Dipole(30,r) - 3*Dipole(31,r) + 2*Dipole(22,r)*sqrt(3)) + z**4*(Dipole(25,r) - Dipole(24,r)*sqrt(3)))))/(x**2 + y**2), \
	(2,1,1,-1,2,0) : lambda r,x,y,z,Dipole: (x*y*z*(-(z**2*(4*Dipole(22,r) + 6*Dipole(24,r) + (2*Dipole(21,r) + Dipole(23,r) - Dipole(25,r) + 2*(Dipole(26,r) + Dipole(30,r) - Dipole(31,r)))*sqrt(3))) + (x**2 + y**2)*(2*Dipole(22,r) + (-2*Dipole(21,r) + Dipole(23,r) + 2*Dipole(30,r) - Dipole(31,r))*sqrt(3))))/2., \
	(2,1,1,-1,2,1) : lambda r,x,y,z,Dipole: (y*(-((x**4 + x**2*y**2 + y**2*z**2)*Dipole(21,r)) + x**6*(-Dipole(25,r) + Dipole(30,r)) + y**2*((y**4 - z**4)*Dipole(25,r) + z**2*(y**2 + z**2)*Dipole(30,r)) + x**2*(y**4*(Dipole(25,r) + Dipole(30,r)) + z**4*(Dipole(23,r) + Dipole(26,r) + Dipole(30,r)) - y**2*z**2*(Dipole(23,r) + Dipole(26,r) + Dipole(30,r) - 3*Dipole(31,r) + 2*(Dipole(22,r) + Dipole(24,r))*sqrt(3))) - x**4*(y**2*(Dipole(25,r) - 2*Dipole(30,r)) + z**2*(Dipole(23,r) + Dipole(26,r) + 2*Dipole(30,r) - 3*Dipole(31,r) + 2*(Dipole(22,r) + Dipole(24,r))*sqrt(3)))))/(x**2 + y**2), \
	(2,1,1,-1,2,2) : lambda r,x,y,z,Dipole: -(x*y*z*(2*(x**2 + y**2)*(x**2 + 3*y**2)*Dipole(21,r) - 2*y**2*z**4*(Dipole(23,r) + Dipole(25,r)) + y**6*(Dipole(23,r) + 4*Dipole(25,r) - 6*Dipole(30,r) + 3*Dipole(31,r) - 2*Dipole(22,r)*sqrt(3)) - x**6*(Dipole(23,r) + 4*Dipole(25,r) - 2*Dipole(30,r) + 3*Dipole(31,r) - 2*(Dipole(22,r) + 2*Dipole(24,r))*sqrt(3)) + y**4*z**2*(Dipole(23,r) + Dipole(25,r) + 2*(Dipole(26,r) - Dipole(30,r) + Dipole(24,r)*sqrt(3))) + x**2*(2*z**4*(Dipole(23,r) + Dipole(25,r)) - 4*y**2*z**2*(Dipole(30,r) - Dipole(24,r)*sqrt(3)) + y**4*(Dipole(23,r) + 4*Dipole(25,r) - 10*Dipole(30,r) + 3*Dipole(31,r) - 2*Dipole(22,r)*sqrt(3) + 4*Dipole(24,r)*sqrt(3))) + x**4*(-(y**2*(Dipole(23,r) + 4*Dipole(25,r) + 2*Dipole(30,r) + 3*Dipole(31,r) - 2*(Dipole(22,r) + 4*Dipole(24,r))*sqrt(3))) - z**2*(Dipole(23,r) + Dipole(25,r) + 2*(Dipole(26,r) + Dipole(30,r) - Dipole(24,r)*sqrt(3))))))/(2.*(x**2 + y**2)**2), \
	(2,1,1,0,0,0) : lambda r,x,y,z,Dipole: x*((x**2 + y**2 - z**2)*Dipole(16,r) + z**2*Dipole(27,r)*sqrt(3)), \
	(2,1,1,0,1,-1) : lambda r,x,y,z,Dipole: x*y*(-Dipole(17,r) + (x**2 + y**2 - z**2)*Dipole(18,r) - z**2*(Dipole(20,r) + 2*Dipole(28,r) + (Dipole(19,r) - Dipole(29,r))*sqrt(3))), \
	(2,1,1,0,1,0) : lambda r,x,y,z,Dipole: x*z*((x**2 + y**2 - z**2)*Dipole(18,r) + (x**2 + y**2)*Dipole(20,r) + (x**2 + y**2 - z**2)*Dipole(28,r) + ((x**2 + y**2)*Dipole(19,r) + z**2*Dipole(29,r))*sqrt(3)), \
	(2,1,1,0,1,1) : lambda r,x,y,z,Dipole: y**2*Dipole(17,r) + x**2*(x**2 + y**2 - z**2)*Dipole(18,r) + z**2*((y**2 + z**2)*Dipole(28,r) - x**2*(Dipole(20,r) + Dipole(28,r) + (Dipole(19,r) - Dipole(29,r))*sqrt(3))), \
	(2,1,1,0,2,-2) : lambda r,x,y,z,Dipole: (y*((-x**4 + y**4)*Dipole(21,r) - x**2*(x**2 + y**2 - z**2)*(x**2 + y**2 + 2*z**2)*Dipole(23,r) + z**4*(-2*y**2*Dipole(25,r) + (x**2 + y**2)*Dipole(30,r)) + x**2*(x**2 + y**2)*(x**2 + y**2 - z**2)*Dipole(22,r)*sqrt(3) - (x**2 + y**2)*z**2*(y**2*(2*Dipole(25,r) - Dipole(30,r)) + x**2*(-Dipole(25,r) + 2*Dipole(26,r) + 3*Dipole(30,r) - 3*Dipole(31,r) + 2*Dipole(24,r)*sqrt(3)))))/(x**2 + y**2), \
	(2,1,1,0,2,-1) : lambda r,x,y,z,Dipole: x*y*z*(-Dipole(21,r) + (x**2 + y**2 - z**2)*Dipole(23,r) - (2*(x**2 + y**2) + z**2)*Dipole(25,r) + x**2*Dipole(26,r) + y**2*Dipole(26,r) - z**2*Dipole(26,r) + x**2*Dipole(30,r) + y**2*Dipole(30,r) - 3*z**2*Dipole(30,r) + 3*z**2*Dipole(31,r) + (x**2 + y**2 - z**2)*(Dipole(22,r) + Dipole(24,r))*sqrt(3)), \
	(2,1,1,0,2,0) : lambda r,x,y,z,Dipole: (x*(-((x**2 + y**2 - 2*z**2)*(x**2 + y**2 - z**2)*Dipole(22,r)) + 6*(x**2 + y**2)*z**2*Dipole(24,r) - ((x**2 + y**2)*(x**2 + y**2 - z**2)*Dipole(23,r) + z**2*((x**2 + y**2)*Dipole(25,r) - 2*(x**2 + y**2)*Dipole(26,r) - 2*(x**2 + y**2 - z**2)*Dipole(30,r) + (x**2 + y**2 - 2*z**2)*Dipole(31,r)))*sqrt(3)))/2., \
	(2,1,1,0,2,1) : lambda r,x,y,z,Dipole: z*(y**2*Dipole(21,r) + (y**2 + z**2)*(2*y**2*Dipole(25,r) + z**2*Dipole(30,r)) + x**4*(Dipole(23,r) + Dipole(26,r) + Dipole(30,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3)) + x**2*(y**2*(Dipole(23,r) + 2*Dipole(25,r) + Dipole(26,r) + Dipole(30,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3)) - z**2*(Dipole(23,r) - Dipole(25,r) + Dipole(26,r) + 2*Dipole(30,r) - 3*Dipole(31,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3)))), \
	(2,1,1,0,2,2) : lambda r,x,y,z,Dipole: (x*(4*y**2*(x**2 + y**2)*Dipole(21,r) - (x - y)*(x + y)*(x**2 + y**2 - z**2)*(x**2 + y**2 + 2*z**2)*Dipole(23,r) + (x**4 - y**4)*(x**2 + y**2 - z**2)*Dipole(22,r)*sqrt(3) + z**2*(-((x**4 + 7*y**4 + 6*y**2*z**2 + 2*x**2*(4*y**2 + z**2))*Dipole(25,r)) + 2*(x**2 + y**2)*(-x**2 + 3*y**2 + z**2)*Dipole(30,r) - (x**4 - y**4)*(2*Dipole(26,r) - 3*Dipole(31,r) + 2*Dipole(24,r)*sqrt(3)))))/(2.*(x**2 + y**2)), \
	(2,1,1,1,0,0) : lambda r,x,y,z,Dipole: z*((-x**2 + y**2 + z**2)*Dipole(16,r) + x**2*Dipole(27,r)*sqrt(3)), \
	(2,1,1,1,1,-1) : lambda r,x,y,z,Dipole: (y*z*(-(y**2*Dipole(17,r)) + (-x**4 + y**4 + (x**2 + y**2)*z**2)*Dipole(18,r) + x**2*(z**2*Dipole(20,r) - (x**2 + y**2)*(2*Dipole(28,r) + (Dipole(19,r) - Dipole(29,r))*sqrt(3)))))/(x**2 + y**2), \
	(2,1,1,1,1,0) : lambda r,x,y,z,Dipole: y**2*Dipole(17,r) + z**2*(-x**2 + y**2 + z**2)*Dipole(18,r) + x**2*((x**2 + y**2)*Dipole(28,r) - z**2*(Dipole(20,r) + Dipole(28,r) + (Dipole(19,r) - Dipole(29,r))*sqrt(3))), \
	(2,1,1,1,1,1) : lambda r,x,y,z,Dipole: (x*z*(-(y**2*Dipole(17,r)) + (-x**4 + y**4 + (x**2 + y**2)*z**2)*Dipole(18,r) + (x**2 + y**2)*((-x**2 + y**2)*Dipole(28,r) + (y**2*Dipole(19,r) + x**2*Dipole(29,r))*sqrt(3)) + z**2*(x**2*Dipole(20,r) + (x**2 + y**2)*(Dipole(28,r) + Dipole(19,r)*sqrt(3)))))/(x**2 + y**2), \
	(2,1,1,1,2,-2) : lambda r,x,y,z,Dipole: (x*y*z*(-2*y**2*(x**2 + y**2)*Dipole(21,r) - 2*x**2*z**4*(Dipole(23,r) + Dipole(25,r)) + (x**2 + y**2)**2*(x**2*(Dipole(23,r) + Dipole(25,r) - 3*Dipole(30,r) + 3*Dipole(31,r) - (Dipole(22,r) + Dipole(24,r))*sqrt(3)) + y**2*(-3*Dipole(25,r) + Dipole(30,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3))) + (x**2 + y**2)*z**2*(y**2*(-3*Dipole(25,r) + Dipole(30,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3)) + x**2*(Dipole(23,r) - 2*Dipole(25,r) + 2*Dipole(26,r) + Dipole(30,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3)))))/(x**2 + y**2)**2, \
	(2,1,1,1,2,-1) : lambda r,x,y,z,Dipole: (y*(y**2*((x**2 + y**2)**2 - z**4)*Dipole(21,r) - x**2*(x**2 + y**2)**2*(2*Dipole(25,r) - Dipole(30,r)) + z**4*(x**2*Dipole(23,r) + (2*x**2 + y**2)*Dipole(25,r) + x**2*Dipole(26,r) + (x**2 + y**2)*Dipole(22,r)*sqrt(3)) + (x**2 + y**2)*z**2*((x**2 + y**2)*Dipole(25,r) + (-x**2 + y**2)*Dipole(22,r)*sqrt(3) - x**2*(Dipole(23,r) + Dipole(26,r) + 3*Dipole(30,r) - 3*Dipole(31,r) + 2*Dipole(24,r)*sqrt(3)))))/(x**2 + y**2), \
	(2,1,1,1,2,0) : lambda r,x,y,z,Dipole: (z*((x**2 + y**2 - 2*z**2)*(x**2 - y**2 - z**2)*Dipole(22,r) - 6*x**2*z**2*Dipole(24,r) + 2*y**2*Dipole(21,r)*sqrt(3) + (x**2*(x**2 + y**2 - z**2)*Dipole(23,r) - ((x**2 + y**2)**2 + y**2*z**2)*Dipole(25,r) + x**2*(-2*z**2*(Dipole(26,r) + Dipole(30,r) - Dipole(31,r)) + (x**2 + y**2)*(2*Dipole(30,r) - Dipole(31,r))))*sqrt(3)))/2., \
	(2,1,1,1,2,1) : lambda r,x,y,z,Dipole: (x*(y**2*((x**2 + y**2)**2 - z**4)*Dipole(21,r) + x**6*Dipole(30,r) + y**2*(y**2 + z**2)*(2*y**2*Dipole(25,r) + z**2*(-Dipole(25,r) + Dipole(30,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3))) + x**2*(y**4*(4*Dipole(25,r) + Dipole(30,r)) - y**2*z**2*(Dipole(23,r) - 2*Dipole(25,r) + Dipole(26,r) + Dipole(30,r) - 3*Dipole(31,r)) + z**4*(Dipole(23,r) + Dipole(26,r) + Dipole(30,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3))) - x**4*(-2*y**2*(Dipole(25,r) + Dipole(30,r)) + z**2*(Dipole(23,r) - Dipole(25,r) + Dipole(26,r) + 2*Dipole(30,r) - 3*Dipole(31,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3)))))/(x**2 + y**2), \
	(2,1,1,1,2,2) : lambda r,x,y,z,Dipole: (z*(2*y**2*(-x**4 + y**4)*Dipole(21,r) + x**8*(Dipole(23,r) - Dipole(25,r) - 2*Dipole(30,r) + 3*Dipole(31,r) - Dipole(22,r)*sqrt(3)) - y**4*(y**2 + z**2)*(-((y**2 + 2*z**2)*Dipole(25,r)) + y**2*Dipole(22,r)*sqrt(3)) + x**6*(y**2*(Dipole(23,r) - 10*Dipole(25,r) + 2*Dipole(30,r) + 3*Dipole(31,r) + 4*Dipole(24,r)*sqrt(3)) + z**2*(Dipole(23,r) - 2*Dipole(25,r) + 2*(Dipole(26,r) + Dipole(30,r)) + (Dipole(22,r) + 2*Dipole(24,r))*sqrt(3))) + x**4*(-((y**4 + 2*z**4)*Dipole(23,r)) + y**2*(-((16*y**2 + 3*z**2)*Dipole(25,r)) + 2*(5*y**2 + 2*z**2)*Dipole(30,r) - 3*y**2*Dipole(31,r) + (2*y**2 + z**2)*(Dipole(22,r) + 4*Dipole(24,r))*sqrt(3))) + x**2*y**2*(2*z**4*(Dipole(23,r) + 3*Dipole(25,r)) - y**4*(Dipole(23,r) + 6*Dipole(25,r) - 6*Dipole(30,r) + 3*Dipole(31,r) - 4*Dipole(24,r)*sqrt(3)) - y**2*z**2*(Dipole(23,r) + Dipole(22,r)*sqrt(3) - 2*(Dipole(25,r) - Dipole(26,r) + Dipole(30,r) + Dipole(24,r)*sqrt(3))))))/(2.*(x**2 + y**2)**2), \
	(2,2,1,-1,0,0) : lambda r,x,y,z,Dipole: -(y*(2*(2*x**2 + z**2)*Dipole(16,r) + (-x**2 + y**2)*Dipole(27,r)*sqrt(3)))/2., \
	(2,2,1,-1,1,-1) : lambda r,x,y,z,Dipole: (-(x**2*(x**4 - y**4 + 2*(x**2 + 3*y**2)*z**2)*Dipole(17,r)) - 2*(x**2*y + y**3)**2*(2*x**2 + z**2)*Dipole(18,r) + y**2*(-x + y)*(x + y)*z**2*(x**2 + y**2 + 2*z**2)*Dipole(20,r) - 2*(x**2*y + y**3)**2*(2*x**2 + z**2)*Dipole(28,r) + (x - y)*(x + y)*(x**2 + y**2)**2*((x**2 + z**2)*Dipole(19,r) + y**2*Dipole(29,r))*sqrt(3))/(2.*(x**2 + y**2)**2), \
	(2,2,1,-1,1,0) : lambda r,x,y,z,Dipole: (y*(4*x**2*z*Dipole(17,r) - 2*(x**2 + y**2)*(2*x**2*z + z**3)*Dipole(18,r) + 2*(x - y)*(x + y)*z**3*Dipole(20,r) + (x**4 - y**4)*z*(Dipole(20,r) - 2*Dipole(28,r) + (-Dipole(19,r) + Dipole(29,r))*sqrt(3))))/(2.*(x**2 + y**2)), \
	(2,2,1,-1,1,1) : lambda r,x,y,z,Dipole: (x*y*((x - y)*(x + y)*(x**2 + y**2 - 2*z**2)*Dipole(17,r) - 2*(x**2 + y**2)**2*(2*x**2 + z**2)*Dipole(18,r) + 2*(x**2 + y**2)**2*(2*y**2 + z**2)*Dipole(28,r) + (-x**2 + y**2)*(z**2*(x**2 + y**2 + 2*z**2)*Dipole(20,r) + (x**2 + y**2)**2*(Dipole(19,r) - Dipole(29,r))*sqrt(3))))/(2.*(x**2 + y**2)**2), \
	(2,2,1,-1,2,-2) : lambda r,x,y,z,Dipole: (x*(-(((x**2 - y**2)**2*(x**2 + y**2) + 2*(x**4 + 4*x**2*y**2 - y**4)*z**2)*Dipole(21,r)) + 2*y**4*z**4*(2*Dipole(23,r) + 3*Dipole(25,r) + 2*Dipole(26,r)) + x**8*Dipole(24,r)*sqrt(3) + y**8*(-Dipole(25,r) + 4*Dipole(30,r) - 3*Dipole(31,r) + Dipole(24,r)*sqrt(3)) + y**6*z**2*(2*(Dipole(23,r) + 3*Dipole(25,r) + Dipole(26,r)) - (2*Dipole(22,r) + Dipole(24,r))*sqrt(3)) + x**6*(y**2*(5*Dipole(25,r) - 4*Dipole(30,r) + 3*Dipole(31,r) - 4*Dipole(22,r)*sqrt(3)) + z**2*(2*Dipole(25,r) + Dipole(24,r)*sqrt(3))) + x**4*(2*z**4*Dipole(25,r) + y**2*z**2*(-2*Dipole(23,r) + 6*Dipole(25,r) - 2*Dipole(26,r) + (-2*Dipole(22,r) + Dipole(24,r))*sqrt(3)) + y**4*(9*Dipole(25,r) - 4*Dipole(30,r) + 3*Dipole(31,r) - 2*(4*Dipole(22,r) + Dipole(24,r))*sqrt(3))) + x**2*y**2*(-4*z**4*(Dipole(23,r) + Dipole(26,r)) + y**4*(3*Dipole(25,r) + 4*Dipole(30,r) - 3*Dipole(31,r) - 4*Dipole(22,r)*sqrt(3)) + y**2*z**2*(10*Dipole(25,r) - (4*Dipole(22,r) + Dipole(24,r))*sqrt(3)))))/(2.*(x**2 + y**2)**2), \
	(2,2,1,-1,2,-1) : lambda r,x,y,z,Dipole: -(z*(x**2*(x**4 - 5*y**4 + 6*y**2*z**2 + 2*x**2*(-2*y**2 + z**2))*Dipole(21,r) - 2*y**4*z**4*Dipole(26,r) + x**8*(2*Dipole(25,r) - Dipole(24,r)*sqrt(3)) - y**8*(Dipole(25,r) - Dipole(26,r) + 2*Dipole(30,r) - 3*Dipole(31,r) + Dipole(24,r)*sqrt(3)) + y**6*z**2*(2*Dipole(23,r) - 2*Dipole(25,r) + Dipole(26,r) + 2*Dipole(30,r) + (2*Dipole(22,r) + Dipole(24,r))*sqrt(3)) + x**6*(-(y**2*(Dipole(25,r) + Dipole(26,r) - 6*Dipole(30,r) + 3*Dipole(31,r) - 4*Dipole(22,r)*sqrt(3))) + z**2*(2*Dipole(25,r) - Dipole(24,r)*sqrt(3))) + x**4*y**2*(2*(4*y**2 + z**2)*Dipole(22,r)*sqrt(3) - y**2*(9*Dipole(25,r) + Dipole(26,r) - 10*Dipole(30,r) + 3*Dipole(31,r) - 2*Dipole(24,r)*sqrt(3)) - z**2*(2*Dipole(23,r) + 2*Dipole(25,r) + Dipole(26,r) - 2*Dipole(30,r) + Dipole(24,r)*sqrt(3))) + x**2*y**2*(2*z**4*Dipole(26,r) + y**4*(-7*Dipole(25,r) + Dipole(26,r) + 2*Dipole(30,r) + 3*Dipole(31,r)) + 4*y**2*(y**2 + z**2)*Dipole(22,r)*sqrt(3) + y**2*z**2*(-6*Dipole(25,r) + 4*Dipole(30,r) + Dipole(24,r)*sqrt(3)))))/(2.*(x**2 + y**2)**2), \
	(2,2,1,-1,2,0) : lambda r,x,y,z,Dipole: (y*(2*(x**2 + y**2)*(x**2 + y**2 - 2*z**2)*(2*x**2 + z**2)*Dipole(22,r) + 8*x**2*z**2*Dipole(21,r)*sqrt(3) - 4*y**2*z**4*Dipole(26,r)*sqrt(3) + x**6*(5*Dipole(25,r) - Dipole(31,r))*sqrt(3) + y**6*(-Dipole(25,r) + Dipole(31,r))*sqrt(3) + x**2*((3*y**4 + 4*y**2*z**2)*Dipole(25,r) + 4*z**4*Dipole(26,r) + y**4*Dipole(31,r))*sqrt(3) - 2*y**4*z**2*(-3*Dipole(24,r) + (-Dipole(23,r) + Dipole(25,r) + Dipole(26,r) - 2*Dipole(30,r) + Dipole(31,r))*sqrt(3)) + x**4*(y**2*(9*Dipole(25,r) - Dipole(31,r))*sqrt(3) + 2*z**2*(-3*Dipole(24,r) + (-Dipole(23,r) + 3*Dipole(25,r) + Dipole(26,r) - 2*Dipole(30,r) + Dipole(31,r))*sqrt(3)))))/(4.*(x**2 + y**2)), \
	(2,2,1,-1,2,1) : lambda r,x,y,z,Dipole: -(x*y*z*((-5*x**4 - 4*x**2*y**2 + y**4 + 2*(x - y)*(x + y)*z**2)*Dipole(21,r) - 2*y**2*z**4*Dipole(26,r) + y**4*z**2*(2*Dipole(23,r) + 8*Dipole(25,r) + Dipole(26,r) - 2*Dipole(30,r) + 2*Dipole(22,r)*sqrt(3)) + y**6*(9*Dipole(25,r) + Dipole(26,r) - 6*Dipole(30,r) + 3*Dipole(31,r) - 2*Dipole(24,r)*sqrt(3)) + x**6*(3*Dipole(25,r) - Dipole(26,r) + 2*Dipole(30,r) - 3*Dipole(31,r) + 2*(2*Dipole(22,r) + Dipole(24,r))*sqrt(3)) + x**2*(2*z**4*Dipole(26,r) + 4*y**2*z**2*(3*Dipole(25,r) - Dipole(30,r)) + 4*y**2*(y**2 + z**2)*Dipole(22,r)*sqrt(3) + y**4*(21*Dipole(25,r) + Dipole(26,r) - 10*Dipole(30,r) + 3*Dipole(31,r) - 2*Dipole(24,r)*sqrt(3))) + x**4*(-(z**2*(2*Dipole(23,r) - 4*Dipole(25,r) + Dipole(26,r) + 2*Dipole(30,r))) + 2*(4*y**2 + z**2)*Dipole(22,r)*sqrt(3) + y**2*(15*Dipole(25,r) - Dipole(26,r) - 2*Dipole(30,r) - 3*Dipole(31,r) + 2*Dipole(24,r)*sqrt(3)))))/(2.*(x**2 + y**2)**2), \
	(2,2,1,-1,2,2) : lambda r,x,y,z,Dipole: (y*(4*x**2*(x**4 - y**4 + 4*y**2*z**2)*Dipole(21,r) - 4*y**4*z**4*(Dipole(23,r) - Dipole(25,r) + Dipole(26,r)) + y**8*(Dipole(25,r) + 3*Dipole(31,r)) + x**8*(5*Dipole(25,r) + 3*Dipole(31,r) - 4*(Dipole(22,r) + Dipole(24,r))*sqrt(3)) + 2*y**6*z**2*(-Dipole(23,r) + 2*Dipole(25,r) - Dipole(26,r) + 2*Dipole(30,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3)) - 2*x**6*(z**2*Dipole(23,r) - 2*(y**2 + 2*z**2)*Dipole(25,r) + z**2*Dipole(26,r) - 8*y**2*Dipole(30,r) - 2*z**2*Dipole(30,r) + (2*y**2 + z**2)*(Dipole(22,r) + Dipole(24,r))*sqrt(3)) + 2*x**2*y**2*(4*z**4*(Dipole(23,r) + 3*Dipole(25,r) + Dipole(26,r)) + 2*y**4*(-Dipole(25,r) + 4*Dipole(30,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3)) + y**2*z**2*(Dipole(23,r) + 12*Dipole(25,r) + Dipole(26,r) + 6*Dipole(30,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3))) + 2*x**4*(-2*z**4*(Dipole(23,r) - Dipole(25,r) + Dipole(26,r)) + y**2*z**2*(Dipole(23,r) + 14*Dipole(25,r) + Dipole(26,r) + 6*Dipole(30,r) - (Dipole(22,r) + Dipole(24,r))*sqrt(3)) + y**4*(-3*Dipole(25,r) + 16*Dipole(30,r) - 3*Dipole(31,r) + 2*(Dipole(22,r) + Dipole(24,r))*sqrt(3)))))/(4.*(x**2 + y**2)**2), \
	(2,2,1,0,0,0) : lambda r,x,y,z,Dipole: -((x - y)*(x + y)*z*(2*Dipole(16,r) - Dipole(27,r)*sqrt(3)))/2., \
	(2,2,1,0,1,-1) : lambda r,x,y,z,Dipole: (y*z*(4*x**2*Dipole(17,r) + 2*(-x**4 + y**4)*Dipole(18,r) - 2*(x**2 + y**2)*(2*x**2 + z**2)*Dipole(28,r) + (x - y)*(x + y)*((x**2 + y**2 + 2*z**2)*Dipole(20,r) - (x**2 + y**2)*(Dipole(19,r) - Dipole(29,r))*sqrt(3))))/(2.*(x**2 + y**2)), \
	(2,2,1,0,1,0) : lambda r,x,y,z,Dipole: -((x - y)*(x + y)*((x**2 + y**2)*(Dipole(20,r) - Dipole(19,r)*sqrt(3)) + z**2*(2*(Dipole(18,r) + Dipole(20,r) + Dipole(28,r)) - Dipole(29,r)*sqrt(3))))/2., \
	(2,2,1,0,1,1) : lambda r,x,y,z,Dipole: (x*z*(-4*y**2*Dipole(17,r) + 2*(-x**4 + y**4)*Dipole(18,r) + 2*(x**2 + y**2)*(2*y**2 + z**2)*Dipole(28,r) + (x - y)*(x + y)*((x**2 + y**2 + 2*z**2)*Dipole(20,r) - (x**2 + y**2)*(Dipole(19,r) - Dipole(29,r))*sqrt(3))))/(2.*(x**2 + y**2)), \
	(2,2,1,0,2,-2) : lambda r,x,y,z,Dipole: (x*(x - y)*y*(x + y)*z*(4*Dipole(21,r) + 2*(x**2 + y**2 + 2*z**2)*Dipole(23,r) + (5*(x**2 + y**2) + 4*z**2)*Dipole(25,r) + 2*x**2*Dipole(26,r) + 2*y**2*Dipole(26,r) + 4*z**2*Dipole(26,r) - 4*x**2*Dipole(30,r) - 4*y**2*Dipole(30,r) + 3*x**2*Dipole(31,r) + 3*y**2*Dipole(31,r) - 2*(x**2 + y**2)*(Dipole(22,r) + Dipole(24,r))*sqrt(3)))/(2.*(x**2 + y**2)), \
	(2,2,1,0,2,-1) : lambda r,x,y,z,Dipole: -(y*(-4*x**2*z**2*Dipole(21,r) + 2*y**2*z**4*(-Dipole(25,r) + Dipole(26,r) + Dipole(30,r)) + x**6*(4*Dipole(25,r) + Dipole(26,r) - Dipole(24,r)*sqrt(3)) + y**6*(-Dipole(26,r) + Dipole(24,r)*sqrt(3)) - y**4*z**2*(2*Dipole(23,r) + Dipole(25,r) + Dipole(26,r) + 2*Dipole(30,r) - 3*Dipole(31,r) + (2*Dipole(22,r) + Dipole(24,r))*sqrt(3)) + x**2*(-2*z**4*(Dipole(25,r) + Dipole(26,r) - Dipole(30,r)) + 4*y**2*z**2*Dipole(30,r) + y**4*(4*Dipole(25,r) - Dipole(26,r) + Dipole(24,r)*sqrt(3))) + x**4*(y**2*(8*Dipole(25,r) + Dipole(26,r) - Dipole(24,r)*sqrt(3)) + z**2*(2*Dipole(23,r) + Dipole(25,r) + Dipole(26,r) + 6*Dipole(30,r) - 3*Dipole(31,r) + (2*Dipole(22,r) + Dipole(24,r))*sqrt(3)))))/(2.*(x**2 + y**2)), \
	(2,2,1,0,2,0) : lambda r,x,y,z,Dipole: ((x - y)*(x + y)*z*(2*(x**2 + y**2 - 2*z**2)*Dipole(22,r) + 6*(x**2 + y**2)*Dipole(24,r) + (2*(x**2 + y**2)*Dipole(23,r) + (x**2 + y**2 + 2*z**2)*Dipole(25,r) - 2*(x**2 + y**2 + 2*z**2)*Dipole(26,r) - 4*z**2*Dipole(30,r) - (x**2 + y**2 - 2*z**2)*Dipole(31,r))*sqrt(3)))/4., \
	(2,2,1,0,2,1) : lambda r,x,y,z,Dipole: (x*(-4*y**2*z**2*Dipole(21,r) - 2*y**2*z**4*(Dipole(25,r) + Dipole(26,r) - Dipole(30,r)) + y**6*(4*Dipole(25,r) + Dipole(26,r) - Dipole(24,r)*sqrt(3)) + x**6*(-Dipole(26,r) + Dipole(24,r)*sqrt(3)) + y**4*z**2*(2*Dipole(23,r) + Dipole(25,r) + Dipole(26,r) + 6*Dipole(30,r) - 3*Dipole(31,r) + (2*Dipole(22,r) + Dipole(24,r))*sqrt(3)) + x**2*(4*y**2*z**2*Dipole(30,r) + 2*z**4*(-Dipole(25,r) + Dipole(26,r) + Dipole(30,r)) + y**4*(8*Dipole(25,r) + Dipole(26,r) - Dipole(24,r)*sqrt(3))) - x**4*(y**2*(-4*Dipole(25,r) + Dipole(26,r) - Dipole(24,r)*sqrt(3)) + z**2*(2*Dipole(23,r) + Dipole(25,r) + Dipole(26,r) + 2*Dipole(30,r) - 3*Dipole(31,r) + (2*Dipole(22,r) + Dipole(24,r))*sqrt(3)))))/(2.*(x**2 + y**2)), \
	(2,2,1,0,2,2) : lambda r,x,y,z,Dipole: (z*(-16*x**2*y**2*Dipole(21,r) + 4*y**2*z**4*Dipole(25,r) + 4*y**4*z**2*(Dipole(23,r) + Dipole(25,r) + Dipole(26,r) + Dipole(30,r)) + x**6*(2*Dipole(23,r) + Dipole(25,r) + 2*Dipole(26,r) + 3*Dipole(31,r) - 2*(Dipole(22,r) + Dipole(24,r))*sqrt(3)) + y**6*(2*Dipole(23,r) + Dipole(25,r) + 2*Dipole(26,r) + 3*Dipole(31,r) - 2*(Dipole(22,r) + Dipole(24,r))*sqrt(3)) + x**4*(4*z**2*(Dipole(23,r) + Dipole(25,r) + Dipole(26,r) + Dipole(30,r)) + y**2*(-2*Dipole(23,r) - 17*Dipole(25,r) - 2*Dipole(26,r) + 16*Dipole(30,r) - 3*Dipole(31,r) + 2*(Dipole(22,r) + Dipole(24,r))*sqrt(3))) + x**2*(4*z**4*Dipole(25,r) - 8*y**2*z**2*(Dipole(23,r) + Dipole(25,r) + Dipole(26,r) - Dipole(30,r)) + y**4*(-2*Dipole(23,r) - 17*Dipole(25,r) - 2*Dipole(26,r) + 16*Dipole(30,r) - 3*Dipole(31,r) + 2*(Dipole(22,r) + Dipole(24,r))*sqrt(3)))))/(4.*(x**2 + y**2)), \
	(2,2,1,1,0,0) : lambda r,x,y,z,Dipole: x*(2*y**2 + z**2)*Dipole(16,r) + (x*(x - y)*(x + y)*Dipole(27,r)*sqrt(3))/2., \
	(2,2,1,1,1,-1) : lambda r,x,y,z,Dipole: (x*y*((x - y)*(x + y)*(x**2 + y**2 - 2*z**2)*Dipole(17,r) + 2*(x**2 + y**2)**2*(2*y**2 + z**2)*Dipole(18,r) - 2*(x**2 + y**2)**2*(2*x**2 + z**2)*Dipole(28,r) + (-x**2 + y**2)*(z**2*(x**2 + y**2 + 2*z**2)*Dipole(20,r) + (x**2 + y**2)**2*(Dipole(19,r) - Dipole(29,r))*sqrt(3))))/(2.*(x**2 + y**2)**2), \
	(2,2,1,1,1,0) : lambda r,x,y,z,Dipole: (x*z*(-4*y**2*Dipole(17,r) + 2*(x**2 + y**2)*(2*y**2 + z**2)*Dipole(18,r) + 2*(x - y)*(x + y)*z**2*Dipole(20,r) + (x**4 - y**4)*(Dipole(20,r) - 2*Dipole(28,r) + (-Dipole(19,r) + Dipole(29,r))*sqrt(3))))/(2.*(x**2 + y**2)), \
	(2,2,1,1,1,1) : lambda r,x,y,z,Dipole: (y**2*(-x**4 + y**4 + 2*(3*x**2 + y**2)*z**2)*Dipole(17,r) + 2*(x**3 + x*y**2)**2*(2*y**2 + z**2)*Dipole(18,r) - x**2*(x - y)*(x + y)*z**2*(x**2 + y**2 + 2*z**2)*Dipole(20,r) + 2*(x**3 + x*y**2)**2*(2*y**2 + z**2)*Dipole(28,r) + (x - y)*(x + y)*(x**2 + y**2)**2*((y**2 + z**2)*Dipole(19,r) + x**2*Dipole(29,r))*sqrt(3))/(2.*(x**2 + y**2)**2), \
	(2,2,1,1,2,-2) : lambda r,x,y,z,Dipole: (y*(((x**2 - y**2)**2*(x**2 + y**2) + 2*(-x**4 + 4*x**2*y**2 + y**4)*z**2)*Dipole(21,r) - 2*z**4*(2*x**2*(x - y)*(x + y)*Dipole(23,r) + (3*x**4 + y**4)*Dipole(25,r) + 2*x**2*(x - y)*(x + y)*Dipole(26,r)) - (x**2 + y**2)**2*(x**2*(-((x**2 - 5*y**2)*Dipole(25,r)) + (x - y)*(x + y)*(4*Dipole(30,r) - 3*Dipole(31,r))) - 4*x**2*y**2*Dipole(22,r)*sqrt(3) + (x**2 - y**2)**2*Dipole(24,r)*sqrt(3)) + (x**2 + y**2)*z**2*(2*x**2*y**2*(Dipole(23,r) - 2*Dipole(25,r) + Dipole(26,r) + Dipole(22,r)*sqrt(3)) - y**4*(2*Dipole(25,r) + Dipole(24,r)*sqrt(3)) + x**4*(-2*(Dipole(23,r) + 3*Dipole(25,r) + Dipole(26,r)) + (2*Dipole(22,r) + Dipole(24,r))*sqrt(3)))))/(2.*(x**2 + y**2)**2), \
	(2,2,1,1,2,-1) : lambda r,x,y,z,Dipole: (x*y*z*((x**4 - 5*y**4 + 2*y**2*z**2 - 2*x**2*(2*y**2 + z**2))*Dipole(21,r) + 2*y**2*z**4*Dipole(26,r) - y**4*z**2*(2*Dipole(23,r) - 4*Dipole(25,r) + Dipole(26,r) + 2*Dipole(30,r) - 2*Dipole(22,r)*sqrt(3)) + x**6*(9*Dipole(25,r) + Dipole(26,r) - 6*Dipole(30,r) + 3*Dipole(31,r) - 2*Dipole(24,r)*sqrt(3)) + y**6*(3*Dipole(25,r) - Dipole(26,r) + 2*Dipole(30,r) - 3*Dipole(31,r) + 2*(2*Dipole(22,r) + Dipole(24,r))*sqrt(3)) + x**4*(z**2*(2*Dipole(23,r) + 8*Dipole(25,r) + Dipole(26,r) - 2*Dipole(30,r)) + 2*(2*y**2 + z**2)*Dipole(22,r)*sqrt(3) + y**2*(21*Dipole(25,r) + Dipole(26,r) - 10*Dipole(30,r) + 3*Dipole(31,r) - 2*Dipole(24,r)*sqrt(3))) + x**2*(-2*z**4*Dipole(26,r) + 4*y**2*z**2*(3*Dipole(25,r) - Dipole(30,r) + Dipole(22,r)*sqrt(3)) + y**4*(15*Dipole(25,r) - Dipole(26,r) - 2*Dipole(30,r) - 3*Dipole(31,r) + 2*(4*Dipole(22,r) + Dipole(24,r))*sqrt(3)))))/(2.*(x**2 + y**2)**2), \
	(2,2,1,1,2,0) : lambda r,x,y,z,Dipole: (x*(-2*(x**2 + y**2)*(x**2 + y**2 - 2*z**2)*(2*y**2 + z**2)*Dipole(22,r) - 8*y**2*z**2*Dipole(21,r)*sqrt(3) - 4*y**2*z**4*Dipole(26,r)*sqrt(3) + x**6*(Dipole(25,r) - Dipole(31,r))*sqrt(3) + y**6*(-5*Dipole(25,r) + Dipole(31,r))*sqrt(3) + x**2*(-((9*y**4 + 4*y**2*z**2)*Dipole(25,r)) + 4*z**4*Dipole(26,r) + y**4*Dipole(31,r))*sqrt(3) + 2*y**4*z**2*(3*Dipole(24,r) + Dipole(23,r)*sqrt(3) - (3*Dipole(25,r) + Dipole(26,r) - 2*Dipole(30,r) + Dipole(31,r))*sqrt(3)) - x**4*(y**2*(3*Dipole(25,r) + Dipole(31,r))*sqrt(3) + 2*z**2*(3*Dipole(24,r) + Dipole(23,r)*sqrt(3) - (Dipole(25,r) + Dipole(26,r) - 2*Dipole(30,r) + Dipole(31,r))*sqrt(3)))))/(4.*(x**2 + y**2)), \
	(2,2,1,1,2,1) : lambda r,x,y,z,Dipole: (z*(y**2*(-5*x**4 - 4*x**2*y**2 + y**4 + 2*(3*x**2 + y**2)*z**2)*Dipole(21,r) + y**6*(y**2 + z**2)*(2*Dipole(25,r) - Dipole(24,r)*sqrt(3)) - x**8*(Dipole(25,r) - Dipole(26,r) + 2*Dipole(30,r) - 3*Dipole(31,r) + Dipole(24,r)*sqrt(3)) + x**2*y**2*(2*z**4*Dipole(26,r) - y**4*(Dipole(25,r) + Dipole(26,r) - 6*Dipole(30,r) + 3*Dipole(31,r) - 4*Dipole(22,r)*sqrt(3)) - y**2*z**2*(2*Dipole(23,r) + 2*Dipole(25,r) + Dipole(26,r) - 2*Dipole(30,r) + (-2*Dipole(22,r) + Dipole(24,r))*sqrt(3))) + x**6*(y**2*(-7*Dipole(25,r) + Dipole(26,r) + 2*Dipole(30,r) + 3*Dipole(31,r) + 4*Dipole(22,r)*sqrt(3)) + z**2*(2*Dipole(23,r) - 2*Dipole(25,r) + Dipole(26,r) + 2*Dipole(30,r) + (2*Dipole(22,r) + Dipole(24,r))*sqrt(3))) + x**4*(-3*(3*y**4 + 2*y**2*z**2)*Dipole(25,r) - (y**4 + 2*z**4)*Dipole(26,r) + y**2*(2*(5*y**2 + 2*z**2)*Dipole(30,r) - 3*y**2*Dipole(31,r) + (2*y**2 + z**2)*(4*Dipole(22,r) + Dipole(24,r))*sqrt(3)))))/(2.*(x**2 + y**2)**2), \
	(2,2,1,1,2,2) : lambda r,x,y,z,Dipole: (x*(4*y**2*(-x**4 + y**4 + 4*x**2*z**2)*Dipole(21,r) - 4*y**4*z**4*(Dipole(23,r) - Dipole(25,r) + Dipole(26,r)) + x**8*(Dipole(25,r) + 3*Dipole(31,r)) + y**8*(5*Dipole(25,r) + 3*Dipole(31,r) - 4*(Dipole(22,r) + Dipole(24,r))*sqrt(3)) - 2*y**6*z**2*(Dipole(23,r) - 4*Dipole(25,r) + Dipole(26,r) - 2*Dipole(30,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3)) + 2*x**6*(-(z**2*Dipole(23,r)) + 2*(-y**2 + z**2)*Dipole(25,r) - z**2*Dipole(26,r) + 8*y**2*Dipole(30,r) + 2*z**2*Dipole(30,r) + (2*y**2 + z**2)*(Dipole(22,r) + Dipole(24,r))*sqrt(3)) + 2*x**2*y**2*(4*z**4*(Dipole(23,r) + 3*Dipole(25,r) + Dipole(26,r)) + 2*y**4*(Dipole(25,r) + 4*Dipole(30,r) - (Dipole(22,r) + Dipole(24,r))*sqrt(3)) + y**2*z**2*(Dipole(23,r) + 14*Dipole(25,r) + Dipole(26,r) + 6*Dipole(30,r) - (Dipole(22,r) + Dipole(24,r))*sqrt(3))) + 2*x**4*(-2*z**4*(Dipole(23,r) - Dipole(25,r) + Dipole(26,r)) + y**2*z**2*(Dipole(23,r) + 12*Dipole(25,r) + Dipole(26,r) + 6*Dipole(30,r) + (Dipole(22,r) + Dipole(24,r))*sqrt(3)) + y**4*(-3*Dipole(25,r) + 16*Dipole(30,r) - 3*Dipole(31,r) + 2*(Dipole(22,r) + Dipole(24,r))*sqrt(3)))))/(4.*(x**2 + y**2)**2), \
}

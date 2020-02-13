# This file has been generated automatically
from numpy import sin,cos,sqrt,pi
phi = [ \
]

angular_phi = {\
	(0,0,0,0,{0, 0, 0, 0}[[5]],{0, 0, 0, 0}[[6]]) : phi[0],	\
	(0,0,1,0,{0, 0, 1, 0}[[5]],{0, 0, 1, 0}[[6]]) : phi[2],	\
	(0,0,2,0,{0, 0, 2, 0}[[5]],{0, 0, 2, 0}[[6]]) : phi[3],	\
	(1,0,0,0,{1, 0, 0, 0}[[5]],{1, 0, 0, 0}[[6]]) : phi[4],	\
	(1,-1,1,-1,{1, -1, 1, -1}[[5]],{1, -1, 1, -1}[[6]]) : phi[5],	\
	(1,0,1,0,{1, 0, 1, 0}[[5]],{1, 0, 1, 0}[[6]]) : phi[6],	\
	(1,1,1,1,{1, 1, 1, 1}[[5]],{1, 1, 1, 1}[[6]]) : phi[5],	\
	(1,-1,2,-1,{1, -1, 2, -1}[[5]],{1, -1, 2, -1}[[6]]) : phi[7],	\
	(1,0,2,0,{1, 0, 2, 0}[[5]],{1, 0, 2, 0}[[6]]) : phi[8],	\
	(1,1,2,1,{1, 1, 2, 1}[[5]],{1, 1, 2, 1}[[6]]) : phi[7],	\
	(2,0,0,0,{2, 0, 0, 0}[[5]],{2, 0, 0, 0}[[6]]) : phi[9],	\
	(2,-1,1,-1,{2, -1, 1, -1}[[5]],{2, -1, 1, -1}[[6]]) : phi[10],	\
	(2,0,1,0,{2, 0, 1, 0}[[5]],{2, 0, 1, 0}[[6]]) : phi[11],	\
	(2,1,1,1,{2, 1, 1, 1}[[5]],{2, 1, 1, 1}[[6]]) : phi[10],	\
	(2,-2,2,-2,{2, -2, 2, -2}[[5]],{2, -2, 2, -2}[[6]]) : phi[12],	\
	(2,-1,2,-1,{2, -1, 2, -1}[[5]],{2, -1, 2, -1}[[6]]) : phi[13],	\
	(2,0,2,0,{2, 0, 2, 0}[[5]],{2, 0, 2, 0}[[6]]) : phi[14],	\
	(2,1,2,1,{2, 1, 2, 1}[[5]],{2, 1, 2, 1}[[6]]) : phi[13],	\
	(2,2,2,2,{2, 2, 2, 2}[[5]],{2, 2, 2, 2}[[6]]) : phi[12],	\
}

tau2index = {\
	(0,0,0,0,{0, 0, 0, 0}[[5]],{0, 0, 0, 0}[[6]]) : 0,\
	(0,0,1,0,{0, 0, 1, 0}[[5]],{0, 0, 1, 0}[[6]]) : 2,\
	(0,0,2,0,{0, 0, 2, 0}[[5]],{0, 0, 2, 0}[[6]]) : 3,\
	(1,0,0,0,{1, 0, 0, 0}[[5]],{1, 0, 0, 0}[[6]]) : 4,\
	(1,-1,1,-1,{1, -1, 1, -1}[[5]],{1, -1, 1, -1}[[6]]) : 5,\
	(1,0,1,0,{1, 0, 1, 0}[[5]],{1, 0, 1, 0}[[6]]) : 6,\
	(1,1,1,1,{1, 1, 1, 1}[[5]],{1, 1, 1, 1}[[6]]) : 5,\
	(1,-1,2,-1,{1, -1, 2, -1}[[5]],{1, -1, 2, -1}[[6]]) : 7,\
	(1,0,2,0,{1, 0, 2, 0}[[5]],{1, 0, 2, 0}[[6]]) : 8,\
	(1,1,2,1,{1, 1, 2, 1}[[5]],{1, 1, 2, 1}[[6]]) : 7,\
	(2,0,0,0,{2, 0, 0, 0}[[5]],{2, 0, 0, 0}[[6]]) : 9,\
	(2,-1,1,-1,{2, -1, 1, -1}[[5]],{2, -1, 1, -1}[[6]]) : 10,\
	(2,0,1,0,{2, 0, 1, 0}[[5]],{2, 0, 1, 0}[[6]]) : 11,\
	(2,1,1,1,{2, 1, 1, 1}[[5]],{2, 1, 1, 1}[[6]]) : 10,\
	(2,-2,2,-2,{2, -2, 2, -2}[[5]],{2, -2, 2, -2}[[6]]) : 12,\
	(2,-1,2,-1,{2, -1, 2, -1}[[5]],{2, -1, 2, -1}[[6]]) : 13,\
	(2,0,2,0,{2, 0, 2, 0}[[5]],{2, 0, 2, 0}[[6]]) : 14,\
	(2,1,2,1,{2, 1, 2, 1}[[5]],{2, 1, 2, 1}[[6]]) : 13,\
	(2,2,2,2,{2, 2, 2, 2}[[5]],{2, 2, 2, 2}[[6]]) : 12,\
}

index2tau = {\
	0 : (0,0,0,0,{0, 0, 0, 0}[[5]],{0, 0, 0, 0}[[6]]),\
	2 : (0,0,1,0,{0, 0, 1, 0}[[5]],{0, 0, 1, 0}[[6]]),\
	3 : (0,0,2,0,{0, 0, 2, 0}[[5]],{0, 0, 2, 0}[[6]]),\
	5 : (1,-1,1,-1,{1, -1, 1, -1}[[5]],{1, -1, 1, -1}[[6]]),\
	7 : (1,-1,2,-1,{1, -1, 2, -1}[[5]],{1, -1, 2, -1}[[6]]),\
	4 : (1,0,0,0,{1, 0, 0, 0}[[5]],{1, 0, 0, 0}[[6]]),\
	6 : (1,0,1,0,{1, 0, 1, 0}[[5]],{1, 0, 1, 0}[[6]]),\
	8 : (1,0,2,0,{1, 0, 2, 0}[[5]],{1, 0, 2, 0}[[6]]),\
	12 : (2,-2,2,-2,{2, -2, 2, -2}[[5]],{2, -2, 2, -2}[[6]]),\
	10 : (2,-1,1,-1,{2, -1, 1, -1}[[5]],{2, -1, 1, -1}[[6]]),\
	13 : (2,-1,2,-1,{2, -1, 2, -1}[[5]],{2, -1, 2, -1}[[6]]),\
	9 : (2,0,0,0,{2, 0, 0, 0}[[5]],{2, 0, 0, 0}[[6]]),\
	11 : (2,0,1,0,{2, 0, 1, 0}[[5]],{2, 0, 1, 0}[[6]]),\
	14 : (2,0,2,0,{2, 0, 2, 0}[[5]],{2, 0, 2, 0}[[6]]),\
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
